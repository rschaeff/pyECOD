"""
Batch processing commands for the ECOD pipeline
"""

import argparse
import logging
import os
from typing import Dict, Any
from datetime import datetime

from .base_command import BaseCommand
from ecod.models import Batch


logger = logging.getLogger("ecod.cli.batch")

# Define commands in this group
COMMANDS = {
    'create': 'Create a new processing batch',
    'status': 'Check batch status',
    'list': 'List all batches',
    'run': 'Run a full pipeline on a batch'
}

class BatchCommand(BaseCommand):
    def __init__(self, config_path: Optional[str] = None):
        # Call parent initializer to set up logger, db, config, etc.
        super().__init__(config_path)

    def setup_parser(self, parser: argparse.ArgumentParser) -> None:
        """Set up the argument parser for batch commands"""
        subparsers = parser.add_subparsers(dest='command', help='Batch command')
        
        # Create command
        create_parser = subparsers.add_parser('create', help=COMMANDS['create'])
        create_parser.add_argument('--limit', type=int, default=10,
                                help='Maximum number of proteins to include')
        create_parser.add_argument('--type', type=str, default='full',
                                help='Type of batch (blast, hhsearch, full)')
        create_parser.add_argument('--output-dir', type=str,
                                help='Output directory override')
        
        # Status command
        status_parser = subparsers.add_parser('status', help=COMMANDS['status'])
        status_parser.add_argument('--batch-id', type=int, required=True,
                                help='Batch ID to check')
        
        # List command
        list_parser = subparsers.add_parser('list', help=COMMANDS['list'])
        list_parser.add_argument('--limit', type=int, default=10,
                              help='Maximum batches to list')
        list_parser.add_argument('--status', type=str,
                              help='Filter by status (created, processing, completed)')
        
        # Run command
        run_parser = subparsers.add_parser('run', help=COMMANDS['run'])
        run_parser.add_argument('--batch-id', type=int, required=True,
                             help='Batch ID to run')
        run_parser.add_argument('--blast-only', action='store_true',
                             help='Run only BLAST pipeline')
        run_parser.add_argument('--hhsearch-only', action='store_true',
                             help='Run only HHSearch pipeline')
        run_parser.add_argument('--domain-only', action='store_true',
                             help='Run only domain analysis')
        run_parser.add_argument('--partitioned', action='store_true',
                     help='Use adaptive processing paths based on BLAST confidence')

    def run_command(self, args: argparse.Namespace) -> int:
        """Run the specified batch command"""
        # Handle different commands
        if args.command == 'create':
            return self._create_batch(args)
        elif args.command == 'status':
            return self._check_status(args)
        elif args.command == 'list':
            return self._list_batches(args)
        elif args.command == 'run':
            return self._run_batch(args)
        else:
            self.logger.error(f"Unknown command: {args.command}")
            return 1

    def _create_batch(self, args: argparse.Namespace) -> int:
        """Create a new processing batch"""
        try:
            config = self.context.config.copy()
            if args.output_dir:
                config['paths'] = config.get('paths', {})
                config['paths']['output_dir'] = args.output_dir
            
            # Get unprocessed proteins
            query = """
            SELECT 
                p.id, p.pdb_id, p.chain_id, p.source_id, p.length, ps.sequence
            FROM 
                ecod_schema.protein p
            JOIN
                ecod_schema.protein_sequence ps ON p.id = ps.protein_id
            LEFT JOIN (
                SELECT DISTINCT protein_id 
                FROM ecod_schema.process_status 
                WHERE status IN ('success', 'completed')
            ) ps_done ON p.id = ps_done.protein_id
            WHERE 
                ps_done.protein_id IS NULL
                AND ps.sequence IS NOT NULL
            ORDER BY 
                p.id
            LIMIT %s
            """
            
            proteins = db.execute_dict_query(query, (args.limit,))
            
            if not proteins:
                logger.error("No unprocessed proteins found")
                return 1
            
            logger.info(f"Found {len(proteins)} unprocessed proteins")
            
            # Generate batch name with timestamp
            timestamp = datetime.now().strftime("%Y%m%d_%H%M")
            batch_name = f"{args.type}_batch_{timestamp}"
            
            # Get base output directory from config
            base_dir = config.get('paths', {}).get('output_dir', './output')
            
            # Create batch directory
            batch_path = os.path.join(base_dir, batch_name)
            os.makedirs(batch_path, exist_ok=True)
            
            # Create subdirectories
            os.makedirs(os.path.join(batch_path, "query_fastas"), exist_ok=True)
            if args.type in ['blast', 'full']:
                os.makedirs(os.path.join(batch_path, "chain_blast_results"), exist_ok=True)
                os.makedirs(os.path.join(batch_path, "domain_blast_results"), exist_ok=True)
            if args.type in ['hhsearch', 'full']:
                os.makedirs(os.path.join(batch_path, "ecod_dump"), exist_ok=True)
            
            # Insert batch record
            batch_id = db.insert(
                "ecod_schema.batch",
                {
                    "batch_name": batch_name,
                    "base_path": batch_path,
                    "type": args.type,
                    "ref_version": config.get('reference', {}).get('current_version', 'develop291'),
                    "total_items": len(proteins),
                    "status": "created"
                },
                "id"
            )
            
            # Register proteins in this batch
            for protein in proteins:
                # Determine relative path for this protein
                pdb_id = protein['pdb_id']
                chain_id = protein['chain_id']
                rel_path = f"{pdb_id}_{chain_id}"
                
                # Register in process_status
                process_id = db.insert(
                    "ecod_schema.process_status",
                    {
                        "protein_id": protein['id'],
                        "batch_id": batch_id,
                        "current_stage": "fasta_generated",
                        "status": "pending",
                        "relative_path": rel_path
                    },
                    "id"
                )
                
                # Generate FASTA file
                fasta_dir = os.path.join(batch_path, "query_fastas")
                fasta_path = os.path.join(fasta_dir, f"{protein['source_id']}.fa")
                with open(fasta_path, 'w') as f:
                    f.write(f">{protein['source_id']}\n{protein['sequence']}\n")
                
                # Register FASTA file
                db.insert(
                    "ecod_schema.process_file",
                    {
                        "process_id": process_id,
                        "file_type": "fasta",
                        "file_path": f"query_fastas/{protein['source_id']}.fa",
                        "file_exists": True,
                        "file_size": os.path.getsize(fasta_path)
                    }
                )
        
            self.logger.info(f"Created batch {batch_name} with ID {batch_id}")
            print(f"Created batch {batch_id}: {batch_name}")
            print(f"  - Type: {args.type}")
            print(f"  - Items: {len(proteins)}")
            print(f"  - Path: {batch_path}")
        
            return 0
        except Exception as e:
            from ecod.error_handlers import log_exception
            from ecod.exceptions import PipelineError, FileOperationError
        
            # Provide context for the error
            context = {
                "command": "create_batch",
                "output_dir": args.output_dir,
                "type": args.type,
                "limit": args.limit
            }
            
            log_exception(self.logger, e, context=context)
        
            # Determine appropriate exception type
            if "path" in str(e).lower() or "file" in str(e).lower() or "directory" in str(e).lower():
                raise FileOperationError(f"Error creating batch: {str(e)}", context) from e
            else:
                raise PipelineError(f"Error creating batch: {str(e)}", context) from e

    def _check_status(self, args: argparse.Namespace) -> int:
        """Check batch status"""
        try:
            # Get batch information
            query = """
            SELECT 
                b.id, b.batch_name, b.type, b.ref_version, b.total_items, 
                b.completed_items, b.status, b.created_at, b.completed_at
            FROM 
                ecod_schema.batch b
            WHERE 
                b.id = %s
            """
            
            rows = db.execute_dict_query(query, (args.batch_id,))
            
            if not rows:
                logger.error(f"Batch {args.batch_id} not found")
                return 1
            
            batch = rows[0]
            
            # Get process status breakdown
            query = """
            SELECT 
                status, current_stage, COUNT(*) as count
            FROM 
                ecod_schema.process_status
            WHERE 
                batch_id = %s
            GROUP BY 
                status, current_stage
            ORDER BY 
                count DESC
            """
            
            status_rows = db.execute_dict_query(query, (args.batch_id,))
            
            # Display batch information
            print(f"Batch {batch['id']}: {batch['batch_name']}")
            print(f"  - Type: {batch['type']}")
            print(f"  - Reference: {batch['ref_version']}")
            print(f"  - Status: {batch['status']}")
            print(f"  - Progress: {batch['completed_items']}/{batch['total_items']} ({batch['completed_items']/batch['total_items']*100:.1f}%)")
            print(f"  - Created: {batch['created_at']}")
            if batch['completed_at']:
                print(f"  - Completed: {batch['completed_at']}")
            
            # Display status breakdown
            print("\nProcess Status Breakdown:")
            for row in status_rows:
                print(f"  - {row['status']} - {row['current_stage']}: {row['count']}")
            
            return 0
        except Exception as e:
            from ecod.error_handlers import log_exception
            from ecod.exceptions import PipelineError, FileOperationError

            # Provide context for the error
            context = {
                "command": "batch_status",
                "output_dir": args.output_dir,
                "type": args.type,
                "limit": args.limit
            }

            #First try at this
            raise PipelineError(f"Error checking batch status: {str(e)}", context) from e
        
        log_exception(self.logger, e, context=context)

    def _list_batches(self, args: argparse.Namespace) -> int:
        """List all batches"""
        # Build query
        query = """
        SELECT 
            b.id, b.batch_name, b.type, b.ref_version, b.total_items, 
            b.completed_items, b.status, b.created_at, b.completed_at
        FROM 
            ecod_schema.batch b
        """
        
        params = []
        
        if args.status:
            query += " WHERE b.status = %s"
            params.append(args.status)
        
        query += " ORDER BY b.id DESC LIMIT %s"
        params.append(args.limit)
        
        rows = db.execute_dict_query(query, tuple(params))
        
        if not rows:
            print("No batches found")
            return 0
        
        # Display batch list
        print(f"Found {len(rows)} batches:")
        for batch in rows:
            progress = f"{batch['completed_items']}/{batch['total_items']}"
            percent = batch['completed_items']/batch['total_items']*100 if batch['total_items'] > 0 else 0
            print(f"  {batch['id']} - {batch['batch_name']} ({batch['type']}): {progress} ({percent:.1f}%), Status: {batch['status']}")
        
        return 0

    def _run_batch(self, args: argparse.Namespace) -> int:
        """Run a batch through the pipeline"""
        try:
            # Import orchestrator
            from ecod.pipelines.orchestrator import PipelineOrchestrator
            
            # Initialize orchestrator using the context's config path
            orchestrator = PipelineOrchestrator(self.context.config_manager.config_path)
            
            # Determine what to run
            if args.blast_only:
                self.logger.info(f"Running BLAST pipeline for batch {args.batch_id}")
                success = orchestrator.run_blast_pipeline(args.batch_id)
            elif args.hhsearch_only:
                self.logger.info(f"Running HHSearch pipeline for batch {args.batch_id}")
                success = orchestrator.run_hhsearch_pipeline(args.batch_id)
            elif args.domain_only:
                self.logger.info(f"Running domain analysis for batch {args.batch_id}")
                success = orchestrator.run_domain_analysis(args.batch_id)
            elif args.partitioned:
                self.logger.info(f"Running partitioned pipeline for batch {args.batch_id}")
                result = orchestrator.run_partitioned_pipeline(args.batch_id)
                success = result.get("status") == "completed"
                
                # Log additional information about partitioning
                if "paths" in result:
                    self.logger.info(f"Processed {result['paths'].get('blast_only', 0)} proteins with BLAST-only path")
                    self.logger.info(f"Processed {result['paths'].get('full_pipeline', 0)} proteins with full pipeline")
            else:
                # Run full pipeline
                self.logger.info(f"Running full pipeline for batch {args.batch_id}")
                success = orchestrator.run_full_pipeline_for_batch(args.batch_id)
            
            if success:
                self.logger.info(f"Successfully processed batch {args.batch_id}")
                return 0
            else:
                # Since this is a logical failure (not an exception), log it and return error code
                self.logger.error(f"Failed to process batch {args.batch_id}")
                return 1
                
        except Exception as e:
            # Log exception with context
            from ecod.error_handlers import log_exception
            from ecod.exceptions import PipelineError
            
            context = {
                "batch_id": args.batch_id,
                "blast_only": args.blast_only,
                "hhsearch_only": args.hhsearch_only,
                "domain_only": args.domain_only
            }
            
            log_exception(self.logger, e, context=context)
            
            # Re-raise as PipelineError
            raise PipelineError(f"Error running batch {args.batch_id}: {str(e)}", context) from e

# Add these functions to maintain compatibility with main.py

def setup_parser(parser: argparse.ArgumentParser) -> None:
    """Set up the argument parser for batch commands"""
    # Create a temporary command object just for setup
    subparsers = parser.add_subparsers(dest='command', help='Batch command')
    
    # Create command
    create_parser = subparsers.add_parser('create', help=COMMANDS['create'])
    create_parser.add_argument('--limit', type=int, default=10,
                            help='Maximum number of proteins to include')
    create_parser.add_argument('--type', type=str, default='full',
                            help='Type of batch (blast, hhsearch, full)')
    create_parser.add_argument('--output-dir', type=str,
                            help='Output directory override')
    
    # Status command
    status_parser = subparsers.add_parser('status', help=COMMANDS['status'])
    status_parser.add_argument('--batch-id', type=int, required=True,
                            help='Batch ID to check')
    
    # List command
    list_parser = subparsers.add_parser('list', help=COMMANDS['list'])
    list_parser.add_argument('--limit', type=int, default=10,
                          help='Maximum batches to list')
    list_parser.add_argument('--status', type=str,
                          help='Filter by status (created, processing, completed)')
    
    # Run command
    run_parser = subparsers.add_parser('run', help=COMMANDS['run'])
    run_parser.add_argument('--batch-id', type=int, required=True,
                         help='Batch ID to run')
    run_parser.add_argument('--blast-only', action='store_true',
                         help='Run only BLAST pipeline')
    run_parser.add_argument('--hhsearch-only', action='store_true',
                         help='Run only HHSearch pipeline')
    run_parser.add_argument('--domain-only', action='store_true',
                         help='Run only domain analysis')
    run_parser.add_argument('--partitioned', action='store_true',
                 help='Use adaptive processing paths based on BLAST confidence')

def run_command(args: argparse.Namespace) -> int:
    """Run the specified batch command"""
    cmd = BatchCommand(args.config)
    return cmd.run_command(args)