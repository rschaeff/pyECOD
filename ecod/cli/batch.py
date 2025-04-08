"""
Batch processing commands for the ECOD pipeline
"""

import argparse
import logging
import os
from typing import Dict, Any
from datetime import datetime

from ecod.config import ConfigManager
from ecod.core.db_manager import DBManager
from ecod.core.models import Batch

logger = logging.getLogger("ecod.cli.batch")

# Define commands in this group
COMMANDS = {
    'create': 'Create a new processing batch',
    'status': 'Check batch status',
    'list': 'List all batches',
    'run': 'Run a full pipeline on a batch'
}

def setup_parser(parser: argparse.ArgumentParser) -> None:
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

def run_command(args: argparse.Namespace) -> int:
    """Run the specified batch command"""
    # Load configuration
    config_manager = ConfigManager(args.config)
    db_config = config_manager.get_db_config()
    db = DBManager(db_config)
    
    # Handle different commands
    if args.command == 'create':
        return _create_batch(args, db, config_manager)
    elif args.command == 'status':
        return _check_status(args, db)
    elif args.command == 'list':
        return _list_batches(args, db)
    elif args.command == 'run':
        return _run_batch(args, config_manager)
    else:
        logger.error(f"Unknown command: {args.command}")
        return 1

def _create_batch(args: argparse.Namespace, db: DBManager, config_manager: ConfigManager) -> int:
    """Create a new processing batch"""
    # Override output directory if specified
    config = config_manager.config
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
    
    logger.info(f"Created batch {batch_name} with ID {batch_id}")
    print(f"Created batch {batch_id}: {batch_name}")
    print(f"  - Type: {args.type}")
    print(f"  - Items: {len(proteins)}")
    print(f"  - Path: {batch_path}")
    
    return 0

def _check_status(args: argparse.Namespace, db: DBManager) -> int:
    """Check batch status"""
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

def _list_batches(args: argparse.Namespace, db: DBManager) -> int:
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

def _run_batch(args: argparse.Namespace, config_manager: ConfigManager) -> int:
    """Run a batch through the pipeline"""
    # Import orchestrator
    from ecod.pipelines.orchestrator import PipelineOrchestrator
    
    # Initialize orchestrator
    orchestrator = PipelineOrchestrator(config_manager.config_path)
    
    # Determine what to run
    if args.blast_only:
        logger.info(f"Running BLAST pipeline for batch {args.batch_id}")
        success = orchestrator.run_blast_pipeline(args.batch_id)
    elif args.hhsearch_only:
        logger.info(f"Running HHSearch pipeline for batch {args.batch_id}")
        success = orchestrator.run_hhsearch_pipeline(args.batch_id)
    elif args.domain_only:
        logger.info(f"Running domain analysis for batch {args.batch_id}")
        success = orchestrator.run_domain_analysis(args.batch_id)
    else:
        # Run full pipeline
        logger.info(f"Running full pipeline for batch {args.batch_id}")
        success = orchestrator.run_full_pipeline_for_batch(args.batch_id)
    
    if success:
        logger.info(f"Successfully processed batch {args.batch_id}")
        return 0
    else:
        logger.error(f"Failed to process batch {args.batch_id}")
        return 1