"""
Domain analysis commands for the ECOD pipeline
"""

import argparse
import logging

from typing import Dict, Any
from ecod.config import ConfigManager
from ecod.pipelines.domain_analysis.summary import DomainSummary
from ecod.pipelines.domain_analysis.partition import DomainPartition
from ecod.db import DBManager

logger = logging.getLogger("ecod.cli.domain")

# Define commands in this group
COMMANDS = {
    'summary': 'Create domain summary from BLAST and HHSearch results',
    'partition': 'Determine domain boundaries and classifications',
    'analyze': 'Run complete domain analysis pipeline',
}

def setup_parser(parser: argparse.ArgumentParser) -> None:
    """Set up the argument parser for domain analysis commands"""
    subparsers = parser.add_subparsers(dest='command', help='Domain command')
    
    # Summary command
    summary_parser = subparsers.add_parser('summary', help=COMMANDS['summary'])
    summary_parser.add_argument('--batch-id', type=int,
                              help='Batch ID to process')
    summary_parser.add_argument('--pdb-id', type=str,
                              help='PDB identifier')
    summary_parser.add_argument('--chain-id', type=str,
                              help='Chain identifier')
    summary_parser.add_argument('--concise', action='store_true',
                              help='Generate concise summary (no HHSearch data)')
    summary_parser.add_argument('--reference', type=str,
                              help='Reference version')
    
    # Partition command
    partition_parser = subparsers.add_parser('partition', help=COMMANDS['partition'])
    partition_parser.add_argument('--batch-id', type=int,
                                help='Batch ID to process')
    partition_parser.add_argument('--pdb-id', type=str,
                                help='PDB identifier')
    partition_parser.add_argument('--chain-id', type=str,
                                help='Chain identifier')
    partition_parser.add_argument('--input-mode', type=str, default='struct_seqid',
                                help='Input mode (struct_seqid, seqid)')
    partition_parser.add_argument('--reference', type=str,
                                help='Reference version')
    partition_parser.add_argument('--concise', action='store_true',
                                help='Use concise summary files')
    partition_parser.add_argument('--assembly', action='store_true',
                                help='Process as assembly (multiple chains)')
    
    # Analyze command (full pipeline)
    analyze_parser = subparsers.add_parser('analyze', help=COMMANDS['analyze'])
    analyze_parser.add_argument('--batch-id', type=int, required=True,
                              help='Batch ID to process')
    analyze_parser.add_argument('--concise', action='store_true',
                              help='Use concise summary')
    analyze_parser.add_argument('--limit', type=int, default=10,
                              help='Maximum proteins to process')

def run_command(args: argparse.Namespace) -> int:
    """Run the specified domain analysis command"""
    # Load configuration
    config_manager = ConfigManager(args.config)
    db_config = config_manager.get_db_config()
    db = DBManager(db_config)
    
    # Handle different commands
    if args.command == 'summary':
        return _run_summary(args, db, config_manager)
    elif args.command == 'partition':
        return _run_partition(args, db, config_manager)
    elif args.command == 'analyze':
        return _run_analysis(args, db, config_manager)
    else:
        logger.error(f"Unknown command: {args.command}")
        return 1

def _run_summary(args: argparse.Namespace, db: DBManager, config_manager: ConfigManager) -> int:
    """Run domain summary creation"""
    domain_summary = DomainSummary(config_manager.config_path)
    
    if args.batch_id:
        # Get batch information
        query = """
        SELECT 
            b.id, b.batch_name, b.base_path, b.type, b.ref_version
        FROM 
            ecod_schema.batch b
        WHERE 
            b.id = %s
        """
        batch_rows = db.execute_dict_query(query, (args.batch_id,))
        
        if not batch_rows:
            logger.error(f"Batch {args.batch_id} not found")
            return 1
        
        batch_info = batch_rows[0]
        reference = args.reference or batch_info["ref_version"]
        
        # Get batch proteins
        query = """
        SELECT 
            ps.id, p.id AS protein_id, p.pdb_id, p.chain_id, p.source_id,
            ps.current_stage, ps.status, ps.relative_path
        FROM 
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        WHERE 
            ps.batch_id = %s
            AND ps.status = 'success'
            AND EXISTS (
                SELECT 1 FROM ecod_schema.process_file pf
                WHERE pf.process_id = ps.id
                AND pf.file_type IN ('chain_blast_result', 'domain_blast_result')
                AND pf.file_exists = TRUE
            )
        LIMIT %s
        """
        limit = args.limit if hasattr(args, 'limit') else 1
        protein_rows = db.execute_dict_query(query, (args.batch_id, limit))
        
        if not protein_rows:
            logger.error(f"No proteins ready for domain analysis in batch {args.batch_id}")
            return 1
        
        # Process each protein
        success = True
        for protein in protein_rows:
            pdb_id = protein["pdb_id"]
            chain_id = protein["chain_id"]
            
            summary_file = domain_summary.create_summary(
                pdb_id,
                chain_id,
                reference,
                batch_info["base_path"],
                args.concise
            )
            
            if not summary_file:
                logger.error(f"Failed to create domain summary for {pdb_id}_{chain_id}")
                success = False
                continue
                
            logger.info(f"Created domain summary for {pdb_id}_{chain_id}")
        
        return 0 if success else 1
    
    elif args.pdb_id and args.chain_id:
        # Process a single protein
        if not args.reference:
            logger.error("Reference version is required when processing a single protein")
            return 1
        
        dump_dir = config_manager.get_path('output_dir', './output')
        
        summary_file = domain_summary.create_summary(
            args.pdb_id,
            args.chain_id,
            args.reference,
            dump_dir,
            args.concise
        )
        
        if not summary_file:
            logger.error(f"Failed to create domain summary for {args.pdb_id}_{args.chain_id}")
            return 1
            
        logger.info(f"Created domain summary for {args.pdb_id}_{args.chain_id}")
        return 0
        
    else:
        logger.error("Either --batch-id or both --pdb-id and --chain-id must be specified")
        return 1

def _run_partition(args: argparse.Namespace, db: DBManager, config_manager: ConfigManager) -> int:
    """Run domain partition"""
    domain_partition = DomainPartition(config_manager.config_path)
    
    if args.batch_id:
        # Get batch information
        query = """
        SELECT 
            b.id, b.batch_name, b.base_path, b.type, b.ref_version
        FROM 
            ecod_schema.batch b
        WHERE 
            b.id = %s
        """
        batch_rows = db.execute_dict_query(query, (args.batch_id,))
        
        if not batch_rows:
            logger.error(f"Batch {args.batch_id} not found")
            return 1
        
        batch_info = batch_rows[0]
        reference = args.reference or batch_info["ref_version"]
        
        # Get batch proteins
        query = """
        SELECT 
            ps.id, p.id AS protein_id, p.pdb_id, p.chain_id, p.source_id,
            ps.current_stage, ps.status, ps.relative_path
        FROM 
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        WHERE 
            ps.batch_id = %s
            AND ps.status = 'success'
            AND EXISTS (
                SELECT 1 FROM ecod_schema.process_file pf
                WHERE pf.process_id = ps.id
                AND pf.file_type = 'domain_summary'
                AND pf.file_exists = TRUE
            )
        LIMIT %s
        """
        limit = args.limit if hasattr(args, 'limit') else 1
        protein_rows = db.execute_dict_query(query, (args.batch_id, limit))
        
        if not protein_rows:
            logger.error(f"No proteins ready for domain partitioning in batch {args.batch_id}")
            return 1
        
        # Process each protein
        success = True
        for protein in protein_rows:
            pdb_id = protein["pdb_id"]
            chain_id = protein["chain_id"]
            
            if args.assembly:
                # For assembly, we need to gather all chains from the same PDB
                chains = [c["chain_id"] for c in protein_rows if c["pdb_id"] == pdb_id]
                
                partition_file = domain_partition.partition_domains_assembly(
                    pdb_id,
                    chains,
                    batch_info["base_path"],
                    args.input_mode,
                    reference,
                    args.concise
                )
            else:
                partition_file = domain_partition.partition_domains(
                    pdb_id,
                    chain_id,
                    batch_info["base_path"],
                    args.input_mode,
                    reference,
                    args.concise
                )
            
            if not partition_file:
                logger.error(f"Failed to partition domains for {pdb_id}_{chain_id}")
                success = False
                continue
                
            logger.info(f"Created domain partition for {pdb_id}_{chain_id}")
        
        return 0 if success else 1
    
    elif args.pdb_id and args.chain_id:
        # Process a single protein
        if not args.reference:
            logger.error("Reference version is required when processing a single protein")
            return 1
        
        dump_dir = config_manager.get_path('output_dir', './output')
        
        if args.assembly:
            # Treat chain_id as comma-separated list
            chains = args.chain_id.split(',')
            
            partition_file = domain_partition.partition_domains_assembly(
                args.pdb_id,
                chains,
                dump_dir,
                args.input_mode,
                args.reference,
                args.concise
            )
        else:
            partition_file = domain_partition.partition_domains(
                args.pdb_id,
                args.chain_id,
                dump_dir,
                args.input_mode,
                args.reference,
                args.concise
            )
        
        if not partition_file:
            logger.error(f"Failed to partition domains for {args.pdb_id}_{args.chain_id}")
            return 1
            
        logger.info(f"Created domain partition for {args.pdb_id}_{args.chain_id}")
        return 0
        
    else:
        logger.error("Either --batch-id or both --pdb-id and --chain-id must be specified")
        return 1

def _run_analysis(args: argparse.Namespace, db: DBManager, config_manager: ConfigManager) -> int:
    """Run the complete domain analysis pipeline"""
    from ecod.pipelines.domain_analysis.pipeline import DomainAnalysisPipeline
    
    # Initialize domain analysis pipeline
    domain_pipeline = DomainAnalysisPipeline(config_manager.config_path)
    
    # Run pipeline
    success = domain_pipeline.run_pipeline(
        args.batch_id,
        args.concise,
        args.limit
    )
    
    return 0 if success else 1