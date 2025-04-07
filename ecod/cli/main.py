# ecod/cli/main.py
import argparse
import logging
import sys
from pathlib import Path
import os

from ..pipelines.orchestrator import PipelineOrchestrator
from ..db.migration_manager import MigrationManager
from ..core.config import ConfigManager

def setup_logging(verbosity, log_file=None):
    """Configure logging based on verbosity level"""
    log_level = {
        0: logging.WARNING,
        1: logging.INFO,
        2: logging.DEBUG
    }.get(verbosity, logging.DEBUG)
    
    handlers = [
        logging.StreamHandler(sys.stdout)
    ]
    
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )

def main():
    # Create the top-level parser
    parser = argparse.ArgumentParser(description='ECOD Domain Partition Pipeline')
    
    # Global options
    parser.add_argument('--config', type=str, help='Path to configuration file')
    parser.add_argument('-v', '--verbose', action='count', default=0, 
                       help='Increase verbosity (can be used multiple times)')
    parser.add_argument('--log-file', type=str, 
                       help='Log to file in addition to stdout')
    
    # Create subparsers for commands
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # DB migrate command
    migrate_parser = subparsers.add_parser('migrate', help='Apply database migrations')
    migrate_parser.add_argument('--migrations-dir', type=str, 
                             default='ecod/db/migrations',
                             help='Directory containing migration files')
    
    # Full pipeline command
    run_parser = subparsers.add_parser('run', help='Run full pipeline')
    run_parser.add_argument('--batch-size', type=int, default=500, 
                         help='Maximum chains per batch')
    run_parser.add_argument('--max-chains', type=int, default=None,
                         help='Maximum chains to process')
    
    # BLAST pipeline command
    blast_parser = subparsers.add_parser('blast', help='Run BLAST pipeline')
    blast_parser.add_argument('--batch-id', type=int, help='Use existing batch')
    blast_parser.add_argument('--create-batch', action='store_true',
                           help='Create new batch')
    blast_parser.add_argument('--max-chains', type=int, default=None,
                           help='Maximum chains for new batch')
    blast_parser.add_argument('--batch-size', type=int, default=100,
                           help='Batch size for BLAST jobs')
    
    # HHsearch pipeline command
    hhsearch_parser = subparsers.add_parser('hhsearch', help='Run HHsearch pipeline')
    hhsearch_parser.add_argument('--batch-id', type=int, help='Use existing batch')
    hhsearch_parser.add_argument('--create-batch', action='store_true',
                              help='Create new batch')
    hhsearch_parser.add_argument('--max-chains', type=int, default=None,
                              help='Maximum chains for new batch')
    hhsearch_parser.add_argument('--threads', type=int, default=8,
                              help='Threads for HHblits/HHsearch')
    hhsearch_parser.add_argument('--memory', type=str, default='16G',
                              help='Memory for HHblits/HHsearch')
    
    # Status check command
    status_parser = subparsers.add_parser('status', help='Check pipeline status')
    status_parser.add_argument('--batch-id', type=int,
                            help='Check specific batch')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose, args.log_file)
    
    # Load configuration
    config_manager = ConfigManager(args.config)
    
    # Execute command
    if args.command == 'migrate':
        # Run database migrations
        db_config = config_manager.get_db_config()
        migrations_dir = args.migrations_dir
        
        migration_manager = MigrationManager(db_config, migrations_dir)
        migration_manager.apply_migrations()
        
    elif args.command == 'run':
        # Run full pipeline
        orchestrator = PipelineOrchestrator(args.config)
        orchestrator.run_full_pipeline(
            batch_size=args.batch_size,
            max_chains=args.max_chains
        )
        
# ecod/cli/main.py (continued)
    elif args.command == 'blast':
        # Run BLAST pipeline
        orchestrator = PipelineOrchestrator(args.config)
        
        if args.create_batch:
            chains = orchestrator.blast.get_unclassified_chains(args.max_chains)
            if chains:
                batch = orchestrator.blast.create_batch(chains, args.batch_size)
                batch_id = batch.id
                print(f"Created batch {batch_id} with {len(chains)} chains")
            else:
                print("No unclassified chains found")
                return
        elif args.batch_id:
            batch_id = args.batch_id
        else:
            print("Error: Either --batch-id or --create-batch must be specified")
            return
            
        # Run BLAST
        chain_job_ids = orchestrator.blast.run_chain_blast(batch_id, args.batch_size)
        domain_job_ids = orchestrator.blast.run_domain_blast(batch_id, args.batch_size)
        
        print(f"Submitted {len(chain_job_ids)} chain BLAST jobs")
        print(f"Submitted {len(domain_job_ids)} domain BLAST jobs")
        
    elif args.command == 'hhsearch':
        # Run HHsearch pipeline
        orchestrator = PipelineOrchestrator(args.config)
        
        if args.create_batch:
            # Get chains from database
            chains = orchestrator.hhsearch.get_unclassified_chains(args.max_chains)
            if chains:
                batch_id = orchestrator.hhsearch.create_batch(chains, args.batch_size)
                print(f"Created batch {batch_id} with {len(chains)} chains")
            else:
                print("No unclassified chains found")
                return
        elif args.batch_id:
            batch_id = args.batch_id
        else:
            print("Error: Either --batch-id or --create-batch must be specified")
            return
            
        # Generate profiles
        profile_job_ids = orchestrator.hhsearch.generate_profiles(
            batch_id, args.threads, args.memory
        )
        
        print(f"Submitted {len(profile_job_ids)} profile generation jobs")
        
    elif args.command == 'status':
        # Check pipeline status
        orchestrator = PipelineOrchestrator(args.config)
        orchestrator.check_status(args.batch_id)
        
    else:
        parser.print_help()

if __name__ == "__main__":
    main()