# ecod/cli/main.py
import argparse
import sys
import os
import logging
from pathlib import Path

from ..pipelines.orchestrator import PipelineOrchestrator
from ..db.migration_manager import MigrationManager
from ..core.config import ConfigManager
from ..core.logging_config import LoggingManager
from ..core.app_errors import handle_exceptions

@handle_exceptions
def main():
    # Create the top-level parser
    parser = argparse.ArgumentParser(description='ECOD Domain Partition Pipeline')
    
    # Global options
    parser.add_argument('--config', type=str, help='Path to configuration file')
    parser.add_argument('-v', '--verbose', action='count', default=0, 
                       help='Increase verbosity (can be used multiple times)')
    parser.add_argument('--log-file', type=str, 
                       help='Log to file in addition to stdout')
    parser.add_argument('--log-dir', type=str,
                       help='Directory for log files')
    parser.add_argument('--json', action='store_true',
                       help='Output results as JSON')
    
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
    
# ecod/cli/main.py (continued)
    # Setup logging with improved logging configuration
    log_level = min(50, 30 - (args.verbose * 10))  # 0=WARNING, 1=INFO, 2=DEBUG
    
    # Initialize logging with the improved LoggingManager
    logger = LoggingManager.configure(
        verbose=(log_level <= logging.DEBUG),
        log_file=args.log_file,
        log_dir=args.log_dir,
        component="ecod"
    )
    
    logger.info(f"ECOD Pipeline starting - log level: {logging.getLevelName(log_level)}")
    
    # Load configuration
    try:
        config_manager = ConfigManager(args.config)
        logger.info(f"Configuration loaded successfully from {args.config}")
    except Exception as e:
        logger.error(f"Failed to load configuration: {str(e)}")
        print(f"Error: Failed to load configuration: {str(e)}", file=sys.stderr)
        return 1
    
    # Execute command
    if args.command == 'migrate':
        # Run database migrations
        try:
            db_config = config_manager.get_db_config()
            migrations_dir = args.migrations_dir
            
            logger.info(f"Running database migrations from {migrations_dir}")
            migration_manager = MigrationManager(db_config, migrations_dir)
            result = migration_manager.apply_migrations()
            
            if result['success']:
                logger.info(f"Applied {result['migrations_applied']} migrations successfully")
                print(f"Applied {result['migrations_applied']} migrations successfully")
            else:
                logger.error(f"Migration failed: {result['error']}")
                print(f"Error: Migration failed: {result['error']}", file=sys.stderr)
                return 1
        except Exception as e:
            logger.error(f"Migration failed with unexpected error: {str(e)}", exc_info=True)
            print(f"Error: Migration failed: {str(e)}", file=sys.stderr)
            return 1
        
    elif args.command == 'run':
        # Run full pipeline
        try:
            logger.info("Initializing pipeline orchestrator")
            orchestrator = PipelineOrchestrator(args.config)
            
            logger.info(f"Running full pipeline - batch size: {args.batch_size}, max chains: {args.max_chains or 'unlimited'}")
            result = orchestrator.run_full_pipeline(
                batch_size=args.batch_size,
                max_chains=args.max_chains
            )
            
            if result and result.get('status') in ('completed', 'completed_with_errors'):
                logger.info(f"Pipeline completed with status: {result['status']}")
                
                # Print summary
                if args.json:
                    import json
                    print(json.dumps(result, indent=2))
                else:
                    print(f"Pipeline completed with status: {result['status']}")
                    print(f"Batch ID: {result.get('batch_id')}")
                    print(f"Steps completed: {', '.join(result.get('steps_completed', []))}")
                    if result.get('steps_failed'):
                        print(f"Steps with errors: {', '.join(result.get('steps_failed', []))}")
            else:
                if not result:
                    logger.error("Pipeline returned no result")
                    print("Error: Pipeline execution failed with no result", file=sys.stderr)
                else:
                    logger.error(f"Pipeline failed with status: {result.get('status')}")
                    error = result.get('error', 'Unknown error')
                    print(f"Error: Pipeline failed: {error}", file=sys.stderr)
                return 1
                
        except Exception as e:
            logger.error(f"Pipeline execution failed: {str(e)}", exc_info=True)
            print(f"Error: Pipeline execution failed: {str(e)}", file=sys.stderr)
            return 1
            
    elif args.command == 'blast':
        # Run BLAST pipeline
        try:
            orchestrator = PipelineOrchestrator(args.config)
            
            if args.create_batch:
                logger.info(f"Creating new BLAST batch - max chains: {args.max_chains or 'unlimited'}")
                chains = orchestrator.blast.get_unclassified_chains(args.max_chains)
                if chains:
                    batch = orchestrator.blast.create_batch(chains, args.batch_size)
                    batch_id = batch.id
                    logger.info(f"Created batch {batch_id} with {len(chains)} chains")
                    print(f"Created batch {batch_id} with {len(chains)} chains")
                else:
                    logger.warning("No unclassified chains found")
                    print("No unclassified chains found")
                    return 0
            elif args.batch_id:
                batch_id = args.batch_id
                logger.info(f"Using existing batch {batch_id}")
            else:
                logger.error("Either --batch-id or --create-batch must be specified")
                print("Error: Either --batch-id or --create-batch must be specified", file=sys.stderr)
                return 1
                
            # Run BLAST
            logger.info(f"Running chain BLAST for batch {batch_id}")
            chain_job_ids = orchestrator.blast.run_chain_blast(batch_id, args.batch_size)
            
            logger.info(f"Running domain BLAST for batch {batch_id}")
            domain_job_ids = orchestrator.blast.run_domain_blast(batch_id, args.batch_size)
            
            logger.info(f"Submitted {len(chain_job_ids)} chain BLAST jobs and {len(domain_job_ids)} domain BLAST jobs")
            print(f"Submitted {len(chain_job_ids)} chain BLAST jobs and {len(domain_job_ids)} domain BLAST jobs")
            
        except Exception as e:
            logger.error(f"BLAST pipeline failed: {str(e)}", exc_info=True)
            print(f"Error: BLAST pipeline failed: {str(e)}", file=sys.stderr)
            return 1
        
    elif args.command == 'hhsearch':
        # Run HHsearch pipeline
        try:
            orchestrator = PipelineOrchestrator(args.config)
            
            if args.create_batch:
                logger.info(f"Creating new HHsearch batch - max chains: {args.max_chains or 'unlimited'}")
                # Get chains from database
                chains = orchestrator.hhsearch.get_unclassified_chains(args.max_chains)
                if chains:
                    batch_id = orchestrator.hhsearch.create_batch(chains, args.batch_size)
                    logger.info(f"Created batch {batch_id} with {len(chains)} chains")
                    print(f"Created batch {batch_id} with {len(chains)} chains")
                else:
                    logger.warning("No unclassified chains found")
                    print("No unclassified chains found")
                    return 0
            elif args.batch_id:
                batch_id = args.batch_id
                logger.info(f"Using existing batch {batch_id}")
            else:
                logger.error("Either --batch-id or --create-batch must be specified")
                print("Error: Either --batch-id or --create-batch must be specified", file=sys.stderr)
                return 1
                
            # Generate profiles
            logger.info(f"Generating HHblits profiles for batch {batch_id}")
            profile_job_ids = orchestrator.hhsearch.generate_profiles(
                batch_id, args.threads, args.memory
            )
            
            logger.info(f"Submitted {len(profile_job_ids)} profile generation jobs")
            print(f"Submitted {len(profile_job_ids)} profile generation jobs")
            
        except Exception as e:
            logger.error(f"HHsearch pipeline failed: {str(e)}", exc_info=True)
            print(f"Error: HHsearch pipeline failed: {str(e)}", file=sys.stderr)
            return 1
        
    elif args.command == 'status':
        # Check pipeline status
        try:
            orchestrator = PipelineOrchestrator(args.config)
            
            logger.info(f"Checking pipeline status for batch {args.batch_id or 'all'}")
            status = orchestrator.check_status(args.batch_id)
            
            if args.json:
                import json
                print(json.dumps(status, indent=2))
            else:
                # Format and print status information
                if args.batch_id:
                    batch_info = status['batch']
                    print(f"Batch {batch_info['id']} - {batch_info['name']}")
                    print(f"Status: {batch_info['status']}")
                    print(f"Progress: {batch_info['completed_items']}/{batch_info['total_items']} items")
                    print(f"Type: {batch_info['type']}")
                    
                    if 'jobs' in status:
                        job_status = status['jobs']
                        print("\nJobs:")
                        print(f"  Completed: {job_status['completed']}")
                        print(f"  Running: {job_status['running']}")
                        print(f"  Failed: {job_status['failed']}")
                else:
                    # Summary of all batches
                    print("Batch Status Summary:")
                    for batch in status['batches']:
                        print(f"Batch {batch['id']} - {batch['name']}: {batch['status']} - "
                              f"{batch['completed_items']}/{batch['total_items']} completed")
            
        except Exception as e:
            logger.error(f"Status check failed: {str(e)}", exc_info=True)
            print(f"Error: Status check failed: {str(e)}", file=sys.stderr)
            return 1
        
    else:
        parser.print_help()
    
    return 0

if __name__ == "__main__":
    sys.exit(main())