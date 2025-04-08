"""
HHSearch-related commands for the ECOD pipeline
"""

import argparse
import logging
from typing import Dict, Any

from ecod.config import ConfigManager
from ecod.pipelines.hhsearch_pipeline import HHSearchPipeline
from ecod.db import DBManager

logger = logging.getLogger("ecod.cli.hhsearch")

# Define commands in this group
COMMANDS = {
    'run': 'Run HHSearch analyses on protein sequences',
    'check': 'Check the status of HHSearch jobs',
    'parse': 'Parse HHSearch results',
    'profiles': 'Generate HHblits profiles'
}

def setup_parser(parser: argparse.ArgumentParser) -> None:
    """Set up the argument parser for HHSearch commands"""
    subparsers = parser.add_subparsers(dest='command', help='HHSearch command')
    
    # Run command
    run_parser = subparsers.add_parser('run', help=COMMANDS['run'])
    run_parser.add_argument('--batch-id', type=int,
                         help='Use existing batch')
    run_parser.add_argument('--create-batch', action='store_true',
                         help='Create new batch')
    run_parser.add_argument('--max-chains', type=int, default=None,
                         help='Maximum chains for new batch')
    run_parser.add_argument('--threads', type=int, default=8,
                         help='Number of threads for HHSearch/HHblits')
    run_parser.add_argument('--memory', type=str, default='16G',
                         help='Memory allocation for jobs')
    
    # Check command
    check_parser = subparsers.add_parser('check', help=COMMANDS['check'])
    check_parser.add_argument('--batch-id', type=int,
                           help='Check specific batch')
    
    # Parse command
    parse_parser = subparsers.add_parser('parse', help=COMMANDS['parse'])
    parse_parser.add_argument('--process-id', type=int, required=True,
                           help='Process ID to parse results for')
    
    # Profiles command
    profiles_parser = subparsers.add_parser('profiles', help=COMMANDS['profiles'])
    profiles_parser.add_argument('--batch-id', type=int, required=True,
                              help='Batch ID to generate profiles for')
    profiles_parser.add_argument('--threads', type=int, default=8,
                              help='Number of threads for HHblits')
    profiles_parser.add_argument('--memory', type=str, default='16G',
                              help='Memory allocation for jobs')

def run_command(args: argparse.Namespace) -> int:
    """Run the specified HHSearch command"""
    # Load configuration
    config_manager = ConfigManager(args.config)
    db_config = config_manager.get_db_config()
    db = DBManager(db_config)
    
    # Initialize HHSearch pipeline
    from ecod.core.job_manager import JobManager
    job_manager = JobManager(config_manager.config)
    hhsearch_pipeline = HHSearchPipeline(db, job_manager, config_manager.config)
    
    # Handle different commands
    if args.command == 'run':
        return _run_hhsearch(args, hhsearch_pipeline)
    elif args.command == 'check':
        return _check_hhsearch(args, hhsearch_pipeline)
    elif args.command == 'parse':
        return _parse_hhsearch(args, hhsearch_pipeline)
    elif args.command == 'profiles':
        return _generate_profiles(args, hhsearch_pipeline)
    else:
        logger.error(f"Unknown command: {args.command}")
        return 1

def _run_hhsearch(args: argparse.Namespace, hhsearch_pipeline: Any) -> int:
    """Run HHSearch analysis"""
    if args.create_batch:
        # Get chains to process
        chains = hhsearch_pipeline.get_unprocessed_chains(args.max_chains)
        if not chains:
            logger.error("No chains found for HHSearch analysis")
            return 1
        
        # Create batch
        batch_id = hhsearch_pipeline.create_batch(chains)
        logger.info(f"Created batch {batch_id} with {len(chains)} chains")
    elif args.batch_id:
        batch_id = args.batch_id
    else:
        logger.error("Either --batch-id or --create-batch must be specified")
        return 1
    
    # Run HHSearch jobs
    logger.info(f"Running HHSearch for batch {batch_id}")
    hhsearch_job_ids = hhsearch_pipeline.run_hhsearch(
        batch_id, 
        threads=args.threads, 
        memory=args.memory
    )
    
    logger.info(f"Submitted {len(hhsearch_job_ids)} HHSearch jobs")
    return 0 if hhsearch_job_ids else 1

def _check_hhsearch(args: argparse.Namespace, hhsearch_pipeline: Any) -> int:
    """Check HHSearch job status"""
    logger.info("Checking HHSearch job status")
    hhsearch_pipeline.check_status(args.batch_id)
    return 0

def _parse_hhsearch(args: argparse.Namespace, hhsearch_pipeline: Any) -> int:
    """Parse HHSearch results"""
    logger.info(f"Parsing HHSearch results for process {args.process_id}")
    # This depends on implementation details of hhsearch_pipeline
    # Assuming there's a method like parse_results similar to BLAST pipeline
    if hasattr(hhsearch_pipeline, 'parse_results'):
        success = hhsearch_pipeline.parse_results(args.process_id)
        return 0 if success else 1
    else:
        logger.error("HHSearch pipeline does not support parsing results directly")
        return 1

def _generate_profiles(args: argparse.Namespace, hhsearch_pipeline: Any) -> int:
    """Generate HHblits profiles"""
    logger.info(f"Generating HHblits profiles for batch {args.batch_id}")
    profile_job_ids = hhsearch_pipeline.generate_profiles(
        args.batch_id,
        threads=args.threads,
        memory=args.memory
    )
    
    logger.info(f"Submitted {len(profile_job_ids)} profile generation jobs")
    return 0 if profile_job_ids else 1