"""
BLAST-related commands for the ECOD pipeline
"""

import argparse
import logging
from typing import Dict, Any

from ecod.cli.base_command import BaseCommand
from ecod.pipelines.blast_pipeline import BlastPipeline
from ecod.core.context import ApplicationContext
from ecod.config import ConfigManager
from ecod.db import DBManager


logger = logging.getLogger("ecod.cli.blast")

# Define commands in this group
COMMANDS = {
    'run': 'Run BLAST analyses on protein sequences',
    'check': 'Check the status of BLAST jobs',
    'parse': 'Parse BLAST results'
}

class BlastCommand(BaseCommand):
    """Command handler for BLAST operations"""
    
    def setup_parser(self, parser: argparse.ArgumentParser) -> None:
        """Set up the argument parser for BLAST commands"""
        subparsers = parser.add_subparsers(dest='command', help='BLAST command')
        
        # Run command
        run_parser = subparsers.add_parser('run', help=COMMANDS['run'])
        run_parser.add_argument('--batch-id', type=int,
                             help='Use existing batch')
        run_parser.add_argument('--create-batch', action='store_true',
                             help='Create new batch')
        run_parser.add_argument('--max-chains', type=int, default=None,
                             help='Maximum chains for new batch')
        run_parser.add_argument('--batch-size', type=int, default=100,
                             help='Batch size for BLAST jobs')
        run_parser.add_argument('--chain', action='store_true',
                             help='Run chain-wise BLAST')
        run_parser.add_argument('--domain', action='store_true',
                             help='Run domain-wise BLAST')
        
        # Check command
        check_parser = subparsers.add_parser('check', help=COMMANDS['check'])
        check_parser.add_argument('--batch-id', type=int,
                               help='Check specific batch')
        
        # Parse command
        parse_parser = subparsers.add_parser('parse', help=COMMANDS['parse'])
        parse_parser.add_argument('--process-id', type=int, required=True,
                               help='Process ID to parse results for')
    
    def run_command(self, args: argparse.Namespace) -> int:
        """Run the specified BLAST command"""
        # Initialize BLAST pipeline
        from ecod.core.job_manager import JobManager
        job_manager = JobManager(self.config)
        blast_pipeline = BlastPipeline(self.db, job_manager, self.config)
        
        # Handle different commands
        if args.command == 'run':
            return self._run_blast(args, blast_pipeline)
        elif args.command == 'check':
            return self._check_blast(args, blast_pipeline)
        elif args.command == 'parse':
            return self._parse_blast(args, blast_pipeline)
        else:
            self.logger.error(f"Unknown command: {args.command}")
            return 1
    
    def _run_blast(self, args: argparse.Namespace, blast_pipeline: Any) -> int:
        """Run BLAST analysis"""
        if args.create_batch:
            # Get chains to process
            chains = blast_pipeline.get_unclassified_chains(args.max_chains)
            if not chains:
                self.logger.error("No unclassified chains found")
                return 1
            
            # Create batch
            batch = blast_pipeline.create_batch(chains, args.batch_size)
            batch_id = batch.id
            self.logger.info(f"Created batch {batch_id} with {len(chains)} chains")
        elif args.batch_id:
            batch_id = args.batch_id
        else:
            self.logger.error("Either --batch-id or --create-batch must be specified")
            return 1
        
        # Run BLAST
        success = True
        if args.chain or not (args.chain or args.domain):  # Default to chain if neither specified
            self.logger.info(f"Running chain-wise BLAST for batch {batch_id}")
            chain_job_ids = blast_pipeline.run_chain_blast(batch_id, args.batch_size)
            self.logger.info(f"Submitted {len(chain_job_ids)} chain BLAST jobs")
            success = len(chain_job_ids) > 0
        
        if args.domain or not (args.chain or args.domain):  # Default to domain if neither specified
            self.logger.info(f"Running domain-wise BLAST for batch {batch_id}")
            domain_job_ids = blast_pipeline.run_domain_blast(batch_id, args.batch_size)
            self.logger.info(f"Submitted {len(domain_job_ids)} domain BLAST jobs")
            success = success and len(domain_job_ids) > 0
        
        return 0 if success else 1
    
    def _check_blast(self, args: argparse.Namespace, blast_pipeline: Any) -> int:
        """Check BLAST job status"""
        self.logger.info("Checking BLAST job status")
        blast_pipeline.check_job_status(args.batch_id)
        return 0
    
    def _parse_blast(self, args: argparse.Namespace, blast_pipeline: Any) -> int:
        """Parse BLAST results"""
        self.logger.info(f"Parsing BLAST results for process {args.process_id}")
        success = blast_pipeline.parse_results(args.process_id)
        return 0 if success else 1


def setup_parser(parser: argparse.ArgumentParser) -> None:
    """Set up the argument parser for BLAST commands"""
    subparsers = parser.add_subparsers(dest='command', help='BLAST command')
    
    # Run command
    run_parser = subparsers.add_parser('run', help=COMMANDS['run'])
    run_parser.add_argument('--batch-id', type=int,
                         help='Use existing batch')
    run_parser.add_argument('--create-batch', action='store_true',
                         help='Create new batch')
    run_parser.add_argument('--max-chains', type=int, default=None,
                         help='Maximum chains for new batch')
    run_parser.add_argument('--batch-size', type=int, default=100,
                         help='Batch size for BLAST jobs')
    run_parser.add_argument('--chain', action='store_true',
                         help='Run chain-wise BLAST')
    run_parser.add_argument('--domain', action='store_true',
                         help='Run domain-wise BLAST')
    
    # Check command
    check_parser = subparsers.add_parser('check', help=COMMANDS['check'])
    check_parser.add_argument('--batch-id', type=int,
                           help='Check specific batch')
    
    # Parse command
    parse_parser = subparsers.add_parser('parse', help=COMMANDS['parse'])
    parse_parser.add_argument('--process-id', type=int, required=True,
                           help='Process ID to parse results for')

def run_command(args: argparse.Namespace) -> int:
    """Run the specified BLAST command"""
    cmd = BlastCommand(args.config)
    return cmd.run_command(args)
