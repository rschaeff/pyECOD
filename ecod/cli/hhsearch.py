# ecod/cli/hhsearch.py
import argparse
import logging
from typing import Dict, Any

from .base_command import BaseCommand
from ecod.pipelines.hhsearch_pipeline import HHSearchPipeline

# Define commands in this group
COMMANDS = {
    'generate': 'Generate HHblits profiles for sequences',
    'search': 'Run HHSearch against ECOD database',
    'parse': 'Parse HHSearch results',
    'check': 'Check the status of HHSearch jobs'
}

class HHSearchCommand(BaseCommand):
    """Command handler for HHSearch operations"""
    
    def setup_parser(self, parser: argparse.ArgumentParser) -> None:
        """Set up the argument parser for HHSearch commands"""
        subparsers = parser.add_subparsers(dest='command', help='HHSearch command')
        
        # Generate command
        generate_parser = subparsers.add_parser('generate', help=COMMANDS['generate'])
        generate_parser.add_argument('--batch-id', type=int, required=True,
                                  help='Batch ID to process')
        generate_parser.add_argument('--threads', type=int, default=8,
                                  help='Threads per job')
        generate_parser.add_argument('--memory', type=str, default='16G',
                                  help='Memory per job')
        
        # Search command
        search_parser = subparsers.add_parser('search', help=COMMANDS['search'])
        search_parser.add_argument('--batch-id', type=int, required=True,
                                help='Batch ID to process')
        search_parser.add_argument('--threads', type=int, default=8,
                                help='Threads per job')
        search_parser.add_argument('--memory', type=str, default='16G',
                                help='Memory per job')
        
        # Parse command
        parse_parser = subparsers.add_parser('parse', help=COMMANDS['parse'])
        parse_parser.add_argument('--batch-id', type=int, required=True,
                               help='Batch ID to process')
        
        # Check command
        check_parser = subparsers.add_parser('check', help=COMMANDS['check'])
        check_parser.add_argument('--batch-id', type=int,
                               help='Check specific batch')
    
    def run_command(self, args: argparse.Namespace) -> int:
        """Run the specified HHSearch command"""
        # Initialize HHSearch pipeline
        from ecod.core.job_manager import JobManager
        job_manager = JobManager(self.config)
        hhsearch_pipeline = HHSearchPipeline(self.db, job_manager, self.config)
        
        # Handle different commands
        if args.command == 'generate':
            return self._generate_profiles(args, hhsearch_pipeline)
        elif args.command == 'search':
            return self._run_search(args, hhsearch_pipeline)
        elif args.command == 'parse':
            return self._parse_results(args, hhsearch_pipeline)
        elif args.command == 'check':
            return self._check_status(args, hhsearch_pipeline)
        else:
            self.logger.error(f"Unknown command: {args.command}")
            return 1
    
    def _generate_profiles(self, args: argparse.Namespace, pipeline: Any) -> int:
        """Generate HHblits profiles"""
        self.logger.info(f"Generating profiles for batch {args.batch_id}")
        job_ids = pipeline.generate_profiles(
            args.batch_id, args.threads, args.memory
        )
        self.logger.info(f"Submitted {len(job_ids)} profile generation jobs")
        return 0 if job_ids else 1
    
    def _run_search(self, args: argparse.Namespace, pipeline: Any) -> int:
        """Run HHSearch"""
        self.logger.info(f"Running HHSearch for batch {args.batch_id}")
        job_ids = pipeline.run_hhsearch(
            args.batch_id, args.threads, args.memory
        )
        self.logger.info(f"Submitted {len(job_ids)} HHSearch jobs")
        return 0 if job_ids else 1
    
    def _parse_results(self, args: argparse.Namespace, pipeline: Any) -> int:
        """Parse HHSearch results"""
        self.logger.info(f"Parsing HHSearch results for batch {args.batch_id}")
        # This method would need to be implemented in the pipeline class
        if hasattr(pipeline, 'parse_results'):
            success = pipeline.parse_results(args.batch_id)
            return 0 if success else 1
        else:
            self.logger.error("parse_results method not implemented in pipeline")
            return 1
    
    def _check_status(self, args: argparse.Namespace, pipeline: Any) -> int:
        """Check HHSearch job status"""
        self.logger.info("Checking HHSearch job status")
        pipeline.check_status(args.batch_id)
        return 0