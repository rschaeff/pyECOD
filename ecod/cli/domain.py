# ecod/cli/domain.py
"""
Domain analysis commands for the ECOD pipeline
"""

import argparse
import logging
import os
import json
from typing import Dict, Any, List, Optional

from ecod.cli.base_command import BaseCommand, handle_command_errors
from ecod.pipelines.domain_analysis.summary import DomainSummary
from ecod.pipelines.domain_analysis.partition import DomainPartition
from ecod.utils.path_utils import (
    get_standardized_paths,
    get_file_db_path,
    resolve_file_path
)

logger = logging.getLogger("ecod.cli.domain")

# Expanded commands in this group
COMMANDS = {
    'summary': 'Create domain summary from BLAST and HHSearch results',
    'partition': 'Determine domain boundaries and classifications',
    'analyze': 'Run complete domain analysis pipeline',
    'validate': 'Validate domain summary inputs',
    'repair': 'Fix missing or problematic domain files',
    'stats': 'Generate statistics about domains in a batch',
}

class DomainCommand(BaseCommand):
    """Command handler for domain analysis operations"""

    def setup_parser(self, parser: argparse.ArgumentParser) -> None:
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
        summary_parser.add_argument('--blast-only', action='store_true',
                                  help='Skip HHSearch data and use only BLAST results')
        summary_parser.add_argument('--reference', type=str,
                                  help='Reference version')
        summary_parser.add_argument('--limit', type=int,
                                  help='Maximum proteins to process')
        summary_parser.add_argument('--reps-only', action='store_true',
                                  help='Process only representative proteins')
        summary_parser.add_argument('--force', action='store_true',
                                  help='Force regeneration of domain files')

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
        partition_parser.add_argument('--blast-only', action='store_true',
                                    help='Use summary files without HHSearch data')
        partition_parser.add_argument('--assembly', action='store_true',
                                    help='Process as assembly (multiple chains)')
        partition_parser.add_argument('--limit', type=int,
                                    help='Maximum proteins to process')
        partition_parser.add_argument('--reps-only', action='store_true',
                                    help='Process only representative proteins')
        partition_parser.add_argument('--force', action='store_true',
                                    help='Force regeneration of domain files')
        partition_parser.add_argument('--process-ids', type=int, nargs='+',
                                    help='Specific process IDs to process')

        # Analyze command (full pipeline)
        analyze_parser = subparsers.add_parser('analyze', help=COMMANDS['analyze'])
        analyze_parser.add_argument('--batch-id', type=int, required=True,
                                  help='Batch ID to process')
        analyze_parser.add_argument('--blast-only', action='store_true',
                                  help='Skip HHSearch and use only BLAST results')
        analyze_parser.add_argument('--limit', type=int, default=10,
                                  help='Maximum proteins to process')
        analyze_parser.add_argument('--reps-only', action='store_true',
                                  help='Process only representative proteins')
        analyze_parser.add_argument('--force', action='store_true',
                                  help='Force regeneration of domain files')
        analyze_parser.add_argument('--partition-only', action='store_true',
                                  help='Run only the partition step on batches with complete summaries')
        analyze_parser.add_argument('--process-ids', type=str,
                                  help='Comma-separated list of process IDs to analyze')
        analyze_parser.add_argument('--validate-only', action='store_true',
                                  help='Only validate database and filesystem structures')
        analyze_parser.add_argument('--adaptive', action='store_true',
                                  help='Use adaptive domain partitioning strategy')

        # Validate command
        validate_parser = subparsers.add_parser('validate', help=COMMANDS['validate'])
        validate_parser.add_argument('--batch-id', type=int, required=True,
                                   help='Batch ID to validate')
        validate_parser.add_argument('--limit', type=int,
                                   help='Limit number of proteins to validate')
        validate_parser.add_argument('--output', type=str,
                                   help='Output JSON file for statistics')
        validate_parser.add_argument('--reps-only', action='store_true',
                                   help='Validate only representative proteins')

        # Repair command
        repair_parser = subparsers.add_parser('repair', help=COMMANDS['repair'])
        repair_subparsers = repair_parser.add_subparsers(dest='repair_action', help='Repair action')

        # Missing partitions subcommand
        missing_parser = repair_subparsers.add_parser('missing', help='Find and fix missing domain partitions')
        missing_parser.add_argument('--batch-id', type=int, required=True,
                                  help='Batch ID to process')
        missing_parser.add_argument('--summary-only', action='store_true',
                                  help='Only show summary, don\'t repair')
        missing_parser.add_argument('--blast-only', action='store_true',
                                  help='Use only BLAST results (no HHSearch)')
        missing_parser.add_argument('--limit', type=int,
                                  help='Limit number of proteins to process')
        missing_parser.add_argument('--force', action='store_true',
                                  help='Force regeneration of domain files')
        missing_parser.add_argument('--reps-only', action='store_true',
                                  help='Process only representative proteins')

        # DB paths subcommand
        db_paths_parser = repair_subparsers.add_parser('db_paths', help='Repair database file paths')
        db_paths_parser.add_argument('--batch-id', type=int, required=True,
                                   help='Batch ID to process')
        db_paths_parser.add_argument('--dry-run', action='store_true',
                                   help='Don\'t actually update anything')

        # Stats command
        stats_parser = subparsers.add_parser('stats', help=COMMANDS['stats'])
        stats_parser.add_argument('--batch-id', type=int, required=True,
                               help='Batch ID to analyze')
        stats_parser.add_argument('--output', type=str,
                               help='Output JSON file for statistics')
        stats_parser.add_argument('--errors', action='store_true',
                               help='Show error messages for failed analyses')

    @handle_command_errors
    def run_command(self, args: argparse.Namespace) -> int:
        """Run the specified domain analysis command"""
        # Set force_overwrite flag if needed
        if hasattr(args, 'force') and args.force:
            self.context.set_force_overwrite(True)
            self.logger.info("Force overwrite enabled")

        # Handle different commands
        if args.command == 'summary':
            return self._run_summary(args)
        elif args.command == 'partition':
            return self._run_partition(args)
        elif args.command == 'analyze':
            return self._run_analysis(args)
        elif args.command == 'validate':
            return self._run_validate(args)
        elif args.command == 'repair':
            return self._run_repair(args)
        elif args.command == 'stats':
            return self._run_stats(args)
        else:
            self.logger.error(f"Unknown command: {args.command}")
            return 1

    def _run_summary(self, args: argparse.Namespace) -> int:
        """Run domain summary creation"""
        # Implementation using DomainSummary from ecod.pipelines.domain_analysis.summary
        domain_summary = DomainSummary(self.context)

        # Process batch or single protein based on args
        # [Implementation details...]

        return 0

    def _run_partition(self, args: argparse.Namespace) -> int:
        """Run domain partition"""
        # Implementation using DomainPartition from ecod.pipelines.domain_analysis.partition
        domain_partition = DomainPartition(self.context)

        if args.batch_id:
            # Get batch info
            batch_info = self._get_batch_info(args.batch_id)
            if not batch_info:
                return 1

            batch_path = batch_info['base_path']
            reference = args.reference or batch_info['ref_version']

            # If process_ids are provided, use specific processing
            if args.process_ids:
                self.logger.info(f"Processing specific process IDs for batch {args.batch_id}")
                success = domain_partition.process_specific_ids(
                    args.batch_id,
                    args.process_ids,
                    batch_path,
                    reference,
                    args.blast_only
                )
                return 0 if success else 1
            else:
                # Process entire batch
                domain_files = domain_partition.process_batch(
                    args.batch_id,
                    batch_path,
                    reference,
                    args.blast_only,
                    args.limit,
                    args.reps_only
                )

                return 0 if domain_files else 1

        elif args.pdb_id and args.chain_id:
            # Process single protein
            # [Implementation details...]
            return 0

        else:
            self.logger.error("Either --batch-id or both --pdb-id and --chain-id must be specified")
            return 1

    def _run_analysis(self, args: argparse.Namespace) -> int:
        """Run the complete domain analysis pipeline"""
        from ecod.pipelines.domain_analysis.pipeline import DomainAnalysisPipeline

        # Initialize domain analysis pipeline
        domain_pipeline = DomainAnalysisPipeline(self.context)

        # Handle different analysis options
        if args.partition_only:
            self.logger.info(f"Running partition only for batch {args.batch_id}")

            # Check if summaries are complete
            summary_status = self._verify_summary_completion(args.batch_id)
            if not summary_status['complete']:
                self.logger.error(f"Cannot run partition: Summaries incomplete ({summary_status['complete_count']}/{summary_status['total_count']} complete)")
                return 1

            # Get batch info
            batch_info = self._get_batch_info(args.batch_id)
            if not batch_info:
                return 1

            # Initialize partition component
            from ecod.pipelines.domain_analysis.partition import DomainPartition
            partition = DomainPartition(self.context)

            # Run partition with process IDs if specified
            if args.process_ids:
                process_ids = [int(pid.strip()) for pid in args.process_ids.split(',')]
                result = partition.process_specific_ids(
                    args.batch_id,
                    process_ids,
                    batch_info['base_path'],
                    batch_info['ref_version'],
                    args.blast_only
                )
            else:
                result = partition.process_batch(
                    args.batch_id,
                    batch_info['base_path'],
                    batch_info['ref_version'],
                    args.blast_only,
                    args.limit,
                    args.reps_only
                )

            return 0 if result else 1

        elif args.validate_only:
            self.logger.info(f"Validating partition inputs for batch {args.batch_id}")

            # Run validation logic
            # [Implementation details...]

            return 0

        elif args.adaptive:
            self.logger.info(f"Running adaptive domain partitioning for batch {args.batch_id}")

            # Import and use the adaptive partitioning strategy
            # [Implementation details...]

            return 0

        else:
            # Run standard pipeline
            success = domain_pipeline.run_pipeline(
                args.batch_id,
                args.blast_only,
                args.limit,
                args.reps_only
            )

            return 0 if success else 1

    def _run_validate(self, args: argparse.Namespace) -> int:
        """Validate domain summary files"""
        # Implementation of validation logic from domain_partition_tools.py
        # [Implementation details...]
        return 0

    def _run_repair(self, args: argparse.Namespace) -> int:
        """Run repair operations"""
        if args.repair_action == 'missing':
            return self._repair_missing_partitions(args)
        elif args.repair_action == 'db_paths':
            return self._repair_db_paths(args)
        else:
            self.logger.error(f"Unknown repair action: {args.repair_action}")
            return 1

    def _repair_missing_partitions(self, args: argparse.Namespace) -> int:
        """Find and fix missing domain partitions"""
        # Implementation from domain_partition_tools.py repair mode
        # [Implementation details...]
        return 0

    def _repair_db_paths(self, args: argparse.Namespace) -> int:
        """Repair database file paths"""
        # Implementation from domain_partition_tools.py repair mode
        # [Implementation details...]
        return 0

    def _run_stats(self, args: argparse.Namespace) -> int:
        """Generate domain statistics"""
        # Implementation from domain_partition_tools.py analyze mode
        # [Implementation details...]
        return 0

    def _get_batch_info(self, batch_id: int) -> Optional[Dict[str, Any]]:
        """Get batch information from database"""
        query = """
        SELECT id, batch_name, base_path, ref_version
        FROM ecod_schema.batch
        WHERE id = %s
        """

        batch_result = self.db.execute_dict_query(query, (batch_id,))

        if not batch_result:
            self.logger.error(f"Batch {batch_id} not found")
            return None

        return batch_result[0]

    def _verify_summary_completion(self, batch_id: int) -> Dict[str, Any]:
        """Verify that domain summaries are complete for all proteins in batch"""
        query = """
        SELECT
            COUNT(*) as total,
            SUM(CASE WHEN EXISTS (
                SELECT 1 FROM ecod_schema.process_file pf
                WHERE pf.process_id = ps.id
                AND pf.file_type = 'domain_summary'
                AND pf.file_exists = TRUE
            ) THEN 1 ELSE 0 END) as complete_count
        FROM
            ecod_schema.process_status ps
        WHERE
            ps.batch_id = %s
        """

        results = self.db.execute_dict_query(query, (batch_id,))[0]
        total = results.get('total', 0)
        complete = results.get('complete_count', 0)
        is_complete = total > 0 and total == complete

        self.logger.info(f"Summary completion: {complete}/{total} ({is_complete})")

        return {
            'total_count': total,
            'complete_count': complete,
            'complete': is_complete
        }

# Add these functions to maintain compatibility with main.py

def setup_parser(parser: argparse.ArgumentParser) -> None:
    """Set up the argument parser for domain analysis commands"""
    # Create temporary command to setup parser
    cmd = DomainCommand()
    cmd.setup_parser(parser)

def run_command(args: argparse.Namespace) -> int:
    """Run the specified domain analysis command"""
    cmd = DomainCommand(args.config)
    return cmd.run_command(args)
