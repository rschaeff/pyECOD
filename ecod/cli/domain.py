# ecod/cli/domain.py
"""
Domain analysis commands for the ECOD pipeline

This module provides comprehensive domain analysis functionality including:
- Domain summary generation from BLAST and HHSearch results
- Domain partition determination and classification
- Batch processing with SLURM support
- Diagnostic and repair tools
- Real-time monitoring
"""

import argparse
import logging
import os
import sys
import json
import time
import glob
from typing import Dict, Any, List, Optional, Tuple
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

from ecod.cli.base_command import BaseCommand, handle_command_errors
from ecod.pipelines.domain_analysis.summary import DomainSummary
from ecod.pipelines.domain_analysis.partition import (
    DomainPartitionService,
    create_service,
    PartitionOptions,
    ValidationLevel,
    ProcessingMode,
    BatchPartitionResults,
    StatusTracker,
    EvidenceAnalyzer
)
from ecod.pipelines.domain_analysis.pipeline import DomainAnalysisPipeline
from ecod.models.pipeline.partition import DomainPartitionResult
from ecod.utils.path_utils import (
    get_standardized_paths,
    get_file_db_path,
    resolve_file_path,
    get_all_evidence_paths,
    find_files_with_legacy_paths
)
from ecod.jobs.factory import create_job_manager

logger = logging.getLogger("ecod.cli.domain")

# Expanded commands in this group
COMMANDS = {
    'summary': 'Create domain summary from BLAST and HHSearch results',
    'partition': 'Determine domain boundaries and classifications',
    'analyze': 'Analyze domain results and statistics',
    'diagnose': 'Diagnose domain processing issues',
    'repair': 'Fix missing or problematic domain files',
    'monitor': 'Monitor real-time domain processing status',
}

class DomainCommand(BaseCommand):
    """Command handler for domain analysis operations"""

    def __init__(self, config_path: Optional[str] = None):
        super().__init__(config_path)
        self._service = None
        self._pipeline = None

    @property
    def service(self) -> DomainPartitionService:
        """Lazy initialization of domain partition service"""
        if self._service is None:
            # Create service with custom configuration
            service_config = {
                'max_workers': 4,
                'use_multiprocessing': False,
                'save_intermediate': True,
                'track_status': True
            }
            self._service = DomainPartitionService(self.context, service_config)
        return self._service

    @property
    def pipeline(self) -> DomainAnalysisPipeline:
        """Lazy initialization of domain analysis pipeline"""
        if self._pipeline is None:
            self._pipeline = DomainAnalysisPipeline(self.context)
        return self._pipeline

    def setup_parser(self, parser: argparse.ArgumentParser) -> None:
        """Set up the argument parser for domain analysis commands"""
        subparsers = parser.add_subparsers(dest='command', help='Domain command')

        # Summary command
        self._setup_summary_parser(subparsers)

        # Partition command (main functionality)
        self._setup_partition_parser(subparsers)

        # Analyze command
        self._setup_analyze_parser(subparsers)

        # Diagnose command
        self._setup_diagnose_parser(subparsers)

        # Repair command
        self._setup_repair_parser(subparsers)

        # Monitor command
        self._setup_monitor_parser(subparsers)

    def _setup_summary_parser(self, subparsers):
        """Set up summary command parser"""
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

    def _setup_partition_parser(self, subparsers):
        """Set up partition command parser with all sub-actions"""
        partition_parser = subparsers.add_parser('partition', help=COMMANDS['partition'])
        partition_subparsers = partition_parser.add_subparsers(dest='action', help='Partition action')

        # Partition single protein
        single_parser = partition_subparsers.add_parser('single', help='Process a single protein')
        single_parser.add_argument('--pdb-id', type=str, required=True,
                                 help='PDB ID')
        single_parser.add_argument('--chain-id', type=str, required=True,
                                 help='Chain ID')
        single_parser.add_argument('--batch-id', type=int,
                                 help='Batch ID (optional)')
        single_parser.add_argument('--batch-path', type=str,
                                 help='Batch path (required if batch ID not provided)')
        single_parser.add_argument('--reference', type=str,
                                 help='Reference version (required if batch ID not provided)')
        single_parser.add_argument('--blast-only', action='store_true',
                                 help='Use only BLAST results')

        # Partition batch
        batch_parser = partition_subparsers.add_parser('batch', help='Process a batch')
        batch_parser.add_argument('--batch-id', type=int, required=True,
                                help='Batch ID to process')
        batch_parser.add_argument('--blast-only', action='store_true',
                                help='Use only BLAST results (no HHSearch)')
        batch_parser.add_argument('--limit', type=int,
                                help='Maximum number of proteins to process')
        batch_parser.add_argument('--reps-only', action='store_true',
                                help='Process only representative proteins')
        batch_parser.add_argument('--force', action='store_true',
                                help='Force processing even if batch is not ready')
        batch_parser.add_argument('--workers', type=int, default=4,
                                help='Number of parallel workers')

        # Partition specific proteins
        specific_parser = partition_subparsers.add_parser('specific', help='Process specific proteins')
        specific_parser.add_argument('--process-ids', type=int, nargs='+', required=True,
                                   help='Process IDs to process')
        specific_parser.add_argument('--batch-id', type=int,
                                   help='Batch ID (optional)')
        specific_parser.add_argument('--blast-only', action='store_true',
                                   help='Use only BLAST results')

        # Partition all batches
        all_parser = partition_subparsers.add_parser('all', help='Process multiple batches')
        all_parser.add_argument('--batch-ids', type=int, nargs='+',
                              help='Batch IDs to process (default: all batches)')
        all_parser.add_argument('--exclude-batch-ids', type=int, nargs='+', default=[],
                              help='Batch IDs to exclude')
        all_parser.add_argument('--reference', type=str,
                              help='Filter batches by reference version')
        all_parser.add_argument('--blast-only', action='store_true',
                              help='Use only BLAST results (no HHSearch)')
        all_parser.add_argument('--limit-per-batch', type=int,
                              help='Maximum number of proteins to process per batch')
        all_parser.add_argument('--reps-only', action='store_true',
                              help='Process only representative proteins')
        all_parser.add_argument('--batch-size', type=int, default=5,
                              help='Number of batches to process simultaneously')
        all_parser.add_argument('--wait-between-groups', type=int, default=30,
                              help='Seconds to wait between batch groups')
        all_parser.add_argument('--force', action='store_true',
                              help='Force processing even if batches are not ready')
        all_parser.add_argument('--workers', type=int, default=1,
                              help='Number of parallel workers per batch')

        # SLURM-specific options for 'all' action
        all_parser.add_argument('--use-slurm', action='store_true',
                              help='Submit jobs to SLURM instead of running directly')
        all_parser.add_argument('--slurm-threads', type=int, default=8,
                              help='Number of threads to request per SLURM job')
        all_parser.add_argument('--slurm-memory', type=str, default='16G',
                              help='Memory to request per SLURM job')
        all_parser.add_argument('--slurm-time', type=str, default='12:00:00',
                              help='Time limit for SLURM jobs')
        all_parser.add_argument('--wait-for-completion', action='store_true',
                              help='Wait for SLURM jobs to complete')
        all_parser.add_argument('--check-interval', type=int, default=60,
                              help='Seconds between job status checks')
        all_parser.add_argument('--timeout', type=int,
                              help='Maximum time to wait for completion (seconds)')

    def _setup_analyze_parser(self, subparsers):
        """Set up analyze command parser"""
        analyze_parser = subparsers.add_parser('analyze', help=COMMANDS['analyze'])
        analyze_subparsers = analyze_parser.add_subparsers(dest='action', help='Analysis action')

        # Analyze batch status
        status_parser = analyze_subparsers.add_parser('status', help='Check batch status')
        status_parser.add_argument('--batch-ids', type=int, nargs='+',
                                 help='Batch IDs to check (default: all batches)')
        status_parser.add_argument('--blast-only', action='store_true',
                                 help='Check blast-only status')

        # Analyze protein status
        protein_parser = analyze_subparsers.add_parser('protein', help='Check protein status')
        protein_parser.add_argument('--pdb-id', type=str, required=True,
                                  help='PDB ID')
        protein_parser.add_argument('--chain-id', type=str, required=True,
                                  help='Chain ID')
        protein_parser.add_argument('--batch-id', type=int,
                                  help='Batch ID (optional)')

        # Analyze domain counts
        counts_parser = analyze_subparsers.add_parser('counts', help='Analyze domain statistics')
        counts_parser.add_argument('--batch-ids', type=int, nargs='+',
                                 help='Batch IDs to analyze (default: all batches)')
        counts_parser.add_argument('--sample-size', type=int, default=50,
                                 help='Number of proteins to sample per batch')

        # Analyze failures
        failures_parser = analyze_subparsers.add_parser('failures', help='Analyze failure patterns')
        failures_parser.add_argument('--batch-id', type=int,
                                   help='Limit to specific batch')
        failures_parser.add_argument('--limit', type=int, default=20,
                                   help='Maximum failures to analyze')

    def _setup_diagnose_parser(self, subparsers):
        """Set up diagnose command parser"""
        diagnose_parser = subparsers.add_parser('diagnose', help=COMMANDS['diagnose'])
        diagnose_subparsers = diagnose_parser.add_subparsers(dest='action', help='Diagnostic action')

        # Diagnose process
        process_parser = diagnose_subparsers.add_parser('process', help='Diagnose a specific process')
        process_parser.add_argument('--process-id', type=int, required=True,
                                  help='Process ID to diagnose')
        process_parser.add_argument('--fix-in-db', action='store_true',
                                  help='Fix file paths in database if legacy paths are found')

        # Diagnose batch
        batch_parser = diagnose_subparsers.add_parser('batch', help='Diagnose an entire batch')
        batch_parser.add_argument('--batch-id', type=int, required=True,
                                help='Batch ID to diagnose')
        batch_parser.add_argument('--sample-size', type=int, default=10,
                                help='Number of problematic proteins to sample')

        # Check readiness
        readiness_parser = diagnose_subparsers.add_parser('readiness', help='Check batch readiness')
        readiness_parser.add_argument('--batch-id', type=int, required=True,
                                    help='Batch ID to check')
        readiness_parser.add_argument('--detailed', action='store_true',
                                    help='Show detailed protein list')

    def _setup_repair_parser(self, subparsers):
        """Set up repair command parser"""
        repair_parser = subparsers.add_parser('repair', help=COMMANDS['repair'])
        repair_subparsers = repair_parser.add_subparsers(dest='action', help='Repair action')

        # Repair failed processes
        failed_parser = repair_subparsers.add_parser('failed', help='Reset and retry failed processes')
        failed_parser.add_argument('--batch-id', type=int, required=True,
                                 help='Batch ID to repair')
        failed_parser.add_argument('--rerun', action='store_true',
                                 help='Rerun partition after resetting')
        failed_parser.add_argument('--blast-only', action='store_true',
                                 help='Use only BLAST results for rerun')

        # Repair missing files
        missing_parser = repair_subparsers.add_parser('missing', help='Regenerate missing files')
        missing_parser.add_argument('--batch-id', type=int, required=True,
                                  help='Batch ID to repair')
        missing_parser.add_argument('--blast-only', action='store_true',
                                  help='Use only BLAST results')
        missing_parser.add_argument('--limit', type=int,
                                  help='Maximum number of files to repair')

    def _setup_monitor_parser(self, subparsers):
        """Set up monitor command parser"""
        monitor_parser = subparsers.add_parser('monitor', help=COMMANDS['monitor'])
        monitor_parser.add_argument('--batch-id', type=int, required=True,
                                  help='Batch ID to monitor')
        monitor_parser.add_argument('--interval', type=int, default=60,
                                  help='Check interval in seconds')
        monitor_parser.add_argument('--timeout', type=int,
                                  help='Timeout in seconds')

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
            return self._run_analyze(args)
        elif args.command == 'diagnose':
            return self._run_diagnose(args)
        elif args.command == 'repair':
            return self._run_repair(args)
        elif args.command == 'monitor':
            return self._run_monitor(args)
        else:
            self.logger.error(f"Unknown command: {args.command}")
            return 1

    # Summary command implementation
    def _run_summary(self, args: argparse.Namespace) -> int:
        """Run domain summary creation"""
        domain_summary = DomainSummary(self.context)

        if args.batch_id:
            # Process batch
            self.logger.info(f"Creating domain summaries for batch {args.batch_id}")

            # Get batch info
            batch_info = self._get_batch_info(args.batch_id)
            if not batch_info:
                return 1

            success = domain_summary.process_batch(
                args.batch_id,
                batch_info['base_path'],
                args.reference or batch_info['ref_version'],
                args.blast_only,
                args.limit,
                args.reps_only
            )

            return 0 if success else 1

        elif args.pdb_id and args.chain_id:
            # Process single protein
            self.logger.info(f"Creating domain summary for {args.pdb_id}_{args.chain_id}")

            # Implementation for single protein
            # [Would need to implement single protein processing in DomainSummary]

            self.logger.error("Single protein summary not yet implemented")
            return 1

        else:
            self.logger.error("Either --batch-id or both --pdb-id and --chain-id must be specified")
            return 1

    # Partition command implementations
    def _run_partition(self, args: argparse.Namespace) -> int:
        """Run domain partition based on action"""
        if args.action == 'single':
            return self._partition_single(args)
        elif args.action == 'batch':
            return self._partition_batch(args)
        elif args.action == 'specific':
            return self._partition_specific(args)
        elif args.action == 'all':
            return self._partition_all(args)
        else:
            self.logger.error(f"Unknown partition action: {args.action}")
            return 1

    def _partition_single(self, args: argparse.Namespace) -> int:
        """Process domain partition for a single protein"""
        if not args.pdb_id or not args.chain_id:
            self.logger.error("PDB ID and chain ID are required for single protein processing")
            return 1

        # Get batch information if batch ID provided
        batch_path = args.batch_path
        reference = args.reference
        process_id = None

        if args.batch_id:
            batch_info = self._get_batch_info(args.batch_id)
            if not batch_info:
                self.logger.error(f"Batch {args.batch_id} not found")
                return 1

            batch_path = batch_info["base_path"]
            reference = batch_info["ref_version"]

            # Look up process ID for status tracking
            process_id = self._get_process_id_for_protein(args.pdb_id, args.chain_id, args.batch_id)

        if not batch_path or not reference:
            self.logger.error("Batch path and reference are required - specify directly or provide batch ID")
            return 1

        # Check if domain summary exists
        if not self._verify_domain_summary(batch_path, args.pdb_id, args.chain_id, reference, args.blast_only):
            self.logger.error(f"Domain summary not found for {args.pdb_id}_{args.chain_id}")
            self.logger.error("Run domain summary first or provide correct paths")
            return 1

        # Process single protein
        try:
            self.logger.info(f"Processing {args.pdb_id}_{args.chain_id} with service-based architecture")

            # Find domain summary path
            summary_path = self._find_domain_summary(
                batch_path, args.pdb_id, args.chain_id, reference, args.blast_only
            )

            if not summary_path:
                self.logger.error(f"Domain summary file not found for {args.pdb_id}_{args.chain_id}")
                return 1

            # Process using service
            result = self.service.partition_protein(
                pdb_id=args.pdb_id,
                chain_id=args.chain_id,
                summary_path=summary_path,
                output_dir=batch_path,
                process_id=process_id,
                blast_only=args.blast_only
            )

            # Report outcome
            if result.success:
                self.logger.info(f"Successfully processed {args.pdb_id}_{args.chain_id}")

                if result.domain_file:
                    self.logger.info(f"Created domain file: {result.domain_file}")

                if result.is_peptide:
                    self.logger.info("Classified as peptide")
                elif result.is_classified:
                    self.logger.info(f"Classified with {len(result.domains)} domains")

                    # Show domain details
                    for i, domain in enumerate(result.domains):
                        self.logger.info(f"  Domain {i+1}: {domain.range} ({domain.get_classification_level()})")
                else:
                    self.logger.info("Unclassified")

                # Show coverage statistics
                if result.sequence_length > 0:
                    self.logger.info(f"Coverage: {result.coverage:.1%} ({result.residues_assigned}/{result.sequence_length} residues)")

            else:
                self.logger.error(f"Failed to process {args.pdb_id}_{args.chain_id}: {result.error}")

            return 0 if result.success else 1

        except Exception as e:
            self.logger.error(f"Error processing {args.pdb_id}_{args.chain_id}: {str(e)}")
            import traceback
            self.logger.error(traceback.format_exc())
            return 1

    def _partition_batch(self, args: argparse.Namespace) -> int:
        """Process domain partition for a batch"""
        self.logger.info(f"Processing batch {args.batch_id} with service-based architecture")

        # Get batch information
        batch_info = self._get_batch_info(args.batch_id)
        if not batch_info:
            self.logger.error(f"Batch {args.batch_id} not found")
            return 1

        # Verify batch readiness if not forcing
        if not args.force and not self._verify_batch_readiness(args.batch_id, args.blast_only, args.reps_only):
            self.logger.error(f"Batch {args.batch_id} is not ready for domain partition")
            self.logger.error("Use --force to override or run domain summary first")
            return 1

        # Process the batch
        batch_path = batch_info["base_path"]

        try:
            self.logger.info(f"Running service-based partition for batch {args.batch_id}")

            # Create partition options
            options = {
                'blast_only': args.blast_only,
                'representatives_only': args.reps_only,
                'force_overwrite': args.force,
                'validation_level': ValidationLevel.NORMAL
            }

            # Process batch
            results = self.service.partition_batch(
                batch_id=args.batch_id,
                batch_path=batch_path,
                limit=args.limit,
                **options
            )

            # Log results
            self.logger.info(f"Domain partition complete: {results.success_count} succeeded, {results.failure_count} failed")

            # Show summary statistics
            if results.success_count > 0:
                self.logger.info("Classification summary:")
                self.logger.info(f"  Classified: {results.proteins_with_domains}")
                self.logger.info(f"  Peptides: {results.peptides_found}")
                self.logger.info(f"  Unclassified: {results.unclassified_proteins}")
                self.logger.info(f"  Total domains: {results.total_domains_found}")

                if results.proteins_with_domains > 0:
                    avg_domains = results.total_domains_found / results.proteins_with_domains
                    self.logger.info(f"  Average domains per classified protein: {avg_domains:.1f}")

            # Log failures if any
            if results.failure_count > 0:
                for i, (pdb_id, chain_id, error) in enumerate(results.failures[:3]):
                    self.logger.error(f"Failure {i+1}: {pdb_id}_{chain_id}: {error}")

                if results.failure_count > 3:
                    self.logger.error(f"... and {results.failure_count - 3} more failures")

            # Get service statistics
            stats = self.service.get_service_statistics()
            self.logger.info(f"Processing statistics: {stats['service']['proteins_processed']} proteins in "
                           f"{stats['service']['runtime_seconds']:.1f}s "
                           f"({stats['service']['proteins_per_minute']:.1f} proteins/min)")

            return 0 if results.success_count > 0 else 1

        except Exception as e:
            self.logger.error(f"Error processing batch: {str(e)}")
            import traceback
            self.logger.error(traceback.format_exc())
            return 1

    def _partition_specific(self, args: argparse.Namespace) -> int:
        """Process domain partition for specific proteins"""
        # Get batch ID if not specified
        batch_id = args.batch_id
        if not batch_id and args.process_ids:
            batch_id = self._get_batch_id_for_process(args.process_ids[0])

            if not batch_id:
                self.logger.error(f"Could not determine batch ID for process {args.process_ids[0]}")
                return 1

            self.logger.info(f"Using batch ID {batch_id} for process {args.process_ids[0]}")

        # Get batch information
        batch_info = self._get_batch_info(batch_id)
        if not batch_info:
            self.logger.error(f"Batch {batch_id} not found")
            return 1

        # Get protein details for each process ID
        proteins = []
        for process_id in args.process_ids:
            protein_info = self._get_protein_for_process(process_id)
            if protein_info:
                proteins.append(protein_info)
            else:
                self.logger.warning(f"Process ID {process_id} not found")

        if not proteins:
            self.logger.error("No valid proteins found for specified process IDs")
            return 1

        # Process each protein
        success_count = 0
        failure_count = 0

        try:
            self.logger.info(f"Processing {len(proteins)} specific proteins from batch {batch_id}")

            for protein in proteins:
                pdb_id = protein['pdb_id']
                chain_id = protein['chain_id']
                process_id = protein['process_id']

                # Find domain summary
                summary_path = self._find_domain_summary(
                    batch_info["base_path"], pdb_id, chain_id,
                    batch_info["ref_version"], args.blast_only
                )

                if not summary_path:
                    self.logger.error(f"Domain summary not found for {pdb_id}_{chain_id}")
                    failure_count += 1
                    continue

                # Process protein
                result = self.service.partition_protein(
                    pdb_id=pdb_id,
                    chain_id=chain_id,
                    summary_path=summary_path,
                    output_dir=batch_info["base_path"],
                    process_id=process_id,
                    blast_only=args.blast_only
                )

                if result.success:
                    success_count += 1
                    self.logger.info(f"Successfully processed {pdb_id}_{chain_id}")
                else:
                    failure_count += 1
                    self.logger.error(f"Failed to process {pdb_id}_{chain_id}: {result.error}")

            self.logger.info(f"Domain partition complete: {success_count} succeeded, {failure_count} failed")
            return 0 if success_count > 0 else 1

        except Exception as e:
            self.logger.error(f"Error processing specific proteins: {str(e)}")
            import traceback
            self.logger.error(traceback.format_exc())
            return 1

    def _partition_all(self, args: argparse.Namespace) -> int:
        """Process domain partition for all batches"""
        # Determine which batches to process
        batch_ids = args.batch_ids
        if not batch_ids:
            # Get all batches if not specified
            all_batches = self._get_all_batch_ids()
            self.logger.info(f"Found {len(all_batches)} total batches")

            # Filter by reference version if specified
            if args.reference:
                reference_batches = []
                for batch_id in all_batches:
                    batch_info = self._get_batch_info(batch_id)
                    if batch_info and batch_info.get("ref_version") == args.reference:
                        reference_batches.append(batch_id)

                batch_ids = reference_batches
                self.logger.info(f"Filtered to {len(batch_ids)} batches with reference {args.reference}")
            else:
                batch_ids = all_batches

        # Exclude specified batch IDs
        if args.exclude_batch_ids:
            batch_ids = [b_id for b_id in batch_ids if b_id not in args.exclude_batch_ids]
            self.logger.info(f"Excluding specified batches, {len(batch_ids)} batches remaining")

        # Sort batch IDs for consistent processing
        batch_ids.sort()

        if not batch_ids:
            self.logger.warning("No batches to process")
            return 0

        self.logger.info(f"Processing {len(batch_ids)} batches with service-based architecture: {batch_ids}")

        # Determine whether to use SLURM or process directly
        if args.use_slurm:
            return self._submit_batches_to_slurm(batch_ids, args)
        else:
            return self._process_batches_directly(batch_ids, args)

    def _submit_batches_to_slurm(self, batch_ids: List[int], args: argparse.Namespace) -> int:
        """Submit batch processing jobs to SLURM"""
        self.logger.info(f"Submitting {len(batch_ids)} batches to SLURM")

        # Get job manager
        try:
            from ecod.jobs import SlurmJobManager
            job_manager = SlurmJobManager(self.context.config_manager.config)
        except ImportError:
            self.logger.error("Failed to import SlurmJobManager. Make sure SLURM integration is available.")
            return 1

        # Create temporary directory for job scripts
        temp_dir = os.path.join(self.context.config_manager.get_path('output_dir', '/tmp'),
                               f"domain_partition_jobs_{int(time.time())}")
        os.makedirs(temp_dir, exist_ok=True)
        self.logger.info(f"Created job directory: {temp_dir}")

        # Get paths
        config_path = os.path.abspath(self.context.config_manager.config_path)

        # Process batches with appropriate batch size
        batch_groups = []
        for i in range(0, len(batch_ids), args.batch_size):
            batch_groups.append(batch_ids[i:i+args.batch_size])

        self.logger.info(f"Split into {len(batch_groups)} groups with max {args.batch_size} batches per group")

        # Submit jobs
        job_ids = []

        for group_idx, group in enumerate(batch_groups):
            group_dir = os.path.join(temp_dir, f"group_{group_idx}")
            os.makedirs(group_dir, exist_ok=True)

            for batch_id in group:
                # Create job name
                job_name = f"domain_partition_batch_{batch_id}"

                # Build command using CLI
                command = f"ecod --config {config_path} domain partition batch --batch-id {batch_id}"

                if args.blast_only:
                    command += " --blast-only"
                if args.limit_per_batch:
                    command += f" --limit {args.limit_per_batch}"
                if args.reps_only:
                    command += " --reps-only"
                if args.force:
                    command += " --force"

                # Create and submit job
                script_path_gen = job_manager.create_job_script(
                    commands=[command],
                    job_name=job_name,
                    output_dir=group_dir,
                    threads=args.slurm_threads,
                    memory=args.slurm_memory,
                    time=args.slurm_time
                )

                job_id = job_manager.submit_job(script_path_gen)

                if job_id:
                    job_ids.append(job_id)
                    self.logger.info(f"Submitted job for batch {batch_id}, SLURM job ID: {job_id}")
                else:
                    self.logger.error(f"Failed to submit job for batch {batch_id}")

            # Wait between batch groups if specified
            if args.wait_between_groups and group_idx < len(batch_groups) - 1:
                self.logger.info(f"Waiting {args.wait_between_groups} seconds before next batch group")
                time.sleep(args.wait_between_groups)

        self.logger.info(f"Submitted {len(job_ids)} jobs to SLURM")

        # Monitor jobs if requested
        if args.wait_for_completion:
            return self._monitor_slurm_jobs(job_manager, job_ids, args)

        return 0

    def _process_batches_directly(self, batch_ids: List[int], args: argparse.Namespace) -> int:
        """Process batches directly using the service"""
        # Process batches with appropriate batch size
        batch_groups = []
        for i in range(0, len(batch_ids), args.batch_size):
            batch_groups.append(batch_ids[i:i+args.batch_size])

        self.logger.info(f"Split into {len(batch_groups)} groups with max {args.batch_size} batches per group")

        # Track overall statistics
        overall_results = BatchPartitionResults()
        success_count = 0
        failed_batches = []

        for group_idx, group in enumerate(batch_groups):
            self.logger.info(f"Processing batch group {group_idx+1}/{len(batch_groups)}: {group}")

            group_success = 0

            for batch_id in group:
                batch_info = self._get_batch_info(batch_id)
                if not batch_info:
                    self.logger.error(f"Batch {batch_id} not found")
                    failed_batches.append(batch_id)
                    continue

                self.logger.info(f"Processing batch {batch_id} ({batch_info.get('batch_name', '')})")

                # Check batch readiness if not forcing
                if not args.force and not self._verify_batch_readiness(batch_id, args.blast_only, args.reps_only):
                    self.logger.warning(f"Batch {batch_id} is not ready for domain partition, skipping")
                    failed_batches.append(batch_id)
                    continue

                # Process the batch
                try:
                    batch_path = batch_info["base_path"]

                    # Create options
                    options = {
                        'blast_only': args.blast_only,
                        'representatives_only': args.reps_only,
                        'force_overwrite': args.force
                    }

                    results = self.service.partition_batch(
                        batch_id=batch_id,
                        batch_path=batch_path,
                        limit=args.limit_per_batch,
                        **options
                    )

                    # Merge results
                    overall_results.total += results.total
                    overall_results.success_count += results.success_count
                    overall_results.failure_count += results.failure_count
                    overall_results.proteins_with_domains += results.proteins_with_domains
                    overall_results.total_domains_found += results.total_domains_found
                    overall_results.peptides_found += results.peptides_found
                    overall_results.unclassified_proteins += results.unclassified_proteins

                    self.logger.info(f"Batch {batch_id} complete: {results.success_count} succeeded, {results.failure_count} failed")

                    if results.success_count > 0:
                        group_success += 1
                        success_count += 1
                    else:
                        failed_batches.append(batch_id)

                except Exception as e:
                    self.logger.error(f"Error processing batch {batch_id}: {str(e)}")
                    import traceback
                    self.logger.error(traceback.format_exc())
                    failed_batches.append(batch_id)

            self.logger.info(f"Batch group {group_idx+1} complete: {group_success}/{len(group)} succeeded")

            # Wait between batch groups if specified
            if args.wait_between_groups and group_idx < len(batch_groups) - 1:
                self.logger.info(f"Waiting {args.wait_between_groups} seconds before next batch group")
                time.sleep(args.wait_between_groups)

        # Log final summary
        overall_results.finalize()
        summary = overall_results.get_summary()

        self.logger.info(f"All batches processed: {success_count}/{len(batch_ids)} succeeded")
        self.logger.info("Overall statistics:")
        self.logger.info(f"  Total proteins processed: {summary['total_proteins']}")
        self.logger.info(f"  Successful: {summary['successful']}")
        self.logger.info(f"  Failed: {summary['failed']}")
        self.logger.info(f"  Classified: {summary['proteins_with_domains']}")
        self.logger.info(f"  Peptides: {summary['peptides']}")
        self.logger.info(f"  Unclassified: {summary['unclassified']}")
        self.logger.info(f"  Total domains: {summary['total_domains']}")
        self.logger.info(f"  Processing time: {summary['processing_time']:.1f}s")

        if summary['proteins_with_domains'] > 0:
            avg_domains = summary['total_domains'] / summary['proteins_with_domains']
            self.logger.info(f"  Average domains per classified protein: {avg_domains:.1f}")

        if failed_batches:
            self.logger.warning(f"Failed batches: {failed_batches}")

        # Clear caches after processing
        self.service.clear_all_caches()

        return 0 if len(failed_batches) == 0 else 1

    # Analyze command implementations
    def _run_analyze(self, args: argparse.Namespace) -> int:
        """Run analysis based on action"""
        if args.action == 'status':
            return self._analyze_batch_status(args)
        elif args.action == 'protein':
            return self._analyze_protein_status(args)
        elif args.action == 'counts':
            return self._analyze_domain_counts(args)
        elif args.action == 'failures':
            return self._analyze_failures(args)
        else:
            self.logger.error(f"Unknown analyze action: {args.action}")
            return 1

    def _analyze_batch_status(self, args: argparse.Namespace) -> int:
        """Analyze batch status for domain partition"""
        # Determine which batches to analyze
        batch_ids = args.batch_ids
        if not batch_ids:
            batch_ids = self._get_all_batch_ids()
            self.logger.info(f"Analyzing status for all {len(batch_ids)} batches")
        else:
            self.logger.info(f"Analyzing status for {len(batch_ids)} specified batches")

        # Create results table
        results = []

        for batch_id in batch_ids:
            batch_info = self._get_batch_info(batch_id)
            if not batch_info:
                self.logger.warning(f"Batch {batch_id} not found")
                continue

            # Get batch progress using service's tracker
            progress = self.service.tracker.get_batch_progress(batch_id)

            # Calculate readiness and completion
            total = progress.get("total", 0)
            ready = progress.get("files_created", 0)  # Files ready for partition
            complete = progress.get("complete", 0)

            readiness = (ready / total * 100) if total > 0 else 0
            completion = progress.get("complete_pct", 0)

            results.append({
                "batch_id": batch_id,
                "batch_name": batch_info.get("batch_name", ""),
                "total": total,
                "ready": ready,
                "complete": complete,
                "readiness": readiness,
                "completion": completion,
                "status": batch_info.get("status", ""),
                "errors": progress.get("errors", 0)
            })

        # Sort by batch ID
        results.sort(key=lambda x: x["batch_id"])

        # Display results
        self.logger.info(f"{'ID':5s} {'Name':20s} {'Total':8s} {'Ready':8s} {'Complete':8s} {'Ready%':8s} {'Complete%':8s} {'Errors':6s} {'Status':12s}")
        self.logger.info("-" * 100)

        for r in results:
            self.logger.info(f"{r['batch_id']:5d} {r['batch_name'][:20]:20s} {r['total']:8d} " +
                            f"{r['ready']:8d} {r['complete']:8d} {r['readiness']:7.1f}% " +
                            f"{r['completion']:8.1f}% {r['errors']:6d} {r['status'][:12]:12s}")

        self.logger.info("-" * 100)

        # Calculate totals
        total_proteins = sum(r["total"] for r in results)
        total_ready = sum(r["ready"] for r in results)
        total_complete = sum(r["complete"] for r in results)
        total_errors = sum(r["errors"] for r in results)

        self.logger.info(f"Total: {len(results)} batches, {total_proteins} proteins, " +
                        f"{total_ready} ready ({total_ready/total_proteins*100:.1f}%), " +
                        f"{total_complete} complete ({total_complete/total_proteins*100:.1f}%), " +
                        f"{total_errors} errors")

        return 0

    def _analyze_protein_status(self, args: argparse.Namespace) -> int:
        """Analyze protein status"""
        if not args.pdb_id or not args.chain_id:
            self.logger.error("PDB ID and chain ID are required for protein status analysis")
            return 1

        # Find all occurrences of this protein
        proteins = self._find_protein_in_database(args.pdb_id, args.chain_id, args.batch_id)

        if not proteins:
            self.logger.error(f"Protein {args.pdb_id}_{args.chain_id} not found")
            return 1

        # Display protein information
        self.logger.info(f"Found {len(proteins)} instances of {args.pdb_id}_{args.chain_id}")

        for i, p in enumerate(proteins):
            batch_info = self._get_batch_info(p["batch_id"])
            batch_name = batch_info.get("batch_name", "unknown") if batch_info else "unknown"

            self.logger.info(f"Instance {i+1}:")
            self.logger.info(f"  Batch: {p['batch_id']} ({batch_name})")
            self.logger.info(f"  Process ID: {p['process_id']}")
            self.logger.info(f"  Status: {p['current_stage']} / {p['status']}")

            if p.get("is_representative"):
                self.logger.info(f"  Is Representative: Yes")

            if p["error_message"]:
                self.logger.info(f"  Error: {p['error_message']}")

            # Show file information
            files = p.get("files", {})
            self.logger.info("  Files:")
            for file_type, file_info in files.items():
                exists = "EXISTS" if file_info.get("file_exists", False) else "MISSING"
                self.logger.info(f"    {file_type:20s}: {exists}")

                if args.verbose and file_info.get("file_exists", False):
                    file_path = file_info.get("file_path", "")
                    if file_path and batch_info:
                        full_path = os.path.join(batch_info.get("base_path", ""), file_path)
                        self.logger.info(f"      Path: {full_path}")

            # Show domain information if available
            if batch_info and p.get("status") == "success" and "domain_partition" in files:
                domain_info = self._check_domain_result(p["process_id"], args.pdb_id, args.chain_id, batch_info)
                if domain_info:
                    self.logger.info("  Domain Information:")
                    self.logger.info(f"    Is Classified: {domain_info.get('is_classified', False)}")
                    self.logger.info(f"    Is Peptide: {domain_info.get('is_peptide', False)}")
                    self.logger.info(f"    Domain Count: {len(domain_info.get('domains', []))}")
                    self.logger.info(f"    Coverage: {domain_info.get('coverage', 0):.1%}")
                    self.logger.info(f"    Processing Time: {domain_info.get('processing_time', 0):.2f}s")

        return 0

    def _analyze_domain_counts(self, args: argparse.Namespace) -> int:
        """Analyze domain count statistics"""
        # Determine which batches to analyze
        batch_ids = args.batch_ids
        if not batch_ids:
            batch_ids = self._get_all_batch_ids()
            self.logger.info(f"Analyzing domains for all {len(batch_ids)} batches")
        else:
            self.logger.info(f"Analyzing domains for {len(batch_ids)} specified batches")

        # Process each batch
        results = []

        for batch_id in batch_ids:
            batch_info = self._get_batch_info(batch_id)
            if not batch_info:
                continue

            # Get sample files for analysis
            sample_query = """
            SELECT p.pdb_id, p.chain_id, pf.file_path
            FROM ecod_schema.process_status ps
            JOIN ecod_schema.protein p ON ps.protein_id = p.id
            JOIN ecod_schema.process_file pf ON ps.id = pf.process_id
            WHERE pf.file_type = 'domain_partition'
              AND ps.batch_id = %s
              AND ps.status = 'success'
              AND pf.file_exists = TRUE
            ORDER BY RANDOM()
            LIMIT %s
            """

            samples = self.db.execute_dict_query(sample_query, (batch_id, args.sample_size))

            # Analyze using models
            domain_stats = {
                "total_domains": 0,
                "proteins_with_domains": 0,
                "multi_domain_proteins": 0,
                "single_domain_proteins": 0,
                "unclassified_proteins": 0,
                "peptide_proteins": 0,
                "domain_counts": {},
                "confidence_stats": [],
                "evidence_stats": []
            }

            for sample in samples:
                full_path = os.path.join(batch_info["base_path"], sample["file_path"])

                if os.path.exists(full_path):
                    try:
                        # Parse using DomainPartitionResult
                        result = DomainPartitionResult.from_xml_file(full_path)

                        if result.is_peptide:
                            domain_stats["peptide_proteins"] += 1
                        elif result.is_unclassified:
                            domain_stats["unclassified_proteins"] += 1
                        elif result.is_classified:
                            domain_count = len(result.domains)

                            if domain_count > 0:
                                domain_stats["proteins_with_domains"] += 1
                                domain_stats["total_domains"] += domain_count

                                # Track distribution
                                if domain_count not in domain_stats["domain_counts"]:
                                    domain_stats["domain_counts"][domain_count] = 0
                                domain_stats["domain_counts"][domain_count] += 1

                                # Track multi-domain proteins
                                if domain_count > 1:
                                    domain_stats["multi_domain_proteins"] += 1
                                else:
                                    domain_stats["single_domain_proteins"] += 1

                                # Collect confidence and evidence statistics
                                for domain in result.domains:
                                    if hasattr(domain, 'confidence'):
                                        domain_stats["confidence_stats"].append(domain.confidence)
                                    if hasattr(domain, 'evidence'):
                                        domain_stats["evidence_stats"].append(len(domain.evidence))

                    except Exception as e:
                        self.logger.debug(f"Error parsing domain file {full_path}: {str(e)}")

            # Calculate averages
            avg_domains = (domain_stats["total_domains"] / domain_stats["proteins_with_domains"]
                          if domain_stats["proteins_with_domains"] > 0 else 0)

            avg_confidence = (sum(domain_stats["confidence_stats"]) / len(domain_stats["confidence_stats"])
                            if domain_stats["confidence_stats"] else 0)

            avg_evidence = (sum(domain_stats["evidence_stats"]) / len(domain_stats["evidence_stats"])
                          if domain_stats["evidence_stats"] else 0)

            # Add to results
            results.append({
                "batch_id": batch_id,
                "batch_name": batch_info.get("batch_name", ""),
                "sample_size": len(samples),
                "proteins_with_domains": domain_stats["proteins_with_domains"],
                "multi_domain_proteins": domain_stats["multi_domain_proteins"],
                "single_domain_proteins": domain_stats["single_domain_proteins"],
                "unclassified_proteins": domain_stats["unclassified_proteins"],
                "peptide_proteins": domain_stats["peptide_proteins"],
                "avg_domains": avg_domains,
                "avg_confidence": avg_confidence,
                "avg_evidence": avg_evidence,
                "domain_counts": domain_stats["domain_counts"]
            })

        # Display results
        self.logger.info(f"{'ID':5s} {'Name':20s} {'Sample':6s} {'w/Domains':9s} {'Multi':5s} {'Peptides':8s} {'Avg Dom':7s} {'Avg Conf':8s} {'Avg Evid':8s}")
        self.logger.info("-" * 95)

        for r in results:
            self.logger.info(f"{r['batch_id']:5d} {r['batch_name'][:20]:20s} {r['sample_size']:6d} " +
                            f"{r['proteins_with_domains']:9d} " +
                            f"{r['multi_domain_proteins']:5d} {r['peptide_proteins']:8d} " +
                            f"{r['avg_domains']:7.1f} {r['avg_confidence']:8.3f} {r['avg_evidence']:8.1f}")

        self.logger.info("-" * 95)

        # Combined statistics
        if args.verbose:
            combined_counts = {}
            for r in results:
                for count, num in r["domain_counts"].items():
                    if count not in combined_counts:
                        combined_counts[count] = 0
                    combined_counts[count] += num

            self.logger.info("Domain count distribution:")
            for count in sorted(combined_counts.keys()):
                self.logger.info(f"  {count} domains: {combined_counts[count]} proteins")

        return 0

    def _analyze_failures(self, args: argparse.Namespace) -> int:
        """Analyze common failure patterns"""
        self.logger.info("Analyzing domain partition failures")

        # Get failed processes
        query = """
        SELECT ps.id as process_id, p.pdb_id, p.chain_id, ps.batch_id,
               ps.current_stage, ps.error_message, b.base_path
        FROM ecod_schema.process_status ps
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        JOIN ecod_schema.batch b ON ps.batch_id = b.id
        WHERE ps.status = 'error'
          AND ps.current_stage LIKE '%domain_partition%'
        """

        params = []
        if args.batch_id:
            query += " AND ps.batch_id = %s"
            params.append(args.batch_id)

        query += f" ORDER BY ps.updated_at DESC LIMIT {args.limit}"

        failures = self.db.execute_dict_query(query, tuple(params))

        if not failures:
            self.logger.info("No domain partition failures found")
            return 0

        # Analyze failure patterns
        failure_patterns = {}

        for failure in failures:
            error_msg = failure['error_message'] or 'Unknown error'

            # Categorize errors
            if 'summary not found' in error_msg.lower():
                category = 'Missing domain summary'
            elif 'no evidence' in error_msg.lower():
                category = 'No evidence found'
            elif 'validation' in error_msg.lower():
                category = 'Validation error'
            elif 'xml' in error_msg.lower() or 'parse' in error_msg.lower():
                category = 'XML parsing error'
            elif 'permission' in error_msg.lower() or 'access' in error_msg.lower():
                category = 'File access error'
            else:
                category = 'Other'

            if category not in failure_patterns:
                failure_patterns[category] = []

            failure_patterns[category].append(failure)

        # Display analysis
        self.logger.info(f"Analyzed {len(failures)} failures:")

        for category, cases in sorted(failure_patterns.items(), key=lambda x: len(x[1]), reverse=True):
            self.logger.info(f"\n{category}: {len(cases)} cases")

            # Show examples
            for i, case in enumerate(cases[:3]):  # Show first 3 examples
                self.logger.info(f"  {i+1}. {case['pdb_id']}_{case['chain_id']} "
                               f"(Batch {case['batch_id']})")
                if case['error_message']:
                    error_preview = case['error_message'][:150]
                    if len(case['error_message']) > 150:
                        error_preview += "..."
                    self.logger.info(f"     {error_preview}")

        # Provide recommendations
        self.logger.info("\nRecommendations:")

        if 'Missing domain summary' in failure_patterns:
            count = len(failure_patterns['Missing domain summary'])
            self.logger.info(f"  - {count} proteins need domain summary generation")
            self.logger.info("    Run: ecod domain summary --batch-id ...")

        if 'No evidence found' in failure_patterns:
            count = len(failure_patterns['No evidence found'])
            self.logger.info(f"  - {count} proteins have no valid evidence")
            self.logger.info("    These may be truly unclassifiable or need BLAST re-run")

        if 'XML parsing error' in failure_patterns:
            count = len(failure_patterns['XML parsing error'])
            self.logger.info(f"  - {count} proteins have corrupted XML files")
            self.logger.info("    Run: ecod domain repair missing --batch-id ...")

        return 0

    # Diagnose command implementations
    def _run_diagnose(self, args: argparse.Namespace) -> int:
        """Run diagnostics based on action"""
        if args.action == 'process':
            return self._diagnose_process(args)
        elif args.action == 'batch':
            return self._diagnose_batch(args)
        elif args.action == 'readiness':
            return self._diagnose_readiness(args)
        else:
            self.logger.error(f"Unknown diagnose action: {args.action}")
            return 1

    def _diagnose_process(self, args: argparse.Namespace) -> int:
        """Diagnose domain partition issues for a specific process"""
        from ecod.cli.diagnostics import DomainPartitionDiagnostics

        diagnostics = DomainPartitionDiagnostics(self.context)
        success = diagnostics.diagnose_process(args.process_id, args.fix_in_db)

        return 0 if success else 1

    def _diagnose_batch(self, args: argparse.Namespace) -> int:
        """Diagnose domain partition issues for an entire batch"""
        from ecod.cli.diagnostics import DomainPartitionDiagnostics

        diagnostics = DomainPartitionDiagnostics(self.context)
        success = diagnostics.diagnose_batch(args.batch_id, args.sample_size)

        return 0 if success else 1

    def _diagnose_readiness(self, args: argparse.Namespace) -> int:
        """Check batch readiness for domain partition"""
        from ecod.cli.diagnostics import DomainPartitionDiagnostics

        diagnostics = DomainPartitionDiagnostics(self.context)
        mode = 'detailed' if args.detailed else 'summary'
        is_ready = diagnostics.check_readiness(args.batch_id, mode)

        return 0 if is_ready else 1

    # Repair command implementations
    def _run_repair(self, args: argparse.Namespace) -> int:
        """Run repair operations based on action"""
        if args.action == 'failed':
            return self._repair_failed(args)
        elif args.action == 'missing':
            return self._repair_missing(args)
        else:
            self.logger.error(f"Unknown repair action: {args.action}")
            return 1

    def _repair_failed(self, args: argparse.Namespace) -> int:
        """Reset and retry failed processes"""
        # Validate batch ID
        batch_info = self._get_batch_info(args.batch_id)
        if not batch_info:
            self.logger.error(f"Batch {args.batch_id} not found")
            return 1

        # Reset failed processes
        summary_reset = self.pipeline.reset_failed_processes(args.batch_id, 'domain_summary_failed')
        partition_reset = self.pipeline.reset_failed_processes(args.batch_id, 'domain_partition_failed')

        total_reset = summary_reset + partition_reset

        self.logger.info(f"Reset {total_reset} failed processes:")
        self.logger.info(f"  Summary failures reset: {summary_reset}")
        self.logger.info(f"  Partition failures reset: {partition_reset}")

        if total_reset == 0:
            self.logger.info("No failed processes to reset")
            return 0

        # Re-run if requested
        if args.rerun:
            self.logger.info("Re-running partition for reset processes using service")

            # Reprocess the batch
            try:
                results = self.service.reprocess_failed(
                    batch_id=args.batch_id,
                    batch_path=batch_info["base_path"],
                    blast_only=args.blast_only
                )

                self.logger.info(f"Re-run complete: {results.success_count} succeeded, {results.failure_count} failed")
                return 0 if results.success_count > 0 else 1

            except Exception as e:
                self.logger.error(f"Error re-running processes: {str(e)}")
                import traceback
                self.logger.error(traceback.format_exc())
                return 1

        return 0

    def _repair_missing(self, args: argparse.Namespace) -> int:
        """Regenerate missing domain files"""
        # Validate batch ID
        batch_info = self._get_batch_info(args.batch_id)
        if not batch_info:
            self.logger.error(f"Batch {args.batch_id} not found")
            return 1

        # Find processes with missing domain files
        query = """
        SELECT ps.id as process_id, p.pdb_id, p.chain_id
        FROM ecod_schema.process_status ps
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        LEFT JOIN ecod_schema.process_file pf_summary ON
            ps.id = pf_summary.process_id AND
            pf_summary.file_type = %s AND
            pf_summary.file_exists = TRUE
        LEFT JOIN ecod_schema.process_file pf_partition ON
            ps.id = pf_partition.process_id AND
            pf_partition.file_type = 'domain_partition' AND
            pf_partition.file_exists = TRUE
        WHERE ps.batch_id = %s
          AND pf_summary.id IS NOT NULL
          AND pf_partition.id IS NULL
        """

        summary_type = "blast_only_summary" if args.blast_only else "domain_summary"
        rows = self.db.execute_dict_query(query, (summary_type, args.batch_id))

        if not rows:
            self.logger.info(f"No missing domain files found for batch {args.batch_id}")
            return 0

        self.logger.info(f"Found {len(rows)} proteins with missing domain files")

        # Limit number of files to repair
        if args.limit and args.limit < len(rows):
            rows = rows[:args.limit]
            self.logger.info(f"Limiting repair to {args.limit} proteins")

        # Process each missing file
        success_count = 0
        for row in rows:
            try:
                # Find summary path
                summary_path = self._find_domain_summary(
                    batch_info["base_path"], row["pdb_id"], row["chain_id"],
                    batch_info["ref_version"], args.blast_only
                )

                if not summary_path:
                    self.logger.warning(f"No summary found for {row['pdb_id']}_{row['chain_id']}")
                    continue

                # Process protein
                result = self.service.partition_protein(
                    pdb_id=row["pdb_id"],
                    chain_id=row["chain_id"],
                    summary_path=summary_path,
                    output_dir=batch_info["base_path"],
                    process_id=row["process_id"],
                    blast_only=args.blast_only
                )

                if result.success:
                    success_count += 1
                    self.logger.info(f"Regenerated domain file for {row['pdb_id']}_{row['chain_id']}")
                else:
                    self.logger.error(f"Failed to regenerate for {row['pdb_id']}_{row['chain_id']}: {result.error}")

            except Exception as e:
                self.logger.error(f"Error processing {row['pdb_id']}_{row['chain_id']}: {e}")

        self.logger.info(f"Repair complete: {success_count}/{len(rows)} files regenerated")
        return 0 if success_count > 0 else 1

    # Monitor command implementation
    def _run_monitor(self, args: argparse.Namespace) -> int:
        """Monitor batch processing status"""
        # Validate batch ID
        batch_info = self._get_batch_info(args.batch_id)
        if not batch_info:
            self.logger.error(f"Batch {args.batch_id} not found")
            return 1

        self.logger.info(f"Monitoring batch {args.batch_id} ({batch_info.get('batch_name', '')}) " +
                       f"every {args.interval} seconds")

        # Monitor until completion or timeout
        start_time = time.time()
        last_progress = None

        while True:
            # Get current progress using service's tracker
            current_progress = self.service.tracker.get_batch_progress(args.batch_id)

            # Display progress if changed
            if current_progress != last_progress:
                total = current_progress.get("total", 0)
                complete = current_progress.get("complete", 0)
                errors = current_progress.get("errors", 0)
                processing = current_progress.get("processing", 0)

                self.logger.info(f"Progress: {complete}/{total} complete ({current_progress.get('complete_pct', 0):.1f}%), " +
                               f"{processing} processing, {errors} errors")

                last_progress = current_progress

                # Check if complete
                if complete + errors >= total and total > 0:
                    self.logger.info(f"Batch {args.batch_id} processing complete")
                    return 0

            # Check timeout
            if args.timeout and time.time() - start_time > args.timeout:
                self.logger.warning(f"Monitoring timeout reached after {args.timeout} seconds")
                return 0

            # Wait for next check
            time.sleep(args.interval)

    # Utility methods
    def _get_batch_info(self, batch_id: int) -> Optional[Dict[str, Any]]:
        """Get batch information from database"""
        query = """
        SELECT id, batch_name, base_path, ref_version, status
        FROM ecod_schema.batch
        WHERE id = %s
        """

        rows = self.db.execute_dict_query(query, (batch_id,))
        return rows[0] if rows else None

    def _get_all_batch_ids(self) -> List[int]:
        """Get all batch IDs from database"""
        query = "SELECT id FROM ecod_schema.batch ORDER BY id"
        rows = self.db.execute_query(query)
        return [row[0] for row in rows]

    def _get_batch_id_for_process(self, process_id: int) -> Optional[int]:
        """Get batch ID for a process"""
        query = "SELECT batch_id FROM ecod_schema.process_status WHERE id = %s"
        rows = self.db.execute_query(query, (process_id,))
        return rows[0][0] if rows else None

    def _get_process_id_for_protein(self, pdb_id: str, chain_id: str, batch_id: int) -> Optional[int]:
        """Get process ID for a specific protein"""
        query = """
        SELECT ps.id
        FROM ecod_schema.process_status ps
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        WHERE p.pdb_id = %s AND p.chain_id = %s AND ps.batch_id = %s
        """

        try:
            rows = self.db.execute_query(query, (pdb_id, chain_id, batch_id))
            return rows[0][0] if rows else None
        except Exception as e:
            self.logger.error(f"Error getting process ID: {e}")
            return None

    def _get_protein_for_process(self, process_id: int) -> Optional[Dict[str, Any]]:
        """Get protein information for a process ID"""
        query = """
        SELECT p.pdb_id, p.chain_id, ps.id as process_id
        FROM ecod_schema.process_status ps
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        WHERE ps.id = %s
        """

        rows = self.db.execute_dict_query(query, (process_id,))
        return rows[0] if rows else None

    def _find_domain_summary(self, batch_path: str, pdb_id: str, chain_id: str,
                           reference: str, blast_only: bool = False) -> Optional[str]:
        """Find domain summary file"""
        file_type = 'blast_only_summary' if blast_only else 'domain_summary'

        try:
            evidence_paths = get_all_evidence_paths(batch_path, pdb_id, chain_id, reference)

            if file_type in evidence_paths and evidence_paths[file_type]['exists_at']:
                return evidence_paths[file_type]['exists_at']

        except Exception as e:
            self.logger.warning(f"Error finding domain summary: {e}")

        return None

    def _verify_batch_readiness(self, batch_id: int, blast_only: bool = False,
                              reps_only: bool = False) -> bool:
        """Verify that a batch has the necessary domain summaries"""
        summary_type = "blast_only_summary" if blast_only else "domain_summary"
        query = """
        SELECT
            COUNT(*) as total,
            SUM(CASE WHEN EXISTS (
                SELECT 1 FROM ecod_schema.process_file pf
                WHERE pf.process_id = ps.id
                AND pf.file_type = %s
                AND pf.file_exists = TRUE
            ) THEN 1 ELSE 0 END) as ready_count
        FROM ecod_schema.process_status ps
        WHERE ps.batch_id = %s
        """

        params = [summary_type, batch_id]

        if reps_only:
            query += " AND ps.is_representative = TRUE"

        results = self.db.execute_dict_query(query, tuple(params))

        if not results:
            return False

        total = results[0]["total"]
        ready = results[0]["ready_count"]

        # Check if at least 90% of proteins are ready
        return total > 0 and ready / total >= 0.9

    def _verify_domain_summary(self, batch_path: str, pdb_id: str, chain_id: str,
                             reference: str, blast_only: bool = False) -> bool:
        """Verify that domain summary file exists for a protein"""
        all_paths = get_all_evidence_paths(batch_path, pdb_id, chain_id, reference)

        summary_type = "blast_only_summary" if blast_only else "domain_summary"

        return (summary_type in all_paths and all_paths[summary_type]["exists_at"])

    def _find_protein_in_database(self, pdb_id: str, chain_id: str,
                                batch_id: Optional[int] = None) -> List[Dict[str, Any]]:
        """Find a protein in the database"""
        query = """
        SELECT ps.id as process_id, ps.batch_id, ps.current_stage, ps.status,
               ps.error_message, ps.is_representative
        FROM ecod_schema.process_status ps
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        WHERE p.pdb_id = %s AND p.chain_id = %s
        """

        params = [pdb_id, chain_id]

        if batch_id is not None:
            query += " AND ps.batch_id = %s"
            params.append(batch_id)

        query += " ORDER BY ps.batch_id DESC"

        rows = self.db.execute_dict_query(query, tuple(params))

        # Enhance with file information
        for row in rows:
            process_id = row["process_id"]

            file_query = """
            SELECT file_type, file_path, file_exists, file_size
            FROM ecod_schema.process_file
            WHERE process_id = %s
            """

            files = self.db.execute_dict_query(file_query, (process_id,))
            row["files"] = {file["file_type"]: file for file in files}

        return rows

    def _check_domain_result(self, process_id: int, pdb_id: str, chain_id: str,
                           batch_info: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Check domain result file using model-based parsing"""
        # Get file path from database
        query = """
        SELECT file_path FROM ecod_schema.process_file
        WHERE process_id = %s AND file_type = 'domain_partition'
        """

        rows = self.db.execute_query(query, (process_id,))
        if not rows:
            return None

        file_path = rows[0][0]
        full_path = os.path.join(batch_info.get("base_path", ""), file_path)

        if not os.path.exists(full_path):
            return None

        # Parse domain file using model
        try:
            result = DomainPartitionResult.from_xml_file(full_path)

            return {
                "is_classified": result.is_classified,
                "is_unclassified": result.is_unclassified,
                "is_peptide": result.is_peptide,
                "coverage": result.coverage,
                "domains": result.domains,
                "sequence_length": result.sequence_length,
                "processing_time": result.processing_time
            }

        except Exception as e:
            self.logger.error(f"Error parsing domain file: {str(e)}")
            return None

    def _monitor_slurm_jobs(self, job_manager, job_ids: List[str], args: argparse.Namespace) -> int:
        """Monitor SLURM jobs until completion"""
        self.logger.info(f"Waiting for jobs to complete, checking every {args.check_interval} seconds")

        start_time = time.time()
        completed_jobs = set()
        failed_jobs = set()

        while job_ids and (not args.timeout or time.time() - start_time < args.timeout):
            time.sleep(args.check_interval)

            # Check each job
            for job_id in list(job_ids):
                status = job_manager.check_job_status(job_id)

                if status in ["COMPLETED", "COMPLETING"]:
                    self.logger.info(f"Job {job_id} completed successfully")
                    job_ids.remove(job_id)
                    completed_jobs.add(job_id)
                elif status in ["FAILED", "TIMEOUT", "CANCELLED", "NODE_FAIL"]:
                    self.logger.error(f"Job {job_id} failed with status {status}")
                    job_ids.remove(job_id)
                    failed_jobs.add(job_id)

            # Log status
            if job_ids:
                self.logger.info(f"Waiting for {len(job_ids)} remaining jobs")

        # Log final status
        if not job_ids:
            self.logger.info("All jobs completed")
        else:
            self.logger.warning(f"{len(job_ids)} jobs still running at end of monitoring period")

        self.logger.info(f"Job completion summary: {len(completed_jobs)} completed, {len(failed_jobs)} failed")

        return 0 if not failed_jobs and not job_ids else 1

# Compatibility functions for main.py
def setup_parser(parser: argparse.ArgumentParser) -> None:
    """Set up the argument parser for domain analysis commands"""
    cmd = DomainCommand()
    cmd.setup_parser(parser)

def run_command(args: argparse.Namespace) -> int:
    """Run the specified domain analysis command"""
    cmd = DomainCommand(args.config)
    return cmd.run_command(args)
