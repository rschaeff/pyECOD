# ecod/cli/hhsearch.py
"""
HHSearch-related commands for the ECOD pipeline

This module provides comprehensive HHSearch functionality including:
- Profile generation with HHblits
- HHSearch execution against ECOD database
- Result parsing and XML conversion
- Collation with BLAST evidence
- Batch processing with SLURM support
- Analysis and repair tools
"""

import argparse
import logging
import os
import sys
import glob
import json
import time
import random
import xml.etree.ElementTree as ET
from xml.dom import minidom
from datetime import datetime
from typing import Dict, Any, List, Optional, Tuple, Union
from concurrent.futures import ThreadPoolExecutor, as_completed

from ecod.cli.base_command import BaseCommand, handle_command_errors
from ecod.pipelines.hhsearch_pipeline import HHSearchPipeline
from ecod.pipelines.hhsearch.processor import HHSearchProcessor, HHRToXMLConverter
from ecod.pipelines.domain_analysis.hhresult_registrar import HHResultRegistrar
from ecod.utils.hhsearch_utils import HHRParser
from ecod.utils.path_utils import (
    get_standardized_paths,
    get_all_evidence_paths,
    get_file_db_path,
    resolve_file_path,
    find_files_with_legacy_paths,
    migrate_file_to_standard_path,
    scan_batch_directory
)
from ecod.jobs.factory import create_job_manager

logger = logging.getLogger("ecod.cli.hhsearch")

# Enhanced commands in this group
COMMANDS = {
    'generate': 'Generate HHblits profiles for sequences',
    'search': 'Run HHSearch against ECOD database',
    'process': 'Process HHR files to XML format',
    'collate': 'Combine HHSearch results with BLAST evidence',
    'analyze': 'Analyze HHSearch results and diagnose issues',
    'repair': 'Fix missing or problematic HHSearch files',
    'monitor': 'Monitor HHSearch job status',
}

class HHSearchCommand(BaseCommand):
    """Command handler for HHSearch operations"""

    def __init__(self, config_path: Optional[str] = None):
        super().__init__(config_path)
        self._pipeline = None
        self._processor = None
        self._registrar = None

    @property
    def pipeline(self) -> HHSearchPipeline:
        """Lazy initialization of HHSearch pipeline"""
        if self._pipeline is None:
            self._pipeline = HHSearchPipeline(self.context)
        return self._pipeline

    @property
    def processor(self) -> HHSearchProcessor:
        """Lazy initialization of HHSearch processor"""
        if self._processor is None:
            self._processor = HHSearchProcessor(self.context)
        return self._processor

    @property
    def registrar(self) -> HHResultRegistrar:
        """Lazy initialization of HHResult registrar"""
        if self._registrar is None:
            self._registrar = HHResultRegistrar(self.context)
        return self._registrar

    def setup_parser(self, parser: argparse.ArgumentParser) -> None:
        """Set up the argument parser for HHSearch commands"""
        subparsers = parser.add_subparsers(dest='command', help='HHSearch command')

        # Generate command (profile generation)
        self._setup_generate_parser(subparsers)

        # Search command (run HHSearch)
        self._setup_search_parser(subparsers)

        # Process command (convert HHR to XML)
        self._setup_process_parser(subparsers)

        # Collate command (combine with BLAST)
        self._setup_collate_parser(subparsers)

        # Analyze command
        self._setup_analyze_parser(subparsers)

        # Repair command
        self._setup_repair_parser(subparsers)

        # Monitor command
        self._setup_monitor_parser(subparsers)

    def _setup_generate_parser(self, subparsers):
        """Set up generate command parser"""
        generate_parser = subparsers.add_parser('generate', help=COMMANDS['generate'])
        generate_parser.add_argument('--batch-id', type=int, required=True,
                                  help='Batch ID to process')
        generate_parser.add_argument('--pdb-id', type=str,
                                  help='Process specific PDB')
        generate_parser.add_argument('--chain-id', type=str,
                                  help='Process specific chain')
        generate_parser.add_argument('--threads', type=int, default=8,
                                  help='Threads per job')
        generate_parser.add_argument('--memory', type=str, default='16G',
                                  help='Memory per job')
        generate_parser.add_argument('--force', action='store_true',
                                  help='Force regeneration of existing profiles')
        generate_parser.add_argument('--wait', action='store_true',
                                  help='Wait for jobs to complete')
        generate_parser.add_argument('--check-interval', type=int, default=60,
                                  help='Check interval in seconds when waiting')

    def _setup_search_parser(self, subparsers):
        """Set up search command parser"""
        search_parser = subparsers.add_parser('search', help=COMMANDS['search'])
        search_parser.add_argument('--batch-id', type=int, required=True,
                                help='Batch ID to process')
        search_parser.add_argument('--pdb-id', type=str,
                                help='Process specific PDB')
        search_parser.add_argument('--chain-id', type=str,
                                help='Process specific chain')
        search_parser.add_argument('--threads', type=int, default=8,
                                help='Threads per job')
        search_parser.add_argument('--memory', type=str, default='16G',
                                help='Memory per job')
        search_parser.add_argument('--force', action='store_true',
                                help='Force rerun of existing searches')
        search_parser.add_argument('--wait', action='store_true',
                                help='Wait for jobs to complete')
        search_parser.add_argument('--check-interval', type=int, default=60,
                                help='Check interval in seconds when waiting')

    def _setup_process_parser(self, subparsers):
        """Set up process command parser with multiple backends"""
        process_parser = subparsers.add_parser('process', help=COMMANDS['process'])
        process_subparsers = process_parser.add_subparsers(dest='backend', help='Processing backend')

        # Database backend
        db_parser = process_subparsers.add_parser('db', help='Use database backend')
        db_parser.add_argument('--batch-id', type=int, required=True,
                            help='Batch ID to process')
        db_parser.add_argument('--limit', type=int,
                            help='Limit the number of files to process')
        db_parser.add_argument('--force', action='store_true',
                            help='Force reprocessing of existing results')

        # Filesystem backend
        fs_parser = process_subparsers.add_parser('fs', help='Use filesystem backend')
        fs_parser.add_argument('--batch-path', type=str, required=True,
                            help='Path to batch directory')
        fs_parser.add_argument('--ref-version', type=str, default='develop291',
                            help='Reference version')
        fs_parser.add_argument('--limit', type=int,
                            help='Limit the number of files to process')
        fs_parser.add_argument('--force', action='store_true',
                            help='Force reprocessing of existing results')

        # Registration backend
        register_parser = process_subparsers.add_parser('register', help='Register HHR files')
        register_parser.add_argument('--batch-id', type=int, required=True,
                                  help='Batch ID to process')
        register_parser.add_argument('--chains', nargs='+',
                                  help='Specific chains to process (format: pdbid_chainid)')
        register_parser.add_argument('--force', action='store_true',
                                  help='Force reprocessing of existing results')

        # Batch processing
        batch_parser = process_subparsers.add_parser('batches', help='Process multiple batches')
        batch_parser.add_argument('--batch-ids', type=int, nargs='+',
                               help='Specific batch IDs to process')
        batch_parser.add_argument('--exclude-batch-ids', type=int, nargs='+', default=[],
                               help='Batch IDs to exclude')
        batch_parser.add_argument('--max-workers', type=int, default=4,
                               help='Maximum number of worker threads')
        batch_parser.add_argument('--force', action='store_true',
                               help='Force reprocessing of existing results')

    def _setup_collate_parser(self, subparsers):
        """Set up collate command parser"""
        collate_parser = subparsers.add_parser('collate', help=COMMANDS['collate'])
        collate_subparsers = collate_parser.add_subparsers(dest='action', help='Collation action')

        # Single batch collation
        single_parser = collate_subparsers.add_parser('batch', help='Collate single batch')
        single_parser.add_argument('--batch-id', type=int, required=True,
                                help='Batch ID to process')
        single_parser.add_argument('--limit', type=int,
                                help='Limit proteins to process')
        single_parser.add_argument('--force', action='store_true',
                                help='Force reprocessing')

        # All batches collation
        all_parser = collate_subparsers.add_parser('all', help='Collate all batches')
        all_parser.add_argument('--exclude-batch-ids', type=int, nargs='+', default=[],
                             help='Batch IDs to exclude')
        all_parser.add_argument('--limit-per-batch', type=int,
                             help='Limit proteins per batch')
        all_parser.add_argument('--force', action='store_true',
                             help='Force reprocessing')

        # Parallel collation with SLURM
        parallel_parser = collate_subparsers.add_parser('parallel', help='Parallel collation using SLURM')
        parallel_parser.add_argument('--batch-ids', type=int, nargs='+',
                                  help='Specific batch IDs to process')
        parallel_parser.add_argument('--exclude-batch-ids', type=int, nargs='+', default=[],
                                  help='Batch IDs to exclude')
        parallel_parser.add_argument('--threads', type=int, default=4,
                                  help='Threads per job')
        parallel_parser.add_argument('--memory', type=str, default='8G',
                                  help='Memory per job')
        parallel_parser.add_argument('--time', type=str, default='12:00:00',
                                  help='Time limit per job')
        parallel_parser.add_argument('--wait', action='store_true',
                                  help='Wait for jobs to complete')
        parallel_parser.add_argument('--check-interval', type=int, default=60,
                                  help='Check interval in seconds')
        parallel_parser.add_argument('--timeout', type=int,
                                  help='Timeout in seconds')
        parallel_parser.add_argument('--force', action='store_true',
                                  help='Force reprocessing')

    def _setup_analyze_parser(self, subparsers):
        """Set up analyze command parser"""
        analyze_parser = subparsers.add_parser('analyze', help=COMMANDS['analyze'])
        analyze_subparsers = analyze_parser.add_subparsers(dest='action', help='Analysis action')

        # Content analysis
        content_parser = analyze_subparsers.add_parser('content', help='Check content of HHR and XML files')
        content_parser.add_argument('--batch-path', type=str, required=True,
                                  help='Path to batch directory')
        content_parser.add_argument('--ref-version', type=str, default='develop291',
                                  help='Reference version')
        content_parser.add_argument('--sample-size', type=int, default=5,
                                  help='Number of files to examine')
        content_parser.add_argument('--output', type=str,
                                  help='Output JSON file for results')

        # Missing files analysis
        missing_parser = analyze_subparsers.add_parser('missing', help='Find missing files')
        missing_parser.add_argument('--batch-path', type=str, required=True,
                                 help='Path to batch directory')
        missing_parser.add_argument('--ref-version', type=str, default='develop291',
                                 help='Reference version')
        missing_parser.add_argument('--output', type=str,
                                 help='Output JSON file for results')

        # Statistics analysis
        stats_parser = analyze_subparsers.add_parser('stats', help='Generate HHSearch statistics')
        stats_parser.add_argument('--batch-id', type=int, required=True,
                               help='Batch ID to analyze')
        stats_parser.add_argument('--output', type=str,
                               help='Output JSON file for statistics')

    def _setup_repair_parser(self, subparsers):
        """Set up repair command parser"""
        repair_parser = subparsers.add_parser('repair', help=COMMANDS['repair'])
        repair_subparsers = repair_parser.add_subparsers(dest='action', help='Repair action')

        # Repair missing files
        missing_parser = repair_subparsers.add_parser('missing', help='Repair missing files')
        missing_parser.add_argument('--batch-path', type=str, required=True,
                                 help='Path to batch directory')
        missing_parser.add_argument('--ref-version', type=str, default='develop291',
                                 help='Reference version')
        missing_parser.add_argument('--summary-only', action='store_true',
                                 help='Only show summary, don\'t repair')

        # Fix domain summaries
        fix_parser = repair_subparsers.add_parser('summaries', help='Fix domain summaries')
        fix_parser.add_argument('--batch-path', type=str, required=True,
                             help='Path to batch directory')
        fix_parser.add_argument('--ref-version', type=str, default='develop291',
                             help='Reference version')
        fix_parser.add_argument('--limit', type=int,
                             help='Limit files to process')

        # Fix database paths
        db_paths_parser = repair_subparsers.add_parser('db_paths', help='Fix database file paths')
        db_paths_parser.add_argument('--batch-id', type=int, required=True,
                                  help='Batch ID to process')
        db_paths_parser.add_argument('--dry-run', action='store_true',
                                  help='Check but don\'t make changes')

        # Create empty summaries
        empty_parser = repair_subparsers.add_parser('empty', help='Create empty summaries')
        empty_parser.add_argument('--batch-id', type=int, required=True,
                               help='Batch ID to process')
        empty_parser.add_argument('--dry-run', action='store_true',
                               help='Check but don\'t make changes')

    def _setup_monitor_parser(self, subparsers):
        """Set up monitor command parser"""
        monitor_parser = subparsers.add_parser('monitor', help=COMMANDS['monitor'])
        monitor_parser.add_argument('--batch-id', type=int,
                                 help='Monitor specific batch')
        monitor_parser.add_argument('--interval', type=int, default=60,
                                 help='Check interval in seconds')
        monitor_parser.add_argument('--timeout', type=int,
                                 help='Timeout in seconds')

    @handle_command_errors
    def run_command(self, args: argparse.Namespace) -> int:
        """Run the specified HHSearch command"""
        # Set force_overwrite flag if needed
        if hasattr(args, 'force') and args.force:
            self.context.set_force_overwrite(True)
            self.logger.info("Force overwrite enabled")

        # Handle different commands
        if args.command == 'generate':
            return self._run_generate(args)
        elif args.command == 'search':
            return self._run_search(args)
        elif args.command == 'process':
            return self._run_process(args)
        elif args.command == 'collate':
            return self._run_collate(args)
        elif args.command == 'analyze':
            return self._run_analyze(args)
        elif args.command == 'repair':
            return self._run_repair(args)
        elif args.command == 'monitor':
            return self._run_monitor(args)
        else:
            self.logger.error(f"Unknown command: {args.command}")
            return 1

    # Generate command implementation
    def _run_generate(self, args: argparse.Namespace) -> int:
        """Generate HHblits profiles"""
        self.logger.info(f"Generating profiles for batch {args.batch_id}")

        if args.pdb_id and args.chain_id:
            # Process specific chain
            self.logger.info(f"Processing specific chain: {args.pdb_id}_{args.chain_id}")

            # Find protein ID
            query = """
            SELECT p.id
            FROM ecod_schema.protein p
            JOIN ecod_schema.process_status ps ON p.id = ps.protein_id
            WHERE p.pdb_id = %s AND p.chain_id = %s AND ps.batch_id = %s
            """

            rows = self.db.execute_query(query, (args.pdb_id, args.chain_id, args.batch_id))

            if not rows:
                self.logger.error(f"Protein {args.pdb_id}_{args.chain_id} not found in batch {args.batch_id}")
                return 1

            protein_ids = [rows[0][0]]
            job_ids = self.pipeline.generate_profiles_for_proteins(
                args.batch_id, protein_ids, args.threads, args.memory, args.force
            )
        else:
            # Process whole batch
            job_ids = self.pipeline.generate_profiles(
                args.batch_id, args.threads, args.memory, args.force
            )

        self.logger.info(f"Submitted {len(job_ids)} profile generation jobs")

        # Wait if requested
        if args.wait and job_ids:
            return self._wait_for_jobs(job_ids, args.check_interval)

        return 0 if job_ids else 1

    # Search command implementation
    def _run_search(self, args: argparse.Namespace) -> int:
        """Run HHSearch"""
        self.logger.info(f"Running HHSearch for batch {args.batch_id}")

        if args.pdb_id and args.chain_id:
            # Process specific chain
            self.logger.info(f"Processing specific chain: {args.pdb_id}_{args.chain_id}")

            # Find protein ID
            query = """
            SELECT p.id
            FROM ecod_schema.protein p
            JOIN ecod_schema.process_status ps ON p.id = ps.protein_id
            WHERE p.pdb_id = %s AND p.chain_id = %s AND ps.batch_id = %s
            """

            rows = self.db.execute_query(query, (args.pdb_id, args.chain_id, args.batch_id))

            if not rows:
                self.logger.error(f"Protein {args.pdb_id}_{args.chain_id} not found in batch {args.batch_id}")
                return 1

            protein_ids = [rows[0][0]]
            job_ids = self.pipeline.run_hhsearch_for_proteins(
                args.batch_id, protein_ids, args.threads, args.memory, args.force
            )
        else:
            # Process whole batch
            job_ids = self.pipeline.run_hhsearch(
                args.batch_id, args.threads, args.memory, args.force
            )

        self.logger.info(f"Submitted {len(job_ids)} HHSearch jobs")

        # Wait if requested
        if args.wait and job_ids:
            return self._wait_for_jobs(job_ids, args.check_interval)

        return 0 if job_ids else 1

    # Process command implementations
    def _run_process(self, args: argparse.Namespace) -> int:
        """Run process with appropriate backend"""
        if args.backend == 'db':
            return self._process_via_database(args)
        elif args.backend == 'fs':
            return self._process_via_filesystem(args)
        elif args.backend == 'register':
            return self._process_via_registration(args)
        elif args.backend == 'batches':
            return self._process_multiple_batches(args)
        else:
            self.logger.error(f"Unknown backend: {args.backend}")
            return 1

    def _process_via_database(self, args: argparse.Namespace) -> int:
        """Process HHSearch results using database backend"""
        # Get batch info
        batch_info = self._get_batch_info(args.batch_id)
        if not batch_info:
            return 1

        self.logger.info(f"Processing batch {args.batch_id} ({batch_info['batch_name']}) with reference {batch_info['ref_version']}")

        # Process batch
        processed_count = self.processor.process_batch(args.batch_id, args.force)

        if processed_count > 0:
            self.logger.info(f"Successfully processed {processed_count} chains")
            return 0
        else:
            self.logger.warning("No chains were processed")
            return 1

    def _process_via_filesystem(self, args: argparse.Namespace) -> int:
        """Process HHSearch results directly from filesystem"""
        batch_path = args.batch_path
        ref_version = args.ref_version

        self.logger.info(f"Processing HHSearch results from filesystem at {batch_path}")

        # Find all HHR files
        hhr_pattern = os.path.join(batch_path, "hhsearch", f"*.{ref_version}.hhr")
        hhr_files = glob.glob(hhr_pattern)

        self.logger.info(f"Found {len(hhr_files)} HHR files")

        if not hhr_files:
            self.logger.warning(f"No HHR files found matching pattern: {hhr_pattern}")
            return 1

        if args.limit:
            hhr_files = hhr_files[:args.limit]
            self.logger.info(f"Limited processing to {args.limit} files")

        # Create necessary directories
        hhsearch_dir = os.path.join(batch_path, "hhsearch")
        domains_dir = os.path.join(batch_path, "domains")

        os.makedirs(hhsearch_dir, exist_ok=True)
        os.makedirs(domains_dir, exist_ok=True)

        # Initialize parser and converter
        parser = HHRParser(self.logger)
        converter = HHRToXMLConverter(self.logger)

        # Process each HHR file
        processed_count = 0
        for i, hhr_file in enumerate(hhr_files):
            success = self._process_single_hhr_file(
                hhr_file, batch_path, ref_version, parser, converter, args.force
            )

            if success:
                processed_count += 1

            if (i + 1) % 10 == 0:
                self.logger.info(f"Processed {i + 1}/{len(hhr_files)} files")

        self.logger.info(f"Successfully processed {processed_count} chains")
        return 0 if processed_count > 0 else 1

    def _process_via_registration(self, args: argparse.Namespace) -> int:
        """Register HHSearch results using HHResultRegistrar"""
        try:
            if hasattr(args, 'chains') and args.chains:
                result = self.registrar.register_specific_chains(
                    args.batch_id, args.chains, args.force
                )
            else:
                result = self.registrar.register_batch_results(
                    args.batch_id, args.force
                )

            self.logger.info(f"Successfully registered {result} HHR files and converted them to XML")
            return 0 if result > 0 else 1

        except Exception as e:
            self.logger.error(f"Error processing batch: {e}")
            return 1

    def _process_multiple_batches(self, args: argparse.Namespace) -> int:
        """Process multiple batches in parallel"""
        self.logger.info(f"Processing multiple batches with {args.max_workers} workers")

        # Get batch IDs
        if args.batch_ids:
            batch_ids = [b_id for b_id in args.batch_ids if b_id not in args.exclude_batch_ids]
        else:
            batch_ids = self._get_all_batch_ids()
            batch_ids = [b_id for b_id in batch_ids if b_id not in args.exclude_batch_ids]

        self.logger.info(f"Will process {len(batch_ids)} batches: {batch_ids}")

        # Process batches using ThreadPoolExecutor
        results = {}
        failed_batches = []

        with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
            future_to_batch = {
                executor.submit(self._process_single_batch, batch_id, args.force): batch_id
                for batch_id in batch_ids
            }

            for future in as_completed(future_to_batch):
                batch_id = future_to_batch[future]
                try:
                    batch_id, result = future.result()
                    results[batch_id] = result

                    if result < 0:
                        failed_batches.append(batch_id)
                        self.logger.error(f"Batch {batch_id} failed processing")
                    else:
                        self.logger.info(f"Completed batch {batch_id}: registered {result} files")
                except Exception as e:
                    failed_batches.append(batch_id)
                    self.logger.error(f"Batch {batch_id} failed with exception: {e}")

        # Print summary
        total_processed = sum(count for count in results.values() if count > 0)
        self.logger.info(f"Completed processing all batches")
        self.logger.info(f"Total files processed: {total_processed}")

        if failed_batches:
            self.logger.error(f"Failed batches: {failed_batches}")
            return 1

        return 0

    # Collate command implementations
    def _run_collate(self, args: argparse.Namespace) -> int:
        """Run collate based on action"""
        if args.action == 'batch':
            return self._collate_batch(args)
        elif args.action == 'all':
            return self._collate_all_batches(args)
        elif args.action == 'parallel':
            return self._collate_parallel(args)
        else:
            self.logger.error(f"Unknown collate action: {args.action}")
            return 1

    def _collate_batch(self, args: argparse.Namespace) -> int:
        """Collate HHSearch results with BLAST for a single batch"""
        # Get batch info
        batch_info = self._get_batch_info(args.batch_id)
        if not batch_info:
            return 1

        self.logger.info(f"Collating HHSearch results with BLAST for batch {args.batch_id} ({batch_info['batch_name']})")

        # Get representative proteins with HHSearch results
        query = """
        SELECT
            p.id as protein_id, p.pdb_id, p.chain_id, ps.id as process_id
        FROM
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        LEFT JOIN
            ecod_schema.process_file pf ON ps.id = pf.process_id AND pf.file_type = 'hhr'
        WHERE
            ps.batch_id = %s
            AND ps.is_representative = TRUE
            AND pf.file_exists = TRUE
        ORDER BY
            p.pdb_id, p.chain_id
        """

        proteins = self.db.execute_dict_query(query, (args.batch_id,))
        self.logger.info(f"Found {len(proteins)} representative proteins with HHSearch results")

        if args.limit and args.limit < len(proteins):
            proteins = proteins[:args.limit]
            self.logger.info(f"Limited to {args.limit} proteins")

        # Process each protein
        success_count = 0
        for protein in proteins:
            result = self.processor._process_chain(
                protein['pdb_id'],
                protein['chain_id'],
                protein['process_id'],
                batch_info,
                batch_info['ref_version'],
                args.force
            )

            if result:
                success_count += 1
                self.logger.info(f"Successfully processed {protein['pdb_id']}_{protein['chain_id']}")
            else:
                self.logger.warning(f"Failed to process {protein['pdb_id']}_{protein['chain_id']}")

        self.logger.info(f"Successfully collated results for {success_count}/{len(proteins)} proteins")
        return 0 if success_count > 0 else 1

    def _collate_all_batches(self, args: argparse.Namespace) -> int:
        """Collate HHSearch results with BLAST for all batches"""
        self.logger.info("Starting collation across all batches")

        # Get all batch IDs
        batch_ids = self._get_all_batch_ids()
        batch_ids = [b_id for b_id in batch_ids if b_id not in args.exclude_batch_ids]

        self.logger.info(f"Processing {len(batch_ids)} batches")

        # Process each batch
        success_count = 0
        failed_batches = []

        for batch_id in batch_ids:
            batch_info = self._get_batch_info(batch_id)
            if not batch_info:
                failed_batches.append(batch_id)
                continue

            self.logger.info(f"Processing batch {batch_id} ({batch_info['batch_name']})")

            # Create args for single batch
            batch_args = argparse.Namespace(
                batch_id=batch_id,
                force=args.force,
                limit=args.limit_per_batch
            )

            # Process batch
            try:
                result = self._collate_batch(batch_args)
                if result == 0:
                    success_count += 1
                else:
                    failed_batches.append(batch_id)
            except Exception as e:
                failed_batches.append(batch_id)
                self.logger.error(f"Error processing batch {batch_id}: {str(e)}")

        # Log final summary
        self.logger.info(f"Collation complete. Processed {len(batch_ids)} batches.")
        self.logger.info(f"Successful: {success_count}, Failed: {len(failed_batches)}")

        if failed_batches:
            self.logger.warning(f"Failed batches: {failed_batches}")

        return 0 if len(failed_batches) == 0 else 1

    def _collate_parallel(self, args: argparse.Namespace) -> int:
        """Collate batches in parallel using SLURM"""
        self.logger.info(f"Parallel collating batches using SLURM")

        # Get batch IDs
        if args.batch_ids:
            batch_ids = [b_id for b_id in args.batch_ids if b_id not in args.exclude_batch_ids]
        else:
            batch_ids = self._get_all_batch_ids()
            batch_ids = [b_id for b_id in batch_ids if b_id not in args.exclude_batch_ids]

        self.logger.info(f"Will collate {len(batch_ids)} batches: {batch_ids}")

        # Get job manager
        try:
            from ecod.jobs import SlurmJobManager
            job_manager = SlurmJobManager(self.context.config_manager.config)
        except ImportError:
            self.logger.error("Failed to import SlurmJobManager")
            return 1

        # Create temporary directory for job scripts
        temp_dir = os.path.join(
            self.context.config_manager.get_path('output_dir', '/tmp'),
            f"collate_jobs_{int(time.time())}"
        )
        os.makedirs(temp_dir, exist_ok=True)

        # Submit jobs
        job_ids = []
        config_path = os.path.abspath(self.context.config_manager.config_path)

        for batch_id in batch_ids:
            # Create command
            cmd = f"ecod --config {config_path} hhsearch collate batch --batch-id {batch_id}"

            if args.force:
                cmd += " --force"

            # Create job script
            job_name = f"collate_batch_{batch_id}"
            output_dir = os.path.join(temp_dir, f"batch_{batch_id}")
            os.makedirs(output_dir, exist_ok=True)

            script_path = job_manager.create_job_script(
                commands=[cmd],
                job_name=job_name,
                output_dir=output_dir,
                threads=args.threads,
                memory=args.memory,
                time=args.time
            )

            # Submit job
            job_id = job_manager.submit_job(script_path)

            if job_id:
                job_ids.append(job_id)
                self.logger.info(f"Submitted collation job for batch {batch_id}, SLURM job ID: {job_id}")
            else:
                self.logger.error(f"Failed to submit collation job for batch {batch_id}")

        self.logger.info(f"Submitted {len(job_ids)} collation jobs to SLURM")

        # Monitor if requested
        if args.wait:
            return self._monitor_slurm_jobs(job_manager, job_ids, args.check_interval, args.timeout)

        return 0

    # Analyze command implementations
    def _run_analyze(self, args: argparse.Namespace) -> int:
        """Run analysis based on action"""
        if args.action == 'content':
            return self._analyze_content(args)
        elif args.action == 'missing':
            return self._analyze_missing(args)
        elif args.action == 'stats':
            return self._analyze_stats(args)
        else:
            self.logger.error(f"Unknown analyze action: {args.action}")
            return 1

    def _analyze_content(self, args: argparse.Namespace) -> int:
        """Check content of HHR and XML files"""
        self.logger.info(f"Analyzing HHSearch content in {args.batch_path}")

        # Scan batch directory
        files = scan_batch_directory(args.batch_path, args.ref_version)

        hhr_files = files['hhr']
        xml_files = files['hh_xml']

        self.logger.info(f"Found {len(hhr_files)} HHR files and {len(xml_files)} XML files")

        # Select sample files
        if len(hhr_files) > args.sample_size:
            sample_hhr_files = random.sample(hhr_files, args.sample_size)
        else:
            sample_hhr_files = hhr_files

        # Analyze each sample
        results = {
            "total_files": {
                "hhr": len(hhr_files),
                "xml": len(xml_files)
            },
            "hhr_xml_match": len(hhr_files) == len(xml_files),
            "samples": []
        }

        for hhr_file in sample_hhr_files:
            filename = os.path.basename(hhr_file)
            pdb_chain = filename.split('.')[0]

            # Find corresponding XML file
            xml_file = None
            for f in xml_files:
                if os.path.basename(f).startswith(pdb_chain):
                    xml_file = f
                    break

            # Count hits
            hhr_hits = self._count_hits_in_hhr(hhr_file)
            xml_hits = 0

            if xml_file:
                xml_hits, _ = self._check_xml_content(xml_file)

            # Add to results
            results["samples"].append({
                "pdb_chain": pdb_chain,
                "hhr_file": hhr_file,
                "xml_file": xml_file,
                "hhr_hits": hhr_hits,
                "xml_hits": xml_hits,
                "match": hhr_hits == xml_hits
            })

            self.logger.info(f"{pdb_chain}: {hhr_hits} HHR hits, {xml_hits} XML hits, match: {hhr_hits == xml_hits}")

        # Calculate summary statistics
        if results["samples"]:
            match_count = sum(1 for s in results["samples"] if s["match"])
            match_percentage = match_count / len(results["samples"]) * 100

            results["summary"] = {
                "match_percentage": match_percentage,
                "avg_hhr_hits": sum(s["hhr_hits"] for s in results["samples"]) / len(results["samples"]),
                "avg_xml_hits": sum(s["xml_hits"] for s in results["samples"]) / len(results["samples"])
            }

            self.logger.info(f"Analysis summary:")
            self.logger.info(f"  - Match percentage: {match_percentage:.1f}%")
            self.logger.info(f"  - Average HHR hits: {results['summary']['avg_hhr_hits']:.1f}")
            self.logger.info(f"  - Average XML hits: {results['summary']['avg_xml_hits']:.1f}")

        # Write results if output specified
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(results, f, indent=2)
            self.logger.info(f"Analysis results written to {args.output}")

        return 0

    def _analyze_missing(self, args: argparse.Namespace) -> int:
        """Find missing HHR or XML files"""
        self.logger.info(f"Finding missing files in {args.batch_path}")

        # Scan batch directory
        files = scan_batch_directory(args.batch_path, args.ref_version)

        # Extract PDB_CHAIN from filenames
        hhr_chains = set()
        xml_chains = set()

        for hhr_file in files['hhr']:
            filename = os.path.basename(hhr_file)
            pdb_chain = filename.split('.')[0]
            hhr_chains.add(pdb_chain)

        for xml_file in files['hh_xml']:
            filename = os.path.basename(xml_file)
            pdb_chain = filename.split('.')[0]
            xml_chains.add(pdb_chain)

        # Find missing files
        missing_xml = hhr_chains - xml_chains
        xml_no_hhr = xml_chains - hhr_chains

        self.logger.info(f"Found {len(missing_xml)} chains with missing XML files")
        self.logger.info(f"Found {len(xml_no_hhr)} chains with XML files but no HHR files")

        results = {
            "missing_xml": sorted(list(missing_xml)),
            "xml_no_hhr": sorted(list(xml_no_hhr))
        }

        # Write results if output specified
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(results, f, indent=2)
            self.logger.info(f"Missing files list written to {args.output}")

        return 0

    def _analyze_stats(self, args: argparse.Namespace) -> int:
        """Generate HHSearch statistics for a batch"""
        self.logger.info(f"Generating HHSearch statistics for batch {args.batch_id}")

        # Get batch info
        batch_info = self._get_batch_info(args.batch_id)
        if not batch_info:
            return 1

        # Query for HHSearch file statistics
        query = """
        SELECT
            COUNT(DISTINCT ps.id) as total_proteins,
            SUM(CASE WHEN pf_a3m.file_exists THEN 1 ELSE 0 END) as profiles_generated,
            SUM(CASE WHEN pf_hhr.file_exists THEN 1 ELSE 0 END) as hhsearch_complete,
            SUM(CASE WHEN pf_xml.file_exists THEN 1 ELSE 0 END) as xml_generated,
            SUM(CASE WHEN pf_summ.file_exists THEN 1 ELSE 0 END) as summaries_created
        FROM
            ecod_schema.process_status ps
        LEFT JOIN
            ecod_schema.process_file pf_a3m ON ps.id = pf_a3m.process_id AND pf_a3m.file_type = 'a3m'
        LEFT JOIN
            ecod_schema.process_file pf_hhr ON ps.id = pf_hhr.process_id AND pf_hhr.file_type = 'hhr'
        LEFT JOIN
            ecod_schema.process_file pf_xml ON ps.id = pf_xml.process_id AND pf_xml.file_type = 'hh_xml'
        LEFT JOIN
            ecod_schema.process_file pf_summ ON ps.id = pf_summ.process_id AND pf_summ.file_type = 'domain_summary'
        WHERE
            ps.batch_id = %s
        """

        result = self.db.execute_dict_query(query, (args.batch_id,))[0]

        stats = {
            "batch_id": args.batch_id,
            "batch_name": batch_info['batch_name'],
            "total_proteins": result['total_proteins'],
            "profiles_generated": result['profiles_generated'],
            "hhsearch_complete": result['hhsearch_complete'],
            "xml_generated": result['xml_generated'],
            "summaries_created": result['summaries_created'],
            "percentages": {
                "profiles": (result['profiles_generated'] / result['total_proteins'] * 100) if result['total_proteins'] > 0 else 0,
                "hhsearch": (result['hhsearch_complete'] / result['total_proteins'] * 100) if result['total_proteins'] > 0 else 0,
                "xml": (result['xml_generated'] / result['total_proteins'] * 100) if result['total_proteins'] > 0 else 0,
                "summaries": (result['summaries_created'] / result['total_proteins'] * 100) if result['total_proteins'] > 0 else 0
            }
        }

        # Display statistics
        self.logger.info(f"HHSearch Statistics for Batch {args.batch_id} ({batch_info['batch_name']}):")
        self.logger.info(f"  Total proteins: {stats['total_proteins']}")
        self.logger.info(f"  Profiles generated: {stats['profiles_generated']} ({stats['percentages']['profiles']:.1f}%)")
        self.logger.info(f"  HHSearch complete: {stats['hhsearch_complete']} ({stats['percentages']['hhsearch']:.1f}%)")
        self.logger.info(f"  XML generated: {stats['xml_generated']} ({stats['percentages']['xml']:.1f}%)")
        self.logger.info(f"  Summaries created: {stats['summaries_created']} ({stats['percentages']['summaries']:.1f}%)")

        # Write results if output specified
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(stats, f, indent=2)
            self.logger.info(f"Statistics written to {args.output}")

        return 0

    # Repair command implementations
    def _run_repair(self, args: argparse.Namespace) -> int:
        """Run repair based on action"""
        if args.action == 'missing':
            return self._repair_missing(args)
        elif args.action == 'summaries':
            return self._repair_summaries(args)
        elif args.action == 'db_paths':
            return self._repair_db_paths(args)
        elif args.action == 'empty':
            return self._repair_empty(args)
        else:
            self.logger.error(f"Unknown repair action: {args.action}")
            return 1

    def _repair_missing(self, args: argparse.Namespace) -> int:
        """Repair missing HHR or XML files"""
        self.logger.info(f"Repairing missing files in {args.batch_path}")

        # Find missing files
        files = scan_batch_directory(args.batch_path, args.ref_version)

        # Extract PDB_CHAIN mappings
        hhr_chains = {}
        xml_chains = set()

        for hhr_file in files['hhr']:
            filename = os.path.basename(hhr_file)
            pdb_chain = filename.split('.')[0]
            hhr_chains[pdb_chain] = hhr_file

        for xml_file in files['hh_xml']:
            filename = os.path.basename(xml_file)
            pdb_chain = filename.split('.')[0]
            xml_chains.add(pdb_chain)

        # Find chains with missing XML files
        missing_xml = set(hhr_chains.keys()) - xml_chains

        self.logger.info(f"Found {len(missing_xml)} chains with missing XML files")

        if args.summary_only:
            # Only show summary
            for pdb_chain in sorted(missing_xml)[:10]:
                self.logger.info(f"  - {pdb_chain}")

            if len(missing_xml) > 10:
                self.logger.info(f"  ... and {len(missing_xml) - 10} more")

            return 0

        # Initialize parser and converter
        parser = HHRParser(self.logger)
        converter = HHRToXMLConverter(self.logger)

        # Repair missing XML files
        success_count = 0
        for pdb_chain in missing_xml:
            hhr_file = hhr_chains[pdb_chain]

            if '_' not in pdb_chain:
                self.logger.warning(f"Invalid PDB chain format: {pdb_chain}")
                continue

            pdb_id, chain_id = pdb_chain.split('_')

            # Get standardized paths
            paths = get_standardized_paths(args.batch_path, pdb_id, chain_id, args.ref_version)

            try:
                # Parse HHR file
                hhr_data = parser.parse(hhr_file)

                if not hhr_data:
                    self.logger.warning(f"Failed to parse HHR file for {pdb_chain}")
                    continue

                # Convert to XML
                xml_string = converter.convert(hhr_data, pdb_id, chain_id, args.ref_version)

                if not xml_string:
                    self.logger.warning(f"Failed to convert HHR data to XML for {pdb_chain}")
                    continue

                # Save XML file
                if converter.save(xml_string, paths['hh_xml']):
                    success_count += 1
                    self.logger.info(f"Repaired XML file for {pdb_chain}")
                else:
                    self.logger.warning(f"Failed to save XML file for {pdb_chain}")

            except Exception as e:
                self.logger.error(f"Error repairing XML for {pdb_chain}: {str(e)}")

        self.logger.info(f"Repaired {success_count}/{len(missing_xml)} missing XML files")
        return 0

    def _repair_summaries(self, args: argparse.Namespace) -> int:
        """Fix domain summaries to include HHSearch hits"""
        self.logger.info(f"Fixing domain summaries in {args.batch_path}")

        # Find all HHSearch XML files
        xml_pattern = os.path.join(args.batch_path, "hhsearch", f"*.{args.ref_version}.hhsearch.xml")
        xml_files = glob.glob(xml_pattern)

        self.logger.info(f"Found {len(xml_files)} HHSearch XML files")

        if not xml_files:
            self.logger.warning("No HHSearch XML files found")
            return 1

        if args.limit and args.limit < len(xml_files):
            xml_files = xml_files[:args.limit]
            self.logger.info(f"Limited processing to {args.limit} files")

        # Process each XML file
        fixed_count = 0
        for xml_file in xml_files:
            filename = os.path.basename(xml_file)
            pdb_chain = filename.split('.')[0]

            success = self._fix_single_domain_summary(args.batch_path, pdb_chain, args.ref_version)

            if success:
                fixed_count += 1

            if fixed_count % 50 == 0 and fixed_count > 0:
                self.logger.info(f"Fixed {fixed_count} domain summaries so far")

        self.logger.info(f"Successfully fixed {fixed_count} domain summaries")
        return 0 if fixed_count > 0 else 1

    def _repair_db_paths(self, args: argparse.Namespace) -> int:
        """Fix database file paths"""
        self.logger.info(f"Fixing database paths for batch {args.batch_id}")

        # Get batch info
        batch_info = self._get_batch_info(args.batch_id)
        if not batch_info:
            return 1

        # Get all process files for this batch
        query = """
        SELECT
            pf.id, pf.process_id, pf.file_type, pf.file_path, pf.file_exists,
            p.pdb_id, p.chain_id
        FROM
            ecod_schema.process_file pf
        JOIN
            ecod_schema.process_status ps ON pf.process_id = ps.id
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        WHERE
            ps.batch_id = %s
            AND pf.file_type IN ('hhr', 'hh_xml', 'a3m', 'hhm')
        """

        process_files = self.db.execute_dict_query(query, (args.batch_id,))
        self.logger.info(f"Found {len(process_files)} HHSearch-related files in database")

        # Track statistics
        stats = {
            'total': len(process_files),
            'updated': 0,
            'already_standard': 0,
            'errors': 0,
            'file_missing': 0
        }

        # Process each file
        for file in process_files:
            try:
                pdb_id = file['pdb_id']
                chain_id = file['chain_id']
                file_type = file['file_type']
                current_path = file['file_path']

                # Resolve current absolute path
                current_abs_path = resolve_file_path(batch_info['base_path'], current_path)

                # Get standardized paths
                paths = get_standardized_paths(
                    batch_info['base_path'], pdb_id, chain_id,
                    batch_info['ref_version'], create_dirs=False
                )

                # Skip if file type not in standardized paths
                if file_type not in paths:
                    self.logger.warning(f"Unknown file type '{file_type}' for {pdb_id}_{chain_id}")
                    stats['errors'] += 1
                    continue

                # Get standard path
                standard_path = paths[file_type]
                standard_rel_path = get_file_db_path(batch_info['base_path'], standard_path)

                # Skip if already using standard path
                if current_path == standard_rel_path:
                    stats['already_standard'] += 1
                    continue

                # Check if file exists
                if not os.path.exists(current_abs_path):
                    if os.path.exists(standard_path):
                        self.logger.info(f"File exists at standard path: {standard_path}")
                    else:
                        self.logger.warning(f"File missing at both paths")
                        stats['file_missing'] += 1
                        continue
                else:
                    # Migrate file if needed
                    if not os.path.exists(standard_path) and not args.dry_run:
                        migrate_file_to_standard_path(current_abs_path, standard_path)

                # Update database
                if not args.dry_run:
                    self.db.update(
                        "ecod_schema.process_file",
                        {
                            "file_path": standard_rel_path,
                            "file_exists": os.path.exists(standard_path),
                            "file_size": os.path.getsize(standard_path) if os.path.exists(standard_path) else 0
                        },
                        "id = %s",
                        (file['id'],)
                    )
                    self.logger.info(f"Updated path: {current_path} -> {standard_rel_path}")
                else:
                    self.logger.info(f"Would update path: {current_path} -> {standard_rel_path}")

                stats['updated'] += 1

            except Exception as e:
                self.logger.error(f"Error processing file {file['id']}: {str(e)}")
                stats['errors'] += 1

        # Log statistics
        self.logger.info("Update Statistics:")
        self.logger.info(f"Total files: {stats['total']}")
        self.logger.info(f"Already standard: {stats['already_standard']}")
        self.logger.info(f"Updated: {stats['updated']}")
        self.logger.info(f"Missing: {stats['file_missing']}")
        self.logger.info(f"Errors: {stats['errors']}")

        return 0

    def _repair_empty(self, args: argparse.Namespace) -> int:
        """Create empty summary files for peptides or no-hits proteins"""
        self.logger.info(f"Creating empty summaries for batch {args.batch_id}")

        # This would integrate with regenerate_missing_summaries functionality
        # For now, provide a placeholder implementation
        self.logger.warning("Empty summary creation not fully implemented yet")

        return 0

    # Monitor command implementation
    def _run_monitor(self, args: argparse.Namespace) -> int:
        """Monitor HHSearch job status"""
        self.logger.info("Monitoring HHSearch job status")

        start_time = time.time()

        while True:
            if args.batch_id:
                completed, failed, running = self.context.job_manager.check_all_jobs(args.batch_id)
            else:
                completed, failed, running = self.context.job_manager.check_all_jobs()

            self.logger.info(f"Job status: {completed} completed, {failed} failed, {running} running")

            if running == 0:
                self.logger.info("All jobs completed")
                break

            # Check timeout
            if args.timeout and time.time() - start_time > args.timeout:
                self.logger.warning(f"Monitoring timeout reached after {args.timeout} seconds")
                break

            time.sleep(args.interval)

        return 0

    # Helper methods
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

    def _get_all_batch_ids(self) -> List[int]:
        """Get all batch IDs from database"""
        query = "SELECT id FROM ecod_schema.batch ORDER BY id"
        rows = self.db.execute_query(query)
        return [row[0] for row in rows]

    def _wait_for_jobs(self, job_ids: List[str], check_interval: int) -> int:
        """Wait for SLURM jobs to complete"""
        self.logger.info(f"Waiting for {len(job_ids)} jobs to complete")

        while job_ids:
            time.sleep(check_interval)

            completed_jobs = []
            for job_id in job_ids:
                status = self.context.job_manager.check_job_status(job_id)

                if status in ["COMPLETED", "COMPLETING"]:
                    self.logger.info(f"Job {job_id} completed successfully")
                    completed_jobs.append(job_id)
                elif status in ["FAILED", "TIMEOUT", "CANCELLED", "NODE_FAIL"]:
                    self.logger.error(f"Job {job_id} failed with status {status}")
                    completed_jobs.append(job_id)

            for job_id in completed_jobs:
                job_ids.remove(job_id)

            if job_ids:
                self.logger.info(f"Waiting for {len(job_ids)} remaining jobs")

        return 0

    def _monitor_slurm_jobs(self, job_manager, job_ids: List[str],
                           check_interval: int, timeout: Optional[int]) -> int:
        """Monitor SLURM jobs with timeout"""
        self.logger.info(f"Waiting for jobs to complete, checking every {check_interval} seconds")

        start_time = time.time()
        completed_jobs = set()
        failed_jobs = set()

        while job_ids:
            time.sleep(check_interval)

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

            if job_ids:
                self.logger.info(f"Waiting for {len(job_ids)} remaining jobs")

            if timeout and time.time() - start_time > timeout:
                self.logger.warning(f"Timeout reached after {timeout} seconds")
                break

        if not job_ids:
            self.logger.info("All jobs completed")
        else:
            self.logger.warning(f"{len(job_ids)} jobs still running at end of monitoring period")

        return 0 if not failed_jobs and not job_ids else 1

    def _process_single_batch(self, batch_id: int, force: bool) -> Tuple[int, int]:
        """Process a single batch for parallel processing"""
        self.logger.info(f"Processing HHSearch results for batch {batch_id}")

        try:
            result = self.registrar.register_batch_results(batch_id, force)
            self.logger.info(f"Successfully registered {result} HHR files for batch {batch_id}")
            return batch_id, result
        except Exception as e:
            self.logger.error(f"Error processing batch {batch_id}: {e}")
            return batch_id, -1

    def _process_single_hhr_file(self, hhr_file: str, batch_path: str, ref_version: str,
                                parser: HHRParser, converter: HHRToXMLConverter,
                                force: bool) -> bool:
        """Process a single HHR file"""
        # Extract PDB and chain ID
        filename = os.path.basename(hhr_file)
        parts = filename.split('.')
        pdb_chain = parts[0]
        pdb_parts = pdb_chain.split('_')

        if len(pdb_parts) != 2:
            self.logger.warning(f"Invalid filename format: {filename}")
            return False

        pdb_id = pdb_parts[0]
        chain_id = pdb_parts[1]

        # Get standardized paths
        paths = get_standardized_paths(batch_path, pdb_id, chain_id, ref_version)

        # Skip if already exists and not forcing
        if os.path.exists(paths['hh_xml']) and not force:
            self.logger.debug(f"XML already exists for {pdb_chain}, skipping")
            return False

        try:
            # Parse HHR file
            hhr_data = parser.parse(hhr_file)

            if not hhr_data:
                self.logger.warning(f"Failed to parse HHR file: {hhr_file}")
                return False

            # Convert to XML
            xml_string = converter.convert(hhr_data, pdb_id, chain_id, ref_version)

            if not xml_string:
                self.logger.warning(f"Failed to convert HHR data to XML for {pdb_chain}")
                return False

            # Save XML
            if not converter.save(xml_string, paths['hh_xml']):
                self.logger.warning(f"Failed to save XML: {paths['hh_xml']}")
                return False

            self.logger.debug(f"Successfully processed {pdb_chain}")
            return True

        except Exception as e:
            self.logger.error(f"Error processing {pdb_chain}: {str(e)}")
            return False

    def _count_hits_in_hhr(self, hhr_file: str) -> int:
        """Count the number of hits in an HHR file"""
        try:
            with open(hhr_file, 'r') as f:
                content = f.read()

            lines = content.split('\n')
            hit_count = 0

            # Find the table header
            table_start = None
            for i, line in enumerate(lines):
                if line.startswith(' No Hit'):
                    table_start = i + 1
                    break

            if not table_start:
                return 0

            # Count hit entries
            for i in range(table_start, len(lines)):
                line = lines[i].strip()
                if line and line[0].isdigit() and not line.startswith('Q ') and not line.startswith('T '):
                    hit_count += 1

            return hit_count
        except Exception as e:
            self.logger.warning(f"Error reading HHR file {hhr_file}: {str(e)}")
            return 0

    def _check_xml_content(self, xml_file: str) -> Tuple[int, List[str]]:
        """Check the content of an XML file and return hit count and info"""
        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()

            # Look for hh_hit elements
            hit_list = root.find(".//hh_hit_list")
            if hit_list is None:
                return 0, []

            hits = hit_list.findall("hh_hit")

            hit_info = []
            for hit in hits:
                hit_num = hit.get("hit_num", "unknown")
                hit_id = hit.get("hit_id", "unknown")
                probability = hit.get("probability", "unknown")
                hit_info.append(f"{hit_num}: {hit_id} (prob: {probability})")

            return len(hits), hit_info
        except Exception as e:
            self.logger.warning(f"Error parsing XML file {xml_file}: {str(e)}")
            return 0, []

    def _fix_single_domain_summary(self, batch_path: str, pdb_chain: str, ref_version: str) -> bool:
        """Fix a single domain summary to include HHSearch hits"""
        hhsearch_dir = os.path.join(batch_path, "hhsearch")
        domains_dir = os.path.join(batch_path, "domains")

        # Define file paths
        xml_path = os.path.join(hhsearch_dir, f"{pdb_chain}.{ref_version}.hhsearch.xml")
        summary_path = os.path.join(domains_dir, f"{pdb_chain}.{ref_version}.domains.xml")

        if not os.path.exists(xml_path):
            self.logger.warning(f"HHSearch XML file not found: {xml_path}")
            return False

        # Try alternative summary paths if not found
        if not os.path.exists(summary_path):
            alt_paths = [
                os.path.join(domains_dir, f"{pdb_chain}.domains.xml"),
                os.path.join(domains_dir, f"{pdb_chain}.domain_summary.xml"),
                os.path.join(domains_dir, f"{pdb_chain}.{ref_version}.domain_summary.xml")
            ]

            for path in alt_paths:
                if os.path.exists(path):
                    summary_path = path
                    break

            if not os.path.exists(summary_path):
                self.logger.warning(f"Domain summary file not found: {summary_path}")
                return False

        try:
            # Parse XML files
            hh_tree = ET.parse(xml_path)
            hh_root = hh_tree.getroot()

            summary_tree = ET.parse(summary_path)
            summary_root = summary_tree.getroot()

            # Find or create hhsearch_evidence section
            hhsearch_elem = summary_root.find(".//hhsearch_evidence")
            if hhsearch_elem is None:
                hhsearch_elem = ET.SubElement(summary_root, "hhsearch_evidence")
            else:
                # Clear existing content
                for child in list(hhsearch_elem):
                    hhsearch_elem.remove(child)

            # Copy hit_list from HHSearch XML
            hit_list = hh_root.find(".//hh_hit_list")
            if hit_list is None:
                self.logger.warning(f"No hit_list found in HHSearch XML for {pdb_chain}")
                return False

            # Create new hit_list in summary
            new_hit_list = ET.SubElement(hhsearch_elem, "hh_hit_list")

            # Copy all hits
            hit_count = 0
            for hit in hit_list.findall("hh_hit"):
                hit_string = ET.tostring(hit)
                new_hit = ET.fromstring(hit_string)
                new_hit_list.append(new_hit)
                hit_count += 1

            self.logger.info(f"Copied {hit_count} hits to domain summary for {pdb_chain}")

            # Save updated summary
            rough_string = ET.tostring(summary_root, 'utf-8')
            reparsed = minidom.parseString(rough_string)
            pretty_xml = reparsed.toprettyxml(indent="  ")

            with open(summary_path, 'w', encoding='utf-8') as f:
                f.write(pretty_xml)

            return True

        except Exception as e:
            self.logger.error(f"Error fixing domain summary for {pdb_chain}: {str(e)}")
            return False

# Compatibility functions for main.py
def setup_parser(parser: argparse.ArgumentParser) -> None:
    """Set up the argument parser for HHSearch commands"""
    cmd = HHSearchCommand()
    cmd.setup_parser(parser)

def run_command(args: argparse.Namespace) -> int:
    """Run the specified HHSearch command"""
    cmd = HHSearchCommand(args.config)
    return cmd.run_command(args)
