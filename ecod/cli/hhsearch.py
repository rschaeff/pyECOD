# ecod/cli/hhsearch.py
"""
HHSearch-related commands for the ECOD pipeline

This module provides commands for running and managing HHSearch analyses,
including profile generation, search execution, result parsing, and analysis.
"""

import argparse
import logging
import os
import glob
import re
import json
import xml.etree.ElementTree as ET
from typing import Dict, Any, List, Optional, Tuple
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

from ecod.cli.base_command import BaseCommand, handle_command_errors
from ecod.pipelines.hhsearch_pipeline import HHSearchPipeline
from ecod.pipelines.hhsearch.processor import HHSearchProcessor, HHRToXMLConverter
from ecod.pipelines.domain_analysis.hhresult_registrar import HHResultRegistrar
from ecod.utils.hhsearch_utils import HHRParser
from ecod.utils.path_utils import (
    get_standardized_paths,
    get_file_db_path,
    resolve_file_path,
    scan_batch_directory
)

logger = logging.getLogger("ecod.cli.hhsearch")

# Enhanced commands in this group
COMMANDS = {
    'generate': 'Generate HHblits profiles for sequences',
    'search': 'Run HHSearch against ECOD database',
    'parse': 'Parse HHSearch results',
    'check': 'Check the status of HHSearch jobs',
    'process': 'Process HHR files to XML format',
    'analyze': 'Analyze HHSearch results',
    'collate': 'Combine HHSearch results with BLAST evidence',
    'repair': 'Fix missing or problematic HHSearch files',
    'batches': 'Process multiple batches in parallel'
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
        generate_parser.add_argument('--force', action='store_true',
                                  help='Force regeneration of existing profiles')
        generate_parser.add_argument('--pdb-id', type=str,
                                  help='Process specific PDB')
        generate_parser.add_argument('--chain-id', type=str,
                                  help='Process specific chain')

        # Search command
        search_parser = subparsers.add_parser('search', help=COMMANDS['search'])
        search_parser.add_argument('--batch-id', type=int, required=True,
                                help='Batch ID to process')
        search_parser.add_argument('--threads', type=int, default=8,
                                help='Threads per job')
        search_parser.add_argument('--memory', type=str, default='16G',
                                help='Memory per job')
        search_parser.add_argument('--force', action='store_true',
                                help='Force rerun of existing searches')
        search_parser.add_argument('--pdb-id', type=str,
                                help='Process specific PDB')
        search_parser.add_argument('--chain-id', type=str,
                                help='Process specific chain')

        # Parse command
        parse_parser = subparsers.add_parser('parse', help=COMMANDS['parse'])
        parse_parser.add_argument('--batch-id', type=int, required=True,
                               help='Batch ID to process')
        parse_parser.add_argument('--hhr-file', type=str,
                               help='Process specific HHR file')
        parse_parser.add_argument('--force', action='store_true',
                               help='Force reprocessing of existing results')

        # Check command
        check_parser = subparsers.add_parser('check', help=COMMANDS['check'])
        check_parser.add_argument('--batch-id', type=int,
                               help='Check specific batch')
        check_parser.add_argument('--interval', type=int, default=300,
                               help='Check interval in seconds')
        check_parser.add_argument('--continuous', action='store_true',
                               help='Check continuously until completion')

        # Process command
        process_parser = subparsers.add_parser('process', help=COMMANDS['process'])
        process_parser.add_argument('--batch-id', type=int,
                                 help='Process batch by ID')
        process_parser.add_argument('--batch-path', type=str,
                                 help='Process batch by path')
        process_parser.add_argument('--ref-version', type=str, default='develop291',
                                 help='Reference version')
        process_parser.add_argument('--limit', type=int,
                                 help='Limit the number of files to process')
        process_parser.add_argument('--force', action='store_true',
                                 help='Force reprocessing of existing results')

        # Analyze command
        analyze_parser = subparsers.add_parser('analyze', help=COMMANDS['analyze'])
        analyze_subparsers = analyze_parser.add_subparsers(dest='analyze_action', help='Analysis action')

        # Content analysis
        content_parser = analyze_subparsers.add_parser('content', help='Check content of HHR and XML files')
        content_parser.add_argument('--batch-path', type=str, required=True,
                                  help='Path to batch directory')
        content_parser.add_argument('--ref-version', type=str, default='develop291',
                                  help='Reference version')
        content_parser.add_argument('--sample-size', type=int, default=5,
                                  help='Number of files to examine')
        content_parser.add_argument('--output', type=str,
                                  help='Output JSON file for analysis results')

        # Missing files analysis
        missing_parser = analyze_subparsers.add_parser('missing', help='Find missing files')
        missing_parser.add_argument('--batch-path', type=str, required=True,
                                 help='Path to batch directory')
        missing_parser.add_argument('--ref-version', type=str, default='develop291',
                                 help='Reference version')
        missing_parser.add_argument('--output', type=str,
                                 help='Output JSON file for missing files list')

        # Collate command
        collate_parser = subparsers.add_parser('collate', help=COMMANDS['collate'])
        collate_parser.add_argument('--batch-id', type=int, required=True,
                                 help='Batch ID to process')
        collate_parser.add_argument('--limit', type=int,
                                 help='Limit the number of proteins to process')
        collate_parser.add_argument('--force', action='store_true',
                                 help='Force reprocessing of existing collations')

        # Repair command
        repair_parser = subparsers.add_parser('repair', help=COMMANDS['repair'])
        repair_subparsers = repair_parser.add_subparsers(dest='repair_action', help='Repair action')

        # Missing files repair
        missing_repair_parser = repair_subparsers.add_parser('missing', help='Repair missing files')
        missing_repair_parser.add_argument('--batch-path', type=str, required=True,
                                        help='Path to batch directory')
        missing_repair_parser.add_argument('--ref-version', type=str, default='develop291',
                                        help='Reference version')
        missing_repair_parser.add_argument('--summary-only', action='store_true',
                                        help='Only show summary, don\'t repair')

        # Fix database paths
        db_paths_parser = repair_subparsers.add_parser('db_paths', help='Fix database file paths')
        db_paths_parser.add_argument('--batch-id', type=int, required=True,
                                  help='Batch ID to process')
        db_paths_parser.add_argument('--dry-run', action='store_true',
                                  help='Check but don\'t make changes')

        # Batches command
        batches_parser = subparsers.add_parser('batches', help=COMMANDS['batches'])
        batches_parser.add_argument('--batch-ids', type=int, nargs='+',
                                 help='Specific batch IDs to process')
        batches_parser.add_argument('--exclude-batch-ids', type=int, nargs='+', default=[],
                                 help='Batch IDs to exclude')
        batches_parser.add_argument('--max-workers', type=int, default=4,
                                 help='Maximum number of worker threads')
        batches_parser.add_argument('--force', action='store_true',
                                 help='Force reprocessing of existing results')
        batches_parser.add_argument('--output', type=str,
                                 help='Output JSON file for detailed results')

    @handle_command_errors
    def run_command(self, args: argparse.Namespace) -> int:
        """Run the specified HHSearch command"""
        # Handle different commands
        if args.command == 'generate':
            return self._generate_profiles(args)
        elif args.command == 'search':
            return self._run_search(args)
        elif args.command == 'parse':
            return self._parse_results(args)
        elif args.command == 'check':
            return self._check_status(args)
        elif args.command == 'process':
            return self._process_hhr(args)
        elif args.command == 'analyze':
            return self._analyze_results(args)
        elif args.command == 'collate':
            return self._collate_results(args)
        elif args.command == 'repair':
            return self._repair_files(args)
        elif args.command == 'batches':
            return self._process_batches(args)
        else:
            self.logger.error(f"Unknown command: {args.command}")
            return 1

    def _generate_profiles(self, args: argparse.Namespace) -> int:
        """Generate HHblits profiles"""
        # Initialize HHSearch pipeline
        hhsearch_pipeline = HHSearchPipeline(self.context)

        self.logger.info(f"Generating profiles for batch {args.batch_id}")

        if args.pdb_id and args.chain_id:
            # Process specific chain
            self.logger.info(f"Processing specific chain: {args.pdb_id}_{args.chain_id}")

            # Find protein ID for specified PDB and chain
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
            job_ids = hhsearch_pipeline.generate_profiles_for_proteins(
                args.batch_id, protein_ids, args.threads, args.memory, args.force
            )
        else:
            # Process whole batch
            job_ids = hhsearch_pipeline.generate_profiles(
                args.batch_id, args.threads, args.memory, args.force
            )

        self.logger.info(f"Submitted {len(job_ids)} profile generation jobs")
        return 0 if job_ids else 1

    def _run_search(self, args: argparse.Namespace) -> int:
        """Run HHSearch"""
        # Initialize HHSearch pipeline
        hhsearch_pipeline = HHSearchPipeline(self.context)

        self.logger.info(f"Running HHSearch for batch {args.batch_id}")

        if args.pdb_id and args.chain_id:
            # Process specific chain
            self.logger.info(f"Processing specific chain: {args.pdb_id}_{args.chain_id}")

            # Find protein ID for specified PDB and chain
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
            job_ids = hhsearch_pipeline.run_hhsearch_for_proteins(
                args.batch_id, protein_ids, args.threads, args.memory, args.force
            )
        else:
            # Process whole batch
            job_ids = hhsearch_pipeline.run_hhsearch(
                args.batch_id, args.threads, args.memory, args.force
            )

        self.logger.info(f"Submitted {len(job_ids)} HHSearch jobs")
        return 0 if job_ids else 1

    def _parse_results(self, args: argparse.Namespace) -> int:
        """Parse HHSearch results"""
        # Initialize HHResult registrar
        registrar = HHResultRegistrar(self.context)

        if args.hhr_file:
            # Parse specific HHR file
            self.logger.info(f"Parsing specific HHR file: {args.hhr_file}")

            if not os.path.exists(args.hhr_file):
                self.logger.error(f"HHR file not found: {args.hhr_file}")
                return 1

            # Extract PDB and chain ID from filename
            filename = os.path.basename(args.hhr_file)
            match = re.match(r'([^_]+)_([^.]+)', filename)

            if not match:
                self.logger.error(f"Could not extract PDB and chain ID from filename: {filename}")
                return 1

            pdb_id, chain_id = match.groups()

            # Find process ID
            query = """
            SELECT ps.id
            FROM ecod_schema.process_status ps
            JOIN ecod_schema.protein p ON ps.protein_id = p.id
            WHERE p.pdb_id = %s AND p.chain_id = %s AND ps.batch_id = %s
            """

            rows = self.db.execute_query(query, (pdb_id, chain_id, args.batch_id))

            if not rows:
                self.logger.error(f"Process for {pdb_id}_{chain_id} not found in batch {args.batch_id}")
                return 1

            process_id = rows[0][0]

            # Parse and register HHR file
            result = registrar.register_hhr_file(process_id, args.hhr_file, args.force)

            if result:
                self.logger.info(f"Successfully parsed and registered HHR file: {args.hhr_file}")
                return 0
            else:
                self.logger.error(f"Failed to parse and register HHR file: {args.hhr_file}")
                return 1
        else:
            # Parse all HHR files for batch
            self.logger.info(f"Parsing all HHR files for batch {args.batch_id}")

            result = registrar.register_batch_results(args.batch_id, args.force)

            self.logger.info(f"Successfully registered {result} HHR files for batch {args.batch_id}")
            return 0 if result > 0 else 1

    def _check_status(self, args: argparse.Namespace) -> int:
        """Check HHSearch job status"""
        # Initialize HHSearch pipeline
        hhsearch_pipeline = HHSearchPipeline(self.context)

        self.logger.info("Checking HHSearch job status")

        if args.continuous:
            # Check continuously until all jobs are done
            self.logger.info(f"Continuous checking every {args.interval} seconds")

            import time
            while True:
                completed, failed, running = self.context.job_manager.check_all_jobs(args.batch_id)

                self.logger.info(f"Job status: {completed} completed, {failed} failed, {running} running")

                if running == 0:
                    break

                time.sleep(args.interval)
        else:
            # Check once
            completed, failed, running = self.context.job_manager.check_all_jobs(args.batch_id)

            self.logger.info(f"Job status: {completed} completed, {failed} failed, {running} running")

        return 0

    def _process_hhr(self, args: argparse.Namespace) -> int:
        """Process HHR files to XML format"""
        if args.batch_id:
            # Process batch by ID
            # Initialize HHResult registrar
            registrar = HHResultRegistrar(self.context)

            self.logger.info(f"Processing HHR files for batch {args.batch_id}")

            result = registrar.register_batch_results(args.batch_id, args.force)

            self.logger.info(f"Successfully processed {result} HHR files for batch {args.batch_id}")
            return 0 if result > 0 else 1
        elif args.batch_path:
            # Process batch by path
            # Initialize processor
            processor = HHSearchProcessor(self.context)

            self.logger.info(f"Processing HHR files in directory: {args.batch_path}")

            # Find all HHR files
            hhr_pattern = os.path.join(args.batch_path, "hhsearch", f"*.{args.ref_version}.hhr")
            hhr_files = glob.glob(hhr_pattern)

            if not hhr_files:
                self.logger.warning(f"No HHR files found matching pattern: {hhr_pattern}")
                return 1

            self.logger.info(f"Found {len(hhr_files)} HHR files")

            if args.limit:
                hhr_files = hhr_files[:args.limit]
                self.logger.info(f"Limited processing to {args.limit} files")

            # Initialize parser and converter
            parser = HHRParser(self.logger)
            converter = HHRToXMLConverter(self.logger)

            # Process each HHR file
            success_count = 0
            for hhr_file in hhr_files:
                # Extract PDB and chain ID from filename
                filename = os.path.basename(hhr_file)
                match = re.match(r'([^_]+)_([^.]+)', filename)

                if not match:
                    self.logger.warning(f"Could not extract PDB and chain ID from filename: {filename}")
                    continue

                pdb_id, chain_id = match.groups()

                # Get standardized file paths
                paths = get_standardized_paths(args.batch_path, pdb_id, chain_id, args.ref_version)

                # Skip if XML already exists and not forcing
                if os.path.exists(paths['hh_xml']) and not args.force:
                    self.logger.debug(f"XML file already exists for {pdb_id}_{chain_id}, skipping")
                    continue

                try:
                    # Parse HHR file
                    hhr_data = parser.parse(hhr_file)

                    if not hhr_data:
                        self.logger.warning(f"Failed to parse HHR file: {hhr_file}")
                        continue

                    # Convert to XML
                    xml_string = converter.convert(hhr_data, pdb_id, chain_id, args.ref_version)

                    if not xml_string:
                        self.logger.warning(f"Failed to convert HHR data to XML for {pdb_id}_{chain_id}")
                        continue

                    # Save XML file
                    if converter.save(xml_string, paths['hh_xml']):
                        success_count += 1
                        self.logger.info(f"Successfully processed HHR file for {pdb_id}_{chain_id}")
                    else:
                        self.logger.warning(f"Failed to save XML file for {pdb_id}_{chain_id}")

                except Exception as e:
                    self.logger.error(f"Error processing HHR file {hhr_file}: {str(e)}")

            self.logger.info(f"Successfully processed {success_count}/{len(hhr_files)} HHR files")
            return 0 if success_count > 0 else 1
        else:
            self.logger.error("Either --batch-id or --batch-path must be specified")
            return 1

    def _analyze_results(self, args: argparse.Namespace) -> int:
        """Analyze HHSearch results"""
        if args.analyze_action == 'content':
            return self._analyze_content(args)
        elif args.analyze_action == 'missing':
            return self._find_missing_files(args)
        else:
            self.logger.error(f"Unknown analyze action: {args.analyze_action}")
            return 1

    def _analyze_content(self, args: argparse.Namespace) -> int:
        """Check content of HHR and XML files"""
        self.logger.info(f"Analyzing HHSearch content in {args.batch_path}")

        # Scan batch directory
        files = scan_batch_directory(args.batch_path, args.ref_version)

        hhr_files = files['hhr']
        xml_files = files['hh_xml']

        self.logger.info(f"Found {len(hhr_files)} HHR files and {len(xml_files)} XML files")

        # Select sample files for detailed analysis
        if len(hhr_files) > args.sample_size:
            import random
            sample_hhr_files = random.sample(hhr_files, args.sample_size)
        else:
            sample_hhr_files = hhr_files

        # Analyze each sample file
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
            pdb_chain = filename.split('.')[0]  # Format: pdbid_chain

            # Find corresponding XML file
            xml_file = None
            for f in xml_files:
                if os.path.basename(f).startswith(pdb_chain):
                    xml_file = f
                    break

            # Count hits in HHR file
            parser = HHRParser(self.logger)
            hhr_data = parser.parse(hhr_file)
            hhr_hits = len(hhr_data.get("hits", [])) if hhr_data else 0

            # Count hits in XML file
            xml_hits = 0
            if xml_file:
                try:
                    tree = ET.parse(xml_file)
                    root = tree.getroot()
                    hits = root.findall(".//hh_hit")
                    xml_hits = len(hits)
                except Exception as e:
                    self.logger.warning(f"Error parsing XML file {xml_file}: {str(e)}")

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
        match_count = sum(1 for s in results["samples"] if s["match"])
        match_percentage = match_count / len(results["samples"]) if results["samples"] else 0

        results["summary"] = {
            "match_percentage": match_percentage * 100,
            "avg_hhr_hits": sum(s["hhr_hits"] for s in results["samples"]) / len(results["samples"]) if results["samples"] else 0,
            "avg_xml_hits": sum(s["xml_hits"] for s in results["samples"]) / len(results["samples"]) if results["samples"] else 0
        }

        self.logger.info(f"Analysis summary:")
        self.logger.info(f"  - Match percentage: {results['summary']['match_percentage']:.1f}%")
        self.logger.info(f"  - Average HHR hits: {results['summary']['avg_hhr_hits']:.1f}")
        self.logger.info(f"  - Average XML hits: {results['summary']['avg_xml_hits']:.1f}")

        # Write results to file if specified
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(results, f, indent=2)
            self.logger.info(f"Analysis results written to {args.output}")

        return 0

    def _find_missing_files(self, args: argparse.Namespace) -> int:
        """Find missing HHR or XML files"""
        self.logger.info(f"Finding missing files in {args.batch_path}")

        # Scan batch directory
        files = scan_batch_directory(args.batch_path, args.ref_version)

        hhr_files = files['hhr']
        xml_files = files['hh_xml']

        self.logger.info(f"Found {len(hhr_files)} HHR files and {len(xml_files)} XML files")

        # Extract PDB_CHAIN from filenames
        hhr_chains = set()
        xml_chains = set()

        for hhr_file in hhr_files:
            filename = os.path.basename(hhr_file)
            parts = filename.split('.')
            pdb_chain = parts[0]  # Format: pdbid_chain
            hhr_chains.add(pdb_chain)

        for xml_file in xml_files:
            filename = os.path.basename(xml_file)
            parts = filename.split('.')
            pdb_chain = parts[0]  # Format: pdbid_chain
            xml_chains.add(pdb_chain)

        # Find chains with missing files
        missing_xml = hhr_chains - xml_chains
        xml_no_hhr = xml_chains - hhr_chains

        self.logger.info(f"Found {len(missing_xml)} chains with missing XML files")
        self.logger.info(f"Found {len(xml_no_hhr)} chains with XML files but no HHR files")

        results = {
            "missing_xml": sorted(list(missing_xml)),
            "xml_no_hhr": sorted(list(xml_no_hhr))
        }

        # Write results to file if specified
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(results, f, indent=2)
            self.logger.info(f"Missing files list written to {args.output}")

        return 0

    def _collate_results(self, args: argparse.Namespace) -> int:
        """Collate HHSearch results with BLAST evidence"""
        # Initialize HHSearch processor
        processor = HHSearchProcessor(self.context)

        self.logger.info(f"Collating HHSearch results for batch {args.batch_id}")

        # Get batch info
        batch_info = self._get_batch_info(args.batch_id)
        if not batch_info:
            return 1

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
            pdb_id = protein['pdb_id']
            chain_id = protein['chain_id']
            process_id = protein['process_id']

            self.logger.info(f"Processing {pdb_id}_{chain_id}")

            # Process the chain
            result = processor._process_chain(
                pdb_id,
                chain_id,
                process_id,
                batch_info,
                batch_info['ref_version'],
                args.force
            )

            if result:
                success_count += 1
                self.logger.info(f"Successfully processed {pdb_id}_{chain_id}")
            else:
                self.logger.warning(f"Failed to process {pdb_id}_{chain_id}")

        self.logger.info(f"Successfully collated results for {success_count}/{len(proteins)} proteins")
        return 0 if success_count > 0 else 1

    def _repair_files(self, args: argparse.Namespace) -> int:
        """Repair missing or problematic HHSearch files"""
        if args.repair_action == 'missing':
            return self._repair_missing_files(args)
        elif args.repair_action == 'db_paths':
            return self._repair_db_paths(args)
        else:
            self.logger.error(f"Unknown repair action: {args.repair_action}")
            return 1

    def _repair_missing_files(self, args: argparse.Namespace) -> int:
        """Repair missing HHR or XML files"""
        self.logger.info(f"Repairing missing files in {args.batch_path}")

        # First find missing files
        files = scan_batch_directory(args.batch_path, args.ref_version)

        hhr_files = files['hhr']
        xml_files = files['hh_xml']

        # Extract PDB_CHAIN from filenames
        hhr_chains = {}
        xml_chains = set()

        for hhr_file in hhr_files:
            filename = os.path.basename(hhr_file)
            parts = filename.split('.')
            pdb_chain = parts[0]  # Format: pdbid_chain
            hhr_chains[pdb_chain] = hhr_file

        for xml_file in xml_files:
            filename = os.path.basename(xml_file)
            parts = filename.split('.')
            pdb_chain = parts[0]  # Format: pdbid_chain
            xml_chains.add(pdb_chain)

        # Find chains with missing XML files
        missing_xml = set(hhr_chains.keys()) - xml_chains

        self.logger.info(f"Found {len(missing_xml)} chains with missing XML files")

        if args.summary_only:
            # Only show summary
            for pdb_chain in sorted(missing_xml)[:10]:  # Show first 10
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

    def _repair_db_paths(self, args: argparse.Namespace) -> int:
        """Fix database file paths"""
        self.logger.info(f"Fixing database paths for batch {args.batch_id}")

        # Get batch info
        batch_info = self._get_batch_info(args.batch_id)
        if not batch_info:
            return 1

        batch_path = batch_info['base_path']
        ref_version = batch_info['ref_version']

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
                if os.path.isabs(current_path):
                    current_abs_path = current_path
                else:
                    current_abs_path = os.path.join(batch_path, current_path)

                # Get standardized paths
                paths = get_standardized_paths(batch_path, pdb_id, chain_id, ref_version)

                # Skip if file type not in standardized paths
                if file_type not in paths:
                    self.logger.warning(f"Unknown file type '{file_type}' for {pdb_id}_{chain_id}")
                    stats['errors'] += 1
                    continue

                # Get standard path
                standard_path = paths[file_type]
                standard_rel_path = os.path.relpath(standard_path, batch_path)

                # Skip if already using standard path
                if current_path == standard_rel_path:
                    self.logger.debug(f"Already using standard path: {current_path}")
                    stats['already_standard'] += 1
                    continue

                # Check if file exists at current path
                if not os.path.exists(current_abs_path):
                    self.logger.warning(f"File does not exist at current path: {current_abs_path}")

                    # Check if file exists at standard path
                    if os.path.exists(standard_path):
                        self.logger.info(f"File exists at standard path: {standard_path}")
                    else:
                        self.logger.warning(f"File missing at both current and standard paths")
                        stats['file_missing'] += 1
                        continue
                else:
                    # Migrate file to standard path if needed
                    if not os.path.exists(standard_path) and not args.dry_run:
                        os.makedirs(os.path.dirname(standard_path), exist_ok=True)
                        import shutil
                        shutil.copy2(current_abs_path, standard_path)
                        self.logger.info(f"Copied file to standard path: {standard_path}")

                # Update database path
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
                    self.logger.info(f"Updated path in database: {current_path} -> {standard_rel_path}")
                else:
                    self.logger.info(f"Would update path: {current_path} -> {standard_rel_path}")

                stats['updated'] += 1

            except Exception as e:
                self.logger.error(f"Error processing file {file['id']}: {str(e)}")
                stats['errors'] += 1

        # Log statistics
        self.logger.info("Update Statistics:")
        self.logger.info(f"Total files processed: {stats['total']}")
        self.logger.info(f"Files already using standard paths: {stats['already_standard']}")
        self.logger.info(f"Files updated: {stats['updated']}")
        self.logger.info(f"Files missing: {stats['file_missing']}")
        self.logger.info(f"Errors: {stats['errors']}")

        return 0

    def _process_batches(self, args: argparse.Namespace) -> int:
        """Process multiple batches in parallel"""
        self.logger.info(f"Processing multiple batches with {args.max_workers} workers")

        # Get batch IDs to process
        if args.batch_ids:
            # Use specified batch IDs
            batch_ids = args.batch_ids
            self.logger.info(f"Using specified batch IDs: {batch_ids}")
        else:
            # Get all batch IDs from database
            query = """
            SELECT id, batch_name
            FROM ecod_schema.batch
            ORDER BY id
            """

            rows = self.db.execute_dict_query(query)

            if not rows:
                self.logger.error("No batches found in database")
                return 1

            batch_ids = [row['id'] for row in rows]
            self.logger.info(f"Found {len(batch_ids)} batches in database")

        # Filter out excluded batch IDs
        if args.exclude_batch_ids:
            batch_ids = [b_id for b_id in batch_ids if b_id not in args.exclude_batch_ids]
            self.logger.info(f"After exclusions: {len(batch_ids)} batches to process")

        # Define worker function for processing a single batch
        def process_batch(batch_id):
            logger = logging.getLogger(f"ecod.hhsearch.batch_{batch_id}")

            try:
                # Create a new context for the worker
                worker_context = ApplicationContext(self.context.config_manager.config_path)

                # Initialize HHResult registrar
                registrar = HHResultRegistrar(worker_context)

                logger.info(f"Processing HHR files for batch {batch_id}")

                # Process batch
                start_time = datetime.now()
                result = registrar.register_batch_results(batch_id, args.force)
                end_time = datetime.now()

                logger.info(f"Successfully processed {result} HHR files for batch {batch_id}")

                return {
                    "batch_id": batch_id,
                    "processed": result,
                    "success": result > 0,
                    "duration": (end_time - start_time).total_seconds()
                }

            except Exception as e:
                logger.error(f"Error processing batch {batch_id}: {str(e)}")

                return {
                    "batch_id": batch_id,
                    "processed": 0,
                    "success": False,
                    "error": str(e)
                }

        # Process batches in parallel
        results = []
        with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
            futures = [executor.submit(process_batch, batch_id) for batch_id in batch_ids]

            for future in as_completed(futures):
                try:
                    result = future.result()
                    results.append(result)

                    if result['success']:
                        self.logger.info(f"Batch {result['batch_id']} completed: {result['processed']} files in {result['duration']:.1f} seconds")
                    else:
                        self.logger.error(f"Batch {result['batch_id']} failed: {result.get('error', 'Unknown error')}")

                except Exception as e:
                    self.logger.error(f"Error in batch processing: {str(e)}")

        # Log summary
        success_count = sum(1 for r in results if r['success'])
        total_processed = sum(r['processed'] for r in results if r['success'])

        self.logger.info(f"Batch processing complete:")
        self.logger.info(f"  - {success_count}/{len(batch_ids)} batches successful")
        self.logger.info(f"  - {total_processed} total files processed")

        # Write detailed results to output file if specified
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(results, f, indent=2)
            self.logger.info(f"Detailed results written to {args.output}")

        return 0 if success_count > 0 else 1

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

# Add these functions to maintain compatibility with main.py

def setup_parser(parser: argparse.ArgumentParser) -> None:
    """Set up the argument parser for HHSearch commands"""
    # Create temporary command to setup parser
    cmd = HHSearchCommand()
    cmd.setup_parser(parser)

def run_command(args: argparse.Namespace) -> int:
    """Run the specified HHSearch command"""
    cmd = HHSearchCommand(args.config)
    return cmd.run_command(args)
