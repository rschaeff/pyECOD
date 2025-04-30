#!/usr/bin/env python3
"""
collate_hhsearch_results.py - Collate HHSearch and BLAST results for domain analysis

This script integrates HHSearch results with BLAST data to create comprehensive
domain summaries for protein chains in a specified batch.
"""

import os
import sys
import logging
import argparse
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.config import ConfigManager
from ecod.core.context import ApplicationContext
from ecod.pipelines.domain_analysis.summary import DomainSummary
from ecod.pipelines.hhsearch.processor import HHSearchProcessor
from ecod.utils.path_utils import (
    get_standardized_paths,
    get_file_db_path,
    resolve_file_path,
    find_files_with_legacy_paths,
    migrate_file_to_standard_path
)


def setup_logging(verbose: bool = False, log_file: Optional[str] = None) -> None:
    """Configure logging with appropriate handlers and format"""
    log_level = logging.DEBUG if verbose else logging.INFO

    handlers = [logging.StreamHandler()]
    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )


class CollationRunner:
    """Runner for collating HHSearch and BLAST results into domain summaries"""

    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize the collation runner with configuration

        Args:
            config_path: Path to configuration file
        """
        self.logger = logging.getLogger("ecod.collation_runner")

        # Initialize application context (newer pattern)
        self.context = ApplicationContext(config_path)
        self.config = self.context.config
        self.db = self.context.get_db()

        # Initialize processors
        self.domain_summary = DomainSummary(self.context)
        self.hhsearch_processor = HHSearchProcessor(self.context)

    def get_db(self):
        """Access method for database - for backward compatibility"""
        return self.db

    def is_force_overwrite(self) -> bool:
        """Check if force overwrite is enabled - for backward compatibility"""
        return self.config.get('pipeline', {}).get('force_overwrite', False)

    def collate_batch_results(self, batch_id: int, force: bool = False,
                            limit: Optional[int] = None) -> bool:
        """
        Collate HHSearch results with BLAST results for representative processes in a batch

        Args:
            batch_id: ID of the batch to process
            force: Whether to force reprocessing of already processed chains
            limit: Maximum number of chains to process

        Returns:
            True if successful, False otherwise
        """
        self.logger.info(f"Collating HHSearch and BLAST results for batch {batch_id}")

        # Get batch info
        batch_query = """
        SELECT id, batch_name, base_path, ref_version
        FROM ecod_schema.batch
        WHERE id = %s
        """

        batch_info = self.db.execute_dict_query(batch_query, (batch_id,))
        if not batch_info:
            self.logger.error(f"Batch {batch_id} not found")
            return False

        base_path = batch_info[0]['base_path']
        ref_version = batch_info[0]['ref_version']
        batch_name = batch_info[0]['batch_name']

        self.logger.info(f"Processing batch: {batch_name} with reference {ref_version}")

        # Get representative proteins with HHSearch results
        protein_query = """
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

        proteins = self.db.execute_dict_query(protein_query, (batch_id,))
        self.logger.info(f"Found {len(proteins)} representative proteins with HHSearch results")

        if limit and limit < len(proteins):
            proteins = proteins[:limit]
            self.logger.info(f"Limited to {limit} proteins")

        success_count = 0
        for protein in proteins:
            pdb_id = protein['pdb_id']
            chain_id = protein['chain_id']
            process_id = protein['process_id']

            self.logger.info(f"Processing {pdb_id}_{chain_id}")

            # Use standardized paths from path_utils
            paths = get_standardized_paths(
                base_path,
                pdb_id,
                chain_id,
                ref_version,
                create_dirs=True
            )

            # Check if full pipeline domain partition exists
            if os.path.exists(paths['domain_partition']) and not force:
                # Verify this is a full pipeline summary (contains HHSearch evidence)
                is_full_summary = self._check_for_hhsearch_evidence(paths['domain_partition'])

                if is_full_summary:
                    self.logger.info(f"Full pipeline domain summary already exists for {pdb_id}_{chain_id}, skipping")
                    success_count += 1
                    continue
                else:
                    self.logger.info(f"Found blast-only summary for {pdb_id}_{chain_id}, replacing with full pipeline version")

            try:
                # Process chain using HHSearchProcessor
                result = self.hhsearch_processor._process_chain(
                    pdb_id,
                    chain_id,
                    process_id,
                    batch_info[0],
                    ref_version,
                    force
                )

                if result:
                    success_count += 1
                    self.logger.info(f"Successfully processed {pdb_id}_{chain_id}")
                else:
                    self.logger.warning(f"Failed to process {pdb_id}_{chain_id}")

            except Exception as e:
                self.logger.error(f"Error processing {pdb_id}_{chain_id}: {str(e)}")

        self.logger.info(f"Successfully collated results for {success_count}/{len(proteins)} proteins")
        return success_count > 0

    def _check_for_hhsearch_evidence(self, summary_file: str) -> bool:
        """
        Check if a domain summary contains HHSearch evidence

        Args:
            summary_file: Path to summary file

        Returns:
            True if contains HHSearch evidence, False otherwise
        """
        try:
            import xml.etree.ElementTree as ET
            tree = ET.parse(summary_file)
            root = tree.getroot()

            # Look for HHSearch evidence section
            hhsearch_elem = root.find(".//hhsearch_evidence")
            if hhsearch_elem is None:
                return False

            # Check if it has any hit elements
            hits = hhsearch_elem.findall(".//hh_hit")
            return len(hits) > 0
        except Exception as e:
            self.logger.error(f"Error checking summary for HHSearch evidence: {str(e)}")
            return False


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Collate ECOD HHSearch Results with BLAST')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to process')
    parser.add_argument('--limit', type=int,
                      help='Limit the number of proteins to process')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    parser.add_argument('--force', action='store_true',
                      help='Force reprocessing of already processed results')

    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)

    logger = logging.getLogger("main")
    logger.info(f"Starting collation of BLAST and HHSearch results for batch {args.batch_id}")

    # Use absolute path for config if provided
    if args.config and not os.path.isabs(args.config):
        args.config = os.path.abspath(args.config)

    runner = CollationRunner(args.config)

    success = runner.collate_batch_results(args.batch_id, args.force, args.limit)

    if success:
        logger.info(f"Successfully collated results for batch {args.batch_id}")
        return 0
    else:
        logger.error(f"Failed to collate results for batch {args.batch_id}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
