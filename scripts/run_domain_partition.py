#!/usr/bin/env python3
"""
Enhanced Domain Partition Runner - A simplified script for running domain partition
with better path handling and model-based processing
"""

import os
import sys
import logging
import argparse
from pathlib import Path
from typing import List, Dict, Any, Optional, Union

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.config import ConfigManager
from ecod.core.context import ApplicationContext
from ecod.db import DBManager
from ecod.pipelines.domain_analysis.pipeline import DomainAnalysisPipeline
from ecod.pipelines.domain_analysis.partition import DomainPartition
from ecod.utils.path_utils import get_standardized_paths, get_all_evidence_paths, resolve_file_path
from ecod.models import ProteinResult, PipelineResult, ProteinProcessingResult


def setup_logging(verbose=False, log_file=None):
    """Configure logging"""
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

    # Return logger for convenience
    return logging.getLogger("domain_partition_runner")


class DomainPartitionRunner:
    """Enhanced runner for domain partition with improved path handling and batch processing"""

    def __init__(self, config_path=None):
        """Initialize with configuration"""
        self.logger = logging.getLogger("domain_partition_runner")

        # Create application context with configuration
        self.context = ApplicationContext(config_path)

        # Initialize domain analysis pipeline
        self.domain_pipeline = DomainAnalysisPipeline(self.context)

        # Direct access to partition component for specific operations
        self.partition = DomainPartition(self.context)

    def process_batch(self, batch_id: int, blast_only: bool = False,
                     limit: int = None, reps_only: bool = False) -> PipelineResult:
        """
        Process domain partition for a batch of proteins

        Args:
            batch_id: Batch ID to process
            blast_only: Whether to use only BLAST results (no HHSearch)
            limit: Maximum number of proteins to process
            reps_only: Whether to process only representative proteins

        Returns:
            PipelineResult object containing processing results
        """
        self.logger.info(f"Processing batch {batch_id} (blast_only={blast_only}, limit={limit}, reps_only={reps_only})")

        # Run partition-only pipeline with the batch_id using the DomainAnalysisPipeline
        result = self.domain_pipeline.run_pipeline(
            batch_id=batch_id,
            blast_only=blast_only,
            limit=limit,
            partition_only=True,  # Only run partition, assuming summaries exist
            reps_only=reps_only
        )

        # Log statistics from result
        self._log_pipeline_result(result)

        return result

    def process_specific_proteins(self, batch_id: int, protein_ids: List[int],
                                 blast_only: bool = False) -> ProteinProcessingResult:
        """
        Process domain partition for specific proteins

        Args:
            batch_id: Batch ID
            protein_ids: List of protein IDs to process
            blast_only: Whether to use only BLAST results (no HHSearch)

        Returns:
            ProteinProcessingResult object containing processing results
        """
        self.logger.info(f"Processing {len(protein_ids)} specific proteins for batch {batch_id}")

        # Use the process_proteins method to handle specific proteins
        result = self.domain_pipeline.process_proteins(
            batch_id=batch_id,
            protein_ids=protein_ids,
            blast_only=blast_only,
            partition_only=True  # Only run partition, assuming summaries exist
        )

        # Log results
        self._log_specific_protein_results(result)

        return result

    def process_single_protein(self, pdb_id: str, chain_id: str, batch_path: str,
                              reference: str, blast_only: bool = False) -> Dict[str, Any]:
        """
        Process domain partition for a single protein without database

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            batch_path: Path to batch directory
            reference: Reference version
            blast_only: Whether to use only BLAST results (no HHSearch)

        Returns:
            Dictionary with processing result
        """
        self.logger.info(f"Processing single protein {pdb_id}_{chain_id}")

        # Verify domain summary exists before proceeding
        summary_exists = self._verify_domain_summary(pdb_id, chain_id, batch_path, reference, blast_only)
        if not summary_exists:
            return {
                "success": False,
                "message": "Domain summary file not found - cannot proceed with partition"
            }

        # Use the analyze_domain method to process with models
        result = self.domain_pipeline.analyze_domain(
            pdb_id=pdb_id,
            chain_id=chain_id,
            output_dir=batch_path,
            reference=reference,
            blast_only=blast_only
        )

        return result

    def verify_batch_readiness(self, batch_id: int, blast_only: bool = False) -> Dict[str, Any]:
        """
        Verify that a batch has the necessary domain summaries to run partition

        Args:
            batch_id: Batch ID to verify
            blast_only: Whether to check for blast-only summaries

        Returns:
            Dictionary with verification results
        """
        self.logger.info(f"Verifying batch {batch_id} readiness for domain partition")

        # Get batch information
        batch_info = self._get_batch_info(batch_id)
        if not batch_info:
            return {"ready": False, "error": "Batch not found"}

        # Use the verify method from domain pipeline
        is_ready = self.domain_pipeline._verify_summary_completion(batch_id, blast_only)

        # Count total proteins and ready proteins
        db = self.context.db

        # Query to count proteins with domain summaries
        summary_type = "blast_only_summary" if blast_only else "domain_summary"
        query = """
        SELECT
            COUNT(*) as total_proteins,
            SUM(CASE WHEN EXISTS (
                SELECT 1 FROM ecod_schema.process_file pf
                WHERE pf.process_id = ps.id
                AND pf.file_type = %s
                AND pf.file_exists = TRUE
            ) THEN 1 ELSE 0 END) as ready_proteins
        FROM
            ecod_schema.process_status ps
        WHERE
            ps.batch_id = %s
        """

        counts = db.execute_dict_query(query, (summary_type, batch_id))

        total_proteins = counts[0]['total_proteins'] if counts else 0
        ready_proteins = counts[0]['ready_proteins'] if counts else 0

        # Create detailed result
        result = {
            "ready": is_ready,
            "batch_id": batch_id,
            "batch_name": batch_info.get('batch_name', ''),
            "total_proteins": total_proteins,
            "ready_proteins": ready_proteins,
            "ready_percentage": (ready_proteins / total_proteins * 100) if total_proteins > 0 else 0,
            "summary_type": "blast_only" if blast_only else "full"
        }

        # Log the readiness status
        self.logger.info(f"Batch {batch_id} readiness: {ready_proteins}/{total_proteins} "
                       f"proteins ready ({result['ready_percentage']:.1f}%)")

        return result

    def check_paths(self, pdb_id: str, chain_id: str, batch_id: int = None,
                   batch_path: str = None, reference: str = None) -> Dict[str, Any]:
        """
        Check file paths for a protein to help diagnose issues

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            batch_id: Optional batch ID (if known)
            batch_path: Optional batch path (if known)
            reference: Optional reference version (if known)

        Returns:
            Dictionary with path information
        """
        self.logger.info(f"Checking paths for {pdb_id}_{chain_id}")

        # If batch_id is provided but not batch_path or reference, get them from database
        if batch_id and (not batch_path or not reference):
            batch_info = self._get_batch_info(batch_id)
            if batch_info:
                batch_path = batch_info.get('base_path', batch_path)
                reference = batch_info.get('ref_version', reference)

        if not batch_path or not reference:
            return {"error": "Missing required batch_path or reference - provide batch_id or specify directly"}

        # Get all evidence paths
        all_paths = get_all_evidence_paths(batch_path, pdb_id, chain_id, reference)

        # Create a structured result with file existence status
        result = {
            "pdb_id": pdb_id,
            "chain_id": chain_id,
            "batch_path": batch_path,
            "reference": reference,
            "files": {}
        }

        # Check each file type
        for file_type, paths in all_paths.items():
            standard_path = paths.get('standard_path', '')
            legacy_path = paths.get('legacy_path', '')
            exists_at = paths.get('exists_at', '')

            result["files"][file_type] = {
                "standard_path": standard_path,
                "legacy_path": legacy_path,
                "exists": exists_at is not None,
                "path": exists_at or standard_path,
                "size": os.path.getsize(exists_at) if exists_at and os.path.exists(exists_at) else 0
            }

        # Check process status in database if batch_id provided
        if batch_id:
            result["process_info"] = self._get_process_status(batch_id, pdb_id, chain_id)

        return result

    def reset_failed_processes(self, batch_id: int) -> Dict[str, int]:
        """
        Reset failed domain partition processes so they can be retried

        Args:
            batch_id: Batch ID

        Returns:
            Dictionary with reset counts
        """
        self.logger.info(f"Resetting failed processes for batch {batch_id}")

        # Use the pipeline's reset method
        summary_reset = self.domain_pipeline.reset_failed_processes(batch_id, 'domain_summary_failed')
        partition_reset = self.domain_pipeline.reset_failed_processes(batch_id, 'domain_partition_failed')

        result = {
            "summary_reset": summary_reset,
            "partition_reset": partition_reset,
            "total_reset": summary_reset + partition_reset
        }

        self.logger.info(f"Reset {result['total_reset']} failed processes")

        return result

    def _verify_domain_summary(self, pdb_id: str, chain_id: str, batch_path: str,
                             reference: str, blast_only: bool = False) -> bool:
        """Verify that domain summary file exists for a protein"""
        summary_type = 'blast_only_summary' if blast_only else 'domain_summary'

        # Use path_utils to get all possible paths
        all_paths = get_all_evidence_paths(batch_path, pdb_id, chain_id, reference)

        if summary_type in all_paths and all_paths[summary_type]['exists_at']:
            path = all_paths[summary_type]['exists_at']
            self.logger.info(f"Found domain summary at: {path}")
            return True

        self.logger.warning(f"No domain summary found for {pdb_id}_{chain_id}")
        return False

    def _get_batch_info(self, batch_id: int) -> Optional[Dict[str, Any]]:
        """Get batch information from database"""
        # Access db using the correct property (db, not get_db())
        db = self.context.db

        query = """
        SELECT id, batch_name, base_path, ref_version
        FROM ecod_schema.batch
        WHERE id = %s
        """

        try:
            rows = db.execute_dict_query(query, (batch_id,))
            if rows:
                return rows[0]
            else:
                self.logger.error(f"Batch {batch_id} not found")
                return None
        except Exception as e:
            self.logger.error(f"Error retrieving batch information: {e}")
            return None

    def _get_process_status(self, batch_id: int, pdb_id: str, chain_id: str) -> Optional[Dict[str, Any]]:
        """Get process status information from database"""
        # Access db using the correct property
        db = self.context.db

        query = """
        SELECT
            ps.id as process_id,
            ps.protein_id,
            ps.current_stage,
            ps.status,
            ps.error_message,
            ps.is_representative
        FROM
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        WHERE
            ps.batch_id = %s
            AND p.pdb_id = %s
            AND p.chain_id = %s
        """

        try:
            rows = db.execute_dict_query(query, (batch_id, pdb_id, chain_id))
            if rows:
                result = rows[0]

                # Add file information
                file_query = """
                SELECT
                    file_type,
                    file_path,
                    file_exists,
                    file_size
                FROM
                    ecod_schema.process_file
                WHERE
                    process_id = %s
                """

                files = db.execute_dict_query(file_query, (result['process_id'],))
                result['files'] = {file['file_type']: file for file in files}

                return result
            else:
                self.logger.warning(f"No process found for {pdb_id}_{chain_id} in batch {batch_id}")
                return None
        except Exception as e:
            self.logger.error(f"Error retrieving process status: {e}")
            return None

    def _log_pipeline_result(self, result: PipelineResult) -> None:
        """Log pipeline result statistics"""
        self.logger.info(f"Batch processing result: success={result.success}")

        if hasattr(result, 'partition_stats'):
            stats = result.partition_stats
            self.logger.info(f"Partition stats: {stats.get('files_created', 0)} files created")

            # Log errors if any
            if 'errors' in stats and stats['errors']:
                for error in stats['errors'][:5]:  # Show first 5 errors
                    self.logger.error(f"Error: {error}")

                if len(stats['errors']) > 5:
                    self.logger.error(f"... and {len(stats['errors']) - 5} more errors")

    def _log_specific_protein_results(self, result: ProteinProcessingResult) -> None:
        """Log statistics for specific protein processing"""
        self.logger.info(f"Processing result: {result.success_count}/{result.total_count} proteins successful")

        # Log detailed results for each protein
        for protein in result.protein_results:
            status = "SUCCESS" if protein.success else "FAILED"
            self.logger.info(f"{protein.pdb_id}_{protein.chain_id}: {status}")

            if not protein.success and protein.error:
                self.logger.error(f"  Error: {protein.error}")


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Run ECOD Domain Partition with Enhanced Features')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')

    # Basic operation modes
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--batch-id', type=int,
                     help='Batch ID to process')
    group.add_argument('--process-id', type=int, nargs='+',
                     help='One or more process IDs to process')
    group.add_argument('--verify', type=int,
                     help='Verify batch readiness (provide batch ID)')
    group.add_argument('--check-paths', action='store_true',
                     help='Check paths for a specific protein (requires --pdb-id and --chain-id)')
    group.add_argument('--reset', type=int,
                     help='Reset failed processes for batch (provide batch ID)')

    # Protein identifiers for single protein operations
    parser.add_argument('--pdb-id', type=str,
                      help='PDB ID for checking paths or processing single protein')
    parser.add_argument('--chain-id', type=str,
                      help='Chain ID for checking paths or processing single protein')

    # Processing options
    parser.add_argument('--blast-only', action='store_true',
                      help='Use only BLAST results (no HHSearch)')
    parser.add_argument('--limit', type=int,
                      help='Limit the number of proteins to process')
    parser.add_argument('--reps-only', action='store_true',
                      help='Process only representative proteins')

    # Output options
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')

    args = parser.parse_args()

    # Setup logging
    logger = setup_logging(args.verbose, args.log_file)

    # Use absolute path for config if provided
    if args.config and not os.path.isabs(args.config):
        args.config = os.path.abspath(args.config)

    # Initialize runner
    runner = DomainPartitionRunner(args.config)

    # Handle different operation modes

    # Check paths for a protein
    if args.check_paths:
        if not args.pdb_id or not args.chain_id:
            logger.error("--check-paths requires --pdb-id and --chain-id")
            return 1

        batch_id = args.batch_id if args.batch_id else None
        result = runner.check_paths(args.pdb_id, args.chain_id, batch_id)

        # Log detailed path information
        logger.info(f"Path check for {args.pdb_id}_{args.chain_id}:")
        for file_type, info in result.get('files', {}).items():
            status = "EXISTS" if info.get('exists', False) else "MISSING"
            path = info.get('path', '')
            logger.info(f"  {file_type:20s}: {status:8s} {path}")

        # Exit with success if domain summary exists
        domain_summary_exists = result.get('files', {}).get('domain_summary', {}).get('exists', False)
        return 0 if domain_summary_exists else 1

    # Verify batch readiness
    if args.verify:
        result = runner.verify_batch_readiness(args.verify, args.blast_only)

        if result.get('ready', False):
            logger.info(f"Batch {args.verify} is READY for domain partition")
            logger.info(f"  {result['ready_proteins']}/{result['total_proteins']} "
                      f"proteins ready ({result['ready_percentage']:.1f}%)")
            return 0
        else:
            logger.error(f"Batch {args.verify} is NOT READY for domain partition")
            logger.error(f"  Only {result['ready_proteins']}/{result['total_proteins']} "
                       f"proteins ready ({result['ready_percentage']:.1f}%)")
            return 1

    # Reset failed processes
    if args.reset:
        result = runner.reset_failed_processes(args.reset)

        logger.info(f"Reset {result['total_reset']} failed processes:")
        logger.info(f"  Domain summary failures: {result['summary_reset']}")
        logger.info(f"  Domain partition failures: {result['partition_reset']}")

        return 0 if result['total_reset'] > 0 else 1

    # Process specific process IDs
    if args.process_id:
        # Use the first value from process_id list to get batch ID if not provided
        process_id = args.process_id[0]

        # Find the batch ID if not explicitly provided
        batch_id = args.batch_id
        if not batch_id:
            # Query database for batch ID of this process
            db = runner.context.get_db()
            query = "SELECT batch_id FROM ecod_schema.process_status WHERE id = %s"
            rows = db.execute_query(query, (process_id,))

            if rows:
                batch_id = rows[0][0]
                logger.info(f"Found batch ID {batch_id} for process {process_id}")
            else:
                logger.error(f"Could not find batch ID for process {process_id}")
                return 1

        # Process the specified process IDs
        result = runner.process_specific_proteins(
            batch_id=batch_id,
            protein_ids=args.process_id,
            blast_only=args.blast_only
        )

        return 0 if result.success else 1

    # Process a batch
    if args.batch_id:
        result = runner.process_batch(
            batch_id=args.batch_id,
            blast_only=args.blast_only,
            limit=args.limit,
            reps_only=args.reps_only
        )

        return 0 if result.success else 1

    # Should never reach here due to required argument group
    return 1


if __name__ == "__main__":
    sys.exit(main())
