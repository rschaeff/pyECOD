#!/usr/bin/env python3
"""
Enhanced Domain Partition Runner - A simplified script for running domain partition
with better path handling and model-based processing
"""

import os
import sys
import logging
import argparse
from datetime import datetime
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
from ecod.models.domain_analysis.domain_model import DomainModel


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
                     limit: int = None, reps_only: bool = False
    ) -> PipelineResult:
        """
        Process domains for batch with proper handling of DomainPartitionResult

        Args:
            batch_id: Batch ID
            blast_only: Whether to use only BLAST results (no HHSearch)
            limit: Maximum number of proteins to process
            reps_only: Whether to process only representative proteins

        Returns:
            PipelineResult with processing results
        """
        self.logger.info(f"Processing batch {batch_id} (blast_only={blast_only}, limit={limit}, reps_only={reps_only})")

        # Create result object
        result = PipelineResult(batch_id=batch_id)

        # Retrieve batch_info internally
        batch_info = self._get_batch_info(batch_id)
        if not batch_info:
            result.success = False
            result.error = f"Batch ID {batch_id} not found"
            return result

        result.set_batch_info(batch_info)

        # Get partition component
        partition = DomainPartition(self.context)

        # Get batch path and reference from batch_info
        batch_path = batch_info.get('base_path')
        reference = batch_info.get('ref_version')

        # Run the domain partition process
        try:
            self.logger.debug(f"Running partition.process_batch with path={batch_path}, reference={reference}")
            partition_results = partition.process_batch(
                batch_id,
                batch_path,
                reference,
                blast_only,
                limit,
                reps_only
            )

            # Process the results correctly based on the type
            if partition_results is None:
                self.logger.warning("Process batch returned no results")
                result.success = False
                result.error = "No results returned from process_batch"
                result.partition_stats = {
                    "files_created": 0,
                    "total_proteins": 0,
                    "errors": ["No results returned from domain partition"]
                }
            elif isinstance(partition_results, list):
                # Handle list of DomainPartitionResult objects
                success_count = sum(1 for r in partition_results if r.success)
                failed_count = sum(1 for r in partition_results if not r.success)

                result.success = success_count > 0
                result.partition_stats = {
                    "files_created": success_count,
                    "total_proteins": len(partition_results),
                    "errors": [r.error for r in partition_results if not r.success and r.error]
                }
                self.logger.info(f"Results: {success_count} succeeded, {failed_count} failed")
            else:
                # Fallback for backward compatibility
                self.logger.debug(f"Unexpected result type from process_batch: {type(partition_results)}")
                result.success = bool(partition_results)
                result.partition_stats = {
                    "files_created": 0 if not result.success else 1,
                    "total_proteins": 0 if not result.success else 1
                }

        except Exception as e:
            self.logger.error(f"Exception in process_batch: {str(e)}")
            import traceback
            self.logger.debug(traceback.format_exc())
            result.success = False
            result.error = f"Exception in process_batch: {str(e)}"
            result.partition_stats = {
                "files_created": 0,
                "total_proteins": 0,
                "errors": [str(e)]
            }

        return result

    def process_specific_proteins(self, batch_id: int, protein_ids: List[int],
                                 blast_only: bool = False
    ) -> ProteinProcessingResult:
        """
        Process domain partition for specific proteins

        Args:
            batch_id: Batch ID
            protein_ids: List of protein IDs to process
            blast_only: Whether to use only BLAST results (no HHSearch)

        Returns:
            ProteinProcessingResult object containing processing results
        """
        self.logger.info(f"Processing {len(protein_ids)} proteins for batch {batch_id}")

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
                              reference: str, blast_only: bool = False
    ) -> Dict[str, Any]:
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
        self.logger.info(f"Processing {pdb_id}_{chain_id}")

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
        self.logger.info(f"Verifying batch {batch_id} readiness")

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

        # Only log the summary information
        self.logger.debug(f"Batch {batch_id} readiness details: {result}")

        return result

    def check_paths(self, pdb_id: str, chain_id: str, batch_id: int = None,
                   batch_path: str = None, reference: str = None
    ) -> Dict[str, Any]:
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

        return result

    def _verify_domain_summary(self, pdb_id: str, chain_id: str, batch_path: str,
                             reference: str, blast_only: bool = False
    ) -> bool:
        """Verify that domain summary file exists for a protein"""
        summary_type = 'blast_only_summary' if blast_only else 'domain_summary'

        # Use path_utils to get all possible paths
        all_paths = get_all_evidence_paths(batch_path, pdb_id, chain_id, reference)

        if summary_type in all_paths and all_paths[summary_type]['exists_at']:
            path = all_paths[summary_type]['exists_at']
            self.logger.debug(f"Found domain summary at: {path}")
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
                self.logger.debug(f"Found batch info: {batch_id} - {rows[0]['batch_name']}")
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

                self.logger.debug(f"Found process status for {pdb_id}_{chain_id}: {result['current_stage']}/{result['status']}")
                return result
            else:
                self.logger.debug(f"No process found for {pdb_id}_{chain_id} in batch {batch_id}")
                return None
        except Exception as e:
            self.logger.error(f"Error retrieving process status: {e}")
            return None

    def _log_pipeline_result(self, result: PipelineResult) -> None:
        """Log pipeline result statistics"""
        self.logger.debug(f"Batch processing result: success={result.success}")

        if hasattr(result, 'partition_stats'):
            stats = result.partition_stats
            self.logger.debug(f"Partition stats: {stats.get('files_created', 0)} files created")

            # Log errors if any
            if 'errors' in stats and stats['errors']:
                error_count = len(stats['errors'])
                if error_count > 2:
                    self.logger.error(f"{error_count} errors occurred, first 2 shown below:")
                    for error in stats['errors'][:2]:
                        self.logger.error(f"Error: {error}")
                else:
                    for error in stats['errors']:
                        self.logger.error(f"Error: {error}")

    def _log_specific_protein_results(self, result: ProteinProcessingResult) -> None:
        """Log statistics for specific protein processing"""
        self.logger.info(f"Results: {result.success_count}/{result.total_count} proteins successful")

        # Only log failures at the regular info level
        if result.success_count < result.total_count:
            failed_proteins = [p for p in result.protein_results if not p.success]
            if len(failed_proteins) <= 3:
                for protein in failed_proteins:
                    self.logger.warning(f"{protein.pdb_id}_{protein.chain_id} failed: {protein.error}")
            else:
                for protein in failed_proteins[:3]:
                    self.logger.warning(f"{protein.pdb_id}_{protein.chain_id} failed: {protein.error}")
                self.logger.warning(f"...and {len(failed_proteins) - 3} more failures")

        # Log all results at debug level
        for protein in result.protein_results:
            status = "SUCCESS" if protein.success else "FAILED"
            self.logger.debug(f"{protein.pdb_id}_{protein.chain_id}: {status}")


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

        # Log path information - condensed for regular users, detailed in verbose mode
        logger.info(f"Path check for {args.pdb_id}_{args.chain_id}:")

        # Just show whether each file exists or not
        for file_type, info in result.get('files', {}).items():
            status = "EXISTS" if info.get('exists', False) else "MISSING"
            if args.verbose:
                path = info.get('path', '')
                logger.info(f"  {file_type:20s}: {status:8s} {path}")
            else:
                logger.info(f"  {file_type:20s}: {status:8s}")

        # Exit with success if domain summary exists
        domain_summary_exists = result.get('files', {}).get('domain_summary', {}).get('exists', False)
        return 0 if domain_summary_exists else 1

    # Verify batch readiness
    if args.verify:
        result = runner.verify_batch_readiness(args.verify, args.blast_only)

        if result.get('ready', False):
            logger.info(f"Batch {args.verify} is READY for domain partition")
            logger.info(f"  {result['ready_proteins']}/{result['total_proteins']} proteins ready ({result['ready_percentage']:.1f}%)")
            return 0
        else:
            logger.error(f"Batch {args.verify} is NOT READY for domain partition")
            logger.info(f"  Only {result['ready_proteins']}/{result['total_proteins']} proteins ready ({result['ready_percentage']:.1f}%)")
            return 1

    # Reset failed processes
    if args.reset:
        result = runner.reset_failed_processes(args.reset)

        if result['total_reset'] > 0:
            logger.info(f"Reset {result['total_reset']} failed processes ({result['summary_reset']} summary, {result['partition_reset']} partition)")
        else:
            logger.info("No failed processes to reset")

        return 0 if result['total_reset'] > 0 else 1

    # Process specific process IDs
    if args.process_id:
        # Use the first value from process_id list to get batch ID if not provided
        process_id = args.process_id[0]

        # Find the batch ID if not explicitly provided
        batch_id = args.batch_id
        if not batch_id:
            # Query database for batch ID of this process
            db = runner.context.db
            query = "SELECT batch_id FROM ecod_schema.process_status WHERE id = %s"
            rows = db.execute_query(query, (process_id,))

            if rows:
                batch_id = rows[0][0]
                logger.debug(f"Found batch ID {batch_id} for process {process_id}")
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
