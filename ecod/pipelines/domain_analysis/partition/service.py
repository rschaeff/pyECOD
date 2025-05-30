#!/usr/bin/env python3
"""
High-level service interface for domain partitioning.

This module provides the main API for domain partitioning,
handling both single protein and batch processing.
"""

import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple, Union
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import os

from ecod.core.context import ApplicationContext
from ecod.db import DBManager
from ecod.models.pipeline.partition import DomainPartitionResult
from ecod.exceptions import PipelineError, ValidationError

from .models import (
    PartitionOptions, PartitionContext, BatchPartitionResults,
    PartitionStage, ProcessingMode, ValidationLevel
)
from .analyzer import EvidenceAnalyzer
from .processor import PartitionProcessor
from .tracker import StatusTracker


class DomainPartitionService:
    """
    High-level service for domain partitioning.

    Provides a clean API for processing proteins and identifying
    domain boundaries from evidence summaries.
    """

    def __init__(self, context: ApplicationContext,
                 service_config: Optional[Dict[str, Any]] = None):
        """
        Initialize the domain partition service.

        Args:
            context: Application context with configuration
            service_config: Optional service-specific configuration
        """
        self.context = context
        self.config = context.config_manager.config
        self.service_config = service_config or {}
        self.logger = logging.getLogger(__name__)

        # Initialize database
        try:
            db_config = context.config_manager.get_db_config()
            self.db = DBManager(db_config)
        except Exception as e:
            self.logger.warning(f"Database initialization failed: {e}")
            self.db = None  # Service can still work without DB for some operations

        # Create partition options from configuration
        self.default_options = self._create_default_options()

        # Initialize components
        self.analyzer = EvidenceAnalyzer(self.default_options)
        self.processor = PartitionProcessor(self.default_options, self.analyzer)
        self.tracker = StatusTracker(self.db)

        # Service configuration
        self.service_settings = {
            'max_workers': self.service_config.get('max_workers', 4),
            'use_multiprocessing': self.service_config.get('use_multiprocessing', False),
            'batch_size': self.service_config.get('batch_size', 100),
            'save_intermediate': self.service_config.get('save_intermediate', True),
            'track_status': self.service_config.get('track_status', True)
        }

        # Service statistics
        self.service_stats = {
            'proteins_processed': 0,
            'domains_found': 0,
            'peptides_found': 0,
            'unclassified': 0,
            'errors': 0,
            'start_time': datetime.now()
        }

        self.logger.info("DomainPartitionService initialized")

    def _create_default_options(self) -> PartitionOptions:
        """Create default partition options from configuration"""
        partition_config = self.config.get('partition', {})

        # Get coverage settings
        coverage_config = partition_config.get('reference_coverage', {})

        options_dict = {
            # ... existing options ...

            # Reference coverage settings
            'min_reference_coverage': coverage_config.get('min_coverage', 0.7),
            'strict_reference_coverage': coverage_config.get('strict_coverage', 0.9),
            'partial_coverage_threshold': coverage_config.get('partial_threshold', 0.3),
            'extend_to_reference_size': coverage_config.get('extend_to_reference', True),
            'reference_size_tolerance': coverage_config.get('size_tolerance', 0.15),
            'max_extension_length': coverage_config.get('max_extension', 50),
            'use_ungapped_coverage': coverage_config.get('ungapped_coverage', True),
            'combine_partial_evidence': coverage_config.get('combine_partial', True),
        }

        # Create and validate options
        options = PartitionOptions(**options_dict)
        options.validate()

        return options

    def partition_protein(self, pdb_id: str, chain_id: str,
                         summary_path: str, output_dir: str,
                         process_id: Optional[int] = None,
                         **options) -> DomainPartitionResult:
        """Fixed partition_protein method with correct API calls"""
        start_time = datetime.now()

        # Create context
        context = PartitionContext(
            pdb_id=pdb_id,
            chain_id=chain_id,
            reference=self._get_reference(),
            output_dir=Path(output_dir),
            process_id=process_id,
            start_time=start_time
        )

        # Create options
        partition_options = self._create_options(**options)

        # Update status if tracking
        if process_id and self.service_settings['track_status']:
            self.tracker.update_process_status(
                process_id,
                stage="domain_partition_processing",
                status="processing"
            )

        try:
            self.logger.info(f"Partitioning protein {context.protein_id}")

            # Stage 1: Parse domain summary
            context.record_stage_time(PartitionStage.LOADING_SUMMARY)
            summary_data = self.analyzer.parse_domain_summary(summary_path)

            if "error" in summary_data:
                raise PipelineError(f"Failed to parse summary: {summary_data['error']}")

            # Update context with sequence info
            if context.sequence_length == 0:
                # Fallback to database lookup
                try:
                    query = """
                    SELECT sequence_length
                    FROM pdb_analysis.protein
                    WHERE pdb_id = %s AND chain_id = %s
                    """
                    rows = self.db.execute_dict_query(query, (pdb_id, chain_id))
                    if rows:
                        context.sequence_length = rows[0]['sequence_length']
                        self.logger.debug(f"Got sequence length from database: {context.sequence_length}")
                    else:
                        self.logger.warning(f"No sequence length found for {pdb_id}_{chain_id}")
                except Exception as e:
                    self.logger.warning(f"Could not get sequence length from database: {e}")
                    context.sequence_length = 0

            # Continue with rest of processing...
            context.record_stage_time(PartitionStage.IDENTIFYING_BOUNDARIES)
            result = self.processor.process_evidence(evidence_list, context)

            # Check for peptide
            if summary_data.get("is_peptide", False):
                self.logger.info(f"{context.protein_id} is marked as peptide")
                result = DomainPartitionResult(
                    pdb_id=pdb_id,
                    chain_id=chain_id,
                    reference=context.reference,
                    is_peptide=True,
                    is_unclassified=True,
                    success=True,
                    domain_file=str(context.output_file)
                )

                self._finalize_result(result, context, process_id)
                self.service_stats['peptides_found'] += 1
                return result

            # Stage 2: Extract and validate evidence - FIXED METHOD CALL
            context.record_stage_time(PartitionStage.EXTRACTING_EVIDENCE)
            evidence_list = self.analyzer.extract_evidence_with_classification(
                summary_data,
                use_cache=partition_options.use_cache,
                db_lookup_func=self._get_domain_classification if partition_options.use_cache else None
            )

            # Continue with rest of processing...
            context.record_stage_time(PartitionStage.IDENTIFYING_BOUNDARIES)
            # Before line 184, add:

            result = self.processor.process_evidence(evidence_list, context)

            # Stage 4: Save results
            if self.service_settings['save_intermediate'] or partition_options.save_intermediate:
                context.record_stage_time(PartitionStage.SAVING_RESULTS)
                try:
                    success = result.save(output_dir=output_dir)
                    if not success:
                        self.logger.warning(f"Failed to save partition results for {context.protein_id}")
                except Exception as e:
                    self.logger.warning(f"Error saving partition results for {context.protein_id}: {e}")

            # Finalize
            context.record_stage_time(PartitionStage.COMPLETE)
            result.processing_time = context.get_total_time()

            self._finalize_result(result, context, process_id)

            # Update service statistics
            self.service_stats['proteins_processed'] += 1
            if result.domains:
                self.service_stats['domains_found'] += len(result.domains)
            elif result.is_unclassified:
                self.service_stats['unclassified'] += 1

            self.logger.info(
                f"Completed partitioning {context.protein_id} in {result.processing_time:.2f}s. "
                f"Found {len(result.domains)} domains"
            )

            return result

        except Exception as e:
            self.logger.error(f"Error partitioning {context.protein_id}: {str(e)}", exc_info=True)
            self.service_stats['errors'] += 1

            # Update status if tracking - FIXED METHOD CALL
            if process_id and self.service_settings['track_status']:
                self.tracker.update_process_status(
                    process_id,
                    stage="domain_partition_failed",
                    status="error",
                    error_message=str(e)  # Now passed via kwargs
                )

            # Return failed result
            return DomainPartitionResult(
                pdb_id=pdb_id,
                chain_id=chain_id,
                reference=self._get_reference(),
                success=False,
                error=str(e)
            )

    def partition_batch(self, batch_id: int, batch_path: str,
                       limit: Optional[int] = None,
                       representatives_only: bool = False,
                       **options) -> BatchPartitionResults:
        """
        Partition a batch of proteins.

        Args:
            batch_id: Batch ID in database
            batch_path: Base path for batch files
            limit: Optional limit on number of proteins to process
            representatives_only: Process only representative proteins
            **options: Additional processing options

        Returns:
            BatchPartitionResults with all results
        """
        self.logger.info(f"Starting batch partition for batch {batch_id}")
        batch_start = datetime.now()

        # Initialize results
        results = BatchPartitionResults(start_time=batch_start)

        # Get proteins to process
        proteins = self._get_proteins_to_process(batch_id, limit, representatives_only)

        if not proteins:
            self.logger.warning(f"No proteins to process in batch {batch_id}")
            results.finalize()
            return results

        self.logger.info(f"Processing {len(proteins)} proteins from batch {batch_id}")

        # Create process map for status tracking
        process_map = {f"{p['pdb_id']}_{p['chain_id']}": p['process_id'] for p in proteins}

        # Update non-representative status if needed
        if representatives_only:
            self.tracker.update_non_representative_status(batch_id)

        # Process proteins
        if self.service_settings['use_multiprocessing'] and len(proteins) > 1:
            self._process_batch_parallel(proteins, batch_path, process_map, results, **options)
        else:
            self._process_batch_sequential(proteins, batch_path, process_map, results, **options)

        # Finalize results
        results.finalize()

        # Update batch status
        self.tracker.update_batch_completion_status(
            batch_id,
            representatives_only=representatives_only
        )

        self.logger.info(
            f"Batch {batch_id} completed in {results.processing_time:.2f}s. "
            f"Success: {results.success_count}/{results.total} ({results.success_rate:.1f}%)"
        )

        return results

    def _process_batch_sequential(self, proteins: List[Dict[str, Any]],
                                 batch_path: str, process_map: Dict[str, int],
                                 results: BatchPartitionResults, **options) -> None:
        """Process batch sequentially"""
        partition_options = self._create_options(**options)

        for i, protein in enumerate(proteins):
            pdb_id = protein['pdb_id']
            chain_id = protein['chain_id']
            process_id = protein.get('process_id')

            try:
                # Find summary file
                summary_path = self._find_domain_summary(
                    batch_path, pdb_id, chain_id,
                    blast_only=partition_options.blast_only
                )

                if not summary_path:
                    raise FileNotFoundError(f"Domain summary not found for {pdb_id}_{chain_id}")

                # Process protein
                result = self.partition_protein(
                    pdb_id, chain_id, summary_path, batch_path,
                    process_id=process_id, **options
                )

                results.add_result(result)

                # Log progress
                if (i + 1) % 10 == 0:
                    self.logger.info(f"Processed {i + 1}/{len(proteins)} proteins")

            except Exception as e:
                self.logger.error(f"Error processing {pdb_id}_{chain_id}: {e}")

                # Create error result
                error_result = DomainPartitionResult(
                    pdb_id=pdb_id,
                    chain_id=chain_id,
                    reference=self._get_reference(),
                    success=False,
                    error=str(e)
                )
                results.add_result(error_result)

    def _process_batch_parallel(self, proteins: List[Dict[str, Any]],
                               batch_path: str, process_map: Dict[str, int],
                               results: BatchPartitionResults, **options) -> None:
        """Process batch in parallel"""
        max_workers = min(self.service_settings['max_workers'], len(proteins))

        # Use ThreadPoolExecutor for I/O bound tasks
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all tasks
            future_to_protein = {
                executor.submit(
                    self._process_single_protein_wrapper,
                    protein, batch_path, **options
                ): protein
                for protein in proteins
            }

            # Collect results
            completed = 0
            for future in as_completed(future_to_protein):
                protein = future_to_protein[future]
                completed += 1

                try:
                    result = future.result()
                    results.add_result(result)
                except Exception as e:
                    self.logger.error(
                        f"Error processing {protein['pdb_id']}_{protein['chain_id']}: {e}"
                    )

                    # Create error result
                    error_result = DomainPartitionResult(
                        pdb_id=protein['pdb_id'],
                        chain_id=protein['chain_id'],
                        reference=self._get_reference(),
                        success=False,
                        error=str(e)
                    )
                    results.add_result(error_result)

                # Log progress
                if completed % 10 == 0:
                    self.logger.info(f"Completed {completed}/{len(proteins)} proteins")

    def _process_single_protein_wrapper(self, protein: Dict[str, Any],
                                       batch_path: str, **options) -> DomainPartitionResult:
        """Wrapper for processing single protein in parallel execution"""
        pdb_id = protein['pdb_id']
        chain_id = protein['chain_id']
        process_id = protein.get('process_id')

        partition_options = self._create_options(**options)

        # Find summary file
        summary_path = self._find_domain_summary(
            batch_path, pdb_id, chain_id,
            blast_only=partition_options.blast_only
        )

        if not summary_path:
            raise FileNotFoundError(f"Domain summary not found for {pdb_id}_{chain_id}")

        # Process protein
        return self.partition_protein(
            pdb_id, chain_id, summary_path, batch_path,
            process_id=process_id, **options
        )

    def reprocess_failed(self, batch_id: int, batch_path: str,
                        **options) -> BatchPartitionResults:
        """
        Reprocess only failed proteins from a batch.

        Args:
            batch_id: Batch ID
            batch_path: Batch base path
            **options: Processing options

        Returns:
            BatchPartitionResults
        """
        # Get failed proteins
        query = """
        SELECT p.pdb_id, p.chain_id, ps.id as process_id
        FROM ecod_schema.process_status ps
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        WHERE ps.batch_id = %s
          AND ps.current_stage IN ('domain_partition_failed', 'domain_partition_error')
          AND ps.status = 'error'
        """

        try:
            rows = self.db.execute_dict_query(query, (batch_id,))
            proteins = list(rows)

            if not proteins:
                self.logger.info(f"No failed proteins to reprocess in batch {batch_id}")
                return BatchPartitionResults()

            self.logger.info(f"Reprocessing {len(proteins)} failed proteins from batch {batch_id}")

            # Process with force_overwrite
            options['force_overwrite'] = True
            return self.partition_batch(batch_id, batch_path, **options)

        except Exception as e:
            self.logger.error(f"Error reprocessing failed proteins: {e}")
            return BatchPartitionResults()

    def _get_proteins_to_process(self, batch_id: int, limit: Optional[int] = None,
                                representatives_only: bool = False) -> List[Dict[str, Any]]:
        """Get proteins to process from batch"""
        query = """
        SELECT p.id as protein_id, p.pdb_id, p.chain_id,
               ps.id as process_id, ps.is_representative
        FROM ecod_schema.process_status ps
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        LEFT JOIN ecod_schema.process_file pf_summ ON (
            pf_summ.process_id = ps.id AND
            pf_summ.file_type = 'domain_summary' AND
            pf_summ.file_exists = TRUE
        )
        LEFT JOIN ecod_schema.process_file pf_part ON (
            pf_part.process_id = ps.id AND
            pf_part.file_type = 'domain_partition' AND
            pf_part.file_exists = TRUE
        )
        WHERE ps.batch_id = %s
          AND pf_summ.id IS NOT NULL
          AND pf_part.id IS NULL
        """

        params = [batch_id]

        if representatives_only:
            query += " AND ps.is_representative = TRUE"

        query += " ORDER BY p.id"

        if limit is not None:
            query += f" LIMIT {limit}"

        try:
            results = self.db.execute_dict_query(query, tuple(params))
            return list(results)
        except Exception as e:
            self.logger.error(f"Error fetching proteins to process: {e}")
            return []

    def _find_domain_summary(self, batch_path: str, pdb_id: str, chain_id: str,
                           blast_only: bool = False) -> Optional[str]:
        """Find domain summary file"""
        from ecod.utils.path_utils import get_all_evidence_paths

        file_type = 'blast_only_summary' if blast_only else 'domain_summary'

        try:
            evidence_paths = get_all_evidence_paths(
                batch_path, pdb_id, chain_id, self._get_reference()
            )

            if file_type in evidence_paths and evidence_paths[file_type]['exists_at']:
                return evidence_paths[file_type]['exists_at']

        except Exception as e:
            self.logger.warning(f"Error finding domain summary: {e}")

        return None

    def _get_domain_classification(self, domain_id: str) -> Optional[Dict[str, Any]]:
        """Get domain classification from database"""
        query = """
        SELECT t_group, h_group, x_group, a_group,
               is_manual_rep, is_f70, is_f40, is_f99
        FROM pdb_analysis.domain
        WHERE domain_id = %s
        """

        try:
            rows = self.db.execute_dict_query(query, (domain_id,))
            if rows:
                return rows[0]
        except Exception as e:
            self.logger.error(f"Error getting classification for {domain_id}: {e}")

        return None

    def _finalize_result(self, result: DomainPartitionResult,
                        context: PartitionContext,
                        process_id: Optional[int]) -> None:
        """Finalize partition result and update status"""
        # Set processing time
        result.processing_time = context.get_total_time()

        # Update process status if tracking
        if process_id and self.service_settings['track_status']:
            if result.success:
                self.tracker.update_process_status(
                    process_id,
                    stage="domain_partition_complete",
                    status="success"
                )

                # Register output file
                if result.domain_file and os.path.exists(result.domain_file):
                    self.tracker.register_domain_file(
                        process_id,
                        result.domain_file,
                        str(context.output_dir)
                    )
            else:
                self.tracker.update_process_status(
                    process_id,
                    stage="domain_partition_failed",
                    status="error",
                    error_message=result.error
                )

    def _create_options(self, **kwargs) -> PartitionOptions:
        """Create partition options with overrides"""
        # Start with defaults
        options_dict = self.default_options.to_dict()

        # Apply overrides
        options_dict.update(kwargs)

        # Handle special flags
        if hasattr(self.context, 'is_force_overwrite'):
            options_dict['force_overwrite'] = self.context.is_force_overwrite()

        # Create and validate
        options = PartitionOptions.from_dict(options_dict)
        options.validate()

        return options

    def _get_reference(self) -> str:
        """Get reference version"""
        return self.config.get('reference', {}).get('current_version', 'develop291')

    def get_service_statistics(self) -> Dict[str, Any]:
        """Get comprehensive service statistics"""
        runtime = (datetime.now() - self.service_stats['start_time']).total_seconds()

        return {
            'service': {
                **self.service_stats,
                'runtime_seconds': runtime,
                'proteins_per_minute': (self.service_stats['proteins_processed'] / runtime * 60) if runtime > 0 else 0
            },
            'processor': self.processor.get_statistics(),
            'analyzer': self.analyzer.get_cache_statistics(),
            'configuration': {
                'reference': self._get_reference(),
                'default_options': self.default_options.to_dict(),
                'service_settings': self.service_settings
            }
        }

    def clear_all_caches(self) -> None:
        """Clear all service caches"""
        self.analyzer.clear_cache()
        self.processor.clear_cache()
        self.logger.info("Cleared all service caches")

    def validate_setup(self) -> Dict[str, bool]:
        """Validate service setup"""
        validations = {}

        # Check database
        try:
            self.db.test_connection()
            validations['database'] = True
        except Exception as e:
            self.logger.error(f"Database connection failed: {e}")
            validations['database'] = False

        # Check configuration
        validations['config_loaded'] = bool(self.config)
        validations['reference_set'] = bool(self._get_reference())

        # Check options validation
        try:
            self.default_options.validate()
            validations['options_valid'] = True
        except Exception as e:
            self.logger.error(f"Options validation failed: {e}")
            validations['options_valid'] = False

        # Log results
        all_valid = all(validations.values())
        if all_valid:
            self.logger.info("Service validation passed")
        else:
            failed = [k for k, v in validations.items() if not v]
            self.logger.warning(f"Service validation failed for: {failed}")

        return validations

    def __enter__(self):
        """Context manager entry"""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit"""
        self.clear_all_caches()
        if hasattr(self.db, 'close'):
            self.db.close()


# Convenience functions

def create_service(config_path: Optional[str] = None,
                  service_config: Optional[Dict[str, Any]] = None) -> DomainPartitionService:
    """
    Create a domain partition service.

    Args:
        config_path: Optional configuration file path
        service_config: Optional service-specific configuration

    Returns:
        Configured DomainPartitionService
    """
    if config_path:
        context = ApplicationContext(config_path)
    else:
        config_path = os.environ.get('ECOD_CONFIG_PATH', 'config/config.yml')
        context = ApplicationContext(config_path)

    return DomainPartitionService(context, service_config)


def partition_single_protein(pdb_id: str, chain_id: str,
                           summary_path: str, output_dir: str,
                           config_path: Optional[str] = None,
                           **options) -> DomainPartitionResult:
    """
    Convenience function to partition a single protein.

    Args:
        pdb_id: PDB identifier
        chain_id: Chain identifier
        summary_path: Path to domain summary XML
        output_dir: Output directory
        config_path: Optional configuration file path
        **options: Processing options

    Returns:
        DomainPartitionResult
    """
    service = create_service(config_path)
    return service.partition_protein(pdb_id, chain_id, summary_path, output_dir, **options)


def partition_batch(batch_id: int, batch_path: str,
                   config_path: Optional[str] = None,
                   **options) -> BatchPartitionResults:
    """
    Convenience function to partition a batch.

    Args:
        batch_id: Batch ID
        batch_path: Batch base path
        config_path: Optional configuration file path
        **options: Processing options

    Returns:
        BatchPartitionResults
    """
    service = create_service(config_path)
    return service.partition_batch(batch_id, batch_path, **options)
