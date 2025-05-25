#!/usr/bin/env python3
"""
High-level service interface for domain summary generation.

This module provides the main API for generating domain summaries,
handling both single protein and batch processing with a clean interface.
"""

import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple, Union
from concurrent.futures import ProcessPoolExecutor, as_completed
import os

from ecod.core.context import ApplicationContext
from ecod.db import DBManager
from ecod.models.pipeline import DomainPartitionResult
from ecod.exceptions import PipelineError, ValidationError

from .generator import DomainSummaryGenerator, GeneratorConfig
from .models import (
    ProteinIdentifier, SummaryOptions, EvidenceSummary,
    BatchResults, ProcessingStats
)


class DomainSummaryService:
    """
    High-level service for domain summary generation.

    Provides a clean API for processing proteins and generating
    evidence summaries that can be used for domain partitioning.
    """

    def __init__(self, context: ApplicationContext,
                 service_config: Optional[Dict[str, Any]] = None):
        """
        Initialize the domain summary service.

        Args:
            context: Application context with configuration
            service_config: Optional service-specific configuration
        """
        self.context = context
        self.config = context.config_manager.config
        self.service_config = service_config or {}
        self.logger = logging.getLogger(__name__)

        # Initialize database
        db_config = context.config_manager.get_db_config()
        self.db = DBManager(db_config)

        # Initialize generator with configuration
        generator_config = self._create_generator_config()
        self.generator = DomainSummaryGenerator(
            config=self.config,
            db_manager=self.db,
            generator_config=generator_config
        )

        # Track service statistics
        self.service_stats = ProcessingStats()

        # Batch processing configuration
        self.batch_config = {
            'max_workers': self.service_config.get('max_batch_workers', 1),
            'use_multiprocessing': self.service_config.get('use_multiprocessing', False),
            'batch_size': self.service_config.get('batch_size', 100)
        }

        self.logger.info("DomainSummaryService initialized")

    def _create_generator_config(self) -> GeneratorConfig:
        """Create generator configuration from service config"""
        gen_config = self.service_config.get('generator', {})

        return GeneratorConfig(
            max_workers=gen_config.get('max_workers', 4),
            process_timeout=gen_config.get('process_timeout', 300),
            cache_enabled=gen_config.get('cache_enabled', True),
            validate_inputs=gen_config.get('validate_inputs', True),
            parallel_evidence_collection=gen_config.get('parallel_evidence_collection', True),
            skip_on_error=gen_config.get('skip_on_error', False),
            collect_detailed_stats=gen_config.get('collect_detailed_stats', True)
        )

    def process_protein(self, pdb_id: str, chain_id: str, job_dump_dir: str,
                       **options) -> DomainPartitionResult:
        """
        Process a single protein and return domain partition results.

        Args:
            pdb_id: PDB identifier (4 characters)
            chain_id: Chain identifier
            job_dump_dir: Directory containing job files
            **options: Additional options passed to SummaryOptions

        Returns:
            DomainPartitionResult with evidence and preliminary analysis

        Raises:
            ValidationError: If protein identifier is invalid
            PipelineError: If processing fails
        """
        # Start timing
        start_time = datetime.now()

        try:
            # Create protein identifier
            protein = ProteinIdentifier(
                pdb_id=pdb_id,
                chain_id=chain_id,
                reference=self._get_reference()
            )

            # Create processing options
            summary_options = self._create_summary_options(**options)

            # Initialize generator for this job if needed
            self._ensure_generator_initialized(job_dump_dir)

            # Generate evidence summary
            self.logger.info(f"Processing protein {protein.source_id}")
            evidence_summary = self.generator.generate_summary(protein, summary_options)

            # Convert to partition result
            partition_result = self._convert_to_partition_result(
                evidence_summary,
                job_dump_dir,
                summary_options
            )

            # Update statistics
            self.service_stats.proteins_processed += 1
            if partition_result.success and len(partition_result.domains) > 0:
                self.service_stats.proteins_with_evidence += 1

            # Log result
            processing_time = (datetime.now() - start_time).total_seconds()
            self.logger.info(
                f"Completed processing {protein.source_id} in {processing_time:.2f}s. "
                f"Status: {'SUCCESS' if partition_result.success else 'FAILED'}"
            )

            return partition_result

        except ValidationError:
            # Re-raise validation errors
            raise
        except Exception as e:
            # Wrap other errors
            error_msg = f"Failed to process {pdb_id}_{chain_id}: {str(e)}"
            self.logger.error(error_msg, exc_info=True)
            self.service_stats.errors += 1

            # Return error result
            return DomainPartitionResult(
                pdb_id=pdb_id,
                chain_id=chain_id,
                reference=self._get_reference(),
                success=False,
                error=error_msg
            )

    def process_batch(self, proteins: List[Tuple[str, str]], job_dump_dir: str,
                     **options) -> BatchResults:
        """
        Process multiple proteins in batch.

        Args:
            proteins: List of (pdb_id, chain_id) tuples
            job_dump_dir: Directory containing job files
            **options: Additional options passed to SummaryOptions

        Returns:
            BatchResults with all processing results
        """
        self.logger.info(f"Starting batch processing of {len(proteins)} proteins")
        batch_start = datetime.now()

        # Initialize generator once for the batch
        self._ensure_generator_initialized(job_dump_dir)

        # Create options once
        summary_options = self._create_summary_options(**options)

        # Initialize results
        results = BatchResults()

        # Reset batch statistics
        batch_stats = ProcessingStats()
        batch_stats.start_tracking()

        # Determine processing method
        if self.batch_config['use_multiprocessing'] and len(proteins) > 1:
            # Process in parallel using multiprocessing
            self._process_batch_parallel(
                proteins, job_dump_dir, summary_options, results, batch_stats
            )
        else:
            # Process sequentially
            self._process_batch_sequential(
                proteins, job_dump_dir, summary_options, results, batch_stats
            )

        # Finalize statistics
        batch_stats.end_tracking()

        # Log summary
        self.logger.info(
            f"Batch processing completed in {batch_stats.total_processing_time:.2f}s. "
            f"Success: {results.success_count}/{results.total} "
            f"({batch_stats.success_rate:.1f}%)"
        )

        return results

    def _process_batch_sequential(self, proteins: List[Tuple[str, str]],
                                 job_dump_dir: str, options: SummaryOptions,
                                 results: BatchResults, stats: ProcessingStats) -> None:
        """Process batch sequentially"""
        for i, (pdb_id, chain_id) in enumerate(proteins):
            try:
                # Create protein identifier
                protein = ProteinIdentifier(pdb_id=pdb_id, chain_id=chain_id)

                # Generate summary
                summary = self.generator.generate_summary(protein, options)

                # Add to results
                results.add_success(summary)

                # Update statistics
                stats.proteins_processed += 1
                if summary.has_evidence():
                    stats.proteins_with_evidence += 1

                # Log progress
                if (i + 1) % 10 == 0:
                    self.logger.info(f"Processed {i + 1}/{len(proteins)} proteins")

            except Exception as e:
                self.logger.error(f"Error processing {pdb_id}_{chain_id}: {e}")
                results.add_failure(pdb_id, chain_id, str(e))
                stats.errors += 1

    def _process_batch_parallel(self, proteins: List[Tuple[str, str]],
                               job_dump_dir: str, options: SummaryOptions,
                               results: BatchResults, stats: ProcessingStats) -> None:
        """Process batch in parallel using multiprocessing"""
        max_workers = min(self.batch_config['max_workers'], len(proteins))

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Submit all tasks
            future_to_protein = {
                executor.submit(
                    self._process_single_protein_worker,
                    pdb_id, chain_id, job_dump_dir, options.to_dict()
                ): (pdb_id, chain_id)
                for pdb_id, chain_id in proteins
            }

            # Collect results as they complete
            completed = 0
            for future in as_completed(future_to_protein):
                pdb_id, chain_id = future_to_protein[future]
                completed += 1

                try:
                    summary_dict = future.result()

                    # Reconstruct summary from dict
                    summary = self._reconstruct_summary(summary_dict)
                    results.add_success(summary)

                    # Update statistics
                    stats.proteins_processed += 1
                    if summary.has_evidence():
                        stats.proteins_with_evidence += 1

                except Exception as e:
                    self.logger.error(f"Error processing {pdb_id}_{chain_id}: {e}")
                    results.add_failure(pdb_id, chain_id, str(e))
                    stats.errors += 1

                # Log progress
                if completed % 10 == 0:
                    self.logger.info(f"Completed {completed}/{len(proteins)} proteins")

    @staticmethod
    def _process_single_protein_worker(pdb_id: str, chain_id: str,
                                      job_dump_dir: str,
                                      options_dict: Dict[str, Any]) -> Dict[str, Any]:
        """
        Worker function for multiprocessing.

        Static method to avoid pickling issues with instance methods.
        """
        # Create a new service instance in the worker process
        from ecod.core.context import ApplicationContext

        # Get config path from environment or use default
        config_path = os.environ.get('ECOD_CONFIG_PATH', 'config/config.yml')
        context = ApplicationContext(config_path)

        service = DomainSummaryService(context)
        service._ensure_generator_initialized(job_dump_dir)

        # Process protein
        protein = ProteinIdentifier(pdb_id=pdb_id, chain_id=chain_id)
        options = SummaryOptions.from_dict(options_dict)

        summary = service.generator.generate_summary(protein, options)

        # Return as dictionary for pickling
        return summary.get_summary_dict()

    def _reconstruct_summary(self, summary_dict: Dict[str, Any]) -> EvidenceSummary:
        """Reconstruct EvidenceSummary from dictionary"""
        # This is a simplified reconstruction - in practice you might need
        # a more complete from_dict method on EvidenceSummary
        from .models import ProcessingStatus

        summary = EvidenceSummary(
            pdb_id=summary_dict['pdb_id'],
            chain_id=summary_dict['chain_id'],
            reference=summary_dict['reference']
        )

        # Set status
        summary.status = ProcessingStatus[summary_dict['status']]

        # Set metadata
        summary.metadata = summary_dict.get('metadata', {})
        summary.processing_time = summary_dict.get('processing_time', 0.0)

        return summary

    def process_proteins_from_file(self, protein_list_file: str,
                                  job_dump_dir: str, **options) -> BatchResults:
        """
        Process proteins listed in a file.

        Args:
            protein_list_file: Path to file with protein IDs (one per line)
            job_dump_dir: Directory containing job files
            **options: Additional processing options

        Returns:
            BatchResults with all processing results
        """
        # Read protein list
        proteins = self._read_protein_list(protein_list_file)

        if not proteins:
            self.logger.warning(f"No proteins found in {protein_list_file}")
            return BatchResults()

        self.logger.info(f"Read {len(proteins)} proteins from {protein_list_file}")

        # Process batch
        return self.process_batch(proteins, job_dump_dir, **options)

    def _read_protein_list(self, file_path: str) -> List[Tuple[str, str]]:
        """Read protein list from file"""
        proteins = []

        try:
            with open(file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue

                    try:
                        # Try to parse as protein identifier
                        protein = ProteinIdentifier.from_string(line)
                        proteins.append((protein.pdb_id, protein.chain_id))
                    except ValidationError as e:
                        self.logger.warning(f"Skipping invalid protein ID '{line}': {e}")

        except Exception as e:
            self.logger.error(f"Error reading protein list: {e}")

        return proteins

    def _ensure_generator_initialized(self, job_dump_dir: str) -> None:
        """Ensure generator is initialized for the job directory"""
        # Check if already initialized for this directory
        current_job_dir = getattr(self.generator, 'job_dump_dir', None)

        if current_job_dir != Path(job_dump_dir):
            self.logger.debug(f"Initializing generator for job directory: {job_dump_dir}")
            self.generator.initialize_for_job(
                job_dump_dir,
                reference=self._get_reference()
            )

    def _get_reference(self) -> str:
        """Get reference version from configuration"""
        return self.config.get('reference', {}).get('current_version', 'develop291')

    def _create_summary_options(self, **kwargs) -> SummaryOptions:
        """Create SummaryOptions with defaults and overrides"""
        # Start with defaults from config
        defaults = self.config.get('summary_options', {})

        # Apply any overrides
        options_dict = {**defaults, **kwargs}

        # Handle special context flags
        if hasattr(self.context, 'is_force_overwrite'):
            options_dict['force_overwrite'] = self.context.is_force_overwrite()

        # Create and validate options
        options = SummaryOptions(**options_dict)
        options.validate()

        return options

    def _convert_to_partition_result(self, summary: EvidenceSummary,
                                    job_dump_dir: str,
                                    options: SummaryOptions) -> DomainPartitionResult:
        """Convert evidence summary to domain partition result"""
        # Use the built-in conversion method
        result = summary.to_partition_result()

        # Set output file path
        if options.save_intermediate or self.service_config.get('save_summaries', True):
            suffix = ".blast_only" if options.blast_only else ""
            filename = f"{summary.protein_id}.{summary.reference}.evidence_summary{suffix}.xml"

            output_dir = Path(job_dump_dir) / "domains"
            output_dir.mkdir(parents=True, exist_ok=True)

            result.domain_file = str(output_dir / filename)

        return result

    def get_service_statistics(self) -> Dict[str, Any]:
        """Get service-level statistics"""
        stats = {
            'service_stats': self.service_stats.to_dict(),
            'generator_stats': self.generator.get_statistics(),
            'configuration': {
                'reference': self._get_reference(),
                'batch_config': self.batch_config,
                'generator_config': vars(self.generator.generator_config)
            }
        }

        return stats

    def clear_all_caches(self) -> None:
        """Clear all caches in the service"""
        self.generator.clear_caches()
        self.logger.info("Cleared all service caches")

    def validate_setup(self) -> Dict[str, bool]:
        """
        Validate service setup and dependencies.

        Returns:
            Dictionary with validation results
        """
        validations = {}

        # Check database connection
        try:
            self.db.test_connection()
            validations['database'] = True
        except Exception as e:
            self.logger.error(f"Database connection failed: {e}")
            validations['database'] = False

        # Check processor initialization
        validations['processors'] = len(self.generator.processors) > 0

        # Check configuration
        validations['config_loaded'] = bool(self.config)
        validations['reference_set'] = bool(self._get_reference())

        # Log validation results
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
        """Context manager exit - cleanup resources"""
        self.clear_all_caches()
        if hasattr(self.db, 'close'):
            self.db.close()


# Convenience functions for common use cases

def create_service(config_path: Optional[str] = None) -> DomainSummaryService:
    """
    Create a domain summary service with default configuration.

    Args:
        config_path: Optional path to configuration file

    Returns:
        Configured DomainSummaryService instance
    """
    if config_path:
        context = ApplicationContext(config_path)
    else:
        # Use default or environment variable
        config_path = os.environ.get('ECOD_CONFIG_PATH', 'config/config.yml')
        context = ApplicationContext(config_path)

    return DomainSummaryService(context)


def process_single_protein(pdb_id: str, chain_id: str, job_dump_dir: str,
                          config_path: Optional[str] = None,
                          **options) -> DomainPartitionResult:
    """
    Convenience function to process a single protein.

    Args:
        pdb_id: PDB identifier
        chain_id: Chain identifier
        job_dump_dir: Directory containing job files
        config_path: Optional configuration file path
        **options: Processing options

    Returns:
        DomainPartitionResult
    """
    service = create_service(config_path)
    return service.process_protein(pdb_id, chain_id, job_dump_dir, **options)
