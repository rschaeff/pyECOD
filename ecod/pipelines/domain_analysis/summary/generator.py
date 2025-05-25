#!/usr/bin/env python3
"""
Core domain summary generation logic.

This module orchestrates the collection and processing of evidence from
multiple sources to generate comprehensive domain summaries.
"""

import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List, Optional, Set, Tuple, Union
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib

from ecod.db import DBManager
from ecod.models.pipeline import Evidence, DomainPartitionResult
from ecod.exceptions import PipelineError, ValidationError

from .file_locator import EvidenceFileLocator, FileType, SearchResult
from .processors.base import EvidenceProcessor, ProcessingResult, ProcessingStatus
from .processors.blast import create_blast_processor
from .processors.hhsearch import HHSearchEvidenceProcessor
from .processors.self_comparison import SelfComparisonProcessor
from .filters import EvidenceQualityFilter
from .models import (
    ProteinIdentifier, SequenceInfo, SummaryOptions,
    EvidenceCollection, EvidenceSummary, ProcessingStats
)


@dataclass
class GeneratorConfig:
    """Configuration for domain summary generator"""
    max_workers: int = 4
    process_timeout: int = 300
    cache_enabled: bool = True
    validate_inputs: bool = True
    parallel_evidence_collection: bool = True
    skip_on_error: bool = False
    collect_detailed_stats: bool = True


class DomainSummaryGenerator:
    """
    Generates domain summaries by orchestrating evidence collection and processing.

    This class coordinates:
    - File discovery through EvidenceFileLocator
    - Evidence extraction through processors
    - Quality filtering through EvidenceQualityFilter
    - Result aggregation into structured summaries
    """

    def __init__(self, config: Dict[str, Any], db_manager: DBManager,
                 generator_config: Optional[GeneratorConfig] = None):
        """
        Initialize generator with configuration.

        Args:
            config: Application configuration dictionary
            db_manager: Database manager instance
            generator_config: Generator-specific configuration
        """
        self.config = config
        self.db = db_manager
        self.generator_config = generator_config or GeneratorConfig()
        self.logger = logging.getLogger(__name__)

        # Initialize components (will be set per job)
        self.file_locator: Optional[EvidenceFileLocator] = None
        self.processors: Dict[str, EvidenceProcessor] = {}
        self.evidence_filter: Optional[EvidenceQualityFilter] = None

        # Statistics tracking
        self.stats = ProcessingStats()

        # Cache for processed results
        self._summary_cache: Dict[str, EvidenceSummary] = {}

    def initialize_for_job(self, job_dump_dir: Union[str, Path],
                          reference: Optional[str] = None) -> None:
        """
        Initialize components for a specific job directory.

        Args:
            job_dump_dir: Base directory for job files
            reference: Reference version (uses config default if not provided)
        """
        if reference is None:
            reference = self.config.get('reference', {}).get('current_version', 'develop291')

        self.job_dump_dir = Path(job_dump_dir)
        self.reference = reference

        # Initialize file locator
        self.file_locator = EvidenceFileLocator(
            db_manager=self.db,
            job_dump_dir=self.job_dump_dir,
            reference=reference
        )

        # Initialize processors
        self._initialize_processors()

        # Initialize evidence filter
        self.evidence_filter = EvidenceQualityFilter(
            min_confidence=self.config.get('evidence', {}).get('min_confidence', 0.3)
        )

        # Clear caches
        self._summary_cache.clear()

        self.logger.info(
            f"Initialized generator for job: {job_dump_dir} "
            f"(reference: {reference})"
        )

    def _initialize_processors(self) -> None:
        """Initialize evidence processors with configuration"""
        # Get processor configuration
        proc_config = self.config.get('processors', {})

        # Initialize BLAST processors
        blast_config = proc_config.get('blast', {})
        self.processors['chain_blast'] = create_blast_processor(
            blast_type='chain_blast',
            hsp_evalue_threshold=blast_config.get('hsp_evalue_threshold', 0.005),
            hit_coverage_threshold=blast_config.get('hit_coverage_threshold', 0.7),
            logger=self.logger
        )

        self.processors['domain_blast'] = create_blast_processor(
            blast_type='domain_blast',
            hsp_evalue_threshold=blast_config.get('hsp_evalue_threshold', 0.005),
            logger=self.logger
        )

        # Initialize HHSearch processor
        hhsearch_config = proc_config.get('hhsearch', {})
        self.processors['hhsearch'] = HHSearchEvidenceProcessor(
            probability_threshold=hhsearch_config.get('probability_threshold', 0.0),
            logger=self.logger
        )

        # Initialize self-comparison processor
        self.processors['self_comparison'] = SelfComparisonProcessor(
            logger=self.logger
        )

        self.logger.debug(f"Initialized {len(self.processors)} processors")

    def generate_summary(self, protein: ProteinIdentifier,
                        options: Optional[SummaryOptions] = None) -> EvidenceSummary:
        """
        Generate evidence summary for a protein.

        Args:
            protein: Protein identifier
            options: Processing options

        Returns:
            EvidenceSummary with collected and processed evidence
        """
        if not self.file_locator:
            raise RuntimeError("Generator not initialized. Call initialize_for_job first.")

        options = options or SummaryOptions()

        # Check cache
        cache_key = self._get_cache_key(protein, options)
        if self.generator_config.cache_enabled and cache_key in self._summary_cache:
            self.logger.debug(f"Using cached summary for {protein.source_id}")
            return self._summary_cache[cache_key]

        # Start timing
        start_time = datetime.now()

        # Create summary
        summary = EvidenceSummary(
            pdb_id=protein.pdb_id,
            chain_id=protein.chain_id,
            reference=protein.reference or self.reference
        )

        try:
            # Validate protein if configured
            if self.generator_config.validate_inputs:
                self._validate_protein(protein)

            # Check if already processed
            if self._already_processed(protein, options):
                summary.mark_as_skipped("Already processed")
                return summary

            # Get sequence information
            sequence_info = self._get_sequence_info(protein)
            summary.set_sequence_info(sequence_info)

            # Check for peptide
            if sequence_info and sequence_info.is_peptide():
                summary.mark_as_peptide()
                self.stats.peptides_processed += 1
                return self._finalize_summary(summary, start_time, cache_key)

            # Collect evidence from all sources
            evidence_collection = self._collect_all_evidence(protein, options)

            # Filter evidence
            if not options.skip_filtering:
                filtered_collection = self.evidence_filter.filter(evidence_collection)
                evidence_collection = filtered_collection

            # Process discontinuous evidence if needed
            if options.stitch_discontinuous:
                evidence_collection = self._process_discontinuous_evidence(
                    evidence_collection,
                    options.max_gap_for_stitching
                )

            # Add to summary
            summary.add_evidence_collection(evidence_collection)

            # Update statistics
            self.stats.proteins_processed += 1
            if summary.has_evidence():
                self.stats.proteins_with_evidence += 1

        except Exception as e:
            error_msg = f"Error generating summary for {protein.source_id}: {str(e)}"
            self.logger.error(error_msg, exc_info=True)
            summary.mark_as_error(error_msg)
            self.stats.errors += 1

            if not self.generator_config.skip_on_error:
                raise PipelineError(error_msg) from e

        return self._finalize_summary(summary, start_time, cache_key)

    def _finalize_summary(self, summary: EvidenceSummary, start_time: datetime,
                         cache_key: str) -> EvidenceSummary:
        """Finalize summary with timing and caching"""
        # Calculate processing time
        summary.processing_time = (datetime.now() - start_time).total_seconds()

        # Cache if enabled
        if self.generator_config.cache_enabled:
            self._summary_cache[cache_key] = summary

        # Log summary
        self.logger.info(
            f"Generated summary for {summary.protein_id}: "
            f"{summary.total_evidence_count} evidence items in {summary.processing_time:.2f}s"
        )

        return summary

    def _validate_protein(self, protein: ProteinIdentifier) -> None:
        """Validate protein identifier"""
        try:
            # ProteinIdentifier validates itself in __post_init__
            pass
        except Exception as e:
            raise ValidationError(f"Invalid protein identifier: {e}")

    def _already_processed(self, protein: ProteinIdentifier,
                          options: SummaryOptions) -> bool:
        """Check if protein was already processed"""
        if options.force_overwrite:
            return False

        # Check for existing summary file
        suffix = ".blast_only" if options.blast_only else ""
        summary_filename = (
            f"{protein.source_id}.{self.reference}.domain_summary{suffix}.xml"
        )

        result = self.file_locator.find_file(
            protein.pdb_id,
            protein.chain_id,
            FileType.DOMAIN_SUMMARY
        )

        # Check if the specific filename exists
        if result.found:
            for location in result.locations:
                if location.path.name == summary_filename and location.is_valid:
                    return True

        return False

    def _get_sequence_info(self, protein: ProteinIdentifier) -> Optional[SequenceInfo]:
        """Get sequence information for protein"""
        # Find FASTA file
        fasta_path = self.file_locator.find_sequence_file(
            protein.pdb_id,
            protein.chain_id
        )

        if not fasta_path:
            # Try to get from database
            sequence = self._get_sequence_from_db(protein)
            if sequence:
                return SequenceInfo(
                    sequence=sequence,
                    length=len(sequence),
                    md5=hashlib.md5(sequence.encode()).hexdigest()
                )

            self.logger.warning(f"No sequence found for {protein.source_id}")
            return None

        try:
            return SequenceInfo.from_fasta(fasta_path)
        except Exception as e:
            self.logger.error(f"Error reading sequence from {fasta_path}: {e}")
            return None

    def _get_sequence_from_db(self, protein: ProteinIdentifier) -> Optional[str]:
        """Get sequence from database"""
        query = """
        SELECT ps.sequence
        FROM ecod_schema.protein p
        JOIN ecod_schema.protein_sequence ps ON p.id = ps.protein_id
        WHERE p.pdb_id = %s AND p.chain_id = %s
        LIMIT 1
        """

        try:
            rows = self.db.execute_query(query, (protein.pdb_id, protein.chain_id))
            if rows and rows[0][0]:
                return rows[0][0]
        except Exception as e:
            self.logger.error(f"Error getting sequence from database: {e}")

        return None

    def _collect_all_evidence(self, protein: ProteinIdentifier,
                             options: SummaryOptions) -> EvidenceCollection:
        """Collect evidence from all configured sources"""
        collection = EvidenceCollection()

        # Define evidence sources based on options
        sources = self._get_evidence_sources(options)

        if self.generator_config.parallel_evidence_collection and len(sources) > 1:
            # Collect evidence in parallel
            collection = self._collect_evidence_parallel(protein, sources)
        else:
            # Collect evidence sequentially
            collection = self._collect_evidence_sequential(protein, sources)

        # Log collection summary
        self.logger.info(
            f"Collected evidence for {protein.source_id}: "
            f"{collection.get_summary_string()}"
        )

        return collection

    def _get_evidence_sources(self, options: SummaryOptions) -> List[str]:
        """Get list of evidence sources based on options"""
        sources = ['self_comparison', 'chain_blast', 'domain_blast']

        if not options.blast_only:
            sources.append('hhsearch')

        if options.additional_sources:
            sources.extend(options.additional_sources)

        return sources

    def _collect_evidence_sequential(self, protein: ProteinIdentifier,
                                    sources: List[str]) -> EvidenceCollection:
        """Collect evidence sequentially from each source"""
        collection = EvidenceCollection()

        for source in sources:
            try:
                evidence = self._collect_evidence_from_source(protein, source)
                collection.add_source(source, evidence)

                # Update statistics
                if evidence:
                    self.stats.evidence_by_source[source] = \
                        self.stats.evidence_by_source.get(source, 0) + len(evidence)

            except Exception as e:
                self.logger.warning(
                    f"Error collecting {source} evidence for {protein.source_id}: {e}"
                )
                collection.add_error(source, str(e))

        return collection

    def _collect_evidence_parallel(self, protein: ProteinIdentifier,
                                  sources: List[str]) -> EvidenceCollection:
        """Collect evidence in parallel from multiple sources"""
        collection = EvidenceCollection()

        with ThreadPoolExecutor(max_workers=self.generator_config.max_workers) as executor:
            # Submit tasks
            future_to_source = {
                executor.submit(
                    self._collect_evidence_from_source, protein, source
                ): source
                for source in sources
            }

            # Collect results
            for future in as_completed(future_to_source):
                source = future_to_source[future]

                try:
                    evidence = future.result(timeout=self.generator_config.process_timeout)
                    collection.add_source(source, evidence)

                    # Update statistics
                    if evidence:
                        self.stats.evidence_by_source[source] = \
                            self.stats.evidence_by_source.get(source, 0) + len(evidence)

                except Exception as e:
                    self.logger.warning(
                        f"Error collecting {source} evidence for {protein.source_id}: {e}"
                    )
                    collection.add_error(source, str(e))

        return collection

    def _collect_evidence_from_source(self, protein: ProteinIdentifier,
                                     source: str) -> List[Evidence]:
        """Collect evidence from a specific source"""
        # Map source to file type
        file_type_map = {
            'chain_blast': FileType.CHAIN_BLAST,
            'domain_blast': FileType.DOMAIN_BLAST,
            'hhsearch': FileType.HHSEARCH,
            'self_comparison': FileType.SELF_COMPARISON
        }

        file_type = file_type_map.get(source)
        if not file_type:
            self.logger.warning(f"Unknown evidence source: {source}")
            return []

        # Find file
        search_result = self.file_locator.find_file(
            protein.pdb_id,
            protein.chain_id,
            file_type
        )

        if not search_result.found:
            self.logger.debug(
                f"No {source} file found for {protein.source_id}"
            )
            return []

        # Get processor
        processor = self.processors.get(source)
        if not processor:
            self.logger.warning(f"No processor configured for {source}")
            return []

        # Process file
        file_path = search_result.best_location.path
        self.logger.debug(f"Processing {source} file: {file_path}")

        processing_result = processor.process(file_path)

        if processing_result.success:
            return processing_result.evidence
        else:
            self.logger.warning(
                f"Failed to process {source} file: {processing_result.error}"
            )
            return []

    def _process_discontinuous_evidence(self, collection: EvidenceCollection,
                                       max_gap: int) -> EvidenceCollection:
        """Process and stitch discontinuous evidence"""
        processed = EvidenceCollection()

        for source, evidence_list in collection.get_all_evidence().items():
            if source in ['domain_blast', 'hhsearch']:
                # These sources may have discontinuous evidence
                stitched = self._stitch_evidence(evidence_list, max_gap)
                processed.add_source(source, stitched)
            else:
                # Keep as-is
                processed.add_source(source, evidence_list)

        # Copy errors
        for source, error in collection.get_errors().items():
            processed.add_error(source, error)

        return processed

    def _stitch_evidence(self, evidence_list: List[Evidence],
                        max_gap: int) -> List[Evidence]:
        """Stitch together evidence that may represent discontinuous domains"""
        if len(evidence_list) <= 1:
            return evidence_list

        # Group by source_id/domain_id
        groups: Dict[str, List[Evidence]] = {}

        for evidence in evidence_list:
            key = evidence.domain_id or evidence.source_id
            if key not in groups:
                groups[key] = []
            groups[key].append(evidence)

        # Process each group
        stitched = []

        for key, group in groups.items():
            if len(group) == 1:
                stitched.extend(group)
            else:
                # Try to stitch
                stitched_group = self._try_stitch_group(group, max_gap)
                stitched.extend(stitched_group)

        return stitched

    def _try_stitch_group(self, evidence_group: List[Evidence],
                         max_gap: int) -> List[Evidence]:
        """Try to stitch a group of evidence"""
        # Sort by query start position
        sorted_evidence = sorted(
            evidence_group,
            key=lambda e: self._get_start_position(e.query_range)
        )

        stitched = []
        current_set = [sorted_evidence[0]]

        for i in range(1, len(sorted_evidence)):
            evidence = sorted_evidence[i]

            # Check if can be stitched with current set
            can_stitch = False

            for existing in current_set:
                gap = self._calculate_evidence_gap(existing, evidence)
                if gap <= max_gap:
                    can_stitch = True
                    break

            if can_stitch:
                current_set.append(evidence)
            else:
                # Finalize current set and start new one
                if len(current_set) > 1:
                    merged = self._merge_evidence_set(current_set)
                    stitched.append(merged)
                else:
                    stitched.extend(current_set)

                current_set = [evidence]

        # Handle final set
        if len(current_set) > 1:
            merged = self._merge_evidence_set(current_set)
            stitched.append(merged)
        else:
            stitched.extend(current_set)

        return stitched

    def _get_start_position(self, range_str: str) -> int:
        """Extract start position from range string"""
        if not range_str:
            return 0

        # Handle single range
        if "," not in range_str:
            if "-" in range_str:
                return int(range_str.split("-")[0])
            return 0

        # Handle multiple ranges - get first start
        first_range = range_str.split(",")[0]
        if "-" in first_range:
            return int(first_range.split("-")[0])

        return 0

    def _calculate_evidence_gap(self, ev1: Evidence, ev2: Evidence) -> int:
        """Calculate gap between two evidence items"""
        # Get end of ev1
        ranges1 = ev1.query_range.split(",")
        last_range1 = ranges1[-1]
        end1 = int(last_range1.split("-")[1]) if "-" in last_range1 else 0

        # Get start of ev2
        start2 = self._get_start_position(ev2.query_range)

        return start2 - end1 - 1

    def _merge_evidence_set(self, evidence_set: List[Evidence]) -> Evidence:
        """Merge a set of evidence into a single stitched evidence"""
        # Use first evidence as template
        first = evidence_set[0]

        # Combine ranges
        all_query_ranges = []
        all_hit_ranges = []

        for ev in evidence_set:
            all_query_ranges.append(ev.query_range)
            all_hit_ranges.append(ev.hit_range)

        # Create merged evidence
        merged = Evidence(
            type=first.type,
            source_id=first.source_id,
            domain_id=first.domain_id,
            query_range=",".join(all_query_ranges),
            hit_range=",".join(all_hit_ranges),
            confidence=None,  # Will be recalculated
            t_group=first.t_group,
            h_group=first.h_group,
            x_group=first.x_group,
            a_group=first.a_group,
            evalue=min(ev.evalue for ev in evidence_set if ev.evalue is not None),
            probability=max(ev.probability for ev in evidence_set if ev.probability is not None),
            score=sum(ev.score for ev in evidence_set if ev.score is not None),
            hsp_count=sum(ev.hsp_count or 1 for ev in evidence_set),
            extra_attributes={
                "stitched": True,
                "segment_count": len(evidence_set),
                "original_evidence_count": len(evidence_set)
            }
        )

        return merged

    def _get_cache_key(self, protein: ProteinIdentifier,
                      options: SummaryOptions) -> str:
        """Generate cache key for a protein and options"""
        options_str = f"blast_only={options.blast_only},force={options.force_overwrite}"
        return f"{protein.source_id}:{self.reference}:{options_str}"

    def get_statistics(self) -> Dict[str, Any]:
        """Get processing statistics"""
        return {
            "proteins_processed": self.stats.proteins_processed,
            "proteins_with_evidence": self.stats.proteins_with_evidence,
            "peptides_processed": self.stats.peptides_processed,
            "errors": self.stats.errors,
            "evidence_by_source": dict(self.stats.evidence_by_source),
            "cache_size": len(self._summary_cache),
            "processors_loaded": list(self.processors.keys())
        }

    def clear_caches(self) -> None:
        """Clear all caches"""
        self._summary_cache.clear()

        if self.file_locator:
            self.file_locator.clear_cache()

        for processor in self.processors.values():
            processor.clear_cache()

        self.logger.info("Cleared all generator caches")
