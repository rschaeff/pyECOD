#!/usr/bin/env python3
"""
Unified domain partitioning models - Synthesis of all model versions.

This module consolidates all partition-related models into a single,
comprehensive model hierarchy that eliminates conflicts.
"""

from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum, auto
from typing import Dict, Any, List, Optional, Tuple, Set, Union
from pathlib import Path

# Import the gold standard pipeline models
from ecod.models.pipeline.evidence import Evidence
from ecod.models.pipeline.domain import DomainModel
from ecod.models.pipeline.partition import DomainPartitionResult


# =============================================================================
# CORE ENUMS - Single source of truth
# =============================================================================

class PartitionStage(Enum):
    """Unified stages in the partition processing pipeline"""
    INITIALIZING = auto()
    LOADING_SUMMARY = auto()
    VALIDATING_EVIDENCE = auto()
    EXTRACTING_EVIDENCE = auto()
    IDENTIFYING_BOUNDARIES = auto()
    RESOLVING_OVERLAPS = auto()
    ASSIGNING_CLASSIFICATIONS = auto()
    SAVING_RESULTS = auto()
    COMPLETE = auto()
    FAILED = auto()


class ValidationLevel(Enum):
    """Evidence validation levels"""
    LENIENT = "lenient"      # Allow minor issues, focus on processing
    NORMAL = "normal"        # Standard validation
    STRICT = "strict"        # Strict validation, fail on any issues


class ProcessingMode(Enum):
    """Processing modes for domain partitioning"""
    STANDARD = "standard"          # Normal processing
    BLAST_ONLY = "blast_only"      # Use only BLAST evidence
    REPRESENTATIVES = "representatives"  # Process only representative proteins
    REPROCESS = "reprocess"        # Force reprocessing of existing results


# =============================================================================
# UNIFIED CONFIGURATION MODEL
# =============================================================================

@dataclass
class PartitionOptions:
    """
    Unified configuration options for domain partitioning.

    Consolidates all configuration from both model versions into a
    comprehensive, backwards-compatible options class.
    """

    # ========== PROCESSING MODE AND BASIC SETTINGS ==========
    mode: ProcessingMode = ProcessingMode.STANDARD
    validation_level: ValidationLevel = ValidationLevel.NORMAL
    blast_only: bool = False
    representatives_only: bool = False
    force_overwrite: bool = False

    # ========== DOMAIN SIZE CONSTRAINTS ==========
    min_domain_size: int = 20
    max_domain_size: Optional[int] = 2000

    # ========== OVERLAP AND BOUNDARY HANDLING ==========
    overlap_threshold: float = 0.3  # Maximum allowed overlap (0-1)
    resolve_overlaps: bool = True
    merge_gap_tolerance: int = 20

    # ========== EVIDENCE PROCESSING ==========
    use_chain_blast: bool = True
    use_domain_blast: bool = True
    use_hhsearch: bool = True
    min_evidence_confidence: float = 0.0
    require_evidence: bool = True

    # ========== VALIDATION SETTINGS ==========
    validate_coordinates: bool = True
    validate_classifications: bool = True

    # ========== CLASSIFICATION SETTINGS ==========
    prefer_hhsearch_classification: bool = True
    require_full_classification: bool = False
    classification_confidence_weight: Dict[str, float] = field(default_factory=lambda: {
        "hhsearch": 3.0,
        "domain_blast": 2.5,
        "chain_blast": 2.0,
        "blast": 1.5
    })

    # ========== REFERENCE COVERAGE SETTINGS ==========
    min_reference_coverage: float = 0.7  # Minimum 70% coverage required
    strict_reference_coverage: float = 0.9  # Above this, trust boundaries exactly
    partial_coverage_threshold: float = 0.3  # Below this, reject evidence
    extend_to_reference_size: bool = True  # Try to extend to match reference
    reference_size_tolerance: float = 0.15  # Allow 15% size difference
    max_extension_length: int = 50  # Don't extend more than 50 residues
    use_ungapped_coverage: bool = True  # Calculate coverage excluding gaps
    combine_partial_evidence: bool = True  # Try to combine partial alignments

    # ========== PERFORMANCE SETTINGS ==========
    use_cache: bool = True
    parallel_processing: bool = False
    max_workers: int = 4
    batch_size: int = 100

    # ========== OUTPUT SETTINGS ==========
    save_intermediate: bool = False
    include_evidence_in_output: bool = True
    include_metadata: bool = True
    include_detailed_stats: bool = False
    output_format: str = "xml"

    def validate(self) -> None:
        """Comprehensive validation of all options"""
        errors = []

        # Domain size validation
        if self.min_domain_size < 1:
            errors.append("min_domain_size must be at least 1")

        if self.max_domain_size and self.max_domain_size < self.min_domain_size:
            errors.append("max_domain_size must be >= min_domain_size")

        # Threshold validation
        if not 0.0 <= self.overlap_threshold <= 1.0:
            errors.append("overlap_threshold must be between 0.0 and 1.0")

        if not 0.0 <= self.min_evidence_confidence <= 1.0:
            errors.append("min_evidence_confidence must be between 0.0 and 1.0")

        # Coverage validation
        if not 0.0 <= self.min_reference_coverage <= 1.0:
            errors.append("min_reference_coverage must be between 0.0 and 1.0")

        if not 0.0 <= self.strict_reference_coverage <= 1.0:
            errors.append("strict_reference_coverage must be between 0.0 and 1.0")

        if self.strict_reference_coverage < self.min_reference_coverage:
            errors.append("strict_reference_coverage must be >= min_reference_coverage")

        # Performance validation
        if self.max_workers < 1:
            errors.append("max_workers must be positive")

        if self.batch_size < 1:
            errors.append("batch_size must be positive")

        if errors:
            raise ValueError(f"Invalid partition options: {'; '.join(errors)}")

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization"""
        return {
            # Processing mode
            'mode': self.mode.value,
            'validation_level': self.validation_level.value,
            'blast_only': self.blast_only,
            'representatives_only': self.representatives_only,
            'force_overwrite': self.force_overwrite,

            # Domain constraints
            'min_domain_size': self.min_domain_size,
            'max_domain_size': self.max_domain_size,

            # Overlap handling
            'overlap_threshold': self.overlap_threshold,
            'resolve_overlaps': self.resolve_overlaps,
            'merge_gap_tolerance': self.merge_gap_tolerance,

            # Evidence processing
            'use_chain_blast': self.use_chain_blast,
            'use_domain_blast': self.use_domain_blast,
            'use_hhsearch': self.use_hhsearch,
            'min_evidence_confidence': self.min_evidence_confidence,
            'require_evidence': self.require_evidence,

            # Validation
            'validate_coordinates': self.validate_coordinates,
            'validate_classifications': self.validate_classifications,

            # Classification
            'prefer_hhsearch_classification': self.prefer_hhsearch_classification,
            'require_full_classification': self.require_full_classification,
            'classification_confidence_weight': self.classification_confidence_weight,

            # Reference coverage
            'min_reference_coverage': self.min_reference_coverage,
            'strict_reference_coverage': self.strict_reference_coverage,
            'partial_coverage_threshold': self.partial_coverage_threshold,
            'extend_to_reference_size': self.extend_to_reference_size,
            'reference_size_tolerance': self.reference_size_tolerance,
            'max_extension_length': self.max_extension_length,
            'use_ungapped_coverage': self.use_ungapped_coverage,
            'combine_partial_evidence': self.combine_partial_evidence,

            # Performance
            'use_cache': self.use_cache,
            'parallel_processing': self.parallel_processing,
            'max_workers': self.max_workers,
            'batch_size': self.batch_size,

            # Output
            'save_intermediate': self.save_intermediate,
            'include_evidence_in_output': self.include_evidence_in_output,
            'include_metadata': self.include_metadata,
            'include_detailed_stats': self.include_detailed_stats,
            'output_format': self.output_format
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'PartitionOptions':
        """Create from dictionary with enum conversion"""
        # Handle enum conversions
        if 'mode' in data and isinstance(data['mode'], str):
            data['mode'] = ProcessingMode(data['mode'])
        if 'validation_level' in data and isinstance(data['validation_level'], str):
            data['validation_level'] = ValidationLevel(data['validation_level'])

        return cls(**data)

    @classmethod
    def create_blast_only(cls) -> 'PartitionOptions':
        """Factory method for BLAST-only processing"""
        return cls(
            mode=ProcessingMode.BLAST_ONLY,
            blast_only=True,
            use_hhsearch=False,
            min_evidence_confidence=0.1,
            validation_level=ValidationLevel.LENIENT
        )

    @classmethod
    def create_strict(cls) -> 'PartitionOptions':
        """Factory method for strict processing"""
        return cls(
            validation_level=ValidationLevel.STRICT,
            require_evidence=True,
            require_full_classification=True,
            min_evidence_confidence=0.8,
            min_reference_coverage=0.9
        )

    # Add this method to the PartitionOptions class in ecod/pipelines/domain_analysis/partition/models.py

    def __post_init__(self):
        """Validate options after initialization"""
        # Validate enum values
        if isinstance(self.validation_level, str):
            try:
                self.validation_level = ValidationLevel(self.validation_level)
            except ValueError:
                raise ValueError(f"Invalid validation_level: {self.validation_level}. "
                               f"Must be one of: {', '.join([v.value for v in ValidationLevel])}")

        if isinstance(self.mode, str):
            try:
                self.mode = ProcessingMode(self.mode)
            except ValueError:
                raise ValueError(f"Invalid mode: {self.mode}. "
                               f"Must be one of: {', '.join([m.value for m in ProcessingMode])}")

        # Validate all other constraints
        self.validate()


# =============================================================================
# PROCESSING RESULT MODELS
# =============================================================================

@dataclass
class ValidationResult:
    """Unified validation result model"""
    is_valid: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    context: str = ""
    validation_level: ValidationLevel = ValidationLevel.NORMAL
    details: Dict[str, Any] = field(default_factory=dict)

    def add_error(self, message: str) -> None:
        """Add an error message"""
        self.errors.append(message)
        self.is_valid = False

    def add_warning(self, message: str) -> None:
        """Add a warning message"""
        self.warnings.append(message)

    def merge(self, other: 'ValidationResult') -> None:
        """Merge another validation result into this one"""
        self.is_valid = self.is_valid and other.is_valid
        self.errors.extend(other.errors)
        self.warnings.extend(other.warnings)
        self.details.update(other.details)

    def get_summary(self) -> str:
        """Get a summary of validation results"""
        parts = []
        if self.context:
            parts.append(f"Context: {self.context}")

        if self.is_valid:
            parts.append("VALID")
        else:
            parts.append("INVALID")

        if self.errors:
            parts.append(f"Errors: {len(self.errors)}")
        if self.warnings:
            parts.append(f"Warnings: {len(self.warnings)}")

        return " | ".join(parts)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'is_valid': self.is_valid,
            'errors': self.errors,
            'warnings': self.warnings,
            'context': self.context,
            'validation_level': self.validation_level.value,
            'details': self.details,
            'summary': self.get_summary()
        }


@dataclass
class EvidenceGroup:
    """Unified evidence group model"""
    group_id: str = ""
    evidence_items: List[Evidence] = field(default_factory=list)
    consensus_start: Optional[int] = None
    consensus_end: Optional[int] = None
    consensus_confidence: float = 0.0
    position_key: int = 0
    group_type: str = "unknown"  # "blast", "hhsearch", "mixed"

    def __init__(self, group_id: str = "", evidence_items: Optional[List[Evidence]] = None):
        """Initialize evidence group"""
        self.group_id = group_id
        self.evidence_items = evidence_items or []
        self.consensus_start = None
        self.consensus_end = None
        self.consensus_confidence = 0.0
        self.position_key = 0
        self.group_type = "unknown"

        # Update stats if we have evidence
        if self.evidence_items:
            self._update_group_stats()

    def add_evidence(self, evidence: Evidence) -> None:
        """Add evidence to the group"""
        self.evidence_items.append(evidence)
        self._update_group_stats()

    def get_best_evidence(self) -> Optional[Evidence]:
        """Get the highest confidence evidence"""
        if not self.evidence_items:
            return None

        return max(self.evidence_items,
                  key=lambda e: e.confidence if e.confidence is not None else 0.0)

    def _update_group_stats(self) -> None:
        """Update group statistics"""
        if not self.evidence_items:
            return

        # Update confidence score (average)
        confidences = [e.confidence for e in self.evidence_items if e.confidence is not None]
        if confidences:
            self.consensus_confidence = sum(confidences) / len(confidences)

        # Update group type
        types = set(e.type for e in self.evidence_items if e.type)
        if len(types) == 1:
            self.group_type = types.pop()
        else:
            self.group_type = "mixed"

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'group_id': self.group_id,
            'evidence_count': len(self.evidence_items),
            'consensus_start': self.consensus_start,
            'consensus_end': self.consensus_end,
            'consensus_confidence': self.consensus_confidence,
            'group_type': self.group_type
        }


    def _update_group_stats(self) -> None:
        """Update group statistics including consensus start/end positions"""
        if not self.evidence_items:
            return

        # Update confidence score (average)
        confidences = [e.confidence for e in self.evidence_items if e.confidence is not None]
        if confidences:
            self.consensus_confidence = sum(confidences) / len(confidences)

        # Update group type
        types = set(e.type for e in self.evidence_items if e.type)
        if len(types) == 1:
            self.group_type = types.pop()
        else:
            self.group_type = "mixed"

        # CRITICAL FIX: Calculate consensus start/end positions
        self._calculate_consensus_positions()

    def _calculate_consensus_positions(self) -> None:
        """Calculate consensus start and end positions from evidence ranges"""
        if not self.evidence_items:
            self.consensus_start = None
            self.consensus_end = None
            return

        valid_ranges = []

        # Parse ranges from evidence items
        for evidence in self.evidence_items:
            if not evidence.query_range:
                continue

            try:
                # Parse range format like "10-185" or "180-355"
                range_str = evidence.query_range.strip()
                if '-' in range_str:
                    parts = range_str.split('-')
                    if len(parts) == 2:
                        start = int(parts[0].strip())
                        end = int(parts[1].strip())
                        if start > 0 and end > 0 and start <= end:
                            valid_ranges.append((start, end))
            except (ValueError, IndexError):
                # Skip malformed ranges
                continue

        if not valid_ranges:
            self.consensus_start = None
            self.consensus_end = None
            return

        # Calculate consensus using union approach (covers all evidence)
        # This gives maximum coverage which is appropriate for domain finding
        all_starts = [start for start, end in valid_ranges]
        all_ends = [end for start, end in valid_ranges]

        self.consensus_start = min(all_starts)
        self.consensus_end = max(all_ends)

        # Update position key for grouping (use start position)
        self.position_key = self.consensus_start

    def get_consensus_range(self) -> Optional[str]:
        """Get consensus range as string"""
        if self.consensus_start is not None and self.consensus_end is not None:
            return f"{self.consensus_start}-{self.consensus_end}"
        return None

    def get_consensus_size(self) -> int:
        """Get consensus domain size"""
        if self.consensus_start is not None and self.consensus_end is not None:
            return self.consensus_end - self.consensus_start + 1
        return 0


@dataclass
class DomainCandidate:
    """Unified domain candidate model"""
    start: int
    end: int
    evidence_group: EvidenceGroup
    source: str = ""
    confidence: float = 0.0
    protected: bool = False

    @property
    def size(self) -> int:
        """Get domain size"""
        return self.end - self.start + 1

    @property
    def range(self) -> str:
        """Get range string"""
        return f"{self.start}-{self.end}"

    def to_domain_model(self, domain_id: str) -> DomainModel:
        """Convert to DomainModel"""
        best_evidence = self.evidence_group.get_best_evidence()

        return DomainModel(
            id=domain_id,
            start=self.start,
            end=self.end,
            range=self.range,
            source=self.source or (best_evidence.type if best_evidence else "unknown"),
            confidence=self.confidence,
            source_id=best_evidence.source_id if best_evidence else "",
            evidence=self.evidence_group.evidence_items,
            protected=self.protected
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'start': self.start,
            'end': self.end,
            'size': self.size,
            'range': self.range,
            'source': self.source,
            'confidence': self.confidence,
            'protected': self.protected,
            'evidence_group': self.evidence_group.to_dict()
        }


@dataclass
class PartitionContext:
    """Unified partition context model"""
    pdb_id: str
    chain_id: str
    reference: str
    sequence_length: int = 0
    output_dir: Path = field(default_factory=Path)
    process_id: Optional[int] = None
    batch_id: Optional[int] = None

    # Timing information
    start_time: datetime = field(default_factory=datetime.now)
    stage_times: Dict[PartitionStage, float] = field(default_factory=dict)

    # Processing state
    current_stage: PartitionStage = PartitionStage.INITIALIZING

    @property
    def protein_id(self) -> str:
        """Get protein identifier"""
        return f"{self.pdb_id}_{self.chain_id}"

    @property
    def output_file(self) -> Path:
        """Get expected output file path"""
        return self.output_dir / "domains" / f"{self.protein_id}.{self.reference}.domains.xml"

    def record_stage_time(self, stage: PartitionStage) -> None:
        """Record completion time for a stage"""
        elapsed = (datetime.now() - self.start_time).total_seconds()
        self.stage_times[stage] = elapsed
        self.current_stage = stage

    def get_total_time(self) -> float:
        """Get total processing time"""
        return (datetime.now() - self.start_time).total_seconds()

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'pdb_id': self.pdb_id,
            'chain_id': self.chain_id,
            'protein_id': self.protein_id,
            'reference': self.reference,
            'sequence_length': self.sequence_length,
            'output_dir': str(self.output_dir),
            'process_id': self.process_id,
            'batch_id': self.batch_id,
            'current_stage': self.current_stage.name,
            'total_time': self.get_total_time()
        }


@dataclass
class BatchPartitionResults:
    """Unified batch partition results model"""
    total: int = 0
    success_count: int = 0
    failure_count: int = 0
    skipped_count: int = 0

    results: List[DomainPartitionResult] = field(default_factory=list)
    failures: List[Tuple[str, str, str]] = field(default_factory=list)  # (pdb_id, chain_id, error)

    # Statistics
    proteins_with_domains: int = 0
    total_domains_found: int = 0
    peptides_found: int = 0
    unclassified_proteins: int = 0

    # Timing
    start_time: datetime = field(default_factory=datetime.now)
    end_time: Optional[datetime] = None

    def add_result(self, result: DomainPartitionResult) -> None:
        """Add a partition result"""
        self.results.append(result)
        self.total += 1

        if result.success:
            self.success_count += 1

            if result.is_peptide:
                self.peptides_found += 1
            elif result.domains:
                self.proteins_with_domains += 1
                self.total_domains_found += len(result.domains)
            else:
                self.unclassified_proteins += 1
        else:
            self.failure_count += 1
            self.failures.append((result.pdb_id, result.chain_id, result.error or "Unknown error"))

    def add_skipped(self, pdb_id: str, chain_id: str, reason: str) -> None:
        """Record a skipped protein"""
        self.skipped_count += 1
        self.total += 1

    def finalize(self) -> None:
        """Finalize results and calculate statistics"""
        self.end_time = datetime.now()

    @property
    def success_rate(self) -> float:
        """Calculate success rate"""
        if self.total == 0:
            return 0.0
        return (self.success_count / self.total) * 100.0

    @property
    def processing_time(self) -> float:
        """Get total processing time in seconds"""
        if self.end_time is None:
            return (datetime.now() - self.start_time).total_seconds()
        return (self.end_time - self.start_time).total_seconds()

    def get_summary(self) -> Dict[str, Any]:
        """Get summary statistics"""
        return {
            'total_proteins': self.total,
            'successful': self.success_count,
            'failed': self.failure_count,
            'skipped': self.skipped_count,
            'success_rate': self.success_rate,
            'proteins_with_domains': self.proteins_with_domains,
            'total_domains': self.total_domains_found,
            'peptides': self.peptides_found,
            'unclassified': self.unclassified_proteins,
            'processing_time': self.processing_time,
            'average_time_per_protein': self.processing_time / self.total if self.total > 0 else 0
        }

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return self.get_summary()


# =============================================================================
# CACHE AND PERFORMANCE MODELS
# =============================================================================

@dataclass
class ClassificationCache:
    """Unified classification cache model"""
    domain_id_cache: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    chain_domain_cache: Dict[str, List[Dict[str, Any]]] = field(default_factory=dict)
    cache_hits: int = 0
    cache_misses: int = 0

    def get_domain_classification(self, domain_id: str) -> Optional[Dict[str, Any]]:
        """Get cached domain classification"""
        if domain_id in self.domain_id_cache:
            self.cache_hits += 1
            return self.domain_id_cache[domain_id]

        self.cache_misses += 1
        return None

    def set_domain_classification(self, domain_id: str, classification: Dict[str, Any]) -> None:
        """Cache domain classification"""
        self.domain_id_cache[domain_id] = classification

    def get_chain_domains(self, source_id: str) -> Optional[List[Dict[str, Any]]]:
        """Get cached chain domains"""
        if source_id in self.chain_domain_cache:
            self.cache_hits += 1
            return self.chain_domain_cache[source_id]

        self.cache_misses += 1
        return None

    def set_chain_domains(self, source_id: str, domains: List[Dict[str, Any]]) -> None:
        """Cache chain domains"""
        self.chain_domain_cache[source_id] = domains

    def clear(self) -> None:
        """Clear all caches"""
        self.domain_id_cache.clear()
        self.chain_domain_cache.clear()
        self.cache_hits = 0
        self.cache_misses = 0

    def get_stats(self) -> Dict[str, Any]:
        """Get cache statistics"""
        total_requests = self.cache_hits + self.cache_misses
        hit_rate = (self.cache_hits / total_requests * 100.0) if total_requests > 0 else 0.0

        return {
            'cache_hits': self.cache_hits,
            'cache_misses': self.cache_misses,
            'hit_rate': hit_rate,
            'domain_classifications_cached': len(self.domain_id_cache),
            'chain_domains_cached': len(self.chain_domain_cache)
        }


# =============================================================================
# REFERENCE COVERAGE MODELS
# =============================================================================

@dataclass
class ReferenceInfo:
    """Information about a reference domain"""
    domain_id: str
    domain_range: str
    domain_length: int
    pdb_id: str
    chain_id: str
    t_group: Optional[str] = None
    is_discontinuous: bool = False
    discontinuous_ranges: List[Tuple[int, int]] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'domain_id': self.domain_id,
            'domain_range': self.domain_range,
            'domain_length': self.domain_length,
            'pdb_id': self.pdb_id,
            'chain_id': self.chain_id,
            't_group': self.t_group,
            'is_discontinuous': self.is_discontinuous,
            'discontinuous_ranges': self.discontinuous_ranges
        }


@dataclass
class EvidenceWithCoverage(Evidence):
    """Extended evidence with coverage information"""
    reference_coverage: float = 0.0
    reference_info: Optional[ReferenceInfo] = None
    hit_length: int = 0
    alignment_gaps: int = 0
    coverage_warning: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary with coverage info"""
        result = super().to_dict()
        result.update({
            'reference_coverage': self.reference_coverage,
            'hit_length': self.hit_length,
            'alignment_gaps': self.alignment_gaps,
            'coverage_warning': self.coverage_warning
        })
        if self.reference_info:
            result['reference_info'] = self.reference_info.to_dict()
        return result


# =============================================================================
# EXPORT ALL MODELS
# =============================================================================

ProcessingStats = BatchPartitionResults


__all__ = [
    # Enums
    'PartitionStage',
    'ValidationLevel',
    'ProcessingMode',

    # Core models
    'PartitionOptions',
    'ValidationResult',
    'EvidenceGroup',
    'DomainCandidate',
    'PartitionContext',
    'BatchPartitionResults',

    # Cache models
    'ClassificationCache',

    # Coverage models
    'ReferenceInfo',
    'EvidenceWithCoverage',

    # Pipeline models (re-exported)
    'Evidence',
    'DomainModel',
    'DomainPartitionResult',

    #Aliases
    'ProcessingStats'
]
