#!/usr/bin/env python3
"""
Data models for domain partitioning.

This module contains all data models, configuration classes, and options
used by the domain partition service.
"""

from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum, auto
from typing import Dict, Any, List, Optional, Tuple, Set
from pathlib import Path

from ecod.models.pipeline.evidence import Evidence
from ecod.models.pipeline.domain import DomainModel
from ecod.models.pipeline.partition import DomainPartitionResult


class PartitionStage(Enum):
    """Stages in the partition processing pipeline"""
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
    STRICT = "strict"      # Reject any invalid evidence
    NORMAL = "normal"      # Log warnings but continue
    LENIENT = "lenient"    # Accept most evidence with minimal validation


class ProcessingMode(Enum):
    """Processing modes for domain partitioning"""
    STANDARD = "standard"          # Normal processing
    BLAST_ONLY = "blast_only"      # Use only BLAST evidence
    REPRESENTATIVES = "representatives"  # Process only representative proteins
    REPROCESS = "reprocess"        # Force reprocessing of existing results


@dataclass
class PartitionOptions:
    """Configuration options for domain partitioning"""

    # Processing mode
    mode: ProcessingMode = ProcessingMode.STANDARD
    blast_only: bool = False
    representatives_only: bool = False
    force_overwrite: bool = False

    # Validation settings
    validation_level: ValidationLevel = ValidationLevel.NORMAL
    validate_coordinates: bool = True
    validate_classifications: bool = True
    require_evidence: bool = True

    # Evidence processing
    use_chain_blast: bool = True
    use_domain_blast: bool = True
    use_hhsearch: bool = True
    min_evidence_confidence: float = 0.0

    # Domain boundary identification
    min_domain_size: int = 20
    max_domain_size: Optional[int] = None
    merge_gap_tolerance: int = 20
    overlap_threshold: float = 0.3

    # Classification settings
    prefer_hhsearch_classification: bool = True
    require_full_classification: bool = False
    classification_confidence_weight: Dict[str, float] = field(default_factory=lambda: {
        "hhsearch": 3.0,
        "domain_blast": 2.5,
        "chain_blast": 2.0,
        "blast": 1.5
    })

    # Performance settings
    use_cache: bool = True
    parallel_processing: bool = False
    batch_size: int = 100

    # Output settings
    save_intermediate: bool = False
    include_evidence_in_output: bool = True
    include_metadata: bool = True
    output_format: str = "xml"

    # Reference coverage settings
    min_reference_coverage: float = 0.7  # Minimum 70% coverage required
    strict_reference_coverage: float = 0.9  # Above this, trust boundaries exactly
    partial_coverage_threshold: float = 0.3  # Below this, reject evidence

    # Extension settings
    extend_to_reference_size: bool = True  # Try to extend to match reference
    reference_size_tolerance: float = 0.15  # Allow 15% size difference
    max_extension_length: int = 50  # Don't extend more than 50 residues

    # Coverage calculation
    use_ungapped_coverage: bool = True  # Calculate coverage excluding gaps
    combine_partial_evidence: bool = True  # Try to combine partial alignments

    def validate(self) -> None:
        """Validate options for consistency"""
        if self.min_domain_size < 1:
            raise ValueError("min_domain_size must be at least 1")

        if self.max_domain_size and self.max_domain_size < self.min_domain_size:
            raise ValueError("max_domain_size must be >= min_domain_size")

        if not 0.0 <= self.overlap_threshold <= 1.0:
            raise ValueError("overlap_threshold must be between 0.0 and 1.0")

        if not 0.0 <= self.min_evidence_confidence <= 1.0:
            raise ValueError("min_evidence_confidence must be between 0.0 and 1.0")

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization"""
        return {
            'mode': self.mode.value,
            'blast_only': self.blast_only,
            'representatives_only': self.representatives_only,
            'force_overwrite': self.force_overwrite,
            'validation_level': self.validation_level.value,
            'validate_coordinates': self.validate_coordinates,
            'validate_classifications': self.validate_classifications,
            'require_evidence': self.require_evidence,
            'use_chain_blast': self.use_chain_blast,
            'use_domain_blast': self.use_domain_blast,
            'use_hhsearch': self.use_hhsearch,
            'min_evidence_confidence': self.min_evidence_confidence,
            'min_domain_size': self.min_domain_size,
            'max_domain_size': self.max_domain_size,
            'merge_gap_tolerance': self.merge_gap_tolerance,
            'overlap_threshold': self.overlap_threshold,
            'prefer_hhsearch_classification': self.prefer_hhsearch_classification,
            'require_full_classification': self.require_full_classification,
            'classification_confidence_weight': self.classification_confidence_weight,
            'use_cache': self.use_cache,
            'parallel_processing': self.parallel_processing,
            'batch_size': self.batch_size,
            'save_intermediate': self.save_intermediate,
            'include_evidence_in_output': self.include_evidence_in_output,
            'include_metadata': self.include_metadata,
            'output_format': self.output_format
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'PartitionOptions':
        """Create from dictionary"""
        # Handle enum conversions
        if 'mode' in data and isinstance(data['mode'], str):
            data['mode'] = ProcessingMode(data['mode'])
        if 'validation_level' in data and isinstance(data['validation_level'], str):
            data['validation_level'] = ValidationLevel(data['validation_level'])

        return cls(**data)


@dataclass
class ValidationResult:
    """Result of evidence or domain validation"""
    is_valid: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    context: str = ""
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


@dataclass
class EvidenceGroup:
    """Group of related evidence items"""
    evidence_items: List[Evidence] = field(default_factory=list)
    consensus_start: Optional[int] = None
    consensus_end: Optional[int] = None
    consensus_confidence: float = 0.0
    position_key: int = 0

    def add_evidence(self, evidence: Evidence, start: int, end: int) -> None:
        """Add evidence to the group"""
        self.evidence_items.append(evidence)

        # Update consensus boundaries
        if self.consensus_start is None:
            self.consensus_start = start
            self.consensus_end = end
        else:
            # For now, use simple average - could be improved
            self.consensus_start = (self.consensus_start + start) // 2
            self.consensus_end = (self.consensus_end + end) // 2

    def calculate_consensus_confidence(self) -> float:
        """Calculate consensus confidence from all evidence"""
        if not self.evidence_items:
            return 0.0

        confidences = [e.confidence for e in self.evidence_items if e.confidence is not None]
        if not confidences:
            return 0.0

        # Weighted average based on evidence count
        return sum(confidences) / len(confidences)

    def get_best_evidence(self) -> Optional[Evidence]:
        """Get the highest confidence evidence"""
        if not self.evidence_items:
            return None

        return max(self.evidence_items,
                  key=lambda e: e.confidence if e.confidence is not None else 0.0)


@dataclass
class DomainCandidate:
    """Candidate domain during processing"""
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


@dataclass
class PartitionContext:
    """Context information for current partitioning operation"""
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


@dataclass
class BatchPartitionResults:
    """Results from batch partition processing"""
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


@dataclass
class ClassificationCache:
    """Cache for domain classifications"""
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

    def clear(self) -> None:
        """Clear all caches"""
        self.domain_id_cache.clear()
        self.chain_domain_cache.clear()
        self.cache_hits = 0
        self.cache_misses = 0

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

@dataclass
class EvidenceWithCoverage(Evidence):
    """Extended evidence with coverage information"""
    reference_coverage: float = 0.0
    reference_info: Optional[ReferenceInfo] = None
    hit_length: int = 0
    alignment_gaps: int = 0
    coverage_warning: Optional[str] = None

