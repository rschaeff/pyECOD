# ecod/pipelines/domain_analysis/partition/models.py
"""
Models and configuration for domain partition processing.
Defines options, validation levels, and result structures.
"""

from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import List, Any, Dict, Optional
from pathlib import Path

from ecod.models.pipeline.evidence import Evidence
from ecod.models.pipeline.domain import DomainModel
from ecod.models.pipeline.partition import DomainPartitionResult


class ValidationLevel(Enum):
    """Validation strictness levels"""
    LENIENT = "lenient"      # Allow minor issues, focus on processing
    NORMAL = "normal"        # Standard validation
    STRICT = "strict"        # Strict validation, fail on any issues


@dataclass
class PartitionOptions:
    """Configuration options for domain partitioning"""

    # Validation settings
    validation_level: ValidationLevel = ValidationLevel.NORMAL

    # Domain size constraints
    min_domain_size: int = 30
    max_domain_size: int = 2000

    # Overlap handling
    overlap_threshold: float = 0.3  # Maximum allowed overlap (0-1)
    resolve_overlaps: bool = True

    # Evidence filtering
    min_evidence_confidence: float = 0.1
    require_evidence: bool = True

    # Processing options
    use_cache: bool = True
    parallel_processing: bool = False
    max_workers: int = 4

    # Output options
    include_evidence: bool = True
    include_metadata: bool = True
    save_intermediate: bool = False

    def validate(self) -> None:
        """Validate option values"""
        errors = []

        if self.min_domain_size < 1:
            errors.append("min_domain_size must be positive")

        if self.max_domain_size <= self.min_domain_size:
            errors.append("max_domain_size must be greater than min_domain_size")

        if not 0.0 <= self.overlap_threshold <= 1.0:
            errors.append("overlap_threshold must be between 0.0 and 1.0")

        if not 0.0 <= self.min_evidence_confidence <= 1.0:
            errors.append("min_evidence_confidence must be between 0.0 and 1.0")

        if self.max_workers < 1:
            errors.append("max_workers must be positive")

        if errors:
            raise ValueError(f"Invalid partition options: {'; '.join(errors)}")

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'validation_level': self.validation_level.value,
            'min_domain_size': self.min_domain_size,
            'max_domain_size': self.max_domain_size,
            'overlap_threshold': self.overlap_threshold,
            'resolve_overlaps': self.resolve_overlaps,
            'min_evidence_confidence': self.min_evidence_confidence,
            'require_evidence': self.require_evidence,
            'use_cache': self.use_cache,
            'parallel_processing': self.parallel_processing,
            'max_workers': self.max_workers,
            'include_evidence': self.include_evidence,
            'include_metadata': self.include_metadata,
            'save_intermediate': self.save_intermediate
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'PartitionOptions':
        """Create from dictionary"""
        # Handle validation level enum
        validation_level = ValidationLevel.NORMAL
        if 'validation_level' in data:
            level_str = data['validation_level']
            if isinstance(level_str, str):
                try:
                    validation_level = ValidationLevel(level_str)
                except ValueError:
                    validation_level = ValidationLevel.NORMAL

        return cls(
            validation_level=validation_level,
            min_domain_size=data.get('min_domain_size', 30),
            max_domain_size=data.get('max_domain_size', 2000),
            overlap_threshold=data.get('overlap_threshold', 0.3),
            resolve_overlaps=data.get('resolve_overlaps', True),
            min_evidence_confidence=data.get('min_evidence_confidence', 0.1),
            require_evidence=data.get('require_evidence', True),
            use_cache=data.get('use_cache', True),
            parallel_processing=data.get('parallel_processing', False),
            max_workers=data.get('max_workers', 4),
            include_evidence=data.get('include_evidence', True),
            include_metadata=data.get('include_metadata', True),
            save_intermediate=data.get('save_intermediate', False)
        )


@dataclass
class ValidationResult:
    """Result of validation operation"""

    is_valid: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    context: str = ""
    validation_level: ValidationLevel = ValidationLevel.NORMAL

    def add_error(self, message: str) -> None:
        """Add an error message"""
        self.errors.append(message)
        self.is_valid = False

    def add_warning(self, message: str) -> None:
        """Add a warning message"""
        self.warnings.append(message)

    def has_errors(self) -> bool:
        """Check if there are any errors"""
        return len(self.errors) > 0

    def has_warnings(self) -> bool:
        """Check if there are any warnings"""
        return len(self.warnings) > 0

    def get_summary(self) -> str:
        """Get summary of validation result"""
        if self.is_valid:
            if self.has_warnings():
                return f"Valid with {len(self.warnings)} warnings"
            else:
                return "Valid"
        else:
            return f"Invalid: {len(self.errors)} errors, {len(self.warnings)} warnings"

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'is_valid': self.is_valid,
            'errors': self.errors,
            'warnings': self.warnings,
            'context': self.context,
            'validation_level': self.validation_level.value,
            'summary': self.get_summary()
        }


@dataclass
class ProcessingStats:
    """Statistics for processing operations"""

    total_items: int = 0
    processed_items: int = 0
    failed_items: int = 0
    skipped_items: int = 0

    start_time: Optional[float] = None
    end_time: Optional[float] = None

    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    @property
    def success_rate(self) -> float:
        """Calculate success rate"""
        if self.total_items == 0:
            return 0.0
        return self.processed_items / self.total_items

    @property
    def processing_time(self) -> Optional[float]:
        """Calculate processing time"""
        if self.start_time is not None and self.end_time is not None:
            return self.end_time - self.start_time
        return None

    def add_error(self, message: str) -> None:
        """Add an error"""
        self.errors.append(message)
        self.failed_items += 1

    def add_warning(self, message: str) -> None:
        """Add a warning"""
        self.warnings.append(message)

    def record_success(self) -> None:
        """Record a successful item"""
        self.processed_items += 1

    def record_failure(self, error_message: str = "") -> None:
        """Record a failed item"""
        self.failed_items += 1
        if error_message:
            self.errors.append(error_message)

    def record_skip(self, reason: str = "") -> None:
        """Record a skipped item"""
        self.skipped_items += 1
        if reason:
            self.warnings.append(f"Skipped: {reason}")

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'total_items': self.total_items,
            'processed_items': self.processed_items,
            'failed_items': self.failed_items,
            'skipped_items': self.skipped_items,
            'success_rate': self.success_rate,
            'processing_time': self.processing_time,
            'start_time': self.start_time,
            'end_time': self.end_time,
            'errors': self.errors,
            'warnings': self.warnings
        }


@dataclass
class EvidenceGroup:
    """Group of evidence items with common characteristics"""

    group_id: str
    evidence_items: List[Any] = field(default_factory=list)
    position_range: Optional[tuple] = None  # (start, end)
    confidence_score: float = 0.0
    group_type: str = "unknown"  # "blast", "hhsearch", "mixed"

    def add_evidence(self, evidence) -> None:
        """Add evidence to group"""
        self.evidence_items.append(evidence)
        self._update_group_stats()

    def _update_group_stats(self) -> None:
        """Update group statistics"""
        if not self.evidence_items:
            return

        # Update confidence score (average)
        confidences = []
        for evidence in self.evidence_items:
            if hasattr(evidence, 'confidence') and evidence.confidence is not None:
                confidences.append(evidence.confidence)

        if confidences:
            self.confidence_score = sum(confidences) / len(confidences)

        # Update group type
        types = set()
        for evidence in self.evidence_items:
            if hasattr(evidence, 'type'):
                types.add(evidence.type)

        if len(types) == 1:
            self.group_type = types.pop()
        else:
            self.group_type = "mixed"

    def get_position_range(self) -> Optional[tuple]:
        """Get the position range covered by this group"""
        positions = []

        for evidence in self.evidence_items:
            if hasattr(evidence, 'query_range') and evidence.query_range:
                try:
                    if "-" in evidence.query_range:
                        start, end = evidence.query_range.split("-")
                        positions.extend([int(start), int(end)])
                except (ValueError, IndexError):
                    continue

        if positions:
            return (min(positions), max(positions))

        return None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'group_id': self.group_id,
            'evidence_count': len(self.evidence_items),
            'position_range': self.get_position_range(),
            'confidence_score': self.confidence_score,
            'group_type': self.group_type
        }


@dataclass
class EvidenceQualityMetrics:
    """Quality metrics for evidence analysis"""
    total_evidence_count: int = 0
    valid_evidence_count: int = 0
    blast_hit_count: int = 0
    hhsearch_hit_count: int = 0
    self_comparison_count: int = 0

    # Coverage metrics
    sequence_coverage: float = 0.0
    gap_count: int = 0
    largest_gap_size: int = 0
    overlap_count: int = 0
    largest_overlap_size: int = 0

    # Quality scores
    average_confidence: float = 0.0
    confidence_std: float = 0.0
    high_confidence_count: int = 0  # > 0.8
    low_confidence_count: int = 0   # < 0.3

    # Classification metrics
    classified_evidence_count: int = 0
    classification_consistency: float = 0.0  # How consistent are classifications

    # Performance metrics
    processing_time_seconds: float = 0.0
    memory_usage_mb: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'total_evidence': self.total_evidence_count,
            'valid_evidence': self.valid_evidence_count,
            'blast_hits': self.blast_hit_count,
            'hhsearch_hits': self.hhsearch_hit_count,
            'self_comparison': self.self_comparison_count,
            'sequence_coverage': self.sequence_coverage,
            'gap_count': self.gap_count,
            'largest_gap': self.largest_gap_size,
            'overlap_count': self.overlap_count,
            'largest_overlap': self.largest_overlap_size,
            'average_confidence': self.average_confidence,
            'confidence_std': self.confidence_std,
            'high_confidence_count': self.high_confidence_count,
            'low_confidence_count': self.low_confidence_count,
            'classified_evidence': self.classified_evidence_count,
            'classification_consistency': self.classification_consistency,
            'processing_time': self.processing_time_seconds,
            'memory_usage_mb': self.memory_usage_mb
        }


@dataclass
class AnalysisResult:
    """Comprehensive analysis result from EvidenceAnalyzer"""
    success: bool
    file_path: str
    protein_id: Optional[str] = None
    sequence_length: int = 0

    # Evidence data
    evidence_count: int = 0
    filtered_evidence_count: int = 0
    evidence_groups: List[EvidenceGroup] = field(default_factory=list)
    resolved_evidence: List[Any] = field(default_factory=list)  # Evidence objects

    # Analysis results
    domain_suggestions: List[Dict[str, Any]] = field(default_factory=list)
    quality_metrics: Optional[EvidenceQualityMetrics] = None
    classification_analysis: Dict[str, Any] = field(default_factory=dict)
    validation_summary: Dict[str, Any] = field(default_factory=dict)

    # Performance data
    processing_time_seconds: float = 0.0
    cache_stats: Optional[Dict[str, Any]] = None

    # Error information
    error: Optional[str] = None
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        result = {
            'success': self.success,
            'file_path': self.file_path,
            'protein_id': self.protein_id,
            'sequence_length': self.sequence_length,
            'evidence_count': self.evidence_count,
            'filtered_evidence_count': self.filtered_evidence_count,
            'domain_suggestions': self.domain_suggestions,
            'classification_analysis': self.classification_analysis,
            'validation_summary': self.validation_summary,
            'processing_time_seconds': self.processing_time_seconds,
            'cache_stats': self.cache_stats,
            'error': self.error,
            'warnings': self.warnings
        }

        # Add quality metrics if available
        if self.quality_metrics:
            result['quality_metrics'] = self.quality_metrics.to_dict()

        # Add evidence groups
        result['evidence_groups'] = [group.to_dict() for group in self.evidence_groups]

        return result

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'AnalysisResult':
        """Create from dictionary"""
        result = cls(
            success=data.get('success', False),
            file_path=data.get('file_path', ''),
            protein_id=data.get('protein_id'),
            sequence_length=data.get('sequence_length', 0),
            evidence_count=data.get('evidence_count', 0),
            filtered_evidence_count=data.get('filtered_evidence_count', 0),
            domain_suggestions=data.get('domain_suggestions', []),
            classification_analysis=data.get('classification_analysis', {}),
            validation_summary=data.get('validation_summary', {}),
            processing_time_seconds=data.get('processing_time_seconds', 0.0),
            cache_stats=data.get('cache_stats'),
            error=data.get('error'),
            warnings=data.get('warnings', [])
        )

        # Create quality metrics if present
        if 'quality_metrics' in data:
            metrics_data = data['quality_metrics']
            result.quality_metrics = EvidenceQualityMetrics(
                total_evidence_count=metrics_data.get('total_evidence', 0),
                valid_evidence_count=metrics_data.get('valid_evidence', 0),
                blast_hit_count=metrics_data.get('blast_hits', 0),
                hhsearch_hit_count=metrics_data.get('hhsearch_hits', 0),
                self_comparison_count=metrics_data.get('self_comparison', 0),
                sequence_coverage=metrics_data.get('sequence_coverage', 0.0),
                gap_count=metrics_data.get('gap_count', 0),
                largest_gap_size=metrics_data.get('largest_gap', 0),
                overlap_count=metrics_data.get('overlap_count', 0),
                largest_overlap_size=metrics_data.get('largest_overlap', 0),
                average_confidence=metrics_data.get('average_confidence', 0.0),
                confidence_std=metrics_data.get('confidence_std', 0.0),
                high_confidence_count=metrics_data.get('high_confidence_count', 0),
                low_confidence_count=metrics_data.get('low_confidence_count', 0),
                classified_evidence_count=metrics_data.get('classified_evidence', 0),
                classification_consistency=metrics_data.get('classification_consistency', 0.0),
                processing_time_seconds=metrics_data.get('processing_time', 0.0),
                memory_usage_mb=metrics_data.get('memory_usage_mb', 0.0)
            )

        return result


@dataclass
class CacheStatistics:
    """Statistics for classification cache performance"""
    size: int = 0
    max_size: int = 0
    hits: int = 0
    misses: int = 0
    evictions: int = 0
    hit_rate_percent: float = 0.0
    ttl_hours: float = 24.0

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'size': self.size,
            'max_size': self.max_size,
            'hits': self.hits,
            'misses': self.misses,
            'evictions': self.evictions,
            'hit_rate_percent': self.hit_rate_percent,
            'ttl_hours': self.ttl_hours
        }

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
class ProcessingMetrics:
    """Metrics for batch processing operations"""
    files_processed: int = 0
    total_processing_time: float = 0.0
    average_processing_time: float = 0.0
    median_processing_time: float = 0.0
    fastest_processing_time: float = 0.0
    slowest_processing_time: float = 0.0

    # Error statistics
    parse_errors: int = 0
    validation_errors: int = 0
    classification_errors: int = 0

    # Cache performance
    cache_hit_rate: float = 0.0
    cache_misses: int = 0

    # Resource usage
    peak_memory_mb: float = 0.0
    average_memory_mb: float = 0.0
    parallel_efficiency: float = 0.0  # For parallel processing

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'files_processed': self.files_processed,
            'total_processing_time': self.total_processing_time,
            'average_processing_time': self.average_processing_time,
            'median_processing_time': self.median_processing_time,
            'fastest_processing_time': self.fastest_processing_time,
            'slowest_processing_time': self.slowest_processing_time,
            'parse_errors': self.parse_errors,
            'validation_errors': self.validation_errors,
            'classification_errors': self.classification_errors,
            'cache_hit_rate': self.cache_hit_rate,
            'cache_misses': self.cache_misses,
            'peak_memory_mb': self.peak_memory_mb,
            'average_memory_mb': self.average_memory_mb,
            'parallel_efficiency': self.parallel_efficiency
        }
