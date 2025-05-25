#!/usr/bin/env python3
"""
Value objects and models for domain summary generation.

This module provides immutable data structures and models that ensure
type safety and validation throughout the summary generation process.
"""

import re
import hashlib
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum, auto
from pathlib import Path
from typing import Dict, Any, List, Optional, Set, Tuple, Union

from ecod.models.pipeline import Evidence, DomainPartitionResult
from ecod.exceptions import ValidationError


class ProcessingStatus(Enum):
    """Status of summary processing"""
    PENDING = auto()
    PROCESSING = auto()
    COMPLETED = auto()
    SKIPPED = auto()
    ERROR = auto()
    PEPTIDE = auto()


@dataclass(frozen=True)
class ProteinIdentifier:
    """
    Immutable protein identifier with validation.

    Ensures PDB ID and chain ID follow standard formats.
    """
    pdb_id: str
    chain_id: str
    reference: Optional[str] = None

    def __post_init__(self):
        """Validate identifiers on creation"""
        # Validate PDB ID format (4 alphanumeric characters)
        if not re.match(r'^[A-Za-z0-9]{4}$', self.pdb_id):
            raise ValidationError(f"Invalid PDB ID format: {self.pdb_id}")

        # Validate chain ID (typically single character, but can be longer)
        if not self.chain_id or len(self.chain_id) > 4:
            raise ValidationError(f"Invalid chain ID: {self.chain_id}")

        # Normalize to uppercase
        object.__setattr__(self, 'pdb_id', self.pdb_id.lower())
        object.__setattr__(self, 'chain_id', self.chain_id.upper())

    @property
    def source_id(self) -> str:
        """Get combined source identifier"""
        return f"{self.pdb_id}_{self.chain_id}"

    @classmethod
    def from_string(cls, protein_str: str, reference: Optional[str] = None) -> 'ProteinIdentifier':
        """
        Create from string like '1abc_A' or '1abc:A'.

        Args:
            protein_str: Protein identifier string
            reference: Optional reference version

        Returns:
            ProteinIdentifier instance
        """
        # Try different separators
        for sep in ['_', ':', '.', '-']:
            if sep in protein_str:
                parts = protein_str.split(sep, 1)
                if len(parts) == 2:
                    return cls(pdb_id=parts[0], chain_id=parts[1], reference=reference)

        # Try to parse without separator (e.g., "1abcA")
        if len(protein_str) >= 5:
            return cls(pdb_id=protein_str[:4], chain_id=protein_str[4:], reference=reference)

        raise ValidationError(f"Cannot parse protein identifier: {protein_str}")

    def __str__(self) -> str:
        """String representation"""
        return self.source_id

    def __repr__(self) -> str:
        """Detailed representation"""
        ref_str = f", reference='{self.reference}'" if self.reference else ""
        return f"ProteinIdentifier(pdb_id='{self.pdb_id}', chain_id='{self.chain_id}'{ref_str})"


@dataclass
class SequenceInfo:
    """
    Sequence information for a protein.

    Provides sequence data and utility methods for analysis.
    """
    sequence: str
    length: int = field(init=False)
    md5: str = field(init=False)

    def __post_init__(self):
        """Calculate derived fields"""
        self.length = len(self.sequence)
        self.md5 = hashlib.md5(self.sequence.encode()).hexdigest()

    def is_peptide(self, threshold: int = 30) -> bool:
        """Check if sequence is peptide length"""
        return self.length < threshold

    def is_valid_protein_sequence(self) -> bool:
        """Check if sequence contains only valid amino acids"""
        valid_aa = set('ACDEFGHIKLMNPQRSTVWYX')  # Include X for unknown
        return all(aa in valid_aa for aa in self.sequence.upper())

    def get_composition(self) -> Dict[str, float]:
        """Get amino acid composition as percentages"""
        if self.length == 0:
            return {}

        composition = {}
        for aa in 'ACDEFGHIKLMNPQRSTVWY':
            count = self.sequence.upper().count(aa)
            composition[aa] = (count / self.length) * 100

        return composition

    @classmethod
    def from_fasta(cls, fasta_path: Path) -> 'SequenceInfo':
        """
        Create from FASTA file.

        Args:
            fasta_path: Path to FASTA file

        Returns:
            SequenceInfo instance
        """
        sequence = cls._read_fasta(fasta_path)
        return cls(sequence=sequence)

    @staticmethod
    def _read_fasta(fasta_path: Path) -> str:
        """Read sequence from FASTA file"""
        if not fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

        sequence_lines = []

        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('>'):
                    sequence_lines.append(line)

        return ''.join(sequence_lines)

    def to_fasta(self, header: str = "sequence") -> str:
        """Convert to FASTA format"""
        # Format sequence in 80-character lines
        formatted_seq = '\n'.join(
            self.sequence[i:i+80]
            for i in range(0, self.length, 80)
        )
        return f">{header}\n{formatted_seq}\n"


@dataclass
class SummaryOptions:
    """
    Options for summary generation.

    Controls various aspects of evidence collection and processing.
    """
    # Processing options
    blast_only: bool = False
    force_overwrite: bool = False
    skip_filtering: bool = False

    # Evidence options
    include_evidence_details: bool = True
    stitch_discontinuous: bool = True
    max_gap_for_stitching: int = 30

    # Quality thresholds
    min_confidence: float = 0.3
    min_coverage: float = 0.0
    max_evalue: float = 10.0

    # Additional sources
    additional_sources: List[str] = field(default_factory=list)

    # Output options
    save_intermediate: bool = False
    compress_output: bool = False

    def validate(self) -> None:
        """Validate options"""
        if not 0.0 <= self.min_confidence <= 1.0:
            raise ValidationError("min_confidence must be between 0 and 1")

        if not 0.0 <= self.min_coverage <= 100.0:
            raise ValidationError("min_coverage must be between 0 and 100")

        if self.max_evalue <= 0:
            raise ValidationError("max_evalue must be positive")

        if self.max_gap_for_stitching < 0:
            raise ValidationError("max_gap_for_stitching must be non-negative")

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization"""
        return {
            'blast_only': self.blast_only,
            'force_overwrite': self.force_overwrite,
            'skip_filtering': self.skip_filtering,
            'include_evidence_details': self.include_evidence_details,
            'stitch_discontinuous': self.stitch_discontinuous,
            'max_gap_for_stitching': self.max_gap_for_stitching,
            'min_confidence': self.min_confidence,
            'min_coverage': self.min_coverage,
            'max_evalue': self.max_evalue,
            'additional_sources': self.additional_sources,
            'save_intermediate': self.save_intermediate,
            'compress_output': self.compress_output
        }

    @classmethod
    def from_dict(cls, options_dict: Dict[str, Any]) -> 'SummaryOptions':
        """Create from dictionary"""
        return cls(**options_dict)


class EvidenceCollection:
    """
    Collection of evidence from multiple sources.

    Manages evidence organization, access, and statistics.
    """

    def __init__(self):
        """Initialize empty collection"""
        self._evidence_by_source: Dict[str, List[Evidence]] = {}
        self._errors_by_source: Dict[str, str] = {}
        self._metadata: Dict[str, Any] = {}
        self._stats_cache: Optional[Dict[str, Any]] = None

    def add_source(self, source: str, evidence: List[Evidence]) -> None:
        """
        Add evidence from a source.

        Args:
            source: Source name (e.g., 'domain_blast')
            evidence: List of Evidence objects
        """
        self._evidence_by_source[source] = evidence
        self._invalidate_cache()

    def add_evidence(self, source: str, evidence: Evidence) -> None:
        """Add single evidence item to a source"""
        if source not in self._evidence_by_source:
            self._evidence_by_source[source] = []

        self._evidence_by_source[source].append(evidence)
        self._invalidate_cache()

    def add_error(self, source: str, error: str) -> None:
        """Record an error for a source"""
        self._errors_by_source[source] = error

    def get_source_evidence(self, source: str) -> List[Evidence]:
        """Get evidence from a specific source"""
        return self._evidence_by_source.get(source, [])

    def get_all_evidence(self) -> Dict[str, List[Evidence]]:
        """Get all evidence by source"""
        return self._evidence_by_source.copy()

    def get_errors(self) -> Dict[str, str]:
        """Get all errors by source"""
        return self._errors_by_source.copy()

    def has_source(self, source: str) -> bool:
        """Check if collection has evidence from a source"""
        return source in self._evidence_by_source and len(self._evidence_by_source[source]) > 0

    def get_sources(self) -> List[str]:
        """Get list of sources with evidence"""
        return [s for s, e in self._evidence_by_source.items() if e]

    def get_total_evidence_count(self) -> int:
        """Get total number of evidence items"""
        return sum(len(evidence) for evidence in self._evidence_by_source.values())

    def get_statistics(self) -> Dict[str, Any]:
        """Get detailed statistics about the collection"""
        if self._stats_cache is not None:
            return self._stats_cache

        stats = {
            'total_evidence': self.get_total_evidence_count(),
            'sources_with_evidence': len(self.get_sources()),
            'sources_with_errors': len(self._errors_by_source),
            'evidence_by_source': {
                source: len(evidence)
                for source, evidence in self._evidence_by_source.items()
            },
            'evidence_by_type': self._count_by_type(),
            'confidence_distribution': self._get_confidence_distribution(),
            'has_high_confidence': self._has_high_confidence_evidence()
        }

        self._stats_cache = stats
        return stats

    def _count_by_type(self) -> Dict[str, int]:
        """Count evidence by type"""
        type_counts = {}

        for evidence_list in self._evidence_by_source.values():
            for evidence in evidence_list:
                evidence_type = evidence.type
                type_counts[evidence_type] = type_counts.get(evidence_type, 0) + 1

        return type_counts

    def _get_confidence_distribution(self) -> Dict[str, int]:
        """Get distribution of evidence by confidence level"""
        distribution = {
            'high': 0,      # >= 0.9
            'medium': 0,    # >= 0.5
            'low': 0        # < 0.5
        }

        for evidence_list in self._evidence_by_source.values():
            for evidence in evidence_list:
                if evidence.confidence >= 0.9:
                    distribution['high'] += 1
                elif evidence.confidence >= 0.5:
                    distribution['medium'] += 1
                else:
                    distribution['low'] += 1

        return distribution

    def _has_high_confidence_evidence(self) -> bool:
        """Check if collection has any high-confidence evidence"""
        for evidence_list in self._evidence_by_source.values():
            for evidence in evidence_list:
                if evidence.confidence >= 0.9:
                    return True
        return False

    def _invalidate_cache(self) -> None:
        """Invalidate statistics cache"""
        self._stats_cache = None

    def merge(self, other: 'EvidenceCollection') -> None:
        """Merge another collection into this one"""
        for source, evidence in other.get_all_evidence().items():
            if source in self._evidence_by_source:
                self._evidence_by_source[source].extend(evidence)
            else:
                self._evidence_by_source[source] = evidence

        self._errors_by_source.update(other.get_errors())
        self._invalidate_cache()

    def filter_by_confidence(self, min_confidence: float) -> 'EvidenceCollection':
        """Create new collection with filtered evidence"""
        filtered = EvidenceCollection()

        for source, evidence_list in self._evidence_by_source.items():
            filtered_evidence = [
                e for e in evidence_list
                if e.confidence >= min_confidence
            ]
            if filtered_evidence:
                filtered.add_source(source, filtered_evidence)

        # Copy errors
        filtered._errors_by_source = self._errors_by_source.copy()

        return filtered

    def get_summary_string(self) -> str:
        """Get summary string for logging"""
        parts = []

        for source in ['self_comparison', 'chain_blast', 'domain_blast', 'hhsearch']:
            if self.has_source(source):
                count = len(self.get_source_evidence(source))
                parts.append(f"{source}={count}")

        return ", ".join(parts) if parts else "no evidence"


@dataclass
class EvidenceSummary:
    """
    Complete evidence summary for a protein.

    Encapsulates all evidence, metadata, and processing information.
    """
    pdb_id: str
    chain_id: str
    reference: str

    # Status
    status: ProcessingStatus = ProcessingStatus.PENDING

    # Sequence information
    sequence_info: Optional[SequenceInfo] = None

    # Evidence collection
    evidence_collection: Optional[EvidenceCollection] = None

    # Processing metadata
    processing_time: float = 0.0
    timestamp: datetime = field(default_factory=datetime.now)
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    # Additional metadata
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def protein_id(self) -> str:
        """Get protein identifier"""
        return f"{self.pdb_id}_{self.chain_id}"

    @property
    def total_evidence_count(self) -> int:
        """Get total evidence count"""
        if self.evidence_collection:
            return self.evidence_collection.get_total_evidence_count()
        return 0

    def has_evidence(self) -> bool:
        """Check if summary has any evidence"""
        return self.total_evidence_count > 0

    def mark_as_completed(self) -> None:
        """Mark processing as completed"""
        self.status = ProcessingStatus.COMPLETED

    def mark_as_skipped(self, reason: str) -> None:
        """Mark as skipped with reason"""
        self.status = ProcessingStatus.SKIPPED
        self.metadata['skip_reason'] = reason

    def mark_as_error(self, error: str) -> None:
        """Mark as error with message"""
        self.status = ProcessingStatus.ERROR
        self.errors.append(error)

    def mark_as_peptide(self) -> None:
        """Mark as peptide"""
        self.status = ProcessingStatus.PEPTIDE
        self.metadata['is_peptide'] = True

    def set_sequence_info(self, sequence_info: Optional[SequenceInfo]) -> None:
        """Set sequence information"""
        self.sequence_info = sequence_info

        if sequence_info:
            self.metadata['sequence_length'] = sequence_info.length
            self.metadata['sequence_md5'] = sequence_info.md5

    def add_evidence_collection(self, collection: EvidenceCollection) -> None:
        """Add evidence collection"""
        self.evidence_collection = collection

        # Update status if successful
        if self.status == ProcessingStatus.PENDING:
            self.status = ProcessingStatus.COMPLETED

    def add_warning(self, warning: str) -> None:
        """Add a warning message"""
        self.warnings.append(warning)

    def to_partition_result(self) -> DomainPartitionResult:
        """
        Convert to DomainPartitionResult for downstream processing.

        Returns:
            DomainPartitionResult with evidence
        """
        result = DomainPartitionResult(
            pdb_id=self.pdb_id,
            chain_id=self.chain_id,
            reference=self.reference,
            sequence_length=self.sequence_info.length if self.sequence_info else 0
        )

        # Set status flags
        if self.status == ProcessingStatus.PEPTIDE:
            result.is_peptide = True
        elif self.status == ProcessingStatus.ERROR:
            result.success = False
            result.error = "; ".join(self.errors)

        # Add evidence to result metadata
        if self.evidence_collection:
            all_evidence = []
            for source, evidence_list in self.evidence_collection.get_all_evidence().items():
                all_evidence.extend(evidence_list)

            result.evidence_summary = {
                "all_evidence": all_evidence,
                "statistics": self.evidence_collection.get_statistics()
            }

        return result

    def get_summary_dict(self) -> Dict[str, Any]:
        """Get summary as dictionary for reporting"""
        summary = {
            'protein_id': self.protein_id,
            'pdb_id': self.pdb_id,
            'chain_id': self.chain_id,
            'reference': self.reference,
            'status': self.status.name,
            'total_evidence': self.total_evidence_count,
            'processing_time': self.processing_time,
            'timestamp': self.timestamp.isoformat()
        }

        # Add sequence info
        if self.sequence_info:
            summary['sequence_length'] = self.sequence_info.length
            summary['is_peptide'] = self.sequence_info.is_peptide()

        # Add evidence statistics
        if self.evidence_collection:
            summary['evidence_stats'] = self.evidence_collection.get_statistics()

        # Add errors/warnings
        if self.errors:
            summary['errors'] = self.errors
        if self.warnings:
            summary['warnings'] = self.warnings

        # Add metadata
        summary['metadata'] = self.metadata

        return summary


@dataclass
class ProcessingStats:
    """Statistics tracking for batch processing"""
    proteins_processed: int = 0
    proteins_with_evidence: int = 0
    peptides_processed: int = 0
    errors: int = 0
    skipped: int = 0

    # Evidence counts by source
    evidence_by_source: Dict[str, int] = field(default_factory=dict)

    # Timing
    total_processing_time: float = 0.0
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None

    def start_tracking(self) -> None:
        """Start time tracking"""
        self.start_time = datetime.now()

    def end_tracking(self) -> None:
        """End time tracking"""
        self.end_time = datetime.now()
        if self.start_time:
            self.total_processing_time = (self.end_time - self.start_time).total_seconds()

    @property
    def success_rate(self) -> float:
        """Calculate success rate"""
        if self.proteins_processed == 0:
            return 0.0
        return (self.proteins_with_evidence / self.proteins_processed) * 100

    @property
    def average_processing_time(self) -> float:
        """Calculate average processing time per protein"""
        if self.proteins_processed == 0:
            return 0.0
        return self.total_processing_time / self.proteins_processed

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for reporting"""
        return {
            'proteins_processed': self.proteins_processed,
            'proteins_with_evidence': self.proteins_with_evidence,
            'peptides_processed': self.peptides_processed,
            'errors': self.errors,
            'skipped': self.skipped,
            'success_rate': self.success_rate,
            'evidence_by_source': dict(self.evidence_by_source),
            'total_processing_time': self.total_processing_time,
            'average_processing_time': self.average_processing_time,
            'start_time': self.start_time.isoformat() if self.start_time else None,
            'end_time': self.end_time.isoformat() if self.end_time else None
        }


@dataclass
class BatchResults:
    """Results from batch processing"""
    successes: List[EvidenceSummary] = field(default_factory=list)
    failures: List[Tuple[str, str, str]] = field(default_factory=list)  # (pdb_id, chain_id, error)

    @property
    def total(self) -> int:
        """Total proteins processed"""
        return len(self.successes) + len(self.failures)

    @property
    def success_count(self) -> int:
        """Number of successful processings"""
        return len(self.successes)

    @property
    def failure_count(self) -> int:
        """Number of failed processings"""
        return len(self.failures)

    def add_success(self, summary: EvidenceSummary) -> None:
        """Add successful result"""
        self.successes.append(summary)

    def add_failure(self, pdb_id: str, chain_id: str, error: str) -> None:
        """Add failed result"""
        self.failures.append((pdb_id, chain_id, error))

    def get_statistics(self) -> ProcessingStats:
        """Generate statistics from results"""
        stats = ProcessingStats()
        stats.proteins_processed = self.total

        # Count successes with evidence
        for summary in self.successes:
            if summary.has_evidence():
                stats.proteins_with_evidence += 1

            if summary.status == ProcessingStatus.PEPTIDE:
                stats.peptides_processed += 1
            elif summary.status == ProcessingStatus.ERROR:
                stats.errors += 1
            elif summary.status == ProcessingStatus.SKIPPED:
                stats.skipped += 1

            # Count evidence by source
            if summary.evidence_collection:
                for source, evidence in summary.evidence_collection.get_all_evidence().items():
                    stats.evidence_by_source[source] = \
                        stats.evidence_by_source.get(source, 0) + len(evidence)

        stats.errors += self.failure_count

        return stats
