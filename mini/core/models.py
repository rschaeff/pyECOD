# mini/core/models.py (enhanced version with reference coverage tracking)
"""Enhanced models for domain analysis with boundary optimization and comprehensive provenance tracking"""

from dataclasses import dataclass, field
from typing import Optional, List, Set, Tuple, Any, Dict
from enum import Enum
from datetime import datetime
from .sequence_range import SequenceRange


class SegmentType(Enum):
    """Types of unassigned segments"""
    NTERM = "nterm"          # Before first domain
    CTERM = "cterm"          # After last domain
    INTERSTITIAL = "inter"   # Between domains
    SINGLETON = "singleton"  # Single residue


class FragmentSize(Enum):
    """Fragment size classification"""
    SMALL = "small"          # < min_domain_size (candidate for merging)
    LARGE = "large"          # >= min_domain_size (candidate for missed domain)


@dataclass
class AlignmentData:
    """Alignment information for decomposition"""
    query_seq: str
    hit_seq: str
    query_start: int
    query_end: int
    hit_start: int
    hit_end: int


@dataclass
class Evidence:
    """Enhanced evidence model with boundary optimization and comprehensive provenance support"""
    type: str  # 'chain_blast', 'domain_blast', 'hhsearch'
    source_pdb: str  # '6dgv', '2ia4', etc.
    query_range: SequenceRange
    confidence: float = 0.0
    evalue: Optional[float] = None
    domain_id: Optional[str] = None  # 'e6dgvA1'

    # Classification for better partitioning
    t_group: Optional[str] = None
    h_group: Optional[str] = None

    # Reference info for coverage calculation
    reference_length: Optional[int] = None
    alignment_coverage: Optional[float] = None

    # ENHANCED: Reference coverage tracking
    reference_coverage: Optional[float] = None  # Fraction of reference domain covered

    # Optional alignment data (for chain BLAST decomposition)
    alignment: Optional[AlignmentData] = None

    # PROVENANCE TRACKING: Enhanced evidence details
    hit_range: Optional[SequenceRange] = None  # Range in the hit/reference sequence
    hsp_count: Optional[int] = None  # Number of high-scoring pairs
    discontinuous: bool = False  # Whether evidence spans discontinuous regions

    # Source identification for provenance
    source_chain_id: Optional[str] = None  # Chain ID from source PDB
    source_domain_uid: Optional[str] = None  # ECOD UID if available

    def get_positions(self) -> Set[int]:
        """Get all sequence positions covered by this evidence"""
        return set(self.query_range.to_positions_simple())

    def get_reference_coverage(self) -> Optional[float]:
        """Get reference coverage, calculating if not already set"""
        if self.reference_coverage is not None:
            return self.reference_coverage

        # Calculate if we have the necessary data
        if self.hit_range and self.reference_length and self.reference_length > 0:
            return self.hit_range.total_length / self.reference_length

        return None

    def get_quality_metrics(self) -> Dict[str, Any]:
        """Get comprehensive quality metrics for this evidence"""
        metrics = {
            'confidence': self.confidence,
            'evalue': self.evalue,
            'reference_coverage': self.get_reference_coverage(),
            'query_length': self.query_range.total_length,
            'hit_length': self.hit_range.total_length if self.hit_range else None,
            'reference_length': self.reference_length,
            'discontinuous': self.discontinuous
        }

        # Type-specific metrics
        if self.type == 'hhsearch':
            # Calculate original probability from confidence
            original_prob = self.confidence / 0.95 * 100  # Reverse type multiplier
            metrics['original_hhsearch_probability'] = original_prob

        return metrics

    def to_provenance_dict(self) -> dict:
        """Export evidence as provenance dictionary for XML/database storage"""
        return {
            'source_type': self.type,
            'source_id': f"{self.source_pdb}_{self.source_chain_id}" if self.source_chain_id else self.source_pdb,
            'domain_id': self.domain_id,
            'evalue': self.evalue,
            'confidence': self.confidence,
            'query_range': str(self.query_range),
            'hit_range': str(self.hit_range) if self.hit_range else None,
            'reference_coverage': self.get_reference_coverage(),
            'hsp_count': self.hsp_count,
            'discontinuous': self.discontinuous
        }


@dataclass
class Domain:
    """Enhanced domain model with spatial tracking and comprehensive provenance"""
    id: str
    range: SequenceRange
    family: str  # PDB or classification
    evidence_count: int
    source: str  # 'blast' or 'hhsearch'
    evidence_items: List[Evidence]

    # ECOD hierarchy fields
    x_group: Optional[str] = None
    h_group: Optional[str] = None
    t_group: Optional[str] = None

    # Spatial tracking for boundary optimization
    assigned_positions: Set[int] = field(default_factory=set)

    # PROVENANCE TRACKING: Enhanced domain provenance
    primary_evidence: Optional[Evidence] = None  # The main evidence that created this domain
    reference_ecod_domain_id: Optional[str] = None  # Final reference domain assigned
    original_range: Optional[SequenceRange] = None  # Range before boundary optimization
    optimization_actions: List[str] = field(default_factory=list)  # Actions taken during optimization

    # Context for provenance
    creation_timestamp: Optional[datetime] = None
    confidence_score: Optional[float] = None  # Overall domain confidence

    def __post_init__(self):
        """Populate assigned_positions and set defaults during domain creation"""
        if not self.assigned_positions:
            self.assigned_positions = set(self.range.to_positions_simple())

        # Set original range to current range if not specified (before optimization)
        if self.original_range is None:
            self.original_range = self.range

        # Set creation timestamp
        if self.creation_timestamp is None:
            self.creation_timestamp = datetime.now()

        # Set primary evidence if not specified but evidence_items exist
        if self.primary_evidence is None and self.evidence_items:
            # Use the evidence with highest confidence as primary
            self.primary_evidence = max(self.evidence_items, key=lambda e: e.confidence)

        # Set domain confidence from primary evidence
        if self.confidence_score is None and self.primary_evidence:
            self.confidence_score = self.primary_evidence.confidence

    @property
    def start_position(self) -> int:
        """Get the first position of this domain"""
        return min(self.assigned_positions) if self.assigned_positions else 0

    @property
    def end_position(self) -> int:
        """Get the last position of this domain"""
        return max(self.assigned_positions) if self.assigned_positions else 0

    @property
    def length(self) -> int:
        """Get the number of assigned positions"""
        return len(self.assigned_positions)

    @property
    def confidence(self) -> float:
        """Get confidence from the best evidence item (backward compatibility)"""
        if self.confidence_score is not None:
            return self.confidence_score
        if not self.evidence_items:
            return 0.0
        return max(ev.confidence for ev in self.evidence_items)

    def get_positions(self) -> Set[int]:
        """Get all sequence positions covered by this domain"""
        return self.assigned_positions.copy()

    def add_positions(self, positions: Set[int]) -> None:
        """Add positions to this domain (for fragment merging)"""
        self.assigned_positions.update(positions)
        # Update the range representation
        all_positions = sorted(self.assigned_positions)
        self.range = SequenceRange.from_positions(all_positions)

    def remove_positions(self, positions: Set[int]) -> None:
        """Remove positions from this domain (for overlap resolution)"""
        self.assigned_positions.difference_update(positions)
        if self.assigned_positions:
            all_positions = sorted(self.assigned_positions)
            self.range = SequenceRange.from_positions(all_positions)

    def overlaps_with(self, other: 'Domain', tolerance: int = 0) -> bool:
        """Check if this domain overlaps with another domain"""
        if tolerance == 0:
            return bool(self.assigned_positions.intersection(other.assigned_positions))
        else:
            # Extended overlap check with tolerance
            extended_positions = set()
            for pos in self.assigned_positions:
                extended_positions.update(range(pos - tolerance, pos + tolerance + 1))
            return bool(extended_positions.intersection(other.assigned_positions))

    def distance_to(self, other: 'Domain') -> int:
        """Get minimum distance between this domain and another domain"""
        if self.overlaps_with(other):
            return 0

        min_distance = float('inf')
        for pos1 in self.assigned_positions:
            for pos2 in other.assigned_positions:
                min_distance = min(min_distance, abs(pos1 - pos2))

        return int(min_distance) if min_distance != float('inf') else -1

    # PROVENANCE TRACKING: Methods for tracking optimization
    def record_optimization_action(self, action: str, details: Optional[str] = None) -> None:
        """Record a boundary optimization action"""
        if details:
            self.optimization_actions.append(f"{action}:{details}")
        else:
            self.optimization_actions.append(action)

    def was_optimized(self) -> bool:
        """Check if domain boundaries were modified from original evidence"""
        return self.original_range != self.range

    def get_optimization_summary(self) -> dict:
        """Get summary of all optimizations performed on this domain"""
        return {
            'was_optimized': self.was_optimized(),
            'original_range': str(self.original_range),
            'final_range': str(self.range),
            'actions': self.optimization_actions.copy(),
            'position_change': len(self.assigned_positions) - len(set(self.original_range.to_positions_simple()))
        }

    def get_quality_assessment(self) -> Dict[str, Any]:
        """Get comprehensive quality assessment for this domain"""
        assessment = {
            'domain_id': self.id,
            'primary_evidence_quality': None,
            'all_evidence_quality': [],
            'overall_assessment': 'unknown'
        }

        if self.primary_evidence:
            assessment['primary_evidence_quality'] = self.primary_evidence.get_quality_metrics()

        for evidence in self.evidence_items:
            assessment['all_evidence_quality'].append(evidence.get_quality_metrics())

        # Overall assessment based on primary evidence
        if self.primary_evidence:
            confidence = self.primary_evidence.confidence
            ref_coverage = self.primary_evidence.get_reference_coverage()

            issues = []
            if confidence < 0.5:
                issues.append('low_confidence')
            if ref_coverage is not None and ref_coverage < 0.5:
                issues.append('poor_reference_coverage')
            if self.primary_evidence.evalue is not None and self.primary_evidence.evalue > 1.0:
                issues.append('poor_evalue')

            if not issues:
                assessment['overall_assessment'] = 'good'
            elif len(issues) == 1:
                assessment['overall_assessment'] = 'questionable'
            else:
                assessment['overall_assessment'] = 'poor'

            assessment['quality_issues'] = issues

        return assessment

    def to_provenance_dict(self) -> dict:
        """Export domain as provenance dictionary for XML/database storage"""
        base_dict = {
            'id': self.id,
            'range': str(self.range),
            'family': self.family,
            'source': self.source,
            'evidence_count': self.evidence_count,
            'confidence': self.confidence,
            't_group': self.t_group,
            'h_group': self.h_group,
            'x_group': self.x_group,
            'reference_ecod_domain_id': self.reference_ecod_domain_id,
            'was_optimized': self.was_optimized(),
            'original_range': str(self.original_range),
            'optimization_actions': self.optimization_actions.copy()
        }

        if self.primary_evidence:
            base_dict['primary_evidence'] = self.primary_evidence.to_provenance_dict()

        return base_dict


@dataclass
class UnassignedSegment:
    """Represents a gap between domains for boundary optimization"""
    start: int
    end: int
    positions: Set[int] = field(default_factory=set)
    segment_type: SegmentType = SegmentType.INTERSTITIAL
    fragment_size: FragmentSize = FragmentSize.SMALL

    # Spatial relationships
    preceding_domain: Optional[Domain] = None  # Domain that ends before this segment
    following_domain: Optional[Domain] = None  # Domain that starts after this segment
    neighboring_domains: List[Domain] = field(default_factory=list)  # All nearby domains

    def __post_init__(self):
        """Initialize positions from start/end"""
        if not self.positions:
            self.positions = set(range(self.start, self.end + 1))

    @property
    def length(self) -> int:
        """Get the length of this segment"""
        return len(self.positions)

    @property
    def is_singleton(self) -> bool:
        """Check if this is a single residue"""
        return self.length == 1

    @property
    def is_small_fragment(self) -> bool:
        """Check if this is a small fragment (based on fragment_size)"""
        return self.fragment_size == FragmentSize.SMALL

    @property
    def range(self) -> SequenceRange:
        """Get SequenceRange representation of this segment"""
        return SequenceRange.from_positions(sorted(self.positions))

    @classmethod
    def from_positions(cls, positions: List[int], min_domain_size: int = 25) -> 'UnassignedSegment':
        """Create UnassignedSegment from list of positions"""
        if not positions:
            raise ValueError("Cannot create segment from empty positions")

        positions_set = set(positions)
        start = min(positions)
        end = max(positions)

        # Determine fragment size
        fragment_size = FragmentSize.SMALL if len(positions) < min_domain_size else FragmentSize.LARGE

        return cls(
            start=start,
            end=end,
            positions=positions_set,
            fragment_size=fragment_size
        )

    def classify_type(self, all_domains: List[Domain], sequence_length: int) -> SegmentType:
        """Classify this segment based on its position relative to domains"""
        if self.is_singleton:
            return SegmentType.SINGLETON

        # Find domains that start before this segment's end
        domains_before = [d for d in all_domains if d.end_position < self.start]
        # Find domains that start after this segment's start
        domains_after = [d for d in all_domains if d.start_position > self.end]

        if not domains_before and domains_after:
            return SegmentType.NTERM
        elif domains_before and not domains_after:
            return SegmentType.CTERM
        elif domains_before and domains_after:
            return SegmentType.INTERSTITIAL
        else:
            # No domains at all - treat as N-terminal
            return SegmentType.NTERM

    def find_neighbors(self, all_domains: List[Domain], tolerance: int = 5) -> List[Domain]:
        """Find domains within tolerance distance of this segment"""
        neighbors = []

        for domain in all_domains:
            # Check if domain is close to segment start or end
            domain_start = domain.start_position
            domain_end = domain.end_position

            # Distance from domain end to segment start (preceding domain)
            if abs(domain_end - (self.start - 1)) <= tolerance:
                neighbors.append(domain)
                if not self.preceding_domain or domain_end > self.preceding_domain.end_position:
                    self.preceding_domain = domain

            # Distance from segment end to domain start (following domain)
            elif abs((self.end + 1) - domain_start) <= tolerance:
                neighbors.append(domain)
                if not self.following_domain or domain_start < self.following_domain.start_position:
                    self.following_domain = domain

        return neighbors


@dataclass
class DomainLayout:
    """Manages spatial relationships between domains and gaps for boundary optimization"""
    sequence_length: int
    domains: List[Domain] = field(default_factory=list)
    unassigned_segments: List[UnassignedSegment] = field(default_factory=list)

    # Tracking for optimization
    used_positions: Set[int] = field(default_factory=set)
    unused_positions: Set[int] = field(default_factory=set)

    def __post_init__(self):
        """Initialize position tracking"""
        self.update_position_tracking()

    @classmethod
    def from_domains(cls, domains: List[Domain], sequence_length: int) -> 'DomainLayout':
        """Create layout from a list of domains"""
        layout = cls(sequence_length=sequence_length, domains=domains.copy())
        layout.analyze_gaps()
        return layout

    def update_position_tracking(self) -> None:
        """Update used/unused position sets based on current domains"""
        self.used_positions = set()
        for domain in self.domains:
            self.used_positions.update(domain.assigned_positions)

        self.unused_positions = set(range(1, self.sequence_length + 1)) - self.used_positions

    def analyze_gaps(self, min_domain_size: int = 25) -> None:
        """Find and classify all unassigned segments"""
        self.update_position_tracking()
        self.unassigned_segments = []

        if not self.unused_positions:
            return  # No gaps to analyze

        # Group consecutive unused positions into segments
        sorted_unused = sorted(self.unused_positions)
        current_segment_start = sorted_unused[0]
        current_segment_end = sorted_unused[0]

        for pos in sorted_unused[1:]:
            if pos == current_segment_end + 1:
                # Consecutive position - extend current segment
                current_segment_end = pos
            else:
                # Gap - finish current segment and start new one
                segment = self._create_segment(current_segment_start, current_segment_end, min_domain_size)
                self.unassigned_segments.append(segment)
                current_segment_start = pos
                current_segment_end = pos

        # Don't forget the last segment
        segment = self._create_segment(current_segment_start, current_segment_end, min_domain_size)
        self.unassigned_segments.append(segment)

        # Classify all segments
        for segment in self.unassigned_segments:
            segment.segment_type = segment.classify_type(self.domains, self.sequence_length)
            segment.neighboring_domains = segment.find_neighbors(self.domains)

    def _create_segment(self, start: int, end: int, min_domain_size: int) -> UnassignedSegment:
        """Create an unassigned segment with proper classification"""
        length = end - start + 1
        fragment_size = FragmentSize.SMALL if length < min_domain_size else FragmentSize.LARGE

        return UnassignedSegment(
            start=start,
            end=end,
            fragment_size=fragment_size
        )

    def add_domain(self, domain: Domain) -> None:
        """Add a domain to the layout"""
        self.domains.append(domain)
        self.update_position_tracking()

    def remove_domain(self, domain: Domain) -> None:
        """Remove a domain from the layout"""
        if domain in self.domains:
            self.domains.remove(domain)
            self.update_position_tracking()

    def merge_segment_with_domain(self, segment: UnassignedSegment, domain: Domain) -> None:
        """Merge an unassigned segment with a domain (with provenance tracking)"""
        # Record the optimization action
        domain.record_optimization_action(
            "merge_segment",
            f"{segment.segment_type.value}_{segment.start}-{segment.end}"
        )

        domain.add_positions(segment.positions)
        self.update_position_tracking()

        # Remove the segment from unassigned list
        if segment in self.unassigned_segments:
            self.unassigned_segments.remove(segment)

    def split_segment_between_domains(self, segment: UnassignedSegment,
                                    domain1: Domain, domain2: Domain,
                                    split_positions: Tuple[Set[int], Set[int]]) -> None:
        """Split a segment between two domains (with provenance tracking)"""
        positions_to_domain1, positions_to_domain2 = split_positions

        # Record optimization actions
        if positions_to_domain1:
            domain1.record_optimization_action(
                "split_merge",
                f"gained_{len(positions_to_domain1)}_from_{segment.start}-{segment.end}"
            )

        if positions_to_domain2:
            domain2.record_optimization_action(
                "split_merge",
                f"gained_{len(positions_to_domain2)}_from_{segment.start}-{segment.end}"
            )

        # Add positions to respective domains
        domain1.add_positions(positions_to_domain1)
        domain2.add_positions(positions_to_domain2)

        self.update_position_tracking()

        # Remove the segment from unassigned list
        if segment in self.unassigned_segments:
            self.unassigned_segments.remove(segment)

    def get_coverage_stats(self) -> dict:
        """Get coverage statistics for this layout"""
        return {
            'total_residues': self.sequence_length,
            'assigned_residues': len(self.used_positions),
            'unassigned_residues': len(self.unused_positions),
            'coverage_percent': len(self.used_positions) / self.sequence_length * 100 if self.sequence_length > 0 else 0,
            'num_domains': len(self.domains),
            'num_gaps': len(self.unassigned_segments),
            'small_fragments': len([s for s in self.unassigned_segments if s.is_small_fragment]),
            'large_gaps': len([s for s in self.unassigned_segments if not s.is_small_fragment])
        }

    def get_overlapping_domains(self, tolerance: int = 0) -> List[Tuple[Domain, Domain]]:
        """Find pairs of domains that overlap"""
        overlapping_pairs = []

        for i, domain1 in enumerate(self.domains):
            for domain2 in self.domains[i+1:]:
                if domain1.overlaps_with(domain2, tolerance):
                    overlapping_pairs.append((domain1, domain2))

        return overlapping_pairs

    def resolve_small_overlaps(self, max_overlap_size: int = 5) -> int:
        """Resolve small overlaps between domains by trimming (with provenance tracking)"""
        overlaps_resolved = 0

        overlapping_pairs = self.get_overlapping_domains()

        for domain1, domain2 in overlapping_pairs:
            overlap_positions = domain1.assigned_positions.intersection(domain2.assigned_positions)

            if len(overlap_positions) <= max_overlap_size:
                # Small overlap - resolve by giving to the domain with better evidence
                if domain1.confidence > domain2.confidence:
                    domain2.remove_positions(overlap_positions)
                    domain2.record_optimization_action(
                        "overlap_trim",
                        f"lost_{len(overlap_positions)}_to_{domain1.id}"
                    )
                elif domain2.confidence > domain1.confidence:
                    domain1.remove_positions(overlap_positions)
                    domain1.record_optimization_action(
                        "overlap_trim",
                        f"lost_{len(overlap_positions)}_to_{domain2.id}"
                    )
                else:
                    # Equal confidence - give to the larger domain
                    if len(domain1.assigned_positions) >= len(domain2.assigned_positions):
                        domain2.remove_positions(overlap_positions)
                        domain2.record_optimization_action(
                            "overlap_trim",
                            f"lost_{len(overlap_positions)}_to_larger_{domain1.id}"
                        )
                    else:
                        domain1.remove_positions(overlap_positions)
                        domain1.record_optimization_action(
                            "overlap_trim",
                            f"lost_{len(overlap_positions)}_to_larger_{domain2.id}"
                        )

                overlaps_resolved += 1

        if overlaps_resolved > 0:
            self.update_position_tracking()

        return overlaps_resolved

    def get_quality_summary(self) -> Dict[str, Any]:
        """Get quality summary for all domains in this layout"""
        if not self.domains:
            return {'total_domains': 0, 'quality_distribution': {}}

        quality_counts = {'good': 0, 'questionable': 0, 'poor': 0, 'unknown': 0}
        domain_assessments = []

        for domain in self.domains:
            assessment = domain.get_quality_assessment()
            domain_assessments.append(assessment)

            overall = assessment.get('overall_assessment', 'unknown')
            quality_counts[overall] = quality_counts.get(overall, 0) + 1

        return {
            'total_domains': len(self.domains),
            'quality_distribution': quality_counts,
            'domain_assessments': domain_assessments,
            'quality_percentage': {
                quality: (count / len(self.domains)) * 100
                for quality, count in quality_counts.items()
            }
        }


# PROVENANCE TRACKING: Partition-level metadata container
@dataclass
class PartitionMetadata:
    """Metadata for an entire domain partition with comprehensive provenance"""
    pdb_id: str
    chain_id: str

    # Algorithm versioning
    algorithm_version: Optional[str] = None
    git_commit_hash: Optional[str] = None
    processing_timestamp: Optional[datetime] = None

    # File provenance
    source_domain_summary_path: Optional[str] = None
    source_domain_summary_hash: Optional[str] = None
    output_xml_path: Optional[str] = None
    output_xml_hash: Optional[str] = None

    # Processing context
    sequence_length: Optional[int] = None
    batch_id: Optional[str] = None
    process_parameters: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.processing_timestamp is None:
            self.processing_timestamp = datetime.now()


# Add convenience property to Evidence for compatibility
@property
def evidence_confidence(self) -> float:
    """Get confidence value for overlap resolution"""
    return self.confidence

# Extend Evidence class
Evidence.confidence_value = evidence_confidence
