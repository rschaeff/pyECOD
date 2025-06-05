# mini/core/sequence_range.py
"""
Core sequence range handling for protein domain analysis.

Copied from ecod.core.sequence_range to achieve mini isolation.
Unified implementation that replaces the scattered range parsing throughout
the codebase. Based on the battle-tested Perl Range.pm module but with
modern Python design principles.
"""

from dataclasses import dataclass
from typing import List, Tuple, Optional, Union, Set
import re
import logging

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class SequenceSegment:
    """
    Single continuous segment of a sequence range.

    Immutable to prevent accidental modification. Can represent either
    single-chain (chain=None) or multi-chain segments.
    """
    start: int
    end: int
    chain: Optional[str] = None

    def __post_init__(self):
        if self.start > self.end:
            raise ValueError(f"Invalid segment: start {self.start} > end {self.end}")
        if self.start < 1:
            raise ValueError(f"Invalid segment: start {self.start} < 1")

    @property
    def length(self) -> int:
        """Number of residues in this segment"""
        return self.end - self.start + 1

    def contains(self, position: int, chain: Optional[str] = None) -> bool:
        """Check if position is within this segment"""
        position_match = self.start <= position <= self.end
        chain_match = (self.chain is None and chain is None) or (self.chain == chain)
        return position_match and chain_match

    def overlaps(self, other: 'SequenceSegment') -> bool:
        """Check if this segment overlaps with another"""
        # Different chains don't overlap
        if self.chain != other.chain:
            return False

        # Check position overlap
        return not (self.end < other.start or other.end < self.start)

    def gap_to(self, other: 'SequenceSegment') -> Optional[int]:
        """
        Calculate gap between this segment and another.
        Returns None if different chains or overlapping.
        Returns gap size (positive = gap, 0 = adjacent, negative = overlap)
        """
        if self.chain != other.chain or self.overlaps(other):
            return None

        # Calculate gap (positive = space between, negative = overlap)
        if self.end < other.start:
            return other.start - self.end - 1
        else:
            return self.start - other.end - 1

    def merge_with(self, other: 'SequenceSegment') -> 'SequenceSegment':
        """Merge this segment with another (must be same chain)"""
        if self.chain != other.chain:
            raise ValueError(f"Cannot merge segments from different chains: {self.chain} vs {other.chain}")

        return SequenceSegment(
            start=min(self.start, other.start),
            end=max(self.end, other.end),
            chain=self.chain
        )

    def __str__(self) -> str:
        """String representation for debugging"""
        if self.chain:
            return f"{self.chain}:{self.start}-{self.end}"
        else:
            return f"{self.start}-{self.end}"


class SequenceRange:
    """
    Unified sequence range handling for protein domain analysis.

    Replaces all the scattered range parsing and arithmetic throughout
    the codebase. Handles single-chain, multi-chain, and discontinuous ranges.
    """

    def __init__(self, segments: List[SequenceSegment]):
        """
        Create a sequence range from segments.

        Args:
            segments: List of SequenceSegment objects
        """
        if not segments:
            raise ValueError("SequenceRange requires at least one segment")

        self.segments = self._validate_and_sort(segments)

    def _validate_and_sort(self, segments: List[SequenceSegment]) -> List[SequenceSegment]:
        """Validate segments and sort by chain then position"""
        if not segments:
            return []

        # Check chain consistency
        chains = {seg.chain for seg in segments}
        if None in chains and len(chains) > 1:
            raise ValueError("Cannot mix single-chain (None) and multi-chain segments")

        # Sort by chain then start position
        return sorted(segments, key=lambda seg: (seg.chain or '', seg.start))

    @classmethod
    def parse(cls, range_str: str) -> 'SequenceRange':
        """
        Parse any range string format into SequenceRange.

        Supports:
        - Single chain: "1-10,15-20"
        - Multi-chain: "A:1-10,B:15-20"
        - Discontinuous: "2-135,207-248,491-517"
        - Mixed: "A:1-10,15-20" (not recommended but handled)

        Args:
            range_str: Range string to parse

        Returns:
            SequenceRange object

        Raises:
            ValueError: If range string is invalid
        """
        if not range_str or not range_str.strip():
            raise ValueError("Empty range string")

        range_str = range_str.strip()
        segments = []

        # Split on commas, handle each segment
        for segment_str in range_str.split(','):
            segment_str = segment_str.strip()
            if not segment_str:
                continue

            # Try multi-chain format first: "A:1-10" or "A:5"
            if ':' in segment_str:
                chain_match = re.match(r'^([A-Za-z0-9_]+):(.+)$', segment_str)
                if not chain_match:
                    raise ValueError(f"Invalid chain format: {segment_str}")

                chain = chain_match.group(1)
                range_part = chain_match.group(2)

                # Parse the range part
                if '-' in range_part:
                    parts = range_part.split('-', 1)  # Only split on first dash
                    if len(parts) != 2:
                        raise ValueError(f"Invalid range segment: {segment_str}")
                    try:
                        start, end = int(parts[0]), int(parts[1])
                        segments.append(SequenceSegment(start, end, chain))
                    except ValueError as e:
                        raise ValueError(f"Invalid range numbers in {segment_str}: {e}")
                else:
                    # Single position: "A:5"
                    try:
                        pos = int(range_part)
                        segments.append(SequenceSegment(pos, pos, chain))
                    except ValueError as e:
                        raise ValueError(f"Invalid position in {segment_str}: {e}")

            # Try single-chain format: "1-10" or "5"
            else:
                if '-' in segment_str:
                    parts = segment_str.split('-', 1)  # Only split on first dash
                    if len(parts) != 2:
                        raise ValueError(f"Invalid range segment: {segment_str}")
                    try:
                        start, end = int(parts[0]), int(parts[1])
                        segments.append(SequenceSegment(start, end, None))
                    except ValueError as e:
                        raise ValueError(f"Invalid range numbers in {segment_str}: {e}")
                else:
                    # Single position: "5"
                    try:
                        pos = int(segment_str)
                        segments.append(SequenceSegment(pos, pos, None))
                    except ValueError as e:
                        raise ValueError(f"Invalid position in {segment_str}: {e}")

        if not segments:
            raise ValueError(f"No valid segments found in range string: {range_str}")

        return cls(segments)

    @classmethod
    def from_positions(cls, positions: List[int], chain: Optional[str] = None) -> 'SequenceRange':
        """
        Create SequenceRange from list of positions (like Perl rangify).

        Args:
            positions: List of sequence positions
            chain: Optional chain identifier

        Returns:
            SequenceRange with segments for continuous runs
        """
        if not positions:
            raise ValueError("Empty positions list")

        # Sort and deduplicate
        unique_positions = sorted(set(positions))
        segments = []

        start = unique_positions[0]
        end = unique_positions[0]

        for pos in unique_positions[1:]:
            if pos == end + 1:
                # Continuous, extend current segment
                end = pos
            else:
                # Gap found, close current segment and start new one
                segments.append(SequenceSegment(start, end, chain))
                start = end = pos

        # Add final segment
        segments.append(SequenceSegment(start, end, chain))

        return cls(segments)

    def to_positions(self) -> List[Tuple[int, Optional[str]]]:
        """
        Expand range to individual positions (like Perl range_expand).

        Returns:
            List of (position, chain) tuples
        """
        positions = []
        for segment in self.segments:
            for pos in range(segment.start, segment.end + 1):
                positions.append((pos, segment.chain))
        return positions

    def to_positions_simple(self) -> List[int]:
        """
        Expand to positions only (for single-chain ranges).

        Returns:
            List of positions

        Raises:
            ValueError: If this is a multi-chain range
        """
        if self.is_multi_chain:
            raise ValueError("Cannot convert multi-chain range to simple positions")

        positions = []
        for segment in self.segments:
            positions.extend(range(segment.start, segment.end + 1))
        return positions

    @property
    def is_multi_chain(self) -> bool:
        """Check if this range involves multiple chains"""
        chains = {seg.chain for seg in self.segments}
        return len(chains) > 1 or (len(chains) == 1 and None not in chains)

    @property
    def is_discontinuous(self) -> bool:
        """Check if this range has gaps (multiple segments)"""
        return len(self.segments) > 1

    @property
    def chains(self) -> Set[Optional[str]]:
        """Get all chains involved in this range"""
        return {seg.chain for seg in self.segments}

    @property
    def total_length(self) -> int:
        """Total number of residues across all segments"""
        return sum(seg.length for seg in self.segments)

    @property
    def span(self) -> Tuple[int, int]:
        """Overall start and end positions (ignoring gaps)"""
        return (self.segments[0].start, self.segments[-1].end)

    @property
    def start_position(self) -> int:
        """First position in range"""
        return self.segments[0].start

    @property
    def end_position(self) -> int:
        """Last position in range"""
        return self.segments[-1].end

    def coverage(self, other: 'SequenceRange') -> float:
        """
        Calculate coverage of this range by another (like Perl region_coverage).

        Args:
            other: Range to calculate coverage with

        Returns:
            Fraction of this range covered by other (0.0 to 1.0)
        """
        if self.total_length == 0:
            return 0.0

        # Get all positions from both ranges
        self_positions = set(self.to_positions())
        other_positions = set(other.to_positions())

        # Calculate intersection
        intersection = self_positions & other_positions

        return len(intersection) / len(self_positions)

    def overlaps(self, other: 'SequenceRange') -> bool:
        """Check if this range overlaps with another"""
        for self_seg in self.segments:
            for other_seg in other.segments:
                if self_seg.overlaps(other_seg):
                    return True
        return False

    def merge_gaps(self, gap_tolerance: int) -> 'SequenceRange':
        """
        Merge segments separated by small gaps (like Perl ungap_range).

        Args:
            gap_tolerance: Maximum gap size to merge

        Returns:
            New SequenceRange with gaps merged
        """
        if gap_tolerance < 0:
            raise ValueError("Gap tolerance must be non-negative")

        if len(self.segments) <= 1:
            return SequenceRange(self.segments[:])  # Return copy

        # Group segments by chain
        chain_segments = {}
        for segment in self.segments:
            chain = segment.chain
            if chain not in chain_segments:
                chain_segments[chain] = []
            chain_segments[chain].append(segment)

        # Merge gaps within each chain
        merged_segments = []
        for chain, segments in chain_segments.items():
            # Sort segments by start position
            sorted_segments = sorted(segments, key=lambda s: s.start)

            if not sorted_segments:
                continue

            current = sorted_segments[0]

            for next_seg in sorted_segments[1:]:
                gap = current.gap_to(next_seg)

                if gap is not None and gap <= gap_tolerance:
                    # Merge segments
                    current = current.merge_with(next_seg)
                    logger.debug(f"Merged segments with gap {gap}: {current}")
                else:
                    # Gap too large, keep segments separate
                    merged_segments.append(current)
                    current = next_seg

            # Add final segment
            merged_segments.append(current)

        return SequenceRange(merged_segments)

    def map_to_structured(self, structured_positions: List[int],
                         chain: Optional[str] = None) -> 'SequenceRange':
        """
        Map range to structured residues only (like Perl struct_region).

        Args:
            structured_positions: List of structured sequence positions
            chain: Chain identifier for multi-chain ranges

        Returns:
            New SequenceRange containing only structured positions
        """
        if self.is_multi_chain and chain is None:
            raise ValueError("Must specify chain for multi-chain range mapping")

        structured_set = set(structured_positions)
        mapped_positions = []

        for pos, pos_chain in self.to_positions():
            # For single-chain ranges, ignore chain
            if not self.is_multi_chain:
                if pos in structured_set:
                    mapped_positions.append(pos)
            # For multi-chain ranges, match chain
            elif pos_chain == chain and pos in structured_set:
                mapped_positions.append(pos)

        if not mapped_positions:
            raise ValueError("No structured positions found in range")

        return SequenceRange.from_positions(mapped_positions, chain)

    def __str__(self) -> str:
        """
        String representation (like Perl rangify/multi_chain_rangify).

        Returns:
            Compact range string representation
        """
        if not self.segments:
            return ""

        # Single segment optimization
        if len(self.segments) == 1:
            return str(self.segments[0])

        # Multiple segments
        return ",".join(str(seg) for seg in self.segments)

    def __repr__(self) -> str:
        """Developer representation"""
        return f"SequenceRange({self.segments!r})"

    def __eq__(self, other) -> bool:
        """Equality comparison"""
        if not isinstance(other, SequenceRange):
            return False
        return self.segments == other.segments

    def __len__(self) -> int:
        """Total length of range"""
        return self.total_length

    def __contains__(self, item) -> bool:
        """Check if position or (position, chain) is in range"""
        if isinstance(item, tuple) and len(item) == 2:
            pos, chain = item
            return any(seg.contains(pos, chain) for seg in self.segments)
        elif isinstance(item, int):
            # For single-chain ranges only
            if self.is_multi_chain:
                raise ValueError("Use (position, chain) tuple for multi-chain ranges")
            return any(seg.contains(item) for seg in self.segments)
        else:
            return False


# Convenience functions for backward compatibility and ease of use

def parse_range(range_str: str) -> SequenceRange:
    """Convenience function for parsing ranges"""
    return SequenceRange.parse(range_str)


def positions_to_range(positions: List[int], chain: Optional[str] = None) -> SequenceRange:
    """Convenience function for creating ranges from positions"""
    return SequenceRange.from_positions(positions, chain)


def calculate_coverage(range1: SequenceRange, range2: SequenceRange) -> float:
    """Convenience function for coverage calculation"""
    return range1.coverage(range2)
