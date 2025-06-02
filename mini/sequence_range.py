# mini_pyecod/sequence_range.py
"""Minimal but correct sequence range handling"""

from dataclasses import dataclass
from typing import List, Tuple, Set

@dataclass
class SequenceRange:
    """Minimal sequence range for domain analysis"""
    segments: List[Tuple[int, int]]
    
    @classmethod
    def parse(cls, range_str: str) -> 'SequenceRange':
        """Parse range string like '2-248,491-517' or '253-499'"""
        segments = []
        for segment in range_str.split(','):
            segment = segment.strip()
            if '-' in segment:
                parts = segment.split('-', 1)
                start, end = int(parts[0]), int(parts[1])
                if start <= end:
                    segments.append((start, end))
        
        return cls(segments=segments)
    
    @property
    def is_discontinuous(self) -> bool:
        """Check if range has multiple segments"""
        return len(self.segments) > 1
    
    @property
    def span(self) -> Tuple[int, int]:
        """Get overall start and end (ignoring gaps)"""
        if not self.segments:
            return (0, 0)
        return (self.segments[0][0], self.segments[-1][1])
    
    @property
    def size(self) -> int:
        """Total residues covered"""
        return sum(end - start + 1 for start, end in self.segments)
    
    def get_positions(self) -> Set[int]:
        """Get all positions as a set"""
        positions = set()
        for start, end in self.segments:
            positions.update(range(start, end + 1))
        return positions
    
    def overlaps(self, other: 'SequenceRange') -> bool:
        """Check if ranges overlap"""
        my_pos = self.get_positions()
        other_pos = other.get_positions()
        return bool(my_pos & other_pos)
    
    def __str__(self) -> str:
        """String representation"""
        return ','.join(f"{s}-{e}" for s, e in self.segments)
