#!/usr/bin/env python3
"""
Domain models for the ECOD pipeline
Defines data models for domains and their classifications
"""
from dataclasses import dataclass, field
from datetime import datetime
from typing import Optional, List, Dict, Any, Set, Tuple

from ecod.core.errors import ValidationError

@dataclass
class DomainSequence:
    """Domain sequence model"""
    domain_id: int
    sequence: str
    sequence_length: int
    sequence_md5: str
    id: Optional[int] = None
    original_range: Optional[str] = None
    created_at: Optional[datetime] = None
    
    @classmethod
    def from_db_row(cls, row: Dict[str, Any]) -> 'DomainSequence':
        """Create instance from database row
        
        Args:
            row: Database row as dictionary
            
        Returns:
            DomainSequence instance
        """
        return cls(
            id=row.get('id'),
            domain_id=row.get('domain_id'),
            sequence=row.get('sequence', ''),
            sequence_length=row.get('sequence_length', 0),
            sequence_md5=row.get('sequence_md5', ''),
            original_range=row.get('original_range'),
            created_at=row.get('created_at')
        )
    
    def validate(self) -> bool:
        """Validate sequence
        
        Returns:
            True if valid
            
        Raises:
            ValidationError: If validation fails
        """
        if not self.sequence:
            raise ValidationError("Domain sequence cannot be empty")
        
        # Check for valid amino acid characters
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        invalid_chars = set(self.sequence.upper()) - valid_aa
        if invalid_chars:
            raise ValidationError(f"Invalid amino acids in sequence: {', '.join(invalid_chars)}")
        
        # Check sequence length
        if self.sequence_length != len(self.sequence):
            raise ValidationError(f"Sequence length ({self.sequence_length}) does not match actual length ({len(self.sequence)})")
        
        return True
    
    def to_fasta(self, header: Optional[str] = None) -> str:
        """Convert to FASTA format
        
        Args:
            header: FASTA header (default: domain_id)
            
        Returns:
            FASTA formatted string
        """
        header_text = header or f"domain_{self.domain_id}"
        return f">{header_text}\n{self.sequence}\n"
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary
        
        Returns:
            Dictionary representation
        """
        return {
            'id': self.id,
            'domain_id': self.domain_id,
            'sequence': self.sequence,
            'sequence_length': self.sequence_length,
            'sequence_md5': self.sequence_md5,
            'original_range': self.original_range,
            'created_at': self.created_at.isoformat() if self.created_at else None
        }

@dataclass
class DomainDSSPDetail:
    """Domain DSSP (secondary structure) detail model"""
    domain_id: int
    id: Optional[int] = None
    asa: Optional[float] = None
    helix_residues: Optional[str] = None
    strand_residues: Optional[str] = None
    helix_count: Optional[int] = None
    strand_count: Optional[int] = None
    secondary_structure_string: Optional[str] = None
    hbonds_total: Optional[int] = None
    created_at: Optional[datetime] = None
    
    @classmethod
    def from_db_row(cls, row: Dict[str, Any]) -> 'DomainDSSPDetail':
        """Create instance from database row
        
        Args:
            row: Database row as dictionary
            
        Returns:
            DomainDSSPDetail instance
        """
        return cls(
            id=row.get('id'),
            domain_id=row.get('domain_id'),
            asa=row.get('asa'),
            helix_residues=row.get('helix_residues'),
            strand_residues=row.get('strand_residues'),
            helix_count=row.get('helix_count'),
            strand_count=row.get('strand_count'),
            secondary_structure_string=row.get('secondary_structure_string'),
            hbonds_total=row.get('hbonds_total'),
            created_at=row.get('created_at')
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary
        
        Returns:
            Dictionary representation
        """
        return {
            'id': self.id,
            'domain_id': self.domain_id,
            'asa': self.asa,
            'helix_residues': self.helix_residues,
            'strand_residues': self.strand_residues,
            'helix_count': self.helix_count,
            'strand_count': self.strand_count,
            'secondary_structure_string': self.secondary_structure_string,
            'hbonds_total': self.hbonds_total,
            'created_at': self.created_at.isoformat() if self.created_at else None
        }

@dataclass
class DomainRangeSegment:
    """Domain range segment model"""
    start: int
    end: int
    
    def __str__(self) -> str:
        """String representation
        
        Returns:
            String representation (start-end)
        """
        return f"{self.start}-{self.end}"
    
    @property
    def length(self) -> int:
        """Get segment length
        
        Returns:
            Segment length
        """
        return self.end - self.start + 1
    
    def overlaps(self, other: 'DomainRangeSegment') -> bool:
        """Check if segments overlap
        
        Args:
            other: Another segment
            
        Returns:
            True if segments overlap
        """
        return max(self.start, other.start) <= min(self.end, other.end)
    
    def get_positions(self) -> Set[int]:
        """Get set of positions in this segment
        
        Returns:
            Set of positions
        """
        return set(range(self.start, self.end + 1))
    
    @classmethod
    def from_string(cls, segment_str: str) -> 'DomainRangeSegment':
        """Create segment from string
        
        Args:
            segment_str: String representation (start-end)
            
        Returns:
            DomainRangeSegment instance
            
        Raises:
            ValueError: If string format is invalid
        """
        try:
            start_str, end_str = segment_str.split('-')
            return cls(int(start_str), int(end_str))
        except (ValueError, IndexError):
            raise ValueError(f"Invalid segment format: {segment_str}. Expected start-end")

@dataclass
class DomainRange:
    """Domain range model (can have multiple segments)"""
    segments: List[DomainRangeSegment] = field(default_factory=list)
    
    def __str__(self) -> str:
        """String representation
        
        Returns:
            String representation (comma-separated segments)
        """
        return ','.join(str(segment) for segment in self.segments)
    
    @property
    def length(self) -> int:
        """Get total range length
        
        Returns:
            Total range length
        """
        positions = self.get_positions()
        return len(positions)
    
    def get_positions(self) -> Set[int]:
        """Get set of all positions in this range
        
        Returns:
            Set of positions
        """
        positions = set()
        for segment in self.segments:
            positions.update(segment.get_positions())
        return positions
    
    def overlaps(self, other: 'DomainRange') -> bool:
        """Check if ranges overlap
        
        Args:
            other: Another range
            
        Returns:
            True if ranges overlap
        """
        positions1 = self.get_positions()
        positions2 = other.get_positions()
        return bool(positions1.intersection(positions2))
    
    def overlap_size(self, other: 'DomainRange') -> int:
        """Calculate size of overlap
        
        Args:
            other: Another range
            
        Returns:
            Size of overlap (number of positions)
        """
        positions1 = self.get_positions()
        positions2 = other.get_positions()
        return len(positions1.intersection(positions2))
    
    @classmethod
    def from_string(cls, range_str: str) -> 'DomainRange':
        """Create range from string
        
        Args:
            range_str: String representation (comma-separated segments)
            
        Returns:
            DomainRange instance
            
        Raises:
            ValueError: If string format is invalid
        """
        domain_range = cls()
        if not range_str:
            return domain_range
            
        # Split by comma for multi-segment ranges
        for segment_str in range_str.split(','):
            segment = DomainRangeSegment.from_string(segment_str.strip())
            domain_range.segments.append(segment)
            
        return domain_range
    
    def merge(self, other: 'DomainRange') -> 'DomainRange':
        """Merge with another range
        
        Args:
            other: Another range
            
        Returns:
            New merged DomainRange
        """
        # Get all positions from both ranges
        all_positions = self.get_positions().union(other.get_positions())
        if not all_positions:
            return DomainRange()
            
        # Convert back to segments
        sorted_positions = sorted(all_positions)
        merged_segments = []
        
        current_start = sorted_positions[0]
        current_end = current_start
        
        for pos in sorted_positions[1:]:
            if pos == current_end + 1:
                # Continue current segment
                current_end = pos
            else:
                # End current segment and start a new one
                merged_segments.append(DomainRangeSegment(current_start, current_end))
                current_start = pos
                current_end = pos
                
        # Add the last segment
        merged_segments.append(DomainRangeSegment(current_start, current_end))
        
        # Create and return new range
        merged_range = DomainRange()
        merged_range.segments = merged_segments
        return merged_range

@dataclass
class Domain:
    """Domain model"""
    domain_id: str  # e.g., "e1abcA1"
    range: str  # e.g., "1-100,120-150"
    protein_id: Optional[int] = None
    id: Optional[int] = None
    ecod_uid: Optional[int] = None
    ecod_domain_id: Optional[str] = None
    t_group: Optional[str] = None
    h_group: Optional[str] = None
    x_group: Optional[str] = None
    a_group: Optional[str] = None
    is_manual_rep: bool = False
    is_f70: bool = False
    is_f40: bool = False
    is_f99: bool = False
    hcount: Optional[int] = None
    scount: Optional[int] = None
    length: Optional[int] = None
    chain_id: Optional[str] = None
    asym_id: Optional[str] = None
    created_at: Optional[datetime] = None
    updated_at: Optional[datetime] = None
    
    # Related data
    sequence: Optional[DomainSequence] = None
    dssp: Optional[DomainDSSPDetail] = None
    parsed_range: Optional[DomainRange] = None
    
    def __post_init__(self):
        """Post-initialization processing"""
        # Parse the range string
        if self.range and not self.parsed_range:
            self.parsed_range = DomainRange.from_string(self.range)
    
    @classmethod
    def from_db_row(cls, row: Dict[str, Any]) -> 'Domain':
        """Create instance from database row
        
        Args:
            row: Database row as dictionary
            
        Returns:
            Domain instance
        """
        domain = cls(
            id=row.get('id'),
            ecod_uid=row.get('ecod_uid'),
            protein_id=row.get('protein_id'),
            domain_id=row.get('domain_id', ''),
            ecod_domain_id=row.get('ecod_domain_id'),
            range=row.get('range', ''),
            t_group=row.get('t_group'),
            h_group=row.get('h_group'),
            x_group=row.get('x_group'),
            a_group=row.get('a_group'),
            is_manual_rep=row.get('is_manual_rep', False),
            is_f70=row.get('is_f70', False),
            is_f40=row.get('is_f40', False),
            is_f99=row.get('is_f99', False),
            hcount=row.get('hcount'),
            scount=row.get('scount'),
            length=row.get('length'),
            chain_id=row.get('chain_id'),
            asym_id=row.get('asym_id'),
            created_at=row.get('created_at'),
            updated_at=row.get('updated_at')
        )
        
        return domain
    
    def validate(self) -> bool:
        """Validate domain data
        
        Returns:
            True if valid
            
        Raises:
            ValidationError: If validation fails
        """
        if not self.domain_id:
            raise ValidationError("Domain ID cannot be empty")
        
        if not self.range:
            raise ValidationError("Domain range cannot be empty")
            
        # Validate range format
        try:
            DomainRange.from_string(self.range)
        except ValueError as e:
            raise ValidationError(f"Invalid domain range format: {str(e)}")
        
        # Validate sequence if present
        if self.sequence:
            self.sequence.validate()
            
            # Validate length consistency if both are set
            if self.length is not None and self.sequence.sequence_length != self.length:
                raise ValidationError(f"Domain length ({self.length}) does not match sequence length ({self.sequence.sequence_length})")
        
        return True
    
    def get_positions(self) -> Set[int]:
        """Get set of positions in this domain
        
        Returns:
            Set of positions
        """
        if not self.parsed_range:
            self.parsed_range = DomainRange.from_string(self.range)
            
        return self.parsed_range.get_positions()
    
    def overlaps(self, other: 'Domain') -> bool:
        """Check if domains overlap
        
        Args:
            other: Another domain
            
        Returns:
            True if domains overlap
        """
        positions1 = self.get_positions()
        positions2 = other.get_positions()
        return bool(positions1.intersection(positions2))
    
    def overlap_size(self, other: 'Domain') -> int:
        """Calculate size of overlap with another domain
        
        Args:
            other: Another domain
            
        Returns:
            Size of overlap (number of positions)
        """
        positions1 = self.get_positions()
        positions2 = other.get_positions()
        return len(positions1.intersection(positions2))
    
    def overlap_percentage(self, other: 'Domain') -> float:
        """Calculate percentage of overlap with another domain
        
        Args:
            other: Another domain
            
        Returns:
            Percentage of overlap (0-100)
        """
        overlap_size = self.overlap_size(other)
        
        # Calculate percentage based on the smaller domain
        my_size = len(self.get_positions())
        other_size = len(other.get_positions())
        
        if my_size == 0 or other_size == 0:
            return 0.0
            
        min_size = min(my_size, other_size)
        return (overlap_size / min_size) * 100.0
    
    def extract_sequence(self, full_sequence: str) -> str:
        """Extract domain sequence from protein sequence
        
        Args:
            full_sequence: Full protein sequence
            
        Returns:
            Domain sequence
            
        Raises:
            ValueError: If extraction fails
        """
        if not full_sequence:
            raise ValueError("Full sequence cannot be empty")
            
        if not self.parsed_range:
            self.parsed_range = DomainRange.from_string(self.range)
            
        # Extract sequence segments
        segments = []
        for segment in self.parsed_range.segments:
            # Adjust indices to 0-based for Python
            start_idx = segment.start - 1
            end_idx = segment.end
            
            # Validate indices
            if start_idx < 0 or end_idx > len(full_sequence):
                raise ValueError(f"Segment {segment} is out of bounds for sequence of length {len(full_sequence)}")
                
            segments.append(full_sequence[start_idx:end_idx])
            
        # Combine segments
        return ''.join(segments)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary
        
        Returns:
            Dictionary representation
        """
        result = {
            'id': self.id,
            'ecod_uid': self.ecod_uid,
            'protein_id': self.protein_id,
            'domain_id': self.domain_id,
            'ecod_domain_id': self.ecod_domain_id,
            'range': self.range,
            't_group': self.t_group,
            'h_group': self.h_group,
            'x_group': self.x_group,
            'a_group': self.a_group,
            'is_manual_rep': self.is_manual_rep,
            'is_f70': self.is_f70,
            'is_f40': self.is_f40,
            'is_f99': self.is_f99,
            'hcount': self.hcount,
            'scount': self.scount,
            'length': self.length,
            'chain_id': self.chain_id,
            'asym_id': self.asym_id,
            'created_at': self.created_at.isoformat() if self.created_at else None,
            'updated_at': self.updated_at.isoformat() if self.updated_at else None
        }
        
        # Include related data if available
        if self.sequence:
            result['sequence'] = self.sequence.to_dict()
            
        if self.dssp:
            result['dssp'] = self.dssp.to_dict()
            
        return result
        
@dataclass
class DomainClassification:
    """Domain classification model"""
    id: Optional[int] = None
    t_id: Optional[str] = None
    h_id: Optional[str] = None
    x_id: Optional[str] = None
    name: Optional[str] = None
    description: Optional[str] = None
    parent_id: Optional[int] = None
    level: Optional[str] = None
    
    @classmethod
    def from_db_row(cls, row: Dict[str, Any]) -> 'DomainClassification':
        """Create instance from database row
        
        Args:
            row: Database row as dictionary
            
        Returns:
            DomainClassification instance
        """
        return cls(
            id=row.get('id'),
            t_id=row.get('t_id'),
            h_id=row.get('h_id'),
            x_id=row.get('x_id'),
            name=row.get('name'),
            description=row.get('description'),
            parent_id=row.get('parent_id'),
            level=row.get('level')
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary
        
        Returns:
            Dictionary representation
        """
        return {
            'id': self.id,
            't_id': self.t_id,
            'h_id': self.h_id,
            'x_id': self.x_id,
            'name': self.name,
            'description': self.description,
            'parent_id': self.parent_id,
            'level': self.level
        }