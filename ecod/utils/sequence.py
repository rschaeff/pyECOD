#!/usr/bin/env python3
"""
Sequence utilities for the ECOD pipeline
Functions for working with protein sequences
"""
import re
import logging
import hashlib
from typing import List, Dict, Any, Optional, Set, Tuple, Union

from ecod.exceptions import ValidationError

logger = logging.getLogger("ecod.utils.sequence")

# Valid amino acid characters (including ambiguous)
VALID_AA = set('ACDEFGHIKLMNPQRSTVWYX')

def validate_sequence(sequence: str) -> bool:
    """Validate protein sequence
    
    Args:
        sequence: Protein sequence
        
    Returns:
        True if valid
        
    Raises:
        ValidationError: If validation fails
    """
    if not sequence:
        raise ValidationError("Protein sequence cannot be empty")
    
    # Check for valid amino acid characters
    invalid_chars = set(sequence.upper()) - VALID_AA
    if invalid_chars:
        raise ValidationError(f"Invalid amino acids in sequence: {', '.join(invalid_chars)}")
    
    return True

def calculate_md5(sequence: str) -> str:
    """Calculate MD5 hash of a sequence
    
    Args:
        sequence: Protein sequence
        
    Returns:
        MD5 hash as hexadecimal string
    """
    # Clean sequence (uppercase, strip whitespace)
    clean_seq = re.sub(r'\s+', '', sequence.upper())
    
    # Calculate MD5
    return hashlib.md5(clean_seq.encode('utf-8')).hexdigest()

def parse_range(range_str: str) -> List[Tuple[int, int]]:
    """Parse range string into list of (start, end) tuples
    
    Args:
        range_str: Range string, e.g. "1-100,120-150"
        
    Returns:
        List of (start, end) tuples
        
    Raises:
        ValueError: If range format is invalid
    """
    segments = []
    
    if not range_str:
        return segments
    
    # Split by comma for multi-segment ranges
    for segment in range_str.split(','):
        segment = segment.strip()
        
        # Skip empty segments
        if not segment:
            continue
        
        # Handle single residue case
        if segment.isdigit():
            pos = int(segment)
            segments.append((pos, pos))
            continue
        
        # Handle range
        if '-' in segment:
            start_str, end_str = segment.split('-')
            
            try:
                start = int(start_str.strip())
                end = int(end_str.strip())
                
                if start <= end:
                    segments.append((start, end))
                else:
                    raise ValueError(f"Invalid range: {segment} (start > end)")
            except ValueError:
                raise ValueError(f"Invalid range format: {segment}")
        else:
            raise ValueError(f"Invalid range format: {segment}")
    
    return segments

def format_range(segments: List[Tuple[int, int]]) -> str:
    """Format list of (start, end) tuples into range string
    
    Args:
        segments: List of (start, end) tuples
        
    Returns:
        Range string, e.g. "1-100,120-150"
    """
    # Sort segments by start position
    sorted_segments = sorted(segments, key=lambda x: x[0])
    
    # Format segments
    formatted = []
    for start, end in sorted_segments:
        if start == end:
            formatted.append(str(start))
        else:
            formatted.append(f"{start}-{end}")
    
    return ','.join(formatted)

def get_positions_from_range(range_str: str) -> Set[int]:
    """Get set of positions from range string
    
    Args:
        range_str: Range string, e.g. "1-100,120-150"
        
    Returns:
        Set of positions
        
    Raises:
        ValueError: If range format is invalid
    """
    positions = set()
    segments = parse_range(range_str)
    
    for start, end in segments:
        positions.update(range(start, end + 1))
    
    return positions

def extract_sequence_by_range(sequence: str, range_str: str) -> str:
    """Extract subsequence using range string
    
    Args:
        sequence: Full sequence
        range_str: Range string, e.g. "1-100,120-150"
        
    Returns:
        Extracted subsequence
        
    Raises:
        ValueError: If range format is invalid or out of bounds
    """
    if not sequence:
        raise ValueError("Sequence cannot be empty")
    
    segments = parse_range(range_str)
    if not segments:
        return ""
    
    # Extract sequence segments
    result = []
    for start, end in segments:
        # Adjust indices to 0-based for Python
        start_idx = start - 1
        end_idx = end
        
        # Validate indices
        if start_idx < 0 or end_idx > len(sequence):
            raise ValueError(f"Range {start}-{end} is out of bounds for sequence of length {len(sequence)}")
        
        result.append(sequence[start_idx:end_idx])
    
    # Join segments
    return ''.join(result)

def calculate_coverage(query_range: str, subject_range: str, normalize_by: str = 'min') -> float:
    """Calculate coverage between two ranges
    
    Args:
        query_range: Query range string
        subject_range: Subject range string
        normalize_by: Normalization method ('min', 'max', 'query', 'subject')
        
    Returns:
        Coverage as fraction (0.0-1.0)
    """
    if not query_range or not subject_range:
        return 0.0
    
    # Get positions
    query_positions = get_positions_from_range(query_range)
    subject_positions = get_positions_from_range(subject_range)
    
    # Calculate overlap
    overlap = len(query_positions.intersection(subject_positions))
    
    # Normalize
    if normalize_by == 'min':
        denominator = min(len(query_positions), len(subject_positions))
    elif normalize_by == 'max':
        denominator = max(len(query_positions), len(subject_positions))
    elif normalize_by == 'query':
        denominator = len(query_positions)
    elif normalize_by == 'subject':
        denominator = len(subject_positions)
    else:
        raise ValueError(f"Invalid normalization method: {normalize_by}")
    
    if denominator == 0:
        return 0.0
    
    return overlap / denominator

def merge_ranges(range1: str, range2: str) -> str:
    """Merge two range strings
    
    Args:
        range1: First range string
        range2: Second range string
        
    Returns:
        Merged range string
    """
    # Get positions
    positions = get_positions_from_range(range1).union(get_positions_from_range(range2))
    if not positions:
        return ""
    
    # Convert back to segments
    sorted_positions = sorted(positions)
    segments = []
    
    start = sorted_positions[0]
    end = start
    
    for pos in sorted_positions[1:]:
        if pos == end + 1:
            # Continue current segment
            end = pos
        else:
            # End current segment and start a new one
            segments.append((start, end))
            start = pos
            end = pos
    
    # Add the last segment
    segments.append((start, end))
    
    # Format segments
    return format_range(segments)

def calculate_sequence_identity(seq1: str, seq2: str) -> float:
    """Calculate sequence identity between two sequences
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        
    Returns:
        Sequence identity as fraction (0.0-1.0)
    """
    if not seq1 or not seq2:
        return 0.0
    
    # Ensure sequences are same length
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be the same length")
    
    # Count matching positions
    matches = sum(a == b for a, b in zip(seq1, seq2))
    
    # Calculate identity
    return matches / len(seq1)

def get_sequence_composition(sequence: str) -> Dict[str, float]:
    """Calculate amino acid composition of a sequence
    
    Args:
        sequence: Protein sequence
        
    Returns:
        Dictionary with amino acid frequencies
    """
    if not sequence:
        return {}
    
    # Count amino acids
    counts = {}
    for aa in sequence:
        aa = aa.upper()
        if aa in VALID_AA:
            counts[aa] = counts.get(aa, 0) + 1
    
    # Calculate frequencies
    total = len(sequence)
    return {aa: count / total for aa, count in counts.items()}

def translate_alignment_to_sequence_range(alignment: str, seq_start: int) -> str:
    """Translate alignment to sequence range
    
    Args:
        alignment: Alignment string (may contain gaps '-')
        seq_start: Sequence start position
        
    Returns:
        Range string representing alignment
    """
    if not alignment:
        return ""
    
    # Find sequence positions
    positions = []
    pos = seq_start
    
    for char in alignment:
        if char != '-':
            positions.append(pos)
            pos += 1
        
    # Convert positions to range
    segments = []
    if not positions:
        return ""
    
    start = positions[0]
    end = start
    
    for pos in positions[1:]:
        if pos == end + 1:
            # Continue current segment
            end = pos
        else:
            # End current segment and start a new one
            segments.append((start, end))
            start = pos
            end = pos
    
    # Add the last segment
    segments.append((start, end))
    
    # Format segments
    return format_range(segments)