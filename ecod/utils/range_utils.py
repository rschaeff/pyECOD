# ecod/utils/range_utils.py

from typing import List, Tuple, Set, Optional

def parse_range(range_str: str) -> List[Tuple[int, int]]:
    """Parse range string (e.g. "1-100,150-200") into list of (start, end) tuples"""
    segments = []
    if not range_str:
        return segments
        
    for segment in range_str.split(','):
        if '-' in segment:
            try:
                start, end = map(int, segment.split('-'))
                segments.append((start, end))
            except ValueError:
                pass
    return segments

def range_to_positions(range_str: str) -> Set[int]:
    """Convert range string to set of positions"""
    positions = set()
    for start, end in parse_range(range_str):
        positions.update(range(start, end + 1))
    return positions

def positions_to_range(positions: Set[int]) -> str:
    """Convert set of positions to range string"""
    if not positions:
        return ""
        
    # Sort positions
    sorted_positions = sorted(positions)
    
    # Convert to range segments
    segments = []
    seg_start = sorted_positions[0]
    seg_end = seg_start
    
    for pos in sorted_positions[1:]:
        if pos == seg_end + 1:
            # Continue segment
            seg_end = pos
        else:
            # End segment and start new one
            segments.append(f"{seg_start}-{seg_end}")
            seg_start = pos
            seg_end = pos
            
    # Add final segment
    segments.append(f"{seg_start}-{seg_end}")
    
    return ",".join(segments)

def extract_domain_sequence(full_sequence: str, range_str: str) -> str:
    """Extract domain sequence from full sequence based on range"""
    if not full_sequence or not range_str:
        return ""
        
    result = ""
    for start, end in parse_range(range_str):
        # Adjust indices (1-based to 0-based)
        start_idx = max(0, start - 1)
        end_idx = min(len(full_sequence), end)
        result += full_sequence[start_idx:end_idx]
        
    return result
