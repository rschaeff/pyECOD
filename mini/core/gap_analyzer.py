# mini/core/gap_analyzer.py (simplified with SequenceRange integration)
"""
Gap analysis for boundary optimization

Simplified implementation using SequenceRange for all position operations.
Pure sequence-based analysis following mini isolation rules.
"""

from typing import List, Set, Tuple
from .models import Domain, UnassignedSegment, DomainLayout, SegmentType, FragmentSize
from .sequence_range import SequenceRange


def find_unassigned_segments(domains: List[Domain], sequence_length: int, 
                           min_domain_size: int = 25) -> List[UnassignedSegment]:
    """
    Find all unassigned segments between domains.
    
    Simplified implementation using SequenceRange functionality.
    
    Args:
        domains: List of assigned domains
        sequence_length: Total protein sequence length
        min_domain_size: Threshold for small vs large segments
        
    Returns:
        List of UnassignedSegment objects
    """
    # Collect all used positions
    used_positions = set()
    for domain in domains:
        used_positions.update(domain.get_positions())
    
    # Find unused positions
    all_positions = set(range(1, sequence_length + 1))
    unused_positions = all_positions - used_positions
    
    if not unused_positions:
        return []
    
    # Create SequenceRange from unused positions, which automatically groups consecutive positions
    unused_range = SequenceRange.from_positions(sorted(unused_positions))
    
    # Convert each segment to UnassignedSegment
    segments = []
    for segment in unused_range.segments:
        positions = list(range(segment.start, segment.end + 1))
        unassigned_segment = UnassignedSegment.from_positions(positions, min_domain_size)
        segments.append(unassigned_segment)
    
    return segments


def get_gap_analysis_summary(segments: List[UnassignedSegment], sequence_length: int) -> dict:
    """
    Get summary statistics for gap analysis.
    
    Args:
        segments: List of unassigned segments
        sequence_length: Total protein sequence length
        
    Returns:
        Dictionary with gap analysis statistics
    """
    if not segments:
        return {
            'total_unassigned': 0,
            'coverage_percent': 100.0,
            'num_gaps': 0,
            'small_fragments': 0,
            'large_gaps': 0,
            'singletons': 0,
            'nterm_gaps': 0,
            'cterm_gaps': 0,
            'interstitial_gaps': 0
        }
    
    total_unassigned = sum(segment.length for segment in segments)
    coverage_percent = ((sequence_length - total_unassigned) / sequence_length * 100) if sequence_length > 0 else 0
    
    # Count by size
    small_fragments = len([s for s in segments if s.is_small_fragment])
    large_gaps = len([s for s in segments if not s.is_small_fragment])
    singletons = len([s for s in segments if s.is_singleton])
    
    # Count by type
    nterm_gaps = len([s for s in segments if s.segment_type == SegmentType.NTERM])
    cterm_gaps = len([s for s in segments if s.segment_type == SegmentType.CTERM])
    interstitial_gaps = len([s for s in segments if s.segment_type == SegmentType.INTERSTITIAL])
    
    return {
        'total_unassigned': total_unassigned,
        'coverage_percent': coverage_percent,
        'num_gaps': len(segments),
        'small_fragments': small_fragments,
        'large_gaps': large_gaps,
        'singletons': singletons,
        'nterm_gaps': nterm_gaps,
        'cterm_gaps': cterm_gaps,
        'interstitial_gaps': interstitial_gaps,
        'largest_gap': max(segment.length for segment in segments),
        'avg_gap_size': total_unassigned / len(segments)
    }


def validate_gap_analysis(segments: List[UnassignedSegment], domains: List[Domain], 
                         sequence_length: int) -> List[str]:
    """
    Validate gap analysis results and return any issues found.
    
    Args:
        segments: List of unassigned segments
        domains: List of assigned domains
        sequence_length: Total protein sequence length
        
    Returns:
        List of validation error messages (empty if no issues)
    """
    issues = []
    
    # Check that segments don't overlap with domains
    for segment in segments:
        for domain in domains:
            if segment.range.overlaps(domain.range):
                overlap_positions = segment.positions.intersection(domain.get_positions())
                issues.append(f"Segment {segment.start}-{segment.end} overlaps with domain {domain.id} at positions {sorted(overlap_positions)}")
    
    # Check that segments don't overlap with each other
    for i, seg1 in enumerate(segments):
        for seg2 in segments[i+1:]:
            if seg1.range.overlaps(seg2.range):
                overlap_positions = seg1.positions.intersection(seg2.positions)
                issues.append(f"Segments {seg1.start}-{seg1.end} and {seg2.start}-{seg2.end} overlap at positions {sorted(overlap_positions)}")
    
    # Check that all positions are accounted for
    all_domain_positions = set()
    for domain in domains:
        all_domain_positions.update(domain.get_positions())
    
    all_segment_positions = set()
    for segment in segments:
        all_segment_positions.update(segment.positions)
    
    all_covered = all_domain_positions.union(all_segment_positions)
    expected_positions = set(range(1, sequence_length + 1))
    
    missing = expected_positions - all_covered
    if missing:
        issues.append(f"Missing positions not covered by domains or segments: {sorted(missing)}")
    
    extra = all_covered - expected_positions
    if extra:
        issues.append(f"Extra positions outside sequence range: {sorted(extra)}")
    
    return issues


# Convenience function for complete gap analysis
def perform_complete_gap_analysis(domains: List[Domain], sequence_length: int,
                                min_domain_size: int = 25, tolerance: int = 5,
                                validate: bool = True, verbose: bool = False) -> Tuple[List[UnassignedSegment], dict]:
    """
    Perform complete gap analysis with all steps.
    
    Args:
        domains: List of assigned domains
        sequence_length: Total protein sequence length
        min_domain_size: Size threshold for fragment classification
        tolerance: Distance tolerance for neighbor detection
        validate: Whether to run validation checks
        verbose: Whether to print detailed information
        
    Returns:
        Tuple of (segments, summary_stats)
    """
    # Step 1: Find unassigned segments
    segments = find_unassigned_segments(domains, sequence_length, min_domain_size)
    
    if verbose:
        print(f"Found {len(segments)} unassigned segments")
    
    # Step 2: Classify segment types and find neighbors (done automatically in DomainLayout)
    if verbose:
        for segment in segments:
            print(f"  Segment {segment.start}-{segment.end}: length={segment.length}, "
                  f"{segment.fragment_size.value}")
    
    # Step 3: Get summary statistics
    summary = get_gap_analysis_summary(segments, sequence_length)
    
    # Step 4: Validation (optional)
    if validate:
        issues = validate_gap_analysis(segments, domains, sequence_length)
        if issues:
            print("WARNING: Gap analysis validation issues:")
            for issue in issues:
                print(f"  {issue}")
        summary['validation_issues'] = len(issues)
    
    return segments, summary
