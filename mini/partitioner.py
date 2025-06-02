# mini_pyecod/partitioner.py
"""Core partitioning algorithm - NO MERGING!"""

from collections import defaultdict
from typing import List
from .models import Evidence, Domain
from .sequence_range import SequenceRange

def partition_domains(evidence_list: List[Evidence], sequence_length: int) -> List[Domain]:
    """
    Simple partitioning: group by protein family, NO MERGING
    """
    # Group evidence by source protein family
    family_groups = defaultdict(list)
    
    for evidence in evidence_list:
        # Use source PDB as family identifier
        family = evidence.source_pdb
        if not family:
            family = "unknown"
        
        family_groups[family].append(evidence)
    
    # Convert each family group to a domain
    domains = []
    domain_num = 1
    
    for family, family_evidence in family_groups.items():
        if not family_evidence:
            continue
        
        # Collect all ranges for this family
        all_segments = []
        for ev in family_evidence:
            all_segments.extend(ev.query_range.segments)
        
        if not all_segments:
            continue
        
        # Merge overlapping segments WITHIN the same family only
        merged_segments = merge_segments(all_segments)
        
        # Create domain
        domain = Domain(
            id=f"d{domain_num}",
            range=SequenceRange(segments=merged_segments),
            family=family,
            evidence_count=len(family_evidence),
            source=family_evidence[0].type,
            evidence_items=family_evidence
        )
        
        domains.append(domain)
        domain_num += 1
    
    return domains

def merge_segments(segments: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge overlapping segments within same family"""
    if not segments:
        return []
    
    # Sort by start position
    sorted_segs = sorted(segments)
    merged = [sorted_segs[0]]
    
    for start, end in sorted_segs[1:]:
        last_start, last_end = merged[-1]
        
        # Check for overlap or adjacency
        if start <= last_end + 1:
            # Merge
            merged[-1] = (last_start, max(last_end, end))
        else:
            # Keep separate
            merged.append((start, end))
    
    return merged
