# New module: ecod/utils/evidence_processing.py

from typing import List, Dict, Any, Optional, Tuple, Set
import xml.etree.ElementTree as ET
from ecod.models.pipeline import BlastHit, HHSearchHit, DomainSummaryModel

def process_blast_xml(xml_path: str, hit_type: str = "domain_blast") -> List[BlastHit]:
    """
    Process BLAST XML into BlastHit models
    
    This consolidates functionality from:
    - DomainSummary._process_blast()
    - DomainSummary._process_chain_blast_to_dict()
    - HHSearchProcessor's collator._extract_boundaries_from_blast()
    """
    hits = []
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        
        # Get query info
        query_len_str = root.findtext(".//BlastOutput_query-len", "0")
        query_length = int(query_len_str) if query_len_str else 0
        
        # Determine correct XML path based on hit type
        hit_path = ".//chain_blast_run/hits/hit" if hit_type == "chain_blast" else ".//blast_run/hits/hit"
        
        for hit_elem in root.findall(hit_path):
            # Create hit and process details
            hit = BlastHit.from_xml(hit_elem)
            hit.hit_type = hit_type
            hits.append(hit)
            
        return hits
    except Exception as e:
        logging.error(f"Error processing BLAST XML {xml_path}: {str(e)}")
        return []

def process_hhsearch_xml(xml_path: str) -> List[HHSearchHit]:
    """
    Process HHSearch XML into HHSearchHit models
    
    This consolidates functionality from:
    - DomainSummary._process_hhsearch_to_dict()
    - HHSearchProcessor's collator._extract_boundaries_from_hhsearch()
    """
    hits = []
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        
        hit_elements = root.findall(".//hh_hit_list/hh_hit")
        
        for hit_elem in hit_elements:
            hit = HHSearchHit.from_xml(hit_elem)
            hits.append(hit)
            
        return hits
    except Exception as e:
        logging.error(f"Error processing HHSearch XML {xml_path}: {str(e)}")
        return []

def merge_domain_boundaries(boundaries: List[Tuple[int, int, List[str]]]) -> List[Tuple[int, int, List[str]]]:
    """
    Merge overlapping domain boundaries
    
    This consolidates identical implementations from both modules
    """
    if not boundaries:
        return []
        
    # Sort boundaries by start position
    sorted_boundaries = sorted(boundaries, key=lambda x: x[0])
    
    merged = []
    current = sorted_boundaries[0]
    
    for next_boundary in sorted_boundaries[1:]:
        # Check if boundaries overlap
        if next_boundary[0] <= current[1]:
            # Merge boundaries
            current = (
                current[0],
                max(current[1], next_boundary[1]),
                list(set(current[2] + next_boundary[2]))
            )
        else:
            # No overlap, add current to merged list and update current
            merged.append(current)
            current = next_boundary
    
    # Add the last boundary
    merged.append(current)
    
    return merged
