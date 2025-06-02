# mini/decomposer.py
"""Chain BLAST decomposition using alignment-based mapping"""

from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
import csv
from .models import Evidence
from ecod.core.sequence_range import SequenceRange, SequenceSegment

@dataclass
class DomainReference:
    """Reference domain information from database dump"""
    domain_id: str
    pdb_id: str
    chain_id: str
    range: SequenceRange
    length: int
    t_group: Optional[str] = None
    h_group: Optional[str] = None

def load_domain_definitions(csv_path: str) -> Dict[Tuple[str, str], List[DomainReference]]:
    """
    Load domain definitions from database dump CSV

    Expected format:
    domain_id,pdb_id,chain_id,range,length,t_group,h_group
    e2ia4A1,2ia4,A,1-93,93,2.60.40,2.1
    e2ia4A2,2ia4,A,94-259,166,2.60.40,2.1

    Returns:
        Dict mapping (pdb_id, chain_id) -> list of domains
    """
    domains_by_chain = {}

    try:
        with open(csv_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                pdb_id = row['pdb_id'].lower()
                chain_id = row['chain_id']

                # Parse range
                try:
                    range_obj = SequenceRange.parse(row['range'])
                except ValueError as e:
                    print(f"Invalid range for {row['domain_id']}: {row['range']} - {e}")
                    continue

                domain = DomainReference(
                    domain_id=row['domain_id'],
                    pdb_id=pdb_id,
                    chain_id=chain_id,
                    range=range_obj,
                    length=int(row['length']),
                    t_group=row.get('t_group'),
                    h_group=row.get('h_group')
                )

                key = (pdb_id, chain_id)
                if key not in domains_by_chain:
                    domains_by_chain[key] = []
                domains_by_chain[key].append(domain)

        # Sort domains by start position
        for key in domains_by_chain:
            domains_by_chain[key].sort(key=lambda d: d.range.segments[0].start)

        print(f"Loaded domain definitions for {len(domains_by_chain)} chains")

    except FileNotFoundError:
        print(f"Warning: Domain definitions file not found: {csv_path}")
    except Exception as e:
        print(f"Error loading domain definitions: {e}")

    return domains_by_chain

def build_alignment_mapping(query_str: str, hit_str: str,
                           query_start: int, hit_start: int) -> Dict[int, int]:
    """
    Build position mapping from query to hit using alignment strings

    Args:
        query_str: Query alignment string (with gaps as '-')
        hit_str: Hit alignment string (with gaps as '-')
        query_start: Starting position in query sequence
        hit_start: Starting position in hit sequence

    Returns:
        Dict mapping query position -> hit position
    """
    if len(query_str) != len(hit_str):
        raise ValueError(f"Alignment strings have different lengths: {len(query_str)} vs {len(hit_str)}")

    mapping = {}
    query_pos = query_start - 1  # Convert to 0-based
    hit_pos = hit_start - 1

    for q_char, h_char in zip(query_str, hit_str):
        # Advance positions for non-gap characters
        if q_char != '-':
            query_pos += 1
        if h_char != '-':
            hit_pos += 1

        # Record mapping when both are non-gaps
        if q_char != '-' and h_char != '-':
            mapping[query_pos] = hit_pos

    return mapping

def decompose_chain_blast_with_mapping(evidence: Evidence,
                                       hit_query_str: str,
                                       hit_hit_str: str,
                                       query_start: int,
                                       hit_start: int,
                                       domain_refs: List[DomainReference]) -> List[Evidence]:
    """
    Decompose chain BLAST hit using alignment mapping to reference domains

    Args:
        evidence: Chain BLAST evidence
        hit_query_str: Query alignment string from XML
        hit_hit_str: Hit alignment string from XML
        query_start: Query alignment start position
        hit_start: Hit alignment start position
        domain_refs: Reference domain definitions for the hit protein

    Returns:
        List of decomposed evidence corresponding to reference domains
    """
    if not domain_refs:
        return [evidence]  # No decomposition possible

    # Build query -> hit position mapping
    try:
        pos_mapping = build_alignment_mapping(hit_query_str, hit_hit_str, query_start, hit_start)
    except ValueError as e:
        print(f"Failed to build alignment mapping: {e}")
        return [evidence]

    # Invert to get hit -> query mapping
    hit_to_query = {hit: query for query, hit in pos_mapping.items()}

    decomposed = []

    for ref_domain in domain_refs:
        # Get all positions in this reference domain
        ref_positions = set(ref_domain.range.to_positions_simple())

        # Map reference positions to query positions
        query_positions = []
        for ref_pos in sorted(ref_positions):
            # Convert to 0-based for mapping
            if ref_pos - 1 in hit_to_query:
                query_positions.append(hit_to_query[ref_pos - 1] + 1)  # Convert back to 1-based

        if len(query_positions) < 20:  # Skip tiny mapped regions
            continue

        # Create range from mapped positions
        try:
            query_range = SequenceRange.from_positions(query_positions)
        except ValueError as e:
            print(f"Failed to create range from positions: {e}")
            continue

        # Calculate coverage of reference domain
        coverage = len(query_positions) / len(ref_positions)

        # Create decomposed evidence
        new_evidence = Evidence(
            type="chain_blast_decomposed",
            source_pdb=evidence.source_pdb,
            query_range=query_range,
            confidence=evidence.confidence * coverage,  # Adjust by coverage
            evalue=evidence.evalue,
            domain_id=ref_domain.domain_id,
            t_group=ref_domain.t_group,
            h_group=ref_domain.h_group,
            reference_length=ref_domain.length,
            alignment_coverage=coverage
        )

        decomposed.append(new_evidence)

        print(f"  Decomposed to {ref_domain.domain_id}: {query_range} (coverage={coverage:.1%})")

    return decomposed if decomposed else [evidence]
