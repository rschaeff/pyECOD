# mini/decomposer.py
"""Chain BLAST decomposition using alignment-based mapping"""

import os
from typing import List, Dict, Tuple, Optional, Set
from dataclasses import dataclass
import csv
from .models import Evidence
from .sequence_range import SequenceRange, SequenceSegment

@dataclass
class BlacklistEntry:
    """Blacklisted reference chain entry"""
    pdb_id: str
    chain_id: str
    reason: str
    date_added: str = ""
    added_by: str = ""

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

def load_domain_definitions(csv_path: str, verbose: bool = False,
                           blacklist_path: str = None) -> Dict[Tuple[str, str], List[DomainReference]]:
    """
    Load domain definitions from database dump CSV with blacklist support

    Args:
        csv_path: Path to CSV file with domain definitions
        verbose: Whether to print warnings about invalid data
        blacklist_path: Path to blacklist CSV file (optional)
    """
    # Load blacklist
    blacklist = set()
    if blacklist_path:
        blacklist = load_reference_blacklist(blacklist_path, verbose)

    domains_by_chain = {}
    invalid_count = 0
    blacklisted_count = 0

    try:
        with open(csv_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                pdb_id = row['pdb_id'].lower()
                chain_id = row['chain_id']

                # Check blacklist
                if (pdb_id, chain_id) in blacklist:
                    blacklisted_count += 1
                    if verbose:
                        print(f"    Blacklisted: {pdb_id}_{chain_id}")
                    continue

                if verbose and pdb_id == '2ia4':
                    print(f"    Processing 2ia4 entry: {pdb_id}_{chain_id}")

                # Parse range (existing code continues unchanged...)
                try:
                    range_obj = SequenceRange.parse(row['range'])
                except ValueError as e:
                    if verbose:
                        print(f"Invalid range for {row['domain_id']}: {row['range']} - {e}")
                    invalid_count += 1
                    continue

                # Validate length
                try:
                    length = int(row['length'])
                    if length <= 0:
                        invalid_count += 1
                        if verbose:
                            print(f"Warning: Invalid length {length} for {row['domain_id']}")
                        continue
                except (ValueError, KeyError):
                    invalid_count += 1
                    if verbose:
                        print(f"Warning: Missing/invalid length for {row['domain_id']}")
                    continue

                domain = DomainReference(
                    domain_id=row['domain_id'],
                    pdb_id=pdb_id,
                    chain_id=chain_id,
                    range=range_obj,
                    length=length,
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

        print(f"Loaded domain definitions for {len(domains_by_chain)} chain entries")
        if blacklisted_count > 0:
            print(f"Excluded {blacklisted_count} domains from {len(blacklist)} blacklisted chains")
        if invalid_count > 0:
            print(f"Warning: Skipped {invalid_count} domains with invalid data")

        if verbose:
            ia4_chains = [k for k in domains_by_chain.keys() if k[0] == "2ia4"]
            if ia4_chains:
                print(f"2ia4 chains available: {ia4_chains}")
            else:
                print("No 2ia4 chains found in domain definitions")


    except FileNotFoundError:
        print(f"Warning: Domain definitions file not found: {csv_path}")
    except Exception as e:
        print(f"Error loading domain definitions: {e}")

    return domains_by_chain

def load_reference_blacklist(blacklist_path: str, verbose: bool = False) -> Set[Tuple[str, str]]:
    """
    Load reference blacklist from CSV file

    Args:
        blacklist_path: Path to blacklist CSV file
        verbose: Whether to print blacklist information

    Returns:
        Set of (pdb_id, chain_id) tuples to exclude
    """
    blacklist = set()
    entries = []

    if not os.path.exists(blacklist_path):
        if verbose:
            print(f"No blacklist file found: {blacklist_path}")
        return blacklist

    try:
        with open(blacklist_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Skip empty rows
                if not row or not any(row.values()):
                    continue
                
                pdb_id = (row.get('pdb_id') or '').lower().strip()
                chain_id = (row.get('chain_id') or '').strip()
                reason = (row.get('reason') or '').strip()
                date_added = (row.get('date_added') or '').strip()
                added_by = (row.get('added_by') or '').strip()

                if pdb_id and chain_id:
                    blacklist.add((pdb_id, chain_id))
                    entries.append(BlacklistEntry(
                        pdb_id=pdb_id,
                        chain_id=chain_id,
                        reason=reason,
                        date_added=date_added,
                        added_by=added_by
                    ))

        if verbose:
            print(f"Loaded {len(blacklist)} blacklisted reference chains:")
            for entry in entries:
                print(f"  {entry.pdb_id}_{entry.chain_id}: {entry.reason}")

    except Exception as e:
        print(f"Warning: Error loading blacklist from {blacklist_path}: {e}")

    return blacklist

def build_alignment_mapping(query_str: str, hit_str: str,
                           query_start: int, hit_start: int) -> Dict[int, int]:
    """
    Build position mapping from query to hit using alignment strings
    FIXED: Ensure consistent 0-based indexing
    """
    if len(query_str) != len(hit_str):
        raise ValueError(f"Alignment strings have different lengths: {len(query_str)} vs {len(hit_str)}")

    mapping = {}
    query_pos = query_start - 1  # Convert to 0-based
    hit_pos = hit_start - 1      # Convert to 0-based

    for q_char, h_char in zip(query_str, hit_str):
        # Record mapping when both are non-gaps (BEFORE advancing positions)
        if q_char != '-' and h_char != '-':
            mapping[query_pos] = hit_pos

        # Advance positions for non-gap characters
        if q_char != '-':
            query_pos += 1
        if h_char != '-':
            hit_pos += 1

    return mapping

def decompose_chain_blast_with_mapping(evidence: Evidence,
                                       hit_query_str: str,
                                       hit_hit_str: str,
                                       query_start: int,
                                       hit_start: int,
                                       domain_refs: List[DomainReference],
                                       verbose: bool = False) -> List[Evidence]:
    """
    Decompose chain BLAST hit using alignment mapping to reference domains

    Args:
        evidence: Chain BLAST evidence
        hit_query_str: Query alignment string from XML
        hit_hit_str: Hit alignment string from XML
        query_start: Query alignment start position
        hit_start: Hit alignment start position
        domain_refs: Reference domain definitions for the hit protein
        verbose: Whether to print detailed information

    Returns:
        List of decomposed evidence corresponding to reference domains
    """
    if not domain_refs:
        if verbose:
            print(f"  Warning: No domain references for {evidence.source_pdb}")
        return [evidence]  # No decomposition possible



    # Build query -> hit position mapping
    try:
        pos_mapping = build_alignment_mapping(hit_query_str, hit_hit_str, query_start, hit_start)
    except ValueError as e:
        if verbose:
            print(f"Failed to build alignment mapping: {e}")
        return [evidence]

    # Invert to get hit -> query mapping
    hit_to_query = {hit: query for query, hit in pos_mapping.items()}

    decomposed = []
    skipped_count = 0

    for ref_domain in domain_refs:
        # Validate reference domain
        if ref_domain.length <= 0:
            skipped_count += 1
            if verbose:
                print(f"  Warning: Invalid reference length for {ref_domain.domain_id}: {ref_domain.length}")
            continue



        # Get all positions in this reference domain
        ref_positions = set(ref_domain.range.to_positions_simple())

        # Map reference positions to query positions
        query_positions = []
        for ref_pos in sorted(ref_positions):
            # Convert to 0-based for mapping
            if ref_pos - 1 in hit_to_query:
                query_positions.append(hit_to_query[ref_pos - 1] + 1)  # Convert back to 1-based

        if len(query_positions) < 20:  # Skip tiny mapped regions
            if verbose:
                print(f"  Skipping {ref_domain.domain_id}: only {len(query_positions)} positions mapped")
            continue

        # Create range from mapped positions
        try:
            query_range = SequenceRange.from_positions(query_positions)
        except ValueError as e:
            if verbose:
                print(f"Failed to create range from positions: {e}")
            continue

        # Calculate coverage of reference domain (only with valid reference length)
        coverage = len(query_positions) / len(ref_positions)

        if verbose:
            t_group_info = f", T-group={ref_domain.t_group}" if ref_domain.t_group else ""
            print(f"  Decomposed to {ref_domain.domain_id}: {query_range} (coverage={coverage:.1%}{t_group_info})")

        # Create decomposed evidence
        new_evidence = Evidence(
            type="chain_blast_decomposed",
            source_pdb=evidence.source_pdb,
            query_range=query_range,
            confidence=evidence.confidence * coverage,
            evalue=evidence.evalue,
            domain_id=ref_domain.domain_id,

            # NEW PROVENANCE FIELDS:
            source_chain_id=evidence.source_chain_id,
            hit_range=SequenceRange.from_positions(
                sorted([ref_pos for ref_pos in ref_positions if ref_pos - 1 in hit_to_query])
            ),
            hsp_count=evidence.hsp_count,
            discontinuous=query_range.is_discontinuous,
            reference_length=ref_domain.length,
            alignment_coverage=coverage,
            t_group=ref_domain.t_group,
            h_group=ref_domain.h_group
        )

        decomposed.append(new_evidence)

        if verbose:
            print(f"  Decomposed to {ref_domain.domain_id}: {query_range} (coverage={coverage:.1%})")

    if skipped_count > 0 and not verbose:
        print(f"  Skipped {skipped_count} domains with invalid reference lengths")

    return decomposed if decomposed else [evidence]

def decompose_chain_blast_discontinuous(evidence: Evidence,
                                       min_domain: int = 50,
                                       verbose: bool = False) -> List[Evidence]:
    """
    Decompose a discontinuous chain BLAST hit by segments.

    This is a fallback when we don't have alignment data but the hit is discontinuous.
    Each continuous segment becomes a potential domain.

    Note: Without reference domain data, we cannot calculate true coverage or
    set accurate reference lengths.

    6/4/2025 - rds - This method is DEPRECATED and should not be used. AI insanity.

    Args:
        evidence: Chain BLAST evidence with discontinuous range
        min_domain: Minimum domain size
        verbose: Whether to print detailed information

    Returns:
        List of decomposed evidence pieces
    """
    if not evidence.query_range.is_discontinuous:
        return [evidence]

    decomposed = []
    for i, segment in enumerate(evidence.query_range.segments):
        if segment.length >= min_domain:
            # Create new evidence for this segment
            new_evidence = Evidence(
                type="chain_blast_decomposed",
                source_pdb=evidence.source_pdb,
                query_range=SequenceRange([segment]),
                confidence=evidence.confidence * 0.95,  # Small penalty for decomposition
                evalue=evidence.evalue,
                domain_id=f"{evidence.domain_id}_seg{i+1}",
                # Don't set reference_length to segment length - that's query length!
                reference_length=None,  # Unknown without reference data
                alignment_coverage=None,  # Cannot calculate without reference
                t_group=evidence.t_group,  # Preserve if available
                h_group=evidence.h_group   # Preserve if available
            )
            decomposed.append(new_evidence)

            if verbose:
                print(f"  Decomposed discontinuous segment {i+1}: {segment}")
                print(f"    Note: No reference length available for coverage calculation")

    if not decomposed and verbose:
        print(f"  No segments >= {min_domain} residues in discontinuous hit")

    return decomposed if decomposed else [evidence]
