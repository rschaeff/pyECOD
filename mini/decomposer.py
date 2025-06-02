"""Core partitioning algorithm with residue blocking and chain BLAST decomposition"""

from typing import List, Set, Tuple
from .models import Evidence, Domain
from ecod.core.sequence_range import SequenceRange, SequenceSegment

def decompose_chain_blast_hit(evidence: Evidence, min_gap: int = 30, min_domain: int = 50) -> List[Evidence]:
    """
    Decompose a chain BLAST hit into potential domain boundaries.
    
    Strategy: Look for gaps in the alignment or natural break points.
    For now, simple heuristic: if hit is discontinuous, treat each segment as potential domain.
    If continuous and large, look for potential break points.
    
    Args:
        evidence: Chain BLAST evidence to decompose
        min_gap: Minimum gap size to consider a break
        min_domain: Minimum domain size to consider
        
    Returns:
        List of decomposed evidence pieces
    """
    if not evidence.type == "chain_blast":
        return [evidence]
    
    # If already discontinuous, each segment could be a domain
    if evidence.query_range.is_discontinuous:
        decomposed = []
        for i, segment in enumerate(evidence.query_range.segments):
            if segment.length >= min_domain:
                # Create new evidence for this segment
                new_evidence = Evidence(
                    type="chain_blast_decomposed",
                    source_pdb=evidence.source_pdb,
                    query_range=SequenceRange([segment]),
                    confidence=evidence.confidence * 0.9,  # Slight penalty for decomposition
                    evalue=evidence.evalue,
                    domain_id=f"{evidence.domain_id}_seg{i+1}"
                )
                decomposed.append(new_evidence)
        return decomposed if decomposed else [evidence]
    
    # For continuous hits, we need more sophisticated decomposition
    # For now, accept large continuous chain hits as-is
    # TODO: Implement gap detection in alignments
    return [evidence]

def partition_domains(evidence_list: List[Evidence], sequence_length: int, 
                     use_precedence: bool = True) -> List[Domain]:
    """
    Partition domains with residue blocking (inspired by Perl implementation).

    Key principles:
    1. Domains are exclusive - each residue belongs to ONE domain
    2. Process evidence in priority order (precedence or confidence)
    3. Use coverage thresholds to decide domain assignment:
       - NEW_COVERAGE > 0.7: Hit must cover >70% unused residues
       - OLD_COVERAGE < 0.1: Hit must overlap <10% with existing domains
    4. Once residues are assigned, they're blocked from reassignment
    
    Args:
        evidence_list: List of evidence
        sequence_length: Total sequence length
        use_precedence: Use evidence type precedence instead of confidence
    """
    # Track used/unused residues
    used_residues = set()
    unused_residues = set(range(1, sequence_length + 1))

    # Thresholds from Perl implementation analysis
    NEW_COVERAGE_THRESHOLD = 0.7   # Hit must be >70% new residues
    OLD_COVERAGE_THRESHOLD = 0.1   # Hit must have <10% overlap with existing
    MIN_DOMAIN_SIZE = 20          # Minimum domain size

    # Decompose chain BLAST hits
    decomposed_evidence = []
    for evidence in evidence_list:
        if evidence.type == "chain_blast":
            decomposed = decompose_chain_blast_hit(evidence)
            decomposed_evidence.extend(decomposed)
        else:
            decomposed_evidence.append(evidence)
    
    # Sort evidence by precedence or confidence
    if use_precedence:
        # Evidence precedence: domain_blast > hhsearch > chain_blast
        precedence_map = {
            'domain_blast': 1,
            'hhsearch': 2,
            'chain_blast': 3,
            'chain_blast_decomposed': 3.5
        }
        sorted_evidence = sorted(decomposed_evidence,
                               key=lambda e: (precedence_map.get(e.type, 99), 
                                            e.evalue if e.evalue else 999))
    else:
        # Sort by confidence
        sorted_evidence = sorted(decomposed_evidence,
                               key=lambda e: (-e.confidence, e.evalue if e.evalue else 999))

    domains = []
    domain_num = 1

    print(f"\nProcessing {len(sorted_evidence)} evidence items (after decomposition)...")

    for evidence in sorted_evidence:
        # Stop if too few unused residues remain
        if len(unused_residues) < MIN_DOMAIN_SIZE:
            break

        # Get positions covered by this evidence
        try:
            evidence_positions = set(evidence.query_range.to_positions_simple())
        except ValueError:
            # Multi-chain range, extract positions manually
            evidence_positions = set()
            for pos, chain in evidence.query_range.to_positions():
                evidence_positions.add(pos)

        # Skip tiny hits
        if len(evidence_positions) < MIN_DOMAIN_SIZE:
            continue

        # Calculate coverage: what fraction of THIS HIT overlaps with used/unused
        positions_in_unused = evidence_positions.intersection(unused_residues)
        positions_in_used = evidence_positions.intersection(used_residues)

        new_coverage = len(positions_in_unused) / len(evidence_positions) if evidence_positions else 0
        used_coverage = len(positions_in_used) / len(evidence_positions) if evidence_positions else 0

        # Domain assignment decision
        if new_coverage > NEW_COVERAGE_THRESHOLD and used_coverage < OLD_COVERAGE_THRESHOLD:
            # Determine family - use classification if available, otherwise PDB
            family = evidence.t_group or evidence.h_group or evidence.source_pdb or "unknown"

            # Accept this as a domain
            domain = Domain(
                id=f"d{domain_num}",
                range=evidence.query_range,
                family=family,
                evidence_count=1,
                source=evidence.type,
                evidence_items=[evidence]
            )

            # Mark residues as used (blocking)
            used_residues.update(evidence_positions)
            unused_residues.difference_update(evidence_positions)

            domains.append(domain)
            domain_num += 1

            print(f"DEFINE domain {domain.id}: {family} @ {domain.range} "
                  f"(new={new_coverage:.1%}, used={used_coverage:.1%}) [{evidence.type}]")
        else:
            # Debug output for rejected evidence
            if new_coverage <= NEW_COVERAGE_THRESHOLD:
                print(f"REJECT {evidence.source_pdb} @ {evidence.query_range}: "
                      f"insufficient new coverage ({new_coverage:.1%} <= {NEW_COVERAGE_THRESHOLD:.1%})")
            else:
                print(f"REJECT {evidence.source_pdb} @ {evidence.query_range}: "
                      f"too much overlap ({used_coverage:.1%} > {OLD_COVERAGE_THRESHOLD:.1%})")

    # Sort domains by start position - FIX: use .start attribute
    return sorted(domains, key=lambda d: d.range.segments[0].start)
