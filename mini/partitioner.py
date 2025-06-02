# mini/partitioner.py
"""Core partitioning algorithm with residue blocking and proper chain BLAST decomposition"""

from typing import List, Set, Dict, Tuple, Optional
from .models import Evidence, Domain
from .decomposer import decompose_chain_blast_with_mapping, decompose_chain_blast_discontinuous, DomainReference
from ecod.core.sequence_range import SequenceRange

def partition_domains(evidence_list: List[Evidence],
                     sequence_length: int,
                     domain_defs: Optional[Dict[Tuple[str, str], List[DomainReference]]] = None,
                     use_precedence: bool = True) -> List[Domain]:
    """
    Partition domains with residue blocking and chain BLAST decomposition.

    Key principles:
    1. Domains are exclusive - each residue belongs to ONE domain
    2. Process evidence in priority order (chain_blast > domain_blast > hhsearch)
    3. Use coverage thresholds to decide domain assignment:
       - NEW_COVERAGE > 0.7: Hit must cover >70% unused residues
       - OLD_COVERAGE < 0.1: Hit must overlap <10% with existing domains
    4. Once residues are assigned, they're blocked from reassignment

    Args:
        evidence_list: List of evidence
        sequence_length: Total sequence length
        domain_defs: Domain definitions for decomposition
        use_precedence: Use evidence type precedence
    """
    # Track used/unused residues
    used_residues = set()
    unused_residues = set(range(1, sequence_length + 1))

    # Thresholds from Perl implementation analysis
    NEW_COVERAGE_THRESHOLD = 0.7   # Hit must be >70% new residues
    OLD_COVERAGE_THRESHOLD = 0.1   # Hit must have <10% overlap with existing
    MIN_DOMAIN_SIZE = 20          # Minimum domain size

    # Process chain BLAST hits with decomposition
    processed_evidence = []
    domain_defs = domain_defs or {}

    for evidence in evidence_list:
        if evidence.type == "chain_blast" and evidence.alignment is not None:
            # Look up domain definitions for this hit
            hit_key = (evidence.source_pdb, evidence.domain_id.split('_')[-1])  # Extract chain

            if hit_key in domain_defs:
                print(f"\nDecomposing chain BLAST hit {evidence.source_pdb}_{hit_key[1]}...")
                decomposed = decompose_chain_blast_with_mapping(
                    evidence,
                    evidence.alignment.query_seq,
                    evidence.alignment.hit_seq,
                    evidence.alignment.query_start,
                    evidence.alignment.hit_start,
                    domain_defs[hit_key]
                )
                processed_evidence.extend(decomposed)
            else:
                print(f"No domain definitions found for {hit_key}, using as-is")
                processed_evidence.append(evidence)
        else:
            processed_evidence.append(evidence)

    # Sort evidence by precedence
    if use_precedence:
        # CORRECT precedence: chain_blast is HIGHEST priority
        precedence_map = {
            'chain_blast': 1.5,           # Slightly lower to prefer decomposed
            'chain_blast_decomposed': 1,  # Highest - these are more specific
            'domain_blast': 2,            # Good but less comprehensive
            'hhsearch': 3                 # Lowest precedence
        }
        sorted_evidence = sorted(processed_evidence,
                               key=lambda e: (precedence_map.get(e.type, 99),
                                            e.evalue if e.evalue else 999))
    else:
        # Sort by confidence
        sorted_evidence = sorted(processed_evidence,
                               key=lambda e: (-e.confidence, e.evalue if e.evalue else 999))

    domains = []
    domain_num = 1

    print(f"\nProcessing {len(sorted_evidence)} evidence items...")
    print(f"Evidence type distribution:")
    type_counts = {}
    for e in sorted_evidence:
        type_counts[e.type] = type_counts.get(e.type, 0) + 1
    for t, c in sorted(type_counts.items()):
        print(f"  {t}: {c}")

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
            # Determine family - prefer classification, then PDB
            family = evidence.t_group or evidence.h_group or evidence.source_pdb or "unknown"

            # For decomposed hits, try to use the domain ID's PDB
            if evidence.type == "chain_blast_decomposed" and evidence.domain_id:
                # Extract PDB from domain ID like 'e2ia4A1'
                if len(evidence.domain_id) > 4 and evidence.domain_id[0] == 'e':
                    family = evidence.domain_id[1:5]

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
            #if new_coverage <= NEW_COVERAGE_THRESHOLD:
                #print(f"REJECT {evidence.source_pdb} @ {evidence.query_range}: "
                 #     f"insufficient new coverage ({new_coverage:.1%} <= {NEW_COVERAGE_THRESHOLD:.1%})")
            #else:
                #print(f"REJECT {evidence.source_pdb} @ {evidence.query_range}: "
                 #     f"too much overlap ({used_coverage:.1%} > {OLD_COVERAGE_THRESHOLD:.1%})")

    # Sort domains by start position
    return sorted(domains, key=lambda d: d.range.segments[0].start)
