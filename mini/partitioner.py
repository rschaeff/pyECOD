# mini/partitioner.py
"""Core partitioning algorithm with residue blocking"""

from typing import List, Set
from .models import Evidence, Domain
from ecod.core.sequence_range import SequenceRange

def partition_domains(evidence_list: List[Evidence], sequence_length: int) -> List[Domain]:
    """
    Partition domains with residue blocking (inspired by Perl implementation):

    Key principles:
    1. Domains are exclusive - each residue belongs to ONE domain
    2. Process evidence in priority order (best evidence first)
    3. Use coverage thresholds to decide domain assignment:
       - NEW_COVERAGE > 0.7: Hit must cover >70% unused residues
       - OLD_COVERAGE < 0.1: Hit must overlap <10% with existing domains
    4. Once residues are assigned, they're blocked from reassignment
    """
    # Track used/unused residues
    used_residues = set()
    unused_residues = set(range(1, sequence_length + 1))

    # Thresholds from Perl implementation analysis
    NEW_COVERAGE_THRESHOLD = 0.7   # Hit must be >70% new residues
    OLD_COVERAGE_THRESHOLD = 0.1   # Hit must have <10% overlap with existing
    MIN_DOMAIN_SIZE = 20          # Minimum domain size

    # Sort evidence by quality (best first)
    sorted_evidence = sorted(evidence_list,
                           key=lambda e: (-e.confidence, e.evalue if e.evalue else 999))

    domains = []
    domain_num = 1

    for evidence in sorted_evidence:
        # Stop if too few unused residues remain
        if len(unused_residues) < MIN_DOMAIN_SIZE:
            break

        # Get positions covered by this evidence
        # FIXED: Use to_positions_simple() for single-chain ranges
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
                  f"(new={new_coverage:.1%}, used={used_coverage:.1%})")
        else:
            # Debug output for rejected evidence
            if new_coverage <= NEW_COVERAGE_THRESHOLD:
                print(f"REJECT {evidence.source_pdb} @ {evidence.query_range}: "
                      f"insufficient new coverage ({new_coverage:.1%} <= {NEW_COVERAGE_THRESHOLD:.1%})")
            else:
                print(f"REJECT {evidence.source_pdb} @ {evidence.query_range}: "
                      f"too much overlap ({used_coverage:.1%} > {OLD_COVERAGE_THRESHOLD:.1%})")

    # Sort domains by start position
    return sorted(domains, key=lambda d: d.range.segments[0][0])
