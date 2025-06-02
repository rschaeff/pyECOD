# mini/partitioner.py - Corrected coverage thresholds
import xml.etree.ElementTree as ET
from typing import List
from .models import Domain

def partition_domains(evidence_list: List[Evidence], sequence_length: int) -> List[Domain]:
    """
    Partition with proper coverage thresholds:
    - High threshold for new coverage (must be mostly new residues)
    - Low threshold for old coverage (allow ~4% overlap)
    """
    # Track used/unused residues
    used_residues = set()
    unused_residues = set(range(1, sequence_length + 1))

    # Correct thresholds based on Perl implementation
    NEW_COVERAGE_THRESHOLD = 0.7   # Hit must be >70% new residues
    OLD_COVERAGE_THRESHOLD = 0.1   # Hit must be <10% overlap with existing

    # Sort evidence by quality
    sorted_evidence = sorted(evidence_list,
                           key=lambda e: (-e.confidence, e.evalue if e.evalue else 999))

    domains = []
    domain_num = 1

    for evidence in sorted_evidence:
        # Skip if too few unused residues
        if len(unused_residues) < 20:
            break

        # Get positions covered by this evidence
        evidence_positions = set(evidence.query_range.get_positions())

        # Calculate how much of THIS HIT overlaps with used/unused
        positions_in_unused = evidence_positions.intersection(unused_residues)
        positions_in_used = evidence_positions.intersection(used_residues)

        # Coverage = fraction of hit that overlaps with set
        new_coverage = len(positions_in_unused) / len(evidence_positions) if evidence_positions else 0
        used_coverage = len(positions_in_used) / len(evidence_positions) if evidence_positions else 0

        # Domain assignment decision
        if new_coverage > NEW_COVERAGE_THRESHOLD and used_coverage < OLD_COVERAGE_THRESHOLD:
            # Accept this as a domain
            domain = Domain(
                id=f"d{domain_num}",
                range=evidence.query_range,
                family=evidence.source_pdb or evidence.t_group or "unknown",
                evidence_count=1,
                source=evidence.type,
                evidence_items=[evidence]
            )

            # Mark residues as used
            used_residues.update(evidence_positions)
            unused_residues.difference_update(evidence_positions)

            domains.append(domain)
            domain_num += 1

            print(f"DEFINE domain {domain.id}: {domain.family} @ {domain.range} "
                  f"(new={new_coverage:.1%}, used={used_coverage:.1%})")
        else:
            # Debug why rejected
            if new_coverage <= NEW_COVERAGE_THRESHOLD:
                print(f"REJECT {evidence.source_pdb} @ {evidence.query_range}: "
                      f"insufficient new coverage ({new_coverage:.1%})")
            else:
                print(f"REJECT {evidence.source_pdb} @ {evidence.query_range}: "
                      f"too much overlap ({used_coverage:.1%})")

    return domains
