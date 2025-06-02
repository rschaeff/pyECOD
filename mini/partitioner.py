"""Core partitioning algorithm with residue blocking"""

from typing import List, Set, TYPE_CHECKING

if TYPE_CHECKING:
    from mini.models import Evidence, Domain

def partition_domains(evidence_list: List['Evidence'], sequence_length: int, verbose: bool = False) -> List['Domain']:
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
    from mini.models import Domain  # Import here to avoid circular imports

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

    # Track rejection statistics
    rejection_stats = {
        'insufficient_new_coverage': 0,
        'too_much_overlap': 0,
        'too_small': 0,
        'insufficient_residues': 0
    }

    print(f"\nPartitioning with {len(evidence_list)} evidence items for {sequence_length} residue protein")
    print(f"Thresholds: NEW_COVERAGE>{NEW_COVERAGE_THRESHOLD:.0%}, OLD_COVERAGE<{OLD_COVERAGE_THRESHOLD:.0%}, MIN_SIZE={MIN_DOMAIN_SIZE}")

    for i, evidence in enumerate(sorted_evidence):
        # Stop if too few unused residues remain
        if len(unused_residues) < MIN_DOMAIN_SIZE:
            rejection_stats['insufficient_residues'] += len(sorted_evidence) - i
            break

        # Get positions covered by this evidence
        evidence_positions = set(evidence.query_range.to_positions_simple())

        # Skip tiny hits
        if len(evidence_positions) < MIN_DOMAIN_SIZE:
            rejection_stats['too_small'] += 1
            continue

        # Calculate coverage: what fraction of THIS HIT overlaps with used/unused
        positions_in_unused = evidence_positions.intersection(unused_residues)
        positions_in_used = evidence_positions.intersection(used_residues)

        new_coverage = len(positions_in_unused) / len(evidence_positions) if evidence_positions else 0
        used_coverage = len(positions_in_used) / len(evidence_positions) if evidence_positions else 0

        # Domain assignment decision
        if new_coverage > NEW_COVERAGE_THRESHOLD and used_coverage < OLD_COVERAGE_THRESHOLD:
            # Accept this as a domain
            domain = Domain(
                id=f"d{domain_num}",
                range=evidence.query_range,
                family=evidence.t_group or evidence.source_pdb or "unknown",
                evidence_count=1,
                source=evidence.type,
                evidence_items=[evidence]
            )

            # Mark residues as used (blocking)
            used_residues.update(evidence_positions)
            unused_residues.difference_update(evidence_positions)

            domains.append(domain)

            print(f"\n✓ DOMAIN {domain_num}: {domain.family} @ {domain.range}")
            print(f"  Source: {evidence.type}, Confidence: {evidence.confidence:.2f}")
            print(f"  Coverage: {new_coverage:.1%} new, {used_coverage:.1%} overlap")
            print(f"  Residues assigned: {len(evidence_positions)}, Remaining: {len(unused_residues)}")

            domain_num += 1
        else:
            # Track rejection reason
            if new_coverage <= NEW_COVERAGE_THRESHOLD:
                rejection_stats['insufficient_new_coverage'] += 1
                if verbose:
                    print(f"  REJECT {evidence.source_pdb} @ {evidence.query_range}: "
                          f"insufficient new coverage ({new_coverage:.1%})")
            else:
                rejection_stats['too_much_overlap'] += 1
                if verbose:
                    print(f"  REJECT {evidence.source_pdb} @ {evidence.query_range}: "
                          f"too much overlap ({used_coverage:.1%})")

    # Summary of rejections
    total_rejected = sum(rejection_stats.values())
    if total_rejected > 0:
        print(f"\nRejection summary ({total_rejected} evidence items rejected):")
        for reason, count in rejection_stats.items():
            if count > 0:
                print(f"  {reason.replace('_', ' ').title()}: {count}")

    # Coverage summary
    coverage = len(used_residues) / sequence_length if sequence_length > 0 else 0
    print(f"\nFinal coverage: {len(used_residues)}/{sequence_length} residues ({coverage:.1%})")

    # Sort domains by start position
    return sorted(domains, key=lambda d: d.range.segments[0].start)
