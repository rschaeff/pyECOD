"""Core partitioning algorithm with post-selection chain BLAST decomposition"""

from typing import List, Set, Dict, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    from mini.models import Evidence, Domain
    from mini.decomposer import DomainReference

def partition_domains(evidence_list: List['Evidence'],
                     sequence_length: int,
                     domain_definitions: Dict[Tuple[str, str], List['DomainReference']] = None,
                     verbose: bool = False) -> List['Domain']:
    """
    Partition domains with post-selection decomposition:

    CORRECT STRATEGY:
    1. Run residue blocking partitioning (normal algorithm)
    2. Select best evidence (including good chain BLAST hits)
    3. For selected chain BLAST domains, decompose them post-hoc
    4. Replace chain BLAST domains with decomposed components

    This preserves high-quality evidence while enabling decomposition.
    """
    from mini.models import Domain

    # STEP 1: STANDARD RESIDUE BLOCKING PARTITIONING
    print(f"\nSTEP 1: RESIDUE BLOCKING PARTITIONING")
    print("=" * 40)

    # Track used/unused residues
    used_residues = set()
    unused_residues = set(range(1, sequence_length + 1))

    # Evidence type precedence
    type_precedence = {
        'chain_blast': 0,
        'domain_blast': 1,
        'hhsearch': 2
    }

    # Thresholds
    NEW_COVERAGE_THRESHOLD = 0.7
    OLD_COVERAGE_THRESHOLD = 0.1
    MIN_DOMAIN_SIZE = 20

    # Sort evidence by quality (best first)
    sorted_evidence = sorted(evidence_list,
                           key=lambda e: (
                               type_precedence.get(e.type, 3),
                               -e.confidence,
                               e.evalue if e.evalue else 999
                           ))

    selected_domains = []
    domain_num = 1
    rejection_stats = {
        'insufficient_new_coverage': 0,
        'too_much_overlap': 0,
        'too_small': 0,
        'insufficient_residues': 0
    }

    print(f"Partitioning {len(evidence_list)} evidence items for {sequence_length} residue protein")
    print(f"Thresholds: NEW_COVERAGE>{NEW_COVERAGE_THRESHOLD:.0%}, OLD_COVERAGE<{OLD_COVERAGE_THRESHOLD:.0%}, MIN_SIZE={MIN_DOMAIN_SIZE}")

    with_ref_length = sum(1 for e in evidence_list if e.reference_length is not None)
    print(f"Evidence with reference lengths: {with_ref_length}/{len(evidence_list)}")

    if with_ref_length == 0:
        print("ERROR: No evidence has reference lengths")
        return []

    # Standard residue blocking algorithm
    for i, evidence in enumerate(sorted_evidence):
        if len(unused_residues) < MIN_DOMAIN_SIZE:
            rejection_stats['insufficient_residues'] += len(sorted_evidence) - i
            break

        evidence_positions = set(evidence.query_range.to_positions_simple())

        if len(evidence_positions) < MIN_DOMAIN_SIZE:
            rejection_stats['too_small'] += 1
            continue

        positions_in_unused = evidence_positions.intersection(unused_residues)
        positions_in_used = evidence_positions.intersection(used_residues)

        new_coverage = len(positions_in_unused) / len(evidence_positions) if evidence_positions else 0
        used_coverage = len(positions_in_used) / len(evidence_positions) if evidence_positions else 0

        if new_coverage > NEW_COVERAGE_THRESHOLD and used_coverage < OLD_COVERAGE_THRESHOLD:
            # Accept this evidence as a domain
            family = evidence.t_group or evidence.source_pdb or evidence.domain_id or "unknown"

            domain = Domain(
                id=f"d{domain_num}",
                range=evidence.query_range,
                family=family,
                evidence_count=1,
                source=evidence.type,
                evidence_items=[evidence]
            )

            # Mark residues as used
            used_residues.update(evidence_positions)
            unused_residues.difference_update(evidence_positions)
            selected_domains.append(domain)

            print(f"✓ SELECTED: {domain.family} @ {domain.range} (source: {evidence.type})")
            print(f"  Coverage: {new_coverage:.1%} new, {used_coverage:.1%} overlap")
            print(f"  Residues assigned: {len(evidence_positions)}, Remaining: {len(unused_residues)}")

            domain_num += 1
        else:
            if new_coverage <= NEW_COVERAGE_THRESHOLD:
                rejection_stats['insufficient_new_coverage'] += 1
            else:
                rejection_stats['too_much_overlap'] += 1

    print(f"\nSelected {len(selected_domains)} domains before decomposition")

    # STEP 2: POST-SELECTION DECOMPOSITION
    print(f"\nSTEP 2: POST-SELECTION DECOMPOSITION")
    print("=" * 40)

    final_domains = []

    if domain_definitions:
        from mini.decomposer import decompose_chain_blast_with_mapping

        decomposition_stats = {
            'chain_blast_domains': 0,
            'decomposed': 0,
            'kept_original': 0
        }

        for domain in selected_domains:
            if domain.source == 'chain_blast':
                decomposition_stats['chain_blast_domains'] += 1

                # Get the original evidence
                evidence = domain.evidence_items[0]

                # ADD DEBUGGING HERE - let's see what we're actually working with
                print(f"\nDEBUG: Decomposing {domain.family}")
                print(f"  evidence.source_pdb = '{evidence.source_pdb}'")
                print(f"  evidence.domain_id = '{evidence.domain_id}'")

                # Parse the hit key
                if '_' in evidence.domain_id:
                    chain_part = evidence.domain_id.split('_')[-1]
                    print(f"  Parsed chain from domain_id: '{chain_part}'")
                else:
                    chain_part = 'A'
                    print(f"  No chain in domain_id, defaulting to: '{chain_part}'")

                hit_key = (evidence.source_pdb, chain_part)
                print(f"  Constructed hit_key: {hit_key}")

                # Check what's available in domain_definitions
                print(f"  Available domain_definitions keys: {list(domain_definitions.keys())[:10]}...")

                # Check specific matches
                exact_match = hit_key in domain_definitions
                print(f"  Exact match found: {exact_match}")

                if exact_match:
                    print(f"  Number of reference domains: {len(domain_definitions[hit_key])}")
                    domain_ids = [d.domain_id for d in domain_definitions[hit_key]]
                    print(f"  Reference domain IDs: {domain_ids}")

                # Check if decomposition is possible and worthwhile
                can_decompose = (
                    evidence.alignment is not None and
                    hit_key in domain_definitions and
                    len(domain_definitions[hit_key]) > 1  # Multi-domain protein
                )

                print(f"Evaluating {domain.family} for decomposition:")
                print(f"  Has alignment: {evidence.alignment is not None}")
                print(f"  Domain definitions available: {hit_key in domain_definitions}")
                if hit_key in domain_definitions:
                    print(f"  Reference domains: {len(domain_definitions[hit_key])}")

                if can_decompose:
                    print(f"  → Attempting decomposition...")

                    try:
                        decomposed_evidence = decompose_chain_blast_with_mapping(
                            evidence,
                            evidence.alignment.query_seq,
                            evidence.alignment.hit_seq,
                            evidence.alignment.query_start,
                            evidence.alignment.hit_start,
                            domain_definitions[hit_key],
                            verbose=verbose
                        )

                        if len(decomposed_evidence) > 1:
                            # Replace with decomposed domains
                            print(f"  ✓ Decomposed into {len(decomposed_evidence)} domains:")

                            for i, decomp_ev in enumerate(decomposed_evidence):
                                decomp_domain = Domain(
                                    id=f"d{len(final_domains) + i + 1}",
                                    range=decomp_ev.query_range,
                                    family=decomp_ev.domain_id,  # Use specific domain ID
                                    evidence_count=1,
                                    source="chain_blast_decomposed",
                                    evidence_items=[decomp_ev]
                                )
                                final_domains.append(decomp_domain)
                                print(f"    {decomp_ev.domain_id}: {decomp_ev.query_range}")

                            decomposition_stats['decomposed'] += 1
                        else:
                            # Keep original if decomposition didn't help
                            final_domains.append(domain)
                            decomposition_stats['kept_original'] += 1
                            print(f"  → Kept original (decomposition produced only 1 domain)")

                    except Exception as e:
                        print(f"  ✗ Decomposition failed: {e}")
                        final_domains.append(domain)
                        decomposition_stats['kept_original'] += 1
                else:
                    # Keep original if can't decompose
                    final_domains.append(domain)
                    decomposition_stats['kept_original'] += 1
                    print(f"  → Kept original (decomposition not possible)")
            else:
                # Non-chain BLAST domains pass through unchanged
                final_domains.append(domain)
                print(f"Keeping non-chain BLAST domain: {domain.family}")

        print(f"\nDecomposition summary:")
        print(f"  Chain BLAST domains: {decomposition_stats['chain_blast_domains']}")
        print(f"  Successfully decomposed: {decomposition_stats['decomposed']}")
        print(f"  Kept as original: {decomposition_stats['kept_original']}")
    else:
        final_domains = selected_domains
        print("No domain definitions - keeping selected domains as-is")

    # STEP 3: FINAL RESULTS
    print(f"\nSTEP 3: FINAL RESULTS")
    print("=" * 40)

    total_rejected = sum(rejection_stats.values())
    if total_rejected > 0:
        print(f"Rejected {total_rejected} evidence items:")
        for reason, count in rejection_stats.items():
            if count > 0:
                print(f"  {reason.replace('_', ' ').title()}: {count}")

    coverage = len(used_residues) / sequence_length if sequence_length > 0 else 0
    print(f"Final coverage: {len(used_residues)}/{sequence_length} residues ({coverage:.1%})")
    print(f"Final domains: {len(final_domains)}")

    for i, domain in enumerate(final_domains, 1):
        print(f"  {i}. {domain.family}: {domain.range} (source: {domain.source})")

    return sorted(final_domains, key=lambda d: d.range.segments[0].start)
