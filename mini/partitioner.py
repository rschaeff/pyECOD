"""Core partitioning algorithm with residue blocking and chain BLAST decomposition"""

from typing import List, Set, Dict, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    from mini.models import Evidence, Domain
    from mini.decomposer import DomainReference

def partition_domains(evidence_list: List['Evidence'],
                     sequence_length: int,
                     domain_definitions: Dict[Tuple[str, str], List['DomainReference']] = None,
                     verbose: bool = False) -> List['Domain']:
    """
    Partition domains with residue blocking and chain BLAST decomposition:

    Key principles:
    1. DECOMPOSITION: Break chain BLAST hits into component domains using alignment data
    2. RESIDUE BLOCKING: Each residue belongs to ONE domain
    3. PRECEDENCE: Process evidence in priority order (best evidence first)
    4. THRESHOLDS: Use coverage thresholds to decide domain assignment:
       - NEW_COVERAGE > 0.7: Hit must cover >70% unused residues
       - OLD_COVERAGE < 0.1: Hit must overlap <10% with existing domains
    5. BLOCKING: Once residues are assigned, they're blocked from reassignment

    Args:
        evidence_list: List of evidence items
        sequence_length: Target protein sequence length
        domain_definitions: Domain architecture definitions for decomposition
        verbose: Enable verbose output
    """
    from mini.models import Domain  # Import here to avoid circular imports

    # STEP 1: DECOMPOSITION - Process chain BLAST hits first
    print(f"\nSTEP 1: CHAIN BLAST DECOMPOSITION")
    print("=" * 40)

    processed_evidence = []

    if domain_definitions:
        print(f"Decomposing chain BLAST hits using {len(domain_definitions)} protein definitions...")

        from mini.decomposer import decompose_chain_blast_with_mapping, decompose_chain_blast_discontinuous

        decomposition_stats = {
            'chain_blast_total': 0,
            'alignment_decomposed': 0,
            'discontinuous_decomposed': 0,
            'kept_original': 0,
            'total_before': len(evidence_list),
            'total_after': 0
        }

        for evidence in evidence_list:
            if evidence.type == 'chain_blast':
                decomposition_stats['chain_blast_total'] += 1

                # Try alignment-based decomposition first
                hit_key = (evidence.source_pdb, evidence.domain_id.split('_')[-1] if '_' in evidence.domain_id else 'A')

                if evidence.alignment is not None and hit_key in domain_definitions:
                    if verbose:
                        print(f"  Decomposing {evidence.source_pdb}_{hit_key[1]} using alignment data...")

                    try:
                        decomposed = decompose_chain_blast_with_mapping(
                            evidence,
                            evidence.alignment.query_seq,
                            evidence.alignment.hit_seq,
                            evidence.alignment.query_start,
                            evidence.alignment.hit_start,
                            domain_definitions[hit_key],
                            verbose=verbose
                        )

                        if len(decomposed) > 1:
                            processed_evidence.extend(decomposed)
                            decomposition_stats['alignment_decomposed'] += 1
                            if verbose:
                                print(f"    ✓ Decomposed into {len(decomposed)} domains")
                            continue
                    except Exception as e:
                        if verbose:
                            print(f"    ✗ Alignment decomposition failed: {e}")

                # Try discontinuous decomposition as fallback
                if evidence.query_range.is_discontinuous:
                    if verbose:
                        print(f"  Decomposing discontinuous {evidence.source_pdb} by segments...")

                    try:
                        decomposed = decompose_chain_blast_discontinuous(evidence, verbose=verbose)

                        if len(decomposed) > 1:
                            processed_evidence.extend(decomposed)
                            decomposition_stats['discontinuous_decomposed'] += 1
                            if verbose:
                                print(f"    ✓ Decomposed into {len(decomposed)} segments")
                            continue
                    except Exception as e:
                        if verbose:
                            print(f"    ✗ Discontinuous decomposition failed: {e}")

                # Keep original if no decomposition worked
                processed_evidence.append(evidence)
                decomposition_stats['kept_original'] += 1
                if verbose:
                    print(f"  Keeping {evidence.source_pdb} as single domain")
            else:
                # Non-chain BLAST evidence passes through unchanged
                processed_evidence.append(evidence)

        decomposition_stats['total_after'] = len(processed_evidence)

        print(f"Decomposition summary:")
        print(f"  Chain BLAST hits found: {decomposition_stats['chain_blast_total']}")
        print(f"  Alignment-based decomposed: {decomposition_stats['alignment_decomposed']}")
        print(f"  Discontinuous decomposed: {decomposition_stats['discontinuous_decomposed']}")
        print(f"  Kept as original: {decomposition_stats['kept_original']}")
        print(f"  Evidence before: {decomposition_stats['total_before']}")
        print(f"  Evidence after: {decomposition_stats['total_after']}")
    else:
        # No domain definitions provided, use evidence as-is
        processed_evidence = evidence_list
        print("⚠️  No domain definitions provided - skipping chain BLAST decomposition")

    # STEP 2: RESIDUE BLOCKING PARTITIONING
    print(f"\nSTEP 2: RESIDUE BLOCKING PARTITIONING")
    print("=" * 40)

    # Track used/unused residues
    used_residues = set()
    unused_residues = set(range(1, sequence_length + 1))

    # Evidence type precedence: decomposed chain_blast has same priority as chain_blast
    type_precedence = {
        'chain_blast': 0,
        'chain_blast_decomposed': 0,
        'domain_blast': 1,
        'hhsearch': 2
    }

    # Thresholds from Perl implementation analysis
    NEW_COVERAGE_THRESHOLD = 0.7   # Hit must be >70% new residues
    OLD_COVERAGE_THRESHOLD = 0.1   # Hit must have <10% overlap with existing
    MIN_DOMAIN_SIZE = 20          # Minimum domain size

    # Sort evidence by quality (best first)
    sorted_evidence = sorted(processed_evidence,
                           key=lambda e: (
                               type_precedence.get(e.type, 3),
                               -e.confidence,
                               e.evalue if e.evalue else 999
                           ))

    domains = []
    domain_num = 1

    # Track rejection statistics
    rejection_stats = {
        'insufficient_new_coverage': 0,
        'too_much_overlap': 0,
        'too_small': 0,
        'insufficient_residues': 0
    }

    print(f"Partitioning with {len(processed_evidence)} evidence items for {sequence_length} residue protein")
    print(f"Thresholds: NEW_COVERAGE>{NEW_COVERAGE_THRESHOLD:.0%}, OLD_COVERAGE<{OLD_COVERAGE_THRESHOLD:.0%}, MIN_SIZE={MIN_DOMAIN_SIZE}")

    with_ref_length = sum(1 for e in processed_evidence if e.reference_length is not None)
    print(f"Evidence with reference lengths: {with_ref_length}/{len(processed_evidence)}")

    if with_ref_length == 0:
        print("\nERROR: No evidence has reference lengths. Cannot proceed.")
        return []

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

            # Choose family name (prioritize domain_id for decomposed evidence)
            if evidence.type == 'chain_blast_decomposed' and evidence.domain_id:
                family = evidence.domain_id
            else:
                family = evidence.t_group or evidence.source_pdb or evidence.domain_id or "unknown"

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

            print(f"\n✓ DOMAIN {domain_num}: {domain.family} @ {domain.range}")
            print(f"  Source: {evidence.type}, Confidence: {evidence.confidence:.2f}")
            print(f"  Coverage: {new_coverage:.1%} new, {used_coverage:.1%} overlap")
            print(f"  Residues assigned: {len(evidence_positions)}, Remaining: {len(unused_residues)}")

            domain_num += 1
        else:
            # Track rejection reason
            if new_coverage <= NEW_COVERAGE_THRESHOLD:
                rejection_stats['insufficient_new_coverage'] += 1
            else:
                rejection_stats['too_much_overlap'] += 1

    # STEP 3: SUMMARY
    print(f"\nSTEP 3: PARTITIONING SUMMARY")
    print("=" * 40)

    # Summary of rejections
    total_rejected = sum(rejection_stats.values())
    if total_rejected > 0:
        print(f"Rejection summary ({total_rejected} evidence items rejected):")
        for reason, count in rejection_stats.items():
            if count > 0:
                print(f"  {reason.replace('_', ' ').title()}: {count}")

    # Coverage summary
    coverage = len(used_residues) / sequence_length if sequence_length > 0 else 0
    print(f"\nFinal coverage: {len(used_residues)}/{sequence_length} residues ({coverage:.1%})")
    print(f"Final domains: {len(domains)}")

    # Sort domains by start position
    return sorted(domains, key=lambda d: d.range.segments[0].start)
