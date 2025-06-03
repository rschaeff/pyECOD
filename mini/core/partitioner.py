"""Core partitioning algorithm with post-selection chain BLAST decomposition"""

from typing import List, Set, Dict, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    from .models import Evidence, Domain
    from .decomposer import DomainReference

# Import Domain class for runtime use
from mini.core.models import Domain

def partition_domains(evidence_list: List['Evidence'],
                     sequence_length: int,
                     domain_definitions: Dict[Tuple[str, str], List['DomainReference']] = None,
                     verbose: bool = False) -> List['Domain']:
    """
    Partition domains with blacklist-aware evidence filtering
    """
    # STEP 1: FILTER OUT BLACKLISTED CHAIN BLAST EVIDENCE
    print(f"\nSTEP 1: RESIDUE BLOCKING PARTITIONING")
    print("=" * 40)

    # Create blacklist lookup for chain BLAST evidence
    blacklisted_chain_keys = set()
    if domain_definitions is not None:
        # Find all PDB chains that should exist but don't (i.e., were blacklisted)
        all_chain_blast_pdbs = set()
        for ev in evidence_list:
            if ev.type == 'chain_blast':
                chain = ev.domain_id.split('_')[-1] if '_' in ev.domain_id else 'A'
                all_chain_blast_pdbs.add((ev.source_pdb, chain))

        # Identify blacklisted ones (PDB chains mentioned in evidence but not in domain_definitions)
        for (pdb, chain) in all_chain_blast_pdbs:
            if (pdb, chain) not in domain_definitions:
                blacklisted_chain_keys.add((pdb, chain))

        if blacklisted_chain_keys and verbose:
            print(f"Blacklisted chain BLAST targets: {blacklisted_chain_keys}")

    # Filter evidence to remove blacklisted chain BLAST hits
    filtered_evidence = []
    blacklist_filtered_count = 0

    for evidence in evidence_list:
        if evidence.type == 'chain_blast':
            # Check if this chain BLAST evidence targets a blacklisted reference
            chain = evidence.domain_id.split('_')[-1] if '_' in evidence.domain_id else 'A'
            target_key = (evidence.source_pdb, chain)

            if target_key in blacklisted_chain_keys:
                blacklist_filtered_count += 1
                if verbose:
                    print(f"  Filtered out blacklisted chain BLAST: {evidence.source_pdb}_{chain}")
                continue

        filtered_evidence.append(evidence)

    if blacklist_filtered_count > 0:
        print(f"Filtered out {blacklist_filtered_count} chain BLAST evidence targeting blacklisted references")

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

    # Create a tie-breaking score for evidence
    def evidence_sort_key(e):
        # Calculate a tie-breaking score based on multiple factors
        tiebreak_score = 0

        # Prefer evidence with alignment data
        if hasattr(e, 'alignment') and e.alignment is not None:
            tiebreak_score += 10

        # Prefer evidence with higher coverage
        if e.alignment_coverage is not None:
            tiebreak_score += e.alignment_coverage * 5

        # Prefer discontinuous ranges for complex domains
        if e.query_range.is_discontinuous:
            tiebreak_score += 2

        # Prefer evidence with reference length
        if e.reference_length is not None:
            tiebreak_score += 1

        return (
            type_precedence.get(e.type, 3),           # 1. Type precedence
            -e.confidence,                             # 2. Higher confidence first
            e.evalue if e.evalue else 999,            # 3. Lower e-value first
            -tiebreak_score,                           # 4. Higher tiebreak score first
            e.source_pdb or "",                       # 5. Alphabetical by PDB (deterministic)
            e.domain_id or "",                        # 6. Alphabetical by domain ID
            str(e.query_range)                        # 7. Range as final tie-breaker
        )

    # Sort evidence with the enhanced key
    sorted_evidence = sorted(filtered_evidence, key=evidence_sort_key)

    selected_domains = []
    domain_num = 1
    rejection_stats = {
        'insufficient_new_coverage': 0,
        'too_much_overlap': 0,
        'too_small': 0,
        'insufficient_residues': 0
    }

    print(f"Partitioning {len(sorted_evidence)} evidence items for {sequence_length} residue protein")
    print(f"Thresholds: NEW_COVERAGE>{NEW_COVERAGE_THRESHOLD:.0%}, OLD_COVERAGE<{OLD_COVERAGE_THRESHOLD:.0%}, MIN_SIZE={MIN_DOMAIN_SIZE}")

    with_ref_length = sum(1 for e in sorted_evidence if e.reference_length is not None)
    print(f"Evidence with reference lengths: {with_ref_length}/{len(sorted_evidence)}")

    if with_ref_length == 0:
        print("ERROR: No evidence has reference lengths")
        return []

    # Standard residue blocking algorithm (unchanged from here)
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
        from mini.core.decomposer import decompose_chain_blast_with_mapping

        decomposition_stats = {
            'chain_blast_domains': 0,
            'decomposed': 0,
            'kept_original': 0,
            'rejected': 0
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
                        # Extract alignment data
                        alignment = evidence.alignment
                        domain_refs = domain_definitions[hit_key]

                        # Call decomposition with proper arguments
                        decomposed_evidence = decompose_chain_blast_with_mapping(
                            evidence=evidence,
                            hit_query_str=alignment.query_seq,
                            hit_hit_str=alignment.hit_seq,
                            query_start=alignment.query_start,
                            hit_start=alignment.hit_start,
                            domain_refs=domain_refs,
                            verbose=verbose
                        )

                        if len(decomposed_evidence) > 1:
                            # Successful multi-domain decomposition
                            print(f"  ✓ Decomposed into {len(decomposed_evidence)} domains:")
                            for dec_ev in decomposed_evidence:
                                new_domain = Domain(
                                    id=f"d{len(final_domains) + 1}",
                                    range=dec_ev.query_range,
                                    family=dec_ev.domain_id or dec_ev.source_pdb or "decomposed",
                                    evidence_count=1,
                                    source="chain_blast_decomposed",
                                    evidence_items=[dec_ev]
                                )
                                final_domains.append(new_domain)
                                print(f"    → {new_domain.family}: {new_domain.range}")
                            decomposition_stats['decomposed'] += 1

                        elif len(decomposed_evidence) == 1 and len(domain_refs) == 1:
                            # Single domain reference, single domain result - this is OK
                            print(f"  → Single domain reference, keeping decomposed result")
                            dec_ev = decomposed_evidence[0]
                            new_domain = Domain(
                                id=f"d{len(final_domains) + 1}",
                                range=dec_ev.query_range,
                                family=dec_ev.domain_id or dec_ev.source_pdb or domain.family,
                                evidence_count=1,
                                source="chain_blast_decomposed",
                                evidence_items=[dec_ev]
                            )
                            final_domains.append(new_domain)
                            decomposition_stats['kept_original'] += 1

                        else:
                            # Multi-domain reference but decomposition failed/incomplete
                            print(f"  ✗ Decomposition insufficient for multi-domain reference")
                            print(f"    Reference has {len(domain_refs)} domains")
                            print(f"    Decomposition produced {len(decomposed_evidence)} domains")
                            print(f"    REJECTING non-decomposed chain BLAST hit")
                            decomposition_stats['rejected'] += 1
                            # DO NOT add to final_domains!

                    except Exception as e:
                        print(f"  ✗ Decomposition failed: {e}")
                        # For multi-domain references, reject on failure
                        if len(domain_definitions[hit_key]) > 1:
                            print(f"    REJECTING non-decomposed multi-domain hit")
                            decomposition_stats['rejected'] += 1
                        else:
                            # Single domain reference - keep original
                            final_domains.append(domain)
                            decomposition_stats['kept_original'] += 1
                else:
                    # No decomposition possible - only keep if it's likely single-domain
                    if domain.range.is_discontinuous:
                        print(f"  ✗ Discontinuous hit without decomposition - REJECTING")
                        decomposition_stats['rejected'] += 1
                    else:
                        final_domains.append(domain)
                        decomposition_stats['kept_original'] += 1
            else:
                # Non-chain BLAST domains pass through unchanged
                final_domains.append(domain)
                print(f"Keeping non-chain BLAST domain: {domain.family}")

        print(f"\nDecomposition summary:")
        print(f"  Chain BLAST domains: {decomposition_stats['chain_blast_domains']}")
        print(f"  Successfully decomposed: {decomposition_stats['decomposed']}")
        print(f"  Kept as original: {decomposition_stats['kept_original']}")
        print(f"  Rejected: {decomposition_stats['rejected']}")
    else:
        final_domains = selected_domains
        print("No domain definitions - keeping selected domains as-is")

    # STEP 2.5: REPROCESS UNBLOCKED RESIDUES
    if domain_definitions and 'decomposition_stats' in locals() and decomposition_stats.get('rejected', 0) > 0:
        print(f"\nSTEP 2.5: REPROCESSING UNBLOCKED RESIDUES")
        print("=" * 40)

        # Recalculate which residues are actually used by final domains
        actual_used_residues = set()
        for domain in final_domains:
            actual_used_residues.update(domain.range.to_positions_simple())

        # Find unblocked residues
        unblocked_residues = used_residues - actual_used_residues
        if unblocked_residues:
            print(f"Unblocked {len(unblocked_residues)} residues from rejected decompositions")

            # Reset residue tracking
            used_residues = actual_used_residues.copy()
            unused_residues = set(range(1, sequence_length + 1)) - used_residues

            # Process remaining evidence (domain_blast and hhsearch)
            remaining_evidence = [e for e in sorted_evidence
                                if e.type in ['domain_blast', 'hhsearch']]

            print(f"Processing {len(remaining_evidence)} remaining evidence items...")

            reselected_count = 0
            for evidence in remaining_evidence:
                # Skip if too few unused residues remain
                if len(unused_residues) < MIN_DOMAIN_SIZE:
                    break

                evidence_positions = set(evidence.query_range.to_positions_simple())

                # Skip tiny hits
                if len(evidence_positions) < MIN_DOMAIN_SIZE:
                    continue

                # Calculate coverage
                positions_in_unused = evidence_positions.intersection(unused_residues)
                positions_in_used = evidence_positions.intersection(used_residues)

                new_coverage = len(positions_in_unused) / len(evidence_positions) if evidence_positions else 0
                used_coverage = len(positions_in_used) / len(evidence_positions) if evidence_positions else 0

                # Apply same selection criteria
                if new_coverage > NEW_COVERAGE_THRESHOLD and used_coverage < OLD_COVERAGE_THRESHOLD:
                    # Accept this evidence
                    family = evidence.t_group or evidence.source_pdb or evidence.domain_id or "unknown"

                    domain = Domain(
                        id=f"d{len(final_domains) + 1}",
                        range=evidence.query_range,
                        family=family,
                        evidence_count=1,
                        source=evidence.type,
                        evidence_items=[evidence]
                    )

                    # Mark residues as used
                    used_residues.update(evidence_positions)
                    unused_residues.difference_update(evidence_positions)
                    final_domains.append(domain)

                    print(f"✓ RESELECTED: {domain.family} @ {domain.range} (source: {evidence.type})")
                    print(f"  Coverage: {new_coverage:.1%} new, {used_coverage:.1%} overlap")

                    reselected_count += 1

            if reselected_count > 0:
                print(f"Reselected {reselected_count} domains from unblocked regions")
            else:
                print("No additional domains selected from unblocked regions")

    # STEP 3: FINAL RESULTS
    print(f"\nSTEP 3: FINAL RESULTS")
    print("=" * 40)

    total_rejected = sum(rejection_stats.values())
    if total_rejected > 0:
        print(f"Rejected {total_rejected} evidence items:")
        for reason, count in rejection_stats.items():
            if count > 0:
                print(f"  {reason.replace('_', ' ').title()}: {count}")

    # Recalculate final coverage
    final_used_residues = set()
    for domain in final_domains:
        final_used_residues.update(domain.range.to_positions_simple())

    coverage = len(final_used_residues) / sequence_length if sequence_length > 0 else 0
    print(f"Final coverage: {len(final_used_residues)}/{sequence_length} residues ({coverage:.1%})")
    print(f"Final domains: {len(final_domains)}")

    for i, domain in enumerate(final_domains, 1):
        print(f"  {i}. {domain.family}: {domain.range} (source: {domain.source})")

    return sorted(final_domains, key=lambda d: d.range.segments[0].start)
