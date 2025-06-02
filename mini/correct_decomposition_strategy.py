#!/usr/bin/env python3
"""
Correct decomposition strategy: Post-selection decomposition

The RIGHT way:
1. Run residue blocking partitioning (as normal)
2. Select best evidence (including good chain BLAST hits) 
3. For selected chain BLAST domains, decompose them post-hoc
4. Replace chain BLAST domains with their decomposed components

This preserves the high-quality chain BLAST evidence while still getting decomposition.
"""

def create_post_selection_partitioner():
    """Create partitioner with post-selection decomposition"""
    
    partitioner_code = '''"""Core partitioning algorithm with post-selection chain BLAST decomposition"""

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
    
    1. Run residue blocking partitioning (normal algorithm)
    2. Select best evidence (including good chain BLAST hits)
    3. For selected chain BLAST domains, decompose them post-hoc
    4. Replace chain BLAST domains with decomposed components
    
    This preserves high-quality evidence while enabling decomposition.
    """
    from mini.models import Domain

    # STEP 1: STANDARD RESIDUE BLOCKING PARTITIONING
    print(f"\\nSTEP 1: RESIDUE BLOCKING PARTITIONING")
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

            domain_num += 1
        else:
            if new_coverage <= NEW_COVERAGE_THRESHOLD:
                rejection_stats['insufficient_new_coverage'] += 1
            else:
                rejection_stats['too_much_overlap'] += 1

    print(f"\\nSelected {len(selected_domains)} domains before decomposition")

    # STEP 2: POST-SELECTION DECOMPOSITION
    print(f"\\nSTEP 2: POST-SELECTION DECOMPOSITION")
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
                hit_key = (evidence.source_pdb, evidence.domain_id.split('_')[-1] if '_' in evidence.domain_id else 'A')
                
                # Check if decomposition is possible and worthwhile
                can_decompose = (
                    evidence.alignment is not None and
                    hit_key in domain_definitions and
                    len(domain_definitions[hit_key]) > 1  # Multi-domain protein
                )
                
                if can_decompose:
                    if verbose:
                        print(f"  Decomposing selected domain: {domain.family}")
                    
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
                            
                            decomposition_stats['decomposed'] += 1
                            print(f"    ✓ Decomposed into {len(decomposed_evidence)} domains")
                        else:
                            # Keep original if decomposition didn't help
                            final_domains.append(domain)
                            decomposition_stats['kept_original'] += 1
                            
                    except Exception as e:
                        if verbose:
                            print(f"    ✗ Decomposition failed: {e}")
                        final_domains.append(domain)
                        decomposition_stats['kept_original'] += 1
                else:
                    # Keep original if can't decompose
                    final_domains.append(domain)
                    decomposition_stats['kept_original'] += 1
                    if verbose:
                        print(f"  Keeping {domain.family} (no decomposition available)")
            else:
                # Non-chain BLAST domains pass through unchanged
                final_domains.append(domain)
        
        print(f"Decomposition summary:")
        print(f"  Chain BLAST domains: {decomposition_stats['chain_blast_domains']}")
        print(f"  Successfully decomposed: {decomposition_stats['decomposed']}")
        print(f"  Kept as original: {decomposition_stats['kept_original']}")
    else:
        final_domains = selected_domains
        print("No domain definitions - keeping selected domains as-is")

    # STEP 3: SUMMARY
    print(f"\\nSTEP 3: FINAL RESULTS")
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

    return sorted(final_domains, key=lambda d: d.range.segments[0].start)
'''

    # Write the corrected partitioner
    with open("mini/partitioner_correct.py", 'w') as f:
        f.write(partitioner_code)
    
    print("Created mini/partitioner_correct.py with post-selection decomposition")

def main():
    """Main function"""
    
    print("CORRECT DECOMPOSITION STRATEGY")
    print("=" * 50)
    
    print("WRONG approach (what I was doing):")
    print("1. Decompose chain BLAST hits into fragments")
    print("2. Compete fragments against whole hits")
    print("3. Fragments lose due to lower quality")
    print()
    
    print("RIGHT approach:")
    print("1. Run normal residue blocking partitioning")
    print("2. Select best evidence (including good chain BLAST hits)")  
    print("3. Post-hoc decompose selected chain BLAST domains")
    print("4. Replace with decomposed components")
    print()
    
    create_post_selection_partitioner()
    
    print("\\nThis preserves high-quality evidence while enabling decomposition!")

if __name__ == "__main__":
    main()
