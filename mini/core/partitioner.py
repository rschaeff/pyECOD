# mini/core/partitioner.py (refactored and simplified with SequenceRange)
"""
Enhanced partitioning algorithm with iterative evidence processing and boundary optimization

Flow: Chain BLAST → Domain BLAST → HHsearch → Gap Analysis → Fragment Merging → Final domains

Simplified implementation using SequenceRange for all position operations.
"""

from typing import List, Set, Dict, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    from .models import Evidence, Domain
    from .decomposer import DomainReference

# Import for runtime use
from .models import Domain, DomainLayout
from .boundary_optimizer import BoundaryOptimizer
from .sequence_range import SequenceRange


def partition_domains(evidence_list: List['Evidence'],
                     sequence_length: int,
                     domain_definitions: Dict[Tuple[str, str], List['DomainReference']] = None,
                     min_domain_size: int = 25,
                     new_coverage_threshold: float = 0.7,
                     old_coverage_threshold: float = 0.1,
                     neighbor_tolerance: int = 5,
                     verbose: bool = False) -> List['Domain']:
    """
    Enhanced partitioning with iterative evidence processing and boundary optimization.

    Args:
        evidence_list: List of all evidence
        sequence_length: Total protein sequence length
        domain_definitions: Domain definitions for chain BLAST decomposition
        min_domain_size: Minimum domain size and fragment threshold
        new_coverage_threshold: Minimum new coverage required for selection
        old_coverage_threshold: Maximum overlap allowed with existing domains
        neighbor_tolerance: Distance tolerance for boundary optimization
        verbose: Whether to print detailed information

    Returns:
        List of optimized domains
    """

    print(f"\nENHANCED DOMAIN PARTITIONING")
    print("=" * 50)
    print(f"Processing {len(evidence_list)} evidence items for {sequence_length} residue protein")
    print(f"Thresholds: NEW_COVERAGE>{new_coverage_threshold:.0%}, OLD_COVERAGE<{old_coverage_threshold:.0%}, MIN_SIZE={min_domain_size}")

    # Initialize residue tracking using sets
    used_positions = set()
    unused_positions = set(range(1, sequence_length + 1))
    all_domains = []

    # Separate evidence by type
    evidence_by_type = _separate_evidence_by_type(evidence_list, domain_definitions, verbose)

    # PHASE 1: Chain BLAST with mandatory decomposition
    print(f"\nPHASE 1: CHAIN BLAST PROCESSING")
    print("=" * 30)

    chain_domains, used_positions, unused_positions = _process_chain_blast_evidence(
        evidence_by_type['chain_blast'], domain_definitions, used_positions, unused_positions,
        min_domain_size, new_coverage_threshold, old_coverage_threshold, verbose)

    all_domains.extend(chain_domains)
    _print_phase_summary("Chain BLAST", chain_domains, used_positions, sequence_length)

    # PHASE 2: Domain BLAST on remaining residues
    print(f"\nPHASE 2: DOMAIN BLAST PROCESSING")
    print("=" * 30)

    domain_blast_domains, used_positions, unused_positions = _process_standard_evidence(
        evidence_by_type['domain_blast'], used_positions, unused_positions,
        min_domain_size, new_coverage_threshold, old_coverage_threshold,
        domain_definitions, verbose)

    all_domains.extend(domain_blast_domains)
    _print_phase_summary("Domain BLAST", domain_blast_domains, used_positions, sequence_length)

    # PHASE 3: HHsearch on remaining residues
    print(f"\nPHASE 3: HHSEARCH PROCESSING")
    print("=" * 30)

    hhsearch_domains, used_positions, unused_positions = _process_standard_evidence(
        evidence_by_type['hhsearch'], used_positions, unused_positions,
        min_domain_size, new_coverage_threshold, old_coverage_threshold,
        domain_definitions, verbose)

    all_domains.extend(hhsearch_domains)
    _print_phase_summary("HHsearch", hhsearch_domains, used_positions, sequence_length)

    # PHASE 4: Boundary optimization
    print(f"\nPHASE 4: BOUNDARY OPTIMIZATION")
    print("=" * 30)

    if all_domains:
        layout = DomainLayout.from_domains(all_domains, sequence_length)
        optimizer = BoundaryOptimizer()
        optimized_layout = optimizer.optimize_boundaries(
            layout, min_domain_size, neighbor_tolerance, verbose)
        final_domains = optimized_layout.domains

        # Print final results
        final_stats = optimized_layout.get_coverage_stats()
        print(f"\nFINAL RESULTS:")
        print(f"  Coverage: {final_stats['assigned_residues']}/{final_stats['total_residues']} "
              f"residues ({final_stats['coverage_percent']:.1f}%)")
        print(f"  Domains: {final_stats['num_domains']}")
        print(f"  Remaining gaps: {final_stats['num_gaps']}")

        # Domain details
        for i, domain in enumerate(sorted(final_domains, key=lambda d: d.start_position), 1):
            print(f"    {i}. {domain.family}: {domain.range} (source: {domain.source})")

    else:
        print("No domains selected - optimization skipped")
        final_domains = []

    return final_domains

def get_domain_family_name(evidence, classification):
    """Determine domain family name with proper fallback logic"""

    # Prefer T-group if available
    if classification['t_group']:
        return classification['t_group']

    # Fallback to source_pdb for better test compatibility
    if evidence.source_pdb:
        return evidence.source_pdb

    # Last resort: domain_id
    if evidence.domain_id:
        return evidence.domain_id

    # Final fallback
    return 'unclassified'

def safe_extract_chain_id(evidence):
    """Safely extract chain ID from evidence, handling None domain_id"""
    if evidence.domain_id and isinstance(evidence.domain_id, str) and '_' in evidence.domain_id:
        return evidence.domain_id.split('_')[-1]
    return 'A'  # Default chain

def _separate_evidence_by_type(evidence_list: List['Evidence'],
                              domain_definitions: Dict = None,
                              verbose: bool = False) -> Dict[str, List['Evidence']]:
    """Separate evidence by type and apply blacklist filtering"""

    evidence_by_type = {
        'chain_blast': [],
        'domain_blast': [],
        'hhsearch': []
    }

    # Create blacklist for chain BLAST evidence
    blacklisted_chain_keys = set()
    if domain_definitions is not None:
        # Find chain BLAST targets that don't have domain definitions (blacklisted)
        for evidence in evidence_list:
            if evidence.type == 'chain_blast':
                chain = safe_extract_chain_id(evidence)
                target_key = (evidence.source_pdb, chain)
                if target_key not in domain_definitions:
                    blacklisted_chain_keys.add(target_key)

    if blacklisted_chain_keys and verbose:
        print(f"Blacklisted chain BLAST targets: {blacklisted_chain_keys}")

    # Separate and filter evidence
    blacklist_filtered_count = 0

    for evidence in evidence_list:
        if evidence.type == 'chain_blast':
            # Check blacklist
            chain = safe_extract_chain_id(evidence)
            target_key = (evidence.source_pdb, chain)

            if target_key in blacklisted_chain_keys:
                blacklist_filtered_count += 1
                if verbose:
                    print(f"  Filtered blacklisted chain BLAST: {evidence.source_pdb}_{chain}")
                continue

            evidence_by_type['chain_blast'].append(evidence)
        elif evidence.type == 'domain_blast':
            evidence_by_type['domain_blast'].append(evidence)
        elif evidence.type == 'hhsearch':
            evidence_by_type['hhsearch'].append(evidence)

    if blacklist_filtered_count > 0:
        print(f"Filtered {blacklist_filtered_count} blacklisted chain BLAST evidence")

    # Sort each type by priority
    for etype in evidence_by_type:
        evidence_by_type[etype] = _sort_evidence_by_priority(evidence_by_type[etype])

    return evidence_by_type


def _sort_evidence_by_priority(evidence_list: List['Evidence']) -> List['Evidence']:
    """Sort evidence by confidence, e-value, and other factors"""

    def evidence_sort_key(e):
        # Calculate tie-breaking score
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
            -e.confidence,                    # Higher confidence first
            e.evalue if e.evalue else 999,   # Lower e-value first
            -tiebreak_score,                 # Higher tiebreak score first
            e.source_pdb or "",              # Alphabetical by PDB
            e.domain_id or "",               # Alphabetical by domain ID
            str(e.query_range)               # Range as final tie-breaker
        )

    return sorted(evidence_list, key=evidence_sort_key)


def _process_chain_blast_evidence(chain_evidence: List['Evidence'],
                                domain_definitions: Dict,
                                used_positions: Set[int],
                                unused_positions: Set[int],
                                min_domain_size: int,
                                new_coverage_threshold: float,
                                old_coverage_threshold: float,
                                verbose: bool = False) -> Tuple[List['Domain'], Set[int], Set[int]]:
    """Process chain BLAST evidence with mandatory decomposition"""

    selected_domains = []
    decomposition_stats = {
        'evaluated': 0,
        'decomposed': 0,
        'rejected_no_decomp': 0,
        'rejected_failed_decomp': 0,
        'rejected_coverage': 0
    }

    for evidence in chain_evidence:
        decomposition_stats['evaluated'] += 1

        # Check if residues are still available
        evidence_positions = evidence.get_positions()
        if len(evidence_positions) < min_domain_size:
            continue

        # Check coverage thresholds
        positions_in_unused = evidence_positions.intersection(unused_positions)
        positions_in_used = evidence_positions.intersection(used_positions)

        new_coverage = len(positions_in_unused) / len(evidence_positions) if evidence_positions else 0
        used_coverage = len(positions_in_used) / len(evidence_positions) if evidence_positions else 0

        if new_coverage <= new_coverage_threshold or used_coverage >= old_coverage_threshold:
            decomposition_stats['rejected_coverage'] += 1
            continue

        # Attempt decomposition
        if not _can_attempt_decomposition(evidence, domain_definitions):
            decomposition_stats['rejected_no_decomp'] += 1
            if verbose:
                print(f"  Rejected {evidence.source_pdb}: no decomposition data")
            continue

        # Try decomposition
        decomposed_evidence = _attempt_decomposition(evidence, domain_definitions, verbose)

        if not decomposed_evidence or len(decomposed_evidence) < 1:
            decomposition_stats['rejected_failed_decomp'] += 1
            if verbose:
                print(f"  Rejected {evidence.source_pdb}: decomposition failed")
            continue

        # Decomposition succeeded - create domains
        decomposition_stats['decomposed'] += 1

        for i, dec_evidence in enumerate(decomposed_evidence):
            classification = get_evidence_classification(dec_evidence, domain_definitions)

            family_name = get_domain_family_name(dec_evidence, classification)

            domain = Domain(
                id=f"d{len(selected_domains) + 1}",
                range=dec_evidence.query_range,
                family=family_name,
                evidence_count=1,
                source=dec_evidence.type,
                evidence_items=[dec_evidence]
            )

            domain.x_group = classification['x_group']
            domain.h_group = classification['h_group']
            domain.t_group = classification['t_group']

            # Block residues
            domain_positions = domain.get_positions()
            used_positions.update(domain_positions)
            unused_positions.difference_update(domain_positions)

            selected_domains.append(domain)

            if verbose:
                print(f"  ✓ Decomposed domain {domain.id}: {domain.family} @ {domain.range}")

    # Print decomposition summary
    print(f"Chain BLAST decomposition summary:")
    for stat, count in decomposition_stats.items():
        if count > 0:
            print(f"  {stat.replace('_', ' ').title()}: {count}")

    return selected_domains, used_positions, unused_positions


def _process_standard_evidence(evidence_list: List['Evidence'],
                             used_positions: Set[int],
                             unused_positions: Set[int],
                             min_domain_size: int,
                             new_coverage_threshold: float,
                             old_coverage_threshold: float,
                             domain_definitions: Dict = None,
                             verbose: bool = False) -> Tuple[List['Domain'], Set[int], Set[int]]:
    """Process domain BLAST or HHsearch evidence"""

    selected_domains = []

    for evidence in evidence_list:
        # Check remaining residues
        if len(unused_positions) < min_domain_size:
            break

        evidence_positions = evidence.get_positions()

        # Skip tiny evidence
        if len(evidence_positions) < min_domain_size:
            continue

        # Check coverage thresholds
        positions_in_unused = evidence_positions.intersection(unused_positions)
        positions_in_used = evidence_positions.intersection(used_positions)

        new_coverage = len(positions_in_unused) / len(evidence_positions) if evidence_positions else 0
        used_coverage = len(positions_in_used) / len(evidence_positions) if evidence_positions else 0

        if new_coverage > new_coverage_threshold and used_coverage < old_coverage_threshold:
            # Accept this evidence
            classification = get_evidence_classification(evidence, domain_definitions)

            family_name = get_domain_family_name(evidence, classification)

            domain = Domain(
                id=f"d{len(selected_domains) + 1}",
                range=evidence.query_range,
                family=family_name,
                evidence_count=1,
                source=evidence.type,
                evidence_items=[evidence]
            )

            domain.x_group = classification['x_group']
            domain.h_group = classification['h_group']
            domain.t_group = classification['t_group']

            # Block residues
            used_positions.update(evidence_positions)
            unused_positions.difference_update(evidence_positions)
            selected_domains.append(domain)

            if verbose:
                print(f"  ✓ Selected {domain.id}: {domain.family} @ {domain.range}")
                print(f"    Coverage: {new_coverage:.1%} new, {used_coverage:.1%} overlap")

    return selected_domains, used_positions, unused_positions


def _can_attempt_decomposition(evidence: 'Evidence', domain_definitions: Dict) -> bool:
    """Check if decomposition can be attempted for chain BLAST evidence"""

    if not evidence.alignment or not domain_definitions:
        return False

    # Parse hit key
    chain = evidence.domain_id.split('_')[-1] if '_' in evidence.domain_id else 'A'
    hit_key = (evidence.source_pdb, chain)

    return hit_key in domain_definitions


def _attempt_decomposition(evidence: 'Evidence', domain_definitions: Dict,
                         verbose: bool = False) -> List['Evidence']:
    """Attempt to decompose chain BLAST evidence"""

    from mini.core.decomposer import decompose_chain_blast_with_mapping

    # Parse hit key
    chain = evidence.domain_id.split('_')[-1] if '_' in evidence.domain_id else 'A'
    hit_key = (evidence.source_pdb, chain)

    if hit_key not in domain_definitions:
        return []

    domain_refs = domain_definitions[hit_key]
    alignment = evidence.alignment

    try:
        decomposed = decompose_chain_blast_with_mapping(
            evidence=evidence,
            hit_query_str=alignment.query_seq,
            hit_hit_str=alignment.hit_seq,
            query_start=alignment.query_start,
            hit_start=alignment.hit_start,
            domain_refs=domain_refs,
            verbose=verbose
        )

        # Filter for successful decomposition
        if len(decomposed) > 1:
            return decomposed  # Multi-domain decomposition
        elif len(decomposed) == 1 and len(domain_refs) == 1:
            # Single domain reference - check if decomposition actually worked
            if decomposed[0].type == "chain_blast_decomposed":
                return decomposed

        return []  # Decomposition failed or insufficient

    except Exception as e:
        if verbose:
            print(f"    Decomposition error: {e}")
        return []


def _print_phase_summary(phase_name: str, domains: List['Domain'],
                        used_positions: Set[int], sequence_length: int) -> None:
    """Print summary for a processing phase"""

    coverage = len(used_positions) / sequence_length * 100 if sequence_length > 0 else 0
    print(f"{phase_name} results: {len(domains)} domains, "
          f"{len(used_positions)}/{sequence_length} residues ({coverage:.1f}% coverage)")


# Classification functions (keeping existing implementation)
def parse_ecod_hierarchy(t_group_str: str) -> tuple:
    """Parse ECOD T-group into hierarchical components"""
    if not t_group_str:
        return None, None, None

    parts = t_group_str.split('.')
    if len(parts) >= 3:
        x_group = parts[0]
        h_group = f"{parts[0]}.{parts[1]}"
        t_group = f"{parts[0]}.{parts[1]}.{parts[2]}"
        return x_group, h_group, t_group

    return None, None, None


def get_evidence_classification(evidence, domain_definitions=None):
    """Get ECOD taxonomic classification for evidence with better fallbacks"""

    # First try: Direct T-group from evidence (if available)
    if evidence.t_group:
        x_group, h_group, t_group = parse_ecod_hierarchy(evidence.t_group)
        return {
            'x_group': x_group,
            'h_group': h_group,
            't_group': t_group
        }

    # Second try: Lookup domain_id in domain_definitions
    if evidence.domain_id and domain_definitions:
        # Look for exact domain match
        for domain_refs in domain_definitions.values():
            for ref in domain_refs:
                if ref.domain_id == evidence.domain_id and ref.t_group:
                    x_group, h_group, t_group = parse_ecod_hierarchy(ref.t_group)
                    return {
                        'x_group': x_group,
                        'h_group': h_group,
                        't_group': t_group
                    }

    # Fallback: unclassified (but we'll use this in family assignment logic)
    return {
        'x_group': None,
        'h_group': None,
        't_group': None
    }
