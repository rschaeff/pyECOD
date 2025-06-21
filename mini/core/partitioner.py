# mini/core/partitioner.py (enhanced with evidence quality thresholds)
"""
Enhanced partitioning algorithm with evidence quality thresholds and iterative processing

Flow: Chain BLAST → Domain BLAST → HHsearch → Gap Analysis → Fragment Merging → Final domains

ENHANCED: Adds confidence and reference coverage thresholds to prevent poor quality assignments
"""

from typing import List, Set, Dict, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    from .models import Evidence, Domain
    from .decomposer import DomainReference

# Import for runtime use
from .models import Domain, DomainLayout
from .boundary_optimizer import BoundaryOptimizer
from .sequence_range import SequenceRange

from .domain_utils import (
    create_domain_with_provenance,
    get_evidence_classification,
    get_domain_family_name,
    validate_domain_provenance
)
from .evidence_utils import (
    standardize_evidence_list,
    get_evidence_coverage_stats
)

# ENHANCED: Evidence quality thresholds
EVIDENCE_THRESHOLDS = {
    'hhsearch': {
        'min_confidence': 0.65,        # ~68% probability minimum
        'min_reference_coverage': 0.6,  # 60% of reference domain
        'description': 'HHsearch (remote homology)'
    },
    'domain_blast': {
        'min_confidence': 0.5,         # More lenient for domain BLAST
        'min_reference_coverage': 0.5,  # 50% of reference domain
        'description': 'Domain BLAST (reliable homology)'
    },
    'chain_blast_decomposed': {
        'min_confidence': 0.5,
        'min_reference_coverage': 0.5,  # 50% of reference domain
        'description': 'Chain BLAST decomposed'
    }
}

def partition_domains(evidence_list: List['Evidence'],
                     sequence_length: int,
                     domain_definitions: Dict[Tuple[str, str], List['DomainReference']] = None,
                     min_domain_size: int = 25,
                     new_coverage_threshold: float = 0.7,
                     old_coverage_threshold: float = 0.1,
                     neighbor_tolerance: int = 5,
                     apply_quality_thresholds: bool = True,
                     verbose: bool = False) -> List['Domain']:
    """
    Enhanced partitioning with evidence quality thresholds and iterative processing.

    Args:
        evidence_list: List of all evidence
        sequence_length: Total protein sequence length
        domain_definitions: Domain definitions for chain BLAST decomposition
        min_domain_size: Minimum domain size and fragment threshold
        new_coverage_threshold: Minimum new coverage required for selection
        old_coverage_threshold: Maximum overlap allowed with existing domains
        neighbor_tolerance: Distance tolerance for boundary optimization
        apply_quality_thresholds: Whether to apply evidence quality thresholds
        verbose: Whether to print detailed information

    Returns:
        List of optimized domains
    """

    print(f"\nENHANCED DOMAIN PARTITIONING")
    print("=" * 50)
    print(f"Processing {len(evidence_list)} evidence items for {sequence_length} residue protein")
    print(f"Coverage thresholds: NEW_COVERAGE>{new_coverage_threshold:.0%}, OLD_COVERAGE<{old_coverage_threshold:.0%}, MIN_SIZE={min_domain_size}")

    if apply_quality_thresholds:
        print(f"Quality thresholds: ENABLED")
        for etype, thresholds in EVIDENCE_THRESHOLDS.items():
            print(f"  {thresholds['description']}: confidence≥{thresholds['min_confidence']:.2f}, reference_coverage≥{thresholds['min_reference_coverage']:.0%}")
    else:
        print(f"Quality thresholds: DISABLED")

    # STANDARDIZED: Apply evidence standardization before processing
    evidence_list = standardize_evidence_list(
        evidence_list=evidence_list,
        reference_lengths=None,  # Already populated during parsing
        protein_lengths=None,    # Already populated during parsing
        domain_definitions=domain_definitions
    )

    # Show evidence quality summary
    if verbose:
        coverage_stats = get_evidence_coverage_stats(evidence_list)
        print(f"Evidence quality: {coverage_stats['with_complete_provenance']}/{coverage_stats['total']} "
              f"({coverage_stats['provenance_percentage']:.1f}%) with complete provenance")

    # Initialize residue tracking using sets
    used_positions = set()
    unused_positions = set(range(1, sequence_length + 1))
    all_domains = []

    # Separate evidence by type
    evidence_by_type = _separate_evidence_by_type(evidence_list, domain_definitions, verbose)

    # ENHANCED: Apply quality thresholds before processing
    if apply_quality_thresholds:
        evidence_by_type = _apply_quality_thresholds(evidence_by_type, verbose)

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

        # FIXED: Assign proper domain IDs
        for i, domain in enumerate(final_domains, 1):
            domain.id = f"d{i}"

        # Sort domains by sequence position for consistent output
        final_domains.sort(key=lambda d: d.start_position)

        # Print final results
        final_stats = optimized_layout.get_coverage_stats()
        print(f"\nFINAL RESULTS:")
        print(f"  Coverage: {final_stats['assigned_residues']}/{final_stats['total_residues']} "
              f"residues ({final_stats['coverage_percent']:.1f}%)")
        print(f"  Domains: {final_stats['num_domains']}")
        print(f"  Remaining gaps: {final_stats['num_gaps']}")

        # Domain details
        for i, domain in enumerate(final_domains, 1):
            print(f"    {i}. {domain.family}: {domain.range} (source: {domain.source})")

    else:
        print("No domains selected - optimization skipped")
        final_domains = []

    # STANDARDIZED: Validate all domains have proper provenance
    if verbose:
        validation_issues = 0
        for domain in final_domains:
            is_valid, issues = validate_domain_provenance(domain)
            if not is_valid:
                validation_issues += 1
                if verbose:
                    print(f"  Warning: Domain {domain.id} provenance issues: {issues}")

        if validation_issues == 0:
            print(f"  ✓ All {len(final_domains)} domains have complete provenance")
        else:
            print(f"  ⚠️  {validation_issues} domains have provenance issues")

    return final_domains


def _apply_quality_thresholds(evidence_by_type: Dict[str, List['Evidence']],
                             verbose: bool = False) -> Dict[str, List['Evidence']]:
    """Apply quality thresholds to filter poor evidence"""

    filtered_evidence = {}
    rejection_stats = {
        'confidence': 0,
        'reference_coverage': 0,
        'missing_reference_data': 0,
        'total_rejected': 0
    }

    for evidence_type, evidence_list in evidence_by_type.items():
        filtered_list = []

        # Get thresholds for this evidence type
        thresholds = EVIDENCE_THRESHOLDS.get(evidence_type, {})
        min_confidence = thresholds.get('min_confidence', 0.0)
        min_ref_coverage = thresholds.get('min_reference_coverage', 0.0)

        for evidence in evidence_list:
            rejected = False
            rejection_reason = None

            # Check confidence threshold (ALWAYS applied)
            if evidence.confidence < min_confidence:
                rejected = True
                rejection_reason = f"confidence {evidence.confidence:.3f} < {min_confidence:.3f}"
                rejection_stats['confidence'] += 1

            # For evidence types with reference coverage thresholds, REQUIRE complete data
            elif evidence_type in EVIDENCE_THRESHOLDS:
                if not evidence.reference_length or not evidence.hit_range:
                    rejected = True
                    rejection_reason = "incomplete reference data (missing hit_range or reference_length)"
                    rejection_stats['missing_reference_data'] += 1
                else:
                    # Calculate reference coverage and apply threshold
                    ref_coverage = evidence.hit_range.total_length / evidence.reference_length
                    if ref_coverage < min_ref_coverage:
                        rejected = True
                        rejection_reason = f"reference_coverage {ref_coverage:.1%} < {min_ref_coverage:.0%}"
                        rejection_stats['reference_coverage'] += 1

            if rejected:
                rejection_stats['total_rejected'] += 1
                if verbose:
                    print(f"  Rejected {evidence.domain_id or evidence.source_pdb}: {rejection_reason}")
            else:
                filtered_list.append(evidence)

        filtered_evidence[evidence_type] = filtered_list

    # Print rejection summary
    if rejection_stats['total_rejected'] > 0:
        print(f"\nQuality threshold filtering:")
        print(f"  Total rejected: {rejection_stats['total_rejected']}")
        if rejection_stats['confidence'] > 0:
            print(f"    Low confidence: {rejection_stats['confidence']}")
        if rejection_stats['reference_coverage'] > 0:
            print(f"    Poor reference coverage: {rejection_stats['reference_coverage']}")
        if rejection_stats['missing_reference_data'] > 0:
            print(f"    Incomplete reference data: {rejection_stats['missing_reference_data']}")

        # Show remaining evidence counts
        remaining_total = sum(len(elist) for elist in filtered_evidence.values())
        print(f"  Remaining evidence: {remaining_total}")
        for etype, elist in filtered_evidence.items():
            if elist:
                print(f"    {etype}: {len(elist)}")

        # Production warning if lots of evidence missing reference data
        if rejection_stats['missing_reference_data'] > 10:
            print(f"  ⚠️  WARNING: {rejection_stats['missing_reference_data']} evidence items missing reference data")
            print(f"     This suggests profile library sync issues - consider updating profile libraries")

    return filtered_evidence


def safe_extract_chain_id(evidence):
    """Safely extract chain ID from evidence, handling None domain_id"""
    if evidence.domain_id and isinstance(evidence.domain_id, str):
        if '_' in evidence.domain_id:
            return evidence.domain_id.split('_')[-1]
        elif len(evidence.domain_id) > 5:
            return evidence.domain_id[5]  # Extract from domain_id like "e6dgvA1" -> "A"
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

            domain = create_domain_with_provenance(
                evidence=dec_evidence,
                domain_id=f"d{len(selected_domains) + 1}",
                domain_definitions=domain_definitions
            )

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

            # ENHANCED: Create domain with COMPLETE provenance tracking
            domain = create_domain_with_provenance(
                evidence=evidence,
                domain_id=f"d{len(selected_domains) + 1}",
                domain_definitions=domain_definitions
            )

            # Block residues
            used_positions.update(evidence_positions)
            unused_positions.difference_update(evidence_positions)
            selected_domains.append(domain)

            if verbose:
                print(f"  ✓ Selected {domain.id}: {domain.family} @ {domain.range}")
                print(f"    Coverage: {new_coverage:.1%} new, {used_coverage:.1%} overlap")
                print(f"    Provenance: {evidence.type} -> {evidence.domain_id}")

                # ENHANCED: Show quality metrics
                if evidence.reference_length and evidence.hit_range:
                    ref_coverage = evidence.hit_range.total_length / evidence.reference_length
                    print(f"    Reference coverage: {ref_coverage:.1%} ({evidence.hit_range.total_length}/{evidence.reference_length})")
                print(f"    Confidence: {evidence.confidence:.3f}")

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


# ENHANCED: Convenience function to partition with quality control
def partition_domains_with_quality_control(evidence_list: List['Evidence'],
                                         sequence_length: int,
                                         domain_definitions: Dict = None,
                                         strict_mode: bool = True,
                                         verbose: bool = False) -> List['Domain']:
    """
    Partition domains with recommended quality thresholds

    Args:
        evidence_list: List of evidence
        sequence_length: Protein sequence length
        domain_definitions: Domain definitions for decomposition
        strict_mode: Whether to apply strict quality thresholds
        verbose: Detailed output

    Returns:
        List of high-quality domains
    """

    return partition_domains(
        evidence_list=evidence_list,
        sequence_length=sequence_length,
        domain_definitions=domain_definitions,
        apply_quality_thresholds=strict_mode,
        verbose=verbose
    )
