# mini/core/evidence_utils.py
"""
Enhanced utility functions for consistent evidence handling and provenance tracking

This module provides standardized functions to ensure consistent
evidence processing, confidence calculation, reference coverage tracking,
and provenance field population across all evidence types.

ENHANCED: Adds comprehensive reference coverage calculation and quality assessment
"""

from typing import Optional, Dict, Tuple, Any, List
from .models import Evidence, AlignmentData
from .sequence_range import SequenceRange


def calculate_evidence_confidence(evalue: Optional[float] = None,
                                 probability: Optional[float] = None,
                                 evidence_type: str = "unknown",
                                 alignment_coverage: Optional[float] = None,
                                 reference_coverage: Optional[float] = None) -> float:
    """
    Standardized confidence calculation for all evidence types.

    This ensures consistent confidence scoring across chain_blast,
    domain_blast, and hhsearch evidence with enhanced quality weighting.

    Args:
        evalue: E-value from BLAST or HHsearch
        probability: Probability score (typically from HHsearch, 0-100 scale)
        evidence_type: Type of evidence for type-specific adjustments
        alignment_coverage: Coverage fraction for additional weighting
        reference_coverage: Coverage of reference domain for quality assessment

    Returns:
        Confidence score between 0.05 and 0.95
    """
    base_confidence = 0.1  # Minimum confidence for any evidence

    # Primary confidence from e-value (for BLAST evidence)
    if evalue is not None and evalue > 0:
        if evalue < 1e-10:
            base_confidence = 0.9
        elif evalue < 1e-8:
            base_confidence = 0.8
        elif evalue < 1e-5:
            base_confidence = 0.7
        elif evalue < 1e-3:
            base_confidence = 0.6
        elif evalue < 0.01:
            base_confidence = 0.5
        elif evalue < 0.1:
            base_confidence = 0.4
        elif evalue < 1.0:
            base_confidence = 0.3
        else:
            base_confidence = 0.2

    # Alternative confidence from probability (HHsearch)
    elif probability is not None:
        if probability > 1.0:
            # HHsearch probability is 0-100 scale, normalize to 0-1
            probability = probability / 100.0
        base_confidence = max(0.1, min(0.9, probability))

    # Type-specific adjustments
    type_multipliers = {
        'domain_blast': 1.0,      # Domain BLAST is most reliable
        'hhsearch': 0.95,         # HHsearch is very good for remote homology
        'chain_blast': 0.9,       # Chain BLAST needs decomposition
        'chain_blast_decomposed': 0.85  # Decomposed chain BLAST has additional uncertainty
    }

    type_multiplier = type_multipliers.get(evidence_type, 0.8)
    adjusted_confidence = base_confidence * type_multiplier

    # ENHANCED: Reference coverage boost for high-quality matches
    if reference_coverage is not None and reference_coverage > 0.7:
        coverage_boost = min(0.1, (reference_coverage - 0.7) * 0.2)
        adjusted_confidence += coverage_boost

    # ENHANCED: Penalty for poor reference coverage
    elif reference_coverage is not None and reference_coverage < 0.3:
        coverage_penalty = min(0.1, (0.3 - reference_coverage) * 0.2)
        adjusted_confidence -= coverage_penalty

    # Query coverage boost for high-coverage alignments (legacy)
    if alignment_coverage is not None and alignment_coverage > 0.7:
        query_coverage_boost = min(0.05, (alignment_coverage - 0.7) * 0.1)
        adjusted_confidence += query_coverage_boost

    # Ensure confidence stays in valid range
    return max(0.05, min(0.95, adjusted_confidence))


def calculate_reference_coverage(evidence: Evidence) -> Optional[float]:
    """
    Calculate reference coverage for evidence with hit range and reference length.

    Args:
        evidence: Evidence object with hit_range and reference_length

    Returns:
        Coverage fraction (0.0-1.0) or None if cannot be calculated
    """
    if not evidence.hit_range or not evidence.reference_length:
        return None

    if evidence.reference_length <= 0:
        return None

    hit_length = evidence.hit_range.total_length
    return hit_length / evidence.reference_length


def populate_evidence_provenance(evidence: Evidence,
                                hit_range: Optional[SequenceRange] = None,
                                alignment: Optional[AlignmentData] = None,
                                reference_length: Optional[int] = None,
                                domain_id: Optional[str] = None,
                                classification: Optional[Dict[str, str]] = None) -> Evidence:
    """
    Populate provenance fields consistently for any evidence type.

    ENHANCED: Includes comprehensive reference coverage calculation and quality assessment.

    Args:
        evidence: Evidence object to populate
        hit_range: Range in the hit/reference sequence
        alignment: Alignment data for chain BLAST
        reference_length: Length of reference domain/protein
        domain_id: Reference domain ID
        classification: ECOD classification dict

    Returns:
        Evidence object with complete provenance fields
    """
    # Set basic provenance fields
    if hit_range is not None:
        evidence.hit_range = hit_range

    if alignment is not None:
        evidence.alignment = alignment
        evidence.hsp_count = 1  # Most alignments have one HSP
    elif evidence.hsp_count is None:
        evidence.hsp_count = 1  # Default assumption

    if reference_length is not None:
        evidence.reference_length = reference_length

    if domain_id is not None:
        evidence.domain_id = domain_id

    # Set discontinuous flag based on query range
    evidence.discontinuous = evidence.query_range.is_discontinuous

    # ENHANCED: Calculate reference coverage if we have both hit range and reference length
    reference_coverage = calculate_reference_coverage(evidence)
    if reference_coverage is not None:
        evidence.reference_coverage = reference_coverage

        # Also set alignment_coverage for backward compatibility
        # (though reference_coverage is more accurate)
        if evidence.alignment_coverage is None:
            evidence.alignment_coverage = reference_coverage

    # Calculate query coverage if we have hit range but not reference coverage
    elif (evidence.hit_range is not None and
          evidence.reference_length is not None and
          evidence.reference_length > 0):

        # Fallback query-based coverage calculation
        hit_coverage = evidence.hit_range.total_length / evidence.reference_length
        evidence.alignment_coverage = hit_coverage

    # Set classification fields
    if classification:
        evidence.t_group = classification.get('t_group')
        evidence.h_group = classification.get('h_group')

    # Extract chain ID if not already set
    if not evidence.source_chain_id:
        evidence.source_chain_id = extract_chain_id_from_evidence(evidence)

    # ENHANCED: Recalculate confidence with all available information including reference coverage
    evidence.confidence = calculate_evidence_confidence(
        evalue=evidence.evalue,
        evidence_type=evidence.type,
        alignment_coverage=evidence.alignment_coverage,
        reference_coverage=reference_coverage
    )

    return evidence


def extract_chain_id_from_evidence(evidence: Evidence) -> Optional[str]:
    """
    Extract chain ID from evidence using multiple fallback methods.

    Args:
        evidence: Evidence object

    Returns:
        Chain ID string or None if not determinable
    """
    # Method 1: Already set
    if evidence.source_chain_id:
        return evidence.source_chain_id

    # Method 2: Parse from domain_id
    if evidence.domain_id:
        # Handle formats like "e6dgvA1" -> "A" or "6dgv_A" -> "A"
        domain_id = evidence.domain_id

        if '_' in domain_id:
            # Format: "6dgv_A" or "pdb_chain"
            parts = domain_id.split('_')
            if len(parts) >= 2:
                return parts[-1]  # Last part is typically chain

        elif len(domain_id) > 5 and domain_id.startswith('e'):
            # Format: "e6dgvA1" -> extract chain from position 5
            potential_chain = domain_id[5]
            if potential_chain.isalpha():
                return potential_chain

    # Method 3: Default fallback
    return 'A'


def assess_evidence_quality(evidence: Evidence) -> Dict[str, Any]:
    """
    Comprehensive quality assessment for evidence.

    ENHANCED: Includes reference coverage, confidence, and E-value quality metrics.

    Args:
        evidence: Evidence object to assess

    Returns:
        Dictionary with quality assessment results
    """
    assessment = {
        'overall_quality': 'unknown',
        'quality_score': 0.0,
        'issues': [],
        'strengths': [],
        'recommendations': []
    }

    # Assess confidence
    if evidence.confidence >= 0.8:
        assessment['strengths'].append('high_confidence')
        assessment['quality_score'] += 3
    elif evidence.confidence >= 0.6:
        assessment['quality_score'] += 2
    elif evidence.confidence < 0.5:
        assessment['issues'].append('low_confidence')
        assessment['quality_score'] -= 2

    # Assess reference coverage
    ref_coverage = calculate_reference_coverage(evidence)
    if ref_coverage is not None:
        if ref_coverage >= 0.8:
            assessment['strengths'].append('excellent_reference_coverage')
            assessment['quality_score'] += 3
        elif ref_coverage >= 0.6:
            assessment['strengths'].append('good_reference_coverage')
            assessment['quality_score'] += 2
        elif ref_coverage >= 0.4:
            assessment['quality_score'] += 1
        else:
            assessment['issues'].append('poor_reference_coverage')
            assessment['quality_score'] -= 3
            assessment['recommendations'].append('Consider requiring higher reference coverage')
    else:
        assessment['issues'].append('missing_reference_coverage')
        assessment['quality_score'] -= 1

    # Assess E-value (for BLAST evidence)
    if evidence.evalue is not None:
        if evidence.evalue < 1e-10:
            assessment['strengths'].append('excellent_evalue')
            assessment['quality_score'] += 2
        elif evidence.evalue < 1e-5:
            assessment['strengths'].append('good_evalue')
            assessment['quality_score'] += 1
        elif evidence.evalue > 1.0:
            assessment['issues'].append('poor_evalue')
            assessment['quality_score'] -= 2

    # Type-specific assessments
    if evidence.type == 'hhsearch':
        # Calculate original probability
        original_prob = evidence.confidence / 0.95 * 100  # Reverse type multiplier
        if original_prob < 50:
            assessment['issues'].append('low_hhsearch_probability')
            assessment['recommendations'].append('Consider minimum 68% HHsearch probability threshold')
        elif original_prob >= 80:
            assessment['strengths'].append('high_hhsearch_probability')

    # Overall quality determination
    if assessment['quality_score'] >= 5:
        assessment['overall_quality'] = 'excellent'
    elif assessment['quality_score'] >= 3:
        assessment['overall_quality'] = 'good'
    elif assessment['quality_score'] >= 1:
        assessment['overall_quality'] = 'acceptable'
    elif assessment['quality_score'] >= -1:
        assessment['overall_quality'] = 'questionable'
    else:
        assessment['overall_quality'] = 'poor'

    return assessment


def validate_evidence_provenance(evidence: Evidence, strict_mode: bool = True) -> Tuple[bool, List[str]]:
    """
    Validate that evidence has complete provenance information.

    Args:
        evidence: Evidence object to validate
        strict_mode: If True (default), require reference coverage for domain_blast/hhsearch
                    If False, allow missing reference data during parsing

    Returns:
        Tuple of (is_valid, list_of_issues)
    """
    issues = []

    # Required fields for all evidence
    if not evidence.source_pdb:
        issues.append("Missing source_pdb")

    if not evidence.query_range:
        issues.append("Missing query_range")

    if evidence.confidence is None or evidence.confidence <= 0:
        issues.append("Missing or invalid confidence score")

    # Provenance fields that should be set
    if not evidence.source_chain_id:
        issues.append("Missing source_chain_id")

    if evidence.hsp_count is None:
        issues.append("Missing hsp_count")

    # Type-specific validation - STRICT by default
    if evidence.type == "domain_blast":
        if not evidence.domain_id:
            issues.append("Domain BLAST evidence missing domain_id")

        if strict_mode:
            if evidence.reference_length is None:
                issues.append("Domain BLAST evidence missing reference_length")

            # STRICT MODE: Reference coverage is REQUIRED
            ref_coverage = calculate_reference_coverage(evidence)
            if ref_coverage is None:
                issues.append("Domain BLAST evidence missing reference coverage data")

    elif evidence.type == "chain_blast":
        if strict_mode and evidence.reference_length is None:
            issues.append("Chain BLAST evidence missing protein length")

    elif evidence.type == "hhsearch":
        if not evidence.domain_id:
            issues.append("HHsearch evidence missing domain_id")

        if strict_mode:
            if evidence.reference_length is None:
                issues.append("HHsearch evidence missing reference_length")

            # STRICT MODE: Reference coverage is REQUIRED
            ref_coverage = calculate_reference_coverage(evidence)
            if ref_coverage is None:
                issues.append("HHsearch evidence missing reference coverage data")

    # Coverage validation (basic range checks)
    if (evidence.alignment_coverage is not None and
        (evidence.alignment_coverage < 0 or evidence.alignment_coverage > 1)):
        issues.append(f"Invalid alignment_coverage: {evidence.alignment_coverage}")

    # Reference coverage validation (basic range checks)
    ref_coverage = calculate_reference_coverage(evidence)
    if ref_coverage is not None and (ref_coverage < 0 or ref_coverage > 1):
        issues.append(f"Invalid reference_coverage: {ref_coverage}")

    return len(issues) == 0, issues


def create_evidence_summary_dict(evidence: Evidence) -> Dict[str, Any]:
    """
    Create a summary dictionary of evidence with all provenance information.

    ENHANCED: Includes reference coverage and quality assessment.

    Args:
        evidence: Evidence object

    Returns:
        Dictionary with all evidence information
    """
    # Calculate reference coverage
    ref_coverage = calculate_reference_coverage(evidence)

    # Get quality assessment
    quality = assess_evidence_quality(evidence)

    return {
        # Core fields
        'type': evidence.type,
        'source_pdb': evidence.source_pdb,
        'query_range': str(evidence.query_range),
        'confidence': evidence.confidence,
        'evalue': evidence.evalue,

        # Provenance fields
        'source_chain_id': evidence.source_chain_id,
        'domain_id': evidence.domain_id,
        'hit_range': str(evidence.hit_range) if evidence.hit_range else None,
        'hsp_count': evidence.hsp_count,
        'discontinuous': evidence.discontinuous,
        'reference_length': evidence.reference_length,
        'alignment_coverage': evidence.alignment_coverage,

        # ENHANCED: Reference coverage
        'reference_coverage': ref_coverage,
        'reference_coverage_percent': f"{ref_coverage:.1%}" if ref_coverage else None,

        # Classification
        't_group': evidence.t_group,
        'h_group': evidence.h_group,

        # ENHANCED: Quality assessment
        'quality_assessment': quality['overall_quality'],
        'quality_score': quality['quality_score'],
        'quality_issues': quality['issues'],
        'quality_strengths': quality['strengths'],

        # Validation
        'is_valid': validate_evidence_provenance(evidence)[0]
    }


def standardize_evidence_list(evidence_list: list,
                             reference_lengths: Dict[str, int] = None,
                             protein_lengths: Dict[Tuple[str, str], int] = None,
                             domain_definitions: Dict = None) -> list:
    """
    Apply standardized provenance population to a list of evidence.

    ENHANCED: Includes reference coverage calculation for all evidence types.

    Args:
        evidence_list: List of Evidence objects
        reference_lengths: Domain lengths lookup
        protein_lengths: Protein lengths lookup
        domain_definitions: Domain definitions for classification

    Returns:
        List of Evidence objects with standardized provenance
    """
    standardized = []

    # PERFORMANCE FIX: Create reverse lookup for domain classifications (O(1) lookup instead of O(n*m))
    domain_id_to_classification = {}
    if domain_definitions:
        for domain_refs in domain_definitions.values():
            for ref in domain_refs:
                if ref.domain_id and ref.t_group:
                    domain_id_to_classification[ref.domain_id] = {
                        't_group': ref.t_group,
                        'h_group': ref.h_group
                    }

    for evidence in evidence_list:
        # Look up reference length if missing
        if evidence.reference_length is None:
            if evidence.type == "domain_blast" and reference_lengths:
                if evidence.domain_id in reference_lengths:
                    evidence.reference_length = reference_lengths[evidence.domain_id]

            elif evidence.type == "chain_blast" and protein_lengths:
                chain_id = extract_chain_id_from_evidence(evidence)
                lookup_key = (evidence.source_pdb, chain_id)
                if lookup_key in protein_lengths:
                    evidence.reference_length = protein_lengths[lookup_key]

        # PERFORMANCE FIX: O(1) classification lookup
        classification = None
        if evidence.domain_id and evidence.domain_id in domain_id_to_classification:
            classification = domain_id_to_classification[evidence.domain_id]

        # Apply standardized provenance population
        evidence = populate_evidence_provenance(
            evidence=evidence,
            reference_length=evidence.reference_length,
            classification=classification
        )

        standardized.append(evidence)

    return standardized


# ENHANCED: Quality-based filtering functions
def filter_evidence_by_quality(evidence_list: List[Evidence],
                               min_confidence: float = 0.5,
                               min_reference_coverage: float = 0.5,
                               require_reference_data: bool = True) -> Tuple[List[Evidence], Dict[str, int]]:
    """
    Filter evidence list by quality thresholds.

    Args:
        evidence_list: List of evidence to filter
        min_confidence: Minimum confidence threshold
        min_reference_coverage: Minimum reference coverage threshold
        require_reference_data: Whether to require reference coverage data

    Returns:
        Tuple of (filtered_evidence, rejection_stats)
    """
    filtered = []
    stats = {'confidence': 0, 'reference_coverage': 0, 'missing_data': 0}

    for evidence in evidence_list:
        # Check confidence
        if evidence.confidence < min_confidence:
            stats['confidence'] += 1
            continue

        # Check reference coverage
        ref_coverage = calculate_reference_coverage(evidence)
        if ref_coverage is not None:
            if ref_coverage < min_reference_coverage:
                stats['reference_coverage'] += 1
                continue
        elif require_reference_data:
            stats['missing_data'] += 1
            continue

        filtered.append(evidence)

    return filtered, stats


# Convenience functions for common operations
def is_high_confidence_evidence(evidence: Evidence, threshold: float = 0.7) -> bool:
    """Check if evidence meets high confidence threshold"""
    return evidence.confidence >= threshold


def is_high_coverage_evidence(evidence: Evidence, threshold: float = 0.6) -> bool:
    """Check if evidence has high reference coverage"""
    ref_coverage = calculate_reference_coverage(evidence)
    return ref_coverage is not None and ref_coverage >= threshold


def get_evidence_coverage_stats(evidence_list: list) -> Dict[str, Any]:
    """
    Get coverage statistics for a list of evidence.

    ENHANCED: Includes reference coverage statistics.
    """
    total = len(evidence_list)
    if total == 0:
        return {'total': 0}

    with_provenance = sum(1 for e in evidence_list if validate_evidence_provenance(e)[0])
    high_confidence = sum(1 for e in evidence_list if is_high_confidence_evidence(e))
    high_coverage = sum(1 for e in evidence_list if is_high_coverage_evidence(e))

    # Reference coverage statistics
    with_ref_coverage = sum(1 for e in evidence_list if calculate_reference_coverage(e) is not None)
    avg_ref_coverage = None
    if with_ref_coverage > 0:
        coverages = [calculate_reference_coverage(e) for e in evidence_list
                    if calculate_reference_coverage(e) is not None]
        avg_ref_coverage = sum(coverages) / len(coverages)

    by_type = {}
    for evidence in evidence_list:
        by_type[evidence.type] = by_type.get(evidence.type, 0) + 1

    return {
        'total': total,
        'with_complete_provenance': with_provenance,
        'provenance_percentage': (with_provenance / total) * 100,
        'high_confidence': high_confidence,
        'high_confidence_percentage': (high_confidence / total) * 100,
        'high_coverage': high_coverage,
        'high_coverage_percentage': (high_coverage / total) * 100,
        'with_reference_coverage': with_ref_coverage,
        'average_reference_coverage': avg_ref_coverage,
        'by_type': by_type
    }
