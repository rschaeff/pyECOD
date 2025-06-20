# mini/core/domain_utils.py
"""
Utility functions for consistent domain creation and provenance tracking

This module provides standardized functions to ensure consistent
domain creation, classification assignment, and provenance tracking
across all evidence processing paths.
"""

from typing import Optional, Dict, List, Any
from datetime import datetime
from .models import Domain, Evidence
from .sequence_range import SequenceRange


def get_evidence_classification(evidence: Evidence, 
                               domain_definitions: Dict = None) -> Dict[str, Optional[str]]:
    """
    Get ECOD taxonomic classification for evidence with standardized fallbacks.
    
    This ensures consistent classification assignment across all evidence types.
    
    Args:
        evidence: Evidence object
        domain_definitions: Domain definitions for lookup
        
    Returns:
        Dict with x_group, h_group, t_group classification
    """
    # Method 1: Direct T-group from evidence (highest priority)
    if evidence.t_group:
        x_group, h_group, t_group = parse_ecod_hierarchy(evidence.t_group)
        return {
            'x_group': x_group,
            'h_group': h_group,
            't_group': t_group
        }
    
    # Method 2: Lookup domain_id in domain_definitions
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
    
    # Method 3: Use H-group from evidence if available
    if evidence.h_group:
        return {
            'x_group': None,
            'h_group': evidence.h_group,
            't_group': None
        }
    
    # Fallback: unclassified
    return {
        'x_group': None,
        'h_group': None,
        't_group': None
    }


def parse_ecod_hierarchy(t_group_str: str) -> tuple:
    """
    Parse ECOD T-group into hierarchical components.
    
    Args:
        t_group_str: T-group string like "1.1.1"
        
    Returns:
        Tuple of (x_group, h_group, t_group)
    """
    if not t_group_str:
        return None, None, None

    parts = t_group_str.split('.')
    if len(parts) >= 3:
        x_group = parts[0]
        h_group = f"{parts[0]}.{parts[1]}"
        t_group = f"{parts[0]}.{parts[1]}.{parts[2]}"
        return x_group, h_group, t_group
    elif len(parts) == 2:
        x_group = parts[0]
        h_group = f"{parts[0]}.{parts[1]}"
        return x_group, h_group, None
    elif len(parts) == 1:
        return parts[0], None, None

    return None, None, None


def get_domain_family_name(evidence: Evidence, 
                          classification: Dict[str, Optional[str]]) -> str:
    """
    Determine domain family name with standardized priority logic.
    
    Priority order:
    1. T-group (most specific classification)
    2. source_pdb (for test compatibility and provenance)  
    3. domain_id (maintains reference tracking)
    4. 'unclassified' (fallback)
    
    Args:
        evidence: Evidence object
        classification: Classification dict from get_evidence_classification()
        
    Returns:
        Family name string
    """
    # Priority 1: T-group if available (best classification)
    if classification.get('t_group'):
        return classification['t_group']
    
    # Priority 2: source_pdb for better test compatibility
    if evidence.source_pdb:
        return evidence.source_pdb
    
    # Priority 3: domain_id if available (maintains reference tracking)
    if evidence.domain_id:
        return evidence.domain_id
    
    # Fallback
    return 'unclassified'


def create_domain_with_provenance(evidence: Evidence,
                                 domain_id: str,
                                 classification: Dict[str, Optional[str]] = None,
                                 domain_definitions: Dict = None) -> Domain:
    """
    Create a Domain object with complete and consistent provenance tracking.
    
    This function ensures all domains are created with the same provenance
    information regardless of which evidence processing path created them.
    
    Args:
        evidence: Primary evidence that created this domain
        domain_id: Domain identifier (e.g., "d1", "d2")
        classification: Optional pre-computed classification
        domain_definitions: Optional domain definitions for classification lookup
        
    Returns:
        Domain object with complete provenance
    """
    # Get classification if not provided
    if classification is None:
        classification = get_evidence_classification(evidence, domain_definitions)
    
    # Determine family name using standardized logic
    family_name = get_domain_family_name(evidence, classification)
    
    # Create domain with comprehensive provenance
    domain = Domain(
        id=domain_id,
        range=evidence.query_range,
        family=family_name,
        evidence_count=1,
        source=evidence.type,
        evidence_items=[evidence],
        
        # ECOD hierarchy classification
        x_group=classification.get('x_group'),
        h_group=classification.get('h_group'),
        t_group=classification.get('t_group'),
        
        # Comprehensive provenance tracking
        primary_evidence=evidence,
        reference_ecod_domain_id=evidence.domain_id,
        original_range=evidence.query_range,  # Before any optimization
        confidence_score=evidence.confidence,
        creation_timestamp=datetime.now()
    )
    
    return domain


def merge_domain_evidence(domain: Domain, additional_evidence: Evidence) -> Domain:
    """
    Merge additional evidence into an existing domain with provenance tracking.
    
    Args:
        domain: Existing domain
        additional_evidence: Evidence to merge in
        
    Returns:
        Updated domain with merged evidence
    """
    # Add evidence to list
    domain.evidence_items.append(additional_evidence)
    domain.evidence_count = len(domain.evidence_items)
    
    # Update primary evidence if this one is better
    if additional_evidence.confidence > domain.primary_evidence.confidence:
        domain.primary_evidence = additional_evidence
        domain.confidence_score = additional_evidence.confidence
        
        # Update reference ID if new evidence has better reference
        if additional_evidence.domain_id:
            domain.reference_ecod_domain_id = additional_evidence.domain_id
    
    # Record the merge action for provenance
    domain.record_optimization_action(
        "evidence_merge", 
        f"added_{additional_evidence.type}_{additional_evidence.confidence:.2f}"
    )
    
    return domain


def validate_domain_provenance(domain: Domain) -> tuple:
    """
    Validate that a domain has complete provenance information.
    
    Args:
        domain: Domain object to validate
        
    Returns:
        Tuple of (is_valid, list_of_issues)
    """
    issues = []
    
    # Core domain fields
    if not domain.id:
        issues.append("Missing domain ID")
    
    if not domain.range:
        issues.append("Missing domain range")
    
    if not domain.family:
        issues.append("Missing family assignment")
    
    if domain.evidence_count <= 0:
        issues.append("No evidence items")
    
    # Provenance fields
    if not domain.primary_evidence:
        issues.append("Missing primary evidence")
    
    if not domain.original_range:
        issues.append("Missing original range")
    
    if domain.confidence_score is None or domain.confidence_score <= 0:
        issues.append("Missing or invalid confidence score")
    
    if not domain.creation_timestamp:
        issues.append("Missing creation timestamp")
    
    # Evidence consistency
    if domain.evidence_items and len(domain.evidence_items) != domain.evidence_count:
        issues.append("Evidence count mismatch")
    
    if (domain.primary_evidence and 
        domain.evidence_items and 
        domain.primary_evidence not in domain.evidence_items):
        issues.append("Primary evidence not in evidence list")
    
    # Range consistency
    if domain.original_range and domain.range:
        original_positions = set(domain.original_range.to_positions_simple())
        current_positions = domain.get_positions()
        
        # Check if optimization tracking is consistent
        if original_positions != current_positions and not domain.was_optimized():
            issues.append("Range changed but no optimization recorded")
    
    return len(issues) == 0, issues


def create_domain_summary_dict(domain: Domain) -> Dict[str, Any]:
    """
    Create a comprehensive summary dictionary of domain with provenance.
    
    Useful for debugging, logging, and test validation.
    
    Args:
        domain: Domain object
        
    Returns:
        Dictionary with all domain information
    """
    # Basic validation
    is_valid, validation_issues = validate_domain_provenance(domain)
    
    summary = {
        # Core fields
        'id': domain.id,
        'range': str(domain.range),
        'family': domain.family,
        'source': domain.source,
        'evidence_count': domain.evidence_count,
        'confidence': domain.confidence_score,
        
        # Classification
        't_group': domain.t_group,
        'h_group': domain.h_group,
        'x_group': domain.x_group,
        
        # Provenance
        'creation_timestamp': domain.creation_timestamp.isoformat() if domain.creation_timestamp else None,
        'reference_ecod_domain_id': domain.reference_ecod_domain_id,
        'original_range': str(domain.original_range) if domain.original_range else None,
        
        # Optimization tracking
        'was_optimized': domain.was_optimized(),
        'optimization_actions': domain.optimization_actions.copy() if domain.optimization_actions else [],
        
        # Statistics
        'length': domain.length,
        'is_discontinuous': domain.range.is_discontinuous if domain.range else False,
        
        # Primary evidence summary
        'primary_evidence_type': domain.primary_evidence.type if domain.primary_evidence else None,
        'primary_evidence_source': domain.primary_evidence.source_pdb if domain.primary_evidence else None,
        
        # Validation
        'is_valid': is_valid,
        'validation_issues': validation_issues
    }
    
    return summary


def compare_domains_for_consistency(domain1: Domain, domain2: Domain) -> Dict[str, Any]:
    """
    Compare two domains for consistency analysis.
    
    Useful for testing and validation to ensure domains created
    from similar evidence have consistent provenance.
    
    Args:
        domain1: First domain
        domain2: Second domain
        
    Returns:
        Dictionary with comparison results
    """
    comparison = {
        'ranges_match': domain1.range == domain2.range,
        'families_match': domain1.family == domain2.family,
        'classifications_match': (
            domain1.t_group == domain2.t_group and
            domain1.h_group == domain2.h_group and
            domain1.x_group == domain2.x_group
        ),
        'confidence_difference': abs((domain1.confidence_score or 0) - (domain2.confidence_score or 0)),
        'both_have_provenance': (
            validate_domain_provenance(domain1)[0] and
            validate_domain_provenance(domain2)[0]
        )
    }
    
    # Detailed differences
    differences = []
    
    if not comparison['ranges_match']:
        differences.append(f"Range: {domain1.range} vs {domain2.range}")
    
    if not comparison['families_match']:
        differences.append(f"Family: {domain1.family} vs {domain2.family}")
    
    if not comparison['classifications_match']:
        differences.append(f"Classification: {domain1.t_group} vs {domain2.t_group}")
    
    if comparison['confidence_difference'] > 0.1:
        differences.append(f"Confidence: {domain1.confidence_score} vs {domain2.confidence_score}")
    
    comparison['differences'] = differences
    comparison['is_consistent'] = len(differences) == 0
    
    return comparison


def standardize_domain_list(domains: List[Domain]) -> List[Domain]:
    """
    Apply consistency checks and standardization to a list of domains.
    
    Args:
        domains: List of Domain objects
        
    Returns:
        List of domains with standardized provenance
    """
    standardized = []
    
    for i, domain in enumerate(domains):
        # Ensure domain has proper ID if missing
        if not domain.id:
            domain.id = f"d{i+1}"
        
        # Ensure creation timestamp if missing
        if not domain.creation_timestamp:
            domain.creation_timestamp = datetime.now()
        
        # Ensure confidence score if missing
        if domain.confidence_score is None and domain.primary_evidence:
            domain.confidence_score = domain.primary_evidence.confidence
        
        # Validate provenance
        is_valid, issues = validate_domain_provenance(domain)
        if not is_valid:
            print(f"Warning: Domain {domain.id} has provenance issues: {issues}")
        
        standardized.append(domain)
    
    return standardized


# Convenience functions for common domain operations
def get_domain_coverage_stats(domains: List[Domain], sequence_length: int) -> Dict[str, Any]:
    """Get coverage statistics for a list of domains"""
    if not domains:
        return {
            'total_domains': 0,
            'total_coverage': 0,
            'coverage_percentage': 0.0
        }
    
    total_coverage = sum(d.length for d in domains)
    optimized_count = sum(1 for d in domains if d.was_optimized())
    
    with_classification = sum(1 for d in domains if d.t_group)
    high_confidence = sum(1 for d in domains if (d.confidence_score or 0) >= 0.7)
    
    return {
        'total_domains': len(domains),
        'total_coverage': total_coverage,
        'coverage_percentage': (total_coverage / sequence_length * 100) if sequence_length > 0 else 0,
        'optimized_domains': optimized_count,
        'optimization_percentage': (optimized_count / len(domains) * 100),
        'with_classification': with_classification,
        'classification_percentage': (with_classification / len(domains) * 100),
        'high_confidence': high_confidence,
        'high_confidence_percentage': (high_confidence / len(domains) * 100)
    }


def find_domain_inconsistencies(domains: List[Domain]) -> Dict[str, Any]:
    """Find potential inconsistencies in a domain list"""
    issues = {
        'duplicate_ranges': [],
        'missing_provenance': [],
        'confidence_outliers': [],
        'classification_gaps': []
    }
    
    # Check for duplicate ranges
    ranges_seen = {}
    for domain in domains:
        range_str = str(domain.range)
        if range_str in ranges_seen:
            issues['duplicate_ranges'].append({
                'range': range_str,
                'domains': [ranges_seen[range_str], domain.id]
            })
        else:
            ranges_seen[range_str] = domain.id
    
    # Check for missing provenance
    for domain in domains:
        is_valid, validation_issues = validate_domain_provenance(domain)
        if not is_valid:
            issues['missing_provenance'].append({
                'domain_id': domain.id,
                'issues': validation_issues
            })
    
    # Check for confidence outliers
    confidences = [d.confidence_score for d in domains if d.confidence_score is not None]
    if confidences:
        avg_confidence = sum(confidences) / len(confidences)
        for domain in domains:
            if (domain.confidence_score is not None and 
                abs(domain.confidence_score - avg_confidence) > 0.3):
                issues['confidence_outliers'].append({
                    'domain_id': domain.id,
                    'confidence': domain.confidence_score,
                    'average': avg_confidence
                })
    
    # Check for classification gaps
    for domain in domains:
        if not domain.t_group and domain.primary_evidence:
            issues['classification_gaps'].append({
                'domain_id': domain.id,
                'evidence_type': domain.primary_evidence.type,
                'source': domain.primary_evidence.source_pdb
            })
    
    return issues
