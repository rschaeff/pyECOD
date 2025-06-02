# mini/parser.py
"""Parse domain summary XML with reference length support"""

import xml.etree.ElementTree as ET
import csv
from typing import List, Dict, Optional
from .models import Evidence
from ecod.core.sequence_range import SequenceRange

def load_reference_lengths(csv_path: str) -> Dict[str, int]:
    """
    Load domain reference lengths from database export CSV

    Expected CSV format:
    domain_id,length
    e7nwgd21,62
    e7nkyO3,102
    ...
    """
    lengths = {}
    try:
        with open(csv_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                domain_id = row['domain_id']
                length = int(row['length'])
                lengths[domain_id] = length
        print(f"Loaded {len(lengths)} reference domain lengths")
    except FileNotFoundError:
        print(f"Warning: Reference file not found: {csv_path}")
    except Exception as e:
        print(f"Warning: Error loading reference lengths: {e}")

    return lengths

def parse_domain_summary(xml_path: str, ref_lengths: Optional[Dict[str, int]] = None) -> List[Evidence]:
    """
    Parse evidence from domain summary XML

    Args:
        xml_path: Path to domain summary XML file
        ref_lengths: Optional dictionary of domain_id -> length mappings

    Returns:
        List of Evidence objects
    """
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
    except Exception as e:
        print(f"Error parsing XML {xml_path}: {e}")
        return []

    evidence_list = []
    ref_lengths = ref_lengths or {}

    # Parse chain BLAST hits
    for hit in root.findall(".//chain_blast_run/hits/hit"):
        evidence = _parse_chain_blast_hit(hit)
        if evidence:
            evidence_list.append(evidence)

    # Parse domain BLAST hits
    for hit in root.findall(".//blast_run/hits/hit"):
        evidence = _parse_domain_blast_hit(hit, ref_lengths)
        if evidence:
            evidence_list.append(evidence)

    # Parse HHSearch hits
    for hit in root.findall(".//hh_run/hits/hit"):
        evidence = _parse_hhsearch_hit(hit, ref_lengths)
        if evidence:
            evidence_list.append(evidence)

    print(f"Parsed {len(evidence_list)} evidence items from {xml_path}")
    return evidence_list

def _parse_chain_blast_hit(hit_elem: ET.Element) -> Optional[Evidence]:
    """Parse a chain BLAST hit"""
    pdb_id = hit_elem.get("pdb_id", "")
    chain_id = hit_elem.get("chain_id", "")

    query_reg = hit_elem.find("query_reg")
    if query_reg is None or not query_reg.text:
        return None

    try:
        query_range = SequenceRange.parse(query_reg.text.strip())
    except ValueError as e:
        print(f"Invalid range in chain BLAST hit: {query_reg.text} - {e}")
        return None

    # Parse e-value
    evalue = _safe_float(hit_elem.get("evalues", "999"), 999.0)

    # Chain BLAST doesn't have reference coverage, use e-value for confidence
    confidence = _evalue_to_confidence(evalue) * 0.8  # Penalty for chain-level

    return Evidence(
        type="chain_blast",
        source_pdb=pdb_id,
        query_range=query_range,
        evalue=evalue,
        confidence=confidence,
        # Chain hits don't have domain IDs
        domain_id=f"{pdb_id}_{chain_id}" if pdb_id and chain_id else ""
    )

def _parse_domain_blast_hit(hit_elem: ET.Element, ref_lengths: Dict[str, int]) -> Optional[Evidence]:
    """Parse a domain BLAST hit with reference coverage"""
    domain_id = hit_elem.get("domain_id", "")
    if not domain_id:
        return None

    # Extract PDB from domain ID (e.g., 'e6dgvA1' -> '6dgv')
    source_pdb = domain_id[1:5] if len(domain_id) > 4 else ""

    query_reg = hit_elem.find("query_reg")
    if query_reg is None or not query_reg.text:
        return None

    try:
        query_range = SequenceRange.parse(query_reg.text.strip())
    except ValueError as e:
        print(f"Invalid range in domain BLAST hit {domain_id}: {query_reg.text} - {e}")
        return None

    # Get hit range for coverage calculation
    hit_reg = hit_elem.find("hit_reg")
    hit_range = None
    if hit_reg is not None and hit_reg.text:
        try:
            hit_range = SequenceRange.parse(hit_reg.text.strip())
        except ValueError:
            pass

    # Get reference length and calculate coverage
    reference_length = ref_lengths.get(domain_id, 0)
    alignment_coverage = 0.0

    if reference_length > 0 and hit_range:
        alignment_coverage = hit_range.size / reference_length

    # Parse numeric values
    evalue = _safe_float(hit_elem.get("evalues", "999"), 999.0)
    identity = _safe_float(hit_elem.get("identity", "0"), 0.0)

    # Calculate confidence based on e-value, coverage, and identity
    confidence = _calculate_blast_confidence(evalue, alignment_coverage, identity)

    return Evidence(
        type="domain_blast",
        source_pdb=source_pdb,
        query_range=query_range,
        domain_id=domain_id,
        evalue=evalue,
        confidence=confidence,
        reference_length=reference_length,
        alignment_coverage=alignment_coverage,
        # Extract classification if available
        t_group=hit_elem.get("t_group"),
        h_group=hit_elem.get("h_group")
    )

def _parse_hhsearch_hit(hit_elem: ET.Element, ref_lengths: Dict[str, int]) -> Optional[Evidence]:
    """Parse an HHSearch hit with reference coverage"""
    hit_id = hit_elem.get("hit_id", "")
    domain_id = hit_elem.get("domain_id", "")

    # Prefer domain_id over hit_id
    effective_id = domain_id or hit_id
    if not effective_id:
        return None

    # Extract PDB from ID
    source_pdb = ""
    if len(effective_id) > 4 and effective_id[0] in 'ed':
        source_pdb = effective_id[1:5]

    query_reg = hit_elem.find("query_reg")
    if query_reg is None or not query_reg.text:
        return None

    try:
        query_range = SequenceRange.parse(query_reg.text.strip())
    except ValueError as e:
        print(f"Invalid range in HHSearch hit {effective_id}: {query_reg.text} - {e}")
        return None

    # Get hit range for coverage calculation
    hit_reg = hit_elem.find("hit_reg")
    hit_range = None
    if hit_reg is not None and hit_reg.text:
        try:
            hit_range = SequenceRange.parse(hit_reg.text.strip())
        except ValueError:
            pass

    # Get reference length and calculate coverage
    reference_length = ref_lengths.get(effective_id, 0)
    alignment_coverage = 0.0

    if reference_length > 0 and hit_range:
        alignment_coverage = hit_range.size / reference_length

    # Parse HHSearch-specific values
    probability = _safe_float(hit_elem.get("probability", "0"), 0.0)
    evalue = _safe_float(hit_elem.get("evalue", "999"), 999.0)
    score = _safe_float(hit_elem.get("score", "0"), 0.0)

    # Calculate confidence - HHSearch probability is most reliable
    confidence = _calculate_hhsearch_confidence(probability, evalue, alignment_coverage)

    return Evidence(
        type="hhsearch",
        source_pdb=source_pdb,
        query_range=query_range,
        domain_id=effective_id,
        evalue=evalue,
        confidence=confidence,
        reference_length=reference_length,
        alignment_coverage=alignment_coverage,
        # Extract classification if available
        t_group=hit_elem.get("t_group"),
        h_group=hit_elem.get("h_group")
    )

def _safe_float(value: str, default: float) -> float:
    """Safely parse float value with default"""
    try:
        return float(value)
    except (ValueError, TypeError):
        return default

def _evalue_to_confidence(evalue: float) -> float:
    """Convert e-value to confidence score (0-1)"""
    if evalue <= 0:
        return 1.0
    elif evalue < 1e-50:
        return 0.99
    elif evalue < 1e-20:
        return 0.95
    elif evalue < 1e-10:
        return 0.90
    elif evalue < 1e-5:
        return 0.80
    elif evalue < 0.001:
        return 0.70
    elif evalue < 0.01:
        return 0.60
    elif evalue < 0.1:
        return 0.50
    else:
        return 0.40

def _calculate_blast_confidence(evalue: float, coverage: float, identity: float) -> float:
    """
    Calculate BLAST confidence from multiple factors

    Args:
        evalue: BLAST E-value
        coverage: Fraction of reference domain covered
        identity: Percent identity

    Returns:
        Confidence score (0-1)
    """
    # Base confidence from e-value
    base_conf = _evalue_to_confidence(evalue)

    # Boost for good coverage
    coverage_boost = 0.0
    if coverage > 0.9:
        coverage_boost = 0.1
    elif coverage > 0.7:
        coverage_boost = 0.05
    elif coverage < 0.3 and coverage > 0:
        coverage_boost = -0.1  # Penalty for poor coverage

    # Boost for high identity
    identity_boost = 0.0
    if identity > 90:
        identity_boost = 0.05
    elif identity > 70:
        identity_boost = 0.02
    elif identity < 30:
        identity_boost = -0.05  # Penalty for low identity

    # Combine factors
    confidence = base_conf + coverage_boost + identity_boost

    # Clamp to valid range
    return max(0.0, min(1.0, confidence))

def _calculate_hhsearch_confidence(probability: float, evalue: float, coverage: float) -> float:
    """
    Calculate HHSearch confidence from probability and coverage

    Args:
        probability: HHSearch probability (0-100 or 0-1)
        evalue: HHSearch E-value
        coverage: Fraction of reference domain covered

    Returns:
        Confidence score (0-1)
    """
    # Normalize probability to 0-1 range
    if probability > 1.0:
        probability = probability / 100.0

    # HHSearch probability is already a good confidence measure
    base_conf = probability

    # Small adjustments for coverage
    if coverage > 0.8:
        base_conf = min(1.0, base_conf * 1.05)  # Small boost
    elif coverage < 0.3 and coverage > 0:
        base_conf = base_conf * 0.9  # Small penalty

    # Cross-check with e-value
    evalue_conf = _evalue_to_confidence(evalue)

    # Use probability as primary, but sanity check with e-value
    if base_conf > 0.8 and evalue_conf < 0.5:
        # Suspicious: high probability but bad e-value
        base_conf = (base_conf + evalue_conf) / 2

    return max(0.0, min(1.0, base_conf))
