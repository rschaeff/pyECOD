"""Parse domain summary XML with robust chain_id extraction and comprehensive provenance tracking"""

import xml.etree.ElementTree as ET
from typing import List, Dict, Optional, Tuple, Any
from .sequence_range import SequenceRange
from .models import Evidence
from .evidence_utils import (
    populate_evidence_provenance,
    calculate_evidence_confidence,
    extract_chain_id_from_evidence,
    validate_evidence_provenance
)



def create_domain_id_lookup(domain_definitions: Dict = None) -> Dict[str, Tuple[str, str]]:
    """
    Create reverse lookup from domain_id -> (pdb_id, chain_id)
    This is much more robust than parsing domain_id strings

    Args:
        domain_definitions: Dict mapping (pdb_id, chain_id) -> List[DomainReference]

    Returns:
        Dict mapping domain_id -> (pdb_id, chain_id)
    """
    lookup = {}
    if domain_definitions:
        for (pdb_id, chain_id), domain_refs in domain_definitions.items():
            for domain_ref in domain_refs:
                lookup[domain_ref.domain_id] = (pdb_id, chain_id)
    return lookup


def extract_pdb_chain_robust(domain_id: str, domain_lookup: Dict[str, Tuple[str, str]]) -> Tuple[str, str]:
    """
    Robustly extract pdb_id and chain_id from domain_id using lookup first, fallback parsing second

    Args:
        domain_id: Domain identifier (e.g., "e6dgvA1", "d2ia4B2")
        domain_lookup: Lookup table from create_domain_id_lookup()

    Returns:
        Tuple of (pdb_id, chain_id)
    """
    # Primary: Use lookup table (most reliable)
    if domain_id in domain_lookup:
        return domain_lookup[domain_id]

    # Fallback: Parse domain_id (handles cases not in lookup)
    source_pdb = ""
    chain_id = ""

    if len(domain_id) > 4:
        # Common format: "e6dgvA1" -> pdb="6dgv", chain="A"
        source_pdb = domain_id[1:5]

        if len(domain_id) > 5:
            # Extract chain_id more carefully - look for first alphabetic char after PDB
            remainder = domain_id[5:]
            for char in remainder:
                if char.isalpha():
                    chain_id = char
                    break

    return source_pdb, chain_id


def parse_domain_summary(xml_path: str,
                        reference_lengths: Dict[str, int] = None,
                        protein_lengths: Dict[Tuple[str, str], int] = None,
                        blast_alignments: Dict[Tuple[str, str], Any] = None,
                        domain_definitions: Dict = None,
                        verbose: bool = False,
                        require_reference_lengths: bool = True) -> List[Evidence]:
    """
    Parse evidence from domain summary XML with robust chain_id extraction.

    Args:
        xml_path: Path to domain summary XML file
        reference_lengths: Optional dict of reference domain lengths
        protein_lengths: Optional dict of protein lengths by (pdb, chain)
        blast_alignments: Optional dict of BLAST alignment data
        domain_definitions: Optional domain definitions for robust chain_id lookup
        verbose: Whether to print detailed parsing information
        require_reference_lengths: If True, skip evidence without reference lengths

    Returns:
        List of Evidence objects
    """
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
    except ET.ParseError as e:
        print(f"ERROR: Failed to parse domain summary XML {xml_path}: {e}")
        return []
    except FileNotFoundError:
        print(f"ERROR: Domain summary file not found: {xml_path}")
        return []
    except Exception as e:
        print(f"ERROR: Unexpected error parsing {xml_path}: {e}")
        return []

    # Initialize default dicts if not provided
    if reference_lengths is None:
        reference_lengths = {}
    if protein_lengths is None:
        protein_lengths = {}
    if blast_alignments is None:
        blast_alignments = {}

    evidence_list = []
    evidence_counts = {
        'chain_blast': 0,
        'domain_blast': 0,
        'hhsearch': 0
    }
    skipped_counts = {
        'no_reference_length': 0,
        'parse_error': 0
    }

    # Create domain_id lookup for robust chain_id extraction
    domain_lookup = create_domain_id_lookup(domain_definitions)
    chain_id_stats = {
        'via_lookup': 0,
        'via_fallback': 0,
        'extraction_failed': 0
    }

    if verbose and domain_lookup:
        print(f"Created domain lookup for {len(domain_lookup)} domain IDs")

    skipped_chain_blast = []

    # Parse chain BLAST hits
    for hit in root.findall(".//chain_blast_run/hits/hit"):
        pdb_id = hit.get("pdb_id", "")
        chain_id = hit.get("chain_id", "")
        query_reg = hit.find("query_reg")
        if query_reg is not None and query_reg.text:
            try:
                evalue = float(hit.get("evalues", "999"))
                confidence = 0.5  # Default
                if evalue < 1e-10:
                    confidence = 0.9
                elif evalue < 1e-5:
                    confidence = 0.7
                elif evalue < 0.001:
                    confidence = 0.6

                # Get protein length for chain BLAST
                reference_length = None
                if protein_lengths:
                    pdb_lower = pdb_id.lower()
                    # Try multiple lookup patterns
                    lookup_keys = [
                        (pdb_lower, chain_id),
                        (pdb_id, chain_id),
                        f"{pdb_lower}_{chain_id}",
                        f"{pdb_id}_{chain_id}"
                    ]

                    for key in lookup_keys:
                        if key in protein_lengths:
                            reference_length = protein_lengths[key]
                            break

                if require_reference_lengths and reference_length is None:
                    skipped_chain_blast.append(f"{pdb_id}_{chain_id}")
                    skipped_counts['no_reference_length'] += 1
                    continue

                # Create domain_id for chain BLAST
                domain_id = f"{pdb_id}_{chain_id}"

                evidence = Evidence(
                    type="chain_blast",
                    source_pdb=pdb_id,
                    query_range=SequenceRange.parse(query_reg.text),
                    evalue=evalue,
                    domain_id=f"{pdb_id}_{chain_id}"
                )

                evidence = populate_evidence_provenance(
                    evidence=evidence,
                    alignment=blast_alignments.get((pdb_id, chain_id)),
                    reference_length=reference_length
                )

                is_valid, validation_issues = validate_evidence_provenance(evidence)
                if not is_valid:
                    skipped_counts['validation_failed'] += 1
                    if verbose:
                        print(f"  Warning: Evidence validation failed: {validation_issues}")
                    continue

                # Add alignment data if available
                if blast_alignments and (pdb_id, chain_id) in blast_alignments:
                    evidence.alignment = blast_alignments[(pdb_id, chain_id)]

                evidence_list.append(evidence)
                evidence_counts['chain_blast'] += 1

            except Exception as e:
                skipped_counts['parse_error'] += 1
                if verbose:
                    print(f"  Warning: Failed to parse chain BLAST hit {pdb_id}: {e}")

    # Parse domain BLAST hits with robust chain_id extraction
    for hit in root.findall(".//blast_run/hits/hit"):
        domain_id = hit.get("domain_id", "")
        if not domain_id:
            if verbose:
                print(f"  Warning: Skipping domain BLAST hit without domain_id")
            skipped_counts['parse_error'] += 1
            continue

        # ROBUST: Use lookup-based approach for pdb_id and chain_id extraction
        source_pdb, chain_id = extract_pdb_chain_robust(domain_id, domain_lookup)

        # Track extraction method for statistics
        if domain_id in domain_lookup:
            chain_id_stats['via_lookup'] += 1
        elif chain_id:
            chain_id_stats['via_fallback'] += 1
        else:
            chain_id_stats['extraction_failed'] += 1
            if verbose:
                print(f"  Warning: Could not determine chain_id for domain {domain_id}")

        query_reg = hit.find("query_reg")
        hit_reg = hit.find("hit_reg")

        if query_reg is not None and query_reg.text:
            try:
                query_range = SequenceRange.parse(query_reg.text)

                # Look up reference length
                reference_length = None
                if reference_lengths:
                    # Try exact domain_id match first
                    if domain_id in reference_lengths:
                        reference_length = reference_lengths[domain_id]
                    # Try source_pdb match
                    elif source_pdb in reference_lengths:
                        reference_length = reference_lengths[source_pdb]
                    # Try without the 'e' prefix if domain_id starts with 'e'
                    elif domain_id.startswith('e') and domain_id[1:] in reference_lengths:
                        reference_length = reference_lengths[domain_id[1:]]

                # Calculate alignment coverage if we have hit range and reference data
                alignment_coverage = None
                hit_range = None
                if hit_reg is not None and hit_reg.text and reference_length:
                    hit_range = SequenceRange.parse(hit_reg.text)
                    alignment_coverage = hit_range.total_length / reference_length

                # Skip if reference lengths are required but missing
                if require_reference_lengths and not reference_length:
                    if verbose:
                        print(f"  Skipping {domain_id}: no reference length available")
                    skipped_counts['no_reference_length'] += 1
                    continue

                # Parse e-value and calculate confidence
                evalue = float(hit.get("evalues", "999"))
                confidence = 0.5  # Default
                if evalue < 1e-10:
                    confidence = 0.9
                elif evalue < 1e-5:
                    confidence = 0.7

                # Create evidence with robust provenance fields
                evidence = Evidence(
                    type="domain_blast",
                    source_pdb=source_pdb,
                    query_range=query_range,
                    domain_id=domain_id,
                    evalue=evalue
                )

                classification = {
                    't_group': hit.get("t_group"),
                    'h_group': hit.get("h_group")
                }

                evidence = populate_evidence_provenance(
                    evidence=evidence,
                    hit_range=hit_range,
                    reference_length=reference_length,
                    classification=classification
                )

                is_valid, validation_issues = validate_evidence_provenance(evidence)
                if not is_valid:
                    skipped_counts['validation_failed'] += 1
                    if verbose:
                        print(f"  Warning: Evidence validation failed: {validation_issues}")
                    continue

                evidence_list.append(evidence)
                evidence_counts['domain_blast'] += 1

            except Exception as e:
                skipped_counts['parse_error'] += 1
                if verbose:
                    print(f"  Warning: Failed to parse domain BLAST hit {domain_id}: {e}")

    # Parse HHSearch hits
    for hit in root.findall(".//hh_run/hits/hit"):
        hit_id = hit.get("hit_id", "")
        domain_id = hit.get("domain_id", "") or hit_id

        if not domain_id:
            if verbose:
                print(f"  Warning: Skipping HHSearch hit without domain_id or hit_id")
            skipped_counts['parse_error'] += 1
            continue

        # Use robust extraction for HHsearch as well
        source_pdb, chain_id = extract_pdb_chain_robust(domain_id, domain_lookup)

        # Fallback for HHsearch hit_id if domain_id extraction failed
        if not source_pdb and hit_id and len(hit_id) > 4:
            source_pdb = hit_id[1:5]

        query_reg = hit.find("query_reg")
        if query_reg is not None and query_reg.text:
            try:
                prob = float(hit.get("probability", "0"))
                confidence = prob / 100.0 if prob > 1.0 else prob

                # Look for reference length
                reference_length = None
                if reference_lengths:
                    lookup_id = domain_id or hit_id

                    # Try exact domain_id/hit_id match first
                    if lookup_id in reference_lengths:
                        reference_length = reference_lengths[lookup_id]
                    # Try source_pdb match
                    elif source_pdb in reference_lengths:
                        reference_length = reference_lengths[source_pdb]
                    # Try without the 'e' prefix if domain_id starts with 'e'
                    elif lookup_id.startswith('e') and lookup_id[1:] in reference_lengths:
                        reference_length = reference_lengths[lookup_id[1:]]

                # Skip if reference lengths are required but missing
                if require_reference_lengths and not reference_length:
                    if verbose:
                        print(f"  Skipping HHSearch {lookup_id}: no reference length available")
                    skipped_counts['no_reference_length'] += 1
                    continue

                evidence = Evidence(
                    type="hhsearch",
                    source_pdb=source_pdb,
                    query_range=SequenceRange.parse(query_reg.text),
                    domain_id=domain_id,
                    evalue=float(hit.get("evalue", "999"))
                )

                evidence = populate_evidence_provenance(
                    evidence=evidence,
                    reference_length=reference_length
                )

                prob = float(hit.get("probability", "0"))
                evidence.confidence = calculate_evidence_confidence(
                    probability=prob,
                    evidence_type="hhsearch",
                    alignment_coverage=evidence.alignment_coverage
                )

                is_valid, validation_issues = validate_evidence_provenance(evidence)
                if not is_valid:
                    skipped_counts['validation_failed'] += 1
                    if verbose:
                        print(f"  Warning: Evidence validation failed: {validation_issues}")
                    continue

                evidence_list.append(evidence)
                evidence_counts['hhsearch'] += 1

            except Exception as e:
                skipped_counts['parse_error'] += 1
                if verbose:
                    print(f"  Warning: Failed to parse HHSearch hit {hit_id}: {e}")

    # Print summary if verbose or if no evidence found
    if verbose or len(evidence_list) == 0:
        if len(evidence_list) == 0:
            print(f"WARNING: No evidence found in {xml_path}")
        elif verbose:
            print(f"Parsed {len(evidence_list)} evidence items from {xml_path}:")
            for etype, count in evidence_counts.items():
                if count > 0:
                    print(f"  {etype}: {count}")

            # Report chain_id extraction statistics
            if any(chain_id_stats.values()):
                print(f"Chain ID extraction:")
                if chain_id_stats['via_lookup'] > 0:
                    print(f"  Via lookup: {chain_id_stats['via_lookup']}")
                if chain_id_stats['via_fallback'] > 0:
                    print(f"  Via fallback: {chain_id_stats['via_fallback']}")
                if chain_id_stats['extraction_failed'] > 0:
                    print(f"  Failed: {chain_id_stats['extraction_failed']}")

            # Report skipped items
            total_skipped = sum(skipped_counts.values())
            if total_skipped > 0:
                print(f"Skipped {total_skipped} items:")
                for reason, count in skipped_counts.items():
                    if count > 0:
                        print(f"  {reason.replace('_', ' ')}: {count}")

                # Show details for chain BLAST skips
                if skipped_chain_blast and verbose:
                    print(f"  Skipped chain BLAST entries (no protein length): {len(skipped_chain_blast)} total")
                    if len(skipped_chain_blast) <= 10:
                        print(f"    {', '.join(skipped_chain_blast)}")
                    else:
                        print(f"    {', '.join(skipped_chain_blast[:10])} ... and {len(skipped_chain_blast)-10} more")

    return evidence_list


def get_evidence_summary(evidence_list: List[Evidence]) -> Dict[str, Any]:
    """
    Get summary statistics for evidence.

    Args:
        evidence_list: List of Evidence objects

    Returns:
        Dictionary with summary statistics
    """
    if not evidence_list:
        return {
            'total': 0,
            'by_type': {},
            'high_confidence': 0,
            'unique_families': 0,
            'with_provenance': 0
        }

    by_type = {}
    high_conf = 0
    families = set()
    with_provenance = 0

    for ev in evidence_list:
        # Count by type
        by_type[ev.type] = by_type.get(ev.type, 0) + 1

        # Count high confidence
        if ev.confidence > 0.7 or (ev.evalue and ev.evalue < 1e-10):
            high_conf += 1

        # Track unique families
        if ev.source_pdb:
            families.add(ev.source_pdb)

        # Count evidence with provenance data
        if ev.source_chain_id or ev.hit_range or ev.reference_length:
            with_provenance += 1

    return {
        'total': len(evidence_list),
        'by_type': by_type,
        'high_confidence': high_conf,
        'unique_families': len(families),
        'with_provenance': with_provenance,
        'provenance_percent': (with_provenance / len(evidence_list)) * 100 if evidence_list else 0
    }


def load_reference_lengths(reference_file: str = None) -> Dict[str, int]:
    """
    Load reference protein lengths from file.

    Args:
        reference_file: Path to reference lengths file (CSV or similar)

    Returns:
        Dictionary mapping protein_id -> length
    """
    reference_lengths = {}

    if reference_file is None:
        return reference_lengths

    try:
        import csv
        with open(reference_file, 'r') as f:
            reader = csv.reader(f)
            # Skip header if present
            first_row = next(reader, None)
            if first_row and not first_row[1].isdigit():
                # This was a header, continue
                pass
            else:
                # This was data, process it
                if first_row and len(first_row) >= 2:
                    reference_lengths[first_row[0]] = int(first_row[1])

            # Process remaining rows
            for row in reader:
                if len(row) >= 2 and row[1].isdigit():
                    reference_lengths[row[0]] = int(row[1])

    except FileNotFoundError:
        print(f"WARNING: Reference lengths file not found: {reference_file}")
    except Exception as e:
        print(f"WARNING: Error loading reference lengths: {e}")

    return reference_lengths


def load_protein_lengths(protein_file: str = None) -> Dict[Tuple[str, str], int]:
    """
    Load protein lengths from file.

    Args:
        protein_file: Path to protein lengths file

    Returns:
        Dictionary mapping (pdb_id, chain_id) -> length
    """
    protein_lengths = {}

    if protein_file is None:
        return protein_lengths

    try:
        import csv
        with open(protein_file, 'r') as f:
            reader = csv.reader(f)
            # Skip header if present
            first_row = next(reader, None)
            if first_row and not first_row[-1].isdigit():
                # This was a header, continue
                pass
            else:
                # This was data, process it
                if first_row:
                    if len(first_row) >= 3:
                        # Format: pdb_id,chain_id,length
                        protein_lengths[(first_row[0], first_row[1])] = int(first_row[2])
                    elif len(first_row) >= 2 and '_' in first_row[0]:
                        # Format: pdb_chain,length
                        parts = first_row[0].split('_')
                        if len(parts) >= 2:
                            protein_lengths[(parts[0], parts[1])] = int(first_row[1])

            # Process remaining rows
            for row in reader:
                if len(row) >= 3 and row[2].isdigit():
                    # Format: pdb_id,chain_id,length
                    protein_lengths[(row[0], row[1])] = int(row[2])
                elif len(row) >= 2 and '_' in row[0] and row[1].isdigit():
                    # Format: pdb_chain,length
                    parts = row[0].split('_')
                    if len(parts) >= 2:
                        protein_lengths[(parts[0], parts[1])] = int(row[1])

    except FileNotFoundError:
        print(f"WARNING: Protein lengths file not found: {protein_file}")
    except Exception as e:
        print(f"WARNING: Error loading protein lengths: {e}")

    return protein_lengths


def validate_evidence_consistency(evidence_list: List[Evidence]) -> Dict[str, Any]:
    """
    Validate evidence for consistency issues and potential problems.

    Args:
        evidence_list: List of Evidence objects to validate

    Returns:
        Dictionary with validation results
    """
    issues = []
    stats = {
        'total_evidence': len(evidence_list),
        'missing_chain_id': 0,
        'missing_reference_length': 0,
        'low_confidence': 0,
        'invalid_ranges': 0,
        'duplicate_domain_ids': 0
    }

    seen_domain_ids = set()

    for ev in evidence_list:
        # Check for missing chain_id
        if not ev.source_chain_id:
            stats['missing_chain_id'] += 1

        # Check for missing reference length
        if not ev.reference_length:
            stats['missing_reference_length'] += 1

        # Check for low confidence
        if ev.confidence < 0.5:
            stats['low_confidence'] += 1

        # Check for invalid ranges
        try:
            positions = ev.query_range.to_positions_simple()
            if len(positions) == 0:
                stats['invalid_ranges'] += 1
                issues.append(f"Empty range for {ev.domain_id}: {ev.query_range}")
        except Exception as e:
            stats['invalid_ranges'] += 1
            issues.append(f"Invalid range for {ev.domain_id}: {e}")

        # Check for duplicate domain_ids
        if ev.domain_id and ev.domain_id in seen_domain_ids:
            stats['duplicate_domain_ids'] += 1
            issues.append(f"Duplicate domain_id: {ev.domain_id}")
        elif ev.domain_id:
            seen_domain_ids.add(ev.domain_id)

    return {
        'stats': stats,
        'issues': issues,
        'is_valid': len(issues) == 0
    }
