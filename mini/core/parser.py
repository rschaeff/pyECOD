"""Parse domain summary XML"""

import xml.etree.ElementTree as ET
from typing import List, Dict, Optional, Tuple, Any
from .sequence_range import SequenceRange
from .models import Evidence

def parse_domain_summary(xml_path: str,
                        reference_lengths: Dict[str, int] = None,
                        protein_lengths: Dict[Tuple[str, str], int] = None,
                        blast_alignments: Dict[Tuple[str, str], Any] = None,
                        verbose: bool = False,
                        require_reference_lengths: bool = True) -> List[Evidence]:
    """
    Parse evidence from domain summary XML.

    Args:
        xml_path: Path to domain summary XML file
        reference_lengths: Optional dict of reference domain lengths
        protein_lengths: Optional dict of protein lengths by (pdb, chain)
        blast_alignments: Optional dict of BLAST alignment data
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

    skipped_chain_blast = []

    # Parse chain BLAST hits - FIX: Always assign domain_id
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
                    if (pdb_lower, chain_id) in protein_lengths:
                        reference_length = protein_lengths[(pdb_lower, chain_id)]
                    elif f"{pdb_lower}_{chain_id}" in protein_lengths:
                        reference_length = protein_lengths[f"{pdb_lower}_{chain_id}"]
                    elif (pdb_id, chain_id) in protein_lengths:
                        reference_length = protein_lengths[(pdb_id, chain_id)]
                    elif f"{pdb_id}_{chain_id}" in protein_lengths:
                        reference_length = protein_lengths[f"{pdb_id}_{chain_id}"]

                if require_reference_lengths and reference_length is None:
                    skipped_chain_blast.append(f"{pdb_id}_{chain_id}")
                    skipped_counts['no_reference_length'] += 1
                    continue

                # FIXED: Always assign domain_id for chain BLAST
                domain_id = f"{pdb_id}_{chain_id}"  # Consistent format

                evidence = Evidence(
                    type="chain_blast",
                    source_pdb=pdb_id,
                    query_range=SequenceRange.parse(query_reg.text),
                    evalue=evalue,
                    confidence=confidence,
                    domain_id=domain_id,  # ALWAYS set this
                    reference_length=reference_length,
                )

                # Add alignment data if available
                if blast_alignments and (pdb_id, chain_id) in blast_alignments:
                    evidence.alignment = blast_alignments[(pdb_id, chain_id)]

                evidence_list.append(evidence)
                evidence_counts['chain_blast'] += 1
            except Exception as e:
                skipped_counts['parse_error'] += 1
                if verbose:
                    print(f"  Warning: Failed to parse chain BLAST hit {pdb_id}: {e}")

    # Parse domain BLAST hits - domain_id should already be present
    for hit in root.findall(".//blast_run/hits/hit"):
        domain_id = hit.get("domain_id", "")
        if not domain_id:
            # FIXED: Skip hits without domain_id rather than proceeding
            if verbose:
                print(f"  Warning: Skipping domain BLAST hit without domain_id")
            skipped_counts['parse_error'] += 1
            continue

        source_pdb = domain_id[1:5] if len(domain_id) > 4 else ""

        query_reg = hit.find("query_reg")
        hit_reg = hit.find("hit_reg")

        if query_reg is not None and query_reg.text:
            try:
                query_range = SequenceRange.parse(query_reg.text)

                # FIXED: Look up reference length regardless of hit_reg presence
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
                if hit_reg is not None and hit_reg.text and reference_length:
                    hit_range = SequenceRange.parse(hit_reg.text)
                    alignment_coverage = hit_range.total_length / reference_length
                elif verbose and not reference_length:
                    print(f"  Warning: No reference length for {domain_id}, skipping coverage calculation")

                # Skip if reference lengths are required but missing
                if require_reference_lengths and not reference_length:
                    if verbose:
                        print(f"  Skipping {domain_id}: no reference length available")
                    skipped_counts['no_reference_length'] += 1
                    continue

                # Parse e-value
                evalue = float(hit.get("evalues", "999"))

                # Calculate confidence based on e-value and coverage
                confidence = 0.5  # Default
                if evalue < 1e-10:
                    confidence = 0.9
                elif evalue < 1e-5:
                    confidence = 0.7

                evidence = Evidence(
                    type="domain_blast",
                    source_pdb=source_pdb,
                    query_range=query_range,
                    domain_id=domain_id,
                    evalue=evalue,
                    confidence=confidence,

                    # NEW PROVENANCE FIELDS:
                    source_chain_id=chain_id,  # Extract from parsing
                    hit_range=hit_range if hit_reg else None,
                    hsp_count=int(hit.get("hsp_count", "1")),
                    discontinuous=query_range.is_discontinuous,
                    reference_length=reference_length,
                    alignment_coverage=alignment_coverage,
                    t_group=hit.get("t_group"),
                    h_group=hit.get("h_group")
                )
                evidence_list.append(evidence)
                evidence_counts['domain_blast'] += 1
            except Exception as e:
                skipped_counts['parse_error'] += 1
                if verbose:
                    print(f"  Warning: Failed to parse domain BLAST hit {domain_id}: {e}")


    # Parse HHSearch hits
    for hit in root.findall(".//hh_run/hits/hit"):
        hit_id = hit.get("hit_id", "")
        domain_id = hit.get("domain_id", "") or hit_id  # Use hit_id as fallback

        if not domain_id:
            # FIXED: Skip hits without any identifier
            if verbose:
                print(f"  Warning: Skipping HHSearch hit without domain_id or hit_id")
            skipped_counts['parse_error'] += 1
            continue

        source_pdb = ""
        if domain_id and len(domain_id) > 4:
            source_pdb = domain_id[1:5]
        elif hit_id and len(hit_id) > 4:
            source_pdb = hit_id[1:5]

        query_reg = hit.find("query_reg")
        if query_reg is not None and query_reg.text:
            try:
                prob = float(hit.get("probability", "0"))
                confidence = prob / 100.0 if prob > 1.0 else prob

                # Look for reference length - try multiple formats (same as domain_blast)
                reference_length = None
                if reference_lengths:
                    # Use the domain_id or hit_id for lookup
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
                    confidence=confidence,
                    domain_id=domain_id,  # ALWAYS set this
                    evalue=float(hit.get("evalue", "999")),
                    reference_length=reference_length
                )
                evidence_list.append(evidence)
                evidence_counts['hhsearch'] += 1
            except Exception as e:
                skipped_counts['parse_error'] += 1
                if verbose:
                    print(f"  Warning: Failed to parse HHSearch hit {hit_id}: {e}")

    # Only print summary if verbose or if no evidence found
    if verbose or len(evidence_list) == 0:
        if len(evidence_list) == 0:
            print(f"WARNING: No evidence found in {xml_path}")
        elif verbose:
            print(f"Parsed {len(evidence_list)} evidence items from {xml_path}:")
            for etype, count in evidence_counts.items():
                if count > 0:
                    print(f"  {etype}: {count}")

            # Report skipped items as summary
            total_skipped = sum(skipped_counts.values())
            if total_skipped > 0:
                print(f"Skipped {total_skipped} items:")
                for reason, count in skipped_counts.items():
                    if count > 0:
                        print(f"  {reason.replace('_', ' ')}: {count}")

                # Show details only if verbose AND there are chain BLAST skips
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
            'unique_families': 0
        }

    by_type = {}
    high_conf = 0
    families = set()

    for ev in evidence_list:
        # Count by type
        by_type[ev.type] = by_type.get(ev.type, 0) + 1

        # Count high confidence
        if ev.confidence > 0.7 or (ev.evalue and ev.evalue < 1e-10):
            high_conf += 1

        # Track unique families
        if ev.source_pdb:
            families.add(ev.source_pdb)

    return {
        'total': len(evidence_list),
        'by_type': by_type,
        'high_confidence': high_conf,
        'unique_families': len(families)
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
        # Return empty dict if no file specified
        return reference_lengths

    try:
        import csv
        with open(reference_file, 'r') as f:
            # Assume CSV format: protein_id,length
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
        # Return empty dict if no file specified
        return protein_lengths

    try:
        import csv
        with open(protein_file, 'r') as f:
            # Assume CSV format: pdb_id,chain_id,length or pdb_chain,length
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
