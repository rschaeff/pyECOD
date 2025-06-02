"""Parse domain summary XML"""

import xml.etree.ElementTree as ET
from typing import List, Dict, Optional
from ecod.core.sequence_range import SequenceRange
from mini.models import Evidence

def parse_domain_summary(xml_path: str, verbose: bool = False) -> List[Evidence]:
    """
    Parse evidence from domain summary XML.

    Args:
        xml_path: Path to domain summary XML file
        verbose: Whether to print detailed parsing information

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

    evidence_list = []
    evidence_counts = {
        'chain_blast': 0,
        'domain_blast': 0,
        'hhsearch': 0
    }

    # Parse chain BLAST hits
    for hit in root.findall(".//chain_blast_run/hits/hit"):
        pdb_id = hit.get("pdb_id", "")
        query_reg = hit.find("query_reg")
        if query_reg is not None and query_reg.text:
            try:
                evidence = Evidence(
                    type="chain_blast",
                    source_pdb=pdb_id,
                    query_range=SequenceRange.parse(query_reg.text),
                    evalue=float(hit.get("evalues", "999"))
                )
                evidence_list.append(evidence)
                evidence_counts['chain_blast'] += 1
            except Exception as e:
                if verbose:
                    print(f"  Warning: Failed to parse chain BLAST hit {pdb_id}: {e}")

    # Parse domain BLAST hits with better metadata extraction
    for hit in root.findall(".//blast_run/hits/hit"):
        domain_id = hit.get("domain_id", "")
        # Extract PDB from domain ID (e.g., 'e6dgvA1' -> '6dgv')
        source_pdb = domain_id[1:5] if len(domain_id) > 4 else ""

        query_reg = hit.find("query_reg")
        hit_reg = hit.find("hit_reg")

        if query_reg is not None and query_reg.text:
            try:
                query_range = SequenceRange.parse(query_reg.text)

                # Calculate alignment coverage if we have hit range
                alignment_coverage = None
                if hit_reg is not None and hit_reg.text:
                    hit_range = SequenceRange.parse(hit_reg.text)
                    # Proxy calculation - would need DB for true reference length
                    alignment_coverage = hit_range.size / 300.0  # Assume avg domain = 300

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
                    alignment_coverage=alignment_coverage,
                    # Extract classification if available
                    t_group=hit.get("t_group"),
                    h_group=hit.get("h_group")
                )
                evidence_list.append(evidence)
                evidence_counts['domain_blast'] += 1
            except Exception as e:
                if verbose:
                    print(f"  Warning: Failed to parse domain BLAST hit {domain_id}: {e}")

    # Parse HHSearch hits
    for hit in root.findall(".//hh_run/hits/hit"):
        hit_id = hit.get("hit_id", "")
        domain_id = hit.get("domain_id", "")

        # Prefer domain_id for source
        source_pdb = ""
        if domain_id and len(domain_id) > 4:
            source_pdb = domain_id[1:5]
        elif hit_id and len(hit_id) > 4:
            source_pdb = hit_id[1:5]

        query_reg = hit.find("query_reg")
        if query_reg is not None and query_reg.text:
            try:
                prob = float(hit.get("probability", "0"))
                # Convert probability to confidence (0-100 scale to 0-1)
                confidence = prob / 100.0 if prob > 1.0 else prob

                evidence = Evidence(
                    type="hhsearch",
                    source_pdb=source_pdb,
                    query_range=SequenceRange.parse(query_reg.text),
                    confidence=confidence,
                    domain_id=domain_id or hit_id,
                    evalue=float(hit.get("evalue", "999"))
                )
                evidence_list.append(evidence)
                evidence_counts['hhsearch'] += 1
            except Exception as e:
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

    return evidence_list

def get_evidence_summary(evidence_list: List[Evidence]) -> Dict[str, any]:
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
