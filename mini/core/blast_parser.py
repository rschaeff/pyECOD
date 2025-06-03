# mini/blast_parser.py
"""Parse alignment data from raw BLAST XML files"""

import xml.etree.ElementTree as ET
from typing import Dict, Tuple, Optional, Any
from dataclasses import dataclass

@dataclass
class BlastAlignment:
    """BLAST alignment data"""
    query_seq: str
    hit_seq: str
    query_start: int
    query_end: int
    hit_start: int
    hit_end: int
    hit_id: str
    evalue: float

def parse_blast_xml(blast_xml_path: str, verbose: bool = False) -> Dict[Tuple[str, str], BlastAlignment]:
    """
    Parse raw BLAST XML to extract alignment data.

    Args:
        blast_xml_path: Path to BLAST XML file
        verbose: Whether to print detailed parsing information

    Returns:
        Dict mapping (pdb_id, chain_id) -> BlastAlignment
    """
    alignments = {}

    try:
        tree = ET.parse(blast_xml_path)
        root = tree.getroot()
    except ET.ParseError as e:
        # Only show serious parsing errors
        print(f"ERROR: Failed to parse BLAST XML {blast_xml_path}: {e}")
        return alignments
    except FileNotFoundError:
        # Silent fail for missing files (handled by caller)
        return alignments
    except Exception as e:
        # Unexpected errors should be reported
        print(f"ERROR: Unexpected error parsing {blast_xml_path}: {e}")
        return alignments

    # Navigate through BLAST XML structure
    iterations_found = 0
    hits_found = 0
    alignments_extracted = 0

    for iteration in root.findall(".//Iteration"):
        iterations_found += 1
        for hit in iteration.findall(".//Hit"):
            hits_found += 1

            # Extract hit info
            hit_def_elem = hit.find("Hit_def")
            if hit_def_elem is None:
                continue

            hit_def = hit_def_elem.text or ""

            # Parse PDB and chain from hit definition (e.g., "6dgv A")
            parts = hit_def.split()
            if len(parts) >= 2:
                pdb_id = parts[0].lower()
                chain_id = parts[1]

                # Get the first HSP (High-scoring Segment Pair)
                hsp = hit.find(".//Hsp")
                if hsp is not None:
                    try:
                        # Extract alignment data
                        query_seq = hsp.find("Hsp_qseq").text if hsp.find("Hsp_qseq") is not None else ""
                        hit_seq = hsp.find("Hsp_hseq").text if hsp.find("Hsp_hseq") is not None else ""

                        query_start = int(hsp.find("Hsp_query-from").text) if hsp.find("Hsp_query-from") is not None else 1
                        query_end = int(hsp.find("Hsp_query-to").text) if hsp.find("Hsp_query-to") is not None else len(query_seq)
                        hit_start = int(hsp.find("Hsp_hit-from").text) if hsp.find("Hsp_hit-from") is not None else 1
                        hit_end = int(hsp.find("Hsp_hit-to").text) if hsp.find("Hsp_hit-to") is not None else len(hit_seq)

                        evalue = float(hsp.find("Hsp_evalue").text) if hsp.find("Hsp_evalue") is not None else 999.0

                        alignment = BlastAlignment(
                            query_seq=query_seq,
                            hit_seq=hit_seq,
                            query_start=query_start,
                            query_end=query_end,
                            hit_start=hit_start,
                            hit_end=hit_end,
                            hit_id=hit_def,
                            evalue=evalue
                        )

                        alignments[(pdb_id, chain_id)] = alignment
                        alignments_extracted += 1

                    except (ValueError, AttributeError) as e:
                        if verbose:
                            print(f"  Warning: Failed to parse HSP for {hit_def}: {e}")
                        continue

    # Only print summary if alignments were found or in verbose mode
    if alignments_extracted > 0 or verbose:
        if verbose:
            print(f"BLAST parsing summary for {blast_xml_path}:")
            print(f"  Iterations: {iterations_found}")
            print(f"  Hits: {hits_found}")
            print(f"  Alignments extracted: {alignments_extracted}")
        # For non-verbose, no output unless there's an issue

    return alignments

def load_chain_blast_alignments(blast_dir: str, pdb_id: str, chain_id: str, verbose: bool = False) -> Dict[Tuple[str, str], BlastAlignment]:
    """
    Load chain BLAST alignments for a specific query.

    Args:
        blast_dir: Directory containing BLAST XML files
        pdb_id: Query PDB ID
        chain_id: Query chain ID
        verbose: Whether to print detailed information

    Returns:
        Dict mapping (hit_pdb, hit_chain) -> BlastAlignment
    """
    import os

    # Construct expected filename
    blast_file = os.path.join(blast_dir, f"{pdb_id}_{chain_id}.develop291.xml")

    if not os.path.exists(blast_file):
        # Try alternate naming conventions
        blast_file_upper = os.path.join(blast_dir, f"{pdb_id.upper()}_{chain_id}.develop291.xml")
        if os.path.exists(blast_file_upper):
            blast_file = blast_file_upper
        else:
            if verbose:
                print(f"BLAST file not found: {blast_file}")
            return {}

    return parse_blast_xml(blast_file, verbose=verbose)

def get_blast_summary(alignments: Dict[Tuple[str, str], BlastAlignment]) -> Dict[str, Any]:
    """
    Get summary statistics for BLAST alignments.

    Args:
        alignments: Dictionary of BLAST alignments

    Returns:
        Dictionary with summary statistics
    """
    if not alignments:
        return {
            'count': 0,
            'avg_evalue': None,
            'best_evalue': None,
            'avg_coverage': None
        }

    evalues = [a.evalue for a in alignments.values() if a.evalue < 100]
    coverages = []

    for alignment in alignments.values():
        query_len = alignment.query_end - alignment.query_start + 1
        hit_len = alignment.hit_end - alignment.hit_start + 1
        coverage = min(query_len, hit_len) / max(query_len, hit_len) if max(query_len, hit_len) > 0 else 0
        coverages.append(coverage)

    return {
        'count': len(alignments),
        'avg_evalue': sum(evalues) / len(evalues) if evalues else None,
        'best_evalue': min(evalues) if evalues else None,
        'avg_coverage': sum(coverages) / len(coverages) if coverages else None
    }
