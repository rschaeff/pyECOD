# mini/blast_parser.py
"""Parse alignment data from raw BLAST XML files"""

import xml.etree.ElementTree as ET
from typing import Dict, Tuple, Optional
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

def parse_blast_xml(blast_xml_path: str) -> Dict[Tuple[str, str], BlastAlignment]:
    """
    Parse raw BLAST XML to extract alignment data.
    
    Args:
        blast_xml_path: Path to BLAST XML file
        
    Returns:
        Dict mapping (pdb_id, chain_id) -> BlastAlignment
    """
    alignments = {}
    
    try:
        tree = ET.parse(blast_xml_path)
        root = tree.getroot()
    except Exception as e:
        print(f"Error parsing BLAST XML {blast_xml_path}: {e}")
        return alignments
    
    # Navigate through BLAST XML structure
    for iteration in root.findall(".//Iteration"):
        for hit in iteration.findall(".//Hit"):
            # Extract hit info
            hit_def = hit.find("Hit_def").text if hit.find("Hit_def") is not None else ""
            
            # Parse PDB and chain from hit definition (e.g., "6dgv A")
            parts = hit_def.split()
            if len(parts) >= 2:
                pdb_id = parts[0].lower()
                chain_id = parts[1]
                
                # Get the first HSP (High-scoring Segment Pair)
                hsp = hit.find(".//Hsp")
                if hsp is not None:
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
                    
    print(f"Parsed {len(alignments)} alignments from {blast_xml_path}")
    return alignments

def load_chain_blast_alignments(blast_dir: str, pdb_id: str, chain_id: str) -> Dict[Tuple[str, str], BlastAlignment]:
    """
    Load chain BLAST alignments for a specific query.
    
    Args:
        blast_dir: Directory containing BLAST XML files
        pdb_id: Query PDB ID
        chain_id: Query chain ID
        
    Returns:
        Dict mapping (hit_pdb, hit_chain) -> BlastAlignment
    """
    import os
    
    # Construct expected filename
    blast_file = os.path.join(blast_dir, f"{pdb_id}_{chain_id}.develop291.xml")
    
    if not os.path.exists(blast_file):
        print(f"BLAST file not found: {blast_file}")
        # Try alternate naming conventions
        blast_file = os.path.join(blast_dir, f"{pdb_id.upper()}_{chain_id}.develop291.xml")
        if not os.path.exists(blast_file):
            return {}
    
    return parse_blast_xml(blast_file)
