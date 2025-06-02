# mini/parser.py
"""Parse domain summary XML"""

import xml.etree.ElementTree as ET
from typing import List
from .models import Evidence
from ecod.core.sequence_range import SequenceRange

def parse_domain_summary(xml_path: str) -> List[Evidence]:
    """Parse evidence from domain summary XML"""
    tree = ET.parse(xml_path)
    root = tree.getroot()

    evidence_list = []

    # Parse chain BLAST hits
    for hit in root.findall(".//chain_blast_run/hits/hit"):
        pdb_id = hit.get("pdb_id", "")
        query_reg = hit.find("query_reg")
        if query_reg is not None and query_reg.text:
            evidence = Evidence(
                type="chain_blast",
                source_pdb=pdb_id,
                query_range=SequenceRange.parse(query_reg.text),
                evalue=float(hit.get("evalues", "999"))
            )
            evidence_list.append(evidence)

    # Parse domain BLAST hits with classification
    for hit in root.findall(".//blast_run/hits/hit"):
        domain_id = hit.get("domain_id", "")
        # Extract PDB from domain ID (e.g., 'e6dgvA1' -> '6dgv')
        source_pdb = domain_id[1:5] if len(domain_id) > 4 else ""

        query_reg = hit.find("query_reg")
        if query_reg is not None and query_reg.text:
            evidence = Evidence(
                type="domain_blast",
                source_pdb=source_pdb,
                query_range=SequenceRange.parse(query_reg.text),
                domain_id=domain_id,
                evalue=float(hit.get("evalues", "999")),
                # Extract classification if available
                t_group=hit.get("t_group"),
                h_group=hit.get("h_group")
            )
            evidence_list.append(evidence)

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
    
    return evidence_list
