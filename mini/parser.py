# mini/parser.py - Extract more info
import xml.etree.ElementTree as ET
from typing import List
from ecod.core.sequence_range import SequenceRange

def parse_domain_summary(xml_path: str) -> List[Evidence]:
    """Parse evidence with reference coverage info"""
    tree = ET.parse(xml_path)
    root = tree.getroot()

    evidence_list = []

    # Parse domain BLAST hits
    for hit in root.findall(".//blast_run/hits/hit"):
        domain_id = hit.get("domain_id", "")
        source_pdb = domain_id[1:5] if len(domain_id) > 4 else ""

        query_reg = hit.find("query_reg")
        hit_reg = hit.find("hit_reg")

        if query_reg is not None and query_reg.text:
            query_range = SequenceRange.parse(query_reg.text)

            # Calculate alignment coverage if we have hit range
            alignment_coverage = None
            if hit_reg is not None and hit_reg.text:
                hit_range = SequenceRange.parse(hit_reg.text)
                # This is a proxy - we'd need DB for true reference length
                alignment_coverage = hit_range.size / 300.0  # Assume avg domain = 300

            evidence = Evidence(
                type="domain_blast",
                source_pdb=source_pdb,
                query_range=query_range,
                domain_id=domain_id,
                evalue=float(hit.get("evalues", "999")),
                alignment_coverage=alignment_coverage
            )

            # Calculate confidence including coverage
            if evidence.evalue < 1e-10 and alignment_coverage and alignment_coverage > 0.7:
                evidence.confidence = 0.9
            elif evidence.evalue < 1e-5:
                evidence.confidence = 0.7
            else:
                evidence.confidence = 0.5

            evidence_list.append(evidence)
    
    return evidence_list
