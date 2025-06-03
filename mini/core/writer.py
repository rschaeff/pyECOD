# mini_pyecod/writer.py
"""Write domain partition results"""

import xml.etree.ElementTree as ET
from typing import List
from .models import Domain

def write_domain_partition(domains: List[Domain], pdb_id: str, chain_id: str, 
                          output_path: str, reference: str = "mini_pyecod"):
    """Write domains to XML file"""
    root = ET.Element("domain_partition")
    root.set("pdb_id", pdb_id)
    root.set("chain_id", chain_id)
    root.set("reference", reference)
    root.set("is_classified", "true" if domains else "false")
    
    domains_elem = ET.SubElement(root, "domains")
    
    for domain in domains:
        d_elem = ET.SubElement(domains_elem, "domain")
        d_elem.set("id", domain.id)
        d_elem.set("range", str(domain.range))
        d_elem.set("family", domain.family)
        d_elem.set("source", domain.source)
        d_elem.set("evidence_count", str(domain.evidence_count))
        d_elem.set("is_discontinuous", str(domain.range.is_discontinuous).lower())
    
    # Pretty print
    tree = ET.ElementTree(root)
    ET.indent(tree, space="  ")
    tree.write(output_path, encoding="utf-8", xml_declaration=True)
