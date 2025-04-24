# ecod/utils/xml_conversion.py
import xml.etree.ElementTree as ET
from typing import Dict, Any, List, Optional

from ecod.models.pipeline import (
    BlastHit, HHSearchHit, Evidence, DomainCandidate,
    DomainResult, DomainSummary, DomainPartitionResult
)

def parse_domain_summary_xml(file_path: str) -> DomainSummary:
    """Parse domain summary XML into model"""
    try:
        tree = ET.parse(file_path)
        root = tree.getroot()
        
        # Create base summary model
        blast_summ = root.find(".//blast_summ")
        if blast_summ is None:
            raise ValueError("No blast_summ element found in XML")
            
        summary = DomainSummary(
            pdb_id=blast_summ.get("pdb", ""),
            chain_id=blast_summ.get("chain", ""),
            reference=blast_summ.get("reference", ""),
            is_peptide=blast_summ.get("is_peptide", "false").lower() == "true"
        )
        
        # Get sequence length
        query_len_elem = root.find(".//query_len")
        if query_len_elem is not None and query_len_elem.text:
            summary.sequence_length = int(query_len_elem.text.strip())
        
        # Process chain BLAST hits
        for hit_elem in root.findall(".//chain_blast_run/hits/hit"):
            hit = BlastHit(
                hit_id=hit_elem.get("num", ""),
                domain_id="",  # Chain hits don't have domain IDs
                pdb_id=hit_elem.get("pdb_id", ""),
                chain_id=hit_elem.get("chain_id", ""),
                evalue=_parse_float_list(hit_elem.get("evalues", ""))[0],
                hit_type="chain_blast"
            )
            
            # Extract query and hit regions
            query_reg_elem = hit_elem.find("query_reg")
            if query_reg_elem is not None and query_reg_elem.text:
                hit.range = query_reg_elem.text.strip()
                hit.parse_ranges()
            
            hit_reg_elem = hit_elem.find("hit_reg")
            if hit_reg_elem is not None and hit_reg_elem.text:
                hit.hit_range = hit_reg_elem.text.strip()
            
            summary.chain_blast_hits.append(hit)
        
        # Process domain BLAST hits
        for hit_elem in root.findall(".//blast_run/hits/hit"):
            hit = BlastHit(
                hit_id=hit_elem.get("num", ""),
                domain_id=hit_elem.get("domain_id", ""),
                pdb_id=hit_elem.get("pdb_id", ""),
                chain_id=hit_elem.get("chain_id", ""),
                evalue=_parse_float_list(hit_elem.get("evalues", ""))[0],
                hit_type="domain_blast"
            )
            
            # Extract query and hit regions
            query_reg_elem = hit_elem.find("query_reg")
            if query_reg_elem is not None and query_reg_elem.text:
                hit.range = query_reg_elem.text.strip()
                hit.parse_ranges()
            
            hit_reg_elem = hit_elem.find("hit_reg")
            if hit_reg_elem is not None and hit_reg_elem.text:
                hit.hit_range = hit_reg_elem.text.strip()
            
            summary.domain_blast_hits.append(hit)
        
        # Process HHSearch hits
        for hit_elem in root.findall(".//hh_run/hits/hit"):
            hit = HHSearchHit(
                hit_id=hit_elem.get("num", ""),
                domain_id=hit_elem.get("domain_id", ""),
                pdb_id=hit_elem.get("pdb_id", ""),
                chain_id=hit_elem.get("chain_id", ""),
                probability=float(hit_elem.get("probability", "0")),
                evalue=float(hit_elem.get("evalue", "999")),
                score=float(hit_elem.get("score", "0")),
                hit_type="hhsearch"
            )
            
            # Extract query and hit regions
            query_reg_elem = hit_elem.find("query_reg")
            if query_reg_elem is not None and query_reg_elem.text:
                hit.range = query_reg_elem.text.strip()
                hit.parse_ranges()
            
            hit_reg_elem = hit_elem.find("hit_reg")
            if hit_reg_elem is not None and hit_reg_elem.text:
                hit.hit_range = hit_reg_elem.text.strip()
            
            summary.hhsearch_hits.append(hit)
        
        # Process error flags
        error_flags = [
            "no_chain_blast", "chain_blast_no_hits", 
            "no_domain_blast", "domain_blast_no_hits",
            "no_hhsearch", "hhsearch_error"
        ]
        
        for flag in error_flags:
            if blast_summ.get(flag, "false").lower() == "true":
                summary.has_errors = True
                summary.error_details[flag] = True
            
        return summary
        
    except Exception as e:
        # Create error summary
        error_summary = DomainSummary(
            pdb_id=file_path.split("/")[-1].split("_")[0],
            chain_id=file_path.split("/")[-1].split("_")[1].split(".")[0],
            reference="unknown",
            has_errors=True
        )
        error_summary.error_details["parse_error"] = True
        return error_summary

def domain_summary_to_xml(summary: DomainSummary) -> ET.Element:
    """Convert domain summary model to XML"""
    root = ET.Element("blast_summ_doc")
    
    # Create summary node
    blast_summ = ET.SubElement(root, "blast_summ")
    blast_summ.set("pdb", summary.pdb_id)
    blast_summ.set("chain", summary.chain_id)
    
    if summary.reference:
        blast_summ.set("reference", summary.reference)
    
    if summary.is_peptide:
        blast_summ.set("is_peptide", "true")
    
    # Add sequence length
    if summary.sequence_length > 0:
        query_len = ET.SubElement(root, "query_len")
        query_len.text = str(summary.sequence_length)
    
    # Add error flags
    for flag, value in summary.error_details.items():
        if value:
            blast_summ.set(flag, "true")
    
    # Rest of XML generation...
    # (Chain BLAST hits, domain BLAST hits, HHSearch hits)
    
    return root

def _parse_float_list(value_str: str) -> List[float]:
    """Parse comma-separated float values"""
    result = []
    if not value_str:
        return [999.0]  # Default high e-value
        
    for val in value_str.split(","):
        try:
            result.append(float(val))
        except ValueError:
            result.append(999.0)
    
    return result if result else [999.0]