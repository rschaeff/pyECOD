#!/usr/bin/env python3
"""
Debug tool for domain summary creation with model integration.

Usage:
    python debug_domain_summary.py --batch-id [BATCH_ID] --pdb [PDB_ID] --chain [CHAIN_ID] --verbose
"""

import os
import sys
import argparse
import logging
from typing import Dict, List, Any
import xml.etree.ElementTree as ET

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

# Import model classes
try:
    from ecod.models.pipeline import DomainSummaryModel, BlastHit, HHSearchHit
except ImportError:
    print("Warning: Could not import model classes. Running in standalone mode.")
    
    # Define stub model classes for standalone mode
    class BlastHit:
        def __init__(self, **kwargs):
            for key, value in kwargs.items():
                setattr(self, key, value)
            self.range_parsed = []
        
        def parse_ranges(self):
            pass
            
    class HHSearchHit:
        def __init__(self, **kwargs):
            for key, value in kwargs.items():
                setattr(self, key, value)
            self.range_parsed = []
            
        def parse_ranges(self):
            pass
    
    class DomainSummaryModel:
        def __init__(self, pdb_id, chain_id, reference, **kwargs):
            self.pdb_id = pdb_id
            self.chain_id = chain_id
            self.reference = reference
            self.chain_blast_hits = []
            self.domain_blast_hits = []
            self.hhsearch_hits = []
            self.errors = {}
            self.output_file_path = None
            self.processed = False
            
            for key, value in kwargs.items():
                setattr(self, key, value)
                
        def to_xml(self):
            root = ET.Element("blast_summ_doc")
            blast_summ = ET.SubElement(root, "blast_summ", pdb=self.pdb_id, chain=self.chain_id)
            
            # Add errors
            for key, value in self.errors.items():
                if value:
                    blast_summ.set(key, "true")
            
            # Add chain blast hits
            if self.chain_blast_hits:
                chain_blast_run = ET.SubElement(root, "chain_blast_run")
                chain_blast_run.set("program", "blastp")
                hits_node = ET.SubElement(chain_blast_run, "hits")
                
                for hit in self.chain_blast_hits:
                    hit_elem = ET.SubElement(hits_node, "hit")
                    hit_elem.set("num", hit.hit_id)
                    hit_elem.set("pdb_id", hit.pdb_id)
                    hit_elem.set("chain_id", hit.chain_id)
                    if hasattr(hit, "evalue"):
                        hit_elem.set("evalues", str(hit.evalue))
                        
                    query_reg = ET.SubElement(hit_elem, "query_reg")
                    query_reg.text = hit.range
            
            # Add domain blast hits
            if self.domain_blast_hits:
                domain_blast_run = ET.SubElement(root, "blast_run")
                domain_blast_run.set("program", "blastp")
                hits_node = ET.SubElement(domain_blast_run, "hits")
                
                for hit in self.domain_blast_hits:
                    hit_elem = ET.SubElement(hits_node, "hit")
                    if hit.domain_id:
                        hit_elem.set("domain_id", hit.domain_id)
                    hit_elem.set("pdb_id", hit.pdb_id)
                    hit_elem.set("chain_id", hit.chain_id)
                    if hasattr(hit, "evalue"):
                        hit_elem.set("evalues", str(hit.evalue))
                        
                    query_reg = ET.SubElement(hit_elem, "query_reg")
                    query_reg.text = hit.range
            
            # Add HHSearch hits
            if self.hhsearch_hits:
                hh_run = ET.SubElement(root, "hh_run")
                hh_run.set("program", "hhsearch")
                hits_node = ET.SubElement(hh_run, "hits")
                
                for i, hit in enumerate(self.hhsearch_hits):
                    hit_elem = ET.SubElement(hits_node, "hit")
                    if hit.domain_id:
                        hit_elem.set("domain_id", hit.domain_id)
                    hit_elem.set("hit_id", hit.hit_id)
                    hit_elem.set("num", str(i+1))
                    hit_elem.set("probability", str(hit.probability))
                    hit_elem.set("evalue", str(hit.evalue))
                    
                    query_reg = ET.SubElement(hit_elem, "query_reg")
                    query_reg.text = hit.range
            
            return root

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger("domain_summary_debug")

def parse_args():
    parser = argparse.ArgumentParser(description="Debug domain summary creation")
    parser.add_argument("--batch-id", required=True, help="Batch ID")
    parser.add_argument("--pdb", required=True, help="PDB ID")
    parser.add_argument("--chain", required=True, help="Chain ID")
    parser.add_argument("--reference", default="develop291", help="Reference version")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose output")
    parser.add_argument("--test-blast-only", action="store_true", help="Test only BLAST processing")
    parser.add_argument("--test-hhsearch-only", action="store_true", help="Test only HHSearch processing")
    parser.add_argument("--output", help="Path to write output files")
    
    return parser.parse_args()

def get_batch_path(batch_id):
    """Get the path to the batch directory."""
    # This would normally use the config to determine the batch path
    # For now, we're using a hardcoded path structure based on the logs
    return f"/data/ecod/pdb_updates/batches/{batch_id}"

def get_file_paths(batch_path, pdb_id, chain_id, reference):
    """Get the paths to all required files."""
    paths = {
        "fasta": {
            "standard_path": os.path.join(batch_path, "fastas", f"{pdb_id}_{chain_id}.fa"),
            "legacy_path": None,
            "exists_at": None
        },
        "a3m": {
            "standard_path": os.path.join(batch_path, "hhsearch", "profiles", f"{pdb_id}_{chain_id}.a3m"),
            "legacy_path": None,
            "exists_at": None
        },
        "hhm": {
            "standard_path": os.path.join(batch_path, "hhsearch", "profiles", f"{pdb_id}_{chain_id}.hhm"),
            "legacy_path": None,
            "exists_at": None
        },
        "hhr": {
            "standard_path": os.path.join(batch_path, "hhsearch", f"{pdb_id}_{chain_id}.{reference}.hhr"),
            "legacy_path": None,
            "exists_at": None
        },
        "hh_xml": {
            "standard_path": os.path.join(batch_path, "hhsearch", f"{pdb_id}_{chain_id}.{reference}.xml"),
            "legacy_path": None,
            "exists_at": None
        },
        "chain_blast": {
            "standard_path": os.path.join(batch_path, "blast", "chain", f"{pdb_id}_{chain_id}.{reference}.xml"),
            "legacy_path": os.path.join(batch_path, "blast", f"{pdb_id}_{chain_id}.{reference}.chain_blast.xml"),
            "exists_at": None
        },
        "domain_blast": {
            "standard_path": os.path.join(batch_path, "blast", "domain", f"{pdb_id}_{chain_id}.{reference}.xml"),
            "legacy_path": os.path.join(batch_path, "blast", f"{pdb_id}_{chain_id}.{reference}.domain_blast.xml"),
            "exists_at": None
        },
        "domain_summary": {
            "standard_path": os.path.join(batch_path, "domains", f"{pdb_id}_{chain_id}.{reference}.domain_summary.xml"),
            "legacy_path": None,
            "exists_at": None
        },
        "domain_partition": {
            "standard_path": os.path.join(batch_path, "domains", f"{pdb_id}_{chain_id}.{reference}.domains.xml"),
            "legacy_path": None,
            "exists_at": None
        }
    }
    
    # Check which files exist and update the exists_at field
    for file_type, file_paths in paths.items():
        # Check standard path
        if file_paths["standard_path"] and os.path.exists(file_paths["standard_path"]):
            file_paths["exists_at"] = file_paths["standard_path"]
        # Check legacy path
        elif file_paths["legacy_path"] and os.path.exists(file_paths["legacy_path"]):
            file_paths["exists_at"] = file_paths["legacy_path"]
    
    return paths

def analyze_xml_structure(file_path, level=0, max_levels=3):
    """Analyze the structure of an XML file for debugging."""
    try:
        tree = ET.parse(file_path)
        root = tree.getroot()
        return analyze_element_structure(root, level, max_levels)
    except Exception as e:
        return f"Error parsing XML: {str(e)}"

def analyze_element_structure(element, level=0, max_levels=3, path=""):
    """Recursively analyze the structure of an XML element."""
    if level > max_levels:
        return f"{' ' * level}...(truncated due to depth)\n"
    
    current_path = f"{path}/{element.tag}" if path else element.tag
    result = f"{' ' * level}{element.tag}\n"
    
    if element.attrib:
        result += f"{' ' * (level+2)}Attributes: {dict(element.attrib)}\n"
    
    if element.text and element.text.strip():
        text_preview = element.text.strip()[:50]
        if len(element.text.strip()) > 50:
            text_preview += "..."
        result += f"{' ' * (level+2)}Text: {text_preview}\n"
    
    children_tags = {}
    for child in element:
        if child.tag in children_tags:
            children_tags[child.tag] += 1
        else:
            children_tags[child.tag] = 1
    
    if children_tags:
        result += f"{' ' * (level+2)}Children: {children_tags}\n"
    
    # Recursively analyze a limited number of each child type
    analyzed_tags = set()
    for child in element:
        if child.tag not in analyzed_tags:
            result += analyze_element_structure(child, level + 4, max_levels, current_path)
            analyzed_tags.add(child.tag)
    
    return result

def count_hits_in_xml(file_path, hit_paths):
    """Count hits in an XML file using multiple possible paths."""
    try:
        tree = ET.parse(file_path)
        root = tree.getroot()
        
        results = {}
        for path_name, xpath in hit_paths.items():
            hits = root.findall(xpath)
            results[path_name] = len(hits)
        
        return results
    except Exception as e:
        return f"Error counting hits: {str(e)}"

def process_blast_xml(file_path, blast_type):
    """Process BLAST XML and return hits in dictionary format."""
    if not file_path or not os.path.exists(file_path):
        logger.warning(f"{blast_type.capitalize()} BLAST file not found at {file_path}")
        return []
    
    logger.info(f"Processing {blast_type} BLAST from {file_path}")
    
    # First, analyze the XML structure
    structure = analyze_xml_structure(file_path)
    logger.debug(f"{blast_type.capitalize()} BLAST XML structure:\n{structure}")
    
    # Count hits using various possible paths
    hit_paths = {
        "direct_hits": ".//hit",
        "chain_blast_hits": ".//chain_blast_run/hits/hit",
        "blast_run_hits": ".//blast_run/hits/hit",
        "iteration_hits": ".//BlastOutput_iterations/Iteration/Iteration_hits/Hit"
    }
    hit_counts = count_hits_in_xml(file_path, hit_paths)
    logger.info(f"{blast_type.capitalize()} BLAST hit counts: {hit_counts}")
    
    # Parse the XML file
    try:
        tree = ET.parse(file_path)
        root = tree.getroot()
    except Exception as e:
        logger.error(f"Error parsing {blast_type} BLAST XML: {str(e)}")
        return []
    
    # Extract hits based on the XML structure
    hits = []
    hit_elements = []
    
    # Try multiple possible structures
    if blast_type == "chain":
        # Look for chainwise blast runs
        chain_blast_runs = root.findall(".//chain_blast_run")
        if chain_blast_runs:
            for run in chain_blast_runs:
                hit_elements.extend(run.findall(".//hit"))
                
        # Alternative structure
        if not hit_elements:
            hit_elements = root.findall(".//chain_blast_run/hits/hit")
        
        # Another alternative structure
        if not hit_elements:
            hit_elements = root.findall(".//BlastOutput_iterations/Iteration/Iteration_hits/Hit")
    else:
        # Look for domain blast runs
        domain_blast_runs = root.findall(".//blast_run")
        if domain_blast_runs:
            for run in domain_blast_runs:
                hit_elements.extend(run.findall(".//hit"))
                
        # Alternative structure
        if not hit_elements:
            hit_elements = root.findall(".//blast_run/hits/hit")
        
        # Another alternative structure
        if not hit_elements:
            hit_elements = root.findall(".//BlastOutput_iterations/Iteration/Iteration_hits/Hit")
    
    # Last resort - try to find any hit elements
    if not hit_elements:
        hit_elements = root.findall(".//hit")
    
    logger.info(f"Found {len(hit_elements)} hit elements in {blast_type} BLAST")
    
    # Process each hit element
    for i, hit_element in enumerate(hit_elements):
        try:
            hit_dict = element_to_dict(hit_element)
            
            # Extract query region if available
            query_reg = hit_element.find("query_reg")
            if query_reg is not None and query_reg.text:
                hit_dict["range"] = query_reg.text.strip()
                hit_dict["range_parsed"] = parse_range(query_reg.text.strip())
                logger.debug(f"Hit {i+1} query region: {query_reg.text.strip()}")
            else:
                logger.debug(f"Hit {i+1} has no query_reg element")
            
            # Extract hit region if available
            hit_reg = hit_element.find("hit_reg")
            if hit_reg is not None and hit_reg.text:
                hit_dict["hit_range"] = hit_reg.text.strip()
            
            # Extract other important fields if available
            for field in ["evalue", "pident", "qcov", "score"]:
                field_element = hit_element.find(field)
                if field_element is not None and field_element.text:
                    hit_dict[field] = field_element.text.strip()
                    logger.debug(f"Hit {i+1} {field}: {field_element.text.strip()}")
            
            # Add blast type to the hit
            hit_dict["source"] = blast_type
            hit_dict["hit_type"] = f"{blast_type}_blast"
            
            # Set PDB and chain IDs if available
            hit_dict["pdb_id"] = hit_dict.get("pdb_id", hit_element.get("pdb_id", ""))
            hit_dict["chain_id"] = hit_dict.get("chain_id", hit_element.get("chain_id", ""))
            hit_dict["hit_id"] = hit_dict.get("hit_id", hit_element.get("num", str(i+1)))
            hit_dict["domain_id"] = hit_dict.get("domain_id", hit_element.get("domain_id", ""))
            
            # Print summary of the hit
            hit_summary = {
                "hit_id": hit_dict.get("hit_id", "unknown"),
                "pdb_id": hit_dict.get("pdb_id", "N/A"),
                "chain_id": hit_dict.get("chain_id", "N/A"),
                "range": hit_dict.get("range", "N/A"),
                "evalue": hit_dict.get("evalue", "N/A")
            }
            logger.info(f"Hit {i+1} summary: {hit_summary}")
            
            # Add the hit to the list
            hits.append(hit_dict)
        except Exception as e:
            logger.error(f"Error processing hit element {i+1}: {str(e)}", exc_info=True)
    
    logger.info(f"Processed {len(hits)}/{len(hit_elements)} hits from {blast_type} BLAST")
    return hits

def process_hhsearch_xml(file_path):
    """Process HHSearch XML and return hits in dictionary format."""
    if not file_path or not os.path.exists(file_path):
        logger.warning(f"HHSearch file not found at {file_path}")
        return []
    
    logger.info(f"Processing HHSearch from {file_path}")
    
    # First, analyze the XML structure
    structure = analyze_xml_structure(file_path)
    logger.debug(f"HHSearch XML structure:\n{structure}")
    
    # Count hits using various possible paths
    hit_paths = {
        "direct_hits": ".//hit",
        "hits_hit": ".//hits/hit",
        "alignments": ".//alignment"
    }
    hit_counts = count_hits_in_xml(file_path, hit_paths)
    logger.info(f"HHSearch hit counts: {hit_counts}")
    
    # Parse the XML file
    try:
        tree = ET.parse(file_path)
        root = tree.getroot()
    except Exception as e:
        logger.error(f"Error parsing HHSearch XML: {str(e)}")
        return []
    
    # Extract hits based on the XML structure
    hits = []
    hit_elements = []
    
    # Try multiple possible structures
    hit_elements = root.findall(".//hit")
    
    if not hit_elements:
        hit_elements = root.findall(".//hits/hit")
    
    if not hit_elements:
        hit_elements = root.findall(".//alignment")
    
    logger.info(f"Found {len(hit_elements)} hit elements in HHSearch")
    
    # Process each hit element
    for i, hit_element in enumerate(hit_elements):
        try:
            hit_dict = element_to_dict(hit_element)
            
            # Extract query range if available
            query_range = None
            
            # Different possible locations for query range
            query_from_elem = hit_element.find(".//query_from")
            query_to_elem = hit_element.find(".//query_to")
            
            if query_from_elem is not None and query_to_elem is not None:
                query_from = query_from_elem.text.strip()
                query_to = query_to_elem.text.strip()
                if query_from and query_to:
                    query_range = f"{query_from}-{query_to}"
                    logger.debug(f"Hit {i+1} query range from query_from/query_to: {query_range}")
            
            # Try alternate query range format
            if not query_range:
                query_range_elem = hit_element.find(".//query_range")
                if query_range_elem is not None and query_range_elem.text:
                    query_range = query_range_elem.text.strip()
                    logger.debug(f"Hit {i+1} query range from query_range: {query_range}")
            
            if query_range:
                hit_dict["range"] = query_range
                hit_dict["range_parsed"] = parse_range(query_range)
            else:
                logger.debug(f"Hit {i+1} has no query range")
            
            # Extract hit range if available
            hit_range = None
            hit_from_elem = hit_element.find(".//template_from")
            hit_to_elem = hit_element.find(".//template_to")
            
            if hit_from_elem is not None and hit_to_elem is not None:
                hit_from = hit_from_elem.text.strip()
                hit_to = hit_to_elem.text.strip()
                if hit_from and hit_to:
                    hit_range = f"{hit_from}-{hit_to}"
                    hit_dict["hit_range"] = hit_range
            
            # Extract probability if available
            prob_elem = hit_element.find(".//probability")
            if prob_elem is not None and prob_elem.text:
                hit_dict["probability"] = prob_elem.text.strip()
                logger.debug(f"Hit {i+1} probability: {prob_elem.text.strip()}")
            
            # Extract E-value if available
            evalue_elem = hit_element.find(".//evalue")
            if evalue_elem is not None and evalue_elem.text:
                hit_dict["evalue"] = evalue_elem.text.strip()
                logger.debug(f"Hit {i+1} evalue: {evalue_elem.text.strip()}")
            
            # Extract score if available
            score_elem = hit_element.find(".//score")
            if score_elem is not None and score_elem.text:
                hit_dict["score"] = score_elem.text.strip()
            
            # Extract domain ID and hit ID
            if "id" in hit_dict:
                hit_dict["hit_id"] = hit_dict["id"]
                
            if "domain_id" not in hit_dict and "hit_id" in hit_dict:
                # Try to extract domain ID from hit ID
                hit_dict["domain_id"] = hit_dict["hit_id"]
            
            # Add source to the hit
            hit_dict["source"] = "hhsearch"
            
            # Print summary of the hit
            hit_summary = {
                "hit_id": hit_dict.get("hit_id", "unknown"),
                "domain_id": hit_dict.get("domain_id", "N/A"),
                "range": hit_dict.get("range", "N/A"),
                "probability": hit_dict.get("probability", "N/A"),
                "evalue": hit_dict.get("evalue", "N/A")
            }
            logger.info(f"Hit {i+1} summary: {hit_summary}")
            
            # Add the hit to the list
            hits.append(hit_dict)
        except Exception as e:
            logger.error(f"Error processing hit element {i+1}: {str(e)}", exc_info=True)
    
    logger.info(f"Processed {len(hits)}/{len(hit_elements)} hits from HHSearch")
    return hits

def element_to_dict(element):
    """Convert an XML element to a dictionary, extracting attributes and child elements."""
    # Start with attributes
    result = dict(element.attrib)
    
    # Add simple child elements (those with just text content)
    for child in element:
        if len(child) == 0 and child.text and child.text.strip():
            result[child.tag] = child.text.strip()
    
    return result

def parse_range(range_str):
    """Parse a range string into a list of (start, end) tuples."""
    if not range_str:
        return []
    
    ranges = []
    parts = range_str.split(',')
    
    for part in parts:
        if '-' in part:
            try:
                start_str, end_str = part.split('-')
                # Remove any non-digit characters
                start = int(''.join(c for c in start_str if c.isdigit()))
                end = int(''.join(c for c in end_str if c.isdigit()))
                if start <= end:
                    ranges.append((start, end))
            except ValueError:
                logger.warning(f"Invalid range format: {part}")
    
    return ranges

def convert_dict_to_blast_hit(hit_dict):
    """Convert a hit dictionary to a BlastHit object."""
    hit = BlastHit(
        hit_id=hit_dict.get("hit_id", ""),
        domain_id=hit_dict.get("domain_id", ""),
        pdb_id=hit_dict.get("pdb_id", ""),
        chain_id=hit_dict.get("chain_id", ""),
        hit_type=hit_dict.get("hit_type", ""),
        range=hit_dict.get("range", ""),
        hit_range=hit_dict.get("hit_range", "")
    )
    
    # Set range_parsed if available
    if "range_parsed" in hit_dict:
        hit.range_parsed = hit_dict["range_parsed"]
    
    # Set numeric values
    if "evalue" in hit_dict:
        try:
            hit.evalue = float(hit_dict["evalue"])
        except (ValueError, TypeError):
            hit.evalue = 999.0
    
    return hit

def convert_dict_to_hhsearch_hit(hit_dict):
    """Convert a hit dictionary to an HHSearchHit object."""
    hit = HHSearchHit(
        hit_id=hit_dict.get("hit_id", ""),
        domain_id=hit_dict.get("domain_id", ""),
        range=hit_dict.get("range", ""),
        hit_range=hit_dict.get("hit_range", "")
    )
    
    # Set range_parsed if available
    if "range_parsed" in hit_dict:
        hit.range_parsed = hit_dict["range_parsed"]
    
    # Set numeric values
    if "probability" in hit_dict:
        try:
            hit.probability = float(hit_dict["probability"])
        except (ValueError, TypeError):
            hit.probability = 0.0
    
    if "evalue" in hit_dict:
        try:
            hit.evalue = float(hit_dict["evalue"])
        except (ValueError, TypeError):
            hit.evalue = 999.0
    
    if "score" in hit_dict:
        try:
            hit.score = float(hit_dict["score"])
        except (ValueError, TypeError):
            hit.score = 0.0
    
    return hit

def create_domain_summary_model(pdb_id, chain_id, reference, file_paths):
    """Create a domain summary model using the model-based approach."""
    logger.info(f"Creating domain summary model for {pdb_id}_{chain_id} with reference {reference}")
    
    # Initialize domain summary model
    domain_summary = DomainSummaryModel(
        pdb_id=pdb_id,
        chain_id=chain_id,
        reference=reference
    )
    
    # Set output file path
    output_path = file_paths.get('domain_summary', {}).get('standard_path')
    domain_summary.output_file_path = output_path
    
    # Track evidence counts for logging
    evidence_counts = {
        "chain_blast": 0,
        "domain_blast": 0,
        "hhsearch": 0,
        "total": 0
    }
    
    # Process chain BLAST results
    chain_blast_path = file_paths.get('chain_blast', {}).get('exists_at')
    if chain_blast_path:
        chain_blast_hits = process_blast_xml(chain_blast_path, "chain")
        logger.info(f"Found {len(chain_blast_hits)} chain BLAST hits")
        
        # Convert dictionary hits to BlastHit objects
        for hit_dict in chain_blast_hits:
            hit_dict["hit_type"] = "chain_blast"  # Ensure hit type is set
            blast_hit = convert_dict_to_blast_hit(hit_dict)
            domain_summary.chain_blast_hits.append(blast_hit)
            evidence_counts["chain_blast"] += 1
            evidence_counts["total"] += 1
    
    # Process domain BLAST results
    domain_blast_path = file_paths.get('domain_blast', {}).get('exists_at')
    if domain_blast_path:
        domain_blast_hits = process_blast_xml(domain_blast_path, "domain")
        logger.info(f"Found {len(domain_blast_hits)} domain BLAST hits")
        
        # Convert dictionary hits to BlastHit objects
        for hit_dict in domain_blast_hits:
            hit_dict["hit_type"] = "domain_blast"  # Ensure hit type is set
            blast_hit = convert_dict_to_blast_hit(hit_dict)
            domain_summary.domain_blast_hits.append(blast_hit)
            evidence_counts["domain_blast"] += 1
            evidence_counts["total"] += 1
    
    # Process HHSearch results
    hhsearch_path = file_paths.get('hh_xml', {}).get('exists_at')
    if hhsearch_path:
        hhsearch_hits = process_hhsearch_xml(hhsearch_path)
        logger.info(f"Found {len(hhsearch_hits)} HHSearch hits")
        
        # Convert dictionary hits to HHSearchHit objects
        for hit_dict in hhsearch_hits:
            hh_hit = convert_dict_to_hhsearch_hit(hit_dict)
            domain_summary.hhsearch_hits.append(hh_hit)
            evidence_counts["hhsearch"] += 1
            evidence_counts["total"] += 1
    
    # Mark as processed
    domain_summary.processed = True
    
    # Log summary statistics
    logger.info(f"Domain summary statistics for {pdb_id}_{chain_id}:")
    logger.info(f"  Chain BLAST evidence: {evidence_counts['chain_blast']}")
    logger.info(f"  Domain BLAST evidence: {evidence_counts['domain_blast']}")
    logger.info(f"  HHSearch evidence: {evidence_counts['hhsearch']}")
    logger.info(f"  Total evidence: {evidence_counts['total']}")
    
    # Write the XML to file
    if output_path:
        try:
            # Convert domain summary model to XML
            root = domain_summary.to_xml()
            
            # Write XML to file
            tree = ET.ElementTree(root)
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            tree.write(output_path, encoding="utf-8", xml_declaration=True)
            logger.info(f"Wrote domain summary to {output_path} with {evidence_counts['total']} evidence items")
        except Exception as e:
            logger.error(f"Error writing domain summary: {str(e)}", exc_info=True)
    else:
        logger.error("No output path provided for domain summary")
    
    return domain_summary

def main():
    args = parse_args()
    
    # Set log level based on verbosity
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    logger.info(f"Debugging domain summary creation for {args.pdb}_{args.chain} in batch {args.batch_id}")
    
    # Get batch path
    batch_path = get_batch_path(args.batch_id)
    logger.info(f"Batch path: {batch_path}")
    
    # Get file paths
    file_paths = get_file_paths(batch_path, args.pdb, args.chain, args.reference)
    
    # Display file status
    logger.info("File status:")
    for file_type, paths in file_paths.items():
        exists = "exists" if paths["exists_at"] else "missing"
        logger.info(f"  {file_type}: {exists} at {paths['exists_at'] or paths['standard_path']}")
    
    # Test specific components if requested
    if args.test_blast_only:
        logger.info("Testing BLAST processing only")
        chain_blast_path = file_paths.get('chain_blast', {}).get('exists_at')
        domain_blast_path = file_paths.get('domain_blast', {}).get('exists_at')
        
        if chain_blast_path:
            logger.info("Processing chain BLAST file")
            chain_hits = process_blast_xml(chain_blast_path, "chain")
            logger.info(f"Found {len(chain_hits)} chain BLAST hits")
        
        if domain_blast_path:
            logger.info("Processing domain BLAST file")
            domain_hits = process_blast_xml(domain_blast_path, "domain")
            logger.info(f"Found {len(domain_hits)} domain BLAST hits")
    
    elif args.test_hhsearch_only:
        logger.info("Testing HHSearch processing only")
        hhsearch_path = file_paths.get('hh_xml', {}).get('exists_at')
        
        if hhsearch_path:
            logger.info("Processing HHSearch file")
            hhsearch_hits = process_hhsearch_xml(hhsearch_path)
            logger.info(f"Found {len(hhsearch_hits)} HHSearch hits")
    
    else:
        # Create domain summary
        output_path = args.output or file_paths.get('domain_summary', {}).get('standard_path')
        
        logger.info(f"Creating domain summary model and writing to {output_path}")
        domain_summary = create_domain_summary_model(args.pdb, args.chain, args.reference, file_paths)
        
        # Print model stats
        logger.info("Domain summary model created:")
        logger.info(f"  Chain BLAST hits: {len(domain_summary.chain_blast_hits)}")
        logger.info(f"  Domain BLAST hits: {len(domain_summary.domain_blast_hits)}")
        logger.info(f"  HHSearch hits: {len(domain_summary.hhsearch_hits)}")
        logger.info(f"  Output file: {domain_summary.output_file_path}")
        logger.info(f"  Processed: {domain_summary.processed}")

if __name__ == "__main__":
    main()
