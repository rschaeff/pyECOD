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
            """Convert to XML Element"""
            root = ET.Element("blast_summ_doc")

            # Create summary node
            blast_summ = ET.SubElement(root, "blast_summ", pdb=self.pdb_id, chain=self.chain_id)
            
            # Add errors
            for error, value in self.errors.items():
                if value:
                    blast_summ.set(error, "true")
            
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
                    hit_elem.set("hsp_count", str(hit.hsp_count))

                    # Format evalues correctly
                    if hasattr(hit, "evalues") and hit.evalues:
                        hit_elem.set("evalues", ",".join(str(e) for e in hit.evalues))
                    elif hasattr(hit, "evalue"):
                        hit_elem.set("evalues", str(hit.evalue))

                    # Add query region
                    query_reg = ET.SubElement(hit_elem, "query_reg")
                    query_reg.text = hit.range

                    # Add hit region
                    hit_reg = ET.SubElement(hit_elem, "hit_reg")
                    hit_reg.text = hit.hit_range
            
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
                    hit_elem.set("hsp_count", str(hit.hsp_count))

                    # Format evalues correctly
                    if hasattr(hit, "evalues") and hit.evalues:
                        hit_elem.set("evalues", ",".join(str(e) for e in hit.evalues))
                    elif hasattr(hit, "evalue"):
                        hit_elem.set("evalues", str(hit.evalue))

                    # Add discontinuous flag if applicable
                    if hasattr(hit, "discontinuous") and hit.discontinuous:
                        hit_elem.set("discontinuous", "true")

                    # Add query region
                    query_reg = ET.SubElement(hit_elem, "query_reg")
                    query_reg.text = hit.range

                    # Add hit region
                    hit_reg = ET.SubElement(hit_elem, "hit_reg")
                    hit_reg.text = hit.hit_range
            
            # Add HHSearch hits
            if self.hhsearch_hits:
                hh_run = ET.SubElement(root, "hh_run")
                hh_run.set("program", "hhsearch")

                hits_node = ET.SubElement(hh_run, "hits")
                for i, hit in enumerate(self.hhsearch_hits):
                    hit_elem = ET.SubElement(hits_node, "hit")

                    # Set required attributes
                    if hit.domain_id:
                        hit_elem.set("domain_id", hit.domain_id)
                    hit_elem.set("hit_id", hit.hit_id)
                    hit_elem.set("num", str(i+1))
                    hit_elem.set("probability", str(hit.probability))
                    hit_elem.set("evalue", str(hit.evalue))
                    hit_elem.set("score", str(hit.score))
                    
                    # Add query region
                    query_reg = ET.SubElement(hit_elem, "query_reg")
                    query_reg.text = hit.range

                    # Add hit region if available
                    if hasattr(hit, "hit_range") and hit.hit_range:
                        hit_reg = ET.SubElement(hit_elem, "hit_reg")
                        hit_reg.text = hit.hit_range
            
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
    return f"/data/ecod/pdb_updates/batches/alt_rep_batch_006_20250424_1622"

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
            "legacy_path": os.path.join(batch_path, "hhsearch", f"{pdb_id}_{chain_id}.{reference}.hhsearch.xml"),
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
        hits = []

        # Check if this is standard NCBI BLAST output format
        if root.tag == "BlastOutput" or root.find(".//BlastOutput_iterations") is not None:
            # Process standard NCBI BLAST output
            logger.debug(f"Processing standard NCBI BLAST output format")
            hit_elements = root.findall(".//BlastOutput_iterations/Iteration/Iteration_hits/Hit")

            for i, hit_elem in enumerate(hit_elements):
                try:
                    hit_data = {}
                    hit_data["hit_id"] = str(i+1)  # Use position as hit_id if not found

                    # Get Hit_num if available
                    hit_num = hit_elem.find("Hit_num")
                    if hit_num is not None and hit_num.text:
                        hit_data["hit_id"] = hit_num.text

                    # Get Hit_def for PDB ID and chain ID
                    hit_def = hit_elem.find("Hit_def")
                    if hit_def is not None and hit_def.text:
                        hit_def_text = hit_def.text.strip()
                        logger.debug(f"Hit {i+1} definition: {hit_def_text}")

                        if blast_type == "chain":
                            # For chain BLAST, format is typically "2lxo A"
                            parts = hit_def_text.split()
                            if len(parts) >= 2:
                                hit_data["pdb_id"] = parts[0].lower()
                                hit_data["chain_id"] = parts[1]
                        else:
                            # For domain BLAST, format is typically "e2lxoA1 A:1-44 001088464"
                            parts = hit_def_text.split()
                            if len(parts) >= 1:
                                hit_data["domain_id"] = parts[0]
                                # Extract PDB ID from domain ID if it starts with e/d/g/x
                                if len(parts[0]) >= 5 and parts[0][0] in "edgx":
                                    hit_data["pdb_id"] = parts[0][1:5].lower()
                                    if len(parts[0]) > 5:
                                        hit_data["chain_id"] = parts[0][5]

                    # Get HSPs
                    hsps = hit_elem.findall("Hit_hsps/Hsp")
                    if hsps:
                        # Count HSPs
                        hit_data["hsp_count"] = len(hsps)

                        # Get E-values from HSPs
                        evalues = []
                        for hsp in hsps:
                            evalue_elem = hsp.find("Hsp_evalue")
                            if evalue_elem is not None and evalue_elem.text:
                                try:
                                    evalues.append(float(evalue_elem.text))
                                except ValueError:
                                    pass

                        if evalues:
                            hit_data["evalues"] = evalues
                            hit_data["evalue"] = min(evalues)  # Use lowest E-value

                        # Extract range from first HSP (best hit)
                        first_hsp = hsps[0]
                        q_from = first_hsp.find("Hsp_query-from")
                        q_to = first_hsp.find("Hsp_query-to")

                        if q_from is not None and q_to is not None and q_from.text and q_to.text:
                            hit_data["range"] = f"{q_from.text}-{q_to.text}"
                            hit_data["range_parsed"] = [(int(q_from.text), int(q_to.text))]

                        # Extract hit range from first HSP
                        h_from = first_hsp.find("Hsp_hit-from")
                        h_to = first_hsp.find("Hsp_hit-to")

                        if h_from is not None and h_to is not None and h_from.text and h_to.text:
                            hit_data["hit_range"] = f"{h_from.text}-{h_to.text}"

                        # Check if discontinuous (multiple non-contiguous HSPs)
                        if len(hsps) > 1:
                            ranges = []
                            for hsp in hsps:
                                q_from = hsp.find("Hsp_query-from")
                                q_to = hsp.find("Hsp_query-to")
                                if q_from is not None and q_to is not None and q_from.text and q_to.text:
                                    ranges.append((int(q_from.text), int(q_to.text)))

                            # Sort ranges by start position
                            ranges.sort(key=lambda x: x[0])

                            # Check if ranges are discontinuous
                            discontinuous = False
                            for j in range(1, len(ranges)):
                                if ranges[j][0] > ranges[j-1][1] + 1:
                                    discontinuous = True
                                    break

                            if discontinuous:
                                hit_data["discontinuous"] = True
                                hit_data["range"] = ",".join(f"{start}-{end}" for start, end in ranges)
                                hit_data["range_parsed"] = ranges

                    # Add metadata
                    hit_data["source"] = blast_type
                    hit_data["hit_type"] = f"{blast_type}_blast"

                    # Set defaults for required fields
                    hit_data.setdefault("pdb_id", "")
                    hit_data.setdefault("chain_id", "")
                    hit_data.setdefault("domain_id", "")
                    hit_data.setdefault("range", "")
                    hit_data.setdefault("hit_range", "")
                    hit_data.setdefault("evalue", 999.0)
                    hit_data.setdefault("hsp_count", 1)

                    # Summary of the hit for logging
                    hit_summary = {
                        "hit_id": hit_data["hit_id"],
                        "pdb_id": hit_data["pdb_id"],
                        "chain_id": hit_data["chain_id"],
                        "range": hit_data["range"],
                        "evalue": hit_data["evalue"]
                    }
                    logger.info(f"Hit {i+1} summary: {hit_summary}")

                    hits.append(hit_data)

                except Exception as e:
                    logger.error(f"Error processing hit element {i+1}: {str(e)}", exc_info=True)
                    continue
        else:
            # Try to process already-formatted domain summary XML
            logger.debug(f"Processing domain summary format XML")

            # Extract hits based on the XML structure
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
            else:
                # Look for domain blast runs
                domain_blast_runs = root.findall(".//blast_run")
                if domain_blast_runs:
                    for run in domain_blast_runs:
                        hit_elements.extend(run.findall(".//hit"))

                # Alternative structure
                if not hit_elements:
                    hit_elements = root.findall(".//blast_run/hits/hit")

            # Last resort - try to find any hit elements
            if not hit_elements:
                hit_elements = root.findall(".//hit")

            logger.info(f"Found {len(hit_elements)} hit elements in {blast_type} BLAST (domain summary format)")

            # Process each hit element
            for i, hit_element in enumerate(hit_elements):
                try:
                    # Extract attributes
                    hit_data = dict(hit_element.attrib)

                    # Extract query region if available
                    query_reg = hit_element.find("query_reg")
                    if query_reg is not None and query_reg.text:
                        hit_data["range"] = query_reg.text.strip()
                        hit_data["range_parsed"] = parse_range(query_reg.text.strip())

                    # Extract hit region if available
                    hit_reg = hit_element.find("hit_reg")
                    if hit_reg is not None and hit_reg.text:
                        hit_data["hit_range"] = hit_reg.text.strip()

                    # Add metadata
                    hit_data["source"] = blast_type
                    hit_data["hit_type"] = f"{blast_type}_blast"

                    # Set values from attributes
                    hit_data["pdb_id"] = hit_data.get("pdb_id", hit_element.get("pdb_id", ""))
                    hit_data["chain_id"] = hit_data.get("chain_id", hit_element.get("chain_id", ""))
                    hit_data["hit_id"] = hit_data.get("hit_id", hit_element.get("num", str(i+1)))
                    hit_data["domain_id"] = hit_data.get("domain_id", hit_element.get("domain_id", ""))

                    # Parse E-value
                    if "evalues" in hit_data:
                        try:
                            if "," in hit_data["evalues"]:
                                evalues = [float(e) for e in hit_data["evalues"].split(",")]
                                hit_data["evalue"] = min(evalues)
                            else:
                                hit_data["evalue"] = float(hit_data["evalues"])
                        except ValueError:
                            hit_data["evalue"] = 999.0

                    # Summary for logging
                    hit_summary = {
                        "hit_id": hit_data.get("hit_id", "unknown"),
                        "pdb_id": hit_data.get("pdb_id", "N/A"),
                        "chain_id": hit_data.get("chain_id", "N/A"),
                        "range": hit_data.get("range", "N/A"),
                        "evalue": hit_data.get("evalue", "N/A")
                    }
                    logger.info(f"Hit {i+1} summary: {hit_summary}")

                    hits.append(hit_data)
                except Exception as e:
                    logger.error(f"Error processing hit element {i+1}: {str(e)}", exc_info=True)
                    continue

        logger.info(f"Processed {len(hits)} hits for {blast_type} BLAST")
        return hits

    except Exception as e:
        logger.error(f"Error parsing {blast_type} BLAST XML: {str(e)}")
        return []

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
        "hh_hits": ".//hh_hit_list/hh_hit",
        "alignments": ".//alignment"
    }
    hit_counts = count_hits_in_xml(file_path, hit_paths)
    logger.info(f"HHSearch hit counts: {hit_counts}")

    # Parse the XML file
    try:
        tree = ET.parse(file_path)
        root = tree.getroot()
        hits = []

        # First, try to find hh_hit elements in hh_hit_list (preprocessed format)
        hit_elements = root.findall(".//hh_hit_list/hh_hit")

        if hit_elements:
            logger.debug(f"Processing preprocessed HHSearch XML format with {len(hit_elements)} hh_hit elements")

            for i, hit_elem in enumerate(hit_elements):
                try:
                    # Extract attributes
                    hit_data = dict(hit_elem.attrib)

                    # Extract query range
                    query_range = hit_elem.find("query_range")
                    if query_range is not None:
                        # Get range from text content
                        hit_data["range"] = query_range.text.strip()

                        # Also grab start/end attributes if needed
                        if "start" in query_range.attrib and "end" in query_range.attrib:
                            start = query_range.get("start")
                            end = query_range.get("end")
                            if not hit_data["range"]:
                                hit_data["range"] = f"{start}-{end}"

                    # Extract template range
                    template_range = hit_elem.find("template_seqid_range")
                    if template_range is not None:
                        hit_data["hit_range"] = template_range.text.strip()

                    # Convert e_value to evalue
                    if "e_value" in hit_data:
                        try:
                            hit_data["evalue"] = float(hit_data["e_value"])
                            # Remove original key to avoid duplication
                            del hit_data["e_value"]
                        except ValueError:
                            hit_data["evalue"] = 999.0

                    # Convert probability to float
                    if "probability" in hit_data:
                        try:
                            hit_data["probability"] = float(hit_data["probability"])
                        except ValueError:
                            hit_data["probability"] = 0.0

                    # Convert score to float
                    if "score" in hit_data:
                        try:
                            hit_data["score"] = float(hit_data["score"])
                        except ValueError:
                            hit_data["score"] = 0.0

                    # Ensure domain_id is set if ecod_domain_id is available
                    if "domain_id" not in hit_data and "ecod_domain_id" in hit_data:
                        hit_data["domain_id"] = hit_data["ecod_domain_id"]

                    # Make sure hit_id is set
                    if "hit_id" not in hit_data:
                        if "ecod_domain_id" in hit_data:
                            hit_data["hit_id"] = hit_data["ecod_domain_id"]
                        elif "hit_num" in hit_data:
                            hit_data["hit_id"] = hit_data["hit_num"]
                        else:
                            hit_data["hit_id"] = f"hit_{i+1}"

                    # Set defaults for missing fields
                    hit_data.setdefault("probability", 0.0)
                    hit_data.setdefault("evalue", 999.0)
                    hit_data.setdefault("score", 0.0)
                    hit_data.setdefault("range", "")
                    hit_data.setdefault("hit_range", "")
                    hit_data.setdefault("domain_id", "")

                    hit_summary = {
                        "hit_id": hit_data.get("hit_id", "unknown"),
                        "domain_id": hit_data.get("domain_id", "N/A"),
                        "range": hit_data.get("range", "N/A"),
                        "probability": hit_data.get("probability", "N/A"),
                        "evalue": hit_data.get("evalue", "N/A")
                    }
                    logger.info(f"Hit {i+1} summary: {hit_summary}")

                    hits.append(hit_data)
                except Exception as e:
                    logger.error(f"Error processing HH hit element {i+1}: {str(e)}", exc_info=True)
                    continue

        else:
            # Try other element types
            hit_elements = root.findall(".//hit")

            if not hit_elements:
                hit_elements = root.findall(".//hits/hit")

            if not hit_elements:
                hit_elements = root.findall(".//alignment")

            logger.info(f"Found {len(hit_elements)} hit elements in alternative HHSearch formats")

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

                    # Extract probability/evalue/score and convert to float
                    for field, default in [("probability", 0.0), ("evalue", 999.0), ("score", 0.0)]:
                        field_elem = hit_element.find(f".//{field}")
                        if field_elem is not None and field_elem.text:
                            try:
                                hit_dict[field] = float(field_elem.text.strip())
                            except ValueError:
                                hit_dict[field] = default

                    # Set hit_id
                    if "hit_id" not in hit_dict:
                        if "id" in hit_dict:
                            hit_dict["hit_id"] = hit_dict["id"]
                        else:
                            hit_dict["hit_id"] = f"hit_{i+1}"

                    # Set domain_id if missing
                    if "domain_id" not in hit_dict and "hit_id" in hit_dict:
                        hit_dict["domain_id"] = hit_dict["hit_id"]

                    hit_summary = {
                        "hit_id": hit_dict.get("hit_id", "unknown"),
                        "domain_id": hit_dict.get("domain_id", "N/A"),
                        "range": hit_dict.get("range", "N/A"),
                        "probability": hit_dict.get("probability", "N/A"),
                        "evalue": hit_dict.get("evalue", "N/A")
                    }
                    logger.info(f"Hit {i+1} summary: {hit_summary}")

                    hits.append(hit_dict)
                except Exception as e:
                    logger.error(f"Error processing hit element {i+1}: {str(e)}", exc_info=True)
                    continue

        logger.info(f"Processed {len(hits)} hits from HHSearch")
        return hits

    except Exception as e:
        logger.error(f"Error processing HHSearch XML: {str(e)}")
        return []

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
    """Convert a hit dictionary to a BlastHit object with all required attributes."""
    # Create basic hit object
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

    # Set list of evalues if available
    if "evalues" in hit_dict and isinstance(hit_dict["evalues"], list):
        hit.evalues = hit_dict["evalues"]

    # Set HSP count
    if "hsp_count" in hit_dict:
        try:
            hit.hsp_count = int(hit_dict["hsp_count"])
        except (ValueError, TypeError):
            hit.hsp_count = 1
    else:
        hit.hsp_count = 1  # Default to 1 if not specified

    # Set discontinuous flag if available
    if "discontinuous" in hit_dict:
        hit.discontinuous = hit_dict["discontinuous"]

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
    # Try both standard and alternative paths
    hhsearch_paths = [
        file_paths.get('hh_xml', {}).get('exists_at'),
        os.path.join(os.path.dirname(file_paths.get('hhr', {}).get('exists_at', '')),
                    f"{pdb_id}_{chain_id}.{reference}.hhsearch.xml")
    ]

    hhsearch_path = None
    for path in hhsearch_paths:
        if path and os.path.exists(path):
            hhsearch_path = path
            break

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
