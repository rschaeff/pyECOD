# ecod/utils/model_conversion.py - Improved Implementation

import os
import logging
from typing import List, Dict, Any, Optional, Union
import xml.etree.ElementTree as ET
from ecod.models.pipeline import BlastHit, HHSearchHit, DomainSummaryModel

def create_domain_summary(pdb_id: str, chain_id: str, ref_version: str,
                         paths: Dict[str, Dict[str, str]]) -> DomainSummaryModel:
    """
    Create a domain summary model from evidence files

    Args:
        pdb_id: PDB ID
        chain_id: Chain ID
        ref_version: Reference version
        paths: Dictionary with file paths

    Returns:
        DomainSummaryModel with all evidence integrated
    """
    logger = logging.getLogger("ecod.model_conversion")

    # Create base model
    summary = DomainSummaryModel(
        pdb_id=pdb_id,
        chain_id=chain_id,
        reference=ref_version
    )

    logger.info(f"Creating domain summary for {pdb_id}_{chain_id} with reference {ref_version}")

    # Process FASTA for sequence length
    fasta_path = paths.get('fasta', {}).get('exists_at')
    if fasta_path and os.path.exists(fasta_path):
        try:
            with open(fasta_path, 'r') as f:
                lines = f.readlines()
            sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
            summary.sequence = sequence
            summary.sequence_length = len(sequence)
            logger.info(f"Found sequence of length {summary.sequence_length}")
        except Exception as e:
            logger.error(f"Error reading FASTA file: {e}")
            summary.errors["fasta_error"] = True

    # Process chain BLAST
    chain_blast_path = paths.get('chain_blast', {}).get('exists_at')
    if chain_blast_path and os.path.exists(chain_blast_path):
        try:
            # Use the robust BLAST XML parsing function from debug_domain_summary.py approach
            chain_hits = process_blast_xml(chain_blast_path, "chain_blast")
            summary.chain_blast_hits = chain_hits
            logger.info(f"Added {len(chain_hits)} chain BLAST hits")
        except Exception as e:
            logger.error(f"Error processing chain BLAST: {e}")
            summary.errors["chain_blast_error"] = True
    else:
        logger.warning(f"No chain BLAST file found")
        summary.errors["no_chain_blast"] = True

    # Process domain BLAST
    domain_blast_path = paths.get('domain_blast', {}).get('exists_at')
    if domain_blast_path and os.path.exists(domain_blast_path):
        try:
            domain_hits = process_blast_xml(domain_blast_path, "domain_blast")
            summary.domain_blast_hits = domain_hits
            logger.info(f"Added {len(domain_hits)} domain BLAST hits")
        except Exception as e:
            logger.error(f"Error processing domain BLAST: {e}")
            summary.errors["domain_blast_error"] = True
    else:
        logger.warning(f"No domain BLAST file found")
        summary.errors["no_domain_blast"] = True

    # Process HHSearch results
    # First try .xml standard path, then try other possible paths
    hh_xml_path = paths.get('hh_xml', {}).get('exists_at')
    if hh_xml_path and os.path.exists(hh_xml_path):
        try:
            hh_hits = process_hhsearch_xml(hh_xml_path)
            summary.hhsearch_hits = hh_hits
            logger.info(f"Added {len(hh_hits)} HHSearch hits")
        except Exception as e:
            logger.error(f"Error processing HHSearch XML: {e}")
            summary.errors["hhsearch_error"] = True
    else:
        logger.warning(f"No HHSearch XML file found")
        summary.errors["no_hhsearch"] = True

    # Set output file path and mark as processed
    summary.output_file_path = paths.get('domain_summary', {}).get('standard_path')
    summary.processed = True

    # Log summary
    logger.info(f"Domain summary created with {len(summary.chain_blast_hits)} chain BLAST hits, "
               f"{len(summary.domain_blast_hits)} domain BLAST hits, and "
               f"{len(summary.hhsearch_hits)} HHSearch hits")

    return summary


def process_blast_xml(xml_path: str, hit_type: str) -> List[BlastHit]:
    """
    Process BLAST XML and convert to list of BlastHit objects

    Based on debug_domain_summary.py approach which successfully parses various XML formats

    Args:
        xml_path: Path to BLAST XML file
        hit_type: Type of BLAST hits ("chain_blast" or "domain_blast")

    Returns:
        List of BlastHit objects
    """
    logger = logging.getLogger("ecod.model_conversion")
    logger.info(f"Processing {hit_type} XML file: {xml_path}")

    try:
        # First, analyze the XML structure to determine format
        tree = ET.parse(xml_path)
        root = tree.getroot()

        hits = []
        hit_elements = []

        # Try different formats based on root structure
        if root.tag == "BlastOutput" or root.find(".//BlastOutput_iterations") is not None:
            # Standard NCBI BLAST output format
            hit_elements = root.findall(".//BlastOutput_iterations/Iteration/Iteration_hits/Hit")
        else:
            # Try domain summary format structure
            if hit_type == "chain_blast":
                # Look for chainwise blast runs
                hit_elements = root.findall(".//chain_blast_run/hits/hit")
                if not hit_elements:
                    hit_elements = root.findall(".//chain_blast_run//hit")
            else:
                # Look for domain blast runs
                hit_elements = root.findall(".//blast_run/hits/hit")
                if not hit_elements:
                    hit_elements = root.findall(".//blast_run//hit")

            # Last resort - try to find any hit elements
            if not hit_elements:
                hit_elements = root.findall(".//hit")

        logger.info(f"Found {len(hit_elements)} hit elements in {hit_type} BLAST XML")

        # Process each hit
        for i, hit_elem in enumerate(hit_elements):
            try:
                blast_hit = BlastHit(hit_type=hit_type)

                # Handle standard NCBI BLAST format
                if root.tag == "BlastOutput" or root.find(".//BlastOutput_iterations") is not None:
                    # Get Hit_num
                    hit_num = hit_elem.find("Hit_num")
                    if hit_num is not None and hit_num.text:
                        blast_hit.hit_id = hit_num.text
                    else:
                        blast_hit.hit_id = str(i+1)

                    # Get Hit_def for PDB ID and chain ID
                    hit_def = hit_elem.find("Hit_def")
                    if hit_def is not None and hit_def.text:
                        hit_def_text = hit_def.text.strip()

                        if hit_type == "chain_blast":
                            # For chain BLAST, format is typically "2lxo A"
                            parts = hit_def_text.split()
                            if len(parts) >= 2:
                                blast_hit.pdb_id = parts[0].lower()
                                blast_hit.chain_id = parts[1]
                        else:
                            # For domain BLAST, format is typically "e2lxoA1 A:1-44 001088464"
                            parts = hit_def_text.split()
                            if len(parts) >= 1:
                                blast_hit.domain_id = parts[0]
                                # Try to extract PDB ID from domain ID
                                if len(parts[0]) >= 5 and parts[0][0] in "edgx":
                                    blast_hit.pdb_id = parts[0][1:5].lower()
                                    if len(parts[0]) > 5:
                                        blast_hit.chain_id = parts[0][5]

                    # Get HSPs
                    hsps = hit_elem.findall("Hit_hsps/Hsp")
                    if hsps:
                        blast_hit.hsp_count = len(hsps)

                        # Get E-value from first HSP
                        evalue_elem = hsps[0].find("Hsp_evalue")
                        if evalue_elem is not None and evalue_elem.text:
                            try:
                                blast_hit.evalue = float(evalue_elem.text)
                            except ValueError:
                                blast_hit.evalue = 999.0

                        # Get query range from first HSP
                        q_from = hsps[0].find("Hsp_query-from")
                        q_to = hsps[0].find("Hsp_query-to")
                        if q_from is not None and q_to is not None and q_from.text and q_to.text:
                            blast_hit.range = f"{q_from.text}-{q_to.text}"

                        # Get hit range from first HSP
                        h_from = hsps[0].find("Hsp_hit-from")
                        h_to = hsps[0].find("Hsp_hit-to")
                        if h_from is not None and h_to is not None and h_from.text and h_to.text:
                            blast_hit.hit_range = f"{h_from.text}-{h_to.text}"

                        # Check if discontinuous
                        if len(hsps) > 1:
                            ranges = []
                            for hsp in hsps:
                                q_from = hsp.find("Hsp_query-from")
                                q_to = hsp.find("Hsp_query-to")
                                if q_from is not None and q_to is not None and q_from.text and q_to.text:
                                    ranges.append((int(q_from.text), int(q_to.text)))

                            # Sort and check for discontinuity
                            if ranges:
                                ranges.sort(key=lambda x: x[0])
                                discontinuous = False
                                for j in range(1, len(ranges)):
                                    if ranges[j][0] > ranges[j-1][1] + 1:
                                        discontinuous = True
                                        break

                                blast_hit.discontinuous = discontinuous

                                # If discontinuous, create comma-separated range
                                if discontinuous:
                                    blast_hit.range = ",".join(f"{start}-{end}" for start, end in ranges)

                # Handle alternative XML formats
                else:
                    # Extract attributes
                    for attr, value in hit_elem.attrib.items():
                        if attr == "num":
                            blast_hit.hit_id = value
                        elif attr == "domain_id":
                            blast_hit.domain_id = value
                        elif attr == "pdb_id":
                            blast_hit.pdb_id = value
                        elif attr == "chain_id":
                            blast_hit.chain_id = value
                        elif attr == "evalues":
                            try:
                                if "," in value:
                                    evalues = [float(e) for e in value.split(",")]
                                    blast_hit.evalue = min(evalues)
                                    blast_hit.evalues = evalues
                                else:
                                    blast_hit.evalue = float(value)
                            except ValueError:
                                blast_hit.evalue = 999.0
                        elif attr == "hsp_count":
                            try:
                                blast_hit.hsp_count = int(value)
                            except ValueError:
                                blast_hit.hsp_count = 1
                        elif attr == "discontinuous":
                            blast_hit.discontinuous = (value.lower() == "true")

                    # Set hit_id if not already set
                    if not blast_hit.hit_id:
                        blast_hit.hit_id = str(i+1)

                    # Extract query region
                    query_reg = hit_elem.find("query_reg")
                    if query_reg is not None and query_reg.text:
                        blast_hit.range = query_reg.text.strip()

                    # Extract hit region
                    hit_reg = hit_elem.find("hit_reg")
                    if hit_reg is not None and hit_reg.text:
                        blast_hit.hit_range = hit_reg.text.strip()

                # Parse the ranges
                blast_hit.parse_ranges()

                # Add to hits list
                hits.append(blast_hit)

            except Exception as e:
                logger.warning(f"Error processing BLAST hit {i+1}: {str(e)}")
                continue

        logger.info(f"Successfully processed {len(hits)} hits for {hit_type} BLAST")
        return hits

    except Exception as e:
        logger.error(f"Error parsing {hit_type} BLAST XML: {str(e)}")
        return []


def process_hhsearch_xml(xml_path: str) -> List[HHSearchHit]:
    """
    Process HHSearch XML and convert to list of HHSearchHit objects

    Based on debug_domain_summary.py approach which successfully parses various XML formats

    Args:
        xml_path: Path to HHSearch XML file

    Returns:
        List of HHSearchHit objects
    """
    logger = logging.getLogger("ecod.model_conversion")
    logger.info(f"Processing HHSearch XML file: {xml_path}")

    try:
        # Parse the XML
        tree = ET.parse(xml_path)
        root = tree.getroot()
        hits = []

        # First try finding hh_hit elements in hh_hit_list
        hit_elements = root.findall(".//hh_hit_list/hh_hit")

        if not hit_elements:
            # Try alternative paths
            hit_elements = root.findall(".//hits/hit")

        if not hit_elements:
            # Last resort
            hit_elements = root.findall(".//hit")

        logger.info(f"Found {len(hit_elements)} HHSearch hit elements")

        for i, hit_elem in enumerate(hit_elements):
            try:
                hh_hit = HHSearchHit()

                # Extract attributes
                for attr, value in hit_elem.attrib.items():
                    if attr == "hit_id" or attr == "id":
                        hh_hit.hit_id = value
                    elif attr == "domain_id" or attr == "ecod_domain_id":
                        hh_hit.domain_id = value
                    elif attr == "probability":
                        try:
                            hh_hit.probability = float(value)
                        except ValueError:
                            hh_hit.probability = 0.0
                    elif attr == "evalue" or attr == "e_value":
                        try:
                            hh_hit.evalue = float(value)
                        except ValueError:
                            hh_hit.evalue = 999.0
                    elif attr == "score":
                        try:
                            hh_hit.score = float(value)
                        except ValueError:
                            hh_hit.score = 0.0

                # Set default hit_id if not set
                if not hh_hit.hit_id:
                    hh_hit.hit_id = f"hit_{i+1}"

                # Set domain_id from hit_id if not set
                if not hh_hit.domain_id and hh_hit.hit_id:
                    hh_hit.domain_id = hh_hit.hit_id

                # Extract query range
                query_range = hit_elem.find(".//query_range")
                if query_range is not None and query_range.text:
                    hh_hit.range = query_range.text.strip()
                else:
                    # Try alternate format
                    query_from = hit_elem.find(".//query_from")
                    query_to = hit_elem.find(".//query_to")
                    if query_from is not None and query_to is not None and query_from.text and query_to.text:
                        hh_hit.range = f"{query_from.text}-{query_to.text}"

                # Extract hit range
                hit_range = hit_elem.find(".//hit_range") or hit_elem.find(".//template_seqid_range")
                if hit_range is not None and hit_range.text:
                    hh_hit.hit_range = hit_range.text.strip()

                # Parse the ranges
                hh_hit.parse_ranges()

                # Add to hits list
                hits.append(hh_hit)

            except Exception as e:
                logger.warning(f"Error processing HHSearch hit {i+1}: {str(e)}")
                continue

        logger.info(f"Successfully processed {len(hits)} HHSearch hits")
        return hits

    except Exception as e:
        logger.error(f"Error parsing HHSearch XML: {str(e)}")
        return []
