#ecod/utils/model_conversion.py

import os
import logging
from typing import List, Dict, Any, Optional, Union
import xml.etree.ElementTree as ET
from ecod.models.pipeline import BlastHit, HHSearchHit, DomainSummaryModel

def element_to_model(element: ET.Element, model_class):
    """Convert XML Element to model instance"""
    if hasattr(model_class, 'from_xml'):
        return model_class.from_xml(element)
    raise ValueError(f"Model class {model_class.__name__} doesn't have from_xml method")

def xml_file_to_models(xml_path: str, element_path: str, model_class) -> List:
    """Convert XML file to list of model instances"""
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        
        elements = root.findall(element_path)
        return [element_to_model(element, model_class) for element in elements]
    except Exception as e:
        logging.error(f"Error parsing XML file {xml_path}: {e}")
        return []

def xml_to_hits(xml_path: str, hit_type: str) -> List[Union[BlastHit, HHSearchHit]]:
    """Parse XML file to hits of specified type"""
    if hit_type == "hhsearch":
        return xml_file_to_models(xml_path, ".//hh_hit_list/hh_hit", HHSearchHit)
    elif hit_type == "chain_blast":
        return xml_file_to_models(xml_path, ".//chain_blast_run/hits/hit", BlastHit)
    elif hit_type == "domain_blast":
        return xml_file_to_models(xml_path, ".//blast_run/hits/hit", BlastHit)
    else:
        raise ValueError(f"Unknown hit type: {hit_type}")

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
    # First try .hhsearch.xml extension, then fall back to .xml
    hh_xml_path = paths.get('hh_xml', {}).get('exists_at')
    if not hh_xml_path:
        base_path = os.path.join(os.path.dirname(paths.get('hhr', {}).get('exists_at', '')),
                                f"{pdb_id}_{chain_id}.{ref_version}.hhsearch.xml")
        if os.path.exists(base_path):
            hh_xml_path = base_path

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
    Process standard NCBI BLAST XML output format and convert to BlastHit objects

    Args:
        xml_path: Path to BLAST XML file
        hit_type: Type of BLAST hits ("chain_blast" or "domain_blast")

    Returns:
        List of BlastHit objects
    """
    logger = logging.getLogger("ecod.model_conversion")
    logger.info(f"Processing {hit_type} XML file: {xml_path}")

    try:
        # Parse the input XML
        tree = ET.parse(xml_path)
        root = tree.getroot()
        hits = []

        # Find all Hit elements under Iteration_hits
        hit_elements = root.findall(".//Iteration_hits/Hit")
        logger.info(f"Found {len(hit_elements)} hit elements in {hit_type} BLAST XML")

        # Process each hit
        for hit_elem in hit_elements:
            try:
                hit_data = {}

                # Get hit number
                hit_num = hit_elem.find("./Hit_num")
                if hit_num is not None and hit_num.text:
                    hit_data["hit_id"] = hit_num.text

                # Get hit definition and parse PDB ID and chain ID
                hit_def = hit_elem.find("./Hit_def")
                if hit_def is not None and hit_def.text:
                    hit_def_text = hit_def.text.strip()
                    logger.debug(f"Hit definition: {hit_def_text}")

                    if hit_type == "chain_blast":
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
                            # Try to extract PDB ID from domain ID
                            if len(parts[0]) >= 5 and parts[0][0] in "edgx":
                                # Format like e2lxoA1 - pdb_id is characters 2-5
                                hit_data["pdb_id"] = parts[0][1:5].lower()
                                # Chain ID is typically the character after PDB ID
                                if len(parts[0]) > 5:
                                    hit_data["chain_id"] = parts[0][5]

                # Get HSPs (High-scoring Segment Pairs)
                hsps = hit_elem.findall("./Hit_hsps/Hsp")
                if hsps:
                    hit_data["hsp_count"] = len(hsps)

                    # Process the first HSP (best alignment)
                    hsp = hsps[0]

                    # Get E-value
                    evalue_elem = hsp.find("./Hsp_evalue")
                    if evalue_elem is not None and evalue_elem.text:
                        try:
                            hit_data["evalue"] = float(evalue_elem.text)
                        except ValueError:
                            hit_data["evalue"] = 999.0

                    # Get query range
                    q_from = hsp.find("./Hsp_query-from")
                    q_to = hsp.find("./Hsp_query-to")
                    if q_from is not None and q_to is not None and q_from.text and q_to.text:
                        hit_data["range"] = f"{q_from.text}-{q_to.text}"

                    # Get hit range
                    h_from = hsp.find("./Hsp_hit-from")
                    h_to = hsp.find("./Hsp_hit-to")
                    if h_from is not None and h_to is not None and h_from.text and h_to.text:
                        hit_data["hit_range"] = f"{h_from.text}-{h_to.text}"

                    # If multiple HSPs, check if discontinuous
                    if len(hsps) > 1:
                        ranges = []
                        for hsp in hsps:
                            q_from = hsp.find("./Hsp_query-from")
                            q_to = hsp.find("./Hsp_query-to")
                            if q_from is not None and q_to is not None and q_from.text and q_to.text:
                                ranges.append((int(q_from.text), int(q_to.text)))

                        # Sort and check for discontinuity
                        if ranges:
                            ranges.sort(key=lambda x: x[0])
                            discontinuous = False
                            for i in range(1, len(ranges)):
                                if ranges[i][0] > ranges[i-1][1] + 1:
                                    discontinuous = True
                                    break

                            hit_data["discontinuous"] = discontinuous

                            # If discontinuous, create comma-separated range
                            if discontinuous:
                                hit_data["range"] = ",".join(f"{start}-{end}" for start, end in ranges)

                # Set hit type
                hit_data["hit_type"] = hit_type

                # Create BlastHit object (fill in defaults for missing fields)
                hit_data.setdefault("pdb_id", "")
                hit_data.setdefault("chain_id", "")
                hit_data.setdefault("domain_id", "")
                hit_data.setdefault("range", "")
                hit_data.setdefault("hit_range", "")
                hit_data.setdefault("evalue", 999.0)
                hit_data.setdefault("hsp_count", 1)

                blast_hit = BlastHit(**hit_data)
                hits.append(blast_hit)

                logger.debug(f"Processed hit {hit_data['hit_id']}: {hit_data['pdb_id']}_{hit_data['chain_id']} range={hit_data['range']} evalue={hit_data['evalue']}")

            except Exception as e:
                logger.warning(f"Error processing BLAST hit: {str(e)}")
                continue

        logger.info(f"Successfully processed {len(hits)} of {len(hit_elements)} BLAST hits")
        return hits

    except Exception as e:
        logger.error(f"Error processing {hit_type} XML: {str(e)}", exc_info=True)
        return []

def process_hhsearch_xml(xml_path: str) -> List[HHSearchHit]:
    """
    Process HHSearch XML file to extract hits

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

        # Find all hh_hit elements in hh_hit_list
        hit_elements = root.findall(".//hh_hit_list/hh_hit")
        logger.info(f"Found {len(hit_elements)} hit elements in HHSearch XML")

        hits = []
        for hit_elem in hit_elements:
            try:
                # Extract hit data from attributes
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

                # Convert numeric fields to proper types
                if "probability" in hit_data:
                    try:
                        hit_data["probability"] = float(hit_data["probability"])
                    except ValueError:
                        hit_data["probability"] = 0.0

                if "e_value" in hit_data:
                    try:
                        hit_data["evalue"] = float(hit_data["e_value"])
                        # Remove the original e_value key to avoid duplication
                        del hit_data["e_value"]
                    except ValueError:
                        hit_data["evalue"] = 999.0

                if "score" in hit_data:
                    try:
                        hit_data["score"] = float(hit_data["score"])
                    except ValueError:
                        hit_data["score"] = 0.0

                # Create HHSearchHit object with properly typed values
                hit_data.setdefault("probability", 0.0)
                hit_data.setdefault("evalue", 999.0)
                hit_data.setdefault("score", 0.0)
                hit_data.setdefault("range", "")
                hit_data.setdefault("hit_range", "")

                # Ensure critical fields are present
                if "hit_id" not in hit_data and "ecod_domain_id" in hit_data:
                    hit_data["hit_id"] = hit_data["ecod_domain_id"]

                if "domain_id" not in hit_data and "ecod_domain_id" in hit_data:
                    hit_data["domain_id"] = hit_data["ecod_domain_id"]

                hh_hit = HHSearchHit(**hit_data)
                hits.append(hh_hit)

                logger.debug(f"Processed HHSearch hit: {hit_data['hit_id']}, probability={hit_data['probability']}, range={hit_data['range']}")

            except Exception as e:
                logger.warning(f"Error processing HHSearch hit: {str(e)}")
                continue

        logger.info(f"Successfully processed {len(hits)} HHSearch hits")
        return hits

    except Exception as e:
        logger.error(f"Error processing HHSearch XML: {str(e)}", exc_info=True)
        return []

def read_fasta_sequence(fasta_path: str) -> Optional[str]:
    """
    Read sequence from FASTA file

    Args:
        fasta_path: Path to FASTA file

    Returns:
        Sequence string or None if error
    """
    logger = logging.getLogger("ecod.model_conversion")

    try:
        with open(fasta_path, 'r') as f:
            lines = f.readlines()

        # Skip header line and concatenate sequence lines
        sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
        return sequence
    except Exception as e:
        logger.error(f"Error reading FASTA file {fasta_path}: {e}")
        return None
