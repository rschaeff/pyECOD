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
        DomainSummaryModel with evidence
    """
    logger = logging.getLogger("ecod.model_conversion")

    # Create base model
    summary = DomainSummaryModel(
        pdb_id=pdb_id,
        chain_id=chain_id,
        reference=ref_version
    )
    logger.info(f"Starting domain summary for {pdb_id} {chain_id} {ref_version}")

    # Get sequence length if FASTA exists
    fasta_path = paths.get('fasta', {}).get('exists_at')
    if fasta_path and os.path.exists(fasta_path):
        logger.info(f"Processing FASTA file: {fasta_path}")
        sequence = read_fasta_sequence(fasta_path)
        summary.sequence_length = len(sequence) if sequence else 0
    else:
        logger.warning(f"No FASTA file found at {fasta_path}")

    # Process chain BLAST
    chain_blast_path = paths.get('chain_blast', {}).get('exists_at')
    if chain_blast_path and os.path.exists(chain_blast_path):
        try:
            logger.info(f"Processing chain BLAST file: {chain_blast_path}")
            summary.chain_blast_hits = process_blast_xml(chain_blast_path, "chain_blast")
            logger.info(f"Added {len(summary.chain_blast_hits)} chain BLAST hits")
        except Exception as e:
            logger.error(f"Error processing chain BLAST: {e}")
            summary.errors["chain_blast_error"] = True
    else:
        logger.warning(f"No chain BLAST file found at {chain_blast_path}")
        summary.errors["no_chain_blast"] = True

    # Process domain BLAST
    domain_blast_path = paths.get('domain_blast', {}).get('exists_at')
    if domain_blast_path and os.path.exists(domain_blast_path):
        try:
            logger.info(f"Processing domain BLAST file: {domain_blast_path}")
            summary.domain_blast_hits = process_blast_xml(domain_blast_path, "domain_blast")
            logger.info(f"Added {len(summary.domain_blast_hits)} domain BLAST hits")
        except Exception as e:
            logger.error(f"Error processing domain BLAST: {e}")
            summary.errors["domain_blast_error"] = True
    else:
        logger.warning(f"No domain BLAST file found at {domain_blast_path}")
        summary.errors["no_domain_blast"] = True

    # Process HHSearch
    hhsearch_path = paths.get('hh_xml', {}).get('exists_at')
    if hhsearch_path and os.path.exists(hhsearch_path):
        try:
            logger.info(f"Processing HHSearch file: {hhsearch_path}")
            summary.hhsearch_hits = process_hhsearch_xml(hhsearch_path)
            logger.info(f"Added {len(summary.hhsearch_hits)} HHSearch hits")
        except Exception as e:
            logger.error(f"Error processing HHSearch: {e}")
            summary.errors["hhsearch_error"] = True
    else:
        logger.warning(f"No HHSearch file found at {hhsearch_path}")
        summary.errors["no_hhsearch"] = True

    logger.info(f"Completed domain summary creation for {pdb_id}_{chain_id}")
    return summary

    def process_blast_xml(xml_path: str, hit_type: str) -> List[BlastHit]:
    """
    Process BLAST XML file to extract hits

    Args:
        xml_path: Path to BLAST XML file
        hit_type: Type of BLAST hits (chain_blast or domain_blast)

    Returns:
        List of BlastHit objects
    """
    logger = logging.getLogger("ecod.model_conversion")
    logger.info(f"Processing {hit_type} XML file: {xml_path}")

    try:
        return xml_to_hits(xml_path, hit_type)
    except Exception as e:
        logger.error(f"Error processing {hit_type} XML: {e}")
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
        return xml_to_hits(xml_path, "hhsearch")
    except Exception as e:
        logger.error(f"Error processing HHSearch XML: {e}")
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
