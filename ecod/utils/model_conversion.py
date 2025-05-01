#ecod/utils/model_conversion.py

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

# Add to ecod/utils/model_conversion.py

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
    # Create base model
    summary = DomainSummaryModel(
        pdb_id=pdb_id,
        chain_id=chain_id,
        reference=ref_version
    )

    # Get sequence length if FASTA exists
    fasta_path = paths.get('fasta', {}).get('exists_at')
    if fasta_path and os.path.exists(fasta_path):
        sequence = read_fasta_sequence(fasta_path)
        summary.sequence_length = len(sequence) if sequence else 0

    # Process chain BLAST
    chain_blast_path = paths.get('chain_blast', {}).get('exists_at')
    if chain_blast_path and os.path.exists(chain_blast_path):
        try:
            summary.chain_blast_hits = process_blast_xml(chain_blast_path, "chain_blast")
        except Exception as e:
            logging.error(f"Error processing chain BLAST: {e}")
            summary.errors["chain_blast_error"] = True
    else:
        summary.errors["no_chain_blast"] = True

    # Process domain BLAST
    domain_blast_path = paths.get('domain_blast', {}).get('exists_at')
    if domain_blast_path and os.path.exists(domain_blast_path):
        try:
            summary.domain_blast_hits = process_blast_xml(domain_blast_path, "domain_blast")
        except Exception as e:
            logging.error(f"Error processing domain BLAST: {e}")
            summary.errors["domain_blast_error"] = True
    else:
        summary.errors["no_domain_blast"] = True

    # Process HHSearch
    hhsearch_path = paths.get('hh_xml', {}).get('exists_at')
    if hhsearch_path and os.path.exists(hhsearch_path):
        try:
            summary.hhsearch_hits = process_hhsearch_xml(hhsearch_path)
        except Exception as e:
            logging.error(f"Error processing HHSearch: {e}")
            summary.errors["hhsearch_error"] = True
    else:
        summary.errors["no_hhsearch"] = True

    return summary
