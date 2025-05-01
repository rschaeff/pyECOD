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
