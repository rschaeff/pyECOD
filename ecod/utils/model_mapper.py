# ecod/utils/model_mapper.py - New consolidated model conversion
import os
import logging
import xml.etree.ElementTree as ET
from typing import Dict, List, Any, Optional, Type, TypeVar, Generic, Callable

from ecod.models.pipeline import BlastHit, HHSearchHit, DomainSummaryModel
from ecod.utils.xml_core import parse_xml_file, element_to_dict

logger = logging.getLogger(__name__)

T = TypeVar('T')

class XmlModelMapper(Generic[T]):
    """Generic XML to model mapper"""
    
    def __init__(self, model_class: Type[T], element_path: str):
        """Initialize mapper
        
        Args:
            model_class: Model class to map to
            element_path: XPath to find elements
        """
        self.model_class = model_class
        self.element_path = element_path
        
    def xml_to_models(self, xml_path: str) -> List[T]:
        """Convert XML file to list of model instances"""
        tree = parse_xml_file(xml_path)
        if tree is None:
            return []
            
        root = tree.getroot()
        elements = root.findall(self.element_path)
        
        if not hasattr(self.model_class, 'from_xml'):
            logger.error(f"Model class {self.model_class.__name__} missing from_xml method")
            return []
            
        return [self.model_class.from_xml(element) for element in elements]

# Pre-configured mappers for common tasks
blast_hit_mapper = XmlModelMapper(BlastHit, ".//blast_run/hits/hit")
chain_blast_mapper = XmlModelMapper(BlastHit, ".//chain_blast_run/hits/hit")
hhsearch_mapper = XmlModelMapper(HHSearchHit, ".//hh_hit_list/hh_hit")

def load_blast_hits(xml_path: str, hit_type: str = "domain_blast") -> List[BlastHit]:
    """Load BLAST hits from XML file with correct mapper"""
    if not os.path.exists(xml_path):
        logger.warning(f"File not found: {xml_path}")
        return []
        
    mapper = chain_blast_mapper if hit_type == "chain_blast" else blast_hit_mapper
    hits = mapper.xml_to_models(xml_path)
    
    # Set hit_type on all hits
    for hit in hits:
        hit.hit_type = hit_type
        
    return hits

def load_hhsearch_hits(xml_path: str) -> List[HHSearchHit]:
    """Load HHSearch hits from XML file"""
    if not os.path.exists(xml_path):
        logger.warning(f"File not found: {xml_path}")
        return []
        
    return hhsearch_mapper.xml_to_models(xml_path)

def create_domain_summary(pdb_id: str, chain_id: str, ref_version: str,
                         paths: Dict[str, Dict[str, str]]) -> DomainSummaryModel:
    """Create domain summary from evidence files"""
    # Create base model
    summary = DomainSummaryModel(
        pdb_id=pdb_id,
        chain_id=chain_id,
        reference=ref_version
    )
    
    # Process evidence files
    # ... [rest of implementation] ...
    
    return summary
