"""
DEPRECATED MODULE

This module is deprecated and will be removed in a future version.
Please use the following replacements:

- element_to_dict() -> xml_core.element_to_dict()
- ensure_dict() -> xml_core.ensure_dict()
- process_xml_with_model() -> model_mapper.XmlModelMapper.xml_to_models()
- dict_to_xml() -> xml_core.dict_to_element()

Migration guide available at: docs/xml_migration.md
"""

import xml.etree.ElementTree as ET
from typing import Dict, Any, List, Optional, Union, TypeVar, Type, Callable
import logging

import warnings

T = TypeVar('T')

def element_to_dict(*args, **kwargs):
    warnings.warn(
        "element_to_dict() is deprecated. Use xml_core.element_to_dict() instead.",
        DeprecationWarning, stacklevel=2
    )
    from ecod.utils.xml_core import element_to_dict as new_func
    return new_func(*args, **kwargs)

def _ensure_hit_dict(hit):
    """Ensure hit is a dictionary with required fields"""
    if isinstance(hit, dict):
        # Already a dictionary
        return hit.copy()
    elif hasattr(hit, 'attrib'):
        # Convert XML Element to dictionary
        hit_dict = dict(hit.attrib)

        # Extract query region
        query_reg = hit.find("query_reg")
        if query_reg is not None and query_reg.text:
            hit_dict["range"] = query_reg.text.strip()
            hit_dict["range_parsed"] = self._parse_range(query_reg.text.strip())

        return hit_dict
    elif isinstance(hit, str):
        # Not expected - log error and return empty dict
        self.logger.error(f"Received string instead of hit object: {hit}")
        return {}
    else:
        # Not a valid hit
        self.logger.error(f"Received invalid hit type: {type(hit)}")
        return {}

def ensure_dict(*args, **kwargs):
    warnings.warn(
        "ensure_dict() is deprecated. Use xml_core.ensure_dict() instead.",
        DeprecationWarning, stacklevel=2
    )
    from ecod.utils.xml_core import ensure_dict as new_func
    return new_func(*args, **kwargs)

def ensure_list_of_dicts(items: List[Any]) -> List[Dict[str, Any]]:
    """Ensure list contains only dictionaries"""
    return [ensure_dict(item) for item in items]

def process_xml_with_model(
    xml_path: str,
    element_path: str,
    model_class: Type[T],
    from_xml_method: Callable[[ET.Element], T] = None
) -> List[T]:
    """
    Process XML file and convert elements to model instances
    
    Args:
        xml_path: Path to XML file
        element_path: XPath to elements
        model_class: Model class to instantiate
        from_xml_method: Method to convert element to model (default: model_class.from_xml)
        
    Returns:
        List of model instances
    """
    logger = logging.getLogger('ecod.utils.xml')
    
    if from_xml_method is None:
        if hasattr(model_class, 'from_xml'):
            from_xml_method = model_class.from_xml
        else:
            logger.error(f"No from_xml method provided and {model_class.__name__} has no from_xml method")
            return []
    
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        
        elements = root.findall(element_path)
        return [from_xml_method(element) for element in elements]
    except ET.ParseError as e:
        logger.error(f"Error parsing XML file {xml_path}: {e}")
        return []
    except Exception as e:
        logger.error(f"Error processing XML file {xml_path}: {e}")
        return []

# Add to ecod/utils/xml_utils.py

def extract_text_from_dict(d, key, default=""):
    """Extract text content from a dictionary field safely"""
    if key not in d:
        return default

    value = d[key]
    if isinstance(value, list) and len(value) > 0:
        # It's a list, try to get text from the first item
        if isinstance(value[0], dict) and '_text' in value[0]:
            return value[0]['_text']
        return str(value[0])
    elif isinstance(value, dict) and '_text' in value:
        # It's a dictionary with text
        return value['_text']
    else:
        # It's a direct value
        return str(value)

def dict_to_xml(data_dict, root_name):
    """Convert a dictionary to XML Element"""
    if isinstance(root_name, str):
        root = ET.Element(root_name)
    else:
        root = root_name  # Already an Element

    for key, value in data_dict.items():
        if key.startswith('@'):
            # This is an attribute
            root.set(key[1:], str(value))
        elif key == '_text':
            # This is text content
            root.text = str(value)
        elif isinstance(value, list):
            # This is a list of child elements
            for item in value:
                child_elem = ET.SubElement(root, key)
                if isinstance(item, dict):
                    dict_to_xml(item, child_elem)
                else:
                    child_elem.text = str(item)
        elif isinstance(value, dict):
            # This is a child element
            child_elem = ET.SubElement(root, key)
            dict_to_xml(value, child_elem)
        else:
            # This is a simple child element with text
            child_elem = ET.SubElement(root, key)
            child_elem.text = str(value)

    return root

def get_hits_from_dict(data_dict, hit_path):
    """Extract hits from dictionary using a path specification

    Args:
        data_dict: Dictionary containing hit data
        hit_path: List of keys to navigate to hit list
                  e.g. ['blast_summ_doc', 'chain_blast_run', 'hits', 'hit']

    Returns:
        List of hit dictionaries
    """
    current = data_dict
    for key in hit_path[:-1]:  # Navigate to container of hits
        if not current or key not in current:
            return []

        value = current[key]
        if isinstance(value, list):
            if not value:  # Empty list
                return []
            current = value[0]  # Take first element
        else:
            current = value

    # Get the hits from the final key
    hit_key = hit_path[-1]
    if not current or hit_key not in current:
        return []

    hits = current[hit_key]
    return hits if isinstance(hits, list) else [hits]
