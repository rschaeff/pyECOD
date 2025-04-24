# In ecod/utils/xml_utils.py

import xml.etree.ElementTree as ET
from typing import Dict, Any, List, Optional, Union, TypeVar, Type, Callable
import logging

T = TypeVar('T')

def element_to_dict(element: ET.Element, text_field: Optional[str] = None) -> Dict[str, Any]:
    """Convert XML Element to dictionary"""
    result = dict(element.attrib)
    
    # Add text content if specified
    if text_field and element.text and element.text.strip():
        result[text_field] = element.text.strip()
    
    # Add child elements as fields
    for child in element:
        # Get text content from child
        if child.text and child.text.strip():
            result[child.tag] = child.text.strip()
    
    return result

def ensure_dict(obj: Any) -> Dict[str, Any]:
    """Ensure object is a dictionary"""
    if hasattr(obj, 'attrib'):  # It's an XML Element
        return element_to_dict(obj)
    elif isinstance(obj, dict):
        return obj
    else:
        return {}

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