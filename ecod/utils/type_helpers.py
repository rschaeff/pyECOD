# ecod/utils/type_helpers.py
from typing import Dict, Any, List, Optional
import xml.etree.ElementTree as ET

def ensure_dict(obj: Any) -> Dict[str, Any]:
    """Ensure object is a dictionary
    
    If obj is an XML Element, converts its attributes to a dictionary.
    Otherwise, returns obj if it's already a dictionary, or an empty dict.
    
    Args:
        obj: Object to convert
        
    Returns:
        Dictionary version of object
    """
    if hasattr(obj, 'attrib'):  # It's an XML Element
        result = dict(obj.attrib)
        return result
    elif isinstance(obj, dict):
        return obj
    else:
        return {}  # Return empty dict for incompatible types

def ensure_list_of_dicts(items: List[Any]) -> List[Dict[str, Any]]:
    """Ensure list contains only dictionaries
    
    Args:
        items: List of items
        
    Returns:
        List with all items converted to dictionaries
    """
    return [ensure_dict(item) for item in items]

def element_to_dict(element: ET.Element, text_field: Optional[str] = None) -> Dict[str, Any]:
    """Convert XML Element to dictionary
    
    Args:
        element: XML Element to convert
        text_field: Field name for element's text (if any)
        
    Returns:
        Dictionary with element's attributes and text
    """
    result = dict(element.attrib)
    
    # Add text content if specified
    if text_field and element.text:
        result[text_field] = element.text.strip()
    
    # Add child elements
    for child in element:
        if child.tag not in result:
            # Get text content
            if child.text and child.text.strip():
                result[child.tag] = child.text.strip()
    
    return result