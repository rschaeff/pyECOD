# ecod/utils/xml_core.py - New consolidated module
import xml.etree.ElementTree as ET
import logging
from typing import Dict, List, Any, Optional, Union, TypeVar, Type, Callable

logger = logging.getLogger(__name__)

# Core XML parsing functions
def parse_xml_file(file_path: str) -> Optional[ET.ElementTree]:
    """Parse XML file with robust error handling"""
    try:
        return ET.parse(file_path)
    except ET.ParseError as e:
        logger.error(f"XML parsing error for {file_path}: {str(e)}")
        return None
    except Exception as e:
        logger.error(f"Error reading XML file {file_path}: {str(e)}")
        return None

def element_to_dict(element: ET.Element, include_children: bool = True) -> Dict[str, Any]:
    """Convert XML Element to dictionary with consistent approach"""
    if element is None:
        return {}
        
    result = dict(element.attrib)
    
    # Include text content if present
    if element.text and element.text.strip():
        result["_text"] = element.text.strip()
    
    # Include child elements if requested
    if include_children:
        for child in element:
            # Handle repeated elements (create a list)
            if child.tag in result:
                if not isinstance(result[child.tag], list):
                    result[child.tag] = [result[child.tag]]
                result[child.tag].append(element_to_dict(child))
            else:
                # First occurrence of this tag
                result[child.tag] = element_to_dict(child)
    
    return result

def dict_to_element(data: Dict[str, Any], tag_name: str) -> ET.Element:
    """Convert dictionary to XML Element with consistent approach"""
    element = ET.Element(tag_name)
    
    # Set attributes
    for key, value in data.items():
        if key.startswith('@'):
            # This is an attribute
            element.set(key[1:], str(value))
        elif key == '_text':
            # This is element text
            element.text = str(value)
        elif isinstance(value, dict):
            # This is a child element
            child = dict_to_element(value, key)
            element.append(child)
        elif isinstance(value, list):
            # This is a list of child elements
            for item in value:
                if isinstance(item, dict):
                    child = dict_to_element(item, key)
                    element.append(child)
                else:
                    child = ET.SubElement(element, key)
                    child.text = str(item)
        else:
            # This is a simple child with text
            child = ET.SubElement(element, key)
            child.text = str(value)
    
    return element

def element_to_pretty_string(element: ET.Element) -> str:
    """Convert Element to pretty-formatted XML string"""
    from xml.dom import minidom
    rough_string = ET.tostring(element, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")
