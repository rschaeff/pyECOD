# Add to ecod/models/base.py

from typing import Dict, Any, ClassVar, Type
import xml.etree.ElementTree as ET

class XmlSerializable:
    """Interface for models that can be serialized to/from XML"""
    
    # Class variable defining XML element path for finding elements
    xml_element_path: ClassVar[str] = ""
    
    @classmethod
    def from_xml(cls, element: ET.Element) -> 'XmlSerializable':
        """Create instance from XML Element"""
        raise NotImplementedError("Subclasses must implement from_xml")
    
    def to_xml(self) -> ET.Element:
        """Convert to XML Element"""
        raise NotImplementedError("Subclasses must implement to_xml")
    
    @classmethod
    def from_xml_file(cls, file_path: str) -> 'XmlSerializable':
        """Create instance from XML file"""
        from ecod.utils.xml_core import parse_xml_file
        
        tree = parse_xml_file(file_path)
        if tree is None:
            return None
            
        root = tree.getroot()
        return cls.from_xml(root)
    
    def to_xml_file(self, file_path: str) -> bool:
        """Save to XML file"""
        from ecod.utils.xml_core import element_to_pretty_string
        
        element = self.to_xml()
        if element is None:
            return False
            
        try:
            xml_string = element_to_pretty_string(element)
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(xml_string)
            return True
        except Exception as e:
            import logging
            logging.getLogger(__name__).error(f"Error saving XML: {str(e)}")
            return False
