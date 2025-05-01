# tests/utils/test_xml_core.py

import unittest
import xml.etree.ElementTree as ET
from ecod.utils.xml_core import (
    element_to_dict, dict_to_element, element_to_pretty_string
)

class TestXmlCore(unittest.TestCase):
    
    def test_element_to_dict_simple(self):
        """Test conversion of simple element to dict"""
        xml = '<hit num="1" evalue="0.001"><query_reg>1-100</query_reg></hit>'
        element = ET.fromstring(xml)
        result = element_to_dict(element)
        
        self.assertEqual(result['@num'], '1')
        self.assertEqual(result['@evalue'], '0.001')
        self.assertEqual(result['query_reg']['_text'], '1-100')
    
    def test_dict_to_element_simple(self):
        """Test conversion of dict to element"""
        data = {
            '@num': '1',
            '@evalue': '0.001',
            'query_reg': {'_text': '1-100'}
        }
        element = dict_to_element(data, 'hit')
        
        self.assertEqual(element.tag, 'hit')
        self.assertEqual(element.get('num'), '1')
        self.assertEqual(element.get('evalue'), '0.001')
        self.assertEqual(element.find('query_reg').text, '1-100')
    
    # More test cases...
