#!/usr/bin/env python3
"""
Unit tests for the HHSearch parser and processor.
"""

import unittest
import os
import tempfile
import logging
import shutil
from pathlib import Path

from ecod.utils.hhsearch_utils import HHRParser
from ecod.pipelines.hhsearch import HHRToXMLConverter

# Configure logging
logging.basicConfig(level=logging.ERROR)


class TestHHRParser(unittest.TestCase):
    """Tests for the HHR parser"""
    
    def setUp(self):
        """Set up test data and parser instance"""
        self.test_dir = tempfile.mkdtemp()
        self.logger = logging.getLogger("test_logger")
        self.parser = HHRParser(self.logger)
        
        # Create test HHR file with known content
        self.test_hhr_path = os.path.join(self.test_dir, "test.hhr")
        with open(self.test_hhr_path, 'w') as f:
            f.write(self.get_test_hhr_content())
        
        # Create test HHR file with no hits
        self.empty_hhr_path = os.path.join(self.test_dir, "empty.hhr")
        with open(self.empty_hhr_path, 'w') as f:
            f.write(self.get_empty_hhr_content())
        
        # Create test HHR file with malformed content
        self.malformed_hhr_path = os.path.join(self.test_dir, "malformed.hhr")
        with open(self.malformed_hhr_path, 'w') as f:
            f.write(self.get_malformed_hhr_content())
    
    def tearDown(self):
        """Clean up test files"""
        shutil.rmtree(self.test_dir)
    
    def get_test_hhr_content(self):
        """Return content for a test HHR file"""
        return """Query         test_protein
Match_columns 168
No_of_seqs    150 out of 1487
Neff          6.8
Searched_HMMs 44289
Date          Mon Apr 21 15:48:32 2025
Command       /sw/apps/hh-suite/bin/hhsearch -i test.hhm -d ecod_v291 -o test.hhr

 No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
  1 e6gcsh1 h:7-136 002393507       99.4 5.3E-16 1.2E-20  112.2   0.0  121    1-150     8-130 (130)
  2 e5lnki1 i:1-144 001930624       99.3 2.5E-15 5.5E-20  108.6   0.0  121    1-151    17-137 (137)

>e6gcsh1 h:7-136 002393507
Probab=99.44  E-value=5.3e-16  Score=112.22  Aligned_cols=121  Identities=33%  Similarity=0.509  Sum_probs=45.4  Template_Neff=7.400

Q test_protein         1 MSWWSGVWRSVWSALSREVREHVGTDHLGNKYYYVAEYKNWRGQTIREKRIVEAANRKEVDYEAGDIPTEWEAWIRRTRK   80 (168)
Q Consensus            1 ms~~~~~~~~~~~~~~~~~~~~VG~D~~GN~Yye~~~~~~~~g~~~~~rR~ve~~~~~~~~y~~~~iP~EW~aWLr~~R~   80 (168)
                         +.|++++|++|+++++.+.++|||+|.+||+||+...... .   .+++|||+|...   .+++++||+|||+|||++++
T Consensus            8 ~~g~~~~~~~l~~~~~~~~~~lvG~D~~GN~Yye~~~~~~-~---~~~~R~V~y~~~---~~d~~~IP~eW~~WL~~~r~   80 (130)
T e6gcsh1_consen       8 RNGLRGSWRQLNSLGDWRKGKLVGTDQFGNKYYENPKPNS-F---GRRRRWVEYAGK---DYDASQIPPEWHAWLHHTRD   80 (130)
T ss_pred               hcCHHHHHHHHHccCcccCCeEEEEcCCCCeeeeCCCCCC-C---CCCeeEEEcCCC---CCChhhCChhHhhhhcccCC
Confidence              3334444444443333334444444444444444331000 1   123444444331   23444444444444444444


Q test_protein        81 TPPTMEEILKNEKY--REEIKIKSQDFYEKDKLGKETSEELLPSPTATQVKGHASAPYFGREEPSVAPTSTG  150 (168)
Q Consensus           81 ~PPT~eEi~~n~~~--~~~~~~~a~~le~k~~~~~~~~eg~~~~~~~~~~~g~~s~~~~~~~~~s~~p~~~g  150 (168)
                         +||++++++.....  .++|.+|++|+.+      +|+      |+     + +++++|+    +|+|.++.
T Consensus           81 ~pPt~~e~~~~~~~~~~~~h~~n~tgt~~------~y~------p~-----~-t~~~k~~----~W~p~~~~  130 (130)
T e6gcsh1_consen      81 DPPTEEEIMKYPRRKWQKPHKPNLTGTPG------AYV------PY-----S-TTRPKIQ----AWEPPVKQ  130 (130)
T ss_pred               CCCCHHHHhhCCCCCCCCCCCcCCCCChh------hhc------cC-----C-CCCCcee----cCCCCCCC
Confidence              44444433211000  1234444444444      444      44     3 3344444    44444333


>e5lnki1 i:1-144 001930624
Probab=99.35  E-value=2.5e-15  Score=108.57  Aligned_cols=121  Identities=34%  Similarity=0.542  Sum_probs=29.5  Template_Neff=7.800

Q test_protein         1 MSWWSGVWRSVWSALSREVREHVGTDHLGNKYYYVAEYKNWRGQTIREKRIVEAANRKEVDYEAGDIPTEWEAWIRRTRK   80 (168)
Q Consensus            1 ms~~~~~~~~~~~~~~~~~~~~VG~D~~GN~Yye~~~~~~~~g~~~~~rR~ve~~~~~~~~y~~~~iP~EW~aWLr~~R~   80 (168)
                         ++|++++|+.|++.+..+.+++||+|.+||+||+...  . .   .+++|||+|...  +.|++++||+|||+|||++++
T Consensus           17 ~~~~~~~~~~l~~~~~~~~~~~vG~D~~GN~yy~~~~--~-~---~g~~R~V~y~~~--~~~~~~~Ip~eW~~WL~~~r~   88 (137)
T e5lnki1_consen      17 MGGLRGSLRQLWRIDDWRKGKLVGTDQFGNKYYENPK--E-F---PGRRRWVEYAGK--WDYDASQIPPEWHAWLHHTRD   88 (137)
T ss_pred               cCcHHHHHHHHHhcCCcCCCeEEEEeCCCCeEEeCCC--C-C---CCCcceEEeCCC--cccccccCCHHHHHhhcccCC
Confidence              1222222222222222222222222222222222221  0 1   112222222210  012222222222222222222


Q test_protein        81 TPPTMEEILKNEKYREEIKIKSQDFYEKDKLGKETSEELLPSPTATQVKGHASAPYFGREEPSVAPTSTGK  151 (168)
Q Consensus           81 ~PPT~eEi~~n~~~~~~~~~~a~~le~k~~~~~~~~eg~~~~~~~~~~~g~~s~~~~~~~~~s~~p~~~g~  151 (168)
                         +||+.++++........|.+|++|+.+      +|+      |+     + ++++++.    +|+|..+.-
T Consensus           89 ~pPt~~e~~~~~~~~~~~~~n~tgt~~------~y~------p~-----~-~~~~~~~----~w~p~~~~~  137 (137)
T e5lnki1_consen      89 DPPTEEELLKRKWYWQPHKPNLTGTPG------AYV------PY-----S-TTRPKYE----AWQPPSAPY  137 (137)
T ss_pred               CCCCHHHHHHcccccCCCCCCCCCCch------hhc------CC-----C-CCCCcee----cCCCCCCCC
Confidence              222222221111111122222222222      222      22     2 2222222    222222221

"""
    
    def get_empty_hhr_content(self):
        """Return content for a test HHR file with no hits"""
        return """Query         test_protein
Match_columns 168
No_of_seqs    150 out of 1487
Neff          6.8
Searched_HMMs 44289
Date          Mon Apr 21 15:48:32 2025
Command       /sw/apps/hh-suite/bin/hhsearch -i test.hhm -d ecod_v291 -o test.hhr

 No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM

No 1
>
No 2
>
"""
    
    def get_malformed_hhr_content(self):
        """Return content for a test HHR file with malformed content"""
        return """Query         test_protein
Match_columns not-a-number
No_of_seqs    text out of more-text
Date          Mon Apr 21 15:48:32 2025

 No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
  1 e6gcsh1 h:7-136 002393507       99.4 invalid 1.2E-20  112.2   0.0  121    1-150     8-130 (130)

>e6gcsh1 h:7-136 002393507
Probab=99.44  E-value=5.3e-16  Score=112.22  Aligned_cols=121  Identities=33%  Similarity=0.509  Sum_probs=45.4  Template_Neff=7.400

Q test_protein         1 SEQUENCE-WITH-INVALID-POSITIONS   XX (168)
Q Consensus            1 INVALID-CONSENSUS-SEQUENCE   XX (168)
                         |||||||||||||||||||||
T Consensus            8 INVALID-TARGET-CONSENSUS   XX (130)
T e6gcsh1_consen       8 TARGET-SEQUENCE-WITH-ERRORS   XX (130)
"""
    
    def test_parse_header(self):
        """Test header parsing functionality"""
        # Parse the test HHR file
        result = self.parser.parse(self.test_hhr_path)
        self.assertIsNotNone(result, "Parser should not return None for valid HHR file")
        
        # Check header values
        header = result['header']
        self.assertEqual(header['query_id'], 'test_protein', "Query ID not correctly parsed")
        self.assertEqual(header['match_columns'], 168, "Match columns not correctly parsed")
        self.assertEqual(header['no_of_seqs'], 150, "Number of sequences not correctly parsed")
        self.assertEqual(header['total_seqs'], 1487, "Total sequences not correctly parsed")
        self.assertEqual(header['neff'], 6.8, "Neff value not correctly parsed")
        self.assertEqual(header['searched_hmms'], 44289, "Number of searched HMMs not correctly parsed")
        self.assertIn('date', header, "Date should be parsed")
        self.assertIn('command', header, "Command should be parsed")
    
    def test_parse_hit_summary(self):
        """Test hit summary parsing functionality"""
        # Parse the test HHR file
        result = self.parser.parse(self.test_hhr_path)
        self.assertIsNotNone(result, "Parser should not return None for valid HHR file")
        
        # Check number of hits
        hits = result['hits']
        self.assertEqual(len(hits), 2, "Parser should find 2 hits in the test file")
        
        # Check first hit details
        hit1 = hits[0]
        self.assertEqual(hit1['hit_num'], 1, "Hit number not correctly parsed")
        self.assertEqual(hit1['hit_id'], 'e6gcsh1', "Hit ID not correctly parsed")
        self.assertEqual(hit1['probability'], 99.4, "Probability not correctly parsed")
        self.assertEqual(hit1['e_value'], 5.3e-16, "E-value not correctly parsed")
        self.assertEqual(hit1['p_value'], 1.2e-20, "P-value not correctly parsed")
        self.assertEqual(hit1['score'], 112.2, "Score not correctly parsed")
        self.assertEqual(hit1['cols'], "121", "Columns not correctly parsed")
        
        # Check template range parsing from description
        self.assertIn('template_range', hit1, "Template range should be parsed")
        self.assertIn('h:7-136', hit1['description'], "Description should contain template range")
    
    def test_parse_alignments(self):
        """Test alignment parsing functionality"""
        # Parse the test HHR file
        result = self.parser.parse(self.test_hhr_path)
        self.assertIsNotNone(result, "Parser should not return None for valid HHR file")
        
        # Check that alignments are parsed
        hit1 = result['hits'][0]
        self.assertIn('alignment_blocks', hit1, "Alignments should be parsed")
        self.assertGreater(len(hit1['alignment_blocks']), 0, "Should have at least one alignment block")
        
        # Check alignment details for first hit
        block = hit1['alignment_blocks'][0]
        self.assertEqual(block['query_id'], 'test_protein', "Query ID not correctly parsed in alignment")
        self.assertEqual(block['query_start'], 1, "Query start position not correctly parsed")
        self.assertEqual(block['template_id'], 'e6gcsh1_consen', "Template ID not correctly parsed")
        self.assertEqual(block['template_start'], 8, "Template start position not correctly parsed")
        
        # Check sequence content
        self.assertTrue(block['query_seq'].startswith('MSWWSGVWRSVWSALSREV'), "Query sequence not correctly parsed")
        self.assertTrue(block['template_seq'].startswith('RNGLRGSWRQLNSLGDWRK'), "Template sequence not correctly parsed")
    
    def test_calculate_ranges(self):
        """Test range calculation from alignments"""
        # Parse the test HHR file
        result = self.parser.parse(self.test_hhr_path)
        self.assertIsNotNone(result, "Parser should not return None for valid HHR file")
        
        # Check calculated ranges
        hit1 = result['hits'][0]
        self.assertIn('query_range', hit1, "Query range should be calculated")
        self.assertIn('template_range', hit1, "Template range should be calculated")
        
        # Verify that ranges are properly formatted
        self.assertIsInstance(hit1['query_range'], str, "Query range should be a string")
        self.assertGreater(len(hit1['query_range']), 0, "Query range should not be empty")
        
        # Check that range contains start-end format
        self.assertRegex(hit1['query_range'], r'\d+-\d+(,\d+-\d+)*', "Query range should be in format start-end")
    
    def test_alignment_statistics(self):
        """Test alignment statistics calculation"""
        # Parse the test HHR file
        result = self.parser.parse(self.test_hhr_path)
        self.assertIsNotNone(result, "Parser should not return None for valid HHR file")
        
        # Check calculated statistics
        hit1 = result['hits'][0]
        self.assertIn('aligned_cols', hit1, "Aligned columns should be calculated")
        self.assertIn('identity_percentage', hit1, "Identity percentage should be calculated")
        self.assertIn('similarity_percentage', hit1, "Similarity percentage should be calculated")
        
        # Check statistics values
        self.assertIsInstance(hit1['identity_percentage'], float, "Identity percentage should be a float")
        self.assertGreaterEqual(hit1['identity_percentage'], 0.0, "Identity percentage should be >= 0")
        self.assertLessEqual(hit1['identity_percentage'], 100.0, "Identity percentage should be <= 100")
    
    def test_empty_file(self):
        """Test parsing of HHR file with no hits"""
        # Parse the empty HHR file
        result = self.parser.parse(self.empty_hhr_path)
        self.assertIsNotNone(result, "Parser should not return None for empty HHR file")
        
        # Check header and hits
        self.assertIn('header', result, "Header should be parsed")
        self.assertIn('hits', result, "Hits list should be present")
        self.assertEqual(len(result['hits']), 0, "Should find 0 hits in empty file")
    
    def test_malformed_file(self):
        """Test parsing of malformed HHR file"""
        # Parse the malformed HHR file
        result = self.parser.parse(self.malformed_hhr_path)
        self.assertIsNotNone(result, "Parser should not return None for malformed HHR file")
        
        # Check that parser handled malformed values gracefully
        header = result['header']
        self.assertEqual(header['query_id'], 'test_protein', "Should parse valid fields correctly")
        
        # Check how it handled non-numeric values
        self.assertIn('match_columns', header, "Should include field even if value is malformed")
        
        # Check that it attempted to parse hits
        self.assertIn('hits', result, "Hits list should be present")
    
    def test_robustness(self):
        """Test parser robustness with various edge cases"""
        # Test with non-existent file
        nonexistent_path = os.path.join(self.test_dir, "nonexistent.hhr")
        result = self.parser.parse(nonexistent_path)
        self.assertIsNone(result, "Parser should return None for non-existent file")
        
        # Test with empty file
        empty_path = os.path.join(self.test_dir, "empty.txt")
        with open(empty_path, 'w') as f:
            f.write("")
        result = self.parser.parse(empty_path)
        self.assertIsNotNone(result, "Parser should not crash on empty file")
        self.assertEqual(len(result['hits']), 0, "Should find 0 hits in completely empty file")


class TestHHRToXMLConverter(unittest.TestCase):
    """Tests for the HHR to XML converter"""
    
    def setUp(self):
        """Set up test data and converter instance"""
        self.test_dir = tempfile.mkdtemp()
        self.logger = logging.getLogger("test_logger")
        self.parser = HHRParser(self.logger)
        self.converter = HHRToXMLConverter(self.logger)
        
        # Create test HHR file
        self.test_hhr_path = os.path.join(self.test_dir, "test.hhr")
        with open(self.test_hhr_path, 'w') as f:
            f.write(TestHHRParser.get_test_hhr_content(self))
        
        # Parse HHR file for use in tests
        self.hhr_data = self.parser.parse(self.test_hhr_path)
    
    def tearDown(self):
        """Clean up test files"""
        shutil.rmtree(self.test_dir)
    
    def test_convert_to_xml(self):
        """Test conversion of parsed HHR data to XML"""
        # Convert parsed data to XML
        xml_string = self.converter.convert(
            self.hhr_data, 
            pdb_id="1abc", 
            chain_id="A", 
            ref_version="develop291"
        )
        
        self.assertIsNotNone(xml_string, "Converter should not return None")
        self.assertGreater(len(xml_string), 0, "XML output should not be empty")
        
        # Check XML structure
        self.assertTrue(xml_string.startswith('<?xml'), "Should start with XML declaration")
        self.assertIn("<hh_summ_doc>", xml_string, "Should have root element")
        self.assertIn("<metadata>", xml_string, "Should have metadata section")
        self.assertIn("<hh_hit_list>", xml_string, "Should have hit list section")
        
        # Check metadata content
        self.assertIn("<pdb_id>1abc</pdb_id>", xml_string, "Should include PDB ID")
        self.assertIn("<chain_id>A</chain_id>", xml_string, "Should include chain ID")
        self.assertIn("<reference>develop291</reference>", xml_string, "Should include reference version")
    
    def test_save_xml(self):
        """Test saving XML to file"""
        # Convert parsed data to XML
        xml_string = self.converter.convert(
            self.hhr_data, 
            pdb_id="1abc", 
            chain_id="A", 
            ref_version="develop291"
        )
        
        # Save XML to file
        output_path = os.path.join(self.test_dir, "output.xml")
        result = self.converter.save(xml_string, output_path)
        
        self.assertTrue(result, "Save method should return True on success")
        self.assertTrue(os.path.exists(output_path), "Output file should exist")
        
        # Check file content
        with open(output_path, 'r') as f:
            content = f.read()
            self.assertEqual(content, xml_string, "Saved content should match XML string")
    
    def test_empty_data(self):
        """Test converter with empty data"""
        # Create empty data structure
        empty_data = {'header': {}, 'hits': []}
        
        # Convert empty data to XML
        xml_string = self.converter.convert(
            empty_data, 
            pdb_id="1abc", 
            chain_id="A", 
            ref_version="develop291"
        )
        
        self.assertIsNotNone(xml_string, "Converter should not return None for empty data")
        self.assertGreater(len(xml_string), 0, "XML output should not be empty")
        self.assertIn("<hh_hit_list/>", xml_string, "Should have empty hit list")
    
    def test_missing_fields(self):
        """Test converter with data missing some fields"""
        # Modify data to remove some fields
        modified_data = self.hhr_data.copy()
        modified_data['hits'][0].pop('alignment_blocks', None)
        
        # Convert modified data to XML
        xml_string = self.converter.convert(
            modified_data, 
            pdb_id="1abc", 
            chain_id="A", 
            ref_version="develop291"
        )
        
        self.assertIsNotNone(xml_string, "Converter should not return None for data with missing fields")
        self.assertGreater(len(xml_string), 0, "XML output should not be empty")


if __name__ == '__main__':
    unittest.main()