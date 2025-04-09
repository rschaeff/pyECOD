#!/usr/bin/env python3
"""
Unit tests for the Domain Summary module
"""
import os
import unittest
from unittest.mock import patch, MagicMock, mock_open
import xml.etree.ElementTree as ET
import tempfile
import shutil

from ecod.pipelines.domain_analysis.summary import DomainSummary
from ecod.config import ConfigManager
from ecod.db.manager import DBManager
from ecod.exceptions import PipelineError, FileOperationError


class TestDomainSummary(unittest.TestCase):
    """Test cases for DomainSummary class"""
    
    def setUp(self):
        """Set up test fixtures"""
        # Create a mock config
        self.mock_config = {
            'database': {
                'host': 'localhost',
                'port': 5432,
                'user': 'test_user',
                'password': 'test_password',
                'database': 'test_db'
            },
            'paths': {
                'output_dir': '/tmp/ecod_test'
            },
            'force_overwrite': True
        }
        
        # Create a temporary directory for test files
        self.temp_dir = tempfile.mkdtemp()
        
        # Patch the ConfigManager
        self.config_patcher = patch('ecod.pipelines.domain_analysis.summary.ConfigManager')
        self.mock_config_manager = self.config_patcher.start()
        self.mock_config_manager.return_value.config = self.mock_config
        self.mock_config_manager.return_value.get_db_config.return_value = self.mock_config['database']
        
        # Patch the DBManager
        self.db_patcher = patch('ecod.pipelines.domain_analysis.summary.DBManager')
        self.mock_db = self.db_patcher.start()
        
        # Initialize the DomainSummary instance
        self.domain_summary = DomainSummary()
    
    def tearDown(self):
        """Tear down test fixtures"""
        # Stop patches
        self.config_patcher.stop()
        self.db_patcher.stop()
        
        # Remove temporary directory
        shutil.rmtree(self.temp_dir)
    
    def test_init(self):
        """Test initialization of DomainSummary"""
        self.assertIsNotNone(self.domain_summary)
        self.assertEqual(self.domain_summary.hsp_evalue_threshold, 0.005)
        self.assertEqual(self.domain_summary.hit_coverage_threshold, 0.7)
    
    @patch('os.path.exists')
    @patch('os.makedirs')
    def test_create_summary_existing_file(self, mock_makedirs, mock_exists):
        """Test create_summary when output file already exists"""
        # Setup
        mock_exists.return_value = True
        self.domain_summary.config['force_overwrite'] = False
        
        # Call the method
        result = self.domain_summary.create_summary('1abc', 'A', 'develop291', '/tmp/ecod_test', False)
        
        # Assertions
        self.assertIsNotNone(result)
        mock_exists.assert_called_once()
        mock_makedirs.assert_called_once()
    
    @patch('os.path.exists')
    @patch('os.makedirs')
    @patch('ecod.pipelines.domain_analysis.summary.ET')
    def test_create_summary_no_self_comparison(self, mock_et, mock_makedirs, mock_exists):
        """Test create_summary without self-comparison results"""
        # Setup
        mock_exists.side_effect = lambda path: False if 'self_comp.xml' in path else True
        mock_db = MagicMock()
        self.domain_summary.db = mock_db
        
        # Mock the ElementTree functionality
        mock_root = MagicMock()
        mock_et.Element.return_value = mock_root
        mock_blast_summ = MagicMock()
        mock_et.SubElement.return_value = mock_blast_summ
        mock_tree = MagicMock()
        mock_et.ElementTree.return_value = mock_tree
        
        # Mock database query for blast files
        mock_db.execute_query.return_value = [['some/path/file.xml']]
        
        # Mock file processing methods
        self.domain_summary._process_chain_blast = MagicMock()
        self.domain_summary._process_blast = MagicMock()
        
        # Call the method with patched file path checks
        with patch('builtins.open', mock_open()):
            result = self.domain_summary.create_summary('1abc', 'A', 'develop291', '/tmp/ecod_test', True)
        
        # Assertions
        self.assertIsNotNone(result)
        mock_blast_summ.set.assert_any_call("no_selfcomp", "true")
        mock_et.ElementTree.assert_called_once_with(mock_root)
        mock_tree.write.assert_called_once()
    
    def test_process_hit_hsps(self):
        """Test processing of HSPs for a hit"""
        # Create a mock hit node with HSPs
        hit_xml = """
        <Hit>
            <Hit_num>1</Hit_num>
            <Hit_len>100</Hit_len>
            <Hit_hsps>
                <Hsp>
                    <Hsp_evalue>0.001</Hsp_evalue>
                    <Hsp_align-len>80</Hsp_align-len>
                    <Hsp_query-from>10</Hsp_query-from>
                    <Hsp_query-to>90</Hsp_query-to>
                    <Hsp_hit-from>5</Hsp_hit-from>
                    <Hsp_hit-to>85</Hsp_hit-to>
                    <Hsp_qseq>SEQUENCE1</Hsp_qseq>
                    <Hsp_hseq>SEQUENCE2</Hsp_hseq>
                </Hsp>
            </Hit_hsps>
        </Hit>
        """
        hit_node = ET.fromstring(hit_xml)
        
        # Process the HSPs
        result = self.domain_summary._process_hit_hsps(hit_node, 100)
        
        # Assertions
        self.assertIsNotNone(result)
        self.assertEqual(len(result['query_regions']), 1)
        self.assertEqual(result['query_regions'][0], '10-90')
        self.assertEqual(result['hit_regions'][0], '5-85')
        self.assertEqual(result['query_seqs'][0], 'SEQUENCE1')
        self.assertEqual(result['hit_seqs'][0], 'SEQUENCE2')
        self.assertEqual(result['evalues'][0], '0.001')
    
    def test_residue_coverage(self):
        """Test residue coverage calculation"""
        set1 = [1, 2, 3, 4, 5]
        set2 = [3, 4, 5, 6, 7]
        
        result = self.domain_summary._residue_coverage(set1, set2)
        
        self.assertEqual(result, 3)  # Overlap: 3, 4, 5
    
    def test_process_self_comparison(self):
        """Test processing of self-comparison results"""
        # Create a sample self-comparison XML
        xml_content = """
        <self_comparison>
            <structural_repeat aligner="dali" zscore="8.5">
                <ref_range>10-50</ref_range>
                <mob_range>60-100</mob_range>
            </structural_repeat>
            <sequence_repeat_set aligner="hhrepid" type="alpha">
                <sequence_repeat>
                    <range>10-30</range>
                </sequence_repeat>
                <sequence_repeat>
                    <range>40-60</range>
                </sequence_repeat>
            </sequence_repeat_set>
        </self_comparison>
        """
        
        # Save to a temporary file
        temp_file = os.path.join(self.temp_dir, "self_comp.xml")
        with open(temp_file, 'w') as f:
            f.write(xml_content)
        
        # Create a parent node to add results to
        parent_node = ET.Element("parent")
        
        # Process the file
        self.domain_summary._process_self_comparison(temp_file, parent_node)
        
        # Assertions
        self.assertEqual(len(parent_node.findall(".//self_comp_run")), 2)
        self.assertEqual(len(parent_node.findall(".//hit")), 1)
        
        hit = parent_node.find(".//hit")
        self.assertEqual(hit.get("aligner"), "dali")
        self.assertEqual(hit.get("z_score"), "8.5")
        self.assertEqual(hit.find("query_reg").text, "10-50")
        self.assertEqual(hit.find("hit_reg").text, "60-100")
        
        repeat_set = parent_node.find(".//repeat_set")
        self.assertEqual(repeat_set.get("aligner"), "hhrepid")
        self.assertEqual(repeat_set.get("type"), "alpha")
        self.assertEqual(repeat_set.find("seqid_range").text, "10-30,40-60")
    
    @patch('ecod.pipelines.domain_analysis.summary.ET')
    def test_parse_hhr_file_error_handling(self, mock_et):
        """Test error handling during HHR file parsing"""
        # Setup mock to raise exception
        mock_et.parse.side_effect = Exception("Test exception")
        
        # Create a mock parent node
        parent_node = ET.Element("parent")
        
        # Call the method with a non-existent file
        self.domain_summary._process_hhsearch("non_existent_file.hhr", parent_node)
        
        # Verify error was logged but no exception was raised
        # Ideally we'd check the log message, but that's more complex
        # Just making sure the method doesn't crash is a start
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()