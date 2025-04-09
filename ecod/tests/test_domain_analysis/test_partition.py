#!/usr/bin/env python3
"""
Unit tests for the Domain Partition module
"""
import os
import unittest
from unittest.mock import patch, MagicMock, mock_open
import xml.etree.ElementTree as ET
import tempfile
import shutil

from ecod.pipelines.domain_analysis.partition import DomainPartition
from ecod.config import ConfigManager
from ecod.db.manager import DBManager
from ecod.exceptions import PipelineError, FileOperationError


class TestDomainPartition(unittest.TestCase):
    """Test cases for DomainPartition class"""
    
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
        self.config_patcher = patch('ecod.pipelines.domain_analysis.partition.ConfigManager')
        self.mock_config_manager = self.config_patcher.start()
        self.mock_config_manager.return_value.config = self.mock_config
        self.mock_config_manager.return_value.get_db_config.return_value = self.mock_config['database']
        
        # Patch the DBManager
        self.db_patcher = patch('ecod.pipelines.domain_analysis.partition.DBManager')
        self.mock_db = self.db_patcher.start()
        
        # Initialize the DomainPartition instance
        self.domain_partition = DomainPartition()
    
    def tearDown(self):
        """Tear down test fixtures"""
        # Stop patches
        self.config_patcher.stop()
        self.db_patcher.stop()
        
        # Remove temporary directory
        shutil.rmtree(self.temp_dir)
    
    def test_init(self):
        """Test initialization of DomainPartition"""
        self.assertIsNotNone(self.domain_partition)
        self.assertEqual(self.domain_partition.new_coverage_threshold, 0.9)
        self.assertEqual(self.domain_partition.old_coverage_threshold, 0.05)
        self.assertEqual(self.domain_partition.dali_significance_threshold, 5.0)
        self.assertEqual(self.domain_partition.hit_coverage_threshold, 0.7)
        self.assertEqual(self.domain_partition.gap_tol, 20)
    
    def test_get_start_position(self):
        """Test getting start position from range string"""
        # Test simple range
        result = self.domain_partition._get_start_position("10-50")
        self.assertEqual(result, 10)
        
        # Test multi-segment range
        result = self.domain_partition._get_start_position("10-20,30-40")
        self.assertEqual(result, 10)
        
        # Test single position
        result = self.domain_partition._get_start_position("15")
        self.assertEqual(result, 15)
        
        # Test invalid range
        result = self.domain_partition._get_start_position("invalid")
        self.assertEqual(result, 0)
    
    def test_get_end_position(self):
        """Test getting end position from range string"""
        # Test simple range
        result = self.domain_partition._get_end_position("10-50")
        self.assertEqual(result, 50)
        
        # Test multi-segment range
        result = self.domain_partition._get_end_position("10-20,30-40")
        self.assertEqual(result, 40)
        
        # Test single position
        result = self.domain_partition._get_end_position("15")
        self.assertEqual(result, 15)
        
        # Test invalid range
        result = self.domain_partition._get_end_position("invalid")
        self.assertEqual(result, 0)
    
    def test_range_to_positions(self):
        """Test converting range string to set of positions"""
        # Test simple range
        result = self.domain_partition._range_to_positions("10-15")
        self.assertEqual(result, {10, 11, 12, 13, 14, 15})
        
        # Test multi-segment range
        result = self.domain_partition._range_to_positions("10-12,14-16")
        self.assertEqual(result, {10, 11, 12, 14, 15, 16})
        
        # Test single position
        result = self.domain_partition._range_to_positions("10")
        self.assertEqual(result, {10})
        
        # Test empty string
        result = self.domain_partition._range_to_positions("")
        self.assertEqual(result, set())
        
        # Test invalid range
        result = self.domain_partition._range_to_positions("invalid")
        self.assertEqual(result, set())
    
    def test_calculate_overlap_percentage(self):
        """Test calculating overlap percentage between ranges"""
        # Test complete overlap
        result = self.domain_partition._calculate_overlap_percentage("10-20", "10-20", 100)
        self.assertEqual(result, 1.0)
        
        # Test partial overlap
        result = self.domain_partition._calculate_overlap_percentage("10-20", "15-25", 100)
        self.assertEqual(result, 0.5)  # 5/10 overlap
        
        # Test no overlap
        result = self.domain_partition._calculate_overlap_percentage("10-20", "30-40", 100)
        self.assertEqual(result, 0.0)
    
    def test_combine_ranges(self):
        """Test combining multiple range strings"""
        # Test simple ranges
        result = self.domain_partition._combine_ranges(["10-20", "30-40"])
        self.assertEqual(result, "10-20,30-40")
        
        # Test overlapping ranges
        result = self.domain_partition._combine_ranges(["10-25", "20-30"])
        self.assertEqual(result, "10-30")
        
        # Test adjacent ranges
        result = self.domain_partition._combine_ranges(["10-20", "21-30"])
        self.assertEqual(result, "10-30")
        
        # Test single positions
        result = self.domain_partition._combine_ranges(["10", "12", "14"])
        self.assertEqual(result, "10-10,12-12,14-14")
        
        # Test empty list
        result = self.domain_partition._combine_ranges([])
        self.assertEqual(result, "")
        
        # Test list with empty string
        result = self.domain_partition._combine_ranges([""])
        self.assertEqual(result, "")
    
    @patch('os.path.exists')
    @patch('os.makedirs')
    @patch('xml.etree.ElementTree.ElementTree.write')
    def test_partition_domains_existing_file(self, mock_write, mock_makedirs, mock_exists):
        """Test partition_domains when output file already exists"""
        # Setup
        mock_exists.return_value = True
        self.domain_partition.config['force_overwrite'] = False
        
        # Call the method
        result = self.domain_partition.partition_domains('1abc', 'A', '/tmp/ecod_test', 'struct_seqid', 'develop291')
        
        # Assertions
        self.assertEqual(result, os.path.join('/tmp/ecod_test', '1abc_A', 'domains_v12.1abc_A.develop291.xml'))
        mock_exists.assert_called_once()
        mock_makedirs.assert_not_called()
        mock_write.assert_not_called()
    
    @patch('builtins.open', new_callable=mock_open, read_data=">1abc_A\nSEQUENCEABC")
    @patch('os.path.exists')
    @patch('os.makedirs')
    def test_read_fasta_sequence(self, mock_makedirs, mock_exists, mock_file):
        """Test reading sequence from FASTA file"""
        # Setup
        mock_exists.return_value = True
        
        # Call the method
        result = self.domain_partition._read_fasta_sequence('/path/to/file.fa')
        
        # Assertions
        self.assertEqual(result, "SEQUENCEABC")
        mock_file.assert_called_once_with('/path/to/file.fa', 'r')
    
    def test_process_blast_summary(self):
        """Test processing BLAST summary file"""
        # Create a sample blast summary XML
        xml_content = """
        <blast_summ_doc>
            <blast_summ>
                <chain_blast_run>
                    <hits>
                        <hit num="1" pdb_id="2xyz" chain_id="B" evalues="0.001,0.002">
                            <query_reg>10-50,60-100</query_reg>
                            <hit_reg>5-45,55-95</hit_reg>
                        </hit>
                    </hits>
                </chain_blast_run>
                <blast_run>
                    <hits>
                        <hit num="1" domain_id="d2xyzB1" pdb_id="2xyz" chain_id="B" evalues="0.001">
                            <query_reg>10-50</query_reg>
                            <hit_reg>5-45</hit_reg>
                        </hit>
                    </hits>
                </blast_run>
                <hh_run>
                    <hits>
                        <hit num="1" domain_id="d2xyzB1" hh_prob="95.5" hh_score="150.2">
                            <query_reg>10-50</query_reg>
                            <hit_reg>5-45</hit_reg>
                        </hit>
                    </hits>
                </hh_run>
                <self_comp_run programs="dali">
                    <hits>
                        <hit aligner="dali" z_score="8.5">
                            <query_reg>10-50</query_reg>
                            <hit_reg>60-100</hit_reg>
                        </hit>
                    </hits>
                </self_comp_run>
            </blast_summ>
        </blast_summ_doc>
        """
        
        # Save to a temporary file
        temp_file = os.path.join(self.temp_dir, "blast_summ.xml")
        with open(temp_file, 'w') as f:
            f.write(xml_content)
        
        # Process the file
        result = self.domain_partition._process_blast_summary(temp_file)
        
        # Assertions
        self.assertIn('chain_blast', result)
        self.assertIn('domain_blast', result)
        self.assertIn('hhsearch', result)
        self.assertIn('self_comparison', result)
        
        self.assertEqual(len(result['chain_blast']), 1)
        self.assertEqual(len(result['domain_blast']), 1)
        self.assertEqual(len(result['hhsearch']), 1)
        self.assertEqual(len(result['self_comparison']), 1)
        
        # Check chain blast data
        chain_hit = result['chain_blast'][0]
        self.assertEqual(chain_hit['pdb_id'], '2xyz')
        self.assertEqual(chain_hit['chain_id'], 'B')
        self.assertEqual(chain_hit['evalues'], ['0.001', '0.002'])
        self.assertEqual(chain_hit['query_regions'], ['10-50', '60-100'])
        
        # Check HHsearch data
        hh_hit = result['hhsearch'][0]
        self.assertEqual(hh_hit['domain_id'], 'd2xyzB1')
        self.assertEqual(float(hh_hit['probability']), 95.5)
        self.assertEqual(float(hh_hit['score']), 150.2)
    
    def test_determine_domain_boundaries(self):
        """Test determining domain boundaries from blast data"""
        # Mock blast data
        blast_data = {
            'chain_blast': [],
            'domain_blast': [
                {
                    'domain_id': 'd1xyzA1',
                    'evalues': ['0.001'],
                    'query_regions': ['10-50'],
                    'hit_regions': ['5-45']
                }
            ],
            'hhsearch': [
                {
                    'domain_id': 'd2abcB1',
                    'probability': 98.5,
                    'query_region': '60-100',
                    'hit_region': '55-95'
                }
            ],
            'self_comparison': [
                {
                    'aligner': 'dali',
                    'z_score': 8.5,
                    'query_region': '110-150',
                    'hit_region': '160-200'
                }
            ]
        }
        
        # Patch methods used by determine_domain_boundaries
        with patch.multiple(
            self.domain_partition,
            _identify_repeats=MagicMock(return_value=[
                {'range': '110-150', 'quality': 8.5, 'type': 'repeat', 'evidence': []},
                {'range': '160-200', 'quality': 8.5, 'type': 'repeat', 'evidence': []}
            ]),
            _identify_domains_from_blast=MagicMock(return_value=[
                {'range': '10-50', 'quality': 50.0, 'type': 'blast', 'evidence': []}
            ]),
            _identify_domains_from_hhsearch=MagicMock(return_value=[
                {'range': '60-100', 'quality': 98.5, 'type': 'hhsearch', 'evidence': []}
            ]),
            _resolve_domain_boundaries=MagicMock(return_value=[
                {'range': '10-50', 'quality': 50.0, 'type': 'blast', 'evidence': []},
                {'range': '60-100', 'quality': 98.5, 'type': 'hhsearch', 'evidence': []},
                {'range': '110-150', 'quality': 8.5, 'type': 'repeat', 'evidence': []},
                {'range': '160-200', 'quality': 8.5, 'type': 'repeat', 'evidence': []}
            ])
        ) as mocks:
            # Run the method
            result = self.domain_partition._determine_domain_boundaries(blast_data, 200, '1abc_A')
            
            # Assertions
            self.assertEqual(len(result), 4)
            mocks['_identify_repeats'].assert_called_once()
            mocks['_identify_domains_from_blast'].assert_called_once()
            mocks['_identify_domains_from_hhsearch'].assert_called_once()
            mocks['_resolve_domain_boundaries'].assert_called_once()
    
    @patch('ecod.pipelines.domain_analysis.partition.DBManager')
    def test_load_reference_data(self, mock_db_manager):
        """Test loading reference domain classifications"""
        # Setup mock database results
        mock_db = MagicMock()
        mock_db_manager.return_value = mock_db
        mock_db.execute_dict_query.return_value = [
            {'ecod_uid': 1, 'domain_id': 'd1abcA1', 'range': '10-50', 'source_id': '1abc_A'},
            {'ecod_uid': 2, 'domain_id': 'd1abcA2', 'range': '60-100', 'source_id': '1abc_A'},
            {'ecod_uid': 3, 'domain_id': 'd2xyzB1', 'range': '5-45', 'source_id': '2xyz_B'}
        ]
        
        # Call the method
        self.domain_partition.load_reference_data('develop291')
        
        # Assertions
        mock_db.execute_dict_query.assert_called_once()
        
        # Check ref_range_cache
        self.assertIn('1abc_A', self.domain_partition.ref_range_cache)
        self.assertIn('2xyz_B', self.domain_partition.ref_range_cache)
        self.assertEqual(len(self.domain_partition.ref_range_cache['1abc_A']), 2)
        self.assertEqual(len(self.domain_partition.ref_range_cache['2xyz_B']), 1)
        
        # Check ref_domain_uid_lookup
        self.assertEqual(self.domain_partition.ref_domain_uid_lookup['d1abcA1'], 1)
        self.assertEqual(self.domain_partition.ref_domain_uid_lookup['d1abcA2'], 2)
        self.assertEqual(self.domain_partition.ref_domain_uid_lookup['d2xyzB1'], 3)
        
        # Check ref_chain_domains
        self.assertIn('1abc_A', self.domain_partition.ref_chain_domains)
        self.assertIn('2xyz_B', self.domain_partition.ref_chain_domains)
        self.assertEqual(len(self.domain_partition.ref_chain_domains['1abc_A']), 2)
        self.assertEqual(len(self.domain_partition.ref_chain_domains['2xyz_B']), 1)

    def test_resolve_domain_boundaries(self):
        """Test resolving domain boundaries"""
        # Test with empty candidates - should return full chain
        result = self.domain_partition._resolve_domain_boundaries([], 100)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]['range'], '1-100')
        self.assertEqual(result[0]['type'], 'full_chain')
        
        # Test with single candidate
        candidates = [{'range': '10-50', 'quality': 50.0, 'type': 'blast'}]
        result = self.domain_partition._resolve_domain_boundaries(candidates, 100)
        self.assertEqual(len(result), 3)  # Original + gap at start + gap at end
        self.assertEqual(result[0]['range'], '10-50')
        
        # Test with multiple overlapping candidates
        candidates = [
            {'range': '10-60', 'quality': 50.0, 'type': 'blast'},
            {'range': '40-90', 'quality': 98.5, 'type': 'hhsearch'},
        ]
        result = self.domain_partition._resolve_domain_boundaries(candidates, 100)
        # Should keep the higher quality one for the overlap and add gaps
        domains = [d for d in result if d['type