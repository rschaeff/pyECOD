#!/usr/bin/env python3
"""
Domain Processing Debug Tests

Converted from detailed_domain_debug.py to proper pytest tests.
These tests help debug domain processing pipeline issues.
"""

import pytest
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path
from unittest.mock import Mock

# Import components to test
from ecod.core.context import ApplicationContext
from ecod.pipelines.domain_analysis.partition.service import DomainPartitionService


class TestDomainProcessingDebug:
    """Debug tests for domain processing pipeline"""
    
    @pytest.fixture
    def mock_context(self):
        """Create mock application context"""
        context = Mock(spec=ApplicationContext)
        context.config_manager = Mock()
        context.config_manager.config = {
            'database': {
                'host': 'localhost',
                'port': 5432,
                'database': 'ecod_test',
                'user': 'test_user',
                'password': 'test_pass'
            },
            'reference': {
                'current_version': 'develop291'
            },
            'partition': {
                'confidence_thresholds': {
                    'high': 0.9,
                    'medium': 0.7,
                    'low': 0.5
                },
                'evidence_weights': {
                    'domain_blast': 3.0,
                    'hhsearch': 2.5,
                    'chain_blast': 2.0,
                    'blast': 1.5
                },
                'overlap_tolerance': 0.15,
                'min_domain_size': 20,
                'peptide_threshold': 50
            }
        }
        context.config_manager.get_db_config.return_value = context.config_manager.config['database']
        return context
    
    @pytest.fixture  
    def test_xml_file(self, tmp_path):
        """Create test XML file with realistic data"""
        xml_content = """<?xml version="1.0" encoding="utf-8"?>
<blast_summ_doc>
    <blast_summ pdb="3hhp" chain="A"/>
    <blast_run program="blastp">
        <hits>
            <hit domain_id="e3hhpA1" evalues="1e-50" hsp_count="1">
                <query_reg>5-180</query_reg>
                <hit_reg>1-175</hit_reg>
            </hit>
            <hit domain_id="e3hhpA2" evalues="1e-45" hsp_count="1">
                <query_reg>175-350</query_reg>
                <hit_reg>1-175</hit_reg>
            </hit>
            <hit domain_id="e3hhpA3" evalues="1e-40" hsp_count="1">
                <query_reg>345-506</query_reg>
                <hit_reg>1-160</hit_reg>
            </hit>
        </hits>
    </blast_run>
    <hh_run program="hhsearch">
        <hits>
            <hit domain_id="e3hhpA1" probability="95.0" evalue="1e-25" score="80.0">
                <query_reg>10-185</query_reg>
                <hit_reg>5-180</hit_reg>
            </hit>
            <hit domain_id="e3hhpA2" probability="92.0" evalue="1e-22" score="75.0">
                <query_reg>180-355</query_reg>
                <hit_reg>3-178</hit_reg>
            </hit>
            <hit domain_id="e3hhpA3" probability="88.0" evalue="1e-20" score="70.0">
                <query_reg>350-500</query_reg>
                <hit_reg>1-150</hit_reg>
            </hit>
        </hits>
    </hh_run>
</blast_summ_doc>"""
        
        xml_file = tmp_path / "3hhp_A.summary.xml"
        xml_file.write_text(xml_content)
        return xml_file
    
    def test_xml_parsing_debug(self, test_xml_file):
        """Test XML parsing stage of pipeline"""
        from ecod.pipelines.domain_analysis.partition.analyzer import EvidenceAnalyzer
        from ecod.pipelines.domain_analysis.partition.models import PartitionOptions
        
        options = PartitionOptions()
        analyzer = EvidenceAnalyzer(options)
        
        # Parse XML
        result = analyzer.parse_domain_summary(str(test_xml_file))
        
        # Verify parsing worked
        assert 'error' not in result, f"XML parsing failed: {result.get('error')}"
        assert 'blast_hits' in result
        assert 'hhsearch_hits' in result
        
        # Check hit counts
        assert len(result['blast_hits']) == 3
        assert len(result['hhsearch_hits']) == 3
    
    def test_evidence_extraction_debug(self, test_xml_file):
        """Test evidence extraction from parsed XML"""
        from ecod.pipelines.domain_analysis.partition.analyzer import EvidenceAnalyzer
        from ecod.pipelines.domain_analysis.partition.models import PartitionOptions
        
        options = PartitionOptions()
        analyzer = EvidenceAnalyzer(options)
        
        # Parse and extract evidence
        summary_data = analyzer.parse_domain_summary(str(test_xml_file))
        evidence_list = analyzer.extract_evidence_with_classification(summary_data)
        
        # Verify evidence extraction
        assert len(evidence_list) > 0, "Should extract evidence items"
        
        # Check evidence properties
        for evidence in evidence_list:
            assert evidence.type in ['domain_blast', 'hhsearch']
            assert evidence.query_range is not None
            assert evidence.confidence > 0
    
    def test_domain_processing_debug(self, test_xml_file, mock_context):
        """Test full domain processing pipeline"""
        # This test helps debug the full pipeline
        service = DomainPartitionService(mock_context)
        
        # Process the test case
        result = service.partition_protein(
            pdb_id="3hhp",
            chain_id="A", 
            summary_path=str(test_xml_file),
            output_dir=str(test_xml_file.parent)
        )
        
        # Basic checks
        assert result is not None
        assert result.pdb_id == "3hhp"
        assert result.chain_id == "A"
        
        # If processing succeeded, check results
        if result.success:
            assert result.sequence_length > 0
            # Don't assert specific domain counts since this is debug
            print(f"Debug: Found {len(result.domains)} domains")
            print(f"Debug: Coverage = {result.coverage:.3f}")
    
    @pytest.mark.slow
    def test_coverage_calculation_debug(self, test_xml_file, mock_context):
        """Debug coverage calculation specifically"""
        service = DomainPartitionService(mock_context)
        result = service.partition_protein(
            pdb_id="3hhp", chain_id="A",
            summary_path=str(test_xml_file),
            output_dir=str(test_xml_file.parent)
        )
        
        if result.success and result.domains:
            # Calculate expected coverage manually
            sequence_length = 506  # From test XML
            covered_positions = set()
            
            for domain in result.domains:
                if hasattr(domain, 'start') and hasattr(domain, 'end'):
                    covered_positions.update(range(domain.start, domain.end + 1))
            
            manual_coverage = len(covered_positions) / sequence_length
            
            print(f"Debug: Pipeline coverage = {result.coverage:.3f}")
            print(f"Debug: Manual coverage = {manual_coverage:.3f}")
            print(f"Debug: Difference = {abs(result.coverage - manual_coverage):.3f}")
