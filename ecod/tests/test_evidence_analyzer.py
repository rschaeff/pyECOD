#!/usr/bin/env python3
"""
Comprehensive tests for EvidenceAnalyzer class.

Tests core evidence analysis functionality including:
- Domain summary parsing
- Evidence extraction and validation
- Evidence grouping and consensus building
- Domain validation
- XML parsing error handling
"""

import pytest
import tempfile
import os
import xml.etree.ElementTree as ET
from unittest.mock import Mock, MagicMock, patch, mock_open
from datetime import datetime
from pathlib import Path

from ecod.models.pipeline.evidence import Evidence
from ecod.models.pipeline.domain import DomainModel
from ecod.pipelines.domain_analysis.partition.models import (
    PartitionOptions, ValidationLevel, ValidationResult, EvidenceGroup
)
from ecod.pipelines.domain_analysis.partition.analyzer import EvidenceAnalyzer


class TestEvidenceAnalyzerInitialization:
    """Test EvidenceAnalyzer initialization and configuration"""
    
    def test_analyzer_initialization(self):
        """Test basic analyzer initialization"""
        options = PartitionOptions()
        analyzer = EvidenceAnalyzer(options)
        
        assert analyzer.options == options
        assert analyzer.classification_cache is not None
        assert analyzer.stats is not None
        assert 'summaries_parsed' in analyzer.stats
        assert 'evidence_extracted' in analyzer.stats
        assert 'validation_errors' in analyzer.stats
    
    def test_analyzer_with_custom_options(self):
        """Test analyzer with custom validation options"""
        options = PartitionOptions(
            validation_level=ValidationLevel.STRICT,
            require_evidence=True,
            validate_coordinates=True,
            min_evidence_confidence=0.5
        )
        
        analyzer = EvidenceAnalyzer(options)
        
        assert analyzer.options.validation_level == ValidationLevel.STRICT
        assert analyzer.options.require_evidence == True
        assert analyzer.options.min_evidence_confidence == 0.5


class TestDomainSummaryParsing:
    """Test domain summary XML parsing"""
    
    @pytest.fixture
    def analyzer(self):
        """Create analyzer for testing"""
        options = PartitionOptions()
        return EvidenceAnalyzer(options)
    
    def test_parse_valid_summary(self, analyzer):
        """Test parsing a valid domain summary"""
        summary_xml = '''<?xml version="1.0"?>
        <domain_summary pdb_id="1abc" chain_id="A" sequence_length="250">
            <metadata>
                <is_peptide>false</is_peptide>
                <is_discontinuous>false</is_discontinuous>
            </metadata>
            <evidence_list>
                <blast_hits>
                    <hit num="1" domain_id="d1abcA1" evalue="1e-50" evalues="1e-50">
                        <query_range>10-100</query_range>
                        <hit_range>5-95</hit_range>
                    </hit>
                </blast_hits>
                <hhsearch_hits>
                    <hit hit_id="h1" domain_id="d2xyzB1" probability="95.5" evalue="1e-30">
                        <query_range>120-200</query_range>
                        <hit_range>10-90</hit_range>
                    </hit>
                </hhsearch_hits>
            </evidence_list>
        </domain_summary>'''
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xml', delete=False) as f:
            f.write(summary_xml)
            f.flush()
            
            try:
                result = analyzer.parse_domain_summary(f.name)
                
                assert result['pdb_id'] == '1abc'
                assert result['chain_id'] == 'A'
                assert result['sequence_length'] == 250
                assert result['is_peptide'] == False
                assert 'blast_hits' in result
                assert 'hhsearch_hits' in result
                assert len(result['blast_hits']) == 1
                assert len(result['hhsearch_hits']) == 1
                
            finally:
                os.unlink(f.name)
    
    def test_parse_peptide_summary(self, analyzer):
        """Test parsing summary for peptide"""
        summary_xml = '''<?xml version="1.0"?>
        <domain_summary pdb_id="1abc" chain_id="A" sequence_length="15">
            <metadata>
                <is_peptide>true</is_peptide>
            </metadata>
        </domain_summary>'''
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xml', delete=False) as f:
            f.write(summary_xml)
            f.flush()
            
            try:
                result = analyzer.parse_domain_summary(f.name)
                
                assert result['is_peptide'] == True
                assert result['sequence_length'] == 15
                
            finally:
                os.unlink(f.name)
    
    def test_parse_missing_file(self, analyzer):
        """Test parsing non-existent file"""
        result = analyzer.parse_domain_summary('/nonexistent/file.xml')
        
        assert 'error' in result
        assert 'not found' in result['error'].lower()
    
    def test_parse_malformed_xml(self, analyzer):
        """Test parsing malformed XML"""
        malformed_xml = '''<?xml version="1.0"?>
        <domain_summary pdb_id="1abc" chain_id="A">
            <unclosed_tag>
        </domain_summary>'''
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xml', delete=False) as f:
            f.write(malformed_xml)
            f.flush()
            
            try:
                result = analyzer.parse_domain_summary(f.name)
                
                assert 'error' in result
                assert 'xml' in result['error'].lower()
                
            finally:
                os.unlink(f.name)
    
    def test_parse_empty_summary(self, analyzer):
        """Test parsing empty summary"""
        empty_xml = '''<?xml version="1.0"?>
        <domain_summary pdb_id="1abc" chain_id="A" sequence_length="100">
        </domain_summary>'''
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xml', delete=False) as f:
            f.write(empty_xml)
            f.flush()
            
            try:
                result = analyzer.parse_domain_summary(f.name)
                
                assert result['pdb_id'] == '1abc'
                assert result['chain_id'] == 'A'
                assert result['sequence_length'] == 100
                assert result.get('blast_hits', []) == []
                assert result.get('hhsearch_hits', []) == []
                
            finally:
                os.unlink(f.name)
    
    def test_parse_discontinuous_summary(self, analyzer):
        """Test parsing summary with discontinuous domains"""
        discontinuous_xml = '''<?xml version="1.0"?>
        <domain_summary pdb_id="1abc" chain_id="A" sequence_length="300">
            <metadata>
                <is_discontinuous>true</is_discontinuous>
            </metadata>
            <evidence_list>
                <blast_hits>
                    <hit num="1" domain_id="d1abcA1" discontinuous="true">
                        <query_range>10-50,100-150,200-250</query_range>
                        <hit_range>5-45,95-145,195-245</hit_range>
                    </hit>
                </blast_hits>
            </evidence_list>
        </domain_summary>'''
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xml', delete=False) as f:
            f.write(discontinuous_xml)
            f.flush()
            
            try:
                result = analyzer.parse_domain_summary(f.name)
                
                assert result['is_discontinuous'] == True
                assert len(result['blast_hits']) == 1
                hit = result['blast_hits'][0]
                assert hit.get('discontinuous') == 'true'
                assert ',' in hit.get('query_range', '')
                
            finally:
                os.unlink(f.name)


class TestEvidenceExtraction:
    """Test evidence extraction and conversion"""
    
    @pytest.fixture
    def analyzer(self):
        """Create analyzer for testing"""
        options = PartitionOptions(validation_level=ValidationLevel.NORMAL)
        return EvidenceAnalyzer(options)
    
    @pytest.fixture
    def mock_db_lookup(self):
        """Mock database lookup function"""
        def lookup_func(domain_id):
            classifications = {
                'd1abcA1': {
                    't_group': '1.1.1.1',
                    'h_group': '1.1.1',
                    'x_group': '1.1',
                    'a_group': 'a.1'
                },
                'd2xyzB1': {
                    't_group': '2.2.2.2',
                    'h_group': '2.2.2',
                    'x_group': '2.2',
                    'a_group': 'a.2'
                }
            }
            return classifications.get(domain_id)
        
        return lookup_func
    
    def test_extract_evidence_blast_only(self, analyzer):
        """Test extracting only BLAST evidence"""
        summary_data = {
            'blast_hits': [
                {
                    'num': '1',
                    'domain_id': 'd1abcA1',
                    'evalues': '1e-50',
                    'query_range': '10-100',
                    'hit_range': '5-95'
                }
            ],
            'hhsearch_hits': [
                {
                    'hit_id': 'h1',
                    'domain_id': 'd2xyzB1',
                    'probability': '95.5',
                    'query_range': '120-200'
                }
            ]
        }
        
        # Set to blast only
        analyzer.options.use_hhsearch = False
        analyzer.options.use_domain_blast = True
        
        evidence_list = analyzer.extract_evidence_with_classification(summary_data)
        
        # Should only get BLAST evidence
        assert len(evidence_list) == 1
        assert evidence_list[0].type == 'domain_blast'
        assert evidence_list[0].domain_id == 'd1abcA1'
        assert evidence_list[0].evalue == 1e-50
    
    def test_extract_evidence_with_classification(self, analyzer, mock_db_lookup):
        """Test extracting evidence with classification lookup"""
        summary_data = {
            'blast_hits': [
                {
                    'num': '1',
                    'domain_id': 'd1abcA1',
                    'evalues': '1e-50',
                    'query_range': '10-100'
                }
            ],
            'hhsearch_hits': [
                {
                    'hit_id': 'h1',
                    'domain_id': 'd2xyzB1',
                    'probability': '95.5',
                    'query_range': '120-200'
                }
            ]
        }
        
        evidence_list = analyzer.extract_evidence_with_classification(
            summary_data, 
            db_lookup_func=mock_db_lookup
        )
        
        assert len(evidence_list) == 2
        
        # Check BLAST evidence
        blast_evidence = next(e for e in evidence_list if e.type == 'domain_blast')
        assert blast_evidence.domain_id == 'd1abcA1'
        assert blast_evidence.t_group == '1.1.1.1'
        assert blast_evidence.h_group == '1.1.1'
        assert blast_evidence.x_group == '1.1'
        assert blast_evidence.a_group == 'a.1'
        
        # Check HHSearch evidence
        hh_evidence = next(e for e in evidence_list if e.type == 'hhsearch')
        assert hh_evidence.domain_id == 'd2xyzB1'
        assert hh_evidence.t_group == '2.2.2.2'
        assert hh_evidence.probability == 95.5
    
    def test_extract_evidence_filtering_low_confidence(self, analyzer):
        """Test filtering out low confidence evidence"""
        analyzer.options.min_evidence_confidence = 0.8
        
        summary_data = {
            'blast_hits': [
                {
                    'num': '1',
                    'domain_id': 'd1abcA1',
                    'evalues': '1e-50',  # High confidence
                    'query_range': '10-100'
                },
                {
                    'num': '2', 
                    'domain_id': 'd2abcA2',
                    'evalues': '1e-2',   # Low confidence
                    'query_range': '150-200'
                }
            ]
        }
        
        evidence_list = analyzer.extract_evidence_with_classification(summary_data)
        
        # Should only get high confidence evidence
        assert len(evidence_list) == 1
        assert evidence_list[0].domain_id == 'd1abcA1'
        assert evidence_list[0].confidence >= 0.8
    
    def test_extract_evidence_empty_summary(self, analyzer):
        """Test extracting from empty summary"""
        summary_data = {}
        
        evidence_list = analyzer.extract_evidence_with_classification(summary_data)
        
        assert evidence_list == []
    
    def test_extract_evidence_malformed_hits(self, analyzer):
        """Test handling malformed hit data"""
        summary_data = {
            'blast_hits': [
                {
                    'num': '1',
                    'domain_id': 'd1abcA1',
                    'evalues': 'invalid_evalue',  # Malformed
                    'query_range': '10-100'
                },
                {
                    'num': '2',
                    'domain_id': 'd2abcA2',
                    'evalues': '1e-30',  # Valid
                    'query_range': '150-200'
                }
            ]
        }
        
        evidence_list = analyzer.extract_evidence_with_classification(summary_data)
        
        # Should get both (malformed handled gracefully)
        assert len(evidence_list) == 2
        
        # Find the evidence with malformed evalue
        malformed_evidence = next(e for e in evidence_list if e.domain_id == 'd1abcA1')
        assert malformed_evidence.evalue == 999.0  # Fallback value


class TestEvidenceValidation:
    """Test evidence validation logic"""
    
    @pytest.fixture
    def analyzer(self):
        """Create analyzer for testing"""
        options = PartitionOptions()
        return EvidenceAnalyzer(options)
    
    def test_validate_valid_evidence(self, analyzer):
        """Test validating good evidence"""
        evidence = Evidence(
            type='hhsearch',
            domain_id='d1abcA1',
            query_range='10-100',
            hit_range='5-95',
            probability=95.5,
            confidence=0.95
        )
        
        result = analyzer.validate_evidence(evidence, 'test_context')
        
        assert result.is_valid == True
        assert len(result.errors) == 0
    
    def test_validate_evidence_missing_domain_id(self, analyzer):
        """Test validating evidence without domain_id"""
        evidence = Evidence(
            type='hhsearch',
            domain_id='',  # Missing
            query_range='10-100',
            probability=95.5
        )
        
        result = analyzer.validate_evidence(evidence, 'test_context')
        
        assert result.is_valid == False
        assert any('domain_id' in error.lower() for error in result.errors)
    
    def test_validate_evidence_invalid_range(self, analyzer):
        """Test validating evidence with invalid range"""
        evidence = Evidence(
            type='hhsearch',
            domain_id='d1abcA1',
            query_range='invalid-range',  # Invalid format
            probability=95.5
        )
        
        result = analyzer.validate_evidence(evidence, 'test_context')
        
        # Should have warnings about range format
        assert len(result.warnings) > 0 or len(result.errors) > 0
    
    def test_validate_evidence_low_confidence(self, analyzer):
        """Test validating low confidence evidence"""
        analyzer.options.min_evidence_confidence = 0.8
        
        evidence = Evidence(
            type='blast',
            domain_id='d1abcA1',
            query_range='10-100',
            evalue=1e-2,  # Low confidence
            confidence=0.3
        )
        
        result = analyzer.validate_evidence(evidence, 'test_context')
        
        # With normal validation, should be valid but with warnings
        if analyzer.options.validation_level != ValidationLevel.STRICT:
            assert result.is_valid == True
            assert len(result.warnings) > 0
    
    def test_validate_evidence_strict_mode(self, analyzer):
        """Test validation in strict mode"""
        analyzer.options.validation_level = ValidationLevel.STRICT
        analyzer.options.min_evidence_confidence = 0.8
        
        evidence = Evidence(
            type='blast',
            domain_id='d1abcA1',
            query_range='10-100',
            evalue=1e-2,  # Low confidence
            confidence=0.3
        )
        
        result = analyzer.validate_evidence(evidence, 'test_context')
        
        # In strict mode, should be invalid
        assert result.is_valid == False
        assert any('confidence' in error.lower() for error in result.errors)
    
    def test_validate_evidence_missing_coordinates(self, analyzer):
        """Test validating evidence missing coordinates"""
        analyzer.options.validate_coordinates = True
        
        evidence = Evidence(
            type='hhsearch',
            domain_id='d1abcA1',
            query_range='',  # Missing
            probability=95.5
        )
        
        result = analyzer.validate_evidence(evidence, 'test_context')
        
        assert len(result.warnings) > 0 or len(result.errors) > 0


class TestEvidenceGrouping:
    """Test evidence grouping by position"""
    
    @pytest.fixture
    def analyzer(self):
        """Create analyzer for testing"""
        options = PartitionOptions()
        return EvidenceAnalyzer(options)
    
    def test_group_evidence_single_region(self, analyzer):
        """Test grouping evidence in same region"""
        evidence_list = [
            Evidence(type='blast', domain_id='d1', query_range='10-50', confidence=0.8),
            Evidence(type='hhsearch', domain_id='d2', query_range='15-55', confidence=0.9),
            Evidence(type='blast', domain_id='d3', query_range='20-60', confidence=0.7)
        ]
        
        groups = analyzer.group_evidence_by_position(evidence_list, window_size=50)
        
        # Should group into one region
        assert len(groups) == 1
        group = list(groups.values())[0]
        assert len(group.evidence_items) == 3
        assert group.consensus_confidence > 0
    
    def test_group_evidence_separate_regions(self, analyzer):
        """Test grouping evidence in separate regions"""
        evidence_list = [
            Evidence(type='blast', domain_id='d1', query_range='10-50', confidence=0.8),
            Evidence(type='hhsearch', domain_id='d2', query_range='100-150', confidence=0.9),
            Evidence(type='blast', domain_id='d3', query_range='200-250', confidence=0.7)
        ]
        
        groups = analyzer.group_evidence_by_position(evidence_list, window_size=30)
        
        # Should create separate groups
        assert len(groups) == 3
        for group in groups.values():
            assert len(group.evidence_items) == 1
    
    def test_group_evidence_overlapping_regions(self, analyzer):
        """Test grouping overlapping evidence"""
        evidence_list = [
            Evidence(type='blast', domain_id='d1', query_range='10-50', confidence=0.8),
            Evidence(type='hhsearch', domain_id='d2', query_range='40-80', confidence=0.9),
            Evidence(type='blast', domain_id='d3', query_range='70-110', confidence=0.7)
        ]
        
        groups = analyzer.group_evidence_by_position(evidence_list, window_size=40)
        
        # Should intelligently group overlapping evidence
        assert len(groups) >= 1
        
        # Check that consensus positions are reasonable
        for group in groups.values():
            assert group.consensus_start is not None
            assert group.consensus_end is not None
            assert group.consensus_start <= group.consensus_end
    
    def test_group_evidence_empty_list(self, analyzer):
        """Test grouping empty evidence list"""
        groups = analyzer.group_evidence_by_position([], window_size=50)
        
        assert len(groups) == 0
    
    def test_group_evidence_invalid_ranges(self, analyzer):
        """Test grouping evidence with invalid ranges"""
        evidence_list = [
            Evidence(type='blast', domain_id='d1', query_range='invalid', confidence=0.8),
            Evidence(type='hhsearch', domain_id='d2', query_range='10-50', confidence=0.9),
            Evidence(type='blast', domain_id='d3', query_range='', confidence=0.7)
        ]
        
        groups = analyzer.group_evidence_by_position(evidence_list, window_size=50)
        
        # Should handle invalid ranges gracefully
        # Only valid evidence should be grouped
        assert len(groups) >= 0
        
        # Find the group with valid evidence
        valid_groups = [g for g in groups.values() if len(g.evidence_items) > 0]
        if valid_groups:
            assert any(e.domain_id == 'd2' for g in valid_groups for e in g.evidence_items)


class TestDomainValidation:
    """Test domain model validation"""
    
    @pytest.fixture
    def analyzer(self):
        """Create analyzer for testing"""
        options = PartitionOptions()
        return EvidenceAnalyzer(options)
    
    def test_validate_valid_domain(self, analyzer):
        """Test validating a good domain"""
        domain = DomainModel(
            id='test_domain',
            start=10,
            end=100,
            range='10-100',
            confidence=0.9,
            source='hhsearch'
        )
        
        result = analyzer.validate_domain(domain, 'test_context', sequence_length=200)
        
        assert result.is_valid == True
        assert len(result.errors) == 0
    
    def test_validate_domain_too_small(self, analyzer):
        """Test validating domain below minimum size"""
        analyzer.options.min_domain_size = 30
        
        domain = DomainModel(
            id='small_domain',
            start=10,
            end=25,  # Only 16 residues
            range='10-25',
            confidence=0.9
        )
        
        result = analyzer.validate_domain(domain, 'test_context', sequence_length=200)
        
        assert result.is_valid == False
        assert any('size' in error.lower() for error in result.errors)
    
    def test_validate_domain_too_large(self, analyzer):
        """Test validating domain above maximum size"""
        analyzer.options.max_domain_size = 100
        
        domain = DomainModel(
            id='large_domain',
            start=10,
            end=150,  # 141 residues
            range='10-150',
            confidence=0.9
        )
        
        result = analyzer.validate_domain(domain, 'test_context', sequence_length=200)
        
        assert result.is_valid == False
        assert any('size' in error.lower() for error in result.errors)
    
    def test_validate_domain_outside_sequence(self, analyzer):
        """Test validating domain extending beyond sequence"""
        domain = DomainModel(
            id='boundary_domain',
            start=180,
            end=250,  # Beyond sequence length of 200
            range='180-250',
            confidence=0.9
        )
        
        result = analyzer.validate_domain(domain, 'test_context', sequence_length=200)
        
        assert result.is_valid == False
        assert any('sequence' in error.lower() or 'boundary' in error.lower() 
                  for error in result.errors)
    
    def test_validate_domain_invalid_coordinates(self, analyzer):
        """Test validating domain with invalid coordinates"""
        domain = DomainModel(
            id='invalid_domain',
            start=100,
            end=50,  # End before start
            range='100-50',
            confidence=0.9
        )
        
        result = analyzer.validate_domain(domain, 'test_context', sequence_length=200)
        
        assert result.is_valid == False
        assert any('coordinate' in error.lower() or 'start' in error.lower() 
                  for error in result.errors)
    
    def test_validate_domain_low_confidence(self, analyzer):
        """Test validating low confidence domain"""
        analyzer.options.min_evidence_confidence = 0.8
        
        domain = DomainModel(
            id='low_conf_domain',
            start=10,
            end=100,
            range='10-100',
            confidence=0.3  # Below threshold
        )
        
        result = analyzer.validate_domain(domain, 'test_context', sequence_length=200)
        
        # Should have warnings in normal mode
        if analyzer.options.validation_level != ValidationLevel.STRICT:
            assert result.is_valid == True
            assert len(result.warnings) > 0
        else:
            assert result.is_valid == False


class TestCacheManagement:
    """Test evidence analyzer caching"""
    
    @pytest.fixture
    def analyzer(self):
        """Create analyzer for testing"""
        options = PartitionOptions(use_cache=True)
        return EvidenceAnalyzer(options)
    
    def test_classification_cache_hit(self, analyzer):
        """Test classification cache hit"""
        # Pre-populate cache
        analyzer.classification_cache.set_domain_classification(
            'd1abcA1',
            {'t_group': '1.1.1.1', 'h_group': '1.1.1'}
        )
        
        # Lookup should hit cache
        result = analyzer.classification_cache.get_domain_classification('d1abcA1')
        
        assert result is not None
        assert result['t_group'] == '1.1.1.1'
        assert analyzer.classification_cache.cache_hits > 0
    
    def test_classification_cache_miss(self, analyzer):
        """Test classification cache miss"""
        initial_misses = analyzer.classification_cache.cache_misses
        
        result = analyzer.classification_cache.get_domain_classification('nonexistent')
        
        assert result is None
        assert analyzer.classification_cache.cache_misses > initial_misses
    
    def test_cache_statistics(self, analyzer):
        """Test cache statistics reporting"""
        # Perform some cache operations
        analyzer.classification_cache.set_domain_classification('d1', {'t_group': '1.1.1.1'})
        analyzer.classification_cache.get_domain_classification('d1')  # Hit
        analyzer.classification_cache.get_domain_classification('d2')  # Miss
        
        stats = analyzer.get_cache_statistics()
        
        assert 'cache_hits' in stats
        assert 'cache_misses' in stats
        assert 'hit_rate' in stats
        assert stats['cache_hits'] >= 1
        assert stats['cache_misses'] >= 1
    
    def test_clear_cache(self, analyzer):
        """Test cache clearing"""
        # Populate cache
        analyzer.classification_cache.set_domain_classification('d1', {'t_group': '1.1.1.1'})
        
        # Clear cache
        analyzer.clear_cache()
        
        # Should be empty
        result = analyzer.classification_cache.get_domain_classification('d1')
        assert result is None
        
        stats = analyzer.get_cache_statistics()
        assert stats['domain_classifications_cached'] == 0


class TestAnalyzerStatistics:
    """Test analyzer statistics tracking"""
    
    @pytest.fixture
    def analyzer(self):
        """Create analyzer for testing"""
        options = PartitionOptions()
        return EvidenceAnalyzer(options)
    
    def test_statistics_initialization(self, analyzer):
        """Test initial statistics state"""
        stats = analyzer.get_statistics()
        
        assert 'summaries_parsed' in stats
        assert 'evidence_extracted' in stats
        assert 'validation_errors' in stats
        assert 'classification_lookups' in stats
        assert all(v == 0 for v in stats.values() if isinstance(v, int))
    
    def test_statistics_tracking(self, analyzer):
        """Test statistics are updated during operations"""
        # Create a simple summary for parsing
        summary_xml = '''<?xml version="1.0"?>
        <domain_summary pdb_id="1abc" chain_id="A" sequence_length="100">
        </domain_summary>'''
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xml', delete=False) as f:
            f.write(summary_xml)
            f.flush()
            
            try:
                # Parse summary
                analyzer.parse_domain_summary(f.name)
                
                # Check statistics updated
                stats = analyzer.get_statistics()
                assert stats['summaries_parsed'] >= 1
                
            finally:
                os.unlink(f.name)
    
    def test_reset_statistics(self, analyzer):
        """Test resetting statistics"""
        # Generate some statistics
        analyzer.stats['summaries_parsed'] = 5
        analyzer.stats['evidence_extracted'] = 10
        
        # Reset
        analyzer.reset_statistics()
        
        # Should be back to zero
        stats = analyzer.get_statistics()
        assert stats['summaries_parsed'] == 0
        assert stats['evidence_extracted'] == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
