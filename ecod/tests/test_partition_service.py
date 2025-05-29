#!/usr/bin/env python3
"""
Critical tests for domain partition service and processor.

Focus on testing the core algorithm logic that's likely to change during retunes,
service orchestration, and error handling that could break the pipeline.
"""

import pytest
import tempfile
import os
from unittest.mock import Mock, MagicMock, patch, mock_open
from datetime import datetime
from pathlib import Path

# Import the models and service components
from ecod.pipelines.domain_analysis.partition.models import (
    PartitionOptions, ValidationLevel, ProcessingMode, PartitionContext,
    ValidationResult, EvidenceGroup, DomainCandidate, BatchPartitionResults
)
from ecod.pipelines.domain_analysis.partition.processor import PartitionProcessor
from ecod.pipelines.domain_analysis.partition.service import DomainPartitionService
from ecod.pipelines.domain_analysis.partition.tracker import StatusTracker

from ecod.models.pipeline.evidence import Evidence
from ecod.models.pipeline.domain import DomainModel
from ecod.models.pipeline.partition import DomainPartitionResult


class TestPartitionOptions:
    """Test configuration validation and options handling"""
    
    def test_default_options_creation(self):
        """Test creating default partition options"""
        options = PartitionOptions()
        
        assert options.mode == ProcessingMode.STANDARD
        assert options.validation_level == ValidationLevel.NORMAL
        assert options.min_domain_size == 20
        assert options.overlap_threshold == 0.3
        assert options.merge_gap_tolerance == 20
        assert options.use_cache == True
        assert options.blast_only == False
    
    def test_options_validation_success(self):
        """Test valid options pass validation"""
        options = PartitionOptions(
            min_domain_size=30,
            max_domain_size=500,
            overlap_threshold=0.4,
            min_evidence_confidence=0.8
        )
        
        # Should not raise exception
        options.validate()
    
    def test_options_validation_failures(self):
        """Test invalid options fail validation"""
        # min_domain_size too small
        with pytest.raises(ValueError, match="min_domain_size must be at least 1"):
            options = PartitionOptions(min_domain_size=0)
            options.validate()
        
        # max_domain_size smaller than min
        with pytest.raises(ValueError, match="max_domain_size must be >= min_domain_size"):
            options = PartitionOptions(min_domain_size=100, max_domain_size=50)
            options.validate()
        
        # Invalid overlap threshold
        with pytest.raises(ValueError, match="overlap_threshold must be between 0.0 and 1.0"):
            options = PartitionOptions(overlap_threshold=1.5)
            options.validate()
        
        # Invalid confidence threshold
        with pytest.raises(ValueError, match="min_evidence_confidence must be between 0.0 and 1.0"):
            options = PartitionOptions(min_evidence_confidence=-0.1)
            options.validate()
    
    def test_options_serialization(self):
        """Test options to/from dict conversion"""
        original = PartitionOptions(
            mode=ProcessingMode.BLAST_ONLY,
            validation_level=ValidationLevel.STRICT,
            min_domain_size=25,
            overlap_threshold=0.25,
            blast_only=True
        )
        
        # Convert to dict and back
        options_dict = original.to_dict()
        reconstructed = PartitionOptions.from_dict(options_dict)
        
        assert reconstructed.mode == ProcessingMode.BLAST_ONLY
        assert reconstructed.validation_level == ValidationLevel.STRICT
        assert reconstructed.min_domain_size == 25
        assert reconstructed.overlap_threshold == 0.25
        assert reconstructed.blast_only == True


class TestPartitionProcessor:
    """Test core domain identification algorithms"""
    
    @pytest.fixture
    def mock_analyzer(self):
        """Create mock evidence analyzer"""
        analyzer = Mock()
        analyzer.group_evidence_by_position.return_value = {}
        analyzer.validate_domain.return_value = ValidationResult(is_valid=True)
        return analyzer
    
    @pytest.fixture
    def processor(self, mock_analyzer):
        """Create processor with mock analyzer"""
        options = PartitionOptions(min_domain_size=20, overlap_threshold=0.3)
        return PartitionProcessor(options, mock_analyzer)
    
    @pytest.fixture
    def sample_context(self):
        """Create sample partition context"""
        return PartitionContext(
            pdb_id="1abc",
            chain_id="A", 
            reference="develop291",
            sequence_length=200
        )
    
    def test_process_evidence_empty_list(self, processor, sample_context):
        """Test processing with no evidence"""
        result = processor.process_evidence([], sample_context)
        
        assert result.success == True
        assert result.is_unclassified == True
        assert len(result.domains) == 0
        assert result.pdb_id == "1abc"
        assert result.chain_id == "A"
    
    def test_process_evidence_peptide(self, processor):
        """Test processing peptide (too short)"""
        context = PartitionContext(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291", 
            sequence_length=15  # Below min_domain_size
        )
        
        evidence = [Evidence(type="blast", confidence=0.8)]
        result = processor.process_evidence(evidence, context)
        
        assert result.success == True
        assert result.is_peptide == True
        assert result.is_unclassified == True
        assert len(result.domains) == 0
    
    def test_process_evidence_with_domains(self, processor, sample_context, mock_analyzer):
        """Test processing evidence that results in domains"""
        # Create mock evidence
        evidence1 = Evidence(
            type="domain_blast",
            query_range="10-50",
            confidence=0.8,
            domain_id="d1abcA1"
        )
        evidence2 = Evidence(
            type="hhsearch", 
            query_range="80-120",
            confidence=0.9,
            domain_id="d2xyzB1"
        )
        
        # Mock evidence grouping
        group1 = EvidenceGroup(evidence_items=[evidence1])
        group1.consensus_start = 10
        group1.consensus_end = 50
        group1.consensus_confidence = 0.8
        
        group2 = EvidenceGroup(evidence_items=[evidence2])
        group2.consensus_start = 80
        group2.consensus_end = 120
        group2.consensus_confidence = 0.9
        
        mock_analyzer.group_evidence_by_position.return_value = {
            100: group1,  # Position key 
            200: group2
        }
        
        result = processor.process_evidence([evidence1, evidence2], sample_context)
        
        assert result.success == True
        assert result.is_classified == True
        assert len(result.domains) == 2
        assert result.domains[0].start == 10
        assert result.domains[0].end == 50
        assert result.domains[1].start == 80 
        assert result.domains[1].end == 120
    
    def test_resolve_domain_overlaps_no_overlap(self, processor):
        """Test overlap resolution with non-overlapping domains"""
        candidates = [
            DomainCandidate(start=10, end=30, evidence_group=EvidenceGroup(), confidence=0.8),
            DomainCandidate(start=50, end=80, evidence_group=EvidenceGroup(), confidence=0.9)
        ]
        
        resolved = processor.resolve_domain_overlaps(candidates, 100)
        
        assert len(resolved) == 2
        assert resolved[0].start == 10
        assert resolved[1].start == 50
    
    def test_resolve_domain_overlaps_with_overlap(self, processor):
        """Test overlap resolution with overlapping domains"""
        candidates = [
            DomainCandidate(start=10, end=40, evidence_group=EvidenceGroup(), 
                          confidence=0.7, protected=False),
            DomainCandidate(start=30, end=60, evidence_group=EvidenceGroup(), 
                          confidence=0.9, protected=False)  # Higher confidence wins
        ]
        
        resolved = processor.resolve_domain_overlaps(candidates, 100)
        
        # Should keep the higher confidence domain
        assert len(resolved) == 1
        assert resolved[0].start == 30
        assert resolved[0].end == 60
        assert resolved[0].confidence == 0.9
    
    def test_resolve_domain_overlaps_protected_domain(self, processor):
        """Test that protected domains are always kept"""
        candidates = [
            DomainCandidate(start=10, end=40, evidence_group=EvidenceGroup(),
                          confidence=0.9, protected=True),  # Protected but lower confidence
            DomainCandidate(start=30, end=60, evidence_group=EvidenceGroup(),
                          confidence=0.95, protected=False)  # Higher confidence
        ]
        
        resolved = processor.resolve_domain_overlaps(candidates, 100)
        
        # Both should be kept because one is protected
        assert len(resolved) == 2
    
    def test_assign_domain_classifications(self, processor):
        """Test domain classification assignment"""
        evidence = Evidence(
            type="hhsearch",
            confidence=0.9,
            t_group="2002.1.1.1",
            h_group="2002.1.1", 
            x_group="2002.1",
            a_group="a.39"
        )
        
        domain = DomainModel(
            id="test_domain",
            start=10,
            end=50,
            range="10-50",
            evidence=[evidence]
        )
        
        processor.assign_domain_classifications([domain])
        
        # Domain should get classification from evidence
        assert domain.t_group == "2002.1.1.1"
        assert domain.h_group == "2002.1.1"
        assert domain.x_group == "2002.1"
        assert domain.a_group == "a.39"
    
    def test_merge_nearby_candidates(self, processor):
        """Test merging candidates that are close together"""
        candidates = [
            DomainCandidate(start=10, end=30, evidence_group=EvidenceGroup(),
                          confidence=0.8, protected=False),
            DomainCandidate(start=35, end=55, evidence_group=EvidenceGroup(),
                          confidence=0.7, protected=False)  # Gap of 4, within tolerance
        ]
        
        # Set merge tolerance to allow this merge
        processor.options.merge_gap_tolerance = 10
        
        merged = processor._merge_nearby_candidates(candidates)
        
        # Should merge into one domain
        assert len(merged) == 1
        assert merged[0].start == 10
        assert merged[0].end == 55
    
    def test_processor_statistics(self, processor):
        """Test processor statistics tracking"""
        stats = processor.get_statistics()
        
        # Should have expected stat keys
        expected_keys = [
            'domains_identified', 'overlaps_resolved', 'classifications_assigned',
            'domains_validated', 'discontinuous_domains', 'cache_stats'
        ]
        
        for key in expected_keys:
            assert key in stats
    
    def test_processor_cache_clearing(self, processor):
        """Test processor cache clearing"""
        # Should not raise exception
        processor.clear_cache()
        processor.reset_statistics()
        
        # Stats should be reset
        stats = processor.get_statistics()
        assert stats['domains_identified'] == 0


class TestDomainPartitionService:
    """Test high-level service orchestration"""
    
    @pytest.fixture
    def mock_context(self):
        """Create mock application context"""
        context = Mock()
        context.config_manager.config = {
            'partition': {
                'validation_level': 'normal',
                'min_domain_size': 20,
                'overlap_threshold': 0.3
            },
            'reference': {
                'current_version': 'develop291'
            }
        }
        context.config_manager.get_db_config.return_value = {}
        return context
    
    @pytest.fixture 
    def mock_db(self):
        """Create mock database manager"""
        db = Mock()
        db.test_connection.return_value = True
        db.execute_dict_query.return_value = []
        return db
    
    @pytest.fixture
    def service(self, mock_context, mock_db):
        """Create service with mocked dependencies"""
        with patch('ecod.pipelines.domain_analysis.partition.service.DBManager', return_value=mock_db):
            service = DomainPartitionService(mock_context)
            return service
    
    def test_service_initialization(self, service):
        """Test service initializes correctly"""
        assert service.default_options is not None
        assert service.analyzer is not None
        assert service.processor is not None
        assert service.tracker is not None
        assert 'proteins_processed' in service.service_stats
    
    def test_create_default_options_from_config(self, service):
        """Test options creation from config"""
        options = service.default_options
        
        assert options.validation_level == ValidationLevel.NORMAL
        assert options.min_domain_size == 20
        assert options.overlap_threshold == 0.3
    
    @patch('ecod.pipelines.domain_analysis.partition.service.os.path.exists')
    def test_partition_protein_peptide(self, mock_exists, service):
        """Test partitioning a peptide protein"""
        mock_exists.return_value = True

        with patch.object(service.analyzer, 'parse_domain_summary') as mock_parse:
            mock_parse.return_value = {
                "is_peptide": True,
                "sequence_length": 15
            }

            with tempfile.TemporaryDirectory() as temp_dir:
                result = service.partition_protein(
                    "1abc", "A", "/fake/summary.xml", temp_dir
                )
        
        assert result.success == True
        assert result.is_peptide == True
        assert result.is_unclassified == True
        assert len(result.domains) == 0
    
    @patch('ecod.pipelines.domain_analysis.partition.service.os.path.exists')
    def test_partition_protein_with_domains(self, mock_exists, service):
        """Test partitioning protein that results in domains"""
        mock_exists.return_value = True

        with patch.object(service.analyzer, 'parse_domain_summary') as mock_parse, \
             patch.object(service.analyzer, 'extract_evidence_with_classification') as mock_extract, \
             patch.object(service.processor, 'process_evidence') as mock_process:

            mock_parse.return_value = {
                "sequence_length": 200,
                "is_peptide": False
            }

            evidence = Evidence(type="hhsearch", confidence=0.9, query_range="10-50")
            mock_extract.return_value = [evidence]

            # Mock processor returns
            domain = DomainModel(id="test_d1", start=10, end=50, range="10-50")
            mock_result = DomainPartitionResult(
                pdb_id="1abc", chain_id="A", reference="develop291",
                domains=[domain], success=True, is_classified=True
            )
            mock_process.return_value = mock_result

            with tempfile.TemporaryDirectory() as temp_dir:
                result = service.partition_protein(
                    "1abc", "A", "/fake/summary.xml", temp_dir
                )

        assert result.success == True
        assert result.is_classified == True
        assert len(result.domains) == 1
    
    def test_partition_protein_error_handling(self, service):
        """Test service handles errors gracefully"""
        with patch.object(service.analyzer, 'parse_domain_summary') as mock_parse:
            mock_parse.side_effect = Exception("Parse error")

            with tempfile.TemporaryDirectory() as temp_dir:
                result = service.partition_protein(
                    "1abc", "A", "/fake/summary.xml", temp_dir
                )

        assert result.success == False
        assert "Parse error" in result.error
    
    def test_get_proteins_to_process(self, service):
        """Test getting proteins to process from database"""
        # Mock database query
        service.db.execute_dict_query.return_value = [
            {
                'protein_id': 1,
                'pdb_id': '1abc',
                'chain_id': 'A',
                'process_id': 101,
                'is_representative': True
            },
            {
                'protein_id': 2,
                'pdb_id': '2xyz', 
                'chain_id': 'B',
                'process_id': 102,
                'is_representative': False
            }
        ]
        
        proteins = service._get_proteins_to_process(batch_id=123, limit=10)
        
        assert len(proteins) == 2
        assert proteins[0]['pdb_id'] == '1abc'
        assert proteins[1]['pdb_id'] == '2xyz'
    
    def test_batch_processing_empty(self, service):
        """Test batch processing with no proteins"""
        service.db.execute_dict_query.return_value = []
        
        results = service.partition_batch(123, "/fake/batch/path")
        
        assert results.total == 0
        assert results.success_count == 0
        assert results.failure_count == 0
    
    def test_service_statistics(self, service):
        """Test service statistics collection"""
        stats = service.get_service_statistics()
        
        assert 'service' in stats
        assert 'processor' in stats
        assert 'analyzer' in stats
        assert 'configuration' in stats
        
        # Service stats
        service_stats = stats['service']
        assert 'proteins_processed' in service_stats
        assert 'runtime_seconds' in service_stats
        assert 'proteins_per_minute' in service_stats
    
    def test_validate_setup(self, service):
        """Test service setup validation"""
        # Mock successful validation
        service.db.test_connection.return_value = True
        
        validation = service.validate_setup()
        
        assert 'database' in validation
        assert 'config_loaded' in validation
        assert 'reference_set' in validation
        assert 'options_valid' in validation
        
        # Should all be True if setup is valid
        assert validation['database'] == True
        assert validation['config_loaded'] == True
        assert validation['reference_set'] == True
        assert validation['options_valid'] == True
    
    def test_validate_setup_failures(self, service):
        """Test service validation with failures"""
        # Mock database failure
        service.db.test_connection.side_effect = Exception("DB Error")
        
        validation = service.validate_setup()
        
        assert validation['database'] == False
    
    def test_clear_all_caches(self, service):
        """Test clearing all service caches"""
        # Should not raise exception
        service.clear_all_caches()
    
    def test_context_manager(self, service):
        """Test service as context manager"""
        with service as s:
            assert s is service
        
        # Should clean up without error


class TestStatusTracker:
    """Test database status tracking"""
    
    @pytest.fixture
    def mock_db(self):
        """Create mock database"""
        db = Mock()
        db.update.return_value = True
        db.insert.return_value = 1
        db.execute_query.return_value = []
        db.execute_dict_query.return_value = []
        return db
    
    @pytest.fixture
    def tracker(self, mock_db):
        """Create status tracker with mock DB"""
        return StatusTracker(mock_db)
    
    def test_update_process_status_success(self, tracker, mock_db):
        """Test successful process status update"""
        success = tracker.update_process_status(
            process_id=123,
            stage="domain_partition_processing",
            status="processing"
        )
        
        assert success == True
        mock_db.update.assert_called_once()
        
        # Check call arguments
        call_args = mock_db.update.call_args
        assert "ecod_schema.process_status" in call_args[0]
        assert call_args[0][1]["current_stage"] == "domain_partition_processing"
        assert call_args[0][1]["status"] == "processing"
    
    def test_update_process_status_with_error(self, tracker, mock_db):
        """Test process status update with error message"""
        error_msg = "Processing failed due to invalid evidence"
        
        success = tracker.update_process_status(
            process_id=123,
            stage="domain_partition_failed", 
            status="error",
            error_message=error_msg
        )
        
        assert success == True
        call_args = mock_db.update.call_args
        assert call_args[0][1]["error_message"] == error_msg
    
    def test_update_process_status_long_error(self, tracker, mock_db):
        """Test error message truncation"""
        long_error = "x" * 1000  # Very long error message
        
        tracker.update_process_status(
            process_id=123,
            stage="domain_partition_failed",
            status="error", 
            error_message=long_error
        )
        
        call_args = mock_db.update.call_args
        stored_error = call_args[0][1]["error_message"]
        assert len(stored_error) == 500  # Truncated
    
    def test_update_process_status_no_process_id(self, tracker):
        """Test handling of missing process ID"""
        success = tracker.update_process_status(
            process_id=None,
            stage="test",
            status="test"
        )
        
        assert success == False
    
    def test_update_process_status_db_error(self, tracker, mock_db):
        """Test handling of database errors"""
        mock_db.update.side_effect = Exception("DB Error")
        
        success = tracker.update_process_status(
            process_id=123,
            stage="test",
            status="test"
        )
        
        assert success == False
    
    @patch('os.path.exists')
    @patch('os.path.getsize')
    def test_register_domain_file_new(self, mock_getsize, mock_exists, tracker, mock_db):
        """Test registering new domain file"""
        mock_exists.return_value = True
        mock_getsize.return_value = 1024
        mock_db.execute_query.return_value = []  # No existing record
        
        success = tracker.register_domain_file(
            process_id=123,
            file_path="/path/to/domains/1abc_A.develop291.domains.xml",
            base_path="/path/to"
        )
        
        assert success == True
        mock_db.insert.assert_called_once()
        
        # Check insert data
        call_args = mock_db.insert.call_args
        assert "ecod_schema.process_file" in call_args[0]
        assert call_args[0][1]["file_type"] == "domain_partition"
        assert call_args[0][1]["file_exists"] == True
        assert call_args[0][1]["file_size"] == 1024
    
    @patch('os.path.exists')
    def test_register_domain_file_update_existing(self, mock_exists, tracker, mock_db):
        """Test updating existing domain file record"""
        mock_exists.return_value = True
        mock_db.execute_query.return_value = [(456,)]  # Existing record ID
        
        with patch('os.path.getsize', return_value=2048):
            success = tracker.register_domain_file(
                process_id=123,
                file_path="/path/to/domains/1abc_A.develop291.domains.xml",
                base_path="/path/to"
            )
        
        assert success == True
        mock_db.update.assert_called_once()
    
    def test_get_batch_progress(self, tracker, mock_db):
        """Test getting batch progress information"""
        mock_db.execute_dict_query.return_value = [{
            'total': 100,
            'complete': 85,
            'errors': 10,
            'processing': 5,
            'skipped': 0,
            'representatives': 60,
            'files_created': 85
        }]
        
        progress = tracker.get_batch_progress(batch_id=123)
        
        assert progress['total'] == 100
        assert progress['complete'] == 85
        assert progress['errors'] == 10
        assert progress['complete_pct'] == 85.0
        assert progress['error_pct'] == 10.0
        assert progress['progress_pct'] == 95.0  # (85 + 10) / 100
    
    def test_get_batch_progress_empty(self, tracker, mock_db):
        """Test getting progress for empty batch"""
        mock_db.execute_dict_query.return_value = []
        
        progress = tracker.get_batch_progress(batch_id=123)
        
        assert progress['total'] == 0
        assert progress['complete_pct'] == 0
        assert progress['error_pct'] == 0
    
    def test_update_batch_completion_status(self, tracker, mock_db):
        """Test updating batch completion status"""
        # Mock query result showing batch completion
        mock_db.execute_dict_query.return_value = [{
            'total': 100,
            'complete': 95,
            'errors': 5
        }]
        
        success = tracker.update_batch_completion_status(batch_id=123)
        
        assert success == True
        mock_db.update.assert_called_once()
        
        # Check status was set correctly
        call_args = mock_db.update.call_args
        assert call_args[0][1]["status"] == "domain_partition_complete_with_errors"
        assert call_args[0][1]["completed_items"] == 95


class TestValidationResult:
    """Test validation result handling"""
    
    def test_validation_result_creation(self):
        """Test creating validation results"""
        result = ValidationResult(is_valid=True, context="test_context")
        
        assert result.is_valid == True
        assert result.context == "test_context"
        assert len(result.errors) == 0
        assert len(result.warnings) == 0
    
    def test_add_error_marks_invalid(self):
        """Test adding error marks result as invalid"""
        result = ValidationResult(is_valid=True)
        result.add_error("Test error message")
        
        assert result.is_valid == False
        assert "Test error message" in result.errors
    
    def test_add_warning_keeps_valid(self):
        """Test adding warning doesn't affect validity"""
        result = ValidationResult(is_valid=True)
        result.add_warning("Test warning")
        
        assert result.is_valid == True
        assert "Test warning" in result.warnings
    
    def test_merge_validation_results(self):
        """Test merging two validation results"""
        result1 = ValidationResult(is_valid=True)
        result1.add_warning("Warning 1")
        
        result2 = ValidationResult(is_valid=False)
        result2.add_error("Error 1")
        result2.add_warning("Warning 2")
        
        result1.merge(result2)
        
        assert result1.is_valid == False  # Should be false due to merge
        assert "Error 1" in result1.errors
        assert "Warning 1" in result1.warnings
        assert "Warning 2" in result1.warnings
    
    def test_get_summary(self):
        """Test validation summary generation"""
        result = ValidationResult(is_valid=False, context="test")
        result.add_error("Error 1")
        result.add_warning("Warning 1")
        
        summary = result.get_summary()
        
        assert "Context: test" in summary
        assert "INVALID" in summary
        assert "Errors: 1" in summary
        assert "Warnings: 1" in summary


class TestBatchPartitionResults:
    """Test batch processing results"""
    
    def test_batch_results_creation(self):
        """Test creating batch results"""
        results = BatchPartitionResults()
        
        assert results.total == 0
        assert results.success_count == 0
        assert results.failure_count == 0
        assert results.success_rate == 0.0
    
    def test_add_successful_result(self):
        """Test adding successful partition result"""
        results = BatchPartitionResults()
        
        domain = DomainModel(id="d1", start=10, end=50, range="10-50")
        result = DomainPartitionResult(
            pdb_id="1abc", chain_id="A", reference="develop291",
            domains=[domain], success=True, is_classified=True
        )
        
        results.add_result(result)
        
        assert results.total == 1
        assert results.success_count == 1
        assert results.failure_count == 0
        assert results.proteins_with_domains == 1
        assert results.total_domains_found == 1
        assert results.success_rate == 100.0
    
    def test_add_failed_result(self):
        """Test adding failed partition result"""
        results = BatchPartitionResults()
        
        result = DomainPartitionResult(
            pdb_id="1abc", chain_id="A", reference="develop291",
            success=False, error="Processing failed"
        )
        
        results.add_result(result)
        
        assert results.total == 1
        assert results.success_count == 0
        assert results.failure_count == 1
        assert results.success_rate == 0.0
        assert len(results.failures) == 1
        assert results.failures[0] == ("1abc", "A", "Processing failed")
    
    def test_add_peptide_result(self):
        """Test adding peptide result"""
        results = BatchPartitionResults()
        
        result = DomainPartitionResult(
            pdb_id="1abc", chain_id="A", reference="develop291",
            success=True, is_peptide=True
        )
        
        results.add_result(result)
        
        assert results.peptides_found == 1
        assert results.proteins_with_domains == 0
    
    def test_batch_results_summary(self):
        """Test batch results summary generation"""
        results = BatchPartitionResults()
        
        # Add various types of results
        results.add_result(DomainPartitionResult(
            pdb_id="1abc", chain_id="A", reference="develop291",
            domains=[DomainModel(id="d1", start=10, end=50, range="10-50")],
            success=True, is_classified=True
        ))
        
        results.add_result(DomainPartitionResult(
            pdb_id="1xyz", chain_id="B", reference="develop291", 
            success=True, is_peptide=True
        ))
        
        results.add_result(DomainPartitionResult(
            pdb_id="2fail", chain_id="C", reference="develop291",
            success=False, error="Error"
        ))
        
        results.finalize()
        summary = results.get_summary()
        
        assert summary['total_proteins'] == 3
        assert summary['successful'] == 2
        assert summary['failed'] == 1
        assert summary['proteins_with_domains'] == 1
        assert summary['peptides'] == 1
        assert summary['total_domains'] == 1


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
