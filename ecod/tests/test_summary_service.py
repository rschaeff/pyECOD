#!/usr/bin/env python3
"""
Comprehensive test coverage for ecod.pipelines.domain_analysis.summary.service module

Tests the DomainSummaryService class functionality including:
- Service initialization and configuration
- Single protein processing
- Batch processing (sequential and parallel)
- File-based processing 
- Error handling and validation
- Statistics and monitoring
- Convenience functions
- Context manager behavior
"""

import pytest
import tempfile
import os
import json
from unittest.mock import Mock, MagicMock, patch, mock_open, call
from datetime import datetime
from pathlib import Path
from concurrent.futures import Future
import threading
import time

# Mock missing modules before importing
import sys
from unittest.mock import MagicMock

# Mock SelfComparisonProcessor
class MockSelfComparisonProcessor:
    def __init__(self, logger=None):
        self.logger = logger

    def process(self, file_path):
        return []

    def validate_file(self, file_path):
        return False

# Mock EvidenceQualityFilter
class MockEvidenceQualityFilter:
    def __init__(self, min_confidence=0.3):
        self.min_confidence = min_confidence

    def filter(self, evidence_collection):
        return evidence_collection

# Mock the missing modules
mock_self_comparison = MagicMock()
mock_self_comparison.SelfComparisonProcessor = MockSelfComparisonProcessor
sys.modules['ecod.pipelines.domain_analysis.summary.processors.self_comparison'] = mock_self_comparison

mock_filters = MagicMock()
mock_filters.EvidenceQualityFilter = MockEvidenceQualityFilter
sys.modules['ecod.pipelines.domain_analysis.summary.filters'] = mock_filters

from ecod.pipelines.domain_analysis.summary.service import (
    DomainSummaryService, create_service, process_single_protein
)
from ecod.pipelines.domain_analysis.summary.models import (
    ProteinIdentifier, SummaryOptions, EvidenceSummary,
    BatchResults, ProcessingStats, ProcessingStatus, EvidenceCollection
)
from ecod.pipelines.domain_analysis.summary.generator import (
    DomainSummaryGenerator, GeneratorConfig
)
from ecod.models.pipeline import DomainPartitionResult, Evidence
from ecod.exceptions import PipelineError, ValidationError


class TestDomainSummaryServiceInitialization:
    """Test service initialization and configuration"""

    @pytest.fixture
    def mock_context(self):
        """Create mock application context"""
        context = Mock()
        context.config_manager.config = {
            'reference': {'current_version': 'develop291'},
            'summary_options': {
                'blast_only': False,
                'min_confidence': 0.3
            },
            'evidence': {'min_confidence': 0.3},
            'processors': {
                'blast': {'hsp_evalue_threshold': 0.005},
                'hhsearch': {'probability_threshold': 0.0}
            }
        }
        context.config_manager.get_db_config.return_value = {
            'host': 'localhost',
            'database': 'test_db'
        }
        return context

    @pytest.fixture
    def mock_db_manager(self):
        """Create mock database manager"""
        db = Mock()
        db.test_connection.return_value = True
        return db

    def test_service_initialization_default(self, mock_context, mock_db_manager):
        """Test basic service initialization with defaults"""
        with patch('ecod.pipelines.domain_analysis.summary.service.DBManager', return_value=mock_db_manager):
            service = DomainSummaryService(mock_context)

            assert service.context == mock_context
            assert service.config == mock_context.config_manager.config
            assert service.db == mock_db_manager
            assert isinstance(service.generator, DomainSummaryGenerator)
            assert isinstance(service.service_stats, ProcessingStats)
            assert service.batch_config['max_workers'] == 1
            assert service.batch_config['use_multiprocessing'] == False

    def test_service_initialization_with_custom_config(self, mock_context, mock_db_manager):
        """Test service initialization with custom configuration"""
        service_config = {
            'max_batch_workers': 4,
            'use_multiprocessing': True,
            'batch_size': 50,
            'save_summaries': True,
            'generator': {
                'max_workers': 8,
                'cache_enabled': False,
                'parallel_evidence_collection': False
            }
        }

        with patch('ecod.pipelines.domain_analysis.summary.service.DBManager', return_value=mock_db_manager):
            service = DomainSummaryService(mock_context, service_config)

            assert service.service_config == service_config
            assert service.batch_config['max_workers'] == 4
            assert service.batch_config['use_multiprocessing'] == True
            assert service.batch_config['batch_size'] == 50

            # Check generator config
            gen_config = service.generator.generator_config
            assert gen_config.max_workers == 8
            assert gen_config.cache_enabled == False
            assert gen_config.parallel_evidence_collection == False

    def test_generator_config_creation(self, mock_context, mock_db_manager):
        """Test generator configuration creation"""
        service_config = {
            'generator': {
                'process_timeout': 600,
                'validate_inputs': False,
                'skip_on_error': True,
                'collect_detailed_stats': False
            }
        }

        with patch('ecod.pipelines.domain_analysis.summary.service.DBManager', return_value=mock_db_manager):
            service = DomainSummaryService(mock_context, service_config)

            gen_config = service.generator.generator_config
            assert gen_config.process_timeout == 600
            assert gen_config.validate_inputs == False
            assert gen_config.skip_on_error == True
            assert gen_config.collect_detailed_stats == False

    def test_database_initialization_error(self, mock_context):
        """Test handling of database initialization errors"""
        with patch('ecod.pipelines.domain_analysis.summary.service.DBManager', side_effect=Exception("DB connection failed")):
            with pytest.raises(Exception, match="DB connection failed"):
                DomainSummaryService(mock_context)


class TestSingleProteinProcessing:
    """Test single protein processing functionality"""

    @pytest.fixture
    def service(self):
        """Create service with mocked dependencies"""
        context = Mock()
        context.config_manager.config = {
            'reference': {'current_version': 'develop291'},
            'summary_options': {'blast_only': False}
        }
        context.config_manager.get_db_config.return_value = {}
        context.is_force_overwrite.return_value = True  # Return actual boolean

        with patch('ecod.pipelines.domain_analysis.summary.service.DBManager'):
            service = DomainSummaryService(context)
            service.generator = Mock()
            service.generator.initialize_for_job = Mock()  # Make it a proper Mock
            return service

    def test_process_protein_success(self, service):
        """Test successful protein processing"""
        # Mock successful evidence summary
        evidence_summary = EvidenceSummary(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291"
        )
        evidence_summary.status = ProcessingStatus.COMPLETED
        evidence_summary.evidence_collection = EvidenceCollection()
        evidence_summary.evidence_collection.add_source("hhsearch", [Evidence(type="hhsearch", confidence=0.9)])

        service.generator.generate_summary.return_value = evidence_summary

        with tempfile.TemporaryDirectory() as temp_dir:
            result = service.process_protein("1abc", "A", temp_dir)

            assert isinstance(result, DomainPartitionResult)
            assert result.pdb_id == "1abc"
            assert result.chain_id == "A"
            assert result.reference == "develop291"
            assert result.success == True

            # Check generator was called correctly
            service.generator.generate_summary.assert_called_once()
            call_args = service.generator.generate_summary.call_args
            protein = call_args[0][0]
            assert protein.pdb_id == "1abc"
            assert protein.chain_id == "A"

    def test_process_protein_with_options(self, service):
        """Test protein processing with custom options"""
        evidence_summary = EvidenceSummary("1abc", "A", "develop291")
        service.generator.generate_summary.return_value = evidence_summary

        with tempfile.TemporaryDirectory() as temp_dir:
            result = service.process_protein(
                "1abc", "A", temp_dir,
                blast_only=True,
                force_overwrite=True,
                min_confidence=0.8
            )

            # Check options were passed correctly
            call_args = service.generator.generate_summary.call_args
            options = call_args[0][1]
            assert options.blast_only == True
            assert options.force_overwrite == True
            assert options.min_confidence == 0.8

    def test_process_protein_validation_error(self, service):
        """Test handling of protein identifier validation errors"""
        # Invalid PDB ID (too short)
        with tempfile.TemporaryDirectory() as temp_dir:
            with pytest.raises(ValidationError):
                service.process_protein("abc", "A", temp_dir)

    def test_process_protein_pipeline_error(self, service):
        """Test handling of pipeline processing errors"""
        service.generator.generate_summary.side_effect = PipelineError("Processing failed")

        with tempfile.TemporaryDirectory() as temp_dir:
            result = service.process_protein("1abc", "A", temp_dir)

            assert result.success == False
            assert "Processing failed" in result.error
            assert service.service_stats.errors == 1

    def test_process_protein_unexpected_error(self, service):
        """Test handling of unexpected errors"""
        service.generator.generate_summary.side_effect = RuntimeError("Unexpected error")

        with tempfile.TemporaryDirectory() as temp_dir:
            result = service.process_protein("1abc", "A", temp_dir)

            assert result.success == False
            assert "Unexpected error" in result.error
            assert service.service_stats.errors == 1

    def test_process_protein_generator_initialization(self, service):
        """Test generator initialization during processing"""
        evidence_summary = EvidenceSummary("1abc", "A", "develop291")
        service.generator.generate_summary.return_value = evidence_summary
        service.generator.job_dump_dir = None  # Not initialized

        with tempfile.TemporaryDirectory() as temp_dir:
            result = service.process_protein("1abc", "A", temp_dir)

            # Should have initialized generator
            service.generator.initialize_for_job.assert_called_once_with(
                temp_dir, reference="develop291"
            )

    def test_process_protein_statistics_update(self, service):
        """Test statistics are updated correctly"""
        evidence_summary = EvidenceSummary("1abc", "A", "develop291")
        evidence_summary.evidence_collection = EvidenceCollection()
        evidence_summary.evidence_collection.add_source("blast", [Evidence(type="blast")])
        service.generator.generate_summary.return_value = evidence_summary

        initial_count = service.service_stats.proteins_processed

        with tempfile.TemporaryDirectory() as temp_dir:
            # Mock the conversion to return a successful result with domains
            with patch.object(service, '_convert_to_partition_result') as mock_convert:
                # Create a more realistic mock domain that won't break DomainPartitionResult
                mock_domain = Mock()
                mock_domain.confidence = 0.8
                mock_domain.evidence = []  # Empty list instead of Mock

                mock_result = DomainPartitionResult(
                    pdb_id="1abc", chain_id="A", reference="develop291",
                    success=True, domains=[mock_domain]
                )
                mock_convert.return_value = mock_result

                result = service.process_protein("1abc", "A", temp_dir)

                assert service.service_stats.proteins_processed == initial_count + 1
                assert service.service_stats.proteins_with_evidence == 1


class TestBatchProcessing:
    """Test batch processing functionality"""

    @pytest.fixture
    def service(self):
        """Create service for batch testing"""
        context = Mock()
        context.config_manager.config = {
            'reference': {'current_version': 'develop291'},
            'summary_options': {}
        }
        context.config_manager.get_db_config.return_value = {}

        with patch('ecod.pipelines.domain_analysis.summary.service.DBManager'):
            service = DomainSummaryService(context)
            service.generator = Mock()
            service.generator.initialize_for_job = Mock()
            return service

    def test_process_batch_sequential_success(self, service):
        """Test successful sequential batch processing"""
        proteins = [("1abc", "A"), ("2xyz", "B"), ("3def", "C")]

        # Mock successful summaries
        def mock_generate_summary(protein, options):
            summary = EvidenceSummary(protein.pdb_id, protein.chain_id, "develop291")
            summary.evidence_collection = EvidenceCollection()
            summary.evidence_collection.add_source("blast", [Evidence(type="blast")])
            return summary

        service.generator.generate_summary.side_effect = mock_generate_summary

        with tempfile.TemporaryDirectory() as temp_dir:
            results = service.process_batch(proteins, temp_dir)

            assert results.total == 3
            assert results.success_count == 3
            assert results.failure_count == 0
            assert len(results.successes) == 3

            # Check all proteins were processed
            assert service.generator.generate_summary.call_count == 3

    def test_process_batch_sequential_mixed_results(self, service):
        """Test sequential batch with mixed success/failure"""
        proteins = [("1abc", "A"), ("2xyz", "B"), ("3def", "C")]

        def mock_generate_summary(protein, options):
            if protein.pdb_id == "2xyz":
                raise Exception("Processing failed")

            summary = EvidenceSummary(protein.pdb_id, protein.chain_id, "develop291")
            summary.evidence_collection = EvidenceCollection()
            return summary

        service.generator.generate_summary.side_effect = mock_generate_summary

        with tempfile.TemporaryDirectory() as temp_dir:
            results = service.process_batch(proteins, temp_dir)

            assert results.total == 3
            assert results.success_count == 2
            assert results.failure_count == 1

            # Check failure was recorded
            assert len(results.failures) == 1
            assert results.failures[0] == ("2xyz", "B", "Processing failed")

    def test_process_batch_parallel_success(self, service):
        """Test successful parallel batch processing"""
        # Configure for multiprocessing
        service.batch_config['use_multiprocessing'] = True
        service.batch_config['max_workers'] = 2

        proteins = [("1abc", "A"), ("2xyz", "B")]

        # Mock ProcessPoolExecutor
        mock_future1 = Mock()
        mock_future1.result.return_value = {
            'pdb_id': '1abc',
            'chain_id': 'A',
            'reference': 'develop291',
            'status': 'COMPLETED',
            'total_evidence': 1,
            'metadata': {}
        }

        mock_future2 = Mock()
        mock_future2.result.return_value = {
            'pdb_id': '2xyz',
            'chain_id': 'B',
            'reference': 'develop291',
            'status': 'COMPLETED',
            'total_evidence': 1,
            'metadata': {}
        }

        mock_executor = MagicMock()
        mock_executor.__enter__ = Mock(return_value=mock_executor)
        mock_executor.__exit__ = Mock(return_value=None)
        mock_executor.submit.side_effect = [mock_future1, mock_future2]

        with patch('ecod.pipelines.domain_analysis.summary.service.ProcessPoolExecutor', return_value=mock_executor):
            with patch('ecod.pipelines.domain_analysis.summary.service.as_completed', return_value=[mock_future1, mock_future2]):
                with tempfile.TemporaryDirectory() as temp_dir:
                    results = service.process_batch(proteins, temp_dir)

                    assert results.total == 2
                    assert results.success_count == 2
                    assert results.failure_count == 0

    def test_process_batch_parallel_with_errors(self, service):
        """Test parallel batch processing with worker errors"""
        service.batch_config['use_multiprocessing'] = True
        service.batch_config['max_workers'] = 2

        proteins = [("1abc", "A"), ("2xyz", "B")]

        # Mock futures - one success, one failure
        mock_future1 = Mock()
        mock_future1.result.return_value = {
            'pdb_id': '1abc',
            'chain_id': 'A',
            'reference': 'develop291',  # Add missing reference field
            'status': 'COMPLETED',
            'total_evidence': 1,
            'processing_time': 1.0,
            'metadata': {}
        }

        mock_future2 = Mock()
        mock_future2.result.side_effect = Exception("Worker failed")

        mock_executor = MagicMock()
        mock_executor.__enter__ = Mock(return_value=mock_executor)
        mock_executor.__exit__ = Mock(return_value=None)
        mock_executor.submit.side_effect = [mock_future1, mock_future2]

        # Mock the mapping of futures to proteins
        future_to_protein = {mock_future1: ("1abc", "A"), mock_future2: ("2xyz", "B")}

        with patch('ecod.pipelines.domain_analysis.summary.service.ProcessPoolExecutor', return_value=mock_executor):
            with patch('ecod.pipelines.domain_analysis.summary.service.as_completed', return_value=[mock_future1, mock_future2]):
                with tempfile.TemporaryDirectory() as temp_dir:
                    results = service.process_batch(proteins, temp_dir)

                    assert results.total == 2
                    assert results.success_count == 1
                    assert results.failure_count == 1

    def test_process_batch_empty_list(self, service):
        """Test batch processing with empty protein list"""
        with tempfile.TemporaryDirectory() as temp_dir:
            results = service.process_batch([], temp_dir)

            assert results.total == 0
            assert results.success_count == 0
            assert results.failure_count == 0

    def test_process_batch_progress_logging(self, service):
        """Test batch processing logs progress correctly"""
        # Create batch of 15 proteins (should log at 10)
        proteins = [(f"{i:04d}", "A") for i in range(15)]

        def mock_generate_summary(protein, options):
            return EvidenceSummary(protein.pdb_id, protein.chain_id, "develop291")

        service.generator.generate_summary.side_effect = mock_generate_summary

        with patch.object(service.logger, 'info') as mock_log:
            with tempfile.TemporaryDirectory() as temp_dir:
                results = service.process_batch(proteins, temp_dir)

                # Should log progress at 10 proteins processed
                progress_calls = [call for call in mock_log.call_args_list
                                if "Processed 10/" in str(call)]
                assert len(progress_calls) >= 1


class TestFileBasedProcessing:
    """Test file-based protein processing"""

    @pytest.fixture
    def service(self):
        """Create service for file testing"""
        context = Mock()
        context.config_manager.config = {
            'reference': {'current_version': 'develop291'}
        }
        context.config_manager.get_db_config.return_value = {}

        with patch('ecod.pipelines.domain_analysis.summary.service.DBManager'):
            service = DomainSummaryService(context)
            service.generator = Mock()
            service.generator.initialize_for_job = Mock()
            return service

    def test_process_proteins_from_file_success(self, service):
        """Test processing proteins from file"""
        protein_list_content = """1abc_A
2xyz_B
3def_C
# This is a comment
4ghi_D
"""

        def mock_generate_summary(protein, options):
            return EvidenceSummary(protein.pdb_id, protein.chain_id, "develop291")

        service.generator.generate_summary.side_effect = mock_generate_summary

        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            f.write(protein_list_content)
            protein_file = f.name

        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                results = service.process_proteins_from_file(protein_file, temp_dir)

                assert results.total == 4  # Should skip comment
                assert results.success_count == 4
                assert results.failure_count == 0
        finally:
            os.unlink(protein_file)

    def test_process_proteins_from_file_with_invalid_entries(self, service):
        """Test file processing with invalid protein identifiers"""
        protein_list_content = """1abc_A
invalid_protein_id
2xyz_B
abc_X
3def_C
"""

        def mock_generate_summary(protein, options):
            return EvidenceSummary(protein.pdb_id, protein.chain_id, "develop291")

        service.generator.generate_summary.side_effect = mock_generate_summary

        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            f.write(protein_list_content)
            protein_file = f.name

        try:
            with patch.object(service.logger, 'warning') as mock_warning:
                with tempfile.TemporaryDirectory() as temp_dir:
                    results = service.process_proteins_from_file(protein_file, temp_dir)

                    # Should process valid entries and warn about invalid ones
                    assert results.total == 3  # 1abc_A, 2xyz_B, 3def_C
                    assert results.success_count == 3

                    # Should log warnings for invalid entries
                    warning_calls = mock_warning.call_args_list
                    assert len(warning_calls) == 2  # invalid_protein_id and abc_X
        finally:
            os.unlink(protein_file)

    def test_process_proteins_from_file_not_found(self, service):
        """Test processing from non-existent file"""
        with tempfile.TemporaryDirectory() as temp_dir:
            results = service.process_proteins_from_file("/nonexistent/file.txt", temp_dir)

            assert results.total == 0

    def test_process_proteins_from_file_empty(self, service):
        """Test processing from empty file"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            f.write("")  # Empty file
            protein_file = f.name

        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                results = service.process_proteins_from_file(protein_file, temp_dir)

                assert results.total == 0
        finally:
            os.unlink(protein_file)

    def test_read_protein_list_various_formats(self, service):
        """Test reading protein list with various identifier formats"""
        protein_list_content = """1abc_A
2xyz:B
3def.C
4ghi-D
5jklE
"""

        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            f.write(protein_list_content)
            protein_file = f.name

        try:
            proteins = service._read_protein_list(protein_file)

            assert len(proteins) == 5
            assert ("1abc", "A") in proteins
            assert ("2xyz", "B") in proteins
            assert ("3def", "C") in proteins
            assert ("4ghi", "D") in proteins
            assert ("5jkl", "E") in proteins
        finally:
            os.unlink(protein_file)


class TestHelperMethods:
    """Test helper methods and utilities"""

    @pytest.fixture
    def service(self):
        """Create service for helper method testing"""
        context = Mock()
        context.config_manager.config = {
            'reference': {'current_version': 'develop291'},
            'summary_options': {
                'blast_only': False,
                'min_confidence': 0.3
            }
        }
        context.config_manager.get_db_config.return_value = {}
        context.is_force_overwrite.return_value = True

        with patch('ecod.pipelines.domain_analysis.summary.service.DBManager'):
            service = DomainSummaryService(context)
            service.generator = Mock()
            service.generator.initialize_for_job = Mock()
            return service

    def test_get_reference(self, service):
        """Test reference version retrieval"""
        reference = service._get_reference()
        assert reference == "develop291"

    def test_create_summary_options_defaults(self, service):
        """Test summary options creation with defaults"""
        options = service._create_summary_options()

        assert options.blast_only == False
        assert options.min_confidence == 0.3
        assert options.force_overwrite == True  # From context

    def test_create_summary_options_with_overrides(self, service):
        """Test summary options creation with overrides"""
        options = service._create_summary_options(
            blast_only=True,
            min_confidence=0.8,
            stitch_discontinuous=False
        )

        assert options.blast_only == True
        assert options.min_confidence == 0.8
        assert options.stitch_discontinuous == False
        assert options.force_overwrite == True  # Still from context

    def test_ensure_generator_initialized_new_dir(self, service):
        """Test generator initialization for new job directory"""
        service.generator.job_dump_dir = None

        with tempfile.TemporaryDirectory() as temp_dir:
            service._ensure_generator_initialized(temp_dir)

            service.generator.initialize_for_job.assert_called_once_with(
                temp_dir, reference="develop291"
            )

    def test_ensure_generator_initialized_same_dir(self, service):
        """Test generator not re-initialized for same directory"""
        with tempfile.TemporaryDirectory() as temp_dir:
            service.generator.job_dump_dir = Path(temp_dir)

            service._ensure_generator_initialized(temp_dir)

            # Should not call initialize again
            service.generator.initialize_for_job.assert_not_called()

    def test_convert_to_partition_result(self, service):
        """Test conversion of evidence summary to partition result"""
        evidence_summary = EvidenceSummary("1abc", "A", "develop291")
        evidence_summary.evidence_collection = EvidenceCollection()

        with tempfile.TemporaryDirectory() as temp_dir:
            options = SummaryOptions(save_intermediate=True)
            result = service._convert_to_partition_result(evidence_summary, temp_dir, options)

            assert isinstance(result, DomainPartitionResult)
            assert result.pdb_id == "1abc"
            assert result.chain_id == "A"
            assert result.reference == "develop291"

            # Should set domain file path
            assert result.domain_file is not None
            assert "1abc_A.develop291.evidence_summary.xml" in result.domain_file

    def test_convert_to_partition_result_blast_only(self, service):
        """Test conversion with blast_only option"""
        evidence_summary = EvidenceSummary("1abc", "A", "develop291")

        with tempfile.TemporaryDirectory() as temp_dir:
            options = SummaryOptions(save_intermediate=True, blast_only=True)
            result = service._convert_to_partition_result(evidence_summary, temp_dir, options)

            # Should include blast_only in filename
            assert ".blast_only." in result.domain_file

    def test_reconstruct_summary(self, service):
        """Test reconstruction of summary from dictionary"""
        summary_dict = {
            'pdb_id': '1abc',
            'chain_id': 'A',
            'reference': 'develop291',
            'status': 'COMPLETED',
            'total_evidence': 5,
            'processing_time': 1.5,
            'metadata': {'sequence_length': 200}
        }

        summary = service._reconstruct_summary(summary_dict)

        assert isinstance(summary, EvidenceSummary)
        assert summary.pdb_id == "1abc"
        assert summary.chain_id == "A"
        assert summary.reference == "develop291"
        assert summary.status == ProcessingStatus.COMPLETED
        assert summary.processing_time == 1.5
        assert summary.metadata['sequence_length'] == 200


class TestStatisticsAndMonitoring:
    """Test statistics collection and monitoring"""

    @pytest.fixture
    def service(self):
        """Create service for statistics testing"""
        context = Mock()
        context.config_manager.config = {
            'reference': {'current_version': 'develop291'}
        }
        context.config_manager.get_db_config.return_value = {}

        with patch('ecod.pipelines.domain_analysis.summary.service.DBManager'):
            service = DomainSummaryService(context)
            service.generator = Mock()
            service.generator.processors = {'blast': Mock(), 'hhsearch': Mock()}  # Mock with length
            service.generator.get_statistics.return_value = {
                'proteins_processed': 10,
                'proteins_with_evidence': 8,
                'evidence_by_source': {'blast': 15, 'hhsearch': 12}
            }
            return service

    def test_get_service_statistics(self, service):
        """Test service statistics collection"""
        stats = service.get_service_statistics()

        assert 'service_stats' in stats
        assert 'generator_stats' in stats
        assert 'configuration' in stats

        # Check service stats structure
        service_stats = stats['service_stats']
        assert 'proteins_processed' in service_stats

        # Check generator stats
        generator_stats = stats['generator_stats']
        assert generator_stats['proteins_processed'] == 10

        # Check configuration
        config = stats['configuration']
        assert config['reference'] == 'develop291'
        assert 'batch_config' in config
        assert 'generator_config' in config

    def test_validate_setup_success(self, service):
        """Test successful setup validation"""
        # Mock successful database connection
        service.db.test_connection.return_value = True

        validation = service.validate_setup()

        assert validation['database'] == True
        assert validation['processors'] == True  # Generator has processors
        assert validation['config_loaded'] == True
        assert validation['reference_set'] == True

    def test_validate_setup_database_failure(self, service):
        """Test setup validation with database failure"""
        service.db.test_connection.side_effect = Exception("Connection failed")

        validation = service.validate_setup()

        assert validation['database'] == False
        # Other validations should still work
        assert validation['config_loaded'] == True

    def test_clear_all_caches(self, service):
        """Test cache clearing"""
        service.generator.clear_caches = Mock()
        service.clear_all_caches()

        # Should call clear on generator
        service.generator.clear_caches.assert_called_once()


class TestWorkerFunction:
    """Test multiprocessing worker function"""

    def test_process_single_protein_worker(self):
        """Test worker function for multiprocessing"""
        # Instead of trying to test the internal implementation,
        # let's test that the worker function can be called without errors
        # and mock the entire chain

        with patch('ecod.pipelines.domain_analysis.summary.service.ApplicationContext') as mock_context_class:
            with patch('ecod.pipelines.domain_analysis.summary.service.DBManager'):
                with patch('ecod.pipelines.domain_analysis.summary.service.DomainSummaryGenerator') as mock_gen_class:
                    # Mock context
                    mock_context = Mock()
                    mock_context.config_manager.config = {'reference': {'current_version': 'develop291'}}
                    mock_context.config_manager.get_db_config.return_value = {}
                    mock_context_class.return_value = mock_context

                    # Mock generator
                    mock_generator = Mock()
                    mock_gen_class.return_value = mock_generator

                    # Mock evidence summary
                    mock_summary = Mock()
                    mock_summary.get_summary_dict.return_value = {
                        'pdb_id': '1abc',
                        'chain_id': 'A',
                        'status': 'COMPLETED'
                    }
                    mock_generator.generate_summary.return_value = mock_summary

                    # Test the worker function
                    from ecod.pipelines.domain_analysis.summary.service import DomainSummaryService
                    result = DomainSummaryService._process_single_protein_worker(
                        "1abc", "A", "/job/dir", {'blast_only': True}
                    )

                    assert result['pdb_id'] == '1abc'
                    assert result['chain_id'] == 'A'
                    assert result['status'] == 'COMPLETED'


class TestContextManager:
    """Test context manager functionality"""

    @pytest.fixture
    def service(self):
        """Create service for context manager testing"""
        context = Mock()
        context.config_manager.config = {'reference': {'current_version': 'develop291'}}
        context.config_manager.get_db_config.return_value = {}

        with patch('ecod.pipelines.domain_analysis.summary.service.DBManager'):
            service = DomainSummaryService(context)
            service.generator = Mock()
            service.generator.clear_caches = Mock()
            return service

    def test_context_manager_enter_exit(self, service):
        """Test context manager entry and exit"""
        with service as s:
            assert s is service

        # Should have called cleanup methods
        service.generator.clear_caches.assert_called_once()

    def test_context_manager_with_db_close(self, service):
        """Test context manager cleanup with database close"""
        service.db.close = Mock()
        service.generator.clear_caches = Mock()

        with service as s:
            pass

        service.db.close.assert_called_once()
        service.generator.clear_caches.assert_called_once()


class TestConvenienceFunctions:
    """Test convenience functions"""

    @patch('ecod.pipelines.domain_analysis.summary.service.ApplicationContext')
    @patch('ecod.pipelines.domain_analysis.summary.service.DomainSummaryService')
    def test_create_service_with_config_path(self, mock_service_class, mock_context_class):
        """Test create_service with explicit config path"""
        mock_context = Mock()
        mock_context_class.return_value = mock_context

        mock_service = Mock()
        mock_service_class.return_value = mock_service

        result = create_service("/path/to/config.yml")

        assert result == mock_service
        mock_context_class.assert_called_once_with("/path/to/config.yml")
        mock_service_class.assert_called_once_with(mock_context)

    @patch.dict(os.environ, {'ECOD_CONFIG_PATH': '/env/config.yml'})
    @patch('ecod.pipelines.domain_analysis.summary.service.ApplicationContext')
    @patch('ecod.pipelines.domain_analysis.summary.service.DomainSummaryService')
    def test_create_service_from_environment(self, mock_service_class, mock_context_class):
        """Test create_service using environment variable"""
        mock_context = Mock()
        mock_context_class.return_value = mock_context

        mock_service = Mock()
        mock_service_class.return_value = mock_service

        result = create_service()

        assert result == mock_service
        mock_context_class.assert_called_once_with("/env/config.yml")

    @patch('ecod.pipelines.domain_analysis.summary.service.create_service')
    def test_process_single_protein_convenience(self, mock_create_service):
        """Test process_single_protein convenience function"""
        mock_service = Mock()
        mock_create_service.return_value = mock_service

        mock_result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291"
        )
        mock_service.process_protein.return_value = mock_result

        result = process_single_protein(
            "1abc", "A", "/job/dir",
            config_path="/config.yml",
            blast_only=True
        )

        assert result == mock_result
        mock_create_service.assert_called_once_with("/config.yml")
        mock_service.process_protein.assert_called_once_with(
            "1abc", "A", "/job/dir", blast_only=True
        )


class TestErrorHandling:
    """Test comprehensive error handling"""

    @pytest.fixture
    def service(self):
        """Create service for error testing"""
        context = Mock()
        context.config_manager.config = {
            'reference': {'current_version': 'develop291'}
        }
        context.config_manager.get_db_config.return_value = {}

        with patch('ecod.pipelines.domain_analysis.summary.service.DBManager'):
            service = DomainSummaryService(context)
            service.generator = Mock()
            service.generator.initialize_for_job = Mock()
            return service

    def test_invalid_protein_identifier_formats(self, service):
        """Test various invalid protein identifier formats"""
        invalid_identifiers = [
            ("", "A"),          # Empty PDB ID
            ("1ab", "A"),       # Too short PDB ID
            ("1abcd", "A"),     # Too long PDB ID
            ("1abc", ""),       # Empty chain ID
            ("1abc", "ABCDE"),  # Too long chain ID
        ]

        for pdb_id, chain_id in invalid_identifiers:
            with tempfile.TemporaryDirectory() as temp_dir:
                with pytest.raises(ValidationError):
                    service.process_protein(pdb_id, chain_id, temp_dir)

    def test_generator_initialization_failure(self, service):
        """Test handling of generator initialization failure"""
        service.generator.initialize_for_job.side_effect = Exception("Init failed")

        with tempfile.TemporaryDirectory() as temp_dir:
            # Service catches the exception and returns error result
            result = service.process_protein("1abc", "A", temp_dir)

            assert result.success == False
            assert "Init failed" in result.error

    def test_batch_processing_timeout(self, service):
        """Test batch processing with timeouts"""
        proteins = [("1abc", "A")]

        # Mock slow generator
        def slow_generator(protein, options):
            time.sleep(0.1)  # Simulate slow processing
            return EvidenceSummary(protein.pdb_id, protein.chain_id, "develop291")

        service.generator.generate_summary.side_effect = slow_generator

        with tempfile.TemporaryDirectory() as temp_dir:
            # Should complete normally (timeout handling is in ProcessPoolExecutor)
            results = service.process_batch(proteins, temp_dir)
            assert results.total == 1

    def test_memory_constraints_large_batch(self, service):
        """Test handling of large batch sizes"""
        # Create a large batch (but not so large it actually causes memory issues in test)
        proteins = [(f"{i:04d}", "A") for i in range(100)]

        def mock_generate_summary(protein, options):
            return EvidenceSummary(protein.pdb_id, protein.chain_id, "develop291")

        service.generator.generate_summary.side_effect = mock_generate_summary

        with tempfile.TemporaryDirectory() as temp_dir:
            results = service.process_batch(proteins, temp_dir)

            assert results.total == 100
            assert results.success_count == 100


class TestIntegrationScenarios:
    """Test integration scenarios that combine multiple features"""

    @pytest.fixture
    def service(self):
        """Create service for integration testing"""
        context = Mock()
        context.config_manager.config = {
            'reference': {'current_version': 'develop291'},
            'summary_options': {
                'blast_only': False,
                'min_confidence': 0.3
            }
        }
        context.config_manager.get_db_config.return_value = {}

        with patch('ecod.pipelines.domain_analysis.summary.service.DBManager'):
            service = DomainSummaryService(context)
            service.generator = Mock()
            service.generator.initialize_for_job = Mock()
            return service

    def test_full_pipeline_blast_only_mode(self, service):
        """Test complete pipeline in blast-only mode"""
        # Mock evidence summary with only BLAST evidence
        evidence_summary = EvidenceSummary("1abc", "A", "develop291")
        evidence_collection = EvidenceCollection()
        evidence_collection.add_source("blast", [
            Evidence(type="blast", evalue=1e-5, confidence=0.8)
        ])
        evidence_summary.evidence_collection = evidence_collection
        service.generator.generate_summary.return_value = evidence_summary

        with tempfile.TemporaryDirectory() as temp_dir:
            result = service.process_protein(
                "1abc", "A", temp_dir,
                blast_only=True,
                save_intermediate=True
            )

            assert result.success == True
            assert result.pdb_id == "1abc"
            
            # Check options were set correctly
            call_args = service.generator.generate_summary.call_args
            options = call_args[0][1]
            assert options.blast_only == True

    def test_full_pipeline_with_evidence_filtering(self, service):
        """Test pipeline with evidence quality filtering"""
        evidence_summary = EvidenceSummary("1abc", "A", "develop291")
        evidence_collection = EvidenceCollection()
        evidence_collection.add_source("hhsearch", [
            Evidence(type="hhsearch", probability=95.0, confidence=0.95),
            Evidence(type="hhsearch", probability=10.0, confidence=0.1)  # Low quality
        ])
        evidence_summary.evidence_collection = evidence_collection
        service.generator.generate_summary.return_value = evidence_summary

        with tempfile.TemporaryDirectory() as temp_dir:
            result = service.process_protein(
                "1abc", "A", temp_dir,
                min_confidence=0.5,  # Should filter out low quality evidence
                skip_filtering=False
            )

            # Check filtering options were applied
            call_args = service.generator.generate_summary.call_args
            options = call_args[0][1]
            assert options.min_confidence == 0.5
            assert options.skip_filtering == False

    def test_batch_processing_with_file_output(self, service):
        """Test batch processing with file output generation"""
        proteins = [("1abc", "A"), ("2xyz", "B")]

        def mock_generate_summary(protein, options):
            summary = EvidenceSummary(protein.pdb_id, protein.chain_id, "develop291")
            summary.evidence_collection = EvidenceCollection()
            return summary

        service.generator.generate_summary.side_effect = mock_generate_summary

        with tempfile.TemporaryDirectory() as temp_dir:
            results = service.process_batch(
                proteins, temp_dir,
                save_intermediate=True,
                include_evidence_details=True
            )

            assert results.total == 2
            assert results.success_count == 2

            # Check that domain file paths would be set (in actual implementation)
            # This is verified through the convert_to_partition_result method

    def test_error_recovery_in_batch_processing(self, service):
        """Test error recovery and continuation in batch processing"""
        proteins = [("1abc", "A"), ("fail", "B"), ("2xyz", "C")]

        def mock_generate_summary(protein, options):
            if protein.pdb_id == "fail":
                raise Exception("Transient failure")
            
            summary = EvidenceSummary(protein.pdb_id, protein.chain_id, "develop291")
            summary.evidence_collection = EvidenceCollection()
            return summary

        service.generator.generate_summary.side_effect = mock_generate_summary

        with tempfile.TemporaryDirectory() as temp_dir:
            results = service.process_batch(proteins, temp_dir)

            # Should continue processing despite one failure
            assert results.total == 3
            assert results.success_count == 2
            assert results.failure_count == 1
            
            # Failed protein should be recorded
            assert len(results.failures) == 1
            assert results.failures[0][0] == "fail"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
