#!/usr/bin/env python3
"""
Comprehensive tests for domain analysis partition components.

Tests the actual implementation: EvidenceAnalyzer, PartitionProcessor,
DomainPartitionService, StatusTracker, and supporting classes.
"""

import pytest
import tempfile
import os
import json
import xml.etree.ElementTree as ET
from unittest.mock import Mock, MagicMock, patch, mock_open
from datetime import datetime, timedelta
from pathlib import Path
from collections import defaultdict

# Import actual classes from implementation
from ecod.pipelines.domain_analysis.partition.analyzer import (
    EvidenceAnalyzer, ClassificationCache, EvidenceQualityMetrics
)
from ecod.pipelines.domain_analysis.partition.processor import (
    PartitionProcessor, DiscontinuousDomainCandidate
)
from ecod.pipelines.domain_analysis.partition.service import (
    DomainPartitionService, create_service, partition_single_protein
)
from ecod.pipelines.domain_analysis.partition.tracker import (
    StatusTracker, ProcessStage, ProcessStatus, ProcessInfo, BatchInfo
)
from ecod.pipelines.domain_analysis.partition.reference_analyzer import (
    ReferenceCoverageAnalyzer, ReferenceInfo
)
from ecod.pipelines.domain_analysis.partition.models import (
    PartitionOptions, ValidationLevel, ProcessingMode, PartitionContext,
    ValidationResult, EvidenceGroup, DomainCandidate, BatchPartitionResults,
    PartitionStage, EvidenceWithCoverage
)

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
        assert 0.0 <= options.overlap_threshold <= 1.0
        assert options.use_cache == True
        assert options.blast_only == False

    def test_options_validation_success(self):
        """Test valid options pass validation"""
        options = PartitionOptions(
            min_domain_size=30,
            overlap_threshold=0.4,
            min_evidence_confidence=0.8,
            min_reference_coverage=0.8
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

    def test_factory_methods(self):
        """Test factory methods for common configurations"""
        blast_only = PartitionOptions.create_blast_only()
        assert blast_only.blast_only == True
        assert blast_only.use_hhsearch == False
        assert blast_only.mode == ProcessingMode.BLAST_ONLY

        strict = PartitionOptions.create_strict()
        assert strict.validation_level == ValidationLevel.STRICT
        assert strict.require_evidence == True
        assert strict.min_evidence_confidence == 0.8

    def test_options_serialization(self):
        """Test options to/from dict conversion"""
        original = PartitionOptions(
            mode=ProcessingMode.BLAST_ONLY,
            validation_level=ValidationLevel.STRICT,
            min_domain_size=25,
            overlap_threshold=0.25
        )

        # Convert to dict and back
        options_dict = original.to_dict()
        reconstructed = PartitionOptions.from_dict(options_dict)

        assert reconstructed.mode == ProcessingMode.BLAST_ONLY
        assert reconstructed.validation_level == ValidationLevel.STRICT
        assert reconstructed.min_domain_size == 25
        assert reconstructed.overlap_threshold == 0.25


class TestClassificationCache:
    """Test the classification cache system"""

    def test_cache_creation(self):
        """Test cache initialization"""
        cache = ClassificationCache(max_size=100, ttl_hours=12)

        assert cache.max_size == 100
        assert cache.ttl == timedelta(hours=12)
        assert cache.hits == 0
        assert cache.misses == 0

    def test_cache_set_get(self):
        """Test basic cache operations"""
        cache = ClassificationCache()

        classification = {
            "t_group": "2004.1.1.1",
            "h_group": "2004.1.1",
            "x_group": "2004",
            "a_group": "a.39"
        }

        # Set classification
        success = cache.set_domain_classification("d1abcA1", classification)
        assert success == True

        # Get classification
        result = cache.get_domain_classification("d1abcA1")
        assert result is not None
        assert result["t_group"] == "2004.1.1.1"
        assert cache.hits == 1
        assert cache.misses == 0

    def test_cache_miss(self):
        """Test cache miss"""
        cache = ClassificationCache()

        result = cache.get_domain_classification("nonexistent")
        assert result is None
        assert cache.misses == 1

    def test_cache_stats(self):
        """Test cache statistics"""
        cache = ClassificationCache()

        # Add some data
        cache.set_domain_classification("d1", {"t_group": "1.1.1.1"})
        cache.get_domain_classification("d1")  # hit
        cache.get_domain_classification("d2")  # miss

        stats = cache.get_stats()
        assert stats['size'] == 1
        assert stats['hits'] == 1
        assert stats['misses'] == 1
        assert stats['hit_rate_percent'] == 50.0


class TestEvidenceAnalyzer:
    """Test the main evidence analyzer"""

    @pytest.fixture
    def options(self):
        """Create test options"""
        return PartitionOptions(
            validation_level=ValidationLevel.NORMAL,
            min_evidence_confidence=0.1
        )

    @pytest.fixture
    def analyzer(self, options):
        """Create analyzer with mock database"""
        mock_db = Mock()
        analyzer = EvidenceAnalyzer(options, mock_db)
        return analyzer

    def test_analyzer_initialization(self, analyzer):
        """Test analyzer initializes correctly"""
        assert analyzer.options is not None
        assert analyzer.logger is not None
        assert analyzer.classification_cache is not None
        assert len(analyzer.range_patterns) > 0

    def test_parse_domain_summary_file_not_found(self, analyzer):
        """Test parsing non-existent file"""
        result = analyzer.parse_domain_summary("/nonexistent/file.xml")

        assert "error" in result
        assert "not found" in result["error"]

    @patch('os.path.exists')
    @patch('os.access')
    @patch('os.path.getsize')
    def test_parse_domain_summary_permission_denied(self, mock_size, mock_access, mock_exists, analyzer):
        """Test parsing file with no read permission"""
        mock_exists.return_value = True
        mock_access.return_value = False
        mock_size.return_value = 1024

        result = analyzer.parse_domain_summary("/no/permission.xml")

        assert "error" in result
        assert "Permission denied" in result["error"]

    def test_parse_range_comprehensive(self, analyzer):
        """Test range parsing with different formats"""
        # Standard range
        ranges = analyzer._parse_range_comprehensive("10-50")
        assert ranges == [(10, 50)]

        # Multiple ranges
        ranges = analyzer._parse_range_comprehensive("10-50,80-120")
        assert ranges == [(10, 50), (80, 120)]

        # Different separators
        ranges = analyzer._parse_range_comprehensive("10:50")
        assert ranges == [(10, 50)]

        # Invalid range
        ranges = analyzer._parse_range_comprehensive("invalid")
        assert ranges == []

        # Empty string
        ranges = analyzer._parse_range_comprehensive("")
        assert ranges == []

    def test_create_evidence_from_blast(self, analyzer):
        """Test creating Evidence from BLAST hit data"""
        hit_data = {
            "type": "domain_blast",
            "domain_id": "d1abcA1",
            "query_range": "10-50",
            "hit_range": "1-40",
            "evalue": 1e-5,
            "pdb_id": "1abc",
            "chain_id": "A"
        }

        evidence = analyzer._create_evidence_from_blast(hit_data)

        assert evidence is not None
        assert evidence.type == "domain_blast"
        assert evidence.domain_id == "d1abcA1"
        assert evidence.query_range == "10-50"
        assert evidence.evalue == 1e-5

    def test_create_evidence_from_hhsearch(self, analyzer):
        """Test creating Evidence from HHSearch hit data"""
        hit_data = {
            "type": "hhsearch",
            "domain_id": "d1abcA1",
            "query_range": "10-50",
            "hit_range": "1-40",
            "probability": 95.5,
            "evalue": 1e-10,
            "score": 150.2
        }

        evidence = analyzer._create_evidence_from_hhsearch(hit_data)

        assert evidence is not None
        assert evidence.type == "hhsearch"
        assert evidence.domain_id == "d1abcA1"
        assert evidence.probability == 95.5
        assert evidence.evalue == 1e-10

    def test_validate_evidence(self, analyzer):
        """Test evidence validation"""
        # Valid HHSearch evidence
        evidence = Evidence(
            type="hhsearch",
            domain_id="d1abcA1",
            probability=95.0,
            evalue=1e-10,
            confidence=0.95
        )

        result = analyzer.validate_evidence(evidence, "test_context")
        assert result.is_valid == True
        assert len(result.errors) == 0

    def test_validate_evidence_failures(self, analyzer):
        """Test evidence validation failures"""
        # Evidence with no type
        evidence = Evidence(type="", domain_id="d1abcA1")

        result = analyzer.validate_evidence(evidence, "test_context")
        assert result.is_valid == False
        assert any("empty" in error.lower() for error in result.errors)

    def test_filter_evidence(self, analyzer):
        """Test evidence filtering"""
        evidence_list = [
            Evidence(type="hhsearch", probability=95.0, confidence=0.9),  # Good
            Evidence(type="hhsearch", probability=5.0, confidence=0.05),   # Low probability
            Evidence(type="blast", evalue=1e-5, confidence=0.8),          # Good
            Evidence(type="blast", evalue=50.0, confidence=0.1),          # High e-value
        ]

        # Create matching validation results
        validation_results = [
            ValidationResult(is_valid=True),
            ValidationResult(is_valid=True),
            ValidationResult(is_valid=True),
            ValidationResult(is_valid=True)
        ]

        filtered = analyzer.filter_evidence(evidence_list, validation_results)

        # Should filter out low probability and high e-value evidence
        assert len(filtered) == 2
        assert all(e.confidence >= 0.8 for e in filtered)

    def test_calculate_quality_metrics(self, analyzer):
        """Test quality metrics calculation"""
        evidence_list = [
            Evidence(type="hhsearch", confidence=0.9),
            Evidence(type="blast", confidence=0.8),
            Evidence(type="self_comparison", confidence=0.7)
        ]

        metrics = analyzer.calculate_quality_metrics(evidence_list, 200, 1.5)

        assert metrics.total_evidence_count == 3
        assert metrics.valid_evidence_count == 3
        assert metrics.hhsearch_hit_count == 1
        assert metrics.blast_hit_count == 1
        assert metrics.processing_time_seconds == 1.5
        assert 0.0 <= metrics.average_confidence <= 1.0


class TestPartitionProcessor:
    """Test core domain partitioning algorithms"""

    @pytest.fixture
    def options(self):
        """Create test options"""
        return PartitionOptions(
            min_domain_size=20,
            overlap_threshold=0.3,
            merge_gap_tolerance=10
        )

    @pytest.fixture
    def analyzer(self, options):
        """Create mock analyzer"""
        analyzer = Mock()
        analyzer.validate_domain.return_value = ValidationResult(is_valid=True)
        return analyzer

    @pytest.fixture
    def processor(self, options, analyzer):
        """Create processor"""
        return PartitionProcessor(options, analyzer)

    @pytest.fixture
    def context(self):
        """Create test context"""
        return PartitionContext(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            sequence_length=200
        )

    def test_processor_initialization(self, processor):
        """Test processor initializes correctly"""
        assert processor.options is not None
        assert processor.analyzer is not None
        assert processor.classification_cache is not None
        assert isinstance(processor.stats, dict)

    def test_process_evidence_empty(self, processor, context):
        """Test processing with no evidence"""
        result = processor.process_evidence([], context)

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

    def test_resolve_domain_overlaps_no_overlap(self, processor):
        """Test overlap resolution with non-overlapping domains"""
        candidates = [
            DomainCandidate(
                start=10, end=30,
                evidence_group=EvidenceGroup(),
                confidence=0.8, protected=False
            ),
            DomainCandidate(
                start=50, end=80,
                evidence_group=EvidenceGroup(),
                confidence=0.9, protected=False
            )
        ]

        resolved = processor.resolve_domain_overlaps(candidates, 100)

        assert len(resolved) == 2
        assert resolved[0].start == 10
        assert resolved[1].start == 50

    def test_resolve_domain_overlaps_with_overlap(self, processor):
        """Test overlap resolution with overlapping domains"""
        candidates = [
            DomainCandidate(
                start=10, end=40,
                evidence_group=EvidenceGroup(),
                confidence=0.7, protected=False
            ),
            DomainCandidate(
                start=30, end=60,
                evidence_group=EvidenceGroup(),
                confidence=0.9, protected=False  # Higher confidence
            )
        ]

        resolved = processor.resolve_domain_overlaps(candidates, 100)

        # Should keep the higher confidence domain
        assert len(resolved) == 1
        assert resolved[0].start == 30
        assert resolved[0].end == 60
        assert resolved[0].confidence == 0.9

    def test_resolve_domain_overlaps_protected(self, processor):
        """Test that protected domains are always kept"""
        candidates = [
            DomainCandidate(
                start=10, end=40,
                evidence_group=EvidenceGroup(),
                confidence=0.8, protected=True  # Protected
            ),
            DomainCandidate(
                start=30, end=60,
                evidence_group=EvidenceGroup(),
                confidence=0.9, protected=False  # Higher confidence but not protected
            )
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
            DomainCandidate(
                start=10, end=30,
                evidence_group=EvidenceGroup(evidence_items=[Evidence(type="blast")]),
                confidence=0.8, protected=False
            ),
            DomainCandidate(
                start=35, end=55,
                evidence_group=EvidenceGroup(evidence_items=[Evidence(type="blast")]),
                confidence=0.7, protected=False  # Gap of 4, within tolerance
            )
        ]

        merged = processor._merge_nearby_candidates(candidates)

        # Should merge into one domain
        assert len(merged) == 1
        assert merged[0].start == 10
        assert merged[0].end == 55

    def test_discontinuous_domain_candidate(self):
        """Test discontinuous domain candidate"""
        segments = [(10, 50), (80, 120), (150, 180)]
        candidate = DiscontinuousDomainCandidate(
            segments=segments,
            evidence_group=EvidenceGroup(),
            confidence=0.8
        )

        assert candidate.start == 10  # First segment start
        assert candidate.end == 180   # Last segment end
        assert candidate.size == 41 + 41 + 31  # Sum of segment sizes
        assert candidate.range == "10-50,80-120,150-180"

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


class TestReferenceCoverageAnalyzer:
    """Test reference coverage analysis"""

    @pytest.fixture
    def mock_db(self):
        """Create mock database"""
        db = Mock()
        return db

    @pytest.fixture
    def analyzer(self, mock_db):
        """Create reference coverage analyzer"""
        return ReferenceCoverageAnalyzer(mock_db, min_coverage=0.7)

    def test_analyzer_initialization(self, analyzer):
        """Test analyzer initializes correctly"""
        assert analyzer.min_coverage == 0.7
        assert analyzer._reference_cache == {}

    def test_analyze_evidence_coverage_no_domain(self, analyzer):
        """Test coverage analysis without domain_id"""
        evidence = Evidence(type="blast", query_range="10-50")

        enhanced = analyzer.analyze_evidence_coverage(evidence)

        assert isinstance(enhanced, EvidenceWithCoverage)
        assert enhanced.reference_coverage == 0.0
        assert "No domain_id" in enhanced.coverage_warning

    def test_get_reference_info_cache_hit(self, analyzer):
        """Test reference info caching"""
        # Populate cache
        ref_info = ReferenceInfo(
            domain_id="d1abcA1",
            domain_range="10-150",
            domain_length=141,
            pdb_id="1abc",
            chain_id="A"
        )
        analyzer._reference_cache["d1abcA1"] = ref_info

        result = analyzer.get_reference_info("d1abcA1")

        assert result == ref_info
        assert analyzer._cache_stats['hits'] == 1

    def test_parse_range_segments(self, analyzer):
        """Test range segment parsing"""
        # Single range
        segments = analyzer._parse_range_segments("10-50")
        assert segments == [(10, 50)]

        # Multiple ranges
        segments = analyzer._parse_range_segments("10-50,80-120")
        assert segments == [(10, 50), (80, 120)]

        # Empty string
        segments = analyzer._parse_range_segments("")
        assert segments == []

    def test_suggest_extended_boundaries(self, analyzer):
        """Test boundary extension suggestions"""
        ref_info = ReferenceInfo(
            domain_id="d1abcA1",
            domain_range="1-200",
            domain_length=200,
            pdb_id="1abc",
            chain_id="A"
        )

        evidence = EvidenceWithCoverage(
            type="hhsearch",
            domain_id="d1abcA1",
            reference_info=ref_info
        )

        # Current assignment is too small
        new_start, new_end = analyzer.suggest_extended_boundaries(
            evidence, 50, 100, 300  # Current: 51 residues, expected: 200
        )

        # Should extend to get closer to expected length
        assert new_start < 50  # Extended left
        assert new_end > 100   # Extended right
        assert (new_end - new_start + 1) > 51  # Larger than original


class TestDomainPartitionService:
    """Test high-level service interface"""

    @pytest.fixture
    def mock_context(self):
        """Create mock application context"""
        context = Mock()
        context.config_manager.config = {
            'partition': {
                'validation_level': 'normal',
                'min_domain_size': 20,
                'overlap_threshold': 0.3,
                'reference_coverage': {
                    'min_coverage': 0.7,
                    'strict_coverage': 0.9
                }
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
        assert options.min_reference_coverage == 0.7

    @patch('os.path.exists')
    def test_partition_protein_file_not_found(self, mock_exists, service):
        """Test partitioning with missing summary file"""
        mock_exists.return_value = False

        with patch.object(service.analyzer, 'parse_domain_summary') as mock_parse:
            mock_parse.return_value = {"error": "File not found"}

            with tempfile.TemporaryDirectory() as temp_dir:
                result = service.partition_protein(
                    "1abc", "A", "/fake/summary.xml", temp_dir
                )

        assert result.success == False
        assert "File not found" in result.error

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

    def test_context_manager(self, service):
        """Test service as context manager"""
        with service as s:
            assert s is service

        # Should clean up without error


class TestStatusTracker:
    """Test process status tracking"""

    @pytest.fixture
    def mock_db(self):
        """Create mock database"""
        db = Mock()
        db.execute_query.return_value = None
        db.execute_dict_query.return_value = []
        return db

    @pytest.fixture
    def tracker(self, mock_db):
        """Create status tracker with mock DB"""
        return StatusTracker(mock_db)

    def test_tracker_initialization(self, tracker):
        """Test tracker initializes correctly"""
        assert tracker.db_manager is not None
        assert isinstance(tracker.processes, dict)
        assert isinstance(tracker.batches, dict)

    def test_start_process(self, tracker):
        """Test starting process tracking"""
        success = tracker.start_process(
            process_id=123,
            pdb_id="1abc",
            chain_id="A",
            batch_id=456
        )

        assert success == True
        assert 123 in tracker.processes

        process_info = tracker.processes[123]
        assert process_info.pdb_id == "1abc"
        assert process_info.chain_id == "A"
        assert process_info.batch_id == 456

    def test_update_process_status_success(self, tracker):
        """Test updating process status"""
        # First start a process
        tracker.start_process(123, "1abc", "A")

        # Update status
        success = tracker.update_process_status(
            process_id=123,
            stage="domain_partition_processing",
            status="processing"
        )

        assert success == True

        # Check the process was updated
        process_info = tracker.processes[123]
        assert process_info.stage == ProcessStage.DOMAIN_PARTITIONING
        assert process_info.status == ProcessStatus.RUNNING

    def test_update_process_status_invalid_values(self, tracker):
        """Test handling invalid stage/status values"""
        tracker.start_process(123, "1abc", "A")

        success = tracker.update_process_status(
            process_id=123,
            stage="invalid_stage",
            status="invalid_status"
        )

        assert success == False

    def test_add_process_error(self, tracker):
        """Test adding error to process"""
        tracker.start_process(123, "1abc", "A")

        success = tracker.add_process_error(
            process_id=123,
            error_message="Processing failed",
            retry=False
        )

        assert success == True

        process_info = tracker.processes[123]
        assert process_info.error_message == "Processing failed"
        assert process_info.status == ProcessStatus.FAILED

    def test_add_process_error_with_retry(self, tracker):
        """Test adding error with retry"""
        tracker.start_process(123, "1abc", "A")

        success = tracker.add_process_error(
            process_id=123,
            error_message="Temporary failure",
            retry=True
        )

        assert success == True

        process_info = tracker.processes[123]
        assert process_info.retry_count == 1
        assert process_info.status == ProcessStatus.RETRYING

    def test_register_domain_file(self, tracker):
        """Test registering domain file"""
        tracker.start_process(123, "1abc", "A")

        success = tracker.register_domain_file(
            process_id=123,
            file_path="/path/to/domains/1abc_A.xml",
            base_path="/path/to"
        )

        assert success == True

        process_info = tracker.processes[123]
        assert "/path/to/domains/1abc_A.xml" in process_info.files_created

    def test_get_process_status(self, tracker):
        """Test getting process status"""
        tracker.start_process(123, "1abc", "A", batch_id=456)

        status = tracker.get_process_status(123)

        assert status is not None
        assert status['process_id'] == 123
        assert status['pdb_id'] == "1abc"
        assert status['chain_id'] == "A"
        assert status['batch_id'] == 456

    def test_get_batch_status(self, tracker):
        """Test getting batch status"""
        # Start several processes in same batch
        tracker.start_process(123, "1abc", "A", batch_id=456)
        tracker.start_process(124, "1abc", "B", batch_id=456)
        tracker.start_process(125, "2xyz", "A", batch_id=456)

        # Complete one process
        tracker.update_process_stage(123, ProcessStage.COMPLETED, ProcessStatus.COMPLETED)

        status = tracker.get_batch_status(456)

        assert status is not None
        assert status['batch_id'] == 456
        assert status['total_processes'] == 3
        assert status['completed'] == 1
        assert status['running'] == 2

    def test_performance_metrics(self, tracker):
        """Test performance metrics collection"""
        metrics = tracker.get_performance_metrics()

        assert isinstance(metrics, dict)
        assert 'processes_per_hour' in metrics
        assert 'error_rate_percent' in metrics

    def test_cleanup_old_processes(self, tracker):
        """Test cleaning up old processes"""
        # Add a process and mark it as old
        tracker.start_process(123, "1abc", "A")
        process_info = tracker.processes[123]
        process_info.status = ProcessStatus.COMPLETED
        process_info.last_updated = datetime.now() - timedelta(days=2)

        # Cleanup processes older than 1 day
        cleaned = tracker.cleanup_old_processes(timedelta(days=1))

        assert cleaned == 1
        assert 123 not in tracker.processes


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
            success=True, is_peptide=True, is_unclassified=True
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
            success=True, is_peptide=True, is_unclassified=True
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


class TestConvenienceFunctions:
    """Test convenience functions"""

    @patch('ecod.pipelines.domain_analysis.partition.service.ApplicationContext')
    def test_create_service(self, mock_context_class):
        """Test create_service function"""
        mock_context = Mock()
        mock_context_class.return_value = mock_context

        with patch('ecod.pipelines.domain_analysis.partition.service.DomainPartitionService') as mock_service_class:
            service = create_service('/path/to/config.yml')

            mock_context_class.assert_called_once_with('/path/to/config.yml')
            mock_service_class.assert_called_once_with(mock_context, None)

    @patch('ecod.pipelines.domain_analysis.partition.service.create_service')
    def test_partition_single_protein(self, mock_create_service):
        """Test partition_single_protein convenience function"""
        mock_service = Mock()
        mock_create_service.return_value = mock_service

        mock_result = DomainPartitionResult(
            pdb_id="1abc", chain_id="A", reference="develop291",
            success=True
        )
        mock_service.partition_protein.return_value = mock_result

        result = partition_single_protein(
            "1abc", "A", "/path/to/summary.xml", "/path/to/output"
        )

        assert result.success == True
        mock_service.partition_protein.assert_called_once()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
