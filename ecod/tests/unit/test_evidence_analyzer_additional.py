#!/usr/bin/env python3
"""
Additional tests to improve coverage for domain analysis partition components.
"""

import pytest
import tempfile
import os
import json
import xml.etree.ElementTree as ET
from unittest.mock import Mock, MagicMock, patch, mock_open
from datetime import datetime, timedelta
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import threading
import time

# Import classes being tested
from ecod.pipelines.domain_analysis.partition.analyzer import EvidenceAnalyzer
from ecod.pipelines.domain_analysis.partition.processor import PartitionProcessor
from ecod.pipelines.domain_analysis.partition.service import DomainPartitionService
from ecod.pipelines.domain_analysis.partition.tracker import StatusTracker
from ecod.pipelines.domain_analysis.partition.models import *
from ecod.models.pipeline.evidence import Evidence
from ecod.models.pipeline.domain import DomainModel


class TestEvidenceAnalyzerAdditional:
    """Additional coverage for EvidenceAnalyzer"""

    @pytest.fixture
    def analyzer(self):
        options = PartitionOptions(validation_level=ValidationLevel.STRICT)
        return EvidenceAnalyzer(options)

    def test_xml_repair_functionality(self, analyzer):
        """Test XML repair for malformed files"""
        malformed_xml = """<?xml version="1.0"?>
        <blast_summ_doc>
            <blast_summ pdb="1abc" chain="A"/>
            <hit domain_id="d1abcA1" evalue="1e-5">
                <query_reg>10-50 & invalid</query_reg>
            </hit>
        </blast_summ_doc>"""
        
        with patch('builtins.open', mock_open(read_data=malformed_xml)):
            result = analyzer._parse_with_repair('/fake/path.xml', 'Parse error')
            
        # Should attempt repair but handle safely
        assert 'error' in result or 'blast_summ' in result

    def test_parallel_processing_enabled(self, analyzer):
        """Test analyzer with parallel processing enabled"""
        options = PartitionOptions(parallel_processing=True, max_workers=2)
        analyzer_parallel = EvidenceAnalyzer(options)
        
        # Should have executor
        assert analyzer_parallel.executor is not None
        
        # Cleanup
        analyzer_parallel.cleanup_resources()

    def test_encoding_fallback_parsing(self, analyzer):
        """Test parsing with different encodings"""
        # Simulate file with encoding issues
        with patch('builtins.open') as mock_file:
            mock_file.side_effect = [
                UnicodeDecodeError('utf-8', b'', 0, 1, 'invalid start byte'),
                mock_open(read_data='<blast_summ_doc></blast_summ_doc>')()
            ]
            
            result = analyzer._parse_with_encoding_fallback('/fake/path.xml')
            
        # Should handle encoding fallback
        assert isinstance(result, dict)

    def test_classification_cache_ttl_expiry(self, analyzer):
        """Test classification cache TTL expiry"""
        if analyzer.classification_cache:
            cache = analyzer.classification_cache
            
            # Add entry
            cache.set_domain_classification('test_domain', {'t_group': '1.1.1.1'})
            
            # Should be in cache
            result = cache.get_domain_classification('test_domain')
            assert result is not None
            
            # Simulate expiry by setting old timestamp
            cache._timestamps['test_domain'] = datetime.now() - timedelta(hours=25)
            
            # Should be expired
            result = cache.get_domain_classification('test_domain')
            assert result is None

    def test_evidence_quality_edge_cases(self, analyzer):
        """Test quality metrics with edge cases"""
        # Empty evidence list
        metrics = analyzer.calculate_quality_metrics([], 100, 1.0)
        assert metrics.total_evidence_count == 0
        assert metrics.sequence_coverage == 0.0
        
        # Single evidence item
        evidence = [Evidence(type="blast", confidence=0.8, query_range="10-50")]
        metrics = analyzer.calculate_quality_metrics(evidence, 100, 1.0)
        assert metrics.total_evidence_count == 1
        assert metrics.sequence_coverage == 41.0  # 41 residues out of 100


class TestPartitionProcessorAdditional:
    """Additional coverage for PartitionProcessor"""

    @pytest.fixture
    def processor_with_db(self):
        options = PartitionOptions()
        analyzer = Mock()
        db_manager = Mock()
        return PartitionProcessor(options, analyzer, db_manager)

    def test_discontinuous_domain_edge_cases(self):
        """Test discontinuous domain handling edge cases"""
        # Empty segments
        candidate = DiscontinuousDomainCandidate(
            segments=[],
            evidence_group=EvidenceGroup(),
            confidence=0.8
        )
        assert candidate.start == 0
        assert candidate.end == 0
        assert candidate.size == 0

        # Single segment (should behave like continuous)
        candidate = DiscontinuousDomainCandidate(
            segments=[(10, 50)],
            evidence_group=EvidenceGroup(),
            confidence=0.8
        )
        assert candidate.start == 10
        assert candidate.end == 50
        assert candidate.size == 41

    def test_complex_overlap_resolution(self, processor_with_db):
        """Test complex overlap scenarios"""
        # Multiple overlapping domains with mixed protection
        candidates = [
            DomainCandidate(start=10, end=50, evidence_group=EvidenceGroup(), 
                          confidence=0.9, protected=True),
            DomainCandidate(start=30, end=70, evidence_group=EvidenceGroup(), 
                          confidence=0.8, protected=False),
            DomainCandidate(start=60, end=100, evidence_group=EvidenceGroup(), 
                          confidence=0.7, protected=True)
        ]
        
        resolved = processor_with_db.resolve_domain_overlaps(candidates, 150)
        
        # All protected domains should be kept
        protected_count = sum(1 for d in resolved if d.protected)
        assert protected_count >= 2

    def test_classification_database_fallback(self, processor_with_db):
        """Test classification assignment with database fallback"""
        # Mock database response
        processor_with_db.db.execute_dict_query.return_value = [{
            't_group': '2004.1.1.1',
            'h_group': '2004.1.1',
            'x_group': '2004.1',
            'a_group': 'a.39',
            'is_manual_rep': True
        }]
        
        # Domain with incomplete classification
        evidence = Evidence(type="blast", domain_id="d1abcA1", confidence=0.8)
        domain = DomainModel(id="test", start=10, end=50, range="10-50", evidence=[evidence])
        
        processor_with_db.assign_domain_classifications([domain])
        
        # Should get classification from database
        assert domain.t_group == '2004.1.1.1'

    def test_processor_statistics_reset(self, processor_with_db):
        """Test processor statistics tracking and reset"""
        # Modify some stats
        processor_with_db.stats['domains_identified'] = 10
        processor_with_db.stats['overlaps_resolved'] = 5
        
        stats = processor_with_db.get_statistics()
        assert stats['domains_identified'] == 10
        
        # Reset stats
        processor_with_db.reset_statistics()
        stats = processor_with_db.get_statistics()
        assert stats['domains_identified'] == 0


class TestStatusTrackerAdvanced:
    """Advanced StatusTracker testing"""

    @pytest.fixture
    def tracker_with_config(self):
        config = {
            'max_retry_attempts': 2,
            'retry_delay_seconds': 1,
            'cleanup_completed_hours': 1
        }
        mock_db = Mock()
        return StatusTracker(mock_db, config)

    def test_database_reconnection_cycle(self, tracker_with_config):
        """Test database reconnection logic"""
        tracker = tracker_with_config
        
        # Initially available
        assert tracker.db_available
        
        # Simulate database failure
        tracker.db_available = False
        tracker.last_db_check = datetime.now() - timedelta(minutes=10)
        
        # Mock successful reconnection
        with patch.object(tracker, '_test_database_connection', return_value=True):
            tracker._maybe_reconnect_database()
            assert tracker.db_available

    def test_retry_mechanism(self, tracker_with_config):
        """Test error retry mechanism"""
        tracker = tracker_with_config
        tracker.start_process(123, "1abc", "A")
        
        # First error with retry
        success = tracker.add_process_error(123, "Temporary failure", retry=True)
        assert success
        
        process_info = tracker.processes[123]
        assert process_info.retry_count == 1
        assert process_info.status == ProcessStatus.RETRYING
        
        # Second error (should still retry)
        tracker.add_process_error(123, "Another failure", retry=True)
        assert process_info.retry_count == 2
        
        # Third error (should fail permanently)
        tracker.add_process_error(123, "Final failure", retry=True)
        assert process_info.status == ProcessStatus.FAILED

    def test_performance_metrics_calculation(self, tracker_with_config):
        """Test performance metrics calculation"""
        tracker = tracker_with_config
        
        # Add some completed processes
        for i in range(5):
            tracker.start_process(i, f"test{i}", "A")
            tracker.processes[i].status = ProcessStatus.COMPLETED
            tracker.processes[i].end_time = datetime.now()
        
        # Add some failed processes
        for i in range(5, 8):
            tracker.start_process(i, f"test{i}", "A")
            tracker.processes[i].status = ProcessStatus.FAILED
            tracker.processes[i].end_time = datetime.now()
        
        # Update metrics
        tracker._update_performance_metrics()
        
        metrics = tracker.get_performance_metrics()
        assert 'processes_per_hour' in metrics
        assert 'error_rate_percent' in metrics

    def test_batch_completion_status_update(self, tracker_with_config):
        """Test batch completion status tracking"""
        tracker = tracker_with_config
        batch_id = 456
        
        # Create batch processes
        for i in range(3):
            tracker.start_process(i, f"test{i}", "A", batch_id=batch_id)
        
        # Complete some processes
        tracker.processes[0].status = ProcessStatus.COMPLETED
        tracker.processes[1].status = ProcessStatus.FAILED
        tracker.processes[2].status = ProcessStatus.RUNNING
        
        success = tracker.update_batch_completion_status(batch_id, representatives_only=False)
        assert success

    def test_non_representative_status_update(self, tracker_with_config):
        """Test non-representative protein status updates"""
        tracker = tracker_with_config
        batch_id = 789
        
        # Create processes with different representative status
        tracker.start_process(1, "rep1", "A", batch_id=batch_id)
        tracker.processes[1].metadata['is_representative'] = True
        
        tracker.start_process(2, "nonrep1", "A", batch_id=batch_id)
        tracker.processes[2].metadata['is_representative'] = False
        
        success = tracker.update_non_representative_status(batch_id)
        assert success
        
        # Non-representative should be marked as skipped
        assert tracker.processes[2].status == ProcessStatus.CANCELLED


class TestServiceIntegration:
    """Integration testing for DomainPartitionService"""

    @pytest.fixture
    def integration_service(self):
        mock_context = Mock()
        mock_context.config_manager.config = {
            'partition': {'validation_level': 'normal'},
            'reference': {'current_version': 'test'}
        }
        mock_context.config_manager.get_db_config.return_value = {}
        
        with patch('ecod.pipelines.domain_analysis.partition.service.DBManager'):
            service = DomainPartitionService(mock_context)
            
        return service

    def test_batch_processing_with_errors(self, integration_service):
        """Test batch processing with mixed success/failure"""
        service = integration_service
        
        # Mock database responses
        service.db.execute_dict_query.return_value = [
            {'protein_id': 1, 'pdb_id': '1abc', 'chain_id': 'A', 'process_id': 101, 'is_representative': True},
            {'protein_id': 2, 'pdb_id': 'fail', 'chain_id': 'B', 'process_id': 102, 'is_representative': True}
        ]
        
        # Mock file finding - one success, one failure
        def mock_find_summary(batch_path, pdb_id, chain_id, blast_only=False):
            if pdb_id == '1abc':
                return '/fake/path/1abc_A.xml'
            else:
                return None  # File not found
        
        with patch.object(service, '_find_domain_summary', side_effect=mock_find_summary):
            with patch.object(service, 'partition_protein') as mock_partition:
                # Mock successful partition for 1abc
                mock_partition.return_value = Mock(success=True, domains=[])
                
                results = service.partition_batch(123, '/fake/batch')
                
        assert results.total >= 1
        assert results.failure_count >= 1  # At least one should fail due to missing file

    def test_service_resource_cleanup(self, integration_service):
        """Test service resource cleanup"""
        service = integration_service
        
        # Service should clean up resources
        with service as s:
            assert s is service
        
        # Should exit cleanly

    def test_configuration_validation_failures(self):
        """Test service setup with invalid configurations"""
        mock_context = Mock()
        mock_context.config_manager.config = {}  # Empty config
        mock_context.config_manager.get_db_config.side_effect = Exception("DB config error")
        
        with patch('ecod.pipelines.domain_analysis.partition.service.DBManager'):
            service = DomainPartitionService(mock_context)
            validation = service.validate_setup()
            
        # Should detect configuration issues
        assert not all(validation.values())


class TestErrorHandlingAndEdgeCases:
    """Test error handling and edge cases"""

    def test_evidence_with_malformed_ranges(self):
        """Test evidence with malformed range strings"""
        options = PartitionOptions()
        analyzer = EvidenceAnalyzer(options)
        
        # Test various malformed ranges
        malformed_ranges = [
            "",
            "invalid",
            "10-",
            "-50",
            "abc-def",
            "50-10",  # Reverse range
            "10-50-90",  # Too many parts
        ]
        
        for range_str in malformed_ranges:
            ranges = analyzer._parse_range_comprehensive(range_str)
            # Should handle gracefully (empty list or valid ranges only)
            assert isinstance(ranges, list)

    def test_domain_validation_edge_cases(self):
        """Test domain validation with edge cases"""
        options = PartitionOptions(min_domain_size=20)
        analyzer = EvidenceAnalyzer(options)
        
        # Domain with start > end
        domain = DomainModel(id="bad1", start=50, end=10, range="50-10")
        result = analyzer.validate_domain(domain, "test")
        assert not result.is_valid
        
        # Domain too small
        domain = DomainModel(id="small", start=10, end=15, range="10-15")
        result = analyzer.validate_domain(domain, "test", sequence_length=100)
        assert len(result.warnings) > 0 or not result.is_valid

    def test_partition_context_timing(self):
        """Test partition context timing functionality"""
        context = PartitionContext(
            pdb_id="1abc",
            chain_id="A", 
            reference="test",
            sequence_length=100
        )
        
        # Record stage times
        context.record_stage_time(PartitionStage.LOADING_SUMMARY)
        time.sleep(0.01)  # Small delay
        context.record_stage_time(PartitionStage.COMPLETE)
        
        total_time = context.get_total_time()
        assert total_time > 0
        assert len(context.stage_times) == 2

    def test_concurrent_status_tracking(self):
        """Test status tracker under concurrent access"""
        mock_db = Mock()
        tracker = StatusTracker(mock_db)
        
        def worker(worker_id):
            for i in range(5):
                process_id = worker_id * 100 + i
                tracker.start_process(process_id, f"test{process_id}", "A")
                tracker.update_process_stage(
                    process_id, 
                    ProcessStage.COMPLETED, 
                    ProcessStatus.COMPLETED
                )
        
        # Run multiple workers concurrently
        threads = []
        for i in range(3):
            thread = threading.Thread(target=worker, args=(i,))
            threads.append(thread)
            thread.start()
        
        for thread in threads:
            thread.join()
        
        # Should have processed all without errors
        assert len(tracker.processes) == 15  # 3 workers * 5 processes each


class TestMemoryAndPerformance:
    """Test memory usage and performance characteristics"""

    def test_large_batch_memory_usage(self):
        """Test memory usage with large batches"""
        results = BatchPartitionResults()
        
        # Add many results
        for i in range(1000):
            result = Mock()
            result.success = True
            result.domains = []
            result.is_peptide = False
            result.is_classified = False
            results.add_result(result)
        
        summary = results.get_summary()
        assert summary['total_proteins'] == 1000

    def test_cache_memory_management(self):
        """Test cache memory management"""
        cache = ClassificationCache(max_size=100, ttl_hours=1)
        
        # Fill beyond capacity
        for i in range(150):
            cache.set_domain_classification(f"domain_{i}", {'t_group': f'group_{i}'})
        
        # Should not exceed max size due to eviction
        stats = cache.get_stats()
        assert stats['size'] <= 100
        assert stats['evictions'] > 0


# Run additional tests
if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
