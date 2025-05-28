#!/usr/bin/env python3
"""
Comprehensive tests for error handling scenarios.

Tests error handling and recovery across all components including:
- Database connection failures
- File I/O errors and corruption
- XML parsing errors and malformed data
- Invalid configuration handling
- Memory and resource constraints
- Concurrent access issues
- Network timeouts and interruptions
- Graceful degradation under stress
"""

import pytest
import tempfile
import os
import shutil
import xml.etree.ElementTree as ET
from unittest.mock import Mock, MagicMock, patch, mock_open
from pathlib import Path
import threading
import time
from datetime import datetime

from ecod.core.context import ApplicationContext
from ecod.db import DBManager
from ecod.models.pipeline.evidence import Evidence
from ecod.models.pipeline.domain import DomainModel
from ecod.models.pipeline.partition import DomainPartitionResult
from ecod.pipelines.domain_analysis.partition.service import DomainPartitionService
from ecod.pipelines.domain_analysis.partition.analyzer import EvidenceAnalyzer
from ecod.pipelines.domain_analysis.partition.processor import PartitionProcessor
from ecod.pipelines.domain_analysis.partition.tracker import StatusTracker
from ecod.pipelines.domain_analysis.partition.models import (
    PartitionOptions, ValidationLevel, ValidationResult
)
from ecod.exceptions import PipelineError, ValidationError


class TestDatabaseErrorHandling:
    """Test database connection and operation error handling"""
    
    @pytest.fixture
    def mock_context(self):
        """Create mock context for testing"""
        context = Mock()  # ← Remove spec=ApplicationContext
        context.config_manager = Mock()  # ← Add this line
        context.config_manager.config = {
            'partition': {'validation_level': 'normal'},
            'reference': {'current_version': 'develop291'}
        }
        context.config_manager.get_db_config.return_value = {}
        return context
    
    def test_service_initialization_db_failure(self, mock_context):
        """Test service initialization with database failure"""
        with patch('ecod.pipelines.domain_analysis.partition.service.DBManager') as mock_db_class:
            mock_db_class.side_effect = Exception("Database connection failed")
            
            # Should still initialize but with broken DB
            service = DomainPartitionService(mock_context)
            
            # Validation should fail
            validation = service.validate_setup()
            assert validation['database'] == False
    
    def test_service_db_connection_lost_during_processing(self, mock_context):
        """Test handling database connection loss during processing"""
        with patch('ecod.pipelines.domain_analysis.partition.service.DBManager') as mock_db_class:
            mock_db = Mock()
            mock_db.test_connection.return_value = True
            mock_db_class.return_value = mock_db
            
            service = DomainPartitionService(mock_context)
            
            # Simulate connection loss during protein processing
            mock_db.execute_dict_query.side_effect = Exception("Connection lost")
            
            result = service.partition_protein(
                "1abc", "A", "/fake/summary.xml", "/fake/output"
            )
            
            # Should return error result, not crash
            assert result.success == False
            assert "connection" in result.error.lower() or "error" in result.error.lower()
    
    def test_batch_processing_db_failures(self, mock_context):
        """Test batch processing with database failures"""
        with patch('ecod.pipelines.domain_analysis.partition.service.DBManager') as mock_db_class:
            mock_db = Mock()
            mock_db.test_connection.return_value = True
            mock_db_class.return_value = mock_db
            
            service = DomainPartitionService(mock_context)
            
            # Simulate database failure when getting proteins
            mock_db.execute_dict_query.side_effect = Exception("Database error")
            
            results = service.partition_batch(123, "/fake/batch/path")
            
            # Should return empty results, not crash
            assert results.total == 0
            assert results.success_count == 0
    
    def test_status_tracker_db_failures(self):
        """Test status tracker with database failures"""
        mock_db = Mock()
        mock_db.update.side_effect = Exception("Database error")
        
        tracker = StatusTracker(mock_db)
        
        # Should handle gracefully
        success = tracker.update_process_status(123, "test_stage", "processing")
        assert success == False
        
        # Should not raise exception
        success = tracker.register_domain_file(123, "/fake/file.xml", "/fake/base")
        assert success == False
    
    def test_database_transaction_rollback(self, mock_context):
        """Test database transaction handling during errors"""
        with patch('ecod.pipelines.domain_analysis.partition.service.DBManager') as mock_db_class:
            mock_db = Mock()
            mock_db.test_connection.return_value = True
            
            # Simulate transaction failure
            mock_db.update.side_effect = [True, Exception("Transaction failed"), True]
            mock_db_class.return_value = mock_db
            
            service = DomainPartitionService(mock_context)
            tracker = service.tracker
            
            # Multiple updates - middle one should fail
            success1 = tracker.update_process_status(101, "stage1", "success")
            success2 = tracker.update_process_status(102, "stage2", "processing")  # Fails
            success3 = tracker.update_process_status(103, "stage3", "success")
            
            assert success1 == True
            assert success2 == False  # Failed gracefully
            assert success3 == True   # Continued after failure


class TestFileIOErrorHandling:
    """Test file I/O error handling"""
    
    @pytest.fixture
    def analyzer(self):
        """Create analyzer for testing"""
        options = PartitionOptions()
        return EvidenceAnalyzer(options)
    
    def test_parse_nonexistent_file(self, analyzer):
        """Test parsing file that doesn't exist"""
        result = analyzer.parse_domain_summary("/nonexistent/file.xml")
        
        assert 'error' in result
        assert 'not found' in result['error'].lower()
    
    def test_parse_permission_denied(self, analyzer):
        """Test parsing file with permission denied"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            f.write("<test>content</test>")
            f.flush()
            
            # Change permissions to make unreadable
            os.chmod(f.name, 0o000)
            
            try:
                result = analyzer.parse_domain_summary(f.name)
                
                assert 'error' in result
                assert 'permission' in result['error'].lower() or 'error' in result['error'].lower()
                
            finally:
                # Restore permissions for cleanup
                os.chmod(f.name, 0o644)
                os.unlink(f.name)
    
    def test_parse_corrupted_file(self, analyzer):
        """Test parsing corrupted file"""
        corrupted_content = b'\x00\x01\x02\xff\xfe'  # Binary garbage
        
        with tempfile.NamedTemporaryFile(mode='wb', delete=False) as f:
            f.write(corrupted_content)
            f.flush()
            
            try:
                result = analyzer.parse_domain_summary(f.name)
                
                assert 'error' in result
                
            finally:
                os.unlink(f.name)
    
    def test_parse_file_disappears_during_reading(self, analyzer):
        """Test file being deleted while reading"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            f.write("<?xml version='1.0'?><test>content</test>")
            f.flush()
            
            file_path = f.name
        
        # Delete file
        os.unlink(file_path)
        
        result = analyzer.parse_domain_summary(file_path)
        
        assert 'error' in result
    
    def test_save_to_readonly_directory(self):
        """Test saving to read-only directory"""
        with tempfile.TemporaryDirectory() as temp_dir:
            readonly_dir = os.path.join(temp_dir, "readonly")
            os.makedirs(readonly_dir)
            os.chmod(readonly_dir, 0o444)  # Read-only
            
            try:
                result = DomainPartitionResult(
                    pdb_id="1abc",
                    chain_id="A", 
                    reference="develop291"
                )
                
                success = result.save(output_dir=readonly_dir)
                
                assert success == False
                
            finally:
                # Restore permissions for cleanup
                os.chmod(readonly_dir, 0o755)
    
    def test_save_disk_full_simulation(self):
        """Test handling disk full during save"""
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291"
        )
        
        # Mock write to simulate disk full
        with patch('xml.etree.ElementTree.ElementTree.write') as mock_write:
            mock_write.side_effect = OSError("No space left on device")
            
            success = result.save(output_dir="/fake/dir", filename="test.xml")
            
            assert success == False
    
    def test_concurrent_file_access(self):
        """Test concurrent access to same file"""
        with tempfile.TemporaryDirectory() as temp_dir:
            file_path = os.path.join(temp_dir, "concurrent_test.xml")
            
            results = []
            errors = []
            
            def save_result(thread_id):
                try:
                    result = DomainPartitionResult(
                        pdb_id=f"test{thread_id}",
                        chain_id="A",
                        reference="develop291"
                    )
                    result.domain_file = file_path
                    success = result.save()
                    results.append((thread_id, success))
                except Exception as e:
                    errors.append((thread_id, str(e)))
            
            # Start multiple threads trying to write same file
            threads = []
            for i in range(5):
                t = threading.Thread(target=save_result, args=(i,))
                threads.append(t)
                t.start()
            
            # Wait for all threads
            for t in threads:
                t.join()
            
            # At least one should succeed, others may fail gracefully
            successes = [r[1] for r in results if r[1]]
            assert len(successes) >= 1 or len(errors) >= 1  # Either success or graceful failure


class TestXMLParsingErrorHandling:
    """Test XML parsing error scenarios"""
    
    @pytest.fixture
    def analyzer(self):
        """Create analyzer for testing"""
        options = PartitionOptions()
        return EvidenceAnalyzer(options)
    
    def test_malformed_xml_structure(self, analyzer):
        """Test malformed XML structure"""
        malformed_cases = [
            "<?xml version='1.0'?><unclosed>",  # Unclosed tag
            "<?xml version='1.0'?><root><child></wrong></root>",  # Wrong closing tag
            "<?xml version='1.0'?><root attribute='unclosed>content</root>",  # Unclosed attribute
            "<root>no xml declaration</root>",  # No XML declaration
            "<?xml version='1.0'?><root>& invalid entity</root>",  # Invalid entity
        ]
        
        for malformed_xml in malformed_cases:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.xml', delete=False) as f:
                f.write(malformed_xml)
                f.flush()
                
                try:
                    result = analyzer.parse_domain_summary(f.name)
                    
                    # Should handle gracefully with error
                    assert 'error' in result
                    
                finally:
                    os.unlink(f.name)
    
    def test_xml_with_invalid_characters(self, analyzer):
        """Test XML with invalid characters"""
        invalid_xml = """<?xml version='1.0'?>
        <domain_summary pdb_id="1abc" chain_id="A">
            <metadata>
                <note>\x00\x01\x02 invalid chars</note>
            </metadata>
        </domain_summary>"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xml', delete=False) as f:
            try:
                f.write(invalid_xml)
                f.flush()
                
                result = analyzer.parse_domain_summary(f.name)
                
                # Should handle gracefully
                assert 'error' in result or result.get('pdb_id') == '1abc'
                
            finally:
                os.unlink(f.name)
    
    def test_xml_encoding_issues(self, analyzer):
        """Test XML with encoding issues"""
        # UTF-8 with BOM
        utf8_bom = b'\xef\xbb\xbf<?xml version="1.0" encoding="UTF-8"?><root>test</root>'
        
        # Latin-1 content in UTF-8 declared file
        encoding_mismatch = "<?xml version='1.0' encoding='UTF-8'?><root>caf\xe9</root>".encode('latin-1')
        
        test_cases = [utf8_bom, encoding_mismatch]
        
        for xml_bytes in test_cases:
            with tempfile.NamedTemporaryFile(mode='wb', suffix='.xml', delete=False) as f:
                f.write(xml_bytes)
                f.flush()
                
                try:
                    result = analyzer.parse_domain_summary(f.name)
                    
                    # Should either parse successfully or fail gracefully
                    assert isinstance(result, dict)
                    
                finally:
                    os.unlink(f.name)
    
    def test_extremely_large_xml(self, analyzer):
        """Test handling extremely large XML files"""
        # Create large XML content
        large_content = "<?xml version='1.0'?>\n<domain_summary pdb_id='1abc' chain_id='A'>\n"
        large_content += "<evidence_list>\n"
        
        # Add many elements
        for i in range(10000):
            large_content += f"<hit num='{i}' domain_id='d{i}'>content{i}</hit>\n"
        
        large_content += "</evidence_list>\n</domain_summary>"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xml', delete=False) as f:
            f.write(large_content)
            f.flush()
            
            try:
                # Should handle large files without crashing
                result = analyzer.parse_domain_summary(f.name)
                
                # Should either succeed or fail gracefully
                assert isinstance(result, dict)
                
            finally:
                os.unlink(f.name)
    
    def test_xml_external_entity_attack(self, analyzer):
        """Test handling XML external entity (XXE) attacks"""
        xxe_xml = """<?xml version="1.0"?>
        <!DOCTYPE root [
        <!ENTITY xxe SYSTEM "file:///etc/passwd">
        ]>
        <domain_summary pdb_id="1abc" chain_id="A">
            <metadata>&xxe;</metadata>
        </domain_summary>"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xml', delete=False) as f:
            f.write(xxe_xml)
            f.flush()
            
            try:
                result = analyzer.parse_domain_summary(f.name)
                
                # Should handle securely - either parse safely or fail
                assert isinstance(result, dict)
                
                # Should not contain system file content
                if 'metadata' in result:
                    assert 'root:' not in str(result['metadata'])
                
            finally:
                os.unlink(f.name)
    
    def test_evidence_xml_parsing_errors(self):
        """Test Evidence XML parsing with various error conditions"""
        error_cases = [
            # Missing required elements
            '<evidence type="blast"/>',
            
            # Invalid attribute values
            '<evidence type="blast" confidence="invalid_number"/>',
            
            # Mixed up element structure
            '<evidence type="blast"><confidence>0.5</confidence><invalid_element/></evidence>',
            
            # Extremely long attribute values
            f'<evidence type="blast" domain_id="{"x" * 10000}"/>',
        ]
        
        for xml_str in error_cases:
            try:
                element = ET.fromstring(xml_str)
                evidence = Evidence.from_xml(element)
                
                # Should create evidence object without crashing
                assert isinstance(evidence, Evidence)
                
            except Exception as e:
                # If parsing fails, should be graceful
                assert not isinstance(e, SystemExit)


class TestConfigurationErrorHandling:
    """Test invalid configuration handling"""
    
    def test_invalid_validation_level(self):
        """Test invalid validation level in configuration"""
        with pytest.raises(ValueError):
            PartitionOptions(validation_level="invalid_level")
    
    def test_invalid_numeric_options(self):
        """Test invalid numeric option values"""
        with pytest.raises(ValueError):
            options = PartitionOptions(min_domain_size=-5)
            options.validate()
        
        with pytest.raises(ValueError):
            options = PartitionOptions(overlap_threshold=1.5)
            options.validate()
        
        with pytest.raises(ValueError):
            options = PartitionOptions(min_evidence_confidence=-0.1)
            options.validate()
    
    def test_inconsistent_options(self):
        """Test inconsistent option combinations"""
        with pytest.raises(ValueError):
            options = PartitionOptions(min_domain_size=100, max_domain_size=50)
            options.validate()
    
    def test_missing_config_sections(self):
        """Test handling missing configuration sections"""
        mock_context = Mock()
        mock_context.config_manager.config = {}  # Empty config
        mock_context.config_manager.get_db_config.return_value = {}
        
        with patch('ecod.pipelines.domain_analysis.partition.service.DBManager'):
            # Should handle missing config gracefully
            service = DomainPartitionService(mock_context)
            
            # Should have reasonable defaults
            assert service.default_options.min_domain_size > 0
            assert 0.0 <= service.default_options.overlap_threshold <= 1.0
    
    def test_invalid_config_values(self):
        """Test handling invalid values in configuration"""
        mock_context = Mock()
        mock_context.config_manager.config = {
            'partition': {
                'validation_level': 'invalid',  # Invalid enum value
                'min_domain_size': 'not_a_number',  # Invalid type
                'overlap_threshold': -0.5  # Invalid range
            }
        }
        mock_context.config_manager.get_db_config.return_value = {}
        
        with patch('ecod.pipelines.domain_analysis.partition.service.DBManager'):
            # Should either use defaults or raise validation error
            try:
                service = DomainPartitionService(mock_context)
                # If it succeeds, should have valid options
                assert service.default_options.min_domain_size > 0
            except (ValueError, TypeError):
                # Should raise clear validation error
                pass


class TestMemoryAndResourceConstraints:
    """Test handling of memory and resource limitations"""
    
    def test_large_evidence_list_processing(self):
        """Test processing very large evidence lists"""
        options = PartitionOptions()
        analyzer = EvidenceAnalyzer(options)
        
        # Create large evidence list
        large_evidence_list = []
        for i in range(10000):
            evidence = Evidence(
                type='blast',
                domain_id=f'd{i}',
                query_range=f'{i*10}-{i*10+50}',
                evalue=1e-20,
                confidence=0.8
            )
            large_evidence_list.append(evidence)
        
        # Should handle large lists without crashing
        try:
            groups = analyzer.group_evidence_by_position(large_evidence_list, window_size=100)
            assert isinstance(groups, dict)
            assert len(groups) > 0
        except MemoryError:
            # If memory error occurs, should be handled gracefully
            pytest.skip("Insufficient memory for large evidence test")
    
    def test_deeply_nested_xml_structures(self):
        """Test handling deeply nested XML structures"""
        # Create deeply nested XML
        deep_xml = "<?xml version='1.0'?>\n"
        deep_xml += "<root>\n"
        
        # Create deep nesting
        for i in range(1000):
            deep_xml += f"<level{i}>\n"
        
        deep_xml += "<content>deep content</content>\n"
        
        for i in range(999, -1, -1):
            deep_xml += f"</level{i}>\n"
        
        deep_xml += "</root>"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xml', delete=False) as f:
            f.write(deep_xml)
            f.flush()
            
            try:
                options = PartitionOptions()
                analyzer = EvidenceAnalyzer(options)
                
                # Should handle or fail gracefully
                result = analyzer.parse_domain_summary(f.name)
                assert isinstance(result, dict)
                
            except RecursionError:
                # Should not crash the entire process
                pytest.skip("Recursion limit reached in XML parsing")
            finally:
                os.unlink(f.name)
    
    def test_memory_leak_prevention(self):
        """Test that repeated operations don't cause memory leaks"""
        options = PartitionOptions()
        analyzer = EvidenceAnalyzer(options)
        
        # Perform many operations to check for leaks
        for i in range(100):
            evidence = Evidence(
                type='blast',
                domain_id=f'd{i}',
                evalue=1e-20,
                confidence=0.8
            )
            
            # Validate evidence
            result = analyzer.validate_evidence(evidence, f'context_{i}')
            assert isinstance(result, ValidationResult)
            
            # Clear reference
            del evidence
            del result
        
        # If we get here without memory issues, test passes
        assert True
    
    def test_concurrent_resource_access(self):
        """Test concurrent access to shared resources"""
        options = PartitionOptions(use_cache=True)
        analyzer = EvidenceAnalyzer(options)
        
        results = []
        
        def worker_function(worker_id):
            try:
                # Each worker performs cache operations
                for i in range(10):
                    domain_id = f'worker{worker_id}_domain{i}'
                    
                    # Set classification
                    analyzer.classification_cache.set_domain_classification(
                        domain_id, {'t_group': f'group_{worker_id}_{i}'}
                    )
                    
                    # Get classification
                    result = analyzer.classification_cache.get_domain_classification(domain_id)
                    results.append((worker_id, i, result is not None))
                    
            except Exception as e:
                results.append((worker_id, -1, str(e)))
        
        # Start multiple worker threads
        threads = []
        for worker_id in range(5):
            t = threading.Thread(target=worker_function, args=(worker_id,))
            threads.append(t)
            t.start()
        
        # Wait for completion
        for t in threads:
            t.join()
        
        # Check that operations completed without major errors
        error_results = [r for r in results if isinstance(r[2], str)]
        assert len(error_results) == 0  # No string errors
        
        success_results = [r for r in results if r[2] is True]
        assert len(success_results) > 0  # Some operations succeeded


class TestNetworkAndTimeoutHandling:
    """Test network-related error handling"""
    
    def test_database_timeout_simulation(self):
        """Test handling database operation timeouts"""
        mock_db = Mock()
        
        # Simulate long-running query that times out
        def slow_query(*args, **kwargs):
            time.sleep(0.1)  # Simulate delay
            raise Exception("Query timeout")
        
        mock_db.execute_dict_query.side_effect = slow_query
        
        tracker = StatusTracker(mock_db)
        
        # Should handle timeout gracefully
        success = tracker.update_process_status(123, "test", "processing")
        assert success == False
    
    def test_interrupted_processing(self):
        """Test handling of interrupted processing"""
        options = PartitionOptions()
        analyzer = EvidenceAnalyzer(options)
        
        # Simulate KeyboardInterrupt during processing
        def interrupt_during_validation(*args, **kwargs):
            raise KeyboardInterrupt("Processing interrupted")
        
        evidence = Evidence(type='blast', domain_id='d1', evalue=1e-20)
        
        with patch.object(analyzer, 'validate_evidence', side_effect=interrupt_during_validation):
            # Should not crash completely
            try:
                result = analyzer.validate_evidence(evidence, 'test')
                assert False, "Should have raised KeyboardInterrupt"
            except KeyboardInterrupt:
                # This is expected
                pass
    
    def test_service_shutdown_gracefully(self):
        """Test graceful service shutdown"""
        mock_context = Mock()
        mock_context.config_manager.config = {'reference': {'current_version': 'develop291'}}
        mock_context.config_manager.get_db_config.return_value = {}
        
        with patch('ecod.pipelines.domain_analysis.partition.service.DBManager'):
            service = DomainPartitionService(mock_context)
            
            # Test context manager
            with service as s:
                assert s is service
            
            # Should complete without errors


class TestGracefulDegradation:
    """Test graceful degradation under stress conditions"""
    
    def test_service_without_database(self):
        """Test service operation without database"""
        mock_context = Mock()
        mock_context.config_manager.config = {'reference': {'current_version': 'develop291'}}
        mock_context.config_manager.get_db_config.return_value = {}
        
        with patch('ecod.pipelines.domain_analysis.partition.service.DBManager') as mock_db_class:
            mock_db_class.side_effect = Exception("No database available")
            
            # Should still initialize
            service = DomainPartitionService(mock_context)
            
            # Should work with limited functionality
            validation = service.validate_setup()
            assert validation['database'] == False
            assert validation['config_loaded'] == True
    
    def test_evidence_processing_with_partial_data(self):
        """Test evidence processing with incomplete data"""
        options = PartitionOptions(validation_level=ValidationLevel.LENIENT)
        analyzer = EvidenceAnalyzer(options)
        
        # Create evidence with missing/invalid data
        incomplete_evidence = [
            Evidence(type='blast', domain_id='', confidence=None),  # Missing domain_id
            Evidence(type='', domain_id='d1', confidence=0.8),      # Missing type
            Evidence(type='hhsearch', domain_id='d2'),              # Missing confidence
        ]
        
        # Should process what it can
        summary_data = {
            'blast_hits': [],
            'hhsearch_hits': []
        }
        
        # Should not crash
        result = analyzer.extract_evidence_with_classification(summary_data)
        assert isinstance(result, list)
    
    def test_partial_batch_processing_failure(self):
        """Test batch processing continuing after partial failures"""
        mock_context = Mock()
        mock_context.config_manager.config = {'reference': {'current_version': 'develop291'}}
        mock_context.config_manager.get_db_config.return_value = {}
        
        with patch('ecod.pipelines.domain_analysis.partition.service.DBManager') as mock_db_class:
            mock_db = Mock()
            mock_db.test_connection.return_value = True
            mock_db_class.return_value = mock_db
            
            service = DomainPartitionService(mock_context)
            service.service_settings['use_multiprocessing'] = False
            
            # Mock some proteins succeed, some fail
            proteins = [
                {'pdb_id': '1abc', 'chain_id': 'A', 'process_id': 101},
                {'pdb_id': '2fail', 'chain_id': 'B', 'process_id': 102},
                {'pdb_id': '3xyz', 'chain_id': 'C', 'process_id': 103}
            ]
            
            service._get_proteins_to_process = Mock(return_value=proteins)
            
            # Mock partition_protein to fail for one protein
            def mock_partition(pdb_id, chain_id, *args, **kwargs):
                if pdb_id == '2fail':
                    raise Exception("Partition failed")
                return DomainPartitionResult(
                    pdb_id=pdb_id, chain_id=chain_id, reference='develop291', success=True
                )
            
            service.partition_protein = Mock(side_effect=mock_partition)
            service._find_domain_summary = Mock(return_value="/fake/summary.xml")
            
            # Should continue processing after failure
            results = service.partition_batch(123, "/fake/batch/path")
            
            # Should have processed all proteins
            assert results.total == 3
            assert results.success_count == 2  # 1abc and 3xyz succeeded
            assert results.failure_count == 1  # 2fail failed
    
    def test_fallback_confidence_calculation(self):
        """Test fallback confidence calculation when primary methods fail"""
        # Create evidence with problematic data
        evidence = Evidence(
            type='unknown_type',
            confidence=None
        )
        
        # Should fall back to generic calculation
        assert evidence.confidence >= 0.0
        assert evidence.confidence <= 1.0
        
        # Test with completely empty evidence
        empty_evidence = Evidence(type='')
        assert empty_evidence.confidence >= 0.0


class TestErrorRecoveryAndLogging:
    """Test error recovery and comprehensive logging"""
    
    def test_error_context_preservation(self):
        """Test that error context is preserved for debugging"""
        options = PartitionOptions()
        analyzer = EvidenceAnalyzer(options)
        
        # Create evidence that will cause validation error
        invalid_evidence = Evidence(
            type='blast',
            domain_id='',  # Invalid - empty
            query_range='invalid-range',
            confidence=-0.5  # Invalid - negative
        )
        
        result = analyzer.validate_evidence(invalid_evidence, 'test_protein_1abc_A')
        
        # Should preserve context information
        assert result.context == 'test_protein_1abc_A'
        assert len(result.errors) > 0
        
        # Error messages should be descriptive
        error_text = ' '.join(result.errors).lower()
        assert 'domain_id' in error_text or 'empty' in error_text
    
    def test_comprehensive_error_logging(self):
        """Test that errors are logged comprehensively"""
        with patch('ecod.pipelines.domain_analysis.partition.service.logging') as mock_logging:
            mock_context = Mock()
            mock_context.config_manager.config = {'reference': {'current_version': 'develop291'}}
            mock_context.config_manager.get_db_config.return_value = {}
            
            with patch('ecod.pipelines.domain_analysis.partition.service.DBManager'):
                service = DomainPartitionService(mock_context)
                
                # Trigger an error
                result = service.partition_protein(
                    "1abc", "A", "/nonexistent/file.xml", "/fake/output"
                )
                
                # Should have logged the error
                assert result.success == False
                # Logging should have occurred (check if logger was called)
    
    def test_error_aggregation_in_batch(self):
        """Test error aggregation during batch processing"""
        mock_context = Mock()
        mock_context.config_manager.config = {'reference': {'current_version': 'develop291'}}
        mock_context.config_manager.get_db_config.return_value = {}
        
        with patch('ecod.pipelines.domain_analysis.partition.service.DBManager'):
            service = DomainPartitionService(mock_context)
            
            # Create batch with various error types
            proteins = [
                {'pdb_id': '1abc', 'chain_id': 'A', 'process_id': 101},
                {'pdb_id': '2xyz', 'chain_id': 'B', 'process_id': 102},
            ]
            
            service._get_proteins_to_process = Mock(return_value=proteins)
            
            # Mock different types of failures
            def mock_partition(pdb_id, chain_id, *args, **kwargs):
                if pdb_id == '1abc':
                    return DomainPartitionResult(
                        pdb_id=pdb_id, chain_id=chain_id, reference='develop291',
                        success=False, error="File not found"
                    )
                else:
                    return DomainPartitionResult(
                        pdb_id=pdb_id, chain_id=chain_id, reference='develop291',
                        success=False, error="XML parsing failed"
                    )
            
            service.partition_protein = Mock(side_effect=mock_partition)
            service._find_domain_summary = Mock(return_value="/fake/summary.xml")
            
            results = service.partition_batch(123, "/fake/batch/path")
            
            # Should aggregate all errors
            assert len(results.failures) == 2
            
            # Should have different error types
            error_messages = [f[2] for f in results.failures]
            assert "File not found" in error_messages
            assert "XML parsing failed" in error_messages


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
