#!/usr/bin/env python3
"""
Comprehensive tests for batch processing functionality.

Tests batch processing integration including:
- Sequential batch processing
- Parallel batch processing  
- Progress tracking and status updates
- Error handling during batch operations
- Reprocessing failed proteins
- File discovery and validation
- Database status management
"""

import pytest
import tempfile
import os
import shutil
from datetime import datetime
from pathlib import Path
from unittest.mock import Mock, MagicMock, patch, call
from concurrent.futures import Future

from ecod.core.context import ApplicationContext
from ecod.db import DBManager
from ecod.models.pipeline.partition import DomainPartitionResult
from ecod.pipelines.domain_analysis.partition.service import DomainPartitionService
from ecod.pipelines.domain_analysis.partition.models import (
    BatchPartitionResults, PartitionOptions, ProcessingMode, ValidationLevel
)


class TestBatchProcessingSetup:
    """Test batch processing setup and initialization"""
    
    @pytest.fixture
    def mock_context(self):
        context = Mock(spec=ApplicationContext)
        context.config_manager = Mock()  # ← Add this line
        context.config_manager.config = {
            'partition': {'validation_level': 'normal'},
            'reference': {'current_version': 'develop291'}
        }
        context.config_manager.get_db_config.return_value = {}
        return context
    
    @pytest.fixture
    def mock_db(self):
        """Create mock database manager"""
        db = Mock()  # ← Remove spec=DBManager
        db.test_connection.return_value = True
        db.execute_dict_query.return_value = []
        db.execute_query.return_value = []
        return db
    
    @pytest.fixture
    def service(self, mock_context, mock_db):
        """Create service with mocked dependencies"""
        with patch('ecod.pipelines.domain_analysis.partition.service.DBManager', return_value=mock_db):
            with patch('ecod.pipelines.domain_analysis.partition.service.PartitionProcessor') as mock_processor:
                # Make sure PartitionProcessor is also mocked to avoid real processing
                service = DomainPartitionService(mock_context)
                return service

    def test_batch_service_initialization(self, service):
        """Test batch service initializes correctly"""
        assert service.service_settings['max_workers'] >= 1
        assert service.service_settings['batch_size'] > 0
        assert 'proteins_processed' in service.service_stats
    
    def test_custom_service_config(self, mock_context, mock_db):
        """Test service with custom configuration"""
        service_config = {
            'max_workers': 8,
            'use_multiprocessing': True,
            'batch_size': 200,
            'save_intermediate': False
        }
        
        with patch('ecod.pipelines.domain_analysis.partition.service.DBManager', return_value=mock_db):
            service = DomainPartitionService(mock_context, service_config)
            
            assert service.service_settings['max_workers'] == 8
            assert service.service_settings['use_multiprocessing'] == True
            assert service.service_settings['batch_size'] == 200
            assert service.service_settings['save_intermediate'] == False


class TestBatchProteinRetrieval:
    """Test retrieving proteins to process from database"""
    
    @pytest.fixture
    def service(self):
        """Create service for testing"""
        mock_context = Mock()
        mock_context.config_manager.config = {'reference': {'current_version': 'develop291'}}
        mock_context.config_manager.get_db_config.return_value = {}
        
        with patch('ecod.pipelines.domain_analysis.partition.service.DBManager') as mock_db_class:
            service = DomainPartitionService(mock_context)
            return service
    
    def test_get_proteins_to_process_basic(self, service):
        """Test basic protein retrieval"""
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
        
        # Check query was called correctly
        service.db.execute_dict_query.assert_called_once()
        query_args = service.db.execute_dict_query.call_args[0]
        assert 'batch_id = %s' in query_args[0]
        assert query_args[1] == (123,)
    
    def test_get_proteins_representatives_only(self, service):
        """Test retrieving only representative proteins"""
        service.db.execute_dict_query.return_value = [
            {
                'protein_id': 1,
                'pdb_id': '1abc',
                'chain_id': 'A',
                'process_id': 101,
                'is_representative': True
            }
        ]
        
        proteins = service._get_proteins_to_process(
            batch_id=123, 
            representatives_only=True
        )
        
        assert len(proteins) == 1
        assert proteins[0]['is_representative'] == True
        
        # Check representatives filter was applied
        query_args = service.db.execute_dict_query.call_args[0]
        assert 'is_representative = TRUE' in query_args[0]
    
    def test_get_proteins_with_limit(self, service):
        """Test protein retrieval with limit"""
        service.db.execute_dict_query.return_value = []
        
        service._get_proteins_to_process(batch_id=123, limit=50)
        
        query_args = service.db.execute_dict_query.call_args[0]
        assert 'LIMIT 50' in query_args[0]
    
    def test_get_proteins_db_error(self, service):
        """Test handling database error during protein retrieval"""
        service.db.execute_dict_query.side_effect = Exception("DB Error")
        
        proteins = service._get_proteins_to_process(batch_id=123)
        
        assert proteins == []
    
    def test_get_proteins_empty_batch(self, service):
        """Test retrieving from empty batch"""
        service.db.execute_dict_query.return_value = []
        
        proteins = service._get_proteins_to_process(batch_id=123)
        
        assert proteins == []


class TestSequentialBatchProcessing:
    """Test sequential batch processing"""
    
    @pytest.fixture
    def temp_batch_dir(self):
        """Create temporary batch directory structure"""
        temp_dir = tempfile.mkdtemp()
        
        # Create domains subdirectory
        domains_dir = Path(temp_dir) / "domains"
        domains_dir.mkdir()
        
        # Create some test summary files
        summaries_dir = Path(temp_dir) / "summaries"
        summaries_dir.mkdir()
        
        summary_content = '''<?xml version="1.0"?>
        <domain_summary pdb_id="1abc" chain_id="A" sequence_length="150">
            <evidence_list>
                <hhsearch_hits>
                    <hit hit_id="h1" domain_id="d1abcA1" probability="95.0">
                        <query_range>10-100</query_range>
                    </hit>
                </hhsearch_hits>
            </evidence_list>
        </domain_summary>'''
        
        (summaries_dir / "1abc_A.develop291.domain_summary.xml").write_text(summary_content)
        (summaries_dir / "2xyz_B.develop291.domain_summary.xml").write_text(summary_content.replace("1abc", "2xyz").replace("A", "B"))
        
        yield temp_dir
        
        # Cleanup
        shutil.rmtree(temp_dir)
    
    @pytest.fixture
    def service_with_mocks(self):
        """Create service with comprehensive mocks"""
        mock_context = Mock()
        mock_context.config_manager.config = {
            'partition': {'validation_level': 'normal'},
            'reference': {'current_version': 'develop291'}
        }
        mock_context.config_manager.get_db_config.return_value = {}
        
        with patch('ecod.pipelines.domain_analysis.partition.service.DBManager') as mock_db_class:
            mock_db = Mock()
            mock_db.test_connection.return_value = True
            mock_db_class.return_value = mock_db
            
            service = DomainPartitionService(mock_context)
            service.service_settings['use_multiprocessing'] = False  # Force sequential
            
            return service, mock_db
    
    def test_sequential_batch_processing_success(self, service_with_mocks, temp_batch_dir):
        """Test successful sequential batch processing"""
        service, mock_db = service_with_mocks
        
        # Mock protein list
        proteins = [
            {'pdb_id': '1abc', 'chain_id': 'A', 'process_id': 101, 'is_representative': True},
            {'pdb_id': '2xyz', 'chain_id': 'B', 'process_id': 102, 'is_representative': True}
        ]
        
        service._get_proteins_to_process = Mock(return_value=proteins)
        
        # Mock file finding
        def mock_find_summary(batch_path, pdb_id, chain_id, blast_only=False):
            return f"{batch_path}/summaries/{pdb_id}_{chain_id}.develop291.domain_summary.xml"
        
        service._find_domain_summary = Mock(side_effect=mock_find_summary)
        
        # Mock individual processing
        def mock_partition_protein(pdb_id, chain_id, summary_path, output_dir, **kwargs):
            return DomainPartitionResult(
                pdb_id=pdb_id,
                chain_id=chain_id,
                reference='develop291',
                success=True,
                is_classified=True,
                domains=[{'id': f'd{pdb_id}{chain_id}1', 'start': 10, 'end': 100}]
            )
        
        service.partition_protein = Mock(side_effect=mock_partition_protein)
        
        # Process batch
        results = service.partition_batch(123, temp_batch_dir)
        
        # Check results
        assert results.total == 2
        assert results.success_count == 2
        assert results.failure_count == 0
        assert len(results.results) == 2
        
        # Check all proteins were processed
        assert service.partition_protein.call_count == 2
    
    def test_sequential_batch_processing_with_failures(self, service_with_mocks, temp_batch_dir):
        """Test sequential processing with some failures"""
        service, mock_db = service_with_mocks
        
        proteins = [
            {'pdb_id': '1abc', 'chain_id': 'A', 'process_id': 101},
            {'pdb_id': '2xyz', 'chain_id': 'B', 'process_id': 102},
            {'pdb_id': '3fail', 'chain_id': 'C', 'process_id': 103}
        ]
        
        service._get_proteins_to_process = Mock(return_value=proteins)
        
        def mock_find_summary(batch_path, pdb_id, chain_id, blast_only=False):
            if pdb_id == '3fail':
                return None  # File not found
            return f"{batch_path}/summaries/{pdb_id}_{chain_id}.develop291.domain_summary.xml"
        
        service._find_domain_summary = Mock(side_effect=mock_find_summary)
        
        def mock_partition_protein(pdb_id, chain_id, summary_path, output_dir, **kwargs):
            if pdb_id == '2xyz':
                return DomainPartitionResult(
                    pdb_id=pdb_id,
                    chain_id=chain_id,
                    reference='develop291',
                    success=False,
                    error="Processing failed"
                )
            return DomainPartitionResult(
                pdb_id=pdb_id,
                chain_id=chain_id,
                reference='develop291',
                success=True,
                is_classified=True
            )
        
        service.partition_protein = Mock(side_effect=mock_partition_protein)
        
        # Process batch
        results = service.partition_batch(123, temp_batch_dir)
        
        # Check results
        assert results.total == 3
        assert results.success_count == 1  # Only 1abc succeeded
        assert results.failure_count == 2  # 2xyz failed processing, 3fail file not found
        assert len(results.failures) == 2
        
        # Check failure reasons
        failure_proteins = [f[0] + '_' + f[1] for f in results.failures]
        assert '2xyz_B' in failure_proteins
        assert '3fail_C' in failure_proteins
    
    def test_sequential_batch_empty(self, service_with_mocks, temp_batch_dir):
        """Test sequential processing of empty batch"""
        service, mock_db = service_with_mocks
        
        service._get_proteins_to_process = Mock(return_value=[])
        
        results = service.partition_batch(123, temp_batch_dir)
        
        assert results.total == 0
        assert results.success_count == 0
        assert results.failure_count == 0
    
    def test_sequential_progress_logging(self, service_with_mocks, temp_batch_dir):
        """Test progress logging during sequential processing"""
        service, mock_db = service_with_mocks
        
        # Create many proteins to trigger progress logging
        proteins = [
            {'pdb_id': f'{i:04d}', 'chain_id': 'A', 'process_id': 100 + i}
            for i in range(25)  # More than progress interval of 10
        ]
        
        service._get_proteins_to_process = Mock(return_value=proteins)
        service._find_domain_summary = Mock(return_value="/fake/summary.xml")
        service.partition_protein = Mock(return_value=DomainPartitionResult(
            pdb_id='test', chain_id='A', reference='develop291', success=True
        ))
        
        with patch('ecod.pipelines.domain_analysis.partition.service.logging') as mock_logging:
            results = service.partition_batch(123, temp_batch_dir)
            
            # Should have logged progress at intervals
            assert results.total == 25
            assert results.success_count == 25


class TestParallelBatchProcessing:
    """Test parallel batch processing"""
    
    @pytest.fixture
    def service_with_parallel(self):
        """Create service configured for parallel processing"""
        mock_context = Mock()
        mock_context.config_manager.config = {
            'partition': {'validation_level': 'normal'},
            'reference': {'current_version': 'develop291'}
        }
        mock_context.config_manager.get_db_config.return_value = {}
        
        service_config = {
            'max_workers': 4,
            'use_multiprocessing': True,
            'batch_size': 10
        }
        
        with patch('ecod.pipelines.domain_analysis.partition.service.DBManager'):
            service = DomainPartitionService(mock_context, service_config)
            return service
    
    def test_parallel_batch_processing_success(self, service_with_parallel):
        """Test successful parallel batch processing"""
        proteins = [
            {'pdb_id': '1abc', 'chain_id': 'A', 'process_id': 101},
            {'pdb_id': '2xyz', 'chain_id': 'B', 'process_id': 102},
            {'pdb_id': '3def', 'chain_id': 'C', 'process_id': 103}
        ]
        
        service_with_parallel._get_proteins_to_process = Mock(return_value=proteins)
        
        # Mock the wrapper function that would be called in parallel
        def mock_wrapper(protein, batch_path, **options):
            return DomainPartitionResult(
                pdb_id=protein['pdb_id'],
                chain_id=protein['chain_id'],
                reference='develop291',
                success=True,
                is_classified=True
            )
        
        service_with_parallel._process_single_protein_wrapper = Mock(side_effect=mock_wrapper)
        
        # Mock ThreadPoolExecutor
        with patch('ecod.pipelines.domain_analysis.partition.service.ThreadPoolExecutor') as mock_executor:
            # Create mock futures
            mock_futures = []
            for protein in proteins:
                mock_future = Mock(spec=Future)
                mock_future.result.return_value = mock_wrapper(protein, "/fake/path")
                mock_futures.append(mock_future)
            
            mock_executor_instance = Mock()
            mock_executor.return_value.__enter__ = Mock(return_value=mock_executor_instance)
            mock_executor.return_value.__exit__ = Mock(return_value=None)
            
            # Setup submit to return our mock futures
            future_iter = iter(mock_futures)
            mock_executor_instance.submit = Mock(side_effect=lambda *args, **kwargs: next(future_iter))
            
            # Mock as_completed to return futures in order
            with patch('ecod.pipelines.domain_analysis.partition.service.as_completed', 
                      return_value=mock_futures):
                results = service_with_parallel.partition_batch(123, "/fake/batch/path")
        
        assert results.success_count == 3
        assert results.failure_count == 0
        assert len(results.results) == 3
    
    def test_parallel_batch_processing_with_failures(self, service_with_parallel):
        """Test parallel processing with some failures"""
        proteins = [
            {'pdb_id': '1abc', 'chain_id': 'A', 'process_id': 101},
            {'pdb_id': '2fail', 'chain_id': 'B', 'process_id': 102}
        ]
        
        service_with_parallel._get_proteins_to_process = Mock(return_value=proteins)
        
        def mock_wrapper(protein, batch_path, **options):
            if protein['pdb_id'] == '2fail':
                raise Exception("Processing failed")
            return DomainPartitionResult(
                pdb_id=protein['pdb_id'],
                chain_id=protein['chain_id'],
                reference='develop291',
                success=True
            )
        
        service_with_parallel._process_single_protein_wrapper = Mock(side_effect=mock_wrapper)
        
        with patch('ecod.pipelines.domain_analysis.partition.service.ThreadPoolExecutor') as mock_executor:
            # Create mock futures - one success, one failure
            success_future = Mock(spec=Future)
            success_future.result.return_value = DomainPartitionResult(
                pdb_id='1abc', chain_id='A', reference='develop291', success=True
            )
            
            failure_future = Mock(spec=Future)
            failure_future.result.side_effect = Exception("Processing failed")
            
            mock_futures = [success_future, failure_future]
            
            mock_executor_instance = Mock()
            mock_executor.return_value.__enter__ = Mock(return_value=mock_executor_instance)
            mock_executor.return_value.__exit__ = Mock(return_value=None)
            
            future_iter = iter(mock_futures)
            mock_executor_instance.submit = Mock(side_effect=lambda *args, **kwargs: next(future_iter))
            
            with patch('ecod.pipelines.domain_analysis.partition.service.as_completed', 
                      return_value=mock_futures):
                results = service_with_parallel.partition_batch(123, "/fake/batch/path")
        
        assert results.success_count == 1
        assert results.failure_count == 1
        assert len(results.failures) == 1
        assert results.failures[0][0] == '2fail'  # Failed protein
    
    def test_parallel_worker_limit(self, service_with_parallel):
        """Test parallel processing respects worker limits"""
        # Test with more proteins than workers
        proteins = [
            {'pdb_id': f'{i:04d}', 'chain_id': 'A', 'process_id': 100 + i}
            for i in range(10)  # More than max_workers (4)
        ]

        service_with_parallel._get_proteins_to_process = Mock(return_value=proteins)
        service_with_parallel._process_single_protein_wrapper = Mock(
            return_value=DomainPartitionResult(
                pdb_id='test', chain_id='A', reference='develop291', success=True
            )
        )

        with patch('ecod.pipelines.domain_analysis.partition.service.ThreadPoolExecutor') as mock_executor:
            # Mock the context manager
            mock_executor_instance = Mock()
            mock_executor.return_value.__enter__ = Mock(return_value=mock_executor_instance)
            mock_executor.return_value.__exit__ = Mock(return_value=None)

            # Create mock futures
            mock_futures = [Mock(spec=Future) for _ in proteins]
            for future in mock_futures:
                future.result.return_value = DomainPartitionResult(
                    pdb_id='test', chain_id='A', reference='develop291', success=True
                )

            # Mock submit to return futures
            future_iter = iter(mock_futures)
            mock_executor_instance.submit = Mock(side_effect=lambda *args, **kwargs: next(future_iter))

            # Mock as_completed to prevent hanging
            with patch('ecod.pipelines.domain_analysis.partition.service.as_completed',
                      return_value=mock_futures):

                service_with_parallel.partition_batch(123, "/fake/batch/path")

                # Should use min(max_workers, len(proteins)) = min(4, 10) = 4
                mock_executor.assert_called_once_with(max_workers=4)


class TestReprocessingFailed:
    """Test reprocessing failed proteins"""
    
    @pytest.fixture
    def service_with_reprocess(self):
        """Create service for reprocessing tests"""
        mock_context = Mock()
        mock_context.config_manager.config = {'reference': {'current_version': 'develop291'}}
        mock_context.config_manager.get_db_config.return_value = {}
        
        with patch('ecod.pipelines.domain_analysis.partition.service.DBManager'):
            service = DomainPartitionService(mock_context)
            return service
    
    def test_reprocess_failed_proteins(self, service_with_reprocess):
        """Test reprocessing only failed proteins"""
        # Mock failed proteins query
        failed_proteins = [
            {'pdb_id': '1fail', 'chain_id': 'A', 'process_id': 201},
            {'pdb_id': '2fail', 'chain_id': 'B', 'process_id': 202}
        ]
        
        service_with_reprocess.db.execute_dict_query.return_value = failed_proteins
        
        # Mock the partition_batch method to track reprocessing
        mock_batch_result = BatchPartitionResults()
        mock_batch_result.total = 2
        mock_batch_result.success_count = 2
        
        service_with_reprocess.partition_batch = Mock(return_value=mock_batch_result)
        
        # Call reprocess_failed
        results = service_with_reprocess.reprocess_failed(123, "/fake/batch/path")
        
        # Should have queried for failed proteins
        service_with_reprocess.db.execute_dict_query.assert_called_once()
        query_args = service_with_reprocess.db.execute_dict_query.call_args[0]
        assert 'domain_partition_failed' in query_args[0]
        assert 'domain_partition_error' in query_args[0]
        assert query_args[1] == (123,)
        
        # Should have called partition_batch with force_overwrite
        service_with_reprocess.partition_batch.assert_called_once_with(
            123, "/fake/batch/path", force_overwrite=True
        )
        
        assert results.success_count == 2
    
    def test_reprocess_failed_no_failures(self, service_with_reprocess):
        """Test reprocessing when no failed proteins exist"""
        service_with_reprocess.db.execute_dict_query.return_value = []
        
        results = service_with_reprocess.reprocess_failed(123, "/fake/batch/path")
        
        assert results.total == 0
    
    def test_reprocess_failed_db_error(self, service_with_reprocess):
        """Test reprocessing with database error"""
        service_with_reprocess.db.execute_dict_query.side_effect = Exception("DB Error")
        
        results = service_with_reprocess.reprocess_failed(123, "/fake/batch/path")
        
        assert results.total == 0  # Empty result on error


class TestFileSummaryDiscovery:
    """Test domain summary file discovery"""
    
    @pytest.fixture
    def temp_batch_structure(self):
        """Create temporary batch directory with realistic structure"""
        temp_dir = tempfile.mkdtemp()
        base_path = Path(temp_dir)
        
        # Create typical batch structure
        summaries_dir = base_path / "summaries"
        summaries_dir.mkdir(parents=True)
        
        blast_summaries_dir = base_path / "blast_summaries"
        blast_summaries_dir.mkdir(parents=True)
        
        # Create some test files
        (summaries_dir / "1abc_A.develop291.domain_summary.xml").touch()
        (summaries_dir / "2xyz_B.develop291.domain_summary.xml").touch()
        (blast_summaries_dir / "1abc_A.develop291.blast_only_summary.xml").touch()
        
        yield temp_dir
        shutil.rmtree(temp_dir)
    
    @pytest.fixture
    def service_for_files(self):
        """Create service for file discovery tests"""
        mock_context = Mock()
        mock_context.config_manager.config = {'reference': {'current_version': 'develop291'}}
        mock_context.config_manager.get_db_config.return_value = {}
        
        with patch('ecod.pipelines.domain_analysis.partition.service.DBManager'):
            service = DomainPartitionService(mock_context)
            return service
    
    def test_find_domain_summary_normal(self, service_for_files, temp_batch_structure):
        """Test finding normal domain summary"""
        with patch('ecod.utils.path_utils.get_all_evidence_paths') as mock_get_paths:
            mock_get_paths.return_value = {
                'domain_summary': {
                    'exists_at': f"{temp_batch_structure}/summaries/1abc_A.develop291.domain_summary.xml"
                }
            }
            
            result = service_for_files._find_domain_summary(
                temp_batch_structure, "1abc", "A"
            )
            
            assert result is not None
            assert "1abc_A.develop291.domain_summary.xml" in result
    
    def test_find_domain_summary_blast_only(self, service_for_files, temp_batch_structure):
        """Test finding BLAST-only summary"""
        with patch('ecod.utils.path_utils.get_all_evidence_paths') as mock_get_paths:
            mock_get_paths.return_value = {
                'blast_only_summary': {
                    'exists_at': f"{temp_batch_structure}/blast_summaries/1abc_A.develop291.blast_only_summary.xml"
                }
            }
            
            result = service_for_files._find_domain_summary(
                temp_batch_structure, "1abc", "A", blast_only=True
            )
            
            assert result is not None
            assert "blast_only_summary.xml" in result
    
    def test_find_domain_summary_not_found(self, service_for_files, temp_batch_structure):
        """Test handling missing summary file"""
        with patch('ecod.utils.path_utils.get_all_evidence_paths') as mock_get_paths:
            mock_get_paths.return_value = {
                'domain_summary': {
                    'exists_at': None
                }
            }
            
            result = service_for_files._find_domain_summary(
                temp_batch_structure, "missing", "X"
            )
            
            assert result is None
    
    def test_find_domain_summary_path_error(self, service_for_files, temp_batch_structure):
        """Test handling path utility error"""
        with patch('ecod.utils.path_utils.get_all_evidence_paths') as mock_get_paths:
            mock_get_paths.side_effect = Exception("Path error")
            
            result = service_for_files._find_domain_summary(
                temp_batch_structure, "1abc", "A"
            )
            
            assert result is None


class TestBatchStatusTracking:
    """Test batch status tracking and database updates"""
    
    @pytest.fixture
    def service_with_tracker(self):
        """Create service with status tracking"""
        mock_context = Mock()
        mock_context.config_manager.config = {'reference': {'current_version': 'develop291'}}
        mock_context.config_manager.get_db_config.return_value = {}
        
        with patch('ecod.pipelines.domain_analysis.partition.service.DBManager'):
            service = DomainPartitionService(mock_context)
            service.service_settings['track_status'] = True
            return service
    
    def test_batch_completion_status_update(self, service_with_tracker):
        """Test batch completion status gets updated"""
        proteins = [
            {'pdb_id': '1abc', 'chain_id': 'A', 'process_id': 101}
        ]
        
        service_with_tracker._get_proteins_to_process = Mock(return_value=proteins)
        service_with_tracker._find_domain_summary = Mock(return_value="/fake/summary.xml")
        service_with_tracker.partition_protein = Mock(return_value=DomainPartitionResult(
            pdb_id='1abc', chain_id='A', reference='develop291', success=True
        ))
        
        # Process batch
        service_with_tracker.partition_batch(123, "/fake/batch/path")
        
        # Should have called tracker to update batch status
        service_with_tracker.tracker.update_batch_completion_status.assert_called_once_with(
            123, representatives_only=False
        )
    
    def test_non_representative_status_update(self, service_with_tracker):
        """Test non-representative status update"""
        proteins = [
            {'pdb_id': '1abc', 'chain_id': 'A', 'process_id': 101, 'is_representative': True}
        ]
        
        service_with_tracker._get_proteins_to_process = Mock(return_value=proteins)
        service_with_tracker._find_domain_summary = Mock(return_value="/fake/summary.xml")
        service_with_tracker.partition_protein = Mock(return_value=DomainPartitionResult(
            pdb_id='1abc', chain_id='A', reference='develop291', success=True
        ))
        
        # Process batch with representatives only
        service_with_tracker.partition_batch(123, "/fake/batch/path", representatives_only=True)
        
        # Should have updated non-representative status
        service_with_tracker.tracker.update_non_representative_status.assert_called_once_with(123)
    
    def test_process_status_tracking_disabled(self, service_with_tracker):
        """Test when status tracking is disabled"""
        service_with_tracker.service_settings['track_status'] = False
        
        proteins = [{'pdb_id': '1abc', 'chain_id': 'A', 'process_id': 101}]
        service_with_tracker._get_proteins_to_process = Mock(return_value=proteins)
        service_with_tracker._find_domain_summary = Mock(return_value="/fake/summary.xml")
        service_with_tracker.partition_protein = Mock(return_value=DomainPartitionResult(
            pdb_id='1abc', chain_id='A', reference='develop291', success=True
        ))
        
        service_with_tracker.partition_batch(123, "/fake/batch/path")
        
        # Should not have called status updates
        service_with_tracker.tracker.update_non_representative_status.assert_not_called()


class TestBatchResultsManagement:
    """Test BatchPartitionResults functionality"""
    
    def test_batch_results_initialization(self):
        """Test batch results initialization"""
        results = BatchPartitionResults()
        
        assert results.total == 0
        assert results.success_count == 0
        assert results.failure_count == 0
        assert results.success_rate == 0.0
        assert results.start_time is not None
        assert results.end_time is None
    
    def test_batch_results_add_success(self):
        """Test adding successful results"""
        results = BatchPartitionResults()
        
        success_result = DomainPartitionResult(
            pdb_id='1abc', chain_id='A', reference='develop291',
            success=True, is_classified=True,
            domains=[{'id': 'd1', 'start': 10, 'end': 100}]
        )
        
        results.add_result(success_result)
        
        assert results.total == 1
        assert results.success_count == 1
        assert results.failure_count == 0
        assert results.success_rate == 100.0
        assert results.proteins_with_domains == 1
        assert results.total_domains_found == 1
    
    def test_batch_results_add_failure(self):
        """Test adding failed results"""
        results = BatchPartitionResults()
        
        failure_result = DomainPartitionResult(
            pdb_id='1fail', chain_id='A', reference='develop291',
            success=False, error="Processing failed"
        )
        
        results.add_result(failure_result)
        
        assert results.total == 1
        assert results.success_count == 0
        assert results.failure_count == 1
        assert results.success_rate == 0.0
        assert len(results.failures) == 1
        assert results.failures[0] == ('1fail', 'A', 'Processing failed')
    
    def test_batch_results_add_peptide(self):
        """Test adding peptide results"""
        results = BatchPartitionResults()
        
        peptide_result = DomainPartitionResult(
            pdb_id='1pep', chain_id='A', reference='develop291',
            success=True, is_peptide=True
        )
        
        results.add_result(peptide_result)
        
        assert results.peptides_found == 1
        assert results.proteins_with_domains == 0
    
    def test_batch_results_mixed(self):
        """Test batch results with mixed outcomes"""
        results = BatchPartitionResults()
        
        # Add various result types
        results.add_result(DomainPartitionResult(
            pdb_id='1abc', chain_id='A', reference='develop291',
            success=True, is_classified=True,
            domains=[{'id': 'd1', 'start': 10, 'end': 100}]
        ))
        
        results.add_result(DomainPartitionResult(
            pdb_id='2pep', chain_id='B', reference='develop291',
            success=True, is_peptide=True
        ))
        
        results.add_result(DomainPartitionResult(
            pdb_id='3unc', chain_id='C', reference='develop291',
            success=True, is_unclassified=True
        ))
        
        results.add_result(DomainPartitionResult(
            pdb_id='4fail', chain_id='D', reference='develop291',
            success=False, error="Error"
        ))
        
        assert results.total == 4
        assert results.success_count == 3
        assert results.failure_count == 1
        assert results.proteins_with_domains == 1
        assert results.peptides_found == 1
        assert results.unclassified_proteins == 1
        assert results.success_rate == 75.0
    
    def test_batch_results_finalization(self):
        """Test batch results finalization"""
        results = BatchPartitionResults()
        
        # Add some results
        results.add_result(DomainPartitionResult(
            pdb_id='1abc', chain_id='A', reference='develop291', success=True
        ))
        
        assert results.end_time is None
        
        results.finalize()
        
        assert results.end_time is not None
        assert results.processing_time > 0
    
    def test_batch_results_summary(self):
        """Test batch results summary generation"""
        results = BatchPartitionResults()
        
        results.add_result(DomainPartitionResult(
            pdb_id='1abc', chain_id='A', reference='develop291',
            success=True, domains=[{'id': 'd1'}]
        ))
        
        results.finalize()
        summary = results.get_summary()
        
        assert summary['total_proteins'] == 1
        assert summary['successful'] == 1
        assert summary['failed'] == 0
        assert summary['success_rate'] == 100.0
        assert 'processing_time' in summary
        assert 'average_time_per_protein' in summary


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
