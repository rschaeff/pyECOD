#!/usr/bin/env python3
"""
Basic import and functionality tests.

Quick smoke tests that verify core components work without requiring
external dependencies like pytest-cov.
"""

def test_evidence_imports():
    """Test that Evidence model imports and creates correctly"""
    from ecod.models.pipeline.evidence import Evidence
    
    # Test basic creation
    evidence = Evidence(type='blast', domain_id='d1abcA1')
    assert evidence.type == 'blast'
    assert evidence.domain_id == 'd1abcA1'
    
    # Test confidence auto-calculation
    evidence_with_evalue = Evidence(type='blast', evalue=1e-20, confidence=None)
    assert evidence_with_evalue.confidence > 0.0
    assert evidence_with_evalue.confidence <= 1.0


def test_domain_model_imports():
    """Test that DomainModel imports and creates correctly"""
    from ecod.models.pipeline.domain import DomainModel
    
    # Test basic creation
    domain = DomainModel(
        id='test_domain',
        start=10,
        end=100,
        range='10-100'
    )
    
    assert domain.id == 'test_domain'
    assert domain.start == 10
    assert domain.end == 100
    assert domain.size == 91


def test_partition_result_imports():
    """Test that DomainPartitionResult imports and creates correctly"""
    from ecod.models.pipeline.partition import DomainPartitionResult
    
    # Test basic creation
    result = DomainPartitionResult(
        pdb_id='1abc',
        chain_id='A',
        reference='develop291'
    )
    
    assert result.pdb_id == '1abc'
    assert result.chain_id == 'A'
    assert result.reference == 'develop291'
    assert result.success == True


def test_partition_options_imports():
    """Test that PartitionOptions imports and validates correctly"""
    from ecod.pipelines.domain_analysis.partition.models import PartitionOptions
    
    # Test basic creation
    options = PartitionOptions()
    assert options.min_domain_size > 0
    assert 0.0 <= options.overlap_threshold <= 1.0
    
    # Test validation
    options.validate()  # Should not raise


def test_evidence_analyzer_imports():
    """Test that EvidenceAnalyzer imports correctly"""
    from ecod.pipelines.domain_analysis.partition.analyzer import EvidenceAnalyzer
    from ecod.pipelines.domain_analysis.partition.models import PartitionOptions
    
    # Test basic creation
    options = PartitionOptions()
    analyzer = EvidenceAnalyzer(options)
    
    assert analyzer.options == options
    assert hasattr(analyzer, 'stats')


def test_partition_processor_imports():
    """Test that PartitionProcessor imports correctly"""
    from ecod.pipelines.domain_analysis.partition.processor import PartitionProcessor
    from ecod.pipelines.domain_analysis.partition.analyzer import EvidenceAnalyzer
    from ecod.pipelines.domain_analysis.partition.models import PartitionOptions
    
    # Test basic creation
    options = PartitionOptions()
    analyzer = EvidenceAnalyzer(options)
    processor = PartitionProcessor(options, analyzer)
    
    assert processor.options == options
    assert processor.analyzer == analyzer


def test_status_tracker_imports():
    """Test that StatusTracker imports correctly"""
    from ecod.pipelines.domain_analysis.partition.tracker import StatusTracker
    from unittest.mock import Mock
    
    # Test basic creation with mock DB
    mock_db = Mock()
    tracker = StatusTracker(mock_db)
    
    assert tracker.db == mock_db


def test_service_imports():
    """Test that DomainPartitionService imports correctly"""
    from ecod.pipelines.domain_analysis.partition.service import DomainPartitionService
    
    # Just test import - actual initialization requires more setup
    assert DomainPartitionService is not None


def test_convenience_functions():
    """Test that convenience functions import correctly"""
    from ecod.pipelines.domain_analysis.partition.service import (
        create_service, partition_single_protein, partition_batch
    )
    
    # Just test imports
    assert create_service is not None
    assert partition_single_protein is not None
    assert partition_batch is not None


if __name__ == '__main__':
    """Run basic tests without pytest"""
    import sys
    
    test_functions = [
        test_evidence_imports,
        test_domain_model_imports,
        test_partition_result_imports,
        test_partition_options_imports,
        test_evidence_analyzer_imports,
        test_partition_processor_imports,
        test_status_tracker_imports,
        test_service_imports,
        test_convenience_functions
    ]
    
    passed = 0
    failed = 0
    
    print("Running basic import and functionality tests...")
    print("=" * 50)
    
    for test_func in test_functions:
        try:
            test_func()
            print(f"✅ {test_func.__name__}")
            passed += 1
        except Exception as e:
            print(f"❌ {test_func.__name__}: {str(e)}")
            failed += 1
    
    print("=" * 50)
    print(f"Results: {passed} passed, {failed} failed")
    
    if failed > 0:
        print("\nThere are import or basic functionality issues.")
        print("Please check your Python path and module structure.")
        sys.exit(1)
    else:
        print("\n🎉 All basic tests passed! Core modules are working.")
        print("\nNext steps:")
        print("1. Install test dependencies: pip install pytest pytest-cov pytest-html")
        print("2. Run comprehensive tests: python test_runner_config.py --coverage")
        sys.exit(0)
