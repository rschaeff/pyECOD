#!/usr/bin/env python3
"""
Domain partition module - New service-oriented architecture

This module provides domain partitioning functionality through a clean,
modular service architecture that replaces the old monolithic approach.
"""

# Import the new service components
from .service import DomainPartitionService, create_service, partition_single_protein, partition_batch
from .processor import PartitionProcessor
from .analyzer import EvidenceAnalyzer
from .tracker import StatusTracker
from .models import (
    PartitionOptions, ValidationLevel, ProcessingMode, PartitionContext,
    ValidationResult, EvidenceGroup, DomainCandidate, BatchPartitionResults,
    PartitionStage, ClassificationCache
)

# Backward compatibility - import the old class but mark as deprecated
from ..partition_legacy import DomainPartition

import warnings

# Main public API
__all__ = [
    # Primary service interface
    'DomainPartitionService',
    
    # Core components
    'PartitionProcessor',
    'EvidenceAnalyzer', 
    'StatusTracker',
    
    # Configuration and models
    'PartitionOptions',
    'ValidationLevel',
    'ProcessingMode',
    'PartitionContext',
    'ValidationResult',
    'EvidenceGroup',
    'DomainCandidate',
    'BatchPartitionResults',
    'PartitionStage',
    'ClassificationCache',
    
    # Convenience functions
    'create_service',
    'partition_single_protein',
    'partition_batch',
    
    # Backward compatibility (deprecated)
    'DomainPartition',
]

# Issue deprecation warning when old class is accessed
def __getattr__(name):
    if name == 'DomainPartition':
        warnings.warn(
            "DomainPartition is deprecated. Use DomainPartitionService instead. "
            "See the migration guide for details.",
            DeprecationWarning,
            stacklevel=2
        )
        return DomainPartition
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

# Version info
__version__ = "2.0.0"
__author__ = "ECOD Team"

# Quick access patterns for common use cases
def quick_partition_protein(pdb_id: str, chain_id: str, summary_path: str, 
                          output_dir: str, config_path: str = None, **options):
    """
    Quick function to partition a single protein with sensible defaults.
    
    This is the recommended way for simple use cases.
    """
    return partition_single_protein(pdb_id, chain_id, summary_path, output_dir, 
                                  config_path, **options)

def quick_partition_batch(batch_id: int, batch_path: str, 
                         config_path: str = None, **options):
    """
    Quick function to partition a batch with sensible defaults.
    
    This is the recommended way for batch processing.
    """
    return partition_batch(batch_id, batch_path, config_path, **options)

# Service factory with common configurations
def create_standard_service(config_path: str = None):
    """Create service with standard production settings"""
    service_config = {
        'max_workers': 4,
        'use_multiprocessing': False,  # Conservative default
        'batch_size': 100,
        'save_intermediate': True,
        'track_status': True
    }
    return create_service(config_path, service_config)

def create_fast_service(config_path: str = None):
    """Create service optimized for speed"""
    service_config = {
        'max_workers': 8,
        'use_multiprocessing': True,
        'batch_size': 50,
        'save_intermediate': False,
        'track_status': False
    }
    return create_service(config_path, service_config)

def create_reliable_service(config_path: str = None):
    """Create service optimized for reliability"""
    service_config = {
        'max_workers': 2,
        'use_multiprocessing': False,
        'batch_size': 25,
        'save_intermediate': True,
        'track_status': True
    }
    return create_service(config_path, service_config)
