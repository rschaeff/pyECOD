# ecod/evaluation/__init__.py
"""
ECOD Evaluation Framework

This module provides tools for algorithm version management, test set creation,
and systematic evaluation of domain partitioning algorithms.
"""

# For now, just provide access to the main classes
try:
    from .algorithm_versions import AlgorithmVersionManager, AlgorithmVersion, AlgorithmStatus
except ImportError:
    # Handle case where algorithm_versions module doesn't exist yet
    pass

__all__ = ['AlgorithmVersionManager', 'AlgorithmVersion', 'AlgorithmStatus']
