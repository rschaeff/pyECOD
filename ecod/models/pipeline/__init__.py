# ecod/models/pipeline/__init__.py
"""
Gold Standard Pipeline Models for pyECOD

This module contains the consolidated, gold standard models for domain analysis
and processing workflows. These models replace all previous implementations.
"""

# Gold Standard Models - Primary Exports
from .evidence import Evidence, DomainEvidence, BlastEvidence, HHSearchEvidence
from .domain import DomainModel, Domain
from .partition import DomainPartitionResult, PartitionResult

__all__ = [
    # Evidence Models (Unified)
    'Evidence',           # Primary evidence model
    'DomainEvidence',     # Alias for Evidence
    'BlastEvidence',      # Alias for Evidence
    'HHSearchEvidence',   # Alias for Evidence

    # Domain Models (Consolidated)
    'DomainModel',        # Primary domain model
    'Domain',             # Alias for DomainModel

    # Partition Results (Enhanced)
    'DomainPartitionResult',  # Primary partition result model
    'PartitionResult',        # Alias for DomainPartitionResult
]

# Model metadata
__version__ = "2.0.0"
__status__ = "gold_standard"

# Quick access to model classes
EVIDENCE_MODEL = Evidence
DOMAIN_MODEL = DomainModel
PARTITION_MODEL = DomainPartitionResult

def get_model_info():
    """Get information about the current pipeline models"""
    return {
        "version": __version__,
        "status": __status__,
        "evidence_model": Evidence.__name__,
        "domain_model": DomainModel.__name__,
        "partition_model": DomainPartitionResult.__name__,
        "available_models": __all__
    }
