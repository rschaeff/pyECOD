#!/usr/bin/env python3
"""
ECOD Pipeline Models Module - Cleaned and Consolidated

This module provides access to the gold standard models for the ECOD pipeline.
All models have been consolidated under the pipeline namespace for consistency.
"""

# Gold Standard Pipeline Models (Primary)
from ecod.models.pipeline.evidence import Evidence, DomainEvidence, BlastEvidence, HHSearchEvidence
from ecod.models.pipeline.domain import DomainModel, Domain
from ecod.models.pipeline.partition import DomainPartitionResult, PartitionResult

# Core Infrastructure Models
from .base import XmlSerializable

# Protein and Structure Models
from .protein import (
    Protein, ProteinSequence, ProteinStructure,
    PDBChain, ChainSequence, PDBEntry
)

# Job and Workflow Models
from .job import (
    Batch, ProcessStatus, ProcessFile,
    Job, JobItem, ECODVersion, ReferenceResource
)

# Legacy Pipeline Models (for existing pipeline code)
from .pipeline import (
    BlastHit, HHSearchHit, DomainSummaryModel, PipelineResult,
    ProteinResult, ProteinProcessingResult
)

__all__ = [
    # Gold Standard Models
    'Evidence', 'DomainEvidence', 'BlastEvidence', 'HHSearchEvidence',
    'DomainModel', 'Domain',
    'DomainPartitionResult', 'PartitionResult',

    # Core Infrastructure
    'XmlSerializable',

    # Protein Models
    'Protein', 'ProteinSequence', 'ProteinStructure',
    'PDBChain', 'ChainSequence', 'PDBEntry',

    # Job Models
    'Batch', 'ProcessStatus', 'ProcessFile',
    'Job', 'JobItem', 'ECODVersion', 'ReferenceResource',

    # Legacy Pipeline Models
    'BlastHit', 'HHSearchHit', 'DomainSummaryModel', 'PipelineResult',
    'ProteinResult', 'ProteinProcessingResult'
]

# Version info for model compatibility
__version__ = "2.0.0"
__model_version__ = "consolidated"

# Migration notes for developers
_MIGRATION_NOTES = """
Model Consolidation - Breaking Changes:

REMOVED MODELS:
- models.evidence.DomainEvidence -> Use models.pipeline.evidence.Evidence
- models.domain.DomainModel -> Use models.pipeline.domain.DomainModel
- models.domain_analysis.* -> All replaced by pipeline models
- utils.evidence_bridge -> No longer needed

NEW GOLD STANDARD:
- Evidence: Unified evidence model (replaces DomainEvidence, BlastEvidence, HHSearchEvidence)
- DomainModel: Consolidated domain model with full functionality
- DomainPartitionResult: Enhanced partition result model

USAGE:
from ecod.models import Evidence, DomainModel, DomainPartitionResult
# OR
from ecod.models.pipeline import Evidence, DomainModel, DomainPartitionResult
"""
