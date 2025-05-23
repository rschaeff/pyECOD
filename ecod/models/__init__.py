#!/usr/bin/env python3
"""
ECOD Pipeline Models Module - Cleaned and Consolidated

This module provides access to the gold standard models for the ECOD pipeline.
All models have been consolidated under the pipeline namespace for consistency.

IMPORTANT: BlastHit and HHSearchHit have been consolidated into Evidence model.
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

# Legacy Pipeline Models - Import from pipeline.py file (not pipeline/ directory)
# Import using importlib to avoid name conflict with pipeline/ directory
import importlib.util
import os

def _import_legacy_pipeline_models():
    """Import legacy models from pipeline.py file"""
    pipeline_file = os.path.join(os.path.dirname(__file__), 'pipeline.py')

    if os.path.exists(pipeline_file):
        spec = importlib.util.spec_from_file_location("pipeline_legacy", pipeline_file)
        if spec and spec.loader:
            pipeline_module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(pipeline_module)
            return pipeline_module

    return None

# Import legacy models
_pipeline_legacy = _import_legacy_pipeline_models()
if _pipeline_legacy:
    DomainSummaryModel = _pipeline_legacy.DomainSummaryModel
    PipelineResult = _pipeline_legacy.PipelineResult
    ProteinResult = _pipeline_legacy.ProteinResult
    ProteinProcessingResult = _pipeline_legacy.ProteinProcessingResult
else:
    # Fallback stubs if file not found
    class DomainSummaryModel:
        pass
    class PipelineResult:
        pass
    class ProteinResult:
        pass
    class ProteinProcessingResult:
        pass

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

    # Legacy Pipeline Models (Note: BlastHit, HHSearchHit consolidated into Evidence)
    'DomainSummaryModel', 'PipelineResult',
    'ProteinResult', 'ProteinProcessingResult'
]

# Version info for model compatibility
__version__ = "2.0.0"
__model_version__ = "consolidated"

# Migration notes for developers
_MIGRATION_NOTES = """
Model Consolidation - Breaking Changes:

CONSOLIDATED MODELS:
- models.pipeline.BlastHit -> Use models.pipeline.evidence.Evidence
- models.pipeline.HHSearchHit -> Use models.pipeline.evidence.Evidence
- models.evidence.DomainEvidence -> Use models.pipeline.evidence.Evidence
- models.domain.DomainModel -> Use models.pipeline.domain.DomainModel
- models.domain_analysis.* -> All replaced by pipeline models
- utils.evidence_bridge -> No longer needed

NEW GOLD STANDARD:
- Evidence: Unified evidence model (replaces BlastHit, HHSearchHit, all evidence types)
- DomainModel: Consolidated domain model with full functionality
- DomainPartitionResult: Enhanced partition result model

USAGE:
from ecod.models import Evidence, DomainModel, DomainPartitionResult
# OR
from ecod.models.pipeline import Evidence, DomainModel, DomainPartitionResult

MIGRATION EXAMPLES:
# Old:
blast_hit = BlastHit.from_xml(element)
evidence = Evidence.from_blast_hit(blast_hit)

# New:
evidence = Evidence.from_blast_xml(element, "domain_blast")

# Old:
hhsearch_hit = HHSearchHit.from_xml(element)
evidence = Evidence.from_hhsearch_hit(hhsearch_hit)

# New:
evidence = Evidence.from_hhsearch_xml(element)
"""
