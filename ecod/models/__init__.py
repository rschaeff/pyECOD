#!/usr/bin/env python3
"""
ECOD Pipeline Models Module
"""

# models/__init__.py - Add temporary imports
from ecod.models.pipeline.evidence import Evidence as DomainEvidence
from ecod.models.pipeline.domain import DomainModel
from ecod.models.pipeline.partition import DomainPartitionResult

# Keep existing imports working
from ecod.models.evidence import DomainEvidence as LegacyDomainEvidence  # deprecate later

from .protein import (
    Protein, ProteinSequence, ProteinStructure, 
    PDBChain, ChainSequence, PDBEntry
)
from .domain import (
    DomainSequence, DomainDSSPDetail,
    DomainRange, DomainRangeSegment, DomainClassification
)
from .job import (
    Batch, ProcessStatus, ProcessFile, 
    Job, JobItem, ECODVersion, ReferenceResource
)

from .pipeline import (
    BlastHit, HHSearchHit, DomainSummaryModel, PipelineResult,
    ProteinResult, ProteinProcessingResult
    )

from .domain_analysis.partition_result import (
    LegacyDomainPartitionResult
    )



__all__ = [
    'Protein', 'ProteinSequence', 'ProteinStructure',
    'PDBChain', 'ChainSequence', 'PDBEntry',
    'Domain', 'DomainSequence', 'DomainDSSPDetail',
    'DomainRange', 'DomainRangeSegment', 'DomainClassification',
    'Batch', 'ProcessStatus', 'ProcessFile',
    'Job', 'JobItem', 'ECODVersion', 'ReferenceResource', 'BlastHHit',
    "HHSearchHit", "DomainSummaryModel", "PipelineResult", "ProteinResult",
    "ProteinProcessingResult", "DomainPartitionResult"
]
