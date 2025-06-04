# mini/models.py
"""Minimal models for domain analysis"""

from dataclasses import dataclass
from typing import Optional, List, Any
from ecod.core.sequence_range import SequenceRange

@dataclass
class AlignmentData:
    """Alignment information for decomposition"""
    query_seq: str
    hit_seq: str
    query_start: int
    query_end: int
    hit_start: int
    hit_end: int


@dataclass
class Evidence:
    """Minimal evidence model"""
    type: str  # 'chain_blast', 'domain_blast', 'hhsearch'
    source_pdb: str  # '6dgv', '2ia4', etc.
    query_range: SequenceRange
    confidence: float = 0.0
    evalue: Optional[float] = None
    domain_id: Optional[str] = None  # 'e6dgvA1'

    # Classification for better partitioning
    t_group: Optional[str] = None
    h_group: Optional[str] = None

    # Reference info for coverage calculation
    reference_length: Optional[int] = None  # Actual length from reference data
    alignment_coverage: Optional[float] = None  # Only set if reference_length available

    # Optional alignment data (for chain BLAST decomposition)
    alignment: Optional[Any] = None

@dataclass  
class Domain:
    """Minimal domain model"""
    id: str
    range: SequenceRange
    family: str  # PDB or classification
    evidence_count: int
    source: str  # 'blast' or 'hhsearch'
    evidence_items: List[Evidence]

    #ECOD hierarchy fields
    x_group: Optional[str] = None
    h_group: Optional[str] = None
    t_group: Optional[str] = None
