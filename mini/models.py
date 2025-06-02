# mini_pyecod/models.py
"""Minimal models for domain analysis"""

from dataclasses import dataclass
from typing import Optional, List
from .sequence_range import SequenceRange

@dataclass
class Evidence:
    """Minimal evidence model"""
    type: str  # 'chain_blast', 'domain_blast', 'hhsearch'
    source_pdb: str  # '6dgv', '2ia4', etc.
    query_range: SequenceRange
    confidence: float = 0.0
    evalue: Optional[float] = None
    domain_id: Optional[str] = None  # 'e6dgvA1'

@dataclass  
class Domain:
    """Minimal domain model"""
    id: str
    range: SequenceRange
    family: str  # PDB or classification
    evidence_count: int
    source: str  # 'blast' or 'hhsearch'
    evidence_items: List[Evidence]
