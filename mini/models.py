# mini/models.py
"""Minimal models for domain analysis"""

from dataclasses import dataclass
from typing import Optional, List
from ecod.core.sequence_range import SequenceRange  # Use the real one!

@dataclass
class Evidence:
    """Minimal evidence model"""
    type: str  # 'chain_blast', 'domain_blast', 'hhsearch'
    source_pdb: str  # '6dgv', '2ia4', etc.
    query_range: SequenceRange
    confidence: float = 0.0
    evalue: Optional[float] = None
    domain_id: Optional[str] = None  # 'e6dgvA1'

    # Add classification for better partitioning
    t_group: Optional[str] = None
    h_group: Optional[str] = None

@dataclass  
class Domain:
    """Minimal domain model"""
    id: str
    range: SequenceRange
    family: str  # PDB or classification
    evidence_count: int
    source: str  # 'blast' or 'hhsearch'
    evidence_items: List[Evidence]
