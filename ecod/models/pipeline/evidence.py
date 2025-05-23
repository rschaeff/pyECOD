# models/pipeline/evidence.py
from dataclasses import dataclass, field
from typing import Dict, Any, Optional, List, Union
import xml.etree.ElementTree as ET
from ecod.models.base import XmlSerializable

@dataclass
class Evidence(XmlSerializable):
    """Unified evidence model - consolidates DomainEvidence and Evidence classes"""
    type: str  # "blast", "hhsearch", "domain_blast", "chain_blast"
    source_id: str = ""
    domain_id: str = ""
    query_range: str = ""
    hit_range: str = ""
    confidence: float = 0.0
    
    # Classification data
    t_group: Optional[str] = None
    h_group: Optional[str] = None
    x_group: Optional[str] = None
    a_group: Optional[str] = None
    
    # Source-specific attributes (replaces both attributes dict and specific fields)
    evalue: Optional[float] = None
    probability: Optional[float] = None
    score: Optional[float] = None
    identity: Optional[float] = None
    coverage: Optional[float] = None
    hsp_count: Optional[int] = None
    
    # Additional attributes for extensibility
    extra_attributes: Dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self):
        """Auto-calculate confidence from available metrics"""
        if self.confidence == 0.0:
            if self.probability is not None and self.probability > 1.0:
                # HHSearch probability (0-100) to confidence (0-1)
                self.confidence = self.probability / 100.0
            elif self.evalue is not None:
                # E-value to confidence
                self.confidence = 1.0 / (1.0 + self.evalue)

    @classmethod
    def from_blast_hit(cls, hit, hit_type="domain_blast") -> 'Evidence':
        """Create from BlastHit - replaces BlastEvidence.from_blast_hit"""
        return cls(
            type=hit_type,
            source_id=getattr(hit, "domain_id", "") or getattr(hit, "hit_id", ""),
            domain_id=getattr(hit, "domain_id", ""),
            query_range=getattr(hit, "range", ""),
            hit_range=getattr(hit, "hit_range", ""),
            evalue=getattr(hit, "evalue", 999.0),
            hsp_count=getattr(hit, "hsp_count", 0),
            extra_attributes={
                "pdb_id": getattr(hit, "pdb_id", ""),
                "chain_id": getattr(hit, "chain_id", ""),
                "discontinuous": getattr(hit, "discontinuous", False)
            }
        )
    
    @classmethod
    def from_hhsearch_hit(cls, hit) -> 'Evidence':
        """Create from HHSearchHit - replaces HHSearchEvidence.from_hhsearch_hit"""
        return cls(
            type="hhsearch",
            source_id=getattr(hit, "domain_id", "") or getattr(hit, "hit_id", ""),
            domain_id=getattr(hit, "domain_id", ""),
            query_range=getattr(hit, "range", ""),
            hit_range=getattr(hit, "hit_range", ""),
            probability=getattr(hit, "probability", 0.0),
            evalue=getattr(hit, "evalue", 999.0),
            score=getattr(hit, "score", 0.0)
        )

# Backward compatibility aliases
DomainEvidence = Evidence
BlastEvidence = Evidence 
HHSearchEvidence = Evidence
