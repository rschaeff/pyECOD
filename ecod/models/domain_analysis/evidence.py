# ecod/models/domain_analysis/evidence.py
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional
import xml.etree.ElementTree as ET
from ecod.models.base import XmlSerializable

@dataclass
class Evidence(XmlSerializable):
    """Base evidence model for domain boundary decisions"""
    type: str  # "blast", "hhsearch", "self_comparison"
    source_id: str = ""  # Domain ID, hit ID, etc.
    query_range: str = ""
    hit_range: str = ""
    confidence: float = 0.0  # Normalized confidence score (0-1)
    
    # Classification info (if available)
    t_group: Optional[str] = None
    h_group: Optional[str] = None
    x_group: Optional[str] = None
    a_group: Optional[str] = None
    
    # Source-specific attributes
    attributes: Dict[str, Any] = field(default_factory=dict)
    
    def to_xml(self) -> ET.Element:
        """Convert to XML Element"""
        element = ET.Element("evidence")
        element.set("type", self.type)
        element.set("source_id", self.source_id)
        
        # Add ranges
        query_range = ET.SubElement(element, "query_range")
        query_range.text = self.query_range
        
        hit_range = ET.SubElement(element, "hit_range")
        hit_range.text = self.hit_range
        
        # Add confidence
        element.set("confidence", f"{self.confidence:.4f}")
        
        # Add classification if available
        for cls_attr in ["t_group", "h_group", "x_group", "a_group"]:
            value = getattr(self, cls_attr)
            if value:
                element.set(cls_attr, value)
        
        # Add source-specific attributes
        for key, value in self.attributes.items():
            if value is not None:
                element.set(key, str(value))
        
        return element
    
    @classmethod
    def from_xml(cls, element: ET.Element) -> 'Evidence':
        """Create from XML Element"""
        evidence = cls(
            type=element.get("type", ""),
            source_id=element.get("source_id", ""),
            confidence=float(element.get("confidence", "0.0"))
        )
        
        # Get ranges
        query_range = element.find("query_range")
        if query_range is not None and query_range.text:
            evidence.query_range = query_range.text
            
        hit_range = element.find("hit_range")
        if hit_range is not None and hit_range.text:
            evidence.hit_range = hit_range.text
        
        # Get classification
        for cls_attr in ["t_group", "h_group", "x_group", "a_group"]:
            value = element.get(cls_attr)
            if value:
                setattr(evidence, cls_attr, value)
        
        # Get other attributes
        for key, value in element.attrib.items():
            if key not in ["type", "source_id", "confidence", "t_group", "h_group", "x_group", "a_group"]:
                evidence.attributes[key] = value
        
        return evidence

@dataclass
class BlastEvidence(Evidence):
    """Evidence from BLAST results"""
    evalue: float = 999.0
    identity: float = 0.0
    coverage: float = 0.0
    
    @classmethod
    def from_blast_hit(cls, hit, hit_type="domain_blast"):
        """Create from BlastHit model"""
        return cls(
            type=hit_type,
            source_id=hit.domain_id or hit.hit_id,
            query_range=hit.range,
            hit_range=hit.hit_range,
            confidence=1.0 / (1.0 + hit.evalue),  # Convert evalue to confidence
            evalue=hit.evalue,
            attributes={
                "hit_id": hit.hit_id,
                "pdb_id": hit.pdb_id,
                "chain_id": hit.chain_id
            }
        )

@dataclass
class HHSearchEvidence(Evidence):
    """Evidence from HHSearch results"""
    probability: float = 0.0
    evalue: float = 999.0
    score: float = 0.0
    
    @classmethod
    def from_hhsearch_hit(cls, hit):
        """Create from HHSearchHit model"""
        return cls(
            type="hhsearch",
            source_id=hit.domain_id or hit.hit_id,
            query_range=hit.range,
            hit_range=hit.hit_range,
            confidence=hit.probability / 100.0,  # Convert probability to confidence
            probability=hit.probability,
            evalue=hit.evalue,
            score=hit.score
        )
