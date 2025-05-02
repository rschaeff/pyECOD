# ecod/models/domain_analysis/domain_model.py
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Set
import xml.etree.ElementTree as ET
from ecod.models.base import XmlSerializable
from ecod.models.domain_analysis.evidence import Evidence

@dataclass
class DomainModel(XmlSerializable):
    """Model for a single domain within a protein chain"""
    id: str
    start: int
    end: int
    range: str
    
    # Classification
    t_group: Optional[str] = None
    h_group: Optional[str] = None
    x_group: Optional[str] = None
    a_group: Optional[str] = None
    
    # Domain properties
    source: str = ""
    confidence: float = 0.0
    source_id: str = ""
    is_manual_rep: bool = False
    is_f70: bool = False
    is_f40: bool = False
    is_f99: bool = False
    
    # Evidence
    evidence: List[Evidence] = field(default_factory=list)
    
    @property
    def size(self) -> int:
        """Get domain size"""
        return self.end - self.start + 1
    
    def get_positions(self) -> Set[int]:
        """Get set of positions covered by this domain"""
        return set(range(self.start, self.end + 1))
    
    def to_xml(self) -> ET.Element:
        """Convert to XML Element"""
        element = ET.Element("domain")
        
        # Basic attributes
        element.set("id", self.id)
        element.set("start", str(self.start))
        element.set("end", str(self.end))
        element.set("range", self.range)
        
        # Classification
        for cls_attr in ["t_group", "h_group", "x_group", "a_group"]:
            value = getattr(self, cls_attr)
            if value:
                element.set(cls_attr, value)
        
        # Domain properties
        if self.source:
            element.set("source", self.source)
        element.set("confidence", f"{self.confidence:.4f}")
        if self.source_id:
            element.set("source_id", self.source_id)
        
        # Representative/Filter flags
        element.set("is_manual_rep", str(self.is_manual_rep).lower())
        element.set("is_f70", str(self.is_f70).lower())
        element.set("is_f40", str(self.is_f40).lower()) 
        element.set("is_f99", str(self.is_f99).lower())
        
        # Add evidence if available
        if self.evidence:
            evidence_elem = ET.SubElement(element, "evidence_list")
            for evidence in self.evidence:
                evidence_elem.append(evidence.to_xml())
        
        return element
    
    @classmethod
    def from_xml(cls, element: ET.Element) -> 'DomainModel':
        """Create from XML Element"""
        from ecod.models.domain_analysis.evidence import Evidence
        
        domain = cls(
            id=element.get("id", ""),
            start=int(element.get("start", "0")),
            end=int(element.get("end", "0")),
            range=element.get("range", ""),
            t_group=element.get("t_group"),
            h_group=element.get("h_group"),
            x_group=element.get("x_group"),
            a_group=element.get("a_group"),
            source=element.get("source", ""),
            confidence=float(element.get("confidence", "0.0")),
            source_id=element.get("source_id", ""),
            is_manual_rep=element.get("is_manual_rep", "false").lower() == "true",
            is_f70=element.get("is_f70", "false").lower() == "true",
            is_f40=element.get("is_f40", "false").lower() == "true",
            is_f99=element.get("is_f99", "false").lower() == "true"
        )
        
        # Get evidence if available
        evidence_list = element.find("evidence_list")
        if evidence_list is not None:
            for ev_elem in evidence_list.findall("evidence"):
                evidence = Evidence.from_xml(ev_elem)
