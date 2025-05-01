# ecod/models/domain_analysis/domain_candidate.py
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Set
import xml.etree.ElementTree as ET
from ecod.models.base import XmlSerializable
from ecod.models.domain_analysis.evidence import Evidence

@dataclass
class DomainCandidate(XmlSerializable):
    """Domain candidate model for boundary determination"""
    start: int
    end: int
    evidence: List[Evidence] = field(default_factory=list)
    confidence: float = 0.0  # Overall confidence score
    source: str = ""  # Primary evidence source
    id: Optional[str] = None  # Candidate ID (assigned during processing)
    
    # Classification info (if available)
    t_group: Optional[str] = None
    h_group: Optional[str] = None
    x_group: Optional[str] = None
    a_group: Optional[str] = None
    
    # Processing metadata
    protected: bool = False  # Flag for high-confidence candidates
    
    @property
    def size(self) -> int:
        """Get domain size"""
        return self.end - self.start + 1
    
    @property
    def range(self) -> str:
        """Get domain range"""
        return f"{self.start}-{self.end}"
    
    def get_positions(self) -> Set[int]:
        """Get set of positions covered by this domain"""
        return set(range(self.start, self.end + 1))
    
    def overlaps(self, other: 'DomainCandidate') -> bool:
        """Check if this candidate overlaps with another"""
        return max(self.start, other.start) <= min(self.end, other.end)
    
    def overlap_size(self, other: 'DomainCandidate') -> int:
        """Calculate overlap size with another candidate"""
        if not self.overlaps(other):
            return 0
        return min(self.end, other.end) - max(self.start, other.start) + 1
    
    def overlap_percentage(self, other: 'DomainCandidate') -> float:
        """Calculate overlap percentage with another candidate"""
        overlap = self.overlap_size(other)
        min_size = min(self.size, other.size)
        return (overlap / min_size) * 100.0 if min_size > 0 else 0.0
    
    def add_evidence(self, evidence: Evidence) -> None:
        """Add evidence to candidate"""
        self.evidence.append(evidence)
        
        # Update confidence based on new evidence
        self._update_confidence()
        
        # Update classification if not already set
        self._update_classification(evidence)
    
    def _update_confidence(self) -> None:
        """Update overall confidence based on evidence"""
        if not self.evidence:
            self.confidence = 0.0
            return
            
        # Calculate weighted average of evidence confidences
        total_weight = 0.0
        weighted_sum = 0.0
        
        for e in self.evidence:
            # Determine weight based on evidence type
            weight = 1.0
            if e.type == "domain_blast":
                weight = 3.0  # Prioritize domain BLAST
            elif e.type == "hhsearch":
                weight = 2.0  # Next priority to HHSearch
                
            weighted_sum += e.confidence * weight
            total_weight += weight
        
        self.confidence = weighted_sum / total_weight if total_weight > 0 else 0.0
    
    def _update_classification(self, evidence: Evidence) -> None:
        """Update classification from evidence if not already set"""
        for cls_attr in ["t_group", "h_group", "x_group", "a_group"]:
            if getattr(self, cls_attr) is None:
                evidence_value = getattr(evidence, cls_attr)
                if evidence_value:
                    setattr(self, cls_attr, evidence_value)
    
    def to_xml(self) -> ET.Element:
        """Convert to XML Element"""
        element = ET.Element("domain_candidate")
        element.set("start", str(self.start))
        element.set("end", str(self.end))
        element.set("size", str(self.size))
        element.set("confidence", f"{self.confidence:.4f}")
        
        if self.id:
            element.set("id", self.id)
        
        if self.source:
            element.set("source", self.source)
        
        if self.protected:
            element.set("protected", "true")
        
        # Add classification if available
        for cls_attr in ["t_group", "h_group", "x_group", "a_group"]:
            value = getattr(self, cls_attr)
            if value:
                element.set(cls_attr, value)
        
        # Add evidence
        if self.evidence:
            evidence_list = ET.SubElement(element, "evidence_list")
            for e in self.evidence:
                evidence_list.append(e.to_xml())
        
        return element
    
    @classmethod
    def from_xml(cls, element: ET.Element) -> 'DomainCandidate':
        """Create from XML Element"""
        from ecod.models.domain_analysis.evidence import Evidence
        
        candidate = cls(
            start=int(element.get("start", "0")),
            end=int(element.get("end", "0")),
            confidence=float(element.get("confidence", "0.0")),
            source=element.get("source", ""),
            id=element.get("id"),
            protected=element.get("protected", "false").lower() == "true"
        )
        
        # Get classification
        for cls_attr in ["t_group", "h_group", "x_group", "a_group"]:
            value = element.get(cls_attr)
            if value:
                setattr(candidate, cls_attr, value)
        
        # Get evidence
        evidence_list = element.find("evidence_list")
        if evidence_list is not None:
            for e_elem in evidence_list.findall("evidence"):
                evidence = Evidence.from_xml(e_elem)
                candidate.evidence.append(evidence)
        
        return candidate
