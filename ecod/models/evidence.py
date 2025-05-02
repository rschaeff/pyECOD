# ecod/models/domain_analysis/evidence.py
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional
from ecod.models.base import XmlSerializable

# Enhance DomainEvidence model to support hierarchical evidence structure
@dataclass
class DomainEvidence(XmlSerializable):
    """Evidence supporting a domain boundary decision"""
    type: str  # "blast", "hhsearch", etc.
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

    # Type-specific details
    details: Dict[str, Any] = field(default_factory=dict)

    def to_xml(self) -> ET.Element:
        """Convert to XML Element with complete details"""
        element = ET.Element("evidence")
        element.set("type", self.type)
        element.set("source_id", self.source_id)

        if self.domain_id:
            element.set("domain_id", self.domain_id)

        element.set("confidence", f"{self.confidence:.4f}")

        # Add range elements
        query_range = ET.SubElement(element, "query_range")
        query_range.text = self.query_range

        hit_range = ET.SubElement(element, "hit_range")
        hit_range.text = self.hit_range

        # Add classification if available
        for cls_attr in ["t_group", "h_group", "x_group", "a_group"]:
            value = getattr(self, cls_attr)
            if value:
                element.set(cls_attr, value)

        # Add type-specific details
        if self.details:
            details_elem = ET.SubElement(element, "details")
            for key, value in self.details.items():
                if isinstance(value, (str, int, float, bool)):
                    details_elem.set(key, str(value))

        return element
