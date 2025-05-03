"""
Enhanced DomainEvidence Class for pyECOD

This module implements a comprehensive DomainEvidence class that bridges
between dictionary-based evidence representation and model-based approaches.
The class follows the design outlined in the Domain Model Serialization Strategy
document and supports both backward and forward compatibility.
"""

from dataclasses import dataclass, field
from typing import Dict, Any, Optional, ClassVar, List, Union
import xml.etree.ElementTree as ET
from ecod.models.base import XmlSerializable


@dataclass
class DomainEvidence(XmlSerializable):
    """Evidence supporting a domain boundary decision"""
    type: str  # "blast", "hhsearch", "domain_blast", "chain_blast", etc.
    source_id: str = ""  # Domain ID, hit ID, etc.
    domain_id: str = ""  # Domain ID, if available
    query_range: str = ""  # Query region in range format
    hit_range: str = ""  # Hit region in range format
    confidence: float = 0.0  # Normalized confidence score (0-1)

    # Classification data (if available)
    t_group: Optional[str] = None
    h_group: Optional[str] = None
    x_group: Optional[str] = None
    a_group: Optional[str] = None

    # Source-specific attributes
    attributes: Dict[str, Any] = field(default_factory=dict)

    # Element path for XML parsing
    xml_element_path: ClassVar[str] = ".//evidence"

    def __post_init__(self):
        """Normalize confidence value to 0-1 range"""
        # Handle probability to confidence conversion for HHSearch
        if self.type == "hhsearch" and "probability" in self.attributes:
            prob = self.attributes.get("probability", 0.0)
            if isinstance(prob, (int, float)) and prob > 1.0:
                # HHSearch probabilities are 0-100, normalize to 0-1
                self.confidence = prob / 100.0

        # Handle e-value to confidence conversion for BLAST
        elif self.type in ["blast", "domain_blast", "chain_blast"] and "evalue" in self.attributes:
            evalue = self.attributes.get("evalue", 999.0)
            if isinstance(evalue, (int, float)):
                # Convert e-value to confidence score (higher e-value = lower confidence)
                self.confidence = 1.0 / (1.0 + evalue)

    @classmethod
    def from_dict(cls, evidence_dict: Dict[str, Any]) -> 'DomainEvidence':
        """Create evidence from dictionary representation

        This method handles conversion from the traditional dictionary
        representation used in the codebase to the structured model.

        Args:
            evidence_dict: Dictionary with evidence data

        Returns:
            DomainEvidence instance
        """
        # Extract common fields
        evidence_type = evidence_dict.get("type", "unknown")
        source_id = evidence_dict.get("source_id", "")
        domain_id = evidence_dict.get("domain_id", "")

        # Get range information - handle different field names
        query_range = evidence_dict.get("query_range", "")
        if not query_range:
            # Try alternate field names
            query_range = evidence_dict.get("range", "")
            if not query_range:
                query_range = evidence_dict.get("query_regions", "")

        hit_range = evidence_dict.get("hit_range", "")

        # Get confidence or calculate from other fields
        confidence = evidence_dict.get("confidence", 0.0)

        # Get classification data
        t_group = evidence_dict.get("t_group")
        h_group = evidence_dict.get("h_group")
        x_group = evidence_dict.get("x_group")
        a_group = evidence_dict.get("a_group")

        # Create attributes dict with remaining fields
        attributes = {}
        for key, value in evidence_dict.items():
            if key not in ["type", "source_id", "domain_id", "query_range", "hit_range",
                          "confidence", "t_group", "h_group", "x_group", "a_group",
                          "range", "query_regions"]:
                attributes[key] = value

        # Create instance
        return cls(
            type=evidence_type,
            source_id=source_id,
            domain_id=domain_id,
            query_range=query_range,
            hit_range=hit_range,
            confidence=confidence,
            t_group=t_group,
            h_group=h_group,
            x_group=x_group,
            a_group=a_group,
            attributes=attributes
        )

    @classmethod
    def from_blast_hit(cls, hit, hit_type="domain_blast") -> 'DomainEvidence':
        """Create evidence from a BlastHit model

        Args:
            hit: BlastHit instance
            hit_type: Type of BLAST hit (domain_blast, chain_blast)

        Returns:
            DomainEvidence instance
        """
        # Extract common fields
        source_id = getattr(hit, "domain_id", "") or getattr(hit, "hit_id", "")
        query_range = getattr(hit, "range", "")
        hit_range = getattr(hit, "hit_range", "")

        # Get e-value for confidence calculation
        evalue = getattr(hit, "evalue", 999.0)

        # Create attributes with remaining fields
        attributes = {
            "hit_id": getattr(hit, "hit_id", ""),
            "pdb_id": getattr(hit, "pdb_id", ""),
            "chain_id": getattr(hit, "chain_id", ""),
            "evalue": evalue,
            "hsp_count": getattr(hit, "hsp_count", 0),
        }

        # Add additional attributes if available
        if hasattr(hit, "discontinuous"):
            attributes["discontinuous"] = hit.discontinuous

        # Create instance
        return cls(
            type=hit_type,
            source_id=source_id,
            domain_id=getattr(hit, "domain_id", ""),
            query_range=query_range,
            hit_range=hit_range,
            confidence=1.0 / (1.0 + evalue),  # Convert e-value to confidence
            attributes=attributes
        )

    @classmethod
    def from_hhsearch_hit(cls, hit) -> 'DomainEvidence':
        """Create evidence from an HHSearchHit model

        Args:
            hit: HHSearchHit instance

        Returns:
            DomainEvidence instance
        """
        # Extract common fields
        source_id = getattr(hit, "domain_id", "") or getattr(hit, "hit_id", "")
        query_range = getattr(hit, "range", "")
        hit_range = getattr(hit, "hit_range", "")

        # Get probability for confidence calculation
        probability = getattr(hit, "probability", 0.0)

        # Create attributes with remaining fields
        attributes = {
            "hit_id": getattr(hit, "hit_id", ""),
            "probability": probability,
            "evalue": getattr(hit, "evalue", 999.0),
            "score": getattr(hit, "score", 0.0)
        }

        # Create instance
        return cls(
            type="hhsearch",
            source_id=source_id,
            domain_id=getattr(hit, "domain_id", ""),
            query_range=query_range,
            hit_range=hit_range,
            confidence=probability / 100.0,  # Convert probability to confidence
            attributes=attributes
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation

        This method provides backward compatibility with code that expects
        dictionary-based evidence.

        Returns:
            Dictionary representation
        """
        # Start with base fields
        result = {
            "type": self.type,
            "source_id": self.source_id,
            "domain_id": self.domain_id,
            "query_range": self.query_range,
            "hit_range": self.hit_range,
            "confidence": self.confidence
        }

        # Add classification fields if available
        for cls_attr in ["t_group", "h_group", "x_group", "a_group"]:
            value = getattr(self, cls_attr)
            if value:
                result[cls_attr] = value

        # Add all attributes
        for key, value in self.attributes.items():
            result[key] = value

        # Add backward compatibility aliases (for code that expects specific field names)
        if "query_range" in result and self.query_range:
            result["range"] = self.query_range

        return result

    def to_xml(self) -> ET.Element:
        """Convert to XML Element

        Returns:
            XML Element representing the evidence
        """
        element = ET.Element("evidence")
        element.set("type", self.type)
        element.set("source_id", self.source_id)

        if self.domain_id:
            element.set("domain_id", self.domain_id)

        # Add confidence
        element.set("confidence", f"{self.confidence:.4f}")

        # Add ranges as child elements
        if self.query_range:
            query_range = ET.SubElement(element, "query_range")
            query_range.text = self.query_range

        if self.hit_range:
            hit_range = ET.SubElement(element, "hit_range")
            hit_range.text = self.hit_range

        # Add classification if available
        for cls_attr in ["t_group", "h_group", "x_group", "a_group"]:
            value = getattr(self, cls_attr)
            if value:
                element.set(cls_attr, value)

        # Add source-specific attributes
        for key, value in self.attributes.items():
            if value is not None:
                if isinstance(value, bool):
                    element.set(key, str(value).lower())
                else:
                    element.set(key, str(value))

        return element

    @classmethod
    def from_xml(cls, element: ET.Element) -> 'DomainEvidence':
        """Create from XML Element

        Args:
            element: XML Element

        Returns:
            DomainEvidence instance
        """
        # Extract basic attributes
        evidence_type = element.get("type", "unknown")
        source_id = element.get("source_id", "")
        domain_id = element.get("domain_id", "")
        confidence = float(element.get("confidence", "0.0"))

        # Extract range elements
        query_range = ""
        hit_range = ""

        query_range_elem = element.find("query_range")
        if query_range_elem is not None and query_range_elem.text:
            query_range = query_range_elem.text.strip()

        hit_range_elem = element.find("hit_range")
        if hit_range_elem is not None and hit_range_elem.text:
            hit_range = hit_range_elem.text.strip()

        # Extract classification
        t_group = element.get("t_group")
        h_group = element.get("h_group")
        x_group = element.get("x_group")
        a_group = element.get("a_group")

        # Create attributes dict with remaining fields
        attributes = {}
        for key, value in element.attrib.items():
            if key not in ["type", "source_id", "domain_id", "confidence",
                          "t_group", "h_group", "x_group", "a_group"]:
                # Convert numeric values
                if value.replace('.', '', 1).isdigit():
                    try:
                        if '.' in value:
                            attributes[key] = float(value)
                        else:
                            attributes[key] = int(value)
                    except ValueError:
                        attributes[key] = value
                # Convert boolean values
                elif value.lower() in ["true", "false"]:
                    attributes[key] = value.lower() == "true"
                else:
                    attributes[key] = value

        # Create instance
        return cls(
            type=evidence_type,
            source_id=source_id,
            domain_id=domain_id,
            query_range=query_range,
            hit_range=hit_range,
            confidence=confidence,
            t_group=t_group,
            h_group=h_group,
            x_group=x_group,
            a_group=a_group,
            attributes=attributes
        )

    @classmethod
    def create_from_multiple(cls, items: List[Dict[str, Any]],
                          consolidate: bool = True) -> List['DomainEvidence']:
        """Create multiple evidence items from a list of dictionaries

        Args:
            items: List of evidence dictionaries
            consolidate: Whether to consolidate duplicates

        Returns:
            List of DomainEvidence instances
        """
        evidence_list = []
        for item in items:
            evidence = cls.from_dict(item)
            evidence_list.append(evidence)

        if consolidate:
            return cls.consolidate_duplicates(evidence_list)

        return evidence_list

    @staticmethod
    def consolidate_duplicates(evidence_list: List['DomainEvidence']) -> List['DomainEvidence']:
        """Consolidate duplicate evidence items

        Args:
            evidence_list: List of evidence items

        Returns:
            Consolidated list with duplicates removed
        """
        # Create a dictionary to track unique evidence
        unique_evidence = {}

        for evidence in evidence_list:
            # Create a key based on source and range
            key = f"{evidence.type}:{evidence.source_id}:{evidence.query_range}"

            # If we already have this evidence, keep the one with higher confidence
            if key in unique_evidence:
                existing = unique_evidence[key]
                if evidence.confidence > existing.confidence:
                    unique_evidence[key] = evidence
            else:
                unique_evidence[key] = evidence

        # Return consolidated list
        return list(unique_evidence.values())
