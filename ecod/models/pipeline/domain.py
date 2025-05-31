# ecod/models/pipeline/domain.py
"""
Consolidated DomainModel for pyECOD

This module consolidates the DomainModel classes from:
- models/domain.py -> DomainModel (analysis-focused)
- models/domain_analysis/domain_model.py -> DomainModel (workflow-focused)

Into a single, comprehensive domain model for pipeline processing.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Set, Union
import xml.etree.ElementTree as ET
from ecod.models.base import XmlSerializable
import math


@dataclass
class DomainModel(XmlSerializable):
    """
    Consolidated domain model for pipeline processing.
    
    Merges functionality from both existing DomainModel implementations
    to provide a single, comprehensive model for domain analysis workflows.
    """
    
    # Core identity fields
    id: str
    start: int
    end: int
    range: str
    
    # Classification hierarchy (ECOD classification)
    t_group: Optional[str] = None
    h_group: Optional[str] = None
    x_group: Optional[str] = None
    a_group: Optional[str] = None
    
    # Domain properties and metadata
    source: str = ""  # Source of domain boundary ("hhsearch", "blast", etc.)
    confidence: float = 0.0  # Overall confidence score (0-1)
    source_id: str = ""  # ID of the source evidence (domain_id, hit_id, etc.)
    
    # ECOD representative and filter flags
    is_manual_rep: bool = False
    is_f70: bool = False
    is_f40: bool = False
    is_f99: bool = False
    
    # Evidence supporting this domain boundary and classification
    evidence: List[Union['Evidence', Dict[str, Any]]] = field(default_factory=list)
    
    # Additional metadata for analysis workflows
    protected: bool = False  # Flag for high-confidence domains that shouldn't be merged
    quality_score: Optional[float] = None  # Additional quality metric
    notes: Optional[str] = None  # Free-text notes about domain determination
    
    def __post_init__(self):
        """Post-initialization processing"""
        # Ensure evidence uses the new Evidence model when available
        self._standardize_evidence()

        # CRITICAL FIX: Apply classifications from evidence AFTER standardization
        for evidence in self.evidence:
            self._update_classification_from_evidence(evidence)

        # CRITICAL FIX: Calculate confidence AFTER evidence processing
        if self.confidence == 0.0 and self.evidence:
            self.confidence = self._calculate_confidence()

        # Set protected flag for very high confidence domains
        if self.confidence >= 0.98:
            self.protected = True
    
    def _standardize_evidence(self):
        """Convert any dictionary evidence to Evidence objects when possible"""

        standardized_evidence = []
        for ev in self.evidence:
            if isinstance(ev, dict):
                # CRITICAL FIX: Only skip dictionaries that are completely invalid
                if 'type' not in ev or ev.get('type') is None or ev.get('type') == '':
                    continue  # Skip malformed evidence (no type)

                try:
                    # Try to convert to Evidence object
                    evidence_obj = Evidence.from_dict(ev)
                    standardized_evidence.append(evidence_obj)
                except Exception:
                    # CRITICAL FIX: If conversion fails, keep the original dictionary
                    # This preserves valid evidence that just can't be converted
                    standardized_evidence.append(ev)
            else:
                # Keep Evidence objects and other types as-is
                standardized_evidence.append(ev)
        
        self.evidence = standardized_evidence
    
    def _calculate_confidence(self) -> float:
        """
        FIXED: Calculate confidence score from evidence with robust edge case handling
        """
        if not self.evidence:
            return 0.0

        # Evidence type weights
        weights = {
            "domain_blast": 3.0,
            "hhsearch": 1.0,
            "chain_blast": 5.0
        }

        total_weight = 0.0
        weighted_sum = 0.0

        for ev in self.evidence:
            # Extract evidence type and confidence with robust handling
            if hasattr(ev, 'type'):
                ev_type = ev.type
                ev_confidence = getattr(ev, 'confidence', 0.0)
            elif isinstance(ev, dict):
                ev_type = ev.get('type')  # FIXED: Don't default to 'unknown'
                ev_confidence = ev.get('confidence', 0.0)
            else:
                continue  # Skip invalid evidence

            # CRITICAL FIX 1: Skip evidence with invalid types or confidence
            if ev_type is None or ev_type == '':
                continue

            if ev_confidence is None:
                continue

            # CRITICAL FIX 2: Handle NaN and infinity values
            if not isinstance(ev_confidence, (int, float)):
                continue

            if math.isnan(ev_confidence) or math.isinf(ev_confidence):
                continue

            # CRITICAL FIX 3: Skip evidence with invalid confidence ranges
            if ev_confidence < 0.0 or ev_confidence > 1.0:
                continue

            weight = weights.get(ev_type, 1.0)
            weighted_sum += ev_confidence * weight
            total_weight += weight

        if total_weight == 0:
            return 0.0

        confidence = weighted_sum / total_weight

        # CRITICAL FIX 4: Ensure final confidence is always in valid range
        return max(0.0, min(1.0, confidence))
    
    @property
    def size(self) -> int:
        """Get domain size in residues"""
        return self.end - self.start + 1
    
    @property
    def length(self) -> int:
        """Alias for size for backward compatibility"""
        return self.size
    
    def get_positions(self) -> Set[int]:
        """Get set of positions covered by this domain"""
        return set(range(self.start, self.end + 1))
    
    def get_range_segments(self) -> List[tuple]:
        """Parse range string into list of (start, end) tuples"""
        segments = []
        if not self.range:
            return [(self.start, self.end)]
        
        for segment in self.range.split(','):
            if '-' in segment:
                try:
                    start, end = map(int, segment.split('-'))
                    segments.append((start, end))
                except ValueError:
                    pass
        
        return segments if segments else [(self.start, self.end)]
    
    def overlaps(self, other: 'DomainModel') -> bool:
        """Check if this domain overlaps with another"""
        return max(self.start, other.start) <= min(self.end, other.end)
    
    def overlap_size(self, other: 'DomainModel') -> int:
        """Calculate overlap size with another domain"""
        if not self.overlaps(other):
            return 0
        return min(self.end, other.end) - max(self.start, other.start) + 1
    
    def overlap_percentage(self, other: 'DomainModel') -> float:
        """Calculate overlap percentage with another domain"""
        overlap = self.overlap_size(other)
        min_size = min(self.size, other.size)
        return (overlap / min_size) * 100.0 if min_size > 0 else 0.0

    def recalculate_confidence(self) -> None:
        """
        FIXED: Force recalculation of domain confidence and update self.confidence
        """
        self.confidence = self._calculate_confidence()
    
    def add_evidence(self, evidence: Union['Evidence', Dict[str, Any]]) -> None:
        """
        FIXED: Add evidence to this domain and update confidence
        """
        self.evidence.append(evidence)

        # Update classification if not already set
        self._update_classification_from_evidence(evidence)

        # CRITICAL FIX: Recalculate confidence after adding evidence
        self.confidence = self._calculate_confidence()
    
    def _update_classification_from_evidence(self, evidence: Union['Evidence', Dict[str, Any]]) -> None:
        """
        FIXED: Update classification from evidence if not already set
        """
        # Extract classification from evidence
        for cls_attr in ["t_group", "h_group", "x_group", "a_group"]:
            # Only update if domain doesn't already have this classification
            if getattr(self, cls_attr) is None:
                # Get value from evidence
                evidence_value = None

                if hasattr(evidence, cls_attr):
                    evidence_value = getattr(evidence, cls_attr)
                elif isinstance(evidence, dict):
                    evidence_value = evidence.get(cls_attr)

                # CRITICAL FIX: Set classification if evidence has a value (including empty strings)
                # Use 'is not None' instead of truthiness check to preserve empty strings
                if evidence_value is not None:
                    setattr(self, cls_attr, evidence_value)
    
    def is_classified(self) -> bool:
        """Check if domain has any classification"""
        return any([self.t_group, self.h_group, self.x_group, self.a_group])
    
    def is_fully_classified(self) -> bool:
        """Check if domain has complete ECOD classification"""
        return all([self.t_group, self.h_group, self.x_group, self.a_group])
    
    def get_classification_level(self) -> str:
        """Get the deepest classification level available"""
        if self.a_group:
            return "A-group"
        elif self.x_group:
            return "X-group"
        elif self.h_group:
            return "H-group"
        elif self.t_group:
            return "T-group"
        else:
            return "Unclassified"
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation"""
        result = {
            "id": self.id,
            "start": self.start,
            "end": self.end,
            "range": self.range,
            "size": self.size,
            "source": self.source,
            "confidence": self.confidence,
            "source_id": self.source_id,
            "t_group": self.t_group,
            "h_group": self.h_group,
            "x_group": self.x_group,
            "a_group": self.a_group,
            "is_manual_rep": self.is_manual_rep,
            "is_f70": self.is_f70,
            "is_f40": self.is_f40,
            "is_f99": self.is_f99,
            "protected": self.protected
        }
        
        # Add evidence as dictionaries
        if self.evidence:
            evidence_dicts = []
            for ev in self.evidence:
                if hasattr(ev, 'to_dict'):
                    evidence_dicts.append(ev.to_dict())
                elif isinstance(ev, dict):
                    evidence_dicts.append(ev)
                else:
                    # Try to convert object to dict
                    evidence_dicts.append(vars(ev) if hasattr(ev, '__dict__') else str(ev))
            result["evidence"] = evidence_dicts
        
        # Add optional fields if they have values
        if self.quality_score is not None:
            result["quality_score"] = self.quality_score
        if self.notes:
            result["notes"] = self.notes
        
        return result
    
    @classmethod
    def from_dict(cls, domain_dict: Dict[str, Any], domain_id: Optional[str] = None) -> 'DomainModel':
        """Create DomainModel from dictionary"""
        # Extract core fields
        domain_id = domain_id or domain_dict.get("id", domain_dict.get("domain_id", "unknown"))
        start = domain_dict.get("start", 0)
        end = domain_dict.get("end", 0)
        range_str = domain_dict.get("range", f"{start}-{end}")
        
        # Create domain
        domain = cls(
            id=domain_id,
            start=start,
            end=end,
            range=range_str,
            t_group=domain_dict.get("t_group"),
            h_group=domain_dict.get("h_group"),
            x_group=domain_dict.get("x_group"),
            a_group=domain_dict.get("a_group"),
            source=domain_dict.get("source", ""),
            confidence=domain_dict.get("confidence", 0.0),
            source_id=domain_dict.get("source_id", ""),
            is_manual_rep=domain_dict.get("is_manual_rep", False),
            is_f70=domain_dict.get("is_f70", False),
            is_f40=domain_dict.get("is_f40", False),
            is_f99=domain_dict.get("is_f99", False),
            protected=domain_dict.get("protected", False),
            quality_score=domain_dict.get("quality_score"),
            notes=domain_dict.get("notes")
        )
        
        # Add evidence if present
        if "evidence" in domain_dict and domain_dict["evidence"]:
            for ev_item in domain_dict["evidence"]:
                domain.add_evidence(ev_item)
        
        return domain
    
    def to_xml(self) -> ET.Element:
        """Convert to XML Element with comprehensive information"""
        element = ET.Element("domain")
        
        # Basic attributes
        element.set("id", self.id)
        element.set("start", str(self.start))
        element.set("end", str(self.end))
        element.set("range", self.range)
        element.set("size", str(self.size))
        
        # Source and confidence
        if self.source:
            element.set("source", self.source)
        element.set("confidence", f"{self.confidence:.4f}")
        if self.source_id:
            element.set("source_id", self.source_id)
        
        # Classification
        for cls_attr in ["t_group", "h_group", "x_group", "a_group"]:
            value = getattr(self, cls_attr)
            if value:
                element.set(cls_attr, value)
        
        # Representative/Filter flags
        element.set("is_manual_rep", str(self.is_manual_rep).lower())
        element.set("is_f70", str(self.is_f70).lower())
        element.set("is_f40", str(self.is_f40).lower())
        element.set("is_f99", str(self.is_f99).lower())
        
        # Additional metadata
        if self.protected:
            element.set("protected", "true")
        if self.quality_score is not None:
            element.set("quality_score", f"{self.quality_score:.4f}")
        if self.notes:
            element.set("notes", self.notes)
        
        # Add evidence if available
        if self.evidence:
            evidence_elem = ET.SubElement(element, "evidence_list")
            evidence_elem.set("count", str(len(self.evidence)))
            
            for i, evidence in enumerate(self.evidence):
                if hasattr(evidence, 'to_xml'):
                    # Use Evidence object's XML method
                    evidence_elem.append(evidence.to_xml())
                else:
                    # Create XML from dictionary or other object
                    ev_elem = ET.SubElement(evidence_elem, "evidence")
                    ev_elem.set("index", str(i))
                    
                    if isinstance(evidence, dict):
                        for key, value in evidence.items():
                            if value is not None:
                                if isinstance(value, bool):
                                    ev_elem.set(key, str(value).lower())
                                else:
                                    ev_elem.set(key, str(value))
                    else:
                        # Convert object to attributes
                        for attr in dir(evidence):
                            if not attr.startswith('_'):
                                value = getattr(evidence, attr)
                                if value is not None and not callable(value):
                                    if isinstance(value, bool):
                                        ev_elem.set(attr, str(value).lower())
                                    else:
                                        ev_elem.set(attr, str(value))
        
        return element
    
    @classmethod
    def from_xml(cls, element: ET.Element) -> 'DomainModel':
        """Create from XML Element"""
        # Extract basic attributes
        domain_id = element.get("id", "unknown")
        start = int(element.get("start", "0"))
        end = int(element.get("end", "0"))
        range_str = element.get("range", f"{start}-{end}")
        
        # Create domain
        domain = cls(
            id=domain_id,
            start=start,
            end=end,
            range=range_str,
            source=element.get("source", ""),
            confidence=float(element.get("confidence", "0.0")),
            source_id=element.get("source_id", ""),
            t_group=element.get("t_group"),
            h_group=element.get("h_group"),
            x_group=element.get("x_group"),
            a_group=element.get("a_group"),
            is_manual_rep=element.get("is_manual_rep", "false").lower() == "true",
            is_f70=element.get("is_f70", "false").lower() == "true",
            is_f40=element.get("is_f40", "false").lower() == "true",
            is_f99=element.get("is_f99", "false").lower() == "true",
            protected=element.get("protected", "false").lower() == "true",
            quality_score=float(element.get("quality_score")) if element.get("quality_score") else None,
            notes=element.get("notes")
        )
        
        # Extract evidence if available
        evidence_list = element.find("evidence_list")
        if evidence_list is not None:
            try:
                from ecod.models.pipeline.evidence import Evidence
                NEW_EVIDENCE_MODEL = True
            except ImportError:
                NEW_EVIDENCE_MODEL = False
            
            for ev_elem in evidence_list.findall("evidence"):
                if NEW_EVIDENCE_MODEL:
                    try:
                        evidence = Evidence.from_xml(ev_elem)
                        domain.evidence.append(evidence)
                    except Exception:
                        # Fall back to dictionary
                        evidence_dict = dict(ev_elem.attrib)
                        domain.evidence.append(evidence_dict)
                else:
                    # Use dictionary approach
                    evidence_dict = dict(ev_elem.attrib)
                    domain.evidence.append(evidence_dict)
        
        return domain
    
    def merge_with(self, other: 'DomainModel') -> 'DomainModel':
        """
        FIXED: Merge this domain with another, keeping the best evidence
        """
        # Determine which domain has higher confidence
        if self.confidence >= other.confidence:
            primary, secondary = self, other
        else:
            primary, secondary = other, self

        # Create merged domain
        merged = DomainModel(
            id=f"{primary.id}_merged",
            start=min(self.start, other.start),
            end=max(self.end, other.end),
            range=f"{min(self.start, other.start)}-{max(self.end, other.end)}",
            source=primary.source,
            confidence=max(self.confidence, other.confidence),
            source_id=primary.source_id,

            # CRITICAL FIX: Properly merge classifications (primary wins, secondary fills gaps)
            t_group=primary.t_group or secondary.t_group,
            h_group=primary.h_group or secondary.h_group,
            x_group=primary.x_group or secondary.x_group,  # FIXED: was using a_group incorrectly
            a_group=primary.a_group or secondary.a_group,

            is_manual_rep=primary.is_manual_rep or secondary.is_manual_rep,
            is_f70=primary.is_f70 or secondary.is_f70,
            is_f40=primary.is_f40 or secondary.is_f40,
            is_f99=primary.is_f99 or secondary.is_f99,
            protected=primary.protected or secondary.protected,
            quality_score=primary.quality_score or secondary.quality_score,
            notes=f"Merged from {self.id} and {other.id}"
        )

        # Merge evidence
        merged.evidence = primary.evidence + secondary.evidence
        
        return merged
    
    def split_at(self, position: int) -> tuple['DomainModel', 'DomainModel']:
        """Split domain at specified position"""
        if position <= self.start or position >= self.end:
            raise ValueError(f"Split position {position} not within domain range {self.start}-{self.end}")
        
        # Create first domain
        domain1 = DomainModel(
            id=f"{self.id}_part1",
            start=self.start,
            end=position - 1,
            range=f"{self.start}-{position - 1}",
            source=self.source,
            confidence=self.confidence * 0.9,  # Reduce confidence for split domains
            source_id=self.source_id,
            t_group=self.t_group,
            h_group=self.h_group,
            x_group=self.x_group,
            a_group=self.a_group,
            is_manual_rep=self.is_manual_rep,
            is_f70=self.is_f70,
            is_f40=self.is_f40,
            is_f99=self.is_f99,
            evidence=self.evidence.copy(),
            notes=f"Split from {self.id} at position {position}"
        )
        
        # Create second domain
        domain2 = DomainModel(
            id=f"{self.id}_part2",
            start=position,
            end=self.end,
            range=f"{position}-{self.end}",
            source=self.source,
            confidence=self.confidence * 0.9,
            source_id=self.source_id,
            t_group=self.t_group,
            h_group=self.h_group,
            x_group=self.x_group,
            a_group=self.a_group,
            is_manual_rep=self.is_manual_rep,
            is_f70=self.is_f70,
            is_f40=self.is_f40,
            is_f99=self.is_f99,
            evidence=self.evidence.copy(),
            notes=f"Split from {self.id} at position {position}"
        )
        
        return domain1, domain2
    
    def __str__(self) -> str:
        """String representation"""
        return f"Domain({self.id}, {self.range}, {self.get_classification_level()})"
    
    def __repr__(self) -> str:
        """Detailed string representation"""
        return (f"DomainModel(id='{self.id}', range='{self.range}', "
                f"source='{self.source}', confidence={self.confidence:.3f}, "
                f"classification='{self.get_classification_level()}')")


# Backward compatibility aliases
Domain = DomainModel  # For any code expecting 'Domain' class
