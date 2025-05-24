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
    confidence: float = None
    
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

    def _calculate_confidence(self) -> float:
        """
        Calculate confidence score from available metrics

        Returns:
            float: Confidence score between 0.0 and 1.0

        Confidence calculation strategy by evidence type:
        - HHSearch: Primary use probability, fallback to evalue/score
        - BLAST: Primary use evalue, consider identity/coverage
        - Chain BLAST: Use evalue with domain mapping penalty
        """
        logger = logging.getLogger(__name__)

        try:
            if self.type == "hhsearch":
                return self._calculate_hhsearch_confidence()
            elif self.type in ("domain_blast", "blast"):
                return self._calculate_blast_confidence()
            elif self.type == "chain_blast":
                return self._calculate_chain_blast_confidence()
            else:
                # Generic calculation for unknown types
                return self._calculate_generic_confidence()

        except Exception as e:
            logger.warning(f"Error calculating confidence for {self.type} evidence: {e}")
            return 0.0

    def _calculate_hhsearch_confidence(self) -> float:
        """Calculate confidence for HHSearch evidence"""
        confidence = 0.0

        # Primary: Use probability (handle both 0-1 and 0-100 scales)
        if self.probability is not None:
            if self.probability > 1.0:
                # Assume 0-100 scale (typical HHSearch format)
                confidence = min(self.probability / 100.0, 1.0)
            else:
                # Assume 0-1 scale
                confidence = min(self.probability, 1.0)

            # HHSearch probabilities above 90% are very reliable
            if confidence >= 0.9:
                confidence = min(confidence * 1.05, 1.0)  # Small boost for high confidence

        # Fallback: Use evalue if probability not available
        elif self.evalue is not None:
            confidence = self._evalue_to_confidence(self.evalue)
            # Penalty for using evalue instead of probability in HHSearch
            confidence *= 0.9

        # Additional boost: Consider score if available
        if self.score is not None and self.score > 0:
            # HHSearch scores above 50 are typically very good
            score_factor = min(self.score / 100.0, 0.1)  # Max 10% boost
            confidence = min(confidence + score_factor, 1.0)

        return max(confidence, 0.0)

    def _calculate_blast_confidence(self) -> float:
        """Calculate confidence for BLAST evidence"""
        confidence = 0.0

        # Primary: Use evalue
        if self.evalue is not None:
            confidence = self._evalue_to_confidence(self.evalue)

            # Adjust based on identity if available
            if self.identity is not None:
                identity_factor = self.identity / 100.0 if self.identity > 1.0 else self.identity
                # High identity boosts confidence, low identity reduces it
                if identity_factor >= 0.5:
                    confidence *= (1.0 + (identity_factor - 0.5))  # Up to 50% boost
                else:
                    confidence *= identity_factor * 2  # Significant penalty for low identity

            # Adjust based on coverage if available
            if self.coverage is not None:
                coverage_factor = self.coverage / 100.0 if self.coverage > 1.0 else self.coverage
                # Good coverage is important for domain BLAST
                if coverage_factor >= 0.7:
                    confidence *= 1.1  # 10% boost for good coverage
                elif coverage_factor < 0.3:
                    confidence *= 0.8  # 20% penalty for poor coverage

        # Fallback: Use probability if available (some BLAST variants provide this)
        elif self.probability is not None:
            confidence = self.probability if self.probability <= 1.0 else self.probability / 100.0

        return max(min(confidence, 1.0), 0.0)

    def _calculate_chain_blast_confidence(self) -> float:
        """Calculate confidence for chain BLAST evidence (typically mapped to domains)"""
        # Start with regular BLAST confidence
        confidence = self._calculate_blast_confidence()

        # Apply penalty for chain-level mapping uncertainty
        confidence *= 0.85  # 15% penalty for mapping uncertainty

        # Additional penalty if this looks like a multi-domain chain
        if self.hsp_count is not None and self.hsp_count > 3:
            confidence *= 0.9  # Additional penalty for complex alignments

        return max(confidence, 0.0)

    def _calculate_generic_confidence(self) -> float:
        """Calculate confidence for unknown evidence types"""
        confidence = 0.0

        # Try each metric in order of preference
        if self.probability is not None:
            confidence = self.probability if self.probability <= 1.0 else self.probability / 100.0
        elif self.evalue is not None:
            confidence = self._evalue_to_confidence(self.evalue)
        elif self.score is not None and self.score > 0:
            # Assume score is normalized to ~100
            confidence = min(self.score / 100.0, 1.0)

        # Conservative penalty for unknown types
        confidence *= 0.8

        return max(min(confidence, 1.0), 0.0)

    def _evalue_to_confidence(self, evalue: float) -> float:
        """
        Convert E-value to confidence score

        Args:
            evalue: E-value (lower is better)

        Returns:
            float: Confidence score (0-1, higher is better)

        Formula based on sigmoid transformation:
        confidence = 1 / (1 + log10(evalue + 1))

        E-value ranges and typical confidence:
        - 1e-10: ~0.91 (excellent)
        - 1e-5:  ~0.83 (very good)
        - 1e-3:  ~0.75 (good)
        - 0.1:   ~0.50 (marginal)
        - 1.0:   ~0.33 (poor)
        - 10.0:  ~0.23 (very poor)
        """
        if evalue <= 0:
            return 1.0

        try:
            # Use log10 transformation with adjustment for very small e-values
            log_evalue = math.log10(evalue + 1e-100)  # Avoid log(0)

            # Sigmoid-like transformation
            # Good e-values (< 1e-3) get high confidence
            # Poor e-values (> 1) get low confidence
            confidence = 1.0 / (1.0 + abs(log_evalue) / 10.0)

            # Special handling for very good e-values
            if evalue < 1e-10:
                confidence = min(confidence * 1.1, 1.0)  # Small boost
            elif evalue > 10.0:
                confidence = confidence * 0.5  # Significant penalty

            return max(min(confidence, 1.0), 0.0)

        except (ValueError, OverflowError):
            # Fallback for edge cases
            return 0.5 if evalue < 1.0 else 0.1

    def set_confidence(self, confidence: float) -> None:
        """
        Explicitly set confidence value

        Args:
            confidence: Confidence value (0-1)

        Raises:
            ValueError: If confidence is not in valid range
        """
        if not 0.0 <= confidence <= 1.0:
            raise ValueError(f"Confidence must be between 0.0 and 1.0, got {confidence}")
        self.confidence = confidence

    def get_confidence_explanation(self) -> str:
        """
        Get human-readable explanation of how confidence was calculated

        Returns:
            str: Explanation of confidence calculation
        """
        if self.confidence is None:
            return "No confidence calculated"

        explanation = f"Confidence: {self.confidence:.3f} "

        if self.type == "hhsearch":
            if self.probability is not None:
                explanation += f"(from HHSearch probability: {self.probability})"
            elif self.evalue is not None:
                explanation += f"(from E-value: {self.evalue})"
        elif self.type in ("domain_blast", "blast"):
            if self.evalue is not None:
                explanation += f"(from BLAST E-value: {self.evalue}"
                if self.identity is not None:
                    explanation += f", identity: {self.identity}%"
                if self.coverage is not None:
                    explanation += f", coverage: {self.coverage}%"
                explanation += ")"
        elif self.type == "chain_blast":
            explanation += f"(chain BLAST with mapping penalty)"

        return explanation

    def recalculate_confidence(self) -> float:
        """
        Force recalculation of confidence even if already set

        Returns:
            float: New confidence value
        """
        old_confidence = self.confidence
        self.confidence = self._calculate_confidence()

        logger = logging.getLogger(__name__)
        logger.debug(f"Recalculated confidence: {old_confidence} -> {self.confidence}")


    @classmethod
    def from_blast_xml(cls, element: ET.Element, blast_type: str = "domain_blast") -> 'Evidence':
        """Create Evidence from BLAST hit XML element"""
        # Extract basic attributes
        hit_id = element.get("num", "")
        domain_id = element.get("domain_id", "")
        pdb_id = element.get("pdb_id", "")
        chain_id = element.get("chain_id", "")

        # Parse e-value from evalues attribute (may be comma-separated)
        evalues_str = element.get("evalues", "999")
        try:
            evalue = float(evalues_str.split(",")[0])
        except (ValueError, IndexError):
            evalue = 999.0

        # Parse HSP count
        try:
            hsp_count = int(element.get("hsp_count", "0"))
        except ValueError:
            hsp_count = 0

        # Extract query and hit regions
        query_range = ""
        hit_range = ""

        query_reg = element.find("query_reg")
        if query_reg is not None and query_reg.text:
            query_range = query_reg.text.strip()

        hit_reg = element.find("hit_reg")
        if hit_reg is not None and hit_reg.text:
            hit_range = hit_reg.text.strip()

        # Create evidence
        evidence = cls(
            type=blast_type,
            source_id=domain_id or hit_id,
            domain_id=domain_id,
            query_range=query_range,
            hit_range=hit_range,
            evalue=evalue,
            hsp_count=hsp_count,
            confidence=None,
            extra_attributes={
                "pdb_id": pdb_id,
                "chain_id": chain_id,
                "hit_id": hit_id,
                "discontinuous": element.get("discontinuous", "false").lower() == "true"
            }
        )

        return evidence

    @classmethod
    def from_hhsearch_xml(cls, element: ET.Element) -> 'Evidence':
        """Create Evidence from HHSearch hit XML element"""
        # Extract basic attributes
        hit_id = element.get("hit_id", "")
        domain_id = element.get("domain_id", "")

        # Parse numeric values
        try:
            probability = float(element.get("probability", "0"))
        except ValueError:
            probability = 0.0

        try:
            evalue = float(element.get("evalue", "999"))
        except ValueError:
            evalue = 999.0

        try:
            score = float(element.get("score", "0"))
        except ValueError:
            score = 0.0

        # Extract query and hit regions
        query_range = ""
        hit_range = ""

        query_reg = element.find("query_reg")
        if query_reg is not None and query_reg.text:
            query_range = query_reg.text.strip()

        hit_reg = element.find("hit_reg")
        if hit_reg is not None and hit_reg.text:
            hit_range = hit_reg.text.strip()

        # Create evidence
        evidence = cls(
            type="hhsearch",
            source_id=domain_id or hit_id,
            domain_id=domain_id,
            query_range=query_range,
            hit_range=hit_range,
            probability=probability,
            evalue=evalue,
            score=score,
            confidence=None,
            extra_attributes={
                "hit_id": hit_id,
                "num": element.get("num", "")
            }
        )

        return evidence

    @classmethod
    def from_xml(cls, element: ET.Element) -> 'Evidence':
        """Create Evidence from XML element - auto-detects type"""
        # Determine evidence type based on attributes
        if element.get("probability") is not None or element.get("score") is not None:
            # HHSearch hit
            return cls.from_hhsearch_xml(element)
        elif element.get("evalues") is not None or element.get("hsp_count") is not None:
            # BLAST hit - determine subtype
            if element.get("domain_id"):
                return cls.from_blast_xml(element, "domain_blast")
            else:
                return cls.from_blast_xml(element, "chain_blast")
        else:
            # Generic evidence - extract what we can
            return cls(
                type="unknown",
                source_id=element.get("source_id", ""),
                domain_id=element.get("domain_id", ""),
                query_range=element.get("query_range", ""),
                hit_range=element.get("hit_range", ""),
                confidence=float(element.get("confidence", "0"))
            )

    def to_xml(self) -> ET.Element:
        """Convert Evidence to XML element"""
        element = ET.Element("evidence")

        # Set basic attributes
        element.set("type", self.type)
        if self.source_id:
            element.set("source_id", self.source_id)
        if self.domain_id:
            element.set("domain_id", self.domain_id)
        element.set("confidence", f"{self.confidence:.4f}")

        # Set type-specific attributes
        if self.evalue is not None:
            element.set("evalue", str(self.evalue))
        if self.probability is not None:
            element.set("probability", str(self.probability))
        if self.score is not None:
            element.set("score", str(self.score))
        if self.hsp_count is not None:
            element.set("hsp_count", str(self.hsp_count))
        if self.identity is not None:
            element.set("identity", str(self.identity))
        if self.coverage is not None:
            element.set("coverage", str(self.coverage))

        # Set classification
        for cls_attr in ["t_group", "h_group", "x_group", "a_group"]:
            value = getattr(self, cls_attr)
            if value:
                element.set(cls_attr, value)

        # Add ranges as child elements
        if self.query_range:
            query_reg = ET.SubElement(element, "query_range")
            query_reg.text = self.query_range

        if self.hit_range:
            hit_reg = ET.SubElement(element, "hit_range")
            hit_reg.text = self.hit_range

        # Add extra attributes
        for key, value in self.extra_attributes.items():
            if value is not None:
                if isinstance(value, bool):
                    element.set(key, str(value).lower())
                else:
                    element.set(key, str(value))

        return element

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

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation"""
        result = {
            "type": self.type,
            "source_id": self.source_id,
            "domain_id": self.domain_id,
            "query_range": self.query_range,
            "hit_range": self.hit_range,
            "confidence": self.confidence
        }

        # Add optional numeric fields
        for field in ["evalue", "probability", "score", "identity", "coverage", "hsp_count"]:
            value = getattr(self, field)
            if value is not None:
                result[field] = value

        # Add classification
        for cls_attr in ["t_group", "h_group", "x_group", "a_group"]:
            value = getattr(self, cls_attr)
            if value:
                result[cls_attr] = value

        # Add extra attributes
        if self.extra_attributes:
            result.update(self.extra_attributes)

        return result

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'Evidence':
        """Create Evidence from dictionary"""
        # Extract core fields
        evidence = cls(
            type=data.get("type", "unknown"),
            source_id=data.get("source_id", ""),
            domain_id=data.get("domain_id", ""),
            query_range=data.get("query_range", ""),
            hit_range=data.get("hit_range", ""),
            confidence=data.get("confidence", 0.0)
        )

        # Set optional numeric fields
        for field in ["evalue", "probability", "score", "identity", "coverage", "hsp_count"]:
            if field in data:
                setattr(evidence, field, data[field])

        # Set classification
        for cls_attr in ["t_group", "h_group", "x_group", "a_group"]:
            if cls_attr in data:
                setattr(evidence, cls_attr, data[cls_attr])

        # Collect extra attributes
        extra_keys = set(data.keys()) - {
            "type", "source_id", "domain_id", "query_range", "hit_range", "confidence",
            "evalue", "probability", "score", "identity", "coverage", "hsp_count",
            "t_group", "h_group", "x_group", "a_group"
        }

        for key in extra_keys:
            evidence.extra_attributes[key] = data[key]

        return evidence

# Backward compatibility aliases
DomainEvidence = Evidence
BlastEvidence = Evidence 
HHSearchEvidence = Evidence
