# models/pipeline/evidence.py - CORRECTED VERSION
from dataclasses import dataclass, field
from typing import Dict, Any, Optional, List, Union
import xml.etree.ElementTree as ET
import math
import logging
from ecod.models.base import XmlSerializable

@dataclass
class Evidence(XmlSerializable):
    """Unified evidence model - consolidates DomainEvidence and Evidence classes"""
    type: str  # "blast", "hhsearch", "domain_blast", "chain_blast"
    source_id: str = ""
    domain_id: str = ""
    query_range: str = ""
    hit_range: str = ""
    confidence: Optional[float] = None  # None = auto-calculate, any other value = preserve

    # Classification data
    t_group: Optional[str] = None
    h_group: Optional[str] = None
    x_group: Optional[str] = None
    a_group: Optional[str] = None

    # Source-specific attributes
    evalue: Optional[float] = None
    probability: Optional[float] = None
    score: Optional[float] = None
    identity: Optional[float] = None
    coverage: Optional[float] = None
    hsp_count: Optional[int] = None

    # Additional attributes for extensibility
    extra_attributes: Dict[str, Any] = field(default_factory=dict)

    # Internal flag to track if confidence was explicitly set
    _confidence_explicitly_set: bool = field(default=False, init=False)

    def __post_init__(self):
        """Auto-calculate confidence from available metrics if not explicitly set"""
        # If confidence is None, it means auto-calculate
        # If confidence is any other value (including 0.0), it was explicitly set
        if self.confidence is None:
            self.confidence = self._calculate_confidence()
            self._confidence_explicitly_set = False
        else:
            self._confidence_explicitly_set = True

    def _calculate_confidence(self) -> float:
        """
        Calculate confidence score from available metrics

        Returns:
            float: Confidence score between 0.0 and 1.0
        """
        try:
            if self.type == "hhsearch":
                return self._calculate_hhsearch_confidence()
            elif self.type in ("domain_blast", "blast"):
                return self._calculate_blast_confidence()
            elif self.type == "chain_blast":
                return self._calculate_chain_blast_confidence()
            else:
                return self._calculate_generic_confidence()
        except Exception as e:
            logging.getLogger(__name__).warning(f"Error calculating confidence for {self.type} evidence: {e}")
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
        """
        Calculate confidence for unknown evidence types

        Priority order (most to least reliable):
        1. Probability (statistically calibrated)
        2. E-value (universally meaningful statistical measure)
        3. Score (method-dependent, least reliable)
        """
        confidence = 0.0

        # Priority 1: Try probability first (most reliable)
        if self.probability is not None:
            confidence = self.probability if self.probability <= 1.0 else self.probability / 100.0

        # Priority 2: Try E-value (statistically meaningful)
        elif self.evalue is not None:
            confidence = self._evalue_to_confidence(self.evalue)
            # Small penalty for not having probability in a generic context
            confidence *= 0.95

        # Priority 3: Fall back to score (least reliable)
        elif self.score is not None and self.score > 0:
            # Assume score is normalized to ~100, but apply significant penalty
            # since scores are often method-specific and not comparable
            confidence = min(self.score / 100.0, 1.0)
            confidence *= 0.7  # Larger penalty for unreliable scores

        # Conservative penalty for unknown types (they haven't been validated)
        confidence *= 0.8

        return max(min(confidence, 1.0), 0.0)

    def _evalue_to_confidence(self, evalue: float) -> float:
        """
        Convert E-value to confidence score

        E-value interpretation for bioinformatics:
        - E-value is expected number of hits by chance
        - Smaller E-values = better matches = higher confidence
        - E-value of 1e-5 means 1 in 100,000 chance of random match

        Args:
            evalue: E-value from BLAST or similar search

        Returns:
            float: Confidence score between 0.0 and 1.0
        """
        if evalue <= 0:
            return 1.0

        try:
            # Use negative log10 so smaller E-values give higher confidence
            neg_log_evalue = -math.log10(evalue + 1e-100)  # Avoid log(0)

            # Map to confidence using biologically meaningful thresholds
            if neg_log_evalue >= 10:  # E-value <= 1e-10 (excellent)
                # Scale from 0.9 to 1.0 for very significant hits
                excess = min(neg_log_evalue - 10, 10)  # Cap at 20 (-log10)
                confidence = 0.9 + 0.1 * (excess / 10)

            elif neg_log_evalue >= 5:  # 1e-10 < E-value <= 1e-5 (very good)
                # Scale from 0.7 to 0.9
                confidence = 0.7 + 0.2 * (neg_log_evalue - 5) / 5

            elif neg_log_evalue >= 2:  # 1e-5 < E-value <= 1e-2 (good)
                # Scale from 0.5 to 0.7
                confidence = 0.5 + 0.2 * (neg_log_evalue - 2) / 3

            elif neg_log_evalue >= 0:  # 1e-2 < E-value <= 1.0 (moderate)
                # Scale from 0.2 to 0.5
                confidence = 0.2 + 0.3 * neg_log_evalue / 2

            else:  # E-value > 1.0 (poor)
                # Exponential decay for poor hits
                confidence = 0.2 * math.exp(neg_log_evalue)  # neg_log_evalue is negative

            # Apply boost for exceptionally good E-values
            if evalue < 1e-10:
                confidence = min(confidence * 1.05, 1.0)  # Small boost
            # Apply penalty for poor E-values
            elif evalue > 10.0:
                confidence = confidence * 0.8  # Significant penalty

            return max(min(confidence, 1.0), 0.0)

        except (ValueError, OverflowError):
            # Fallback for edge cases
            return 0.5 if evalue < 1.0 else 0.1

    def set_confidence(self, confidence: float) -> None:
        """Explicitly set confidence value"""
        if not 0.0 <= confidence <= 1.0:
            raise ValueError(f"Confidence must be between 0.0 and 1.0, got {confidence}")
        self.confidence = confidence
        self._confidence_explicitly_set = True

    def get_confidence_explanation(self) -> str:
        """Get human-readable explanation of how confidence was calculated"""
        if self.confidence is None:
            return "No confidence calculated"

        explanation = f"Confidence: {self.confidence:.3f} "

        if self._confidence_explicitly_set:
            explanation += "(explicitly set)"
        elif self.type == "hhsearch":
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
        """Force recalculation of confidence even if already set"""
        old_confidence = self.confidence
        self.confidence = self._calculate_confidence()
        self._confidence_explicitly_set = False

        logger = logging.getLogger(__name__)
        logger.debug(f"Recalculated confidence: {old_confidence} -> {self.confidence}")

        return self.confidence

    @classmethod
    def from_blast_xml(cls, element: ET.Element, blast_type: str = "domain_blast") -> 'Evidence':
        """
        Create Evidence from BLAST hit XML element

        FIXED: Better error handling for invalid numeric values
        """
        # Extract basic attributes
        hit_id = element.get("num", "")
        domain_id = element.get("domain_id", "")
        pdb_id = element.get("pdb_id", "")
        chain_id = element.get("chain_id", "")

        # FIXED: Parse e-value with proper error handling
        evalues_str = element.get("evalues", "999")
        try:
            evalue = float(evalues_str.split(",")[0])
            # CRITICAL FIX: Handle inf and nan values
            if not math.isfinite(evalue):
                evalue = 999.0
        except (ValueError, IndexError):
            evalue = 999.0

        # Parse HSP count safely
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

        # Create evidence (confidence will be auto-calculated since confidence=None)
        evidence = cls(
            type=blast_type,
            source_id=domain_id or hit_id,  # Prefer domain_id for source_id
            domain_id=domain_id,
            query_range=query_range,
            hit_range=hit_range,
            evalue=evalue,
            hsp_count=hsp_count,
            confidence=None,  # Auto-calculate
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
        """Create Evidence from HHSearch hit XML element - FIXED source_id preservation"""
        # Extract basic attributes
        hit_id = element.get("hit_id", "")
        domain_id = element.get("domain_id", "")

        # CRITICAL FIX: Preserve original source_id if present, don't overwrite with domain_id
        source_id = element.get("source_id", "") or domain_id or hit_id

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

        # CRITICAL FIX: Check for explicit confidence in XML
        confidence_str = element.get("confidence")
        if confidence_str:
            try:
                confidence = float(confidence_str)
            except ValueError:
                confidence = None  # Will auto-calculate
        else:
            confidence = None  # Will auto-calculate

        # Create evidence
        evidence = cls(
            type="hhsearch",
            source_id=source_id,  # FIXED: Don't overwrite with domain_id
            domain_id=domain_id,
            query_range=query_range,
            hit_range=hit_range,
            probability=probability,
            evalue=evalue,
            score=score,
            confidence=confidence,  # FIXED: Use explicit confidence if present
            extra_attributes={
                "hit_id": hit_id,
                "num": element.get("num", "")
            }
        )

        # CRITICAL FIX: If confidence was in XML, mark it as explicitly set
        if confidence_str:
            evidence._confidence_explicitly_set = True

        return evidence

    @classmethod
    def from_xml(cls, element: ET.Element) -> 'Evidence':
        """Create Evidence from XML element - FIXED to preserve confidence and source_id"""

        # Trust the stored type if present
        stored_type = element.get("type")

        if stored_type == "hhsearch":
            return cls.from_hhsearch_xml(element)
        elif stored_type in ("domain_blast", "blast", "chain_blast"):
            return cls.from_blast_xml(element, stored_type)
        elif stored_type:
            # Parse as generic evidence with preserved confidence state
            return cls._from_xml_preserve_state(element, stored_type)
        else:
            # Auto-detect for legacy support
            return cls._from_xml_auto_detect(element)

    @classmethod
    def _from_xml_preserve_state(cls, element: ET.Element, evidence_type: str) -> 'Evidence':
        """Parse XML preserving original source_id and confidence state"""

        # CRITICAL FIX: Preserve original source_id, don't auto-substitute
        source_id = element.get("source_id", "")
        domain_id = element.get("domain_id", "")

        # Parse numeric values safely
        evalue = cls._parse_float_safe(element.get("evalue"), None)
        probability = cls._parse_float_safe(element.get("probability"), None)
        score = cls._parse_float_safe(element.get("score"), None)
        identity = cls._parse_float_safe(element.get("identity"), None)
        coverage = cls._parse_float_safe(element.get("coverage"), None)
        hsp_count = cls._parse_int_safe(element.get("hsp_count"), None)

        # CRITICAL FIX: Parse confidence preserving explicit state
        confidence_str = element.get("confidence")
        if confidence_str:
            confidence = cls._parse_float_safe(confidence_str, None)
            explicitly_set = True
        else:
            confidence = None  # Will auto-calculate
            explicitly_set = False

        # Extract ranges
        query_range = ""
        hit_range = ""

        for range_element in ["query_range", "query_reg"]:
            query_reg = element.find(range_element)
            if query_reg is not None and query_reg.text:
                query_range = query_reg.text.strip()
                break

        for range_element in ["hit_range", "hit_reg"]:
            hit_reg = element.find(range_element)
            if hit_reg is not None and hit_reg.text:
                hit_range = hit_reg.text.strip()
                break

        # Create evidence
        evidence = cls(
            type=evidence_type,
            source_id=source_id,  # FIXED: Preserve original source_id
            domain_id=domain_id,
            query_range=query_range,
            hit_range=hit_range,
            confidence=confidence,  # FIXED: Use explicit confidence if present
            evalue=evalue,
            probability=probability,
            score=score,
            identity=identity,
            coverage=coverage,
            hsp_count=hsp_count
        )

        # CRITICAL FIX: Set confidence state correctly
        evidence._confidence_explicitly_set = explicitly_set

        return evidence

    @classmethod
    def _parse_float_safe(cls, value_str, default=None):
        """Safely parse float with proper error handling"""
        if value_str is None:
            return default

        try:
            result = float(value_str)
            # Check for inf and nan
            if not (result == result and result != float('inf') and result != float('-inf')):
                return default
            return result
        except (ValueError, TypeError):
            return default

    @classmethod
    def _parse_int_safe(cls, value_str, default=None):
        """Safely parse int with proper error handling"""
        if value_str is None:
            return default

        try:
            return int(value_str)
        except (ValueError, TypeError):
            return default

    # TARGETED FIX 3: Ensure auto-detection fallback works correctly
    @classmethod
    def _from_xml_auto_detect(cls, element: ET.Element) -> 'Evidence':
        """Auto-detect evidence type (fallback for legacy XML)"""

        has_probability = element.get("probability") is not None
        has_score = element.get("score") is not None
        has_evalues = element.get("evalues") is not None
        has_evalue = element.get("evalue") is not None
        has_domain_id = element.get("domain_id") is not None
        has_hsp_count = element.get("hsp_count") is not None

        # HHSearch detection
        if has_probability:
            return cls.from_hhsearch_xml(element)

        # BLAST detection
        elif has_evalues or has_hsp_count:
            if has_domain_id:
                return cls.from_blast_xml(element, "domain_blast")
            else:
                return cls.from_blast_xml(element, "chain_blast")

        # Generic BLAST
        elif has_evalue:
            return cls.from_blast_xml(element, "blast")

        else:
            # Unknown - parse generically
            return cls._from_xml_preserve_state(element, "unknown")

    @classmethod
    def _from_xml_generic(cls, element: ET.Element, evidence_type: str) -> 'Evidence':
        """Parse XML as generic evidence with known type"""
        # Extract basic attributes - PRESERVE original source_id
        source_id = element.get("source_id", "")
        domain_id = element.get("domain_id", "")

        # Parse numeric values safely
        evalue = cls._parse_float_safe(element.get("evalue"), None)
        probability = cls._parse_float_safe(element.get("probability"), None)
        score = cls._parse_float_safe(element.get("score"), None)
        identity = cls._parse_float_safe(element.get("identity"), None)
        coverage = cls._parse_float_safe(element.get("coverage"), None)
        hsp_count = cls._parse_int_safe(element.get("hsp_count"), None)

        # CRITICAL FIX: Parse confidence but preserve auto-calculation state
        confidence_str = element.get("confidence")
        if confidence_str:
            confidence = cls._parse_float_safe(confidence_str, None)
            # If we loaded confidence from XML, it was explicitly stored
            explicitly_set = True
        else:
            confidence = None  # Will auto-calculate
            explicitly_set = False

        # Extract ranges
        query_range = ""
        hit_range = ""

        query_reg = element.find("query_range") or element.find("query_reg")
        if query_reg is not None and query_reg.text:
            query_range = query_reg.text.strip()

        hit_reg = element.find("hit_range") or element.find("hit_reg")
        if hit_reg is not None and hit_reg.text:
            hit_range = hit_reg.text.strip()

        # Extract classification
        t_group = element.get("t_group")
        h_group = element.get("h_group")
        x_group = element.get("x_group")
        a_group = element.get("a_group")

        # Create evidence
        evidence = cls(
            type=evidence_type,
            source_id=source_id,
            domain_id=domain_id,
            query_range=query_range,
            hit_range=hit_range,
            confidence=confidence,
            evalue=evalue,
            probability=probability,
            score=score,
            identity=identity,
            coverage=coverage,
            hsp_count=hsp_count,
            t_group=t_group,
            h_group=h_group,
            x_group=x_group,
            a_group=a_group
        )

        # CRITICAL FIX: Set the confidence state correctly
        evidence._confidence_explicitly_set = explicitly_set

        # Extract extra attributes (anything not in core fields)
        core_attrs = {
            "type", "source_id", "domain_id", "confidence", "evalue",
            "probability", "score", "identity", "coverage", "hsp_count",
            "t_group", "h_group", "x_group", "a_group"
        }

        for key, value in element.attrib.items():
            if key not in core_attrs:
                evidence.extra_attributes[key] = value

        return evidence

    @classmethod
    def _from_xml_auto_detect(cls, element: ET.Element) -> 'Evidence':
        """Auto-detect evidence type from XML (fallback for legacy data)"""
        # IMPROVED: Better auto-detection logic
        has_probability = element.get("probability") is not None
        has_score = element.get("score") is not None
        has_evalues = element.get("evalues") is not None
        has_evalue = element.get("evalue") is not None
        has_domain_id = element.get("domain_id") is not None
        has_hsp_count = element.get("hsp_count") is not None

        # HHSearch: has probability OR (has score AND no evalues)
        if has_probability or (has_score and not has_evalues and not has_evalue):
            return cls.from_hhsearch_xml(element)

        # BLAST: has evalues OR evalue OR hsp_count
        elif has_evalues or has_hsp_count:
            if has_domain_id:
                return cls.from_blast_xml(element, "domain_blast")
            else:
                return cls.from_blast_xml(element, "chain_blast")

        # Generic BLAST: has evalue but no other BLAST indicators
        elif has_evalue:
            return cls.from_blast_xml(element, "blast")

        else:
            # Completely unknown - parse as generic
            return cls._from_xml_generic(element, "unknown")

    @classmethod
    def _parse_float_safe(cls, value_str: str, default=None):
        """Safely parse float with proper error handling"""
        if value_str is None:
            return default

        try:
            result = float(value_str)
            # CRITICAL FIX: Check for inf and nan
            if not math.isfinite(result):
                return default
            return result
        except (ValueError, TypeError):
            return default

    @classmethod
    def _parse_int_safe(cls, value_str: str, default=None):
        """Safely parse int with proper error handling"""
        if value_str is None:
            return default

        try:
            return int(value_str)
        except (ValueError, TypeError):
            return default

    def to_xml(self) -> ET.Element:
        """Convert Evidence to XML element"""
        element = ET.Element("evidence")

        # Set basic attributes
        element.set("type", self.type)
        if self.source_id:
            element.set("source_id", self.source_id)
        if self.domain_id:
            element.set("domain_id", self.domain_id)
        if self.confidence is not None:
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
            confidence=None,  # Auto-calculate
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
            score=getattr(hit, "score", 0.0),
            confidence=None  # Auto-calculate
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
    def from_dict(cls, data: dict) -> 'Evidence':
        """
        Create Evidence from dictionary

        FIXED: Preserve confidence state correctly
        """
        # Extract core fields
        evidence_type = data.get("type", "unknown")
        source_id = data.get("source_id", "")
        domain_id = data.get("domain_id", "")
        query_range = data.get("query_range", "")
        hit_range = data.get("hit_range", "")

        # CRITICAL FIX: Handle confidence state properly
        confidence = data.get("confidence")
        # If confidence is explicitly in the dict, it was explicitly set
        confidence_was_explicit = "confidence" in data and confidence is not None

        evidence = cls(
            type=evidence_type,
            source_id=source_id,
            domain_id=domain_id,
            query_range=query_range,
            hit_range=hit_range,
            confidence=confidence  # None will auto-calculate, value will be preserved
        )

        # Set optional numeric fields
        for field in ["evalue", "probability", "score", "identity", "coverage", "hsp_count"]:
            if field in data and data[field] is not None:
                setattr(evidence, field, data[field])

        # Set classification
        for cls_attr in ["t_group", "h_group", "x_group", "a_group"]:
            if cls_attr in data and data[cls_attr]:
                setattr(evidence, cls_attr, data[cls_attr])

        # CRITICAL FIX: Set confidence state correctly
        if confidence_was_explicit:
            evidence._confidence_explicitly_set = True

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
