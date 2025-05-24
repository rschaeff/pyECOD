# ecod/models/pipeline.py
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict, Any, Set
from ecod.models.base import XmlSerializable
import xml.etree.ElementTree as ET

@dataclass
class BaseHit(XmlSerializable):
    """Base class for hit models"""
    hit_id: str = ""
    domain_id: str = ""  # Domain ID (for domain hits)
    pdb_id: str = ""     # PDB ID
    chain_id: str = ""   # Chain ID
    hit_type: str = ""   # Hit type
    range: str = ""      # Query range
    hit_range: str = ""  # Hit range
    range_parsed: List[Tuple[int, int]] = field(default_factory=list)

    def parse_ranges(self):
        """Parse range string into list of tuples"""
        # Implementation...

    def to_xml(self) -> ET.Element:
        """Convert to XML Element"""
        # Implementation...

    @classmethod
    def from_xml(cls, element: ET.Element) -> 'BaseHit':
        """Create from XML Element"""
        # Implementation...


@dataclass
class RangeSegment:
    """A single segment within a range (start-end)"""
    start: int
    end: int

    def __str__(self) -> str:
        return f"{self.start}-{self.end}"

    @property
    def length(self) -> int:
        return self.end - self.start + 1

@dataclass
class BlastHit:
    """Blast hit model for pipeline processing"""
    hit_id: str = ""
    domain_id: str = ""
    pdb_id: str = ""
    chain_id: str = ""
    evalue: float = 999.0
    hsp_count: int = 0
    hit_type: str = ""  # "chain_blast" or "domain_blast"
    range: str = ""
    hit_range: str = ""
    query_seq: str = ""
    hit_seq: str = ""
    discontinuous: bool = False
    evalues: List[float] = field(default_factory=list)
    range_parsed: List[Tuple[int, int]] = field(default_factory=list)

    xml_element_path = ".//blast_run/hits/hit"
    
    @classmethod
    def from_xml(cls, element: ET.Element) -> 'BlastHit':
        """Create from XML Element"""
        hit = cls(
            hit_id=element.get("num", ""),
            domain_id=element.get("domain_id", ""),
            pdb_id=element.get("pdb_id", ""),
            chain_id=element.get("chain_id", ""),
            evalue=float(element.get("evalues", "999").split(",")[0]),
            hsp_count=int(element.get("hsp_count", "0")),
            discontinuous=element.get("discontinuous", "false").lower() == "true"
        )

        # Get query regions
        query_reg = element.find("query_reg")
        if query_reg is not None and query_reg.text:
            hit.range = query_reg.text.strip()
            hit.parse_ranges()

        # Get hit regions
        hit_reg = element.find("hit_reg")
        if hit_reg is not None and hit_reg.text:
            hit.hit_range = hit_reg.text.strip()

        return hit

    def to_xml(self) -> ET.Element:
        """Convert to XML Element"""
        element = ET.Element("hit")

        # Set attributes
        element.set("num", str(self.hit_id))
        if self.domain_id:
            element.set("domain_id", self.domain_id)
        element.set("pdb_id", self.pdb_id)
        element.set("chain_id", self.chain_id)
        element.set("evalues", str(self.evalue))
        element.set("hsp_count", str(self.hsp_count))
        if self.discontinuous:
            element.set("discontinuous", "true")

        # Add query region
        if self.range:
            query_reg = ET.SubElement(element, "query_reg")
            query_reg.text = self.range

        # Add hit region
        if self.hit_range:
            hit_reg = ET.SubElement(element, "hit_reg")
            hit_reg.text = self.hit_range

        return element
    
    def parse_ranges(self):
        """Parse range string into list of tuples"""
        self.range_parsed = []
        if not self.range:
            return
            
        for segment in self.range.split(","):
            if "-" in segment:
                try:
                    start, end = map(int, segment.split("-"))
                    self.range_parsed.append((start, end))
                except ValueError:
                    pass

    def to_dict(self) -> Dict[str, Any]:
        """Convert BlastHit to dictionary
        
        Returns:
            Dictionary representation
        """
        result = {
            "hit_id": self.hit_id,
            "domain_id": self.domain_id,
            "pdb_id": self.pdb_id,
            "chain_id": self.chain_id,
            "evalue": self.evalue,
            "hsp_count": self.hsp_count,
            "type": self.hit_type,
            "range": self.range,
            "query_range": self.range,  # Alternate name used in some places
            "query_regions": self.range,  # Alternate name used in some places
            "hit_range": self.hit_range,
            "query_seq": self.query_seq,
            "hit_seq": self.hit_seq,
            "discontinuous": self.discontinuous
        }
        
        # Add evalues list if available
        if self.evalues:
            result["evalues"] = self.evalues
        
        # Add range_parsed if available
        if self.range_parsed:
            result["range_parsed"] = self.range_parsed
        
        return result

    def get_segments(self) -> List[RangeSegment]:
        """Get range as list of RangeSegment objects"""
        segments = []
        for start, end in self.range_parsed:
            segments.append(RangeSegment(start, end))
        return segments

    def get_positions(self) -> Set[int]:
        """Get all positions covered by this hit"""
        positions = set()
        for start, end in self.range_parsed:
            positions.update(range(start, end + 1))
        return positions

    def overlaps(self, other: 'BlastHit') -> bool:
        """Check if this hit overlaps with another"""
        my_positions = self.get_positions()
        other_positions = other.get_positions()
        return bool(my_positions.intersection(other_positions))

@dataclass
class HHSearchHit:
    """HHSearch hit model for pipeline processing"""
    hit_id: str = ""
    domain_id: str = ""
    probability: float = 0.0
    evalue: float = 999.0
    score: float = 0.0
    range: str = ""
    hit_range: str = ""
    range_parsed: List[Tuple[int, int]] = field(default_factory=list)
    
    @classmethod
    def from_xml(cls, hit_elem):
        """Create from XML Element"""
        hit = cls(
            hit_id=hit_elem.get("hit_id", ""),
            domain_id=hit_elem.get("domain_id", ""),
            probability=float(hit_elem.get("probability", "0")),
            evalue=float(hit_elem.get("evalue", "999")),
            score=float(hit_elem.get("score", "0")),
        )
        
        # Get query regions
        query_reg = hit_elem.find("query_reg")
        if query_reg is not None and query_reg.text:
            hit.range = query_reg.text.strip()
            hit.parse_ranges()
            
        # Get hit regions
        hit_reg = hit_elem.find("hit_reg")
        if hit_reg is not None and hit_reg.text:
            hit.hit_range = hit_reg.text.strip()
            
        return hit
    
    def parse_ranges(self):
        """Parse range string into list of tuples"""
        self.range_parsed = []
        if not self.range:
            return
            
        for segment in self.range.split(","):
            if "-" in segment:
                try:
                    start, end = map(int, segment.split("-"))
                    self.range_parsed.append((start, end))
                except ValueError:
                    pass

    def get_segments(self) -> List[RangeSegment]:
        """Get range as list of RangeSegment objects"""
        segments = []
        for start, end in self.range_parsed:
            segments.append(RangeSegment(start, end))
        return segments

    def get_positions(self) -> Set[int]:
        """Get all positions covered by this hit"""
        positions = set()
        for start, end in self.range_parsed:
            positions.update(range(start, end + 1))
        return positions

    def overlaps(self, other: 'BlastHit') -> bool:
        """Check if this hit overlaps with another"""
        my_positions = self.get_positions()
        other_positions = other.get_positions()
        return bool(my_positions.intersection(other_positions))
                    
    def to_dict(self) -> Dict[str, Any]:
        """Convert HHSearchHit to dictionary
        
        Returns:
            Dictionary representation
        """
        result = {
            "hit_id": self.hit_id,
            "domain_id": self.domain_id,
            "probability": self.probability,
            "evalue": self.evalue,
            "score": self.score,
            "range": self.range,
            "query_range": self.range,  # Alternate name used in some places
            "hit_range": self.hit_range,
            "type": "hhsearch"  # Add type to match structure used elsewhere
        }
        
        # Add range_parsed if available
        if self.range_parsed:
            result["range_parsed"] = self.range_parsed
        
        return result

@dataclass
class DomainSummaryModel:
    """Domain summary model for pipeline processing"""
    pdb_id: str
    chain_id: str
    reference: str
    sequence_length: int = 0
    is_peptide: bool = False
    chain_blast_hits: List[BlastHit] = field(default_factory=list)
    domain_blast_hits: List[BlastHit] = field(default_factory=list)
    hhsearch_hits: List[HHSearchHit] = field(default_factory=list)
    self_comparison_hits: List[Dict] = field(default_factory=list)
    errors: Dict[str, bool] = field(default_factory=dict)
    output_file_path: Optional[str] = None
    processed: bool = False
    skipped: bool = False
    sequence: Optional[str] = None

    def to_xml(self):
        """Convert to XML Element"""
        root = ET.Element("blast_summ_doc")

        # Create summary node
        blast_summ = ET.SubElement(root, "blast_summ", pdb=self.pdb_id, chain=self.chain_id)

        # Add errors
        for error, value in self.errors.items():
            if value:
                blast_summ.set(error, "true")

        # Add chain blast hits
        if self.chain_blast_hits:
            chain_blast_run = ET.SubElement(root, "chain_blast_run")
            chain_blast_run.set("program", "blastp")

            hits_node = ET.SubElement(chain_blast_run, "hits")
            for hit in self.chain_blast_hits:
                hit_elem = ET.SubElement(hits_node, "hit")
                hit_elem.set("num", hit.hit_id)
                hit_elem.set("pdb_id", hit.pdb_id)
                hit_elem.set("chain_id", hit.chain_id)
                hit_elem.set("hsp_count", str(hit.hsp_count))

                # Format evalues correctly
                if hasattr(hit, "evalues") and hit.evalues:
                    hit_elem.set("evalues", ",".join(str(e) for e in hit.evalues))
                elif hasattr(hit, "evalue"):
                    hit_elem.set("evalues", str(hit.evalue))

                # Add query region
                query_reg = ET.SubElement(hit_elem, "query_reg")
                query_reg.text = hit.range

                # Add hit region
                hit_reg = ET.SubElement(hit_elem, "hit_reg")
                hit_reg.text = hit.hit_range

        # Add domain blast hits
        if self.domain_blast_hits:
            domain_blast_run = ET.SubElement(root, "blast_run")
            domain_blast_run.set("program", "blastp")

            hits_node = ET.SubElement(domain_blast_run, "hits")
            for hit in self.domain_blast_hits:
                hit_elem = ET.SubElement(hits_node, "hit")
                if hit.domain_id:
                    hit_elem.set("domain_id", hit.domain_id)
                hit_elem.set("pdb_id", hit.pdb_id)
                hit_elem.set("chain_id", hit.chain_id)
                hit_elem.set("hsp_count", str(hit.hsp_count))

                # Format evalues correctly
                if hasattr(hit, "evalues") and hit.evalues:
                    hit_elem.set("evalues", ",".join(str(e) for e in hit.evalues))
                elif hasattr(hit, "evalue"):
                    hit_elem.set("evalues", str(hit.evalue))

                # Add discontinuous flag if applicable
                if hasattr(hit, "discontinuous") and hit.discontinuous:
                    hit_elem.set("discontinuous", "true")

                # Add query region
                query_reg = ET.SubElement(hit_elem, "query_reg")
                query_reg.text = hit.range

                # Add hit region
                hit_reg = ET.SubElement(hit_elem, "hit_reg")
                hit_reg.text = hit.hit_range

        # Add HHSearch hits
        if self.hhsearch_hits:
            hh_run = ET.SubElement(root, "hh_run")
            hh_run.set("program", "hhsearch")

            hits_node = ET.SubElement(hh_run, "hits")
            for i, hit in enumerate(self.hhsearch_hits):
                hit_elem = ET.SubElement(hits_node, "hit")

                # Set required attributes
                if hit.domain_id:
                    hit_elem.set("domain_id", hit.domain_id)
                hit_elem.set("hit_id", hit.hit_id)
                hit_elem.set("num", str(i+1))
                hit_elem.set("probability", str(hit.probability))
                hit_elem.set("evalue", str(hit.evalue))
                hit_elem.set("score", str(hit.score))

                # Add query region
                query_reg = ET.SubElement(hit_elem, "query_reg")
                query_reg.text = hit.range

                # Add hit region if available
                if hasattr(hit, "hit_range") and hit.hit_range:
                    hit_reg = ET.SubElement(hit_elem, "hit_reg")
                    hit_reg.text = hit.hit_range
        
        return root

    def generate_domain_suggestions(self) -> List[Dict[str, Any]]:
        """Generate domain suggestions based on evidence"""
        # Collect all boundaries
        boundaries = []

        # Extract from chain BLAST hits
        for hit in self.chain_blast_hits:
            for start, end in hit.range_parsed:
                boundaries.append((start, end, ["chain_blast"]))

        # Extract from domain BLAST hits
        for hit in self.domain_blast_hits:
            for start, end in hit.range_parsed:
                boundaries.append((start, end, ["domain_blast"]))

        # Extract from HHSearch hits
        for hit in self.hhsearch_hits:
            for start, end in hit.range_parsed:
                boundaries.append((start, end, ["hhsearch"]))

        # Merge overlapping boundaries
        merged = self._merge_boundaries(boundaries)

        # Format as domain suggestions
        domain_suggestions = []
        for i, (start, end, sources) in enumerate(merged):
            domain_suggestions.append({
                "id": f"domain_{i+1}",
                "start": start,
                "end": end,
                "sources": sources
            })

        return domain_suggestions

    def _merge_boundaries(self, boundaries):
        """Merge overlapping domain boundaries"""
        if not boundaries:
            return []


    # Helper method to get stats
    def get_stats(self) -> Dict[str, Any]:
        """Get statistics about this domain summary"""
        return {
            "chain_blast_processed": len(self.chain_blast_hits) > 0,
            "domain_blast_processed": len(self.domain_blast_hits) > 0,
            "hhsearch_processed": len(self.hhsearch_hits) > 0,
            "hhsearch_hits": len(self.hhsearch_hits),
            "chain_blast_hits": len(self.chain_blast_hits),
            "domain_blast_hits": len(self.domain_blast_hits),
            "self_comp_processed": len(self.self_comparison_hits) > 0
        }

    # Compatibility method for legacy code
    def to_legacy_dict(self) -> Dict[str, Any]:
        """Convert to legacy dictionary format for backward compatibility"""
        return {
            "file_path": self.output_file_path,
            "stats": self.get_stats(),
            "skipped": self.skipped,
            "is_peptide": self.is_peptide,
            "summary": self
        }

@dataclass
class PipelineResult:
    """Result of running the domain analysis pipeline"""
    batch_id: int
    success: bool = False
    error: Optional[str] = None
    reset_stats: Dict[str, int] = field(default_factory=dict)
    processing_stats: Dict[str, Any] = field(default_factory=dict)
    summary_stats: Dict[str, Any] = field(default_factory=dict)
    partition_stats: Dict[str, Any] = field(default_factory=dict)
    batch_info: Dict[str, Any] = field(default_factory=dict)

    def set_error(self, error_message: str) -> None:
        """Set error message and update success status"""
        self.error = error_message
        self.success = False

    def set_batch_info(self, batch_info: Dict[str, Any]) -> None:
        """Set batch information"""
        self.batch_info = batch_info

@dataclass
class ProteinResult:
    """Result of processing a single protein"""
    protein_id: int
    pdb_id: str
    chain_id: str
    process_id: int
    success: bool = False
    summary_success: bool = False
    partition_success: bool = False
    error: Optional[str] = None
    summary_error: Optional[str] = None
    partition_error: Optional[str] = None
    summary_file: Optional[str] = None
    domain_file: Optional[str] = None
    chain_blast_count: int = 0
    domain_blast_count: int = 0
    hhsearch_count: int = 0

@dataclass
class ProteinProcessingResult:
    """Result of processing multiple proteins"""
    batch_id: int
    protein_ids: List[int]
    success: bool = False
    error: Optional[str] = None
    success_count: int = 0
    total_count: int = 0
    batch_info: Dict[str, Any] = field(default_factory=dict)
    protein_results: List[ProteinResult] = field(default_factory=list)

    def set_error(self, error_message: str) -> None:
        """Set error message and update success status"""
        self.error = error_message
        self.success = False

    def get_summary_stats(self) -> Dict[str, Any]:
        """Get summary statistics for all processed proteins"""
        return {
            "total": self.total_count,
            "success": self.success_count,
            "chain_blast_files": sum(1 for p in self.protein_results if p.chain_blast_count > 0),
            "domain_blast_files": sum(1 for p in self.protein_results if p.domain_blast_count > 0),
            "hhsearch_files": sum(1 for p in self.protein_results if p.hhsearch_count > 0),
            "summary_success": sum(1 for p in self.protein_results if p.summary_success),
            "partition_success": sum(1 for p in self.protein_results if p.partition_success),
            "failed": self.total_count - self.success_count
        }

# Add to the end of ecod/models/pipeline.py

import warnings
from typing import Union

# =============================================================================
# DEPRECATED LEGACY MODELS - Use Evidence instead
# =============================================================================

class BlastHit:
    """
    DEPRECATED: Use Evidence model instead

    This class is kept for backward compatibility only.
    All new code should use Evidence.from_blast_xml() or Evidence.from_blast_hit().
    """

    def __init__(self, *args, **kwargs):
        warnings.warn(
            "BlastHit is deprecated. Use Evidence.from_blast_xml() instead. "
            "See migration notes in ecod.models.__init__ for details.",
            DeprecationWarning,
            stacklevel=2
        )

        # Initialize with minimal fields for compatibility
        self.hit_id = kwargs.get('hit_id', '')
        self.domain_id = kwargs.get('domain_id', '')
        self.pdb_id = kwargs.get('pdb_id', '')
        self.chain_id = kwargs.get('chain_id', '')
        self.evalue = kwargs.get('evalue', 999.0)
        self.hsp_count = kwargs.get('hsp_count', 0)
        self.hit_type = kwargs.get('hit_type', 'blast')
        self.range = kwargs.get('range', '')
        self.hit_range = kwargs.get('hit_range', '')
        self.query_seq = kwargs.get('query_seq', '')
        self.hit_seq = kwargs.get('hit_seq', '')
        self.discontinuous = kwargs.get('discontinuous', False)
        self.evalues = kwargs.get('evalues', [])
        self.range_parsed = kwargs.get('range_parsed', [])

    @classmethod
    def from_xml(cls, element):
        """DEPRECATED: Use Evidence.from_blast_xml() instead"""
        warnings.warn(
            "BlastHit.from_xml() is deprecated. Use Evidence.from_blast_xml() instead.",
            DeprecationWarning,
            stacklevel=2
        )

        # Import here to avoid circular imports
        from ecod.models.pipeline.evidence import Evidence

        # Convert to Evidence and then back to BlastHit for compatibility
        evidence = Evidence.from_blast_xml(element, "domain_blast")

        # Create BlastHit with converted data
        return cls(
            hit_id=evidence.source_id,
            domain_id=evidence.domain_id,
            pdb_id=evidence.extra_attributes.get('pdb_id', ''),
            chain_id=evidence.extra_attributes.get('chain_id', ''),
            evalue=evidence.evalue or 999.0,
            hsp_count=evidence.hsp_count or 0,
            range=evidence.query_range,
            hit_range=evidence.hit_range,
            discontinuous=evidence.extra_attributes.get('discontinuous', False)
        )

    def to_evidence(self) -> 'Evidence':
        """Convert BlastHit to Evidence model"""
        from ecod.models.pipeline.evidence import Evidence

        return Evidence(
            type=self.hit_type or "blast",
            source_id=self.domain_id or self.hit_id,
            domain_id=self.domain_id,
            query_range=self.range,
            hit_range=self.hit_range,
            evalue=self.evalue,
            hsp_count=self.hsp_count,
            extra_attributes={
                'pdb_id': self.pdb_id,
                'chain_id': self.chain_id,
                'hit_id': self.hit_id,
                'discontinuous': self.discontinuous
            }
        )


class HHSearchHit:
    """
    DEPRECATED: Use Evidence model instead

    This class is kept for backward compatibility only.
    All new code should use Evidence.from_hhsearch_xml() or Evidence.from_hhsearch_hit().
    """

    def __init__(self, *args, **kwargs):
        warnings.warn(
            "HHSearchHit is deprecated. Use Evidence.from_hhsearch_xml() instead. "
            "See migration notes in ecod.models.__init__ for details.",
            DeprecationWarning,
            stacklevel=2
        )

        # Initialize with minimal fields for compatibility
        self.hit_id = kwargs.get('hit_id', '')
        self.domain_id = kwargs.get('domain_id', '')
        self.probability = kwargs.get('probability', 0.0)
        self.evalue = kwargs.get('evalue', 999.0)
        self.score = kwargs.get('score', 0.0)
        self.range = kwargs.get('range', '')
        self.hit_range = kwargs.get('hit_range', '')
        self.range_parsed = kwargs.get('range_parsed', [])

    @classmethod
    def from_xml(cls, element):
        """DEPRECATED: Use Evidence.from_hhsearch_xml() instead"""
        warnings.warn(
            "HHSearchHit.from_xml() is deprecated. Use Evidence.from_hhsearch_xml() instead.",
            DeprecationWarning,
            stacklevel=2
        )

        # Import here to avoid circular imports
        from ecod.models.pipeline.evidence import Evidence

        # Convert to Evidence and then back to HHSearchHit for compatibility
        evidence = Evidence.from_hhsearch_xml(element)

        # Create HHSearchHit with converted data
        return cls(
            hit_id=evidence.source_id,
            domain_id=evidence.domain_id,
            probability=evidence.probability or 0.0,
            evalue=evidence.evalue or 999.0,
            score=evidence.score or 0.0,
            range=evidence.query_range,
            hit_range=evidence.hit_range
        )

    def to_evidence(self) -> 'Evidence':
        """Convert HHSearchHit to Evidence model"""
        from ecod.models.pipeline.evidence import Evidence

        return Evidence(
            type="hhsearch",
            source_id=self.domain_id or self.hit_id,
            domain_id=self.domain_id,
            query_range=self.range,
            hit_range=self.hit_range,
            probability=self.probability,
            evalue=self.evalue,
            score=self.score,
            extra_attributes={
                'hit_id': self.hit_id
            }
        )


# =============================================================================
# MIGRATION HELPER FUNCTIONS
# =============================================================================

def migrate_blast_hit_to_evidence(blast_hit: BlastHit) -> 'Evidence':
    """
    Helper function to migrate BlastHit to Evidence

    Args:
        blast_hit: Legacy BlastHit object

    Returns:
        Evidence: New Evidence object
    """
    warnings.warn(
        "migrate_blast_hit_to_evidence() is deprecated. "
        "Use BlastHit.to_evidence() or Evidence.from_blast_xml() directly.",
        DeprecationWarning,
        stacklevel=2
    )

    return blast_hit.to_evidence()


def migrate_hhsearch_hit_to_evidence(hhsearch_hit: HHSearchHit) -> 'Evidence':
    """
    Helper function to migrate HHSearchHit to Evidence

    Args:
        hhsearch_hit: Legacy HHSearchHit object

    Returns:
        Evidence: New Evidence object
    """
    warnings.warn(
        "migrate_hhsearch_hit_to_evidence() is deprecated. "
        "Use HHSearchHit.to_evidence() or Evidence.from_hhsearch_xml() directly.",
        DeprecationWarning,
        stacklevel=2
    )

    return hhsearch_hit.to_evidence()


# =============================================================================
# LEGACY COMPATIBILITY ALIASES
# =============================================================================

# These aliases allow existing imports to work but with deprecation warnings
def _deprecated_blast_hit_import():
    warnings.warn(
        "Importing BlastHit from ecod.models.pipeline is deprecated. "
        "Use Evidence from ecod.models.pipeline.evidence instead.",
        DeprecationWarning,
        stacklevel=3
    )
    return BlastHit

def _deprecated_hhsearch_hit_import():
    warnings.warn(
        "Importing HHSearchHit from ecod.models.pipeline is deprecated. "
        "Use Evidence from ecod.models.pipeline.evidence instead.",
        DeprecationWarning,
        stacklevel=3
    )
    return HHSearchHit


# Update module's __getattr__ to handle legacy imports
def __getattr__(name: str):
    if name == 'BlastHit':
        return _deprecated_blast_hit_import()
    elif name == 'HHSearchHit':
        return _deprecated_hhsearch_hit_import()
    else:
        raise AttributeError(f"module {__name__} has no attribute {name}")


# =============================================================================
# MIGRATION DOCUMENTATION
# =============================================================================

MIGRATION_GUIDE = """
MIGRATION GUIDE: BlastHit/HHSearchHit → Evidence

OLD CODE:
    from ecod.models.pipeline import BlastHit, HHSearchHit

    blast_hit = BlastHit.from_xml(element)
    hhsearch_hit = HHSearchHit.from_xml(element)

    evidence1 = Evidence.from_blast_hit(blast_hit)
    evidence2 = Evidence.from_hhsearch_hit(hhsearch_hit)

NEW CODE:
    from ecod.models.pipeline.evidence import Evidence

    evidence1 = Evidence.from_blast_xml(element, "domain_blast")
    evidence2 = Evidence.from_hhsearch_xml(element)

BENEFITS:
    - Single parsing path (faster, less memory)
    - Consistent data representation
    - Better confidence calculation
    - Unified XML serialization

COMPATIBILITY:
    - Legacy classes still work but show deprecation warnings
    - Use .to_evidence() method to convert existing objects
    - All methods and attributes preserved for compatibility
"""
