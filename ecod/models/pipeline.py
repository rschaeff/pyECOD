# ecod/models/pipeline.py
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict, Any, Set
from ecod.models.base import XmlSerializable
from ecod.models.pipeline.evidence import Evidence  # ADD THIS LINE
import xml.etree.ElementTree as ET
import warnings

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
    
    def __post_init__(self):
        warnings.warn(
            "BlastHit is deprecated. Use Evidence from ecod.models.pipeline.evidence instead.",
            DeprecationWarning,
            stacklevel=3  # Adjust to point to caller's caller
        )

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

    def __post_init__(self):
        warnings.warn(
            "HHSearchHit is deprecated. Use Evidence from ecod.models.pipeline.evidence instead.",
            DeprecationWarning,
            stacklevel=3  # Adjust to point to caller's caller
        )
    
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

  # In DomainSummaryModel class
    def to_xml(self):
        """Convert to XML Element - updated to handle Evidence objects"""
        root = ET.Element("blast_summ_doc")

        # Create summary node
        blast_summ = ET.SubElement(root, "blast_summ", pdb=self.pdb_id, chain=self.chain_id)

        # Add errors
        for error, value in self.errors.items():
            if value:
                blast_summ.set(error, "true")

        # Add chain blast hits (now Evidence objects)
        if self.chain_blast_hits:
            chain_blast_run = ET.SubElement(root, "chain_blast_run")
            chain_blast_run.set("program", "blastp")

            hits_node = ET.SubElement(chain_blast_run, "hits")
            for i, evidence in enumerate(self.chain_blast_hits):
                # Convert Evidence to BLAST hit XML format
                hit_elem = self._evidence_to_blast_xml(evidence, i + 1)
                hits_node.append(hit_elem)

        # Add domain blast hits (now Evidence objects)
        if self.domain_blast_hits:
            domain_blast_run = ET.SubElement(root, "blast_run")
            domain_blast_run.set("program", "blastp")

            hits_node = ET.SubElement(domain_blast_run, "hits")
            for i, evidence in enumerate(self.domain_blast_hits):
                hit_elem = self._evidence_to_blast_xml(evidence, i + 1)
                hits_node.append(hit_elem)

        # Add HHSearch hits (now Evidence objects)
        if self.hhsearch_hits:
            hh_run = ET.SubElement(root, "hh_run")
            hh_run.set("program", "hhsearch")

            hits_node = ET.SubElement(hh_run, "hits")
            for i, evidence in enumerate(self.hhsearch_hits):
                hit_elem = self._evidence_to_hhsearch_xml(evidence, i + 1)
                hits_node.append(hit_elem)

        return root

    def _evidence_to_blast_xml(self, evidence: Evidence, hit_num: int) -> ET.Element:
        """Convert Evidence object to BLAST hit XML format"""
        hit_elem = ET.Element("hit")
        hit_elem.set("num", str(hit_num))

        # Extract attributes from Evidence
        if evidence.domain_id:
            hit_elem.set("domain_id", evidence.domain_id)

        # Get PDB/chain from extra_attributes
        pdb_id = evidence.extra_attributes.get("pdb_id", "")
        chain_id = evidence.extra_attributes.get("chain_id", "")
        if pdb_id:
            hit_elem.set("pdb_id", pdb_id)
        if chain_id:
            hit_elem.set("chain_id", chain_id)

        # Set e-value and HSP count
        if evidence.evalue is not None:
            hit_elem.set("evalues", str(evidence.evalue))
        if evidence.hsp_count is not None:
            hit_elem.set("hsp_count", str(evidence.hsp_count))

        # Add ranges
        if evidence.query_range:
            query_reg = ET.SubElement(hit_elem, "query_reg")
            query_reg.text = evidence.query_range

        if evidence.hit_range:
            hit_reg = ET.SubElement(hit_elem, "hit_reg")
            hit_reg.text = evidence.hit_range

        return hit_elem

    def _evidence_to_hhsearch_xml(self, evidence: Evidence, hit_num: int) -> ET.Element:
        """Convert Evidence object to HHSearch hit XML format"""
        hit_elem = ET.Element("hit")

        # Set attributes
        hit_elem.set("num", str(hit_num))
        if evidence.domain_id:
            hit_elem.set("domain_id", evidence.domain_id)
        hit_elem.set("hit_id", evidence.source_id or evidence.domain_id or "")

        # Set HHSearch-specific attributes
        if evidence.probability is not None:
            hit_elem.set("probability", str(evidence.probability))
        if evidence.evalue is not None:
            hit_elem.set("evalue", str(evidence.evalue))
        if evidence.score is not None:
            hit_elem.set("score", str(evidence.score))

        # Add ranges
        if evidence.query_range:
            query_reg = ET.SubElement(hit_elem, "query_reg")
            query_reg.text = evidence.query_range

        if evidence.hit_range:
            hit_reg = ET.SubElement(hit_elem, "hit_reg")
            hit_reg.text = evidence.hit_range

        return hit_elem



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
