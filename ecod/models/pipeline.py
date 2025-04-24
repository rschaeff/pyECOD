# ecod/models/pipeline.py
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict, Any

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
    
    @classmethod
    def from_xml(cls, hit_elem):
        """Create from XML Element"""
        hit = cls(
            hit_id=hit_elem.get("num", ""),
            domain_id=hit_elem.get("domain_id", ""),
            pdb_id=hit_elem.get("pdb_id", ""),
            chain_id=hit_elem.get("chain_id", ""),
            evalue=float(hit_elem.get("evalues", "999").split(",")[0]),
            hsp_count=int(hit_elem.get("hsp_count", "0")),
            discontinuous=hit_elem.get("discontinuous", "false").lower() == "true"
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

@dataclass
class DomainSummary:
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
    
    def to_xml(self):
        """Convert to XML Element"""
        import xml.etree.ElementTree as ET
        
        root = ET.Element("blast_summ_doc")
        
        # Create summary node
        blast_summ = ET.SubElement(root, "blast_summ")
        blast_summ.set("pdb", self.pdb_id)
        blast_summ.set("chain", self.chain_id)
        
        # Add errors
        for error, value in self.errors.items():
            if value:
                blast_summ.set(error, "true")
        
        # Add length
        if self.sequence_length > 0:
            query_len = ET.SubElement(root, "query_len")
            query_len.text = str(self.sequence_length)
        
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
                hit_elem.set("evalues", ",".join(str(e) for e in hit.evalues) if hit.evalues else str(hit.evalue))
                
                query_reg = ET.SubElement(hit_elem, "query_reg")
                query_reg.text = hit.range
                
                hit_reg = ET.SubElement(hit_elem, "hit_reg")
                hit_reg.text = hit.hit_range
        
        # Add domain blast hits
        if self.domain_blast_hits:
            domain_blast_run = ET.SubElement(root, "blast_run")
            domain_blast_run.set("program", "blastp")
            
            hits_node = ET.SubElement(domain_blast_run, "hits")
            for hit in self.domain_blast_hits:
                hit_elem = ET.SubElement(hits_node, "hit")
                hit_elem.set("domain_id", hit.domain_id)
                hit_elem.set("pdb_id", hit.pdb_id)
                hit_elem.set("chain_id", hit.chain_id)
                hit_elem.set("hsp_count", str(hit.hsp_count))
                hit_elem.set("evalues", ",".join(str(e) for e in hit.evalues) if hit.evalues else str(hit.evalue))
                
                if hit.discontinuous:
                    hit_elem.set("discontinuous", "true")
                
                query_reg = ET.SubElement(hit_elem, "query_reg")
                query_reg.text = hit.range
                
                hit_reg = ET.SubElement(hit_elem, "hit_reg")
                hit_reg.text = hit.hit_range
        
        # Add HHSearch hits
        if self.hhsearch_hits:
            hh_run = ET.SubElement(root, "hh_run")
            hh_run.set("program", "hhsearch")
            
            hits_node = ET.SubElement(hh_run, "hits")
            for i, hit in enumerate(self.hhsearch_hits):
                hit_elem = ET.SubElement(hits_node, "hit")
                hit_elem.set("domain_id", hit.domain_id)
                hit_elem.set("hit_id", hit.hit_id)
                hit_elem.set("num", str(i+1))
                hit_elem.set("probability", str(hit.probability))
                hit_elem.set("evalue", str(hit.evalue))
                hit_elem.set("score", str(hit.score))
                
                query_reg = ET.SubElement(hit_elem, "query_reg")
                query_reg.text = hit.range
                
                hit_reg = ET.SubElement(hit_elem, "hit_reg")
                hit_reg.text = hit.hit_range
        
        return root