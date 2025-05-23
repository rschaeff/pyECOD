# models/pipeline/partition.py
from dataclasses import dataclass, field
from datetime import datetime
from typing import List, Dict, Any, Optional, Union
import xml.etree.ElementTree as ET
from ecod.models.base import XmlSerializable
from ecod.models.pipeline.domain import DomainModel

@dataclass
class DomainPartitionResult(XmlSerializable):
    """Pipeline output for domain partitioning results"""
    pdb_id: str
    chain_id: str
    reference: str
    domains: List[DomainModel] = field(default_factory=list)
    
    # Status flags
    success: bool = True
    error: Optional[str] = None
    is_classified: bool = False
    is_unclassified: bool = False
    is_peptide: bool = False
    
    # Metadata
    sequence_length: int = 0
    coverage: float = 0.0
    residues_assigned: int = 0
    timestamp: Optional[datetime] = None
    
    # File path
    domain_file: Optional[str] = None
    
    # Serialization control
    include_evidence: bool = True
    include_metadata: bool = True

    def __post_init__(self):
        """Initialize dynamic fields"""
        if not self.timestamp:
            self.timestamp = datetime.now()
        
        # Ensure domains are DomainModel objects
        self._standardize_domains()
        
        # Update classification status
        self.is_classified = len(self.domains) > 0
        self.is_unclassified = len(self.domains) == 0
        
        # Calculate coverage
        if self.coverage == 0.0:
            self.calculate_coverage()

    def _standardize_domains(self):
        """Convert any dictionary domains to DomainModel objects"""
        standardized = []
        for i, domain in enumerate(self.domains):
            if isinstance(domain, dict):
                # Convert dict to DomainModel
                domain_model = DomainModel(
                    id=domain.get("id", f"{self.pdb_id}_{self.chain_id}_d{i+1}"),
                    start=domain.get("start", 0),
                    end=domain.get("end", 0),
                    range=domain.get("range", ""),
                    t_group=domain.get("t_group"),
                    h_group=domain.get("h_group"),
                    x_group=domain.get("x_group"),
                    a_group=domain.get("a_group"),
                    source=domain.get("source", ""),
                    confidence=domain.get("confidence", 0.0),
                    source_id=domain.get("source_id", ""),
                    is_manual_rep=domain.get("is_manual_rep", False),
                    is_f70=domain.get("is_f70", False),
                    is_f40=domain.get("is_f40", False),
                    is_f99=domain.get("is_f99", False)
                )
                
                # Add evidence if present
                if "evidence" in domain and domain["evidence"]:
                    for ev_dict in domain["evidence"]:
                        from ecod.models.pipeline.evidence import Evidence
                        evidence = Evidence.from_dict(ev_dict) if isinstance(ev_dict, dict) else ev_dict
                        domain_model.evidence.append(evidence)
                
                standardized.append(domain_model)
            else:
                standardized.append(domain)
        
        self.domains = standardized

    def calculate_coverage(self) -> None:
        """Calculate sequence coverage from domains"""
        if self.sequence_length == 0:
            self.coverage = 0.0
            self.residues_assigned = 0
            return
        
        positions = set()
        for domain in self.domains:
            # Parse domain range to get positions
            domain_range = domain.range if hasattr(domain, 'range') else ""
            for segment in domain_range.split(","):
                if "-" in segment:
                    try:
                        start, end = map(int, segment.split("-"))
                        positions.update(range(start, end + 1))
                    except ValueError:
                        pass
        
        self.residues_assigned = len(positions)
        self.coverage = len(positions) / self.sequence_length

    def to_xml(self) -> ET.Element:
        """Convert to XML with configurable detail level"""
        root = ET.Element("domain_partition")
        root.set("pdb_id", self.pdb_id)
        root.set("chain_id", self.chain_id)
        root.set("reference", self.reference)
        
        # Status flags
        if self.is_classified:
            root.set("is_classified", "true")
        if self.is_unclassified:
            root.set("is_unclassified", "true")
        if self.is_peptide:
            root.set("is_peptide", "true")
        
        # Metadata section
        if self.include_metadata:
            metadata = ET.SubElement(root, "metadata")
            ET.SubElement(metadata, "sequence_length").text = str(self.sequence_length)
            ET.SubElement(metadata, "coverage").text = f"{self.coverage:.4f}"
            ET.SubElement(metadata, "residues_assigned").text = str(self.residues_assigned)
            ET.SubElement(metadata, "domain_count").text = str(len(self.domains))
            if self.timestamp:
                ET.SubElement(metadata, "timestamp").text = self.timestamp.isoformat()
        
        # Domains section
        domains_elem = ET.SubElement(root, "domains")
        for domain in self.domains:
            domain_xml = domain.to_xml()
            
            # Control evidence inclusion
            if not self.include_evidence:
                evidence_list = domain_xml.find("evidence_list")
                if evidence_list is not None:
                    domain_xml.remove(evidence_list)
            
            domains_elem.append(domain_xml)
        
        return root

    @classmethod
    def from_domains(cls, pdb_id: str, chain_id: str, reference: str,
                    domains: List[Union[Dict[str, Any], DomainModel]],
                    sequence_length: int = 0) -> 'DomainPartitionResult':
        """Create result from domain list"""
        return cls(
            pdb_id=pdb_id,
            chain_id=chain_id,
            reference=reference,
            domains=domains,
            sequence_length=sequence_length
        )
