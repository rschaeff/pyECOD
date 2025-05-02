# ecod/models/domain_analysis/partition_result.py
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Set
import xml.etree.ElementTree as ET
from pathlib import Path
from ecod.models.base import XmlSerializable
from ecod.models.domain_analysis.domain_candidate import DomainCandidate
from ecod.models.domain_analysis.domain_model import DomainModel

@dataclass
class DomainPartitionResult(XmlSerializable):
    """Enhanced model for domain partition results with complete information"""
    pdb_id: str
    chain_id: str
    reference: str
    domains: List[DomainModel] = field(default_factory=list)
    success: bool = True
    error: Optional[str] = None
    domain_file: Optional[str] = None
    is_classified: bool = False
    is_unclassified: bool = False
    is_peptide: bool = False

    # New metadata fields
    sequence_length: int = 0
    coverage: float = 0.0
    residues_assigned: int = 0
    timestamp: Optional[datetime] = None

    # Serialization control
    include_evidence: bool = True
    include_metadata: bool = True

    def __post_init__(self):
        """Initialize any dynamic fields"""
        if not self.timestamp:
            self.timestamp = datetime.now()

    def to_xml(self) -> ET.Element:
        """Convert to XML Element with configurable detail level"""
        root = ET.Element("domain_partition")
        root.set("pdb_id", self.pdb_id)
        root.set("chain_id", self.chain_id)
        root.set("reference", self.reference)

        # Set status attributes
        if self.is_classified:
            root.set("is_classified", "true")
        if self.is_unclassified:
            root.set("is_unclassified", "true")
        if self.is_peptide:
            root.set("is_peptide", "true")

        # Add metadata section if requested
        if self.include_metadata:
            metadata = ET.SubElement(root, "metadata")
            ET.SubElement(metadata, "sequence_length").text = str(self.sequence_length)
            ET.SubElement(metadata, "timestamp").text = self.timestamp.isoformat()
            ET.SubElement(metadata, "domain_count").text = str(len(self.domains))
            ET.SubElement(metadata, "coverage").text = f"{self.coverage:.4f}"
            ET.SubElement(metadata, "residues_assigned").text = str(self.residues_assigned)

            # Add domain statistics
            ET.SubElement(metadata, "domains_with_evidence").text = str(sum(1 for d in self.domains if d.evidence))
            ET.SubElement(metadata, "fully_classified_domains").text = str(sum(
                1 for d in self.domains if d.t_group and d.h_group and d.x_group and d.a_group
            ))

        # Add domains
        domains_elem = ET.SubElement(root, "domains")
        for domain in self.domains:
            # Convert DomainModel to XML with evidence if requested
            domain_xml = domain.to_xml()

            # Handle evidence inclusion/exclusion
            if not self.include_evidence:
                evidence_list = domain_xml.find("evidence_list")
                if evidence_list is not None:
                    domain_xml.remove(evidence_list)

            domains_elem.append(domain_xml)

        return root

    def save(self, output_dir: str = None) -> bool:
        """Save to file with appropriate serialization options

        Args:
            output_dir: Optional output directory (if not using existing path)

        Returns:
            True if saved successfully
        """
        if self.domain_file and not output_dir:
            return self.to_xml_file(self.domain_file)

        # Create path if output_dir is provided
        if output_dir:
            import os
            path = os.path.join(
                output_dir,
                "domains",
                f"{self.pdb_id}_{self.chain_id}.{self.reference}.domains.xml"
            )
            self.domain_file = path
            return self.to_xml_file(path)

        return False

    def calculate_coverage(self) -> None:
        """Calculate coverage from domains"""
        if self.sequence_length == 0:
            self.coverage = 0.0
            self.residues_assigned = 0
            return

        # Get all positions covered by domains
        positions = set()
        for domain in self.domains:
            for segment in domain.range.split(","):
                if "-" in segment:
                    try:
                        start, end = map(int, segment.split("-"))
                        positions.update(range(start, end + 1))
                    except ValueError:
                        pass

        self.residues_assigned = len(positions)
        self.coverage = len(positions) / self.sequence_length if self.sequence_length > 0 else 0.0
