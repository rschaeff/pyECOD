# ecod/models/domain_analysis/partition_result.py - Enhanced Version
from dataclasses import dataclass, field
from datetime import datetime
from typing import List, Dict, Any, Optional, Set, Union
import xml.etree.ElementTree as ET
from pathlib import Path
import logging
import os

from ecod.models.base import XmlSerializable
from ecod.models.domain_analysis.domain_model import DomainModel
from ecod.utils.evidence_bridge import EvidenceBridge  # Import the bridge

@dataclass
class DomainPartitionResult(XmlSerializable):
    """Enhanced model for domain partition results with complete information"""
    pdb_id: str
    chain_id: str
    reference: str
    domains: List[Union[DomainModel, Dict[str, Any]]] = field(default_factory=list)
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
        """Initialize any dynamic fields and standardize domains"""
        if not self.timestamp:
            self.timestamp = datetime.now()

        # Standardize domains to ensure they are all DomainModel objects
        self._standardize_domains()

        # Set classification flags
        self.is_classified = len(self.domains) > 0
        self.is_unclassified = len(self.domains) == 0

        # Calculate coverage if not already set
        if self.coverage == 0.0:
            self.calculate_coverage()

    def _standardize_domains(self):
        """Ensure all domains are DomainModel objects"""
        standardized_domains = []
        for i, domain in enumerate(self.domains):
            if isinstance(domain, dict):
                try:
                    # Use the EvidenceBridge to convert dictionaries to DomainModel
                    domain_model = EvidenceBridge.ensure_domain_model(
                        domain, self.pdb_id, self.chain_id, i
                    )
                    standardized_domains.append(domain_model)
                except Exception as e:
                    logging.getLogger(__name__).warning(
                        f"Error converting domain to DomainModel: {str(e)}"
                    )
                    # Keep the original dictionary in this case
                    standardized_domains.append(domain)
            else:
                # Already a DomainModel or other object
                standardized_domains.append(domain)

        # Replace with standardized domains
        self.domains = standardized_domains

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
            ET.SubElement(metadata, "domains_with_evidence").text = str(
                sum(1 for d in self.domains
                    if (isinstance(d, DomainModel) and d.evidence) or
                       (isinstance(d, dict) and d.get('evidence')))
            )
            ET.SubElement(metadata, "fully_classified_domains").text = str(
                sum(1 for d in self.domains
                    if (isinstance(d, DomainModel) and
                        d.t_group and d.h_group and d.x_group and d.a_group) or
                       (isinstance(d, dict) and
                        d.get('t_group') and d.get('h_group') and
                        d.get('x_group') and d.get('a_group')))
            )

        # Add domains
        domains_elem = ET.SubElement(root, "domains")
        for domain in self.domains:
            # Handle DomainModel and dictionary domains
            if isinstance(domain, DomainModel):
                # Use the DomainModel's to_xml method
                domain_xml = domain.to_xml()

                # Handle evidence inclusion/exclusion
                if not self.include_evidence:
                    evidence_list = domain_xml.find("evidence_list")
                    if evidence_list is not None:
                        domain_xml.remove(evidence_list)

                domains_elem.append(domain_xml)
            else:
                # For dictionary domains, create element manually
                domain_elem = ET.SubElement(domains_elem, "domain")

                # Extract basic info
                for attr in ["id", "start", "end", "range", "source", "source_id",
                           "t_group", "h_group", "x_group", "a_group"]:
                    if attr in domain and domain[attr] is not None:
                        domain_elem.set(attr, str(domain[attr]))

                # Add confidence if available
                if "confidence" in domain:
                    domain_elem.set("confidence", f"{domain['confidence']:.4f}")

                # Add representative/filter flags
                for flag in ["is_manual_rep", "is_f70", "is_f40", "is_f99"]:
                    domain_elem.set(flag, str(domain.get(flag, False)).lower())

                # Add evidence if available and requested
                if self.include_evidence and "evidence" in domain and domain["evidence"]:
                    evidence_list = ET.SubElement(domain_elem, "evidence_list")
                    for ev in domain["evidence"]:
                        # Use EvidenceBridge to standardize evidence
                        std_ev = EvidenceBridge.standardize_evidence(ev)
                        evidence_list.append(std_ev.to_xml())

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
            path = os.path.join(
                output_dir,
                "domains",
                f"{self.pdb_id}_{self.chain_id}.{self.reference}.domains.xml"
            )
            # Create directory if it doesn't exist
            os.makedirs(os.path.dirname(path), exist_ok=True)

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
            if isinstance(domain, DomainModel):
                domain_range = domain.range
            else:
                domain_range = domain.get('range', '')

            for segment in domain_range.split(","):
                if "-" in segment:
                    try:
                        start, end = map(int, segment.split("-"))
                        positions.update(range(start, end + 1))
                    except ValueError:
                        pass

        self.residues_assigned = len(positions)
        self.coverage = len(positions) / self.sequence_length if self.sequence_length > 0 else 0.0

    @classmethod
    def from_domains(cls, pdb_id: str, chain_id: str, reference: str,
                    domains: List[Union[Dict[str, Any], DomainModel]],
                    sequence_length: int = 0) -> 'DomainPartitionResult':
        """
        Create DomainPartitionResult from a list of domains

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            reference: Reference version
            domains: List of domains (dictionary or DomainModel)
            sequence_length: Sequence length

        Returns:
            DomainPartitionResult instance
        """
        result = cls(
            pdb_id=pdb_id,
            chain_id=chain_id,
            reference=reference,
            domains=domains,
            sequence_length=sequence_length,
            is_classified=len(domains) > 0,
            is_unclassified=len(domains) == 0
        )

        # Calculate coverage
        result.calculate_coverage()

        return result

    @classmethod
    def from_xml(cls, element: ET.Element) -> 'DomainPartitionResult':
        """Create from XML Element

        Args:
            element: XML Element

        Returns:
            DomainPartitionResult instance
        """
        # Get basic attributes
        pdb_id = element.get("pdb_id", "")
        chain_id = element.get("chain_id", "")
        reference = element.get("reference", "")
        is_classified = element.get("is_classified", "false").lower() == "true"
        is_unclassified = element.get("is_unclassified", "false").lower() == "true"
        is_peptide = element.get("is_peptide", "false").lower() == "true"

        # Create result
        result = cls(
            pdb_id=pdb_id,
            chain_id=chain_id,
            reference=reference,
            is_classified=is_classified,
            is_unclassified=is_unclassified,
            is_peptide=is_peptide
        )

        # Get metadata
        metadata = element.find("metadata")
        if metadata is not None:
            sequence_length_elem = metadata.find("sequence_length")
            if sequence_length_elem is not None and sequence_length_elem.text:
                result.sequence_length = int(sequence_length_elem.text)

            timestamp_elem = metadata.find("timestamp")
            if timestamp_elem is not None and timestamp_elem.text:
                try:
                    result.timestamp = datetime.fromisoformat(timestamp_elem.text)
                except ValueError:
                    pass

            coverage_elem = metadata.find("coverage")
            if coverage_elem is not None and coverage_elem.text:
                try:
                    result.coverage = float(coverage_elem.text)
                except ValueError:
                    pass

            residues_assigned_elem = metadata.find("residues_assigned")
            if residues_assigned_elem is not None and residues_assigned_elem.text:
                try:
                    result.residues_assigned = int(residues_assigned_elem.text)
                except ValueError:
                    pass

        # Get domains
        domains_elem = element.find("domains")
        if domains_elem is not None:
            for domain_elem in domains_elem.findall("domain"):
                try:
                    domain = DomainModel.from_xml(domain_elem)
                    result.domains.append(domain)
                except Exception as e:
                    logging.getLogger(__name__).warning(
                        f"Error parsing domain from XML: {str(e)}"
                    )

        return result
