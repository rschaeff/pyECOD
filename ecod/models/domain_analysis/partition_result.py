# ecod/models/domain_analysis/partition_result.py
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Set
import xml.etree.ElementTree as ET
from pathlib import Path
from ecod.models.base import XmlSerializable
from ecod.models.domain_analysis.domain_candidate import DomainCandidate

@dataclass
class DomainPartitionResult(XmlSerializable):
    """Model for domain partition results"""
    pdb_id: str
    chain_id: str
    reference: str
    domains: List[Dict[str, Any]] = field(default_factory=list)
    success: bool = True
    error: Optional[str] = None
    domain_file: Optional[str] = None
    is_classified: bool = False
    is_unclassified: bool = False
    is_peptide: bool = False
    
    # Classification information
    classification_info: Dict[str, Any] = field(default_factory=dict)
    
    # Stats and metadata
    processing_metadata: Dict[str, Any] = field(default_factory=dict)
    
    def get_path(self) -> Optional[str]:
        """Get path to domain partition file if available"""
        return self.domain_file
    
    def to_xml(self) -> ET.Element:
        """Convert to XML Element"""
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
        
        # Add domains
        domains_elem = ET.SubElement(root, "domains")
        for domain in self.domains:
            dom_elem = ET.SubElement(domains_elem, "domain")
            # Add domain attributes
            for key, value in domain.items():
                if value is not None and not isinstance(value, (dict, list)):
                    dom_elem.set(str(key), str(value))
            
            # Add classification if available
            if "classification" in domain and domain["classification"]:
                cls_elem = ET.SubElement(dom_elem, "classification")
                for key, value in domain["classification"].items():
                    if value is not None:
                        cls_elem.set(str(key), str(value))
        
        return root
    
    @classmethod
    def from_xml(cls, element: ET.Element) -> 'DomainPartitionResult':
        """Create from XML Element"""
        # Extract basic attributes
        result = cls(
            pdb_id=element.get("pdb_id", ""),
            chain_id=element.get("chain_id", ""),
            reference=element.get("reference", ""),
            is_classified=element.get("is_classified", "false").lower() == "true",
            is_unclassified=element.get("is_unclassified", "false").lower() == "true",
            is_peptide=element.get("is_peptide", "false").lower() == "true"
        )
        
        # Extract domains
        domains_elem = element.find("domains")
        if domains_elem is not None:
            for dom_elem in domains_elem.findall("domain"):
                domain = {key: value for key, value in dom_elem.attrib.items()}
                
                # Extract classification if available
                cls_elem = dom_elem.find("classification")
                if cls_elem is not None:
                    domain["classification"] = {
                        key: value for key, value in cls_elem.attrib.items()
                    }
                
                result.domains.append(domain)
        
        return result
    
    @classmethod
    def from_file(cls, file_path: str) -> Optional['DomainPartitionResult']:
        """Create instance from XML file"""
        try:
            result = cls.from_xml_file(file_path)
            if result:
                result.domain_file = file_path
            return result
        except Exception as e:
            import logging
            logging.getLogger(__name__).error(f"Error loading domain partition: {str(e)}")
            return None
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary
        
        Returns:
            Dictionary representation
        """
        return {
            "pdb_id": self.pdb_id,
            "chain_id": self.chain_id,
            "reference": self.reference,
            "success": self.success,
            "error": self.error,
            "domain_file": self.domain_file,
            "is_classified": self.is_classified,
            "is_unclassified": self.is_unclassified,
            "is_peptide": self.is_peptide,
            "domains": self.domains,
            "classification_info": self.classification_info,
            "processing_metadata": self.processing_metadata
        }
    
    def save(self, output_dir: str = None) -> bool:
        """Save to file
        
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
    
    def __str__(self) -> str:
        """String representation
        
        Returns:
            String representation
        """
        domain_count = len(self.domains)
        domain_str = f"{domain_count} domain{'s' if domain_count != 1 else ''}"
        status = "classified" if self.is_classified else "unclassified"
        return f"{self.pdb_id}_{self.chain_id} ({status}): {domain_str}"
