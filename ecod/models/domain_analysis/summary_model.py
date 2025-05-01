# ecod/models/domain_analysis/summary_model.py
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional
import xml.etree.ElementTree as ET
from ecod.models.base import XmlSerializable
from ecod.models.pipeline import BlastHit, HHSearchHit

@dataclass
class EnhancedDomainSummaryModel(XmlSerializable):
    """Enhanced domain summary model"""
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
    
    # Additional fields for enhanced functionality
    domain_suggestions: List[Dict[str, Any]] = field(default_factory=list)
    classification_sources: Dict[str, Any] = field(default_factory=dict)
    processing_metadata: Dict[str, Any] = field(default_factory=dict)
    
    def to_xml(self) -> ET.Element:
        """Convert to XML Element"""
        # Implementation similar to current DomainSummaryModel.to_xml()
        # But with enhancements for new fields
        
        # ... implementation ...
    
    @classmethod
    def from_xml(cls, element: ET.Element) -> 'EnhancedDomainSummaryModel':
        """Create from XML Element"""
        # Implementation similar to current DomainSummaryModel parsing
        # But with enhancements for new fields
        
        # ... implementation ...
    
    @classmethod
    def from_legacy_model(cls, legacy_model):
        """Convert from legacy DomainSummaryModel"""
        # Create new model with data from legacy model
        
        # ... implementation ...
