# ecod/models/pipeline/partition.py
"""
Consolidated DomainPartitionResult for pyECOD

This module consolidates and enhances the DomainPartitionResult class to work
seamlessly with the new Evidence and DomainModel implementations while maintaining
full backward compatibility with existing code.
"""

from dataclasses import dataclass, field
from datetime import datetime
from typing import List, Dict, Any, Optional, Union
import xml.etree.ElementTree as ET
import os
import logging
from pathlib import Path

from ecod.models.base import XmlSerializable


@dataclass
class DomainPartitionResult(XmlSerializable):
    """
    Enhanced domain partition result model for pipeline processing.

    Consolidates the functionality from the existing DomainPartitionResult while
    integrating seamlessly with the new Evidence and DomainModel implementations.
    """

    # Core identification
    pdb_id: str
    chain_id: str
    reference: str

    # Domain results (flexible to handle both models and dictionaries)
    domains: List[Union['DomainModel', Dict[str, Any]]] = field(default_factory=list)

    # Processing status
    success: bool = True
    error: Optional[str] = None

    # Classification status flags
    is_classified: bool = False
    is_unclassified: bool = False
    is_peptide: bool = False

    # Sequence and coverage metadata
    sequence_length: int = 0
    coverage: float = 0.0
    residues_assigned: int = 0
    residues_unassigned: int = 0

    # File management
    domain_file: Optional[str] = None

    # Processing metadata
    timestamp: Optional[datetime] = None
    processing_time: Optional[float] = None

    # Enhanced analysis metrics
    domain_quality_stats: Dict[str, Any] = field(default_factory=dict)
    evidence_summary: Dict[str, Any] = field(default_factory=dict)

    # Serialization control
    include_evidence: bool = True
    include_metadata: bool = True
    include_detailed_stats: bool = False

    def __post_init__(self):
        """Post-initialization processing and validation"""
        # Set timestamp if not provided
        if not self.timestamp:
            self.timestamp = datetime.now()

        # Standardize domains to use new models when available
        self._standardize_domains()

        # Update classification flags based on domains
        self._update_classification_status()

        # Calculate coverage and statistics
        if self.coverage == 0.0 and self.sequence_length > 0:
            self.calculate_coverage()

        # Generate analysis statistics
        self._update_analysis_stats()

    def _standardize_domains(self):
        """Ensure all domains use the new DomainModel when available"""
        try:
            from ecod.models.pipeline.domain import DomainModel
            NEW_DOMAIN_MODEL = True
        except ImportError:
            NEW_DOMAIN_MODEL = False
            return  # Keep domains as-is if new model not available

        standardized_domains = []
        for i, domain in enumerate(self.domains):
            if isinstance(domain, dict):
                try:
                    # Convert dictionary to DomainModel
                    domain_id = domain.get("id", domain.get("domain_id",
                                          f"{self.pdb_id}_{self.chain_id}_d{i+1}"))
                    domain_model = DomainModel.from_dict(domain, domain_id)
                    standardized_domains.append(domain_model)
                except Exception as e:
                    logging.getLogger(__name__).warning(
                        f"Error converting domain {i} to DomainModel: {str(e)}"
                    )
                    # Keep original dictionary if conversion fails
                    standardized_domains.append(domain)
            else:
                # Already a DomainModel or other object type
                standardized_domains.append(domain)

        self.domains = standardized_domains

    def _update_classification_status(self):
        """Update classification status based on domains"""
        # Handle peptides - they are successfully classified (as peptides) to prevent reprocessing
        if self.is_peptide:
            self.is_classified = True     # Successfully processed/identified as peptide
            self.is_unclassified = False  # Not unclassified - we know what it is
            return

        # For non-peptides, set is_classified to True if we have domains
        if len(self.domains) > 0:
            self.is_classified = True

    def _update_analysis_stats(self):
        """Update detailed analysis statistics"""
        if not self.domains:
            return

        # Domain quality statistics
        confidences = []
        evidence_counts = []
        fully_classified = 0
        protected_count = 0
        source_types = {}

        for domain in self.domains:
            # Extract confidence
            if hasattr(domain, 'confidence'):
                confidences.append(domain.confidence)
            elif isinstance(domain, dict):
                confidences.append(domain.get('confidence', 0.0))

            # Count evidence
            evidence_list = []
            if hasattr(domain, 'evidence'):
                evidence_list = domain.evidence
            elif isinstance(domain, dict):
                evidence_list = domain.get('evidence', [])

            evidence_counts.append(len(evidence_list))

            # Check classification completeness
            if hasattr(domain, 'is_fully_classified') and domain.is_fully_classified():
                fully_classified += 1
            elif isinstance(domain, dict):
                if all(domain.get(cls) for cls in ['t_group', 'h_group', 'x_group', 'a_group']):
                    fully_classified += 1

            # Check if protected
            if hasattr(domain, 'protected') and domain.protected:
                protected_count += 1
            elif isinstance(domain, dict) and domain.get('protected'):
                protected_count += 1

            # Count source types
            source = ""
            if hasattr(domain, 'source'):
                source = domain.source
            elif isinstance(domain, dict):
                source = domain.get('source', 'unknown')

            source_types[source] = source_types.get(source, 0) + 1

        # Calculate statistics
        self.domain_quality_stats = {
            "total_domains": len(self.domains),
            "fully_classified": fully_classified,
            "protected_domains": protected_count,
            "average_confidence": sum(confidences) / len(confidences) if confidences else 0.0,
            "min_confidence": min(confidences) if confidences else 0.0,
            "max_confidence": max(confidences) if confidences else 0.0,
            "average_evidence_count": sum(evidence_counts) / len(evidence_counts) if evidence_counts else 0.0,
            "source_distribution": source_types
        }

        # Evidence summary across all domains
        evidence_types = {}
        total_evidence = 0

        for domain in self.domains:
            evidence_list = []
            if hasattr(domain, 'evidence'):
                evidence_list = domain.evidence
            elif isinstance(domain, dict):
                evidence_list = domain.get('evidence', [])

            for evidence in evidence_list:
                total_evidence += 1
                ev_type = ""
                if hasattr(evidence, 'type'):
                    ev_type = evidence.type
                elif isinstance(evidence, dict):
                    ev_type = evidence.get('type', 'unknown')

                evidence_types[ev_type] = evidence_types.get(ev_type, 0) + 1

        self.evidence_summary = {
            "total_evidence_items": total_evidence,
            "evidence_type_distribution": evidence_types,
            "domains_with_evidence": sum(1 for d in self.domains
                                       if (hasattr(d, 'evidence') and d.evidence) or
                                          (isinstance(d, dict) and d.get('evidence')))
        }

    def calculate_coverage(self) -> None:
        """Calculate sequence coverage from domains"""
        if self.sequence_length == 0:
            self.coverage = 0.0
            self.residues_assigned = 0
            self.residues_unassigned = 0
            return

        # Collect all positions covered by domains
        covered_positions = set()

        for domain in self.domains:
            # Get domain range
            domain_range = ""
            if hasattr(domain, 'range'):
                domain_range = domain.range
            elif isinstance(domain, dict):
                domain_range = domain.get('range', '')

            # Parse range and add positions
            for segment in domain_range.split(","):
                if "-" in segment:
                    try:
                        start, end = map(int, segment.split("-"))
                        covered_positions.update(range(start, end + 1))
                    except ValueError:
                        pass  # Skip invalid ranges

        # Calculate coverage statistics
        self.residues_assigned = len(covered_positions)
        self.residues_unassigned = self.sequence_length - self.residues_assigned
        self.coverage = self.residues_assigned / self.sequence_length if self.sequence_length > 0 else 0.0

    def add_domain(self, domain: Union['DomainModel', Dict[str, Any]]) -> None:
        """Add a domain to the result"""
        self.domains.append(domain)

        # Update status and statistics
        self._update_classification_status()
        self.calculate_coverage()
        self._update_analysis_stats()

    def get_domain_by_id(self, domain_id: str) -> Optional[Union['DomainModel', Dict[str, Any]]]:
        """Get domain by ID"""
        for domain in self.domains:
            domain_id_field = ""
            if hasattr(domain, 'id'):
                domain_id_field = domain.id
            elif isinstance(domain, dict):
                domain_id_field = domain.get('id', domain.get('domain_id', ''))

            if domain_id_field == domain_id:
                return domain

        return None

    def get_domains_by_source(self, source: str) -> List[Union['DomainModel', Dict[str, Any]]]:
        """Get all domains from a specific source"""
        matching_domains = []

        for domain in self.domains:
            domain_source = ""
            if hasattr(domain, 'source'):
                domain_source = domain.source
            elif isinstance(domain, dict):
                domain_source = domain.get('source', '')

            if domain_source == source:
                matching_domains.append(domain)

        return matching_domains

    def get_overlapping_domains(self) -> List[tuple]:
        """Find pairs of overlapping domains"""
        overlaps = []

        try:
            from ecod.models.pipeline.domain import DomainModel
            NEW_DOMAIN_MODEL = True
        except ImportError:
            NEW_DOMAIN_MODEL = False

        for i, domain1 in enumerate(self.domains):
            for j, domain2 in enumerate(self.domains[i+1:], i+1):
                # Check overlap using DomainModel methods if available
                if NEW_DOMAIN_MODEL and hasattr(domain1, 'overlaps') and hasattr(domain2, 'overlaps'):
                    if domain1.overlaps(domain2):
                        overlaps.append((i, j, domain1, domain2))
                else:
                    # Fallback to manual overlap detection
                    start1 = getattr(domain1, 'start', domain1.get('start', 0) if isinstance(domain1, dict) else 0)
                    end1 = getattr(domain1, 'end', domain1.get('end', 0) if isinstance(domain1, dict) else 0)
                    start2 = getattr(domain2, 'start', domain2.get('start', 0) if isinstance(domain2, dict) else 0)
                    end2 = getattr(domain2, 'end', domain2.get('end', 0) if isinstance(domain2, dict) else 0)

                    if max(start1, start2) <= min(end1, end2):
                        overlaps.append((i, j, domain1, domain2))

        return overlaps

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation for backward compatibility"""
        result = {
            "pdb_id": self.pdb_id,
            "chain_id": self.chain_id,
            "reference": self.reference,
            "success": self.success,
            "error": self.error,
            "is_classified": self.is_classified,
            "is_unclassified": self.is_unclassified,
            "is_peptide": self.is_peptide,
            "sequence_length": self.sequence_length,
            "coverage": self.coverage,
            "residues_assigned": self.residues_assigned,
            "residues_unassigned": self.residues_unassigned,
            "domain_file": self.domain_file,
            "timestamp": self.timestamp.isoformat() if self.timestamp else None,
            "domain_count": len(self.domains)
        }

        # Convert domains to dictionaries
        domain_dicts = []
        for domain in self.domains:
            if hasattr(domain, 'to_dict'):
                domain_dicts.append(domain.to_dict())
            elif isinstance(domain, dict):
                domain_dicts.append(domain.copy())
            else:
                # Try to convert object to dictionary
                try:
                    domain_dicts.append(vars(domain))
                except:
                    domain_dicts.append({"error": "Could not serialize domain"})

        result["domains"] = domain_dicts

        # Add optional statistics if requested
        if self.include_detailed_stats:
            result["domain_quality_stats"] = self.domain_quality_stats
            result["evidence_summary"] = self.evidence_summary

        return result

    @classmethod
    def from_dict(cls, result_dict: Dict[str, Any]) -> 'DomainPartitionResult':
        """Create DomainPartitionResult from dictionary"""
        # Extract basic fields
        result = cls(
            pdb_id=result_dict.get("pdb_id", ""),
            chain_id=result_dict.get("chain_id", ""),
            reference=result_dict.get("reference", ""),
            success=result_dict.get("success", True),
            error=result_dict.get("error"),
            is_classified=result_dict.get("is_classified", False),
            is_unclassified=result_dict.get("is_unclassified", False),
            is_peptide=result_dict.get("is_peptide", False),
            sequence_length=result_dict.get("sequence_length", 0),
            coverage=result_dict.get("coverage", 0.0),
            residues_assigned=result_dict.get("residues_assigned", 0),
            residues_unassigned=result_dict.get("residues_unassigned", 0),
            domain_file=result_dict.get("domain_file"),
            include_detailed_stats=result_dict.get("include_detailed_stats", False)
        )

        # Handle timestamp
        if "timestamp" in result_dict and result_dict["timestamp"]:
            try:
                result.timestamp = datetime.fromisoformat(result_dict["timestamp"])
            except (ValueError, TypeError):
                pass

        # Add domains
        if "domains" in result_dict:
            result.domains = result_dict["domains"]  # Will be standardized in __post_init__

        return result

    @classmethod
    def from_domains(cls, pdb_id: str, chain_id: str, reference: str,
                    domains: List[Union['DomainModel', Dict[str, Any]]],
                    sequence_length: int = 0, **kwargs) -> 'DomainPartitionResult':
        """
        Create DomainPartitionResult from a list of domains.

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            reference: Reference version
            domains: List of domains (DomainModel objects or dictionaries)
            sequence_length: Sequence length
            **kwargs: Additional parameters

        Returns:
            DomainPartitionResult instance
        """
        result = cls(
            pdb_id=pdb_id,
            chain_id=chain_id,
            reference=reference,
            domains=domains,
            sequence_length=sequence_length,
            **kwargs
        )

        # Status is determined in __post_init__ based on domains
        return result

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

        # Add success/error information
        root.set("success", str(self.success).lower())
        if self.error:
            root.set("error", self.error)

        # Add metadata section if requested
        if self.include_metadata:
            metadata = ET.SubElement(root, "metadata")
            ET.SubElement(metadata, "sequence_length").text = str(self.sequence_length)
            ET.SubElement(metadata, "coverage").text = f"{self.coverage:.4f}"
            ET.SubElement(metadata, "residues_assigned").text = str(self.residues_assigned)
            ET.SubElement(metadata, "residues_unassigned").text = str(self.residues_unassigned)
            ET.SubElement(metadata, "domain_count").text = str(len(self.domains))

            if self.timestamp:
                ET.SubElement(metadata, "timestamp").text = self.timestamp.isoformat()

            if self.processing_time:
                ET.SubElement(metadata, "processing_time").text = f"{self.processing_time:.3f}"

            # Add detailed statistics if enabled
            if self.include_detailed_stats and self.domain_quality_stats:
                stats_elem = ET.SubElement(metadata, "quality_stats")
                for key, value in self.domain_quality_stats.items():
                    if isinstance(value, dict):
                        # Handle nested dictionaries (like source_distribution)
                        nested_elem = ET.SubElement(stats_elem, key)
                        for nested_key, nested_value in value.items():
                            nested_elem.set(nested_key, str(nested_value))
                    else:
                        stats_elem.set(key, str(value))

                if self.evidence_summary:
                    evidence_elem = ET.SubElement(metadata, "evidence_summary")
                    for key, value in self.evidence_summary.items():
                        if isinstance(value, dict):
                            nested_elem = ET.SubElement(evidence_elem, key)
                            for nested_key, nested_value in value.items():
                                nested_elem.set(nested_key, str(nested_value))
                        else:
                            evidence_elem.set(key, str(value))

        # Add domains section
        domains_elem = ET.SubElement(root, "domains")
        domains_elem.set("count", str(len(self.domains)))

        for domain in self.domains:
            # Use domain's XML serialization if available
            if hasattr(domain, 'to_xml'):
                domain_xml = domain.to_xml()

                # Control evidence inclusion
                if not self.include_evidence:
                    evidence_list = domain_xml.find("evidence_list")
                    if evidence_list is not None:
                        domain_xml.remove(evidence_list)

                domains_elem.append(domain_xml)
            else:
                # Manual XML creation for dictionary domains
                domain_elem = ET.SubElement(domains_elem, "domain")

                if isinstance(domain, dict):
                    # Set basic attributes
                    for attr in ["id", "start", "end", "range", "source", "source_id"]:
                        if attr in domain and domain[attr] is not None:
                            domain_elem.set(attr, str(domain[attr]))

                    # Set classification
                    for cls_attr in ["t_group", "h_group", "x_group", "a_group"]:
                        if cls_attr in domain and domain[cls_attr]:
                            domain_elem.set(cls_attr, domain[cls_attr])

                    # Set confidence
                    if "confidence" in domain:
                        domain_elem.set("confidence", f"{domain['confidence']:.4f}")

                    # Set flags
                    for flag in ["is_manual_rep", "is_f70", "is_f40", "is_f99", "protected"]:
                        if flag in domain:
                            domain_elem.set(flag, str(domain[flag]).lower())

                    # Add evidence if requested
                    if self.include_evidence and "evidence" in domain and domain["evidence"]:
                        evidence_list = ET.SubElement(domain_elem, "evidence_list")
                        evidence_list.set("count", str(len(domain["evidence"])))

                        for ev in domain["evidence"]:
                            if hasattr(ev, 'to_xml'):
                                evidence_list.append(ev.to_xml())
                            else:
                                # Manual evidence XML creation
                                ev_elem = ET.SubElement(evidence_list, "evidence")
                                if isinstance(ev, dict):
                                    for key, value in ev.items():
                                        if value is not None:
                                            ev_elem.set(key, str(value))

        return root

    @classmethod
    def from_xml(cls, element: ET.Element) -> 'DomainPartitionResult':
        """Create from XML Element"""
        # Get basic attributes
        pdb_id = element.get("pdb_id", "")
        chain_id = element.get("chain_id", "")
        reference = element.get("reference", "")
        success = element.get("success", "true").lower() == "true"
        error = element.get("error")
        is_classified = element.get("is_classified", "false").lower() == "true"
        is_unclassified = element.get("is_unclassified", "false").lower() == "true"
        is_peptide = element.get("is_peptide", "false").lower() == "true"

        # Create result
        result = cls(
            pdb_id=pdb_id,
            chain_id=chain_id,
            reference=reference,
            success=success,
            error=error,
            is_classified=is_classified,
            is_unclassified=is_unclassified,
            is_peptide=is_peptide
        )

        # Get metadata
        metadata = element.find("metadata")
        if metadata is not None:
            # Extract metadata values
            for field in ["sequence_length", "coverage", "residues_assigned",
                         "residues_unassigned", "processing_time"]:
                elem = metadata.find(field)
                if elem is not None and elem.text:
                    try:
                        if field in ["sequence_length", "residues_assigned", "residues_unassigned"]:
                            setattr(result, field, int(elem.text))
                        else:
                            setattr(result, field, float(elem.text))
                    except ValueError:
                        pass

            # Extract timestamp
            timestamp_elem = metadata.find("timestamp")
            if timestamp_elem is not None and timestamp_elem.text:
                try:
                    result.timestamp = datetime.fromisoformat(timestamp_elem.text)
                except ValueError:
                    pass

        # Get domains
        domains_elem = element.find("domains")
        if domains_elem is not None:
            try:
                from ecod.models.pipeline.domain import DomainModel
                NEW_DOMAIN_MODEL = True
            except ImportError:
                NEW_DOMAIN_MODEL = False

            for domain_elem in domains_elem.findall("domain"):
                if NEW_DOMAIN_MODEL:
                    try:
                        domain = DomainModel.from_xml(domain_elem)
                        result.domains.append(domain)
                    except Exception as e:
                        logging.getLogger(__name__).warning(
                            f"Error parsing domain from XML: {str(e)}"
                        )
                        # Fall back to dictionary
                        domain_dict = dict(domain_elem.attrib)
                        result.domains.append(domain_dict)
                else:
                    # Use dictionary approach
                    domain_dict = dict(domain_elem.attrib)

                    # Convert numeric fields
                    for field in ["start", "end", "confidence"]:
                        if field in domain_dict:
                            try:
                                if field in ["start", "end"]:
                                    domain_dict[field] = int(domain_dict[field])
                                else:
                                    domain_dict[field] = float(domain_dict[field])
                            except ValueError:
                                pass

                    # Convert boolean fields
                    for field in ["is_manual_rep", "is_f70", "is_f40", "is_f99", "protected"]:
                        if field in domain_dict:
                            domain_dict[field] = domain_dict[field].lower() == "true"

                    result.domains.append(domain_dict)

        return result

    def save(self, output_dir: str = None, filename: str = None) -> bool:
        """
        Save domain partition result to XML file.

        Args:
            output_dir: Output directory (if not using existing domain_file path)
            filename: Custom filename (if not using default naming)

        Returns:
            True if saved successfully, False otherwise
        """
        try:
            # Determine output path
            if filename and output_dir:
                output_path = os.path.join(output_dir, filename)
            elif self.domain_file:
                output_path = self.domain_file
            elif output_dir:
                # Generate default filename
                default_filename = f"{self.pdb_id}_{self.chain_id}.{self.reference}.domains.xml"
                domains_dir = os.path.join(output_dir, "domains")

                # Create domains directory with error handling
                try:
                    os.makedirs(domains_dir, exist_ok=True)
                except (PermissionError, OSError) as e:
                    logging.getLogger(__name__).error(f"Cannot create domains directory {domains_dir}: {str(e)}")
                    return False

                output_path = os.path.join(domains_dir, default_filename)
            else:
                logging.getLogger(__name__).error("No output path specified for saving")
                return False

            # Update domain_file path
            self.domain_file = output_path

            # Create parent directories with error handling
            try:
                parent_dir = os.path.dirname(output_path)
                if parent_dir:  # Only create if there's actually a parent directory
                    os.makedirs(parent_dir, exist_ok=True)
            except (PermissionError, OSError) as e:
                logging.getLogger(__name__).error(f"Cannot create parent directory for {output_path}: {str(e)}")
                return False

            # Save to XML file
            try:
                xml_element = self.to_xml()
                tree = ET.ElementTree(xml_element)
                tree.write(output_path, encoding='utf-8', xml_declaration=True)

                logging.getLogger(__name__).info(f"Saved domain partition: {output_path}")
                return True

            except (PermissionError, OSError, IOError) as e:
                logging.getLogger(__name__).error(f"Error writing domain partition file {output_path}: {str(e)}")
                return False

        except Exception as e:
            # Catch any other unexpected errors
            logging.getLogger(__name__).error(f"Unexpected error saving domain partition: {str(e)}")
            return False

    def get_summary_stats(self) -> Dict[str, Any]:
        """Get summary statistics for reporting"""
        return {
            "pdb_id": self.pdb_id,
            "chain_id": self.chain_id,
            "reference": self.reference,
            "success": self.success,
            "domain_count": len(self.domains),
            "is_classified": self.is_classified,
            "is_peptide": self.is_peptide,
            "sequence_length": self.sequence_length,
            "coverage": self.coverage,
            "residues_assigned": self.residues_assigned,
            "processing_time": self.processing_time,
            "fully_classified_domains": self.domain_quality_stats.get("fully_classified", 0),
            "average_confidence": self.domain_quality_stats.get("average_confidence", 0.0),
            "total_evidence_items": self.evidence_summary.get("total_evidence_items", 0)
        }

    def __str__(self) -> str:
        """String representation"""
        status = "SUCCESS" if self.success else "FAILED"
        if self.is_peptide:
            classification = "PEPTIDE"
        elif self.is_classified:
            classification = f"CLASSIFIED ({len(self.domains)} domains)"
        else:
            classification = "UNCLASSIFIED"

        return f"DomainPartitionResult({self.pdb_id}_{self.chain_id}, {status}, {classification})"

    def __repr__(self) -> str:
        """Detailed string representation"""
        return (f"DomainPartitionResult(pdb_id='{self.pdb_id}', chain_id='{self.chain_id}', "
                f"reference='{self.reference}', domains={len(self.domains)}, "
                f"success={self.success}, coverage={self.coverage:.3f})")


# Backward compatibility alias
PartitionResult = DomainPartitionResult
