# ecod/utils/evidence_bridge.py
"""
Bridge module to facilitate consistent evidence handling during the transition
from dictionary-based to model-based domain analysis.
"""

import logging
from typing import List, Dict, Any, Optional, Union, TypeVar, Tuple
import xml.etree.ElementTree as ET

from ecod.models.domain_analysis.domain_model import DomainModel
from ecod.models.domain_analysis.evidence import Evidence
from ecod.models.evidence import DomainEvidence

# Type variable for domains
DomainType = TypeVar('DomainType', Dict[str, Any], DomainModel)
# Type variable for evidence
EvidenceType = TypeVar('EvidenceType', Dict[str, Any], Evidence, DomainEvidence)


class EvidenceBridge:
    """
    Utilities for converting and standardizing evidence handling between
    dictionary-based and model-based approaches.
    """
    
    @staticmethod
    def standardize_evidence(evidence: EvidenceType) -> DomainEvidence:
        """
        Convert any evidence representation to a standardized DomainEvidence object.
        
        Args:
            evidence: Evidence in any supported format (dict, Evidence, DomainEvidence)
            
        Returns:
            Standardized DomainEvidence object
        """
        logger = logging.getLogger(__name__)
        
        # If already a DomainEvidence, return as is
        if isinstance(evidence, DomainEvidence):
            return evidence
            
        # If it's an Evidence object from domain_analysis module
        if hasattr(evidence, 'type') and hasattr(evidence, 'to_xml'):
            # Convert domain_analysis.Evidence to DomainEvidence
            try:
                return DomainEvidence(
                    type=evidence.type,
                    source_id=getattr(evidence, 'source_id', ''),
                    domain_id=getattr(evidence, 'domain_id', ''),
                    query_range=getattr(evidence, 'query_range', ''),
                    hit_range=getattr(evidence, 'hit_range', ''),
                    confidence=getattr(evidence, 'confidence', 0.0),
                    t_group=getattr(evidence, 't_group', None),
                    h_group=getattr(evidence, 'h_group', None),
                    x_group=getattr(evidence, 'x_group', None),
                    a_group=getattr(evidence, 'a_group', None),
                    attributes={k: v for k, v in evidence.__dict__.items() 
                               if not k.startswith('_') and k not in [
                                   'type', 'source_id', 'domain_id', 'query_range', 
                                   'hit_range', 'confidence', 't_group', 'h_group', 
                                   'x_group', 'a_group'
                               ]}
                )
            except Exception as e:
                logger.warning(f"Error converting Evidence to DomainEvidence: {str(e)}")
                # Fall through to dictionary conversion as backup
        
        # If it's a dictionary
        if isinstance(evidence, dict):
            return DomainEvidence.from_dict(evidence)
            
        # Last resort - try to extract basic information from unknown object
        try:
            evidence_dict = {}
            for attr in ['type', 'source_id', 'domain_id', 'query_range', 'hit_range', 
                       'confidence', 't_group', 'h_group', 'x_group', 'a_group']:
                if hasattr(evidence, attr):
                    evidence_dict[attr] = getattr(evidence, attr)
                    
            # Add remaining attributes to attributes dict
            attributes = {}
            for attr in dir(evidence):
                if not attr.startswith('_') and attr not in evidence_dict and not callable(getattr(evidence, attr)):
                    attributes[attr] = getattr(evidence, attr)
            evidence_dict['attributes'] = attributes
            
            return DomainEvidence.from_dict(evidence_dict)
        except Exception as e:
            logger.error(f"Failed to convert evidence of type {type(evidence)}: {str(e)}")
            # Return a minimal evidence object
            return DomainEvidence(
                type="unknown",
                source_id=str(evidence) if evidence else ""
            )

    @staticmethod
    def ensure_domain_model(domain: DomainType, pdb_id: str = "", chain_id: str = "", 
                          index: int = 0) -> DomainModel:
        """
        Ensure a domain is represented as a DomainModel, converting if necessary.
        
        Args:
            domain: Domain in dictionary or DomainModel format
            pdb_id: PDB identifier (for creating domain ID if needed)
            chain_id: Chain identifier (for creating domain ID if needed)
            index: Domain index (for creating domain ID if needed)
            
        Returns:
            DomainModel representation
        """
        logger = logging.getLogger(__name__)
        
        # If already a DomainModel, return as is
        if isinstance(domain, DomainModel):
            return domain
            
        # If it's a dictionary, convert to DomainModel
        if isinstance(domain, dict):
            # Create domain ID if not present
            domain_id = domain.get('id', f"{pdb_id}_{chain_id}_d{index+1}")
            
            # Get basic fields
            start = domain.get('start', 0)
            end = domain.get('end', 0)
            range_text = domain.get('range', f"{start}-{end}")
            
            # Convert evidence list
            evidence_list = []
            for evidence in domain.get('evidence', []):
                try:
                    standardized_evidence = EvidenceBridge.standardize_evidence(evidence)
                    evidence_list.append(standardized_evidence)
                except Exception as e:
                    logger.warning(f"Error standardizing evidence: {str(e)}")
            
            # Create DomainModel
            return DomainModel(
                id=domain_id,
                start=start,
                end=end,
                range=range_text,
                t_group=domain.get('t_group'),
                h_group=domain.get('h_group'),
                x_group=domain.get('x_group'),
                a_group=domain.get('a_group'),
                source=domain.get('source', ''),
                confidence=domain.get('confidence', 0.0),
                source_id=domain.get('source_id', ''),
                is_manual_rep=domain.get('is_manual_rep', False),
                is_f70=domain.get('is_f70', False),
                is_f40=domain.get('is_f40', False),
                is_f99=domain.get('is_f99', False),
                evidence=evidence_list
            )
        
        # Unsupported type
        logger.error(f"Unsupported domain type: {type(domain)}")
        raise TypeError(f"Unsupported domain type: {type(domain)}")

    @staticmethod
    def ensure_evidence_list(evidence_items: List[Any]) -> List[DomainEvidence]:
        """
        Ensure a list of evidence items are all represented as DomainEvidence objects.
        
        Args:
            evidence_items: List of evidence items in any supported format
            
        Returns:
            List of standardized DomainEvidence objects
        """
        return [EvidenceBridge.standardize_evidence(item) for item in evidence_items]
    
    @staticmethod
    def domain_list_to_models(domains: List[Dict[str, Any]], pdb_id: str, 
                            chain_id: str) -> List[DomainModel]:
        """
        Convert a list of domain dictionaries to DomainModel objects.
        
        Args:
            domains: List of domain dictionaries
            pdb_id: PDB identifier
            chain_id: Chain identifier
            
        Returns:
            List of DomainModel objects
        """
        return [EvidenceBridge.ensure_domain_model(domain, pdb_id, chain_id, i) 
                for i, domain in enumerate(domains)]
    
    @staticmethod
    def add_evidence_to_domain(domain: DomainType, evidence: EvidenceType) -> DomainType:
        """
        Add evidence to a domain, handling both dictionary and model formats.
        
        Args:
            domain: Domain in dictionary or model format
            evidence: Evidence in any supported format
            
        Returns:
            Domain with evidence added (same type as input domain)
        """
        standardized_evidence = EvidenceBridge.standardize_evidence(evidence)
        
        if isinstance(domain, DomainModel):
            domain.evidence.append(standardized_evidence)
            return domain
        elif isinstance(domain, dict):
            if 'evidence' not in domain:
                domain['evidence'] = []
            domain['evidence'].append(standardized_evidence)
            return domain
        else:
            raise TypeError(f"Unsupported domain type: {type(domain)}")
    
    @staticmethod
    def update_classification_from_evidence(domain: DomainType) -> DomainType:
        """
        Update domain classification fields from its evidence, if available.
        
        Args:
            domain: Domain in dictionary or model format
            
        Returns:
            Domain with updated classification (same type as input domain)
        """
        # Get evidence list based on domain type
        evidence_list = domain.evidence if isinstance(domain, DomainModel) else domain.get('evidence', [])
        
        # Find best evidence for classification
        best_evidence = None
        best_score = 0
        
        for evidence in evidence_list:
            # Standardize evidence
            std_evidence = EvidenceBridge.standardize_evidence(evidence)
            
            # Skip evidence without classification
            if not (std_evidence.t_group or std_evidence.h_group or std_evidence.x_group or std_evidence.a_group):
                continue
            
            # Calculate score based on confidence and completeness
            confidence = std_evidence.confidence
            completeness = sum(1 for group in [std_evidence.t_group, std_evidence.h_group, 
                                             std_evidence.x_group, std_evidence.a_group] if group)
            score = confidence * (completeness / 4.0)
            
            if score > best_score:
                best_score = score
                best_evidence = std_evidence
        
        # Update domain classification if best evidence found
        if best_evidence:
            # For DomainModel
            if isinstance(domain, DomainModel):
                if best_evidence.t_group:
                    domain.t_group = best_evidence.t_group
                if best_evidence.h_group:
                    domain.h_group = best_evidence.h_group
                if best_evidence.x_group:
                    domain.x_group = best_evidence.x_group
                if best_evidence.a_group:
                    domain.a_group = best_evidence.a_group
            # For dictionary
            elif isinstance(domain, dict):
                if best_evidence.t_group:
                    domain['t_group'] = best_evidence.t_group
                if best_evidence.h_group:
                    domain['h_group'] = best_evidence.h_group
                if best_evidence.x_group:
                    domain['x_group'] = best_evidence.x_group
                if best_evidence.a_group:
                    domain['a_group'] = best_evidence.a_group
        
        return domain
