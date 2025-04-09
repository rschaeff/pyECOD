#!/usr/bin/env python3
"""
Domain classifier module for ECOD pipeline
Responsible for domain classification and hierarchy assignment
"""

import logging
from typing import Dict, Any, List, Optional, Set, Tuple

from ecod.db.manager import DBManager
from ecod.core.exceptions import ClassificationError


class DomainClassifier:
    """
    Handles domain classification and hierarchy assignment
    
    This class is responsible for retrieving and caching domain classifications
    from the database, and assigning them to new domains based on evidence.
    """
    
    def __init__(self, db_manager: DBManager):
        """
        Initialize with database connection
        
        Args:
            db_manager: Database manager instance
        """
        self.db = db_manager
        self.logger = logging.getLogger("ecod.pipelines.domain_analysis.classifier")
        
        # Initialize classification caches
        self._uid_classification_cache = {}
        self._domain_id_classification_cache = {}
        self._reference_uid_lookup = {}
        
    def get_classification_by_uid(self, ecod_uid: int) -> Optional[Dict[str, Any]]:
        """
        Get domain classification by ECOD UID with caching
        
        Args:
            ecod_uid: ECOD unique identifier
            
        Returns:
            Dictionary with classification data or None if not found
            
        Raises:
            ClassificationError: If database error occurs
        """
        # Check cache first
        if ecod_uid in self._uid_classification_cache:
            return self._uid_classification_cache[ecod_uid]
            
        # Query database
        query = """
        SELECT 
            d.t_group, d.h_group, d.x_group, d.a_group,
            d.is_manual_rep, d.is_f70, d.is_f40, d.is_f99
        FROM 
            pdb_analysis.domain d
        WHERE 
            d.ecod_uid = %s
        """
        
        try:
            rows = self.db.execute_dict_query(query, (ecod_uid,))
            if rows:
                classification = {
                    "t_group": rows[0].get("t_group"),
                    "h_group": rows[0].get("h_group"),
                    "x_group": rows[0].get("x_group"),
                    "a_group": rows[0].get("a_group"),
                    "is_manual_rep": rows[0].get("is_manual_rep", False),
                    "is_f70": rows[0].get("is_f70", False),
                    "is_f40": rows[0].get("is_f40", False),
                    "is_f99": rows[0].get("is_f99", False)
                }
                
                # Cache the result
                self._uid_classification_cache[ecod_uid] = classification
                return classification
                
        except Exception as e:
            self.logger.error(f"Error getting domain classification for UID {ecod_uid}: {e}")
            raise ClassificationError(f"Failed to retrieve classification for UID {ecod_uid}")
            
        return None

    def get_classification_by_domain_id(self, domain_id: str) -> Optional[Dict[str, Any]]:
        """
        Get domain classification from database by domain ID with caching
        
        Args:
            domain_id: ECOD domain identifier (e.g., d2xyzB1)
            
        Returns:
            Dictionary with classification data or None if not found
            
        Raises:
            ClassificationError: If database error occurs
        """
        # Check cache first
        if domain_id in self._domain_id_classification_cache:
            return self._domain_id_classification_cache[domain_id]
        
        # Try to get UID from lookup table and use that
        if domain_id in self._reference_uid_lookup:
            uid = self._reference_uid_lookup[domain_id]
            classification = self.get_classification_by_uid(uid)
            if classification:
                # Cache by domain ID too
                self._domain_id_classification_cache[domain_id] = classification
                return classification
            
        # Query database directly
        query = """
        SELECT 
            d.t_group, d.h_group, d.x_group, d.a_group,
            d.is_manual_rep, d.is_f70, d.is_f40, d.is_f99,
            d.ecod_uid
        FROM 
            pdb_analysis.domain d
        WHERE 
            d.domain_id = %s
        """
        
        try:
            rows = self.db.execute_dict_query(query, (domain_id,))
            if rows:
                classification = {
                    "t_group": rows[0].get("t_group"),
                    "h_group": rows[0].get("h_group"),
                    "x_group": rows[0].get("x_group"),
                    "a_group": rows[0].get("a_group"),
                    "is_manual_rep": False,  # New domains are not manual reps
                    "is_f70": False,
                    "is_f40": False,
                    "is_f99": False
                }
                
                # Cache the result
                self._domain_id_classification_cache[domain_id] = classification
                
                # Update UID lookup and cache by UID too
                uid = rows[0].get("ecod_uid")
                if uid:
                    self._reference_uid_lookup[domain_id] = uid
                    self._uid_classification_cache[uid] = classification
                    
                return classification
                
        except Exception as e:
            self.logger.error(f"Error getting classification for domain {domain_id}: {e}")
            raise ClassificationError(f"Failed to retrieve classification for domain {domain_id}")
            
        return None
    
    def assign_domain_classification(self, domain: Dict[str, Any], evidence: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Assign classification to a domain based on evidence
        
        Args:
            domain: Domain dictionary to classify
            evidence: List of evidence items
            
        Returns:
            Updated domain dictionary with classification
            
        Raises:
            ClassificationError: If classification fails
        """
        # Case 1: Domain already has reference classification
        if "reference" in domain and domain["reference"] and "uid" in domain:
            uid = domain["uid"]
            classification = self.get_classification_by_uid(uid)
            if classification:
                domain.update(classification)
                return domain
        
        # Case 2: Use best evidence to determine classification
        if not evidence:
            return domain
            
        # Find best evidence (highest probability or lowest e-value)
        best_evidence = None
        best_score = 0
        
        for item in evidence:
            domain_id = item.get("domain_id", "")
            if not domain_id or domain_id == "NA":
                continue
            
            # Calculate score based on evidence type
            if item["type"] == "hhsearch":
                score = item.get("probability", 0)
            elif item["type"] == "blast":
                e_value = item.get("evalue", 999)
                score = 100.0 / (1.0 + e_value) if e_value < 10 else 0
            else:
                score = 0
            
            if score > best_score:
                best_score = score
                best_evidence = item
        
        if best_evidence:
            domain_id = best_evidence.get("domain_id", "")
            
            # Get classification for this domain
            classification = self.get_classification_by_domain_id(domain_id)
            if classification:
                domain.update(classification)
                
                # Add evidence source to domain
                if "classification_source" not in domain:
                    domain["classification_source"] = {
                        "domain_id": domain_id,
                        "evidence_type": best_evidence["type"],
                        "score": best_score
                    }
        
        return domain
    
    def load_reference_domains(self) -> None:
        """
        Load reference domain data from database
        
        This preloads domain IDs, UIDs, and classifications for faster lookup
        
        Raises:
            ClassificationError: If loading fails
        """
        query = """
        SELECT 
            d.domain_id, d.ecod_uid, 
            d.t_group, d.h_group, d.x_group, d.a_group,
            d.is_manual_rep, d.is_f70, d.is_f40, d.is_f99
        FROM 
            pdb_analysis.domain d
        WHERE 
            d.is_f70 = TRUE
        LIMIT 10000
        """
        
        try:
            rows = self.db.execute_dict_query(query)
            
            for row in rows:
                domain_id = row.get("domain_id")
                uid = row.get("ecod_uid")
                
                if domain_id and uid:
                    self._reference_uid_lookup[domain_id] = uid
                    
                    # Cache classification
                    classification = {
                        "t_group": row.get("t_group"),
                        "h_group": row.get("h_group"),
                        "x_group": row.get("x_group"),
                        "a_group": row.get("a_group"),
                        "is_manual_rep": row.get("is_manual_rep", False),
                        "is_f70": row.get("is_f70", False),
                        "is_f40": row.get("is_f40", False),
                        "is_f99": row.get("is_f99", False)
                    }
                    
                    self._uid_classification_cache[uid] = classification
                    self._domain_id_classification_cache[domain_id] = classification
                    
            self.logger.info(f"Loaded {len(self._reference_uid_lookup)} reference domains")
                
        except Exception as e:
            self.logger.error(f"Error loading reference domains: {e}")
            raise ClassificationError("Failed to load reference domain data")
    
    def get_domain_statistics(self) -> Dict[str, Any]:
        """
        Get statistics about domain classifications
        
        Returns:
            Dictionary with statistics
        """
        stats = {
            "total_domains": len(self._domain_id_classification_cache),
            "toplogy_groups": set(),
            "homology_groups": set()
        }
        
        for domain_id, classification in self._domain_id_classification_cache.items():
            t_group = classification.get("t_group")
            h_group = classification.get("h_group")
            
            if t_group:
                stats["toplogy_groups"].add(t_group)
            
            if h_group:
                stats["homology_groups"].add(h_group)
        
        stats["toplogy_group_count"] = len(stats["toplogy_groups"])
        stats["homology_group_count"] = len(stats["homology_groups"])
        
        # Remove sets for serialization
        del stats["toplogy_groups"]
        del stats["homology_groups"]
        
        return stats