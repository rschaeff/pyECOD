#!/usr/bin/env python3
"""
Domain repository for the ECOD pipeline
Handles database operations for domains and their classifications
"""
import logging
from typing import List, Dict, Any, Optional, Tuple, Set

from ....errors import DatabaseError, QueryError, ValidationError
from ...db import DBManager
from ...models.domain import Domain, DomainSequence, DomainDSSPDetail, DomainClassification

class DomainRepository:
    """Repository for domain data"""
    
    def __init__(self, db_manager: DBManager):
        """Initialize repository
        
        Args:
            db_manager: Database manager instance
        """
        self.db = db_manager
        self.logger = logging.getLogger("ecod.db.domain_repository")
        
    def get_by_id(self, domain_id: int) -> Optional[Domain]:
        """Get domain by ID
        
        Args:
            domain_id: Domain ID
            
        Returns:
            Domain if found, None otherwise
        """
        query = """
        SELECT d.*
        FROM ecod_schema.domain d
        WHERE d.id = %s
        """
        rows = self.db.execute_dict_query(query, (domain_id,))
        if not rows:
            return None
        
        domain = Domain.from_db_row(rows[0])
        
        # Add sequence if available
        domain.sequence = self.get_sequence(domain_id)
        
        # Add DSSP details if available
        domain.dssp = self.get_dssp_details(domain_id)
        
        return domain
        
    def get_by_ecod_uid(self, ecod_uid: int) -> Optional[Domain]:
        """Get domain by ECOD UID
        
        Args:
            ecod_uid: ECOD UID
            
        Returns:
            Domain if found, None otherwise
        """
        query = """
        SELECT d.*
        FROM ecod_schema.domain d
        WHERE d.ecod_uid = %s
        """
        rows = self.db.execute_dict_query(query, (ecod_uid,))
        if not rows:
            return None
        
        domain = Domain.from_db_row(rows[0])
        
        # Add sequence if available
        domain.sequence = self.get_sequence(domain.id)
        
        # Add DSSP details if available
        domain.dssp = self.get_dssp_details(domain.id)
        
        return domain
        
    def get_by_domain_id(self, domain_str_id: str) -> Optional[Domain]:
        """Get domain by string domain ID (e.g., "e1abcA1")
        
        Args:
            domain_str_id: Domain string ID
            
        Returns:
            Domain if found, None otherwise
        """
        query = """
        SELECT d.*
        FROM ecod_schema.domain d
        WHERE d.domain_id = %s
        """
        rows = self.db.execute_dict_query(query, (domain_str_id,))
        if not rows:
            return None
        
        domain = Domain.from_db_row(rows[0])
        
        # Add sequence if available
        domain.sequence = self.get_sequence(domain.id)
        
        # Add DSSP details if available
        domain.dssp = self.get_dssp_details(domain.id)
        
        return domain
    
    def get_by_protein_id(self, protein_id: int) -> List[Domain]:
        """Get domains by protein ID
        
        Args:
            protein_id: Protein ID
            
        Returns:
            List of domains
        """
        query = """
        SELECT d.*
        FROM ecod_schema.domain d
        WHERE d.protein_id = %s
        ORDER BY d.id
        """
        rows = self.db.execute_dict_query(query, (protein_id,))
        
        domains = []
        for row in rows:
            domain = Domain.from_db_row(row)
            domain.sequence = self.get_sequence(domain.id)
            domain.dssp = self.get_dssp_details(domain.id)
            domains.append(domain)
            
        return domains
    
    def get_by_classification(self, level: str, value: str) -> List[Domain]:
        """Get domains by classification level and value
        
        Args:
            level: Classification level ('t_group', 'h_group', 'x_group', 'a_group')
            value: Classification value
            
        Returns:
            List of domains
        """
        valid_levels = ['t_group', 'h_group', 'x_group', 'a_group']
        if level not in valid_levels:
            raise ValidationError(f"Invalid classification level: {level}. Must be one of {valid_levels}")
            
        query = f"""
        SELECT d.*
        FROM ecod_schema.domain d
        WHERE d.{level} = %s
        ORDER BY d.id
        """
        rows = self.db.execute_dict_query(query, (value,))
        
        domains = []
        for row in rows:
            domain = Domain.from_db_row(row)
            domains.append(domain)
            
        return domains
    
    def search(self, criteria: Dict[str, Any], limit: Optional[int] = None) -> List[Domain]:
        """Search domains by criteria
        
        Args:
            criteria: Search criteria (key-value pairs)
            limit: Maximum number of domains to return
            
        Returns:
            List of matching domains
        """
        # Build query conditions
        conditions = []
        params = []
        
        for key, value in criteria.items():
            if key == 'length_min':
                conditions.append("d.length >= %s")
                params.append(value)
            elif key == 'length_max':
                conditions.append("d.length <= %s")
                params.append(value)
            elif key == 'representative':
                if value:
                    conditions.append("d.is_manual_rep = TRUE")
                else:
                    conditions.append("d.is_manual_rep = FALSE")
            elif key == 'has_structure':
                if value:
                    conditions.append("EXISTS (SELECT 1 FROM ecod_schema.domain_structure ds WHERE ds.domain_id = d.id)")
                else:
                    conditions.append("NOT EXISTS (SELECT 1 FROM ecod_schema.domain_structure ds WHERE ds.domain_id = d.id)")
            else:
                # Default: exact match on domain field
                conditions.append(f"d.{key} = %s")
                params.append(value)
        
        # Build query
        query = """
        SELECT d.*
        FROM ecod_schema.domain d
        """
        
        if conditions:
            query += " WHERE " + " AND ".join(conditions)
            
        query += " ORDER BY d.id"
        
        if limit is not None:
            query += f" LIMIT {limit}"
            
        rows = self.db.execute_dict_query(query, tuple(params))
        
        domains = []
        for row in rows:
            domain = Domain.from_db_row(row)
            domains.append(domain)
            
        return domains
    
    def get_sequence(self, domain_id: int) -> Optional[DomainSequence]:
        """Get domain sequence
        
        Args:
            domain_id: Domain ID
            
        Returns:
            DomainSequence if found, None otherwise
        """
        query = """
        SELECT *
        FROM ecod_schema.domain_sequence
        WHERE domain_id = %s
        """
        rows = self.db.execute_dict_query(query, (domain_id,))
        if not rows:
            return None
        return DomainSequence.from_db_row(rows[0])
    
    def get_dssp_details(self, domain_id: int) -> Optional[DomainDSSPDetail]:
        """Get domain DSSP details
        
        Args:
            domain_id: Domain ID
            
        Returns:
            DomainDSSPDetail if found, None otherwise
        """
        query = """
        SELECT *
        FROM ecod_schema.domain_dssp_detail
        WHERE domain_id = %s
        """
        rows = self.db.execute_dict_query(query, (domain_id,))
        if not rows:
            return None
        return DomainDSSPDetail.from_db_row(rows[0])
    
    def create(self, domain: Domain) -> int:
        """Create a new domain
        
        Args:
            domain: Domain to create
            
        Returns:
            New domain ID
            
        Raises:
            ValidationError: If domain validation fails
        """
        # Validate domain
        domain.validate()
        
        # Check if domain already exists
        if domain.domain_id:
            existing = self.get_by_domain_id(domain.domain_id)
            if existing:
                raise ValidationError(f"Domain with domain_id {domain.domain_id} already exists",
                                    {"existing_id": existing.id})
        
        # Insert domain
        domain_id = self.db.insert(
            "ecod_schema.domain",
            {
                "ecod_uid": domain.ecod_uid,
                "protein_id": domain.protein_id,
                "domain_id": domain.domain_id,
                "ecod_domain_id": domain.ecod_domain_id,
                "range": domain.range,
                "t_group": domain.t_group,
                "h_group": domain.h_group,
                "x_group": domain.x_group,
                "a_group": domain.a_group,
                "is_manual_rep": domain.is_manual_rep,
                "is_f70": domain.is_f70,
                "is_f40": domain.is_f40,
                "is_f99": domain.is_f99,
                "hcount": domain.hcount,
                "scount": domain.scount,
                "length": domain.length,
                "chain_id": domain.chain_id,
                "asym_id": domain.asym_id
            },
            "id"
        )
        
        # Insert sequence if provided
        if domain.sequence:
            self.db.insert(
                "ecod_schema.domain_sequence",
                {
                    "domain_id": domain_id,
                    "sequence": domain.sequence.sequence,
                    "sequence_length": domain.sequence.sequence_length,
                    "sequence_md5": domain.sequence.sequence_md5,
                    "original_range": domain.sequence.original_range
                }
            )
            
        # Insert DSSP details if provided
        if domain.dssp:
            self.db.insert(
                "ecod_schema.domain_dssp_detail",
                {
                    "domain_id": domain_id,
                    "asa": domain.dssp.asa,
                    "helix_residues": domain.dssp.helix_residues,
                    "strand_residues": domain.dssp.strand_residues,
                    "helix_count": domain.dssp.helix_count,
                    "strand_count": domain.dssp.strand_count,
                    "secondary_structure_string": domain.dssp.secondary_structure_string,
                    "hbonds_total": domain.dssp.hbonds_total
                }
            )
            
        self.logger.info(f"Created domain {domain.domain_id} with ID {domain_id}")
        return domain_id
    
    def update(self, domain: Domain) -> bool:
        """Update an existing domain
        
        Args:
            domain: Domain to update
            
        Returns:
            True if update was successful
            
        Raises:
            ValidationError: If domain validation fails
        """
        # Validate domain
        domain.validate()
        
        # Ensure domain ID is set
        if domain.id is None:
            raise ValidationError("Domain ID is required for update")
        
        # Check if domain exists
        existing = self.get_by_id(domain.id)
        if not existing:
            raise ValidationError(f"Domain with ID {domain.id} does not exist")
        
        # Update domain
        rows_updated = self.db.update(
            "ecod_schema.domain",
            {
                "ecod_uid": domain.ecod_uid,
                "protein_id": domain.protein_id,
                "domain_id": domain.domain_id,
                "ecod_domain_id": domain.ecod_domain_id,
                "range": domain.range,
                "t_group": domain.t_group,
                "h_group": domain.h_group,
                "x_group": domain.x_group,
                "a_group": domain.a_group,
                "is_manual_rep": domain.is_manual_rep,
                "is_f70": domain.is_f70,
                "is_f40": domain.is_f40,
                "is_f99": domain.is_f99,
                "hcount": domain.hcount,
                "scount": domain.scount,
                "length": domain.length,
                "chain_id": domain.chain_id,
                "asym_id": domain.asym_id,
                "updated_at": "CURRENT_TIMESTAMP"
            },
            "id = %s",
            (domain.id,)
        )
        
        # Update sequence if provided
        if domain.sequence:
            # Check if sequence exists
            sequence = self.get_sequence(domain.id)
            
            if sequence:
                # Update existing sequence
                self.db.update(
                    "ecod_schema.domain_sequence",
                    {
                        "sequence": domain.sequence.sequence,
                        "sequence_length": domain.sequence.sequence_length,
                        "sequence_md5": domain.sequence.sequence_md5,
                        "original_range": domain.sequence.original_range
                    },
                    "domain_id = %s",
                    (domain.id,)
                )
            else:
                # Insert new sequence
                self.db.insert(
                    "ecod_schema.domain_sequence",
                    {
                        "domain_id": domain.id,
                        "sequence": domain.sequence.sequence,
                        "sequence_length": domain.sequence.sequence_length,
                        "sequence_md5": domain.sequence.sequence_md5,
                        "original_range": domain.sequence.original_range
                    }
                )
                
        # Update DSSP details if provided
        if domain.dssp:
            # Check if DSSP details exist
            dssp = self.get_dssp_details(domain.id)
            
            if dssp:
                # Update existing DSSP details
                self.db.update(
                    "ecod_schema.domain_dssp_detail",
                    {
                        "asa": domain.dssp.asa,
                        "helix_residues": domain.dssp.helix_residues,
                        "strand_residues": domain.dssp.strand_residues,
                        "helix_count": domain.dssp.helix_count,
                        "strand_count": domain.dssp.strand_count,
                        "secondary_structure_string": domain.dssp.secondary_structure_string,
                        "hbonds_total": domain.dssp.hbonds_total
                    },
                    "domain_id = %s",
                    (domain.id,)
                )
            else:
                # Insert new DSSP details
                self.db.insert(
                    "ecod_schema.domain_dssp_detail",
                    {
                        "domain_id": domain.id,
                        "asa