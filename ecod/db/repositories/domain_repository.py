# ecod/db/repositories/domain_repository.py
#!/usr/bin/env python3
"""
Domain repository for the ECOD pipeline
Handles database operations for domains and their classifications
"""
import logging
from typing import List, Dict, Any, Optional, Tuple, Set

# Use absolute imports
from ecod.exceptions import DatabaseError, QueryError, ValidationError
from ecod.db.manager import DBManager

# These models might need to be created or adjusted based on your actual models
try:
    from ecod.models.pipeline.domain import DomainModel as Domain
    from ecod.models.pipeline.domain import DomainSequence, DomainDSSPDetail
except ImportError:
    # Fallback if models don't exist yet
    Domain = dict
    DomainSequence = dict
    DomainDSSPDetail = dict

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
        """Get domain by ID"""
        query = """
        SELECT d.*
        FROM ecod_schema.domain d
        WHERE d.id = %s
        """
        rows = self.db.execute_dict_query(query, (domain_id,))
        if not rows:
            return None

        # Return dict if Domain model not available
        if Domain == dict:
            return rows[0]

        domain = Domain.from_db_row(rows[0])
        return domain

    def get_by_protein_id(self, protein_id: int) -> List[Domain]:
        """Get domains by protein ID"""
        query = """
        SELECT d.*
        FROM ecod_schema.domain d
        WHERE d.protein_id = %s
        ORDER BY d.id
        """
        rows = self.db.execute_dict_query(query, (protein_id,))

        domains = []
        for row in rows:
            if Domain == dict:
                domains.append(row)
            else:
                domain = Domain.from_db_row(row)
                domains.append(domain)

        return domains

    def create(self, domain: Domain) -> int:
        """Create a new domain"""
        # Handle both dict and object types
        if isinstance(domain, dict):
            domain_data = domain
        else:
            domain_data = {
                "ecod_uid": getattr(domain, 'ecod_uid', None),
                "protein_id": getattr(domain, 'protein_id', None),
                "domain_id": getattr(domain, 'domain_id', ''),
                "ecod_domain_id": getattr(domain, 'ecod_domain_id', None),
                "range": getattr(domain, 'range', ''),
                "t_group": getattr(domain, 't_group', None),
                "h_group": getattr(domain, 'h_group', None),
                "x_group": getattr(domain, 'x_group', None),
                "a_group": getattr(domain, 'a_group', None),
                "is_manual_rep": getattr(domain, 'is_manual_rep', False),
                "is_f70": getattr(domain, 'is_f70', False),
                "is_f40": getattr(domain, 'is_f40', False),
                "is_f99": getattr(domain, 'is_f99', False),
                "length": getattr(domain, 'length', None),
            }

        domain_id = self.db.insert(
            "ecod_schema.domain",
            domain_data,
            "id"
        )

        self.logger.info(f"Created domain {domain_data.get('domain_id', 'unknown')} with ID {domain_id}")
        return domain_id
