# ecod/db/repositories/protein_repository.py
#!/usr/bin/env python3
"""
Protein repository for the ECOD pipeline
Handles database operations for proteins and sequences
"""
import logging
from typing import List, Dict, Any, Optional, Tuple

# Use absolute imports - these should work with your project structure
from ecod.exceptions import DatabaseError, QueryError, ValidationError
from ecod.db.manager import DBManager
from ecod.models.protein import Protein, ProteinSequence, ProteinStructure

class ProteinRepository:
    """Repository for protein data"""

    def __init__(self, db_manager: DBManager):
        """Initialize repository

        Args:
            db_manager: Database manager instance
        """
        self.db = db_manager
        self.logger = logging.getLogger("ecod.db.protein_repository")

    def get_by_id(self, protein_id: int) -> Optional[Protein]:
        """Get protein by ID

        Args:
            protein_id: Protein ID

        Returns:
            Protein if found, None otherwise
        """
        query = """
        SELECT p.*, ps.sequence, ps.sequence_md5
        FROM ecod_schema.protein p
        LEFT JOIN ecod_schema.protein_sequence ps ON p.id = ps.protein_id
        WHERE p.id = %s
        """
        rows = self.db.execute_dict_query(query, (protein_id,))
        if not rows:
            return None
        return Protein.from_db_row(rows[0], include_sequence=True)

    def get_by_pdb_chain(self, pdb_id: str, chain_id: str) -> Optional[Protein]:
        """Get protein by PDB ID and chain ID"""
        query = """
        SELECT p.*, ps.sequence, ps.sequence_md5
        FROM ecod_schema.protein p
        LEFT JOIN ecod_schema.protein_sequence ps ON p.id = ps.protein_id
        WHERE p.pdb_id = %s AND p.chain_id = %s
        """
        rows = self.db.execute_dict_query(query, (pdb_id, chain_id))
        if not rows:
            return None
        return Protein.from_db_row(rows[0], include_sequence=True)

    def create(self, protein: Protein) -> int:
        """Create a new protein"""
        # Basic implementation for testing
        protein_id = self.db.insert(
            "ecod_schema.protein",
            {
                "pdb_id": protein.pdb_id,
                "chain_id": protein.chain_id,
                "source_id": protein.source_id,
                "length": protein.length,
                "unp_acc": getattr(protein, 'unp_acc', None),
                "name": getattr(protein, 'name', None),
                "type": getattr(protein, 'type', None),
                "tax_id": getattr(protein, 'tax_id', None)
            },
            "id"
        )

        # Insert sequence if provided
        if hasattr(protein, 'sequence') and protein.sequence:
            self.db.insert(
                "ecod_schema.protein_sequence",
                {
                    "protein_id": protein_id,
                    "sequence": protein.sequence.sequence,
                    "sequence_md5": protein.sequence.sequence_md5 if hasattr(protein.sequence, 'sequence_md5') else None,
                    "fragment_start": getattr(protein.sequence, 'fragment_start', None),
                    "fragment_end": getattr(protein.sequence, 'fragment_end', None)
                }
            )

        self.logger.info(f"Created protein {protein.source_id} with ID {protein_id}")
        return protein_id
