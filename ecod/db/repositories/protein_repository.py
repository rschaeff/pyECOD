#!/usr/bin/env python3
"""
Protein repository for the ECOD pipeline
Handles database operations for proteins and sequences
"""
import logging
from typing import List, Dict, Any, Optional, Tuple

from ...core.errors import DatabaseError, QueryError, ValidationError
from ...db.manager import DBManager
from ...models.protein import Protein, ProteinSequence, ProteinStructure

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
        """Get protein by PDB ID and chain ID
        
        Args:
            pdb_id: PDB ID
            chain_id: Chain ID
            
        Returns:
            Protein if found, None otherwise
        """
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
    
    def get_by_source_id(self, source_id: str) -> Optional[Protein]:
        """Get protein by source ID
        
        Args:
            source_id: Source ID (e.g., "1abc_A")
            
        Returns:
            Protein if found, None otherwise
        """
        query = """
        SELECT p.*, ps.sequence, ps.sequence_md5
        FROM ecod_schema.protein p
        LEFT JOIN ecod_schema.protein_sequence ps ON p.id = ps.protein_id
        WHERE p.source_id = %s
        """
        rows = self.db.execute_dict_query(query, (source_id,))
        if not rows:
            return None
        return Protein.from_db_row(rows[0], include_sequence=True)
    
    def get_all(self, limit: Optional[int] = None, offset: int = 0) -> List[Protein]:
        """Get all proteins
        
        Args:
            limit: Maximum number of proteins to return
            offset: Offset for pagination
            
        Returns:
            List of proteins
        """
        query = """
        SELECT p.*, ps.sequence, ps.sequence_md5
        FROM ecod_schema.protein p
        LEFT JOIN ecod_schema.protein_sequence ps ON p.id = ps.protein_id
        ORDER BY p.id
        """
        
        if limit is not None:
            query += f" LIMIT {limit} OFFSET {offset}"
            
        rows = self.db.execute_dict_query(query)
        return [Protein.from_db_row(row, include_sequence=True) for row in rows]
    
    def search(self, criteria: Dict[str, Any], limit: Optional[int] = None) -> List[Protein]:
        """Search proteins by criteria
        
        Args:
            criteria: Search criteria (key-value pairs)
            limit: Maximum number of proteins to return
            
        Returns:
            List of matching proteins
        """
        # Build query conditions
        conditions = []
        params = []
        
        for key, value in criteria.items():
            if key == 'sequence_contains':
                conditions.append("ps.sequence LIKE %s")
                params.append(f"%{value}%")
            elif key == 'length_min':
                conditions.append("p.length >= %s")
                params.append(value)
            elif key == 'length_max':
                conditions.append("p.length <= %s")
                params.append(value)
            elif key == 'has_domains':
                if value:
                    conditions.append("EXISTS (SELECT 1 FROM ecod_schema.domain d WHERE d.protein_id = p.id)")
                else:
                    conditions.append("NOT EXISTS (SELECT 1 FROM ecod_schema.domain d WHERE d.protein_id = p.id)")
            else:
                # Default: exact match on protein field
                conditions.append(f"p.{key} = %s")
                params.append(value)
        
        # Build query
        query = """
        SELECT p.*, ps.sequence, ps.sequence_md5
        FROM ecod_schema.protein p
        LEFT JOIN ecod_schema.protein_sequence ps ON p.id = ps.protein_id
        """
        
        if conditions:
            query += " WHERE " + " AND ".join(conditions)
            
        query += " ORDER BY p.id"
        
        if limit is not None:
            query += f" LIMIT {limit}"
            
        rows = self.db.execute_dict_query(query, tuple(params))
        return [Protein.from_db_row(row, include_sequence=True) for row in rows]
    
    def create(self, protein: Protein) -> int:
        """Create a new protein
        
        Args:
            protein: Protein to create
            
        Returns:
            New protein ID
            
        Raises:
            ValidationError: If protein validation fails
        """
        # Validate protein
        protein.validate()
        
        # Check if protein already exists
        existing = self.get_by_source_id(protein.source_id)
        if existing:
            raise ValidationError(f"Protein with source ID {protein.source_id} already exists",
                                 {"existing_id": existing.id})
        
        # Insert protein
        protein_id = self.db.insert(
            "ecod_schema.protein",
            {
                "pdb_id": protein.pdb_id,
                "chain_id": protein.chain_id,
                "source_id": protein.source_id,
                "length": protein.length,
                "unp_acc": protein.unp_acc,
                "name": protein.name,
                "type": protein.type,
                "tax_id": protein.tax_id
            },
            "id"
        )
        
        # Insert sequence if provided
        if protein.sequence:
            self.db.insert(
                "ecod_schema.protein_sequence",
                {
                    "protein_id": protein_id,
                    "sequence": protein.sequence.sequence,
                    "sequence_md5": protein.sequence.sequence_md5,
                    "fragment_start": protein.sequence.fragment_start,
                    "fragment_end": protein.sequence.fragment_end
                }
            )
            
        self.logger.info(f"Created protein {protein.source_id} with ID {protein_id}")
        return protein_id
    
    def update(self, protein: Protein) -> bool:
        """Update an existing protein
        
        Args:
            protein: Protein to update
            
        Returns:
            True if update was successful
            
        Raises:
            ValidationError: If protein validation fails
        """
        # Validate protein
        protein.validate()
        
        # Ensure protein ID is set
        if protein.id is None:
            raise ValidationError("Protein ID is required for update")
        
        # Check if protein exists
        existing = self.get_by_id(protein.id)
        if not existing:
            raise ValidationError(f"Protein with ID {protein.id} does not exist")
        
        # Update protein
        rows_updated = self.db.update(
            "ecod_schema.protein",
            {
                "pdb_id": protein.pdb_id,
                "chain_id": protein.chain_id,
                "source_id": protein.source_id,
                "length": protein.length,
                "unp_acc": protein.unp_acc,
                "name": protein.name,
                "type": protein.type,
                "tax_id": protein.tax_id
            },
            "id = %s",
            (protein.id,)
        )
        
        # Update sequence if provided
        if protein.sequence:
            # Check if sequence exists
            query = "SELECT id FROM ecod_schema.protein_sequence WHERE protein_id = %s"
            sequence_rows = self.db.execute_query(query, (protein.id,))
            
            if sequence_rows:
                # Update existing sequence
                self.db.update(
                    "ecod_schema.protein_sequence",
                    {
                        "sequence": protein.sequence.sequence,
                        "sequence_md5": protein.sequence.sequence_md5,
                        "fragment_start": protein.sequence.fragment_start,
                        "fragment_end": protein.sequence.fragment_end
                    },
                    "protein_id = %s",
                    (protein.id,)
                )
            else:
                # Insert new sequence
                self.db.insert(
                    "ecod_schema.protein_sequence",
                    {
                        "protein_id": protein.id,
                        "sequence": protein.sequence.sequence,
                        "sequence_md5": protein.sequence.sequence_md5,
                        "fragment_start": protein.sequence.fragment_start,
                        "fragment_end": protein.sequence.fragment_end
                    }
                )
                
        self.logger.info(f"Updated protein with ID {protein.id}")
        return rows_updated > 0
    
    def delete(self, protein_id: int) -> bool:
        """Delete a protein
        
        Args:
            protein_id: ID of protein to delete
            
        Returns:
            True if delete was successful
        """
        # Check if protein exists
        existing = self.get_by_id(protein_id)
        if not existing:
            return False
            
        try:
            # Delete related records first
            self.db.delete("ecod_schema.protein_sequence", "protein_id = %s", (protein_id,))
            
            # Now delete the protein
            rows_deleted = self.db.delete("ecod_schema.protein", "id = %s", (protein_id,))
            
            self.logger.info(f"Deleted protein with ID {protein_id}")
            return rows_deleted > 0
        except (DatabaseError, QueryError) as e:
            self.logger.error(f"Error deleting protein {protein_id}: {str(e)}")
            return False
    
    def get_structure(self, protein_id: int) -> Optional[ProteinStructure]:
        """Get protein structure metadata
        
        Args:
            protein_id: Protein ID
            
        Returns:
            ProteinStructure if found, None otherwise
        """
        query = """
        SELECT *
        FROM ecod_schema.protein_structure
        WHERE protein_id = %s
        """
        rows = self.db.execute_dict_query(query, (protein_id,))
        if not rows:
            return None
        return ProteinStructure.from_db_row(rows[0])
    
    def save_structure(self, structure: ProteinStructure) -> bool:
        """Save protein structure metadata
        
        Args:
            structure: ProteinStructure to save
            
        Returns:
            True if save was successful
        """
        # Check if structure exists
        query = "SELECT id FROM ecod_schema.protein_structure WHERE protein_id = %s"
        structure_rows = self.db.execute_query(query, (structure.protein_id,))
        
        if structure_rows:
            # Update existing structure
            rows_updated = self.db.update(
                "ecod_schema.protein_structure",
                {
                    "resolution": structure.resolution,
                    "experimental_method": structure.experimental_method,
                    "deposition_date": structure.deposition_date,
                    "release_date": structure.release_date,
                    "r_factor": structure.r_factor
                },
                "protein_id = %s",
                (structure.protein_id,)
            )
            return rows_updated > 0
        else:
            # Insert new structure
            self.db.insert(
                "ecod_schema.protein_structure",
                {
                    "protein_id": structure.protein_id,
                    "resolution": structure.resolution,
                    "experimental_method": structure.experimental_method,
                    "deposition_date": structure.deposition_date,
                    "release_date": structure.release_date,
                    "r_factor": structure.r_factor
                }
            )
            return True
    
    def count_proteins(self, criteria: Optional[Dict[str, Any]] = None) -> int:
        """Count proteins, optionally filtered by criteria
        
        Args:
            criteria: Filter criteria (key-value pairs)
            
        Returns:
            Number of matching proteins
        """
        query = "SELECT COUNT(*) FROM ecod_schema.protein p"
        params = []
        
        # Add criteria if provided
        if criteria:
            conditions = []
            for key, value in criteria.items():
                conditions.append(f"p.{key} = %s")
                params.append(value)
                
            if conditions:
                query += " WHERE " + " AND ".join(conditions)
                
        result = self.db.execute_query(query, tuple(params) if params else None)
        return result[0][0] if result else 0
    
    def get_unprocessed_proteins(self, limit: Optional[int] = None) -> List[Protein]:
        """Get proteins that haven't been processed yet
        
        Args:
            limit: Maximum number of proteins to return
            
        Returns:
            List of unprocessed proteins
        """
        query = """
        SELECT p.*, ps.sequence, ps.sequence_md5
        FROM ecod_schema.protein p
        JOIN ecod_schema.protein_sequence ps ON p.id = ps.protein_id
        LEFT JOIN (
            SELECT DISTINCT protein_id 
            FROM ecod_schema.process_status 
            WHERE status IN ('success', 'completed')
        ) ps_done ON p.id = ps_done.protein_id
        WHERE 
            ps_done.protein_id IS NULL
            AND ps.sequence IS NOT NULL
        ORDER BY 
            p.id
        """
        
        if limit is not None:
            query += f" LIMIT {limit}"
            
        rows = self.db.execute_dict_query(query)
        return [Protein.from_db_row(row, include_sequence=True) for row in rows]