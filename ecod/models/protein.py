#!/usr/bin/env python3
"""
Protein models for the ECOD pipeline
Defines data models for proteins and their sequences
"""
from dataclasses import dataclass, field
from datetime import datetime
from typing import Optional, List, Dict, Any

from ecod.exceptions import ValidationError

@dataclass
class ProteinSequence:
    """Protein sequence model"""
    sequence: str
    sequence_md5: str
    protein_id: Optional[int] = None
    id: Optional[int] = None
    fragment_start: Optional[int] = None
    fragment_end: Optional[int] = None
    created_at: Optional[datetime] = None
    
    @property
    def length(self) -> int:
        """Get sequence length"""
        return len(self.sequence)
    
    @classmethod
    def from_db_row(cls, row: Dict[str, Any]) -> 'ProteinSequence':
        """Create instance from database row
        
        Args:
            row: Database row as dictionary
            
        Returns:
            ProteinSequence instance
        """
        return cls(
            id=row.get('id'),
            protein_id=row.get('protein_id'),
            sequence=row.get('sequence', ''),
            sequence_md5=row.get('sequence_md5', ''),
            fragment_start=row.get('fragment_start'),
            fragment_end=row.get('fragment_end'),
            created_at=row.get('created_at')
        )
    
    def validate(self) -> bool:
        """Validate sequence
        
        Returns:
            True if valid
            
        Raises:
            ValidationError: If validation fails
        """
        if not self.sequence:
            raise ValidationError("Protein sequence cannot be empty")
        
        # Check for valid amino acid characters
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        invalid_chars = set(self.sequence.upper()) - valid_aa
        if invalid_chars:
            raise ValidationError(f"Invalid amino acids in sequence: {', '.join(invalid_chars)}")
        
        return True
    
    def to_fasta(self, header: Optional[str] = None) -> str:
        """Convert to FASTA format
        
        Args:
            header: FASTA header (default: protein_id)
            
        Returns:
            FASTA formatted string
        """
        header_text = header or f"protein_{self.protein_id}"
        return f">{header_text}\n{self.sequence}\n"
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary
        
        Returns:
            Dictionary representation
        """
        return {
            'id': self.id,
            'protein_id': self.protein_id,
            'sequence': self.sequence,
            'sequence_md5': self.sequence_md5,
            'fragment_start': self.fragment_start,
            'fragment_end': self.fragment_end,
            'length': self.length,
            'created_at': self.created_at.isoformat() if self.created_at else None
        }

@dataclass
class ProteinStructure:
    """Protein structure metadata model"""
    protein_id: int
    id: Optional[int] = None
    resolution: Optional[float] = None
    experimental_method: Optional[str] = None
    deposition_date: Optional[datetime] = None
    release_date: Optional[datetime] = None
    r_factor: Optional[float] = None
    created_at: Optional[datetime] = None
    
    @classmethod
    def from_db_row(cls, row: Dict[str, Any]) -> 'ProteinStructure':
        """Create instance from database row
        
        Args:
            row: Database row as dictionary
            
        Returns:
            ProteinStructure instance
        """
        return cls(
            id=row.get('id'),
            protein_id=row.get('protein_id'),
            resolution=row.get('resolution'),
            experimental_method=row.get('experimental_method'),
            deposition_date=row.get('deposition_date'),
            release_date=row.get('release_date'),
            r_factor=row.get('r_factor'),
            created_at=row.get('created_at')
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary
        
        Returns:
            Dictionary representation
        """
        return {
            'id': self.id,
            'protein_id': self.protein_id,
            'resolution': self.resolution,
            'experimental_method': self.experimental_method,
            'deposition_date': self.deposition_date.isoformat() if self.deposition_date else None,
            'release_date': self.release_date.isoformat() if self.release_date else None,
            'r_factor': self.r_factor,
            'created_at': self.created_at.isoformat() if self.created_at else None
        }

@dataclass
class Protein:
    """Protein model"""
    pdb_id: str
    chain_id: str
    source_id: str
    length: int
    id: Optional[int] = None
    unp_acc: Optional[str] = None
    name: Optional[str] = None
    type: Optional[str] = None
    tax_id: Optional[int] = None
    created_at: Optional[datetime] = None
    updated_at: Optional[datetime] = None
    
    # Related data
    sequence: Optional[ProteinSequence] = None
    structure: Optional[ProteinStructure] = None
    
    @classmethod
    def from_db_row(cls, row: Dict[str, Any], include_sequence: bool = False) -> 'Protein':
        """Create instance from database row
        
        Args:
            row: Database row as dictionary
            include_sequence: Whether to include sequence data
            
        Returns:
            Protein instance
        """
        protein = cls(
            id=row.get('id'),
            pdb_id=row.get('pdb_id', ''),
            chain_id=row.get('chain_id', ''),
            source_id=row.get('source_id', ''),
            length=row.get('length', 0),
            unp_acc=row.get('unp_acc'),
            name=row.get('name'),
            type=row.get('type'),
            tax_id=row.get('tax_id'),
            created_at=row.get('created_at'),
            updated_at=row.get('updated_at')
        )
        
        # Add sequence if included in row
        if include_sequence and 'sequence' in row and row['sequence']:
            protein.sequence = ProteinSequence(
                protein_id=row.get('id'),
                sequence=row['sequence'],
                sequence_md5=row.get('sequence_md5', ''),
                fragment_start=row.get('fragment_start'),
                fragment_end=row.get('fragment_end')
            )
            
        return protein
    
    def validate(self) -> bool:
        """Validate protein data
        
        Returns:
            True if valid
            
        Raises:
            ValidationError: If validation fails
        """
        if not self.pdb_id:
            raise ValidationError("PDB ID cannot be empty")
        
        if not self.chain_id:
            raise ValidationError("Chain ID cannot be empty")
        
        if not self.source_id:
            raise ValidationError("Source ID cannot be empty")
        
        # Validate PDB ID format (typically 4 characters)
        import re
        if not re.match(r'^[A-Za-z0-9]{4}$', self.pdb_id):
            raise ValidationError(f"Invalid PDB ID format: {self.pdb_id}")
        
        # Validate sequence if present
        if self.sequence:
            self.sequence.validate()
            
            # Validate length consistency
            if self.length != self.sequence.length:
                raise ValidationError(f"Protein length ({self.length}) does not match sequence length ({self.sequence.length})")
        
        return True
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary
        
        Returns:
            Dictionary representation
        """
        result = {
            'id': self.id,
            'pdb_id': self.pdb_id,
            'chain_id': self.chain_id,
            'source_id': self.source_id,
            'length': self.length,
            'unp_acc': self.unp_acc,
            'name': self.name,
            'type': self.type,
            'tax_id': self.tax_id,
            'created_at': self.created_at.isoformat() if self.created_at else None,
            'updated_at': self.updated_at.isoformat() if self.updated_at else None
        }
        
        # Include related data if available
        if self.sequence:
            result['sequence'] = self.sequence.to_dict()
            
        if self.structure:
            result['structure'] = self.structure.to_dict()
            
        return result

@dataclass
class ChainSequence:
    """Chain sequence model"""
    chain_id: int
    sequence: str
    sequence_length: int
    sequence_md5: str
    id: Optional[int] = None
    created_at: Optional[datetime] = None
    
    @classmethod
    def from_db_row(cls, row: Dict[str, Any]) -> 'ChainSequence':
        """Create instance from database row
        
        Args:
            row: Database row as dictionary
            
        Returns:
            ChainSequence instance
        """
        return cls(
            id=row.get('id'),
            chain_id=row.get('chain_id'),
            sequence=row.get('sequence', ''),
            sequence_length=row.get('sequence_length', 0),
            sequence_md5=row.get('sequence_md5', ''),
            created_at=row.get('created_at')
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary
        
        Returns:
            Dictionary representation
        """
        return {
            'id': self.id,
            'chain_id': self.chain_id,
            'sequence': self.sequence,
            'sequence_length': self.sequence_length,
            'sequence_md5': self.sequence_md5,
            'created_at': self.created_at.isoformat() if self.created_at else None
        }

@dataclass
class PDBChain:
    """PDB chain model"""
    source_id: str
    pdb_id: str
    chain_id: str
    id: Optional[int] = None
    pdb_entry_id: Optional[int] = None
    asym_id: Optional[str] = None
    entity_id: Optional[int] = None
    polymer_type: Optional[str] = None
    unp_acc: Optional[str] = None
    length: Optional[int] = None
    tax_id: Optional[int] = None
    created_at: Optional[datetime] = None
    updated_at: Optional[datetime] = None
    
    # Related data
    sequence: Optional[ChainSequence] = None
    
    @classmethod
    def from_db_row(cls, row: Dict[str, Any], include_sequence: bool = False) -> 'PDBChain':
        """Create instance from database row
        
        Args:
            row: Database row as dictionary
            include_sequence: Whether to include sequence data
            
        Returns:
            PDBChain instance
        """
        chain = cls(
            id=row.get('id'),
            pdb_entry_id=row.get('pdb_entry_id'),
            source_id=row.get('source_id', ''),
            pdb_id=row.get('pdb_id', ''),
            chain_id=row.get('chain_id', ''),
            asym_id=row.get('asym_id'),
            entity_id=row.get('entity_id'),
            polymer_type=row.get('polymer_type'),
            unp_acc=row.get('unp_acc'),
            length=row.get('length'),
            tax_id=row.get('tax_id'),
            created_at=row.get('created_at'),
            updated_at=row.get('updated_at')
        )
        
        # Add sequence if included in row
        if include_sequence and 'sequence' in row and row['sequence']:
            chain.sequence = ChainSequence(
                chain_id=row.get('id'),
                sequence=row['sequence'],
                sequence_length=row.get('sequence_length', len(row['sequence'])),
                sequence_md5=row.get('sequence_md5', '')
            )
            
        return chain
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary
        
        Returns:
            Dictionary representation
        """
        result = {
            'id': self.id,
            'pdb_entry_id': self.pdb_entry_id,
            'source_id': self.source_id,
            'pdb_id': self.pdb_id,
            'chain_id': self.chain_id,
            'asym_id': self.asym_id,
            'entity_id': self.entity_id,
            'polymer_type': self.polymer_type,
            'unp_acc': self.unp_acc,
            'length': self.length,
            'tax_id': self.tax_id,
            'created_at': self.created_at.isoformat() if self.created_at else None,
            'updated_at': self.updated_at.isoformat() if self.updated_at else None
        }
        
        # Include related data if available
        if self.sequence:
            result['sequence'] = self.sequence.to_dict()
            
        return result

@dataclass
class PDBEntry:
    """PDB entry model"""
    pdb_id: str
    id: Optional[int] = None
    title: Optional[str] = None
    experimental_method: Optional[str] = None
    resolution: Optional[float] = None
    r_factor: Optional[float] = None
    deposition_date: Optional[datetime] = None
    release_date: Optional[datetime] = None
    revision_date: Optional[datetime] = None
    version: int = 1
    status: str = 'CURRENT'
    last_updated: Optional[datetime] = None
    
    # Related data
    chains: List[PDBChain] = field(default_factory=list)
    
    @classmethod
    def from_db_row(cls, row: Dict[str, Any]) -> 'PDBEntry':
        """Create instance from database row
        
        Args:
            row: Database row as dictionary
            
        Returns:
            PDBEntry instance
        """
        return cls(
            id=row.get('id'),
            pdb_id=row.get('pdb_id', ''),
            title=row.get('title'),
            experimental_method=row.get('experimental_method'),
            resolution=row.get('resolution'),
            r_factor=row.get('r_factor'),
            deposition_date=row.get('deposition_date'),
            release_date=row.get('release_date'),
            revision_date=row.get('revision_date'),
            version=row.get('version', 1),
            status=row.get('status', 'CURRENT'),
            last_updated=row.get('last_updated')
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary
        
        Returns:
            Dictionary representation
        """
        result = {
            'id': self.id,
            'pdb_id': self.pdb_id,
            'title': self.title,
            'experimental_method': self.experimental_method,
            'resolution': self.resolution,
            'r_factor': self.r_factor,
            'deposition_date': self.deposition_date.isoformat() if self.deposition_date else None,
            'release_date': self.release_date.isoformat() if self.release_date else None,
            'revision_date': self.revision_date.isoformat() if self.revision_date else None,
            'version': self.version,
            'status': self.status,
            'last_updated': self.last_updated.isoformat() if self.last_updated else None
        }
        
        # Include related data if available
        if self.chains:
            result['chains'] = [chain.to_dict() for chain in self.chains]
            
        return result