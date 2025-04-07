# ecod/core/models.py
from dataclasses import dataclass, field
from datetime import datetime
from typing import List, Optional, Dict, Any

@dataclass
class Protein:
    pdb_id: str
    chain_id: str
    source_id: str
    length: int
    id: Optional[int] = None
    sequence: Optional[str] = None
    
    @classmethod
    def from_db_row(cls, row: Dict[str, Any]) -> 'Protein':
        return cls(
            id=row.get('id'),
            pdb_id=row.get('pdb_id'),
            chain_id=row.get('chain_id'),
            source_id=row.get('source_id'),
            length=row.get('length')
        )

@dataclass
class Batch:
    batch_name: str
    base_path: str
    type: str  # 'blast', 'hhsearch', 'classification'
    ref_version: str
    total_items: int
    id: Optional[int] = None
    status: str = "created"
    completed_items: int = 0
    created_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    
    @classmethod
    def from_db_row(cls, row: Dict[str, Any]) -> 'Batch':
        return cls(
            id=row.get('id'),
            batch_name=row.get('batch_name'),
            base_path=row.get('base_path'),
            type=row.get('type'),
            ref_version=row.get('ref_version'),
            total_items=row.get('total_items'),
            status=row.get('status'),
            completed_items=row.get('completed_items'),
            created_at=row.get('created_at'),
            completed_at=row.get('completed_at')
        )

@dataclass
class ProcessStatus:
    protein_id: int
    batch_id: int
    current_stage: str
    id: Optional[int] = None
    status: str = "pending"
    is_representative: bool = False
    relative_path: Optional[str] = None
    error_message: Optional[str] = None
    created_at: Optional[datetime] = None
    updated_at: Optional[datetime] = None
    
    # Additional fields to link related data
    protein: Optional[Protein] = None
    files: List['ProcessFile'] = field(default_factory=list)
    
    @classmethod
    def from_db_row(cls, row: Dict[str, Any]) -> 'ProcessStatus':
        return cls(
            id=row.get('id'),
            protein_id=row.get('protein_id'),
            batch_id=row.get('batch_id'),
            current_stage=row.get('current_stage'),
            status=row.get('status'),
            is_representative=row.get('is_representative'),
            relative_path=row.get('relative_path'),
            error_message=row.get('error_message'),
            created_at=row.get('created_at'),
            updated_at=row.get('updated_at')
        )

@dataclass
class ProcessFile:
    process_id: int
    file_type: str
    file_path: str
    id: Optional[int] = None
    file_exists: bool = False
    file_size: Optional[int] = None
    last_checked: Optional[datetime] = None
    
    @classmethod
    def from_db_row(cls, row: Dict[str, Any]) -> 'ProcessFile':
        return cls(
            id=row.get('id'),
            process_id=row.get('process_id'),
            file_type=row.get('file_type'),
            file_path=row.get('file_path'),
            file_exists=row.get('file_exists'),
            file_size=row.get('file_size'),
            last_checked=row.get('last_checked')
        )