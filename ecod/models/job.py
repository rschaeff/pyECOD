#!/usr/bin/env python3
"""
Job and workflow models for the ECOD pipeline
Defines data models for processing jobs, batches, and workflows.
"""
from dataclasses import dataclass, field
from datetime import datetime
from typing import Optional, List, Dict, Any, Set

from ecod.exceptions import ValidationError

@dataclass
class Batch:
    """Batch processing model"""
    batch_name: str
    base_path: str
    type: str  # 'blast', 'hhsearch', 'classification', etc.
    ref_version: str
    total_items: int
    id: Optional[int] = None
    status: str = "created"
    completed_items: int = 0
    created_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    
    @classmethod
    def from_db_row(cls, row: Dict[str, Any]) -> 'Batch':
        """Create instance from database row
        
        Args:
            row: Database row as dictionary
            
        Returns:
            Batch instance
        """
        return cls(
            id=row.get('id'),
            batch_name=row.get('batch_name', ''),
            base_path=row.get('base_path', ''),
            type=row.get('type', ''),
            ref_version=row.get('ref_version', ''),
            total_items=row.get('total_items', 0),
            status=row.get('status', 'created'),
            completed_items=row.get('completed_items', 0),
            created_at=row.get('created_at'),
            completed_at=row.get('completed_at')
        )
    
    def is_complete(self) -> bool:
        """Check if batch is complete
        
        Returns:
            True if complete
        """
        return self.status == 'completed' or self.completed_items >= self.total_items
    
    def progress_percentage(self) -> float:
        """Calculate progress percentage
        
        Returns:
            Progress percentage (0-100)
        """
        if self.total_items == 0:
            return 0.0
        return (self.completed_items / self.total_items) * 100.0
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary
        
        Returns:
            Dictionary representation
        """
        return {
            'id': self.id,
            'batch_name': self.batch_name,
            'base_path': self.base_path,
            'type': self.type,
            'ref_version': self.ref_version,
            'total_items': self.total_items,
            'status': self.status,
            'completed_items': self.completed_items,
            'progress_percentage': self.progress_percentage(),
            'created_at': self.created_at.isoformat() if self.created_at else None,
            'completed_at': self.completed_at.isoformat() if self.completed_at else None
        }

@dataclass
class ProcessFile:
    """Process file model"""
    process_id: int
    file_type: str
    file_path: str
    id: Optional[int] = None
    file_exists: bool = False
    file_size: Optional[int] = None
    last_checked: Optional[datetime] = None
    
    @classmethod
    def from_db_row(cls, row: Dict[str, Any]) -> 'ProcessFile':
        """Create instance from database row
        
        Args:
            row: Database row as dictionary
            
        Returns:
            ProcessFile instance
        """
        return cls(
            id=row.get('id'),
            process_id=row.get('process_id'),
            file_type=row.get('file_type', ''),
            file_path=row.get('file_path', ''),
            file_exists=row.get('file_exists', False),
            file_size=row.get('file_size'),
            last_checked=row.get('last_checked')
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary
        
        Returns:
            Dictionary representation
        """
        return {
            'id': self.id,
            'process_id': self.process_id,
            'file_type': self.file_type,
            'file_path': self.file_path,
            'file_exists': self.file_exists,
            'file_size': self.file_size,
            'last_checked': self.last_checked.isoformat() if self.last_checked else None
        }

@dataclass
class ProcessStatus:
    """Process status model"""
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
    
    # Related data
    files: List[ProcessFile] = field(default_factory=list)
    
    @classmethod
    def from_db_row(cls, row: Dict[str, Any]) -> 'ProcessStatus':
        """Create instance from database row
        
        Args:
            row: Database row as dictionary
            
        Returns:
            ProcessStatus instance
        """
        return cls(
            id=row.get('id'),
            protein_id=row.get('protein_id'),
            batch_id=row.get('batch_id'),
            current_stage=row.get('current_stage', ''),
            status=row.get('status', 'pending'),
            is_representative=row.get('is_representative', False),
            relative_path=row.get('relative_path'),
            error_message=row.get('error_message'),
            created_at=row.get('created_at'),
            updated_at=row.get('updated_at')
        )
    
    def is_completed(self) -> bool:
        """Check if process is completed
        
        Returns:
            True if completed
        """
        return self.status in ('success', 'completed')
    
    def is_failed(self) -> bool:
        """Check if process failed
        
        Returns:
            True if failed
        """
        return self.status in ('error', 'failed')
    
    def has_file_type(self, file_type: str) -> bool:
        """Check if process has a file of specified type
        
        Args:
            file_type: File type to check
            
        Returns:
            True if process has file of specified type
        """
        return any(f.file_type == file_type for f in self.files)
    
    def get_file(self, file_type: str) -> Optional[ProcessFile]:
        """Get file of specified type
        
        Args:
            file_type: File type to get
            
        Returns:
            ProcessFile if found, None otherwise
        """
        for file in self.files:
            if file.file_type == file_type:
                return file
        return None
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary
        
        Returns:
            Dictionary representation
        """
        result = {
            'id': self.id,
            'protein_id': self.protein_id,
            'batch_id': self.batch_id,
            'current_stage': self.current_stage,
            'status': self.status,
            'is_representative': self.is_representative,
            'relative_path': self.relative_path,
            'error_message': self.error_message,
            'created_at': self.created_at.isoformat() if self.created_at else None,
            'updated_at': self.updated_at.isoformat() if self.updated_at else None
        }
        
        # Include related data if available
        if self.files:
            result['files'] = [f.to_dict() for f in self.files]
            
        return result

@dataclass
class Job:
    """Job model"""
    batch_id: int
    job_type: str
    id: Optional[int] = None
    slurm_job_id: Optional[str] = None
    status: str = "submitted"
    items_count: Optional[int] = None
    submission_time: Optional[datetime] = None
    completion_time: Optional[datetime] = None
    
    # Related data
    items: List[int] = field(default_factory=list)  # List of process_ids
    
    @classmethod
    def from_db_row(cls, row: Dict[str, Any]) -> 'Job':
        """Create instance from database row
        
        Args:
            row: Database row as dictionary
            
        Returns:
            Job instance
        """
        return cls(
            id=row.get('id'),
            batch_id=row.get('batch_id'),
            job_type=row.get('job_type', ''),
            slurm_job_id=row.get('slurm_job_id'),
            status=row.get('status', 'submitted'),
            items_count=row.get('items_count'),
            submission_time=row.get('submission_time'),
            completion_time=row.get('completion_time')
        )
    
    def is_completed(self) -> bool:
        """Check if job is completed
        
        Returns:
            True if completed
        """
        return self.status in ('completed', 'success')
    
    def is_failed(self) -> bool:
        """Check if job failed
        
        Returns:
            True if failed
        """
        return self.status in ('failed', 'error')
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary
        
        Returns:
            Dictionary representation
        """
        result = {
            'id': self.id,
            'batch_id': self.batch_id,
            'job_type': self.job_type,
            'slurm_job_id': self.slurm_job_id,
            'status': self.status,
            'items_count': self.items_count,
            'submission_time': self.submission_time.isoformat() if self.submission_time else None,
            'completion_time': self.completion_time.isoformat() if self.completion_time else None
        }
        
        # Include items if available
        if self.items:
            result['items'] = self.items
            
        return result

@dataclass
class JobItem:
    """Job item model"""
    job_id: int
    process_id: int
    id: Optional[int] = None
    status: str = "pending"
    
    @classmethod
    def from_db_row(cls, row: Dict[str, Any]) -> 'JobItem':
        """Create instance from database row
        
        Args:
            row: Database row as dictionary
            
        Returns:
            JobItem instance
        """
        return cls(
            id=row.get('id'),
            job_id=row.get('job_id'),
            process_id=row.get('process_id'),
            status=row.get('status', 'pending')
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary
        
        Returns:
            Dictionary representation
        """
        return {
            'id': self.id,
            'job_id': self.job_id,
            'process_id': self.process_id,
            'status': self.status
        }

@dataclass
class ECODVersion:
    """ECOD version model"""
    version_name: str
    id: Optional[int] = None
    release_date: Optional[datetime] = None
    is_current: bool = False
    
    @classmethod
    def from_db_row(cls, row: Dict[str, Any]) -> 'ECODVersion':
        """Create instance from database row
        
        Args:
            row: Database row as dictionary
            
        Returns:
            ECODVersion instance
        """
        return cls(
            id=row.get('id'),
            version_name=row.get('version_name', ''),
            release_date=row.get('release_date'),
            is_current=row.get('is_current', False)
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary
        
        Returns:
            Dictionary representation
        """
        return {
            'id': self.id,
            'version_name': self.version_name,
            'release_date': self.release_date.isoformat() if self.release_date else None,
            'is_current': self.is_current
        }

@dataclass
class ReferenceResource:
    """Reference resource model"""
    version_id: int
    resource_type: str
    resource_path: str
    id: Optional[int] = None
    
    @classmethod
    def from_db_row(cls, row: Dict[str, Any]) -> 'ReferenceResource':
        """Create instance from database row
        
        Args:
            row: Database row as dictionary
            
        Returns:
            ReferenceResource instance
        """
        return cls(
            id=row.get('id'),
            version_id=row.get('version_id'),
            resource_type=row.get('resource_type', ''),
            resource_path=row.get('resource_path', '')
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary
        
        Returns:
            Dictionary representation
        """
        return {
            'id': self.id,
            'version_id': self.version_id,
            'resource_type': self.resource_type,
            'resource_path': self.resource_path
        }