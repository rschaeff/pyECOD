#!/usr/bin/env python3
"""
Job and job item models for the ECOD pipeline.
"""
from dataclasses import dataclass, field
from datetime import datetime
from typing import Optional, List, Dict, Any

@dataclass
class JobItem:
    """Job item model"""
    job_id: int
    process_id: int
    id: Optional[int] = None
    status: str = "pending"
    
    @classmethod
    def from_db_row(cls, row: Dict[str, Any]) -> 'JobItem':
        """Create instance from database row"""
        return cls(
            id=row.get('id'),
            job_id=row.get('job_id'),
            process_id=row.get('process_id'),
            status=row.get('status', 'pending')
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'id': self.id,
            'job_id': self.job_id,
            'process_id': self.process_id,
            'status': self.status
        }

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
        """Create instance from database row"""
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
        """Check if job is completed"""
        return self.status in ('completed', 'success')
    
    def is_failed(self) -> bool:
        """Check if job failed"""
        return self.status in ('failed', 'error')
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
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