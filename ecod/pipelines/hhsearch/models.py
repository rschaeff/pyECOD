#!/usr/bin/env python3
"""
Models for HHSearch registration service
"""
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional, List, Dict, Any
from enum import Enum


class RegistrationStatus(Enum):
    """Status of HHSearch result registration"""
    PENDING = "pending"
    FOUND = "found"
    CONVERTING = "converting"
    REGISTERED = "registered"
    FAILED = "failed"
    SKIPPED = "skipped"


@dataclass
class HHSearchFile:
    """Represents an HHSearch result file"""
    process_id: int
    pdb_id: str
    chain_id: str
    file_type: str  # 'hhr' or 'hh_xml'
    file_path: Path
    exists: bool = False
    size: int = 0
    last_modified: Optional[datetime] = None
    
    @property
    def protein_id(self) -> str:
        return f"{self.pdb_id}_{self.chain_id}"
    
    def validate(self) -> bool:
        """Validate file exists and has content"""
        if not self.exists:
            return False
        return self.size > 0


@dataclass
class RegistrationResult:
    """Result of registering HHSearch results for a chain"""
    process_id: int
    pdb_id: str
    chain_id: str
    status: RegistrationStatus
    hhr_registered: bool = False
    xml_generated: bool = False
    xml_registered: bool = False
    error: Optional[str] = None
    processing_time: float = 0.0
    
    @property
    def success(self) -> bool:
        return self.status == RegistrationStatus.REGISTERED
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'process_id': self.process_id,
            'protein_id': f"{self.pdb_id}_{self.chain_id}",
            'status': self.status.value,
            'hhr_registered': self.hhr_registered,
            'xml_generated': self.xml_generated,
            'xml_registered': self.xml_registered,
            'error': self.error,
            'processing_time': self.processing_time
        }


@dataclass
class BatchRegistrationResult:
    """Results from registering a batch of HHSearch results"""
    batch_id: int
    total_chains: int = 0
    registered: int = 0
    skipped: int = 0
    failed: int = 0
    results: List[RegistrationResult] = field(default_factory=list)
    start_time: datetime = field(default_factory=datetime.now)
    end_time: Optional[datetime] = None
    
    @property
    def success_rate(self) -> float:
        if self.total_chains == 0:
            return 0.0
        return (self.registered / self.total_chains) * 100
    
    @property
    def processing_time(self) -> float:
        if not self.end_time:
            return 0.0
        return (self.end_time - self.start_time).total_seconds()
    
    def add_result(self, result: RegistrationResult):
        """Add a registration result"""
        self.results.append(result)
        if result.success:
            self.registered += 1
        elif result.status == RegistrationStatus.SKIPPED:
            self.skipped += 1
        else:
            self.failed += 1
    
    def finalize(self):
        """Mark batch processing as complete"""
        self.end_time = datetime.now()
        self.total_chains = len(self.results)


@dataclass
class ServiceConfig:
    """Configuration for HHSearch registration service"""
    force_regenerate: bool = False
    min_probability: float = 20.0  # Minimum HHSearch probability to include
    validate_hhsearch: bool = True  # Validate HHR files are from HHSearch not HHblits
    migrate_legacy_files: bool = True
    parallel_processing: bool = True
    max_workers: int = 4
    process_timeout: int = 300  # seconds
