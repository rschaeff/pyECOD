#!/usr/bin/env python3
"""
Models for pipeline orchestration
"""
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import Optional, List, Dict, Any


class PipelineMode(Enum):
    """Pipeline execution modes"""
    FULL = "full"                    # Complete pipeline: BLAST + HHSearch + Domain Analysis
    BLAST_ONLY = "blast_only"        # BLAST + Domain Analysis (skip HHSearch)
    ANALYSIS_ONLY = "analysis_only"  # Domain Analysis only (use existing evidence)


class PipelineStage(Enum):
    """Pipeline stages in execution order"""
    INITIALIZATION = "initialization"
    BLAST_CHAIN = "blast_chain"
    BLAST_DOMAIN = "blast_domain"
    HHSEARCH_PROFILE = "hhsearch_profile"
    HHSEARCH_SEARCH = "hhsearch_search"
    HHSEARCH_REGISTRATION = "hhsearch_registration"
    DOMAIN_SUMMARY = "domain_summary"
    DOMAIN_PARTITION = "domain_partition"
    CLASSIFICATION = "classification"
    COMPLETE = "complete"


@dataclass
class StageResult:
    """Result from executing a pipeline stage"""
    stage: PipelineStage
    success: bool
    start_time: datetime
    end_time: Optional[datetime] = None
    items_processed: int = 0
    items_failed: int = 0
    error: Optional[str] = None
    details: Dict[str, Any] = field(default_factory=dict)
    
    @property
    def duration(self) -> float:
        if not self.end_time:
            return 0.0
        return (self.end_time - self.start_time).total_seconds()
    
    def finalize(self, success: bool = True, error: Optional[str] = None):
        """Mark stage as complete"""
        self.end_time = datetime.now()
        self.success = success
        if error:
            self.error = error


@dataclass
class PipelineRun:
    """Represents a complete pipeline execution"""
    run_id: str
    batch_id: int
    mode: PipelineMode
    total_proteins: int
    start_time: datetime = field(default_factory=datetime.now)
    end_time: Optional[datetime] = None
    stages: List[StageResult] = field(default_factory=list)
    
    @property
    def is_complete(self) -> bool:
        return self.end_time is not None
    
    @property
    def duration(self) -> float:
        if not self.end_time:
            return (datetime.now() - self.start_time).total_seconds()
        return (self.end_time - self.start_time).total_seconds()
    
    @property
    def success(self) -> bool:
        return all(stage.success for stage in self.stages)
    
    @property
    def current_stage(self) -> Optional[PipelineStage]:
        """Get the currently executing stage"""
        for stage in self.stages:
            if not stage.end_time:
                return stage.stage
        return None
    
    def add_stage_result(self, result: StageResult):
        """Add a stage result"""
        self.stages.append(result)
    
    def finalize(self):
        """Mark pipeline run as complete"""
        self.end_time = datetime.now()
    
    def get_summary(self) -> Dict[str, Any]:
        """Get run summary"""
        return {
            'run_id': self.run_id,
            'batch_id': self.batch_id,
            'mode': self.mode.value,
            'total_proteins': self.total_proteins,
            'duration': self.duration,
            'success': self.success,
            'stages_completed': len([s for s in self.stages if s.success]),
            'stages_failed': len([s for s in self.stages if not s.success]),
            'current_stage': self.current_stage.value if self.current_stage else None,
            'is_complete': self.is_complete
        }


@dataclass 
class OrchestrationConfig:
    """Configuration for orchestration service"""
    # Execution settings
    checkpoint_enabled: bool = True  # Save progress between stages
    continue_on_error: bool = True   # Continue pipeline if stage fails
    parallel_stages: bool = False    # Run independent stages in parallel
    
    # Stage-specific settings
    blast_batch_size: int = 100
    hhsearch_threads: int = 8
    hhsearch_memory: str = "16G"
    domain_analysis_batch_size: int = 100
    
    # Timeouts (seconds)
    stage_timeout: int = 3600  # 1 hour default
    total_timeout: int = 86400  # 24 hours default
