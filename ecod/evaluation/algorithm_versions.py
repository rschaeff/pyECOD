# ecod/evaluation/algorithm_versions/manager.py
from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, List, Any, Optional
from enum import Enum

class AlgorithmStatus(Enum):
    DEVELOPMENT = "development"
    TESTING = "testing" 
    PRODUCTION = "production"
    DEPRECATED = "deprecated"

@dataclass
class AlgorithmVersion:
    """Comprehensive algorithm version tracking"""
    version_id: str
    name: str
    description: str
    parent_version: Optional[str] = None
    status: AlgorithmStatus = AlgorithmStatus.DEVELOPMENT
    
    # Algorithm configuration
    partition_config: Dict[str, Any] = field(default_factory=dict)
    evidence_weights: Dict[str, float] = field(default_factory=dict)
    coverage_thresholds: Dict[str, float] = field(default_factory=dict)
    
    # Tracking
    created_at: datetime = field(default_factory=datetime.now)
    created_by: str = ""
    test_results: List[Dict[str, Any]] = field(default_factory=list)
    
    def to_config_dict(self) -> Dict[str, Any]:
        """Convert to configuration dictionary for pipeline use"""
        return {
            'version_id': self.version_id,
            'domain_analysis': {
                'partition': self.partition_config,
                'evidence_weights': self.evidence_weights,
                'coverage_thresholds': self.coverage_thresholds
            }
        }

class AlgorithmVersionManager:
    """Manages algorithm versions with database integration"""
    
    def __init__(self, context: ApplicationContext):
        self.context = context
        self.db = context.db
        self._ensure_tables()
    
    def register_version(self, algorithm: AlgorithmVersion) -> int:
        """Register new algorithm version"""
        
    def get_version(self, version_id: str) -> Optional[AlgorithmVersion]:
        """Get algorithm version by ID"""
        
    def list_versions(self, status: Optional[AlgorithmStatus] = None) -> List[AlgorithmVersion]:
        """List algorithm versions"""
        
    def promote_version(self, version_id: str, new_status: AlgorithmStatus) -> bool:
        """Promote algorithm version to new status"""
