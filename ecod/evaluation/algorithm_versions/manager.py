#!/usr/bin/env python3
"""
Complete Algorithm Version Manager Implementation

This module provides full functionality for managing algorithm versions
with database integration.
"""

import json
import logging
from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, List, Any, Optional
from enum import Enum
from pathlib import Path

from ecod.core.context import ApplicationContext
from ecod.exceptions import ValidationError, ConfigurationError


class AlgorithmStatus(Enum):
    """Algorithm deployment status"""
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
    behavioral_flags: Dict[str, bool] = field(default_factory=dict)
    
    # Tracking
    created_at: datetime = field(default_factory=datetime.now)
    created_by: str = ""
    notes: str = ""
    
    # Database tracking (not part of config)
    database_id: Optional[int] = field(default=None, init=False)
    
    def __post_init__(self):
        """Validate algorithm version after creation"""
        if not self.version_id:
            raise ValidationError("Algorithm version_id is required")
        if not self.name:
            raise ValidationError("Algorithm name is required")
    
    def to_config_dict(self) -> Dict[str, Any]:
        """Convert to configuration dictionary for pipeline use"""
        return {
            'version_id': self.version_id,
            'domain_analysis': {
                'partition': self.partition_config,
                'evidence_weights': self.evidence_weights,
                'coverage_thresholds': self.coverage_thresholds,
                'behavioral_flags': self.behavioral_flags
            }
        }
    
    def to_database_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for database storage"""
        config_data = {
            'partition_config': self.partition_config,
            'evidence_weights': self.evidence_weights,
            'coverage_thresholds': self.coverage_thresholds,
            'behavioral_flags': self.behavioral_flags
        }
        
        return {
            'version_id': self.version_id,
            'name': self.name,
            'description': self.description,
            'status': self.status.value,
            'config_data': json.dumps(config_data),
            'created_by': self.created_by,
            'notes': self.notes
        }
    
    @classmethod
    def from_database_row(cls, row: Dict[str, Any]) -> 'AlgorithmVersion':
        """Create AlgorithmVersion from database row"""
        config_data = json.loads(row['config_data']) if isinstance(row['config_data'], str) else row['config_data']
        
        instance = cls(
            version_id=row['version_id'],
            name=row['name'],
            description=row['description'],
            parent_version=row.get('parent_version_id'),
            status=AlgorithmStatus(row['status']),
            partition_config=config_data.get('partition_config', {}),
            evidence_weights=config_data.get('evidence_weights', {}),
            coverage_thresholds=config_data.get('coverage_thresholds', {}),
            behavioral_flags=config_data.get('behavioral_flags', {}),
            created_at=row.get('created_at', datetime.now()),
            created_by=row.get('created_by', ''),
            notes=row.get('notes', '')
        )
        instance.database_id = row.get('id')
        return instance
    
    @classmethod
    def from_config_file(cls, config_path: str) -> 'AlgorithmVersion':
        """Load algorithm version from YAML/JSON configuration file"""
        import yaml
        
        config_path = Path(config_path)
        if not config_path.exists():
            raise ConfigurationError(f"Algorithm config file not found: {config_path}")
        
        with open(config_path, 'r') as f:
            if config_path.suffix.lower() == '.json':
                config_data = json.load(f)
            else:
                config_data = yaml.safe_load(f)
        
        # Extract algorithm metadata
        metadata = config_data.get('algorithm', {})
        domain_analysis = config_data.get('domain_analysis', {})
        
        return cls(
            version_id=metadata.get('version_id', ''),
            name=metadata.get('name', ''),
            description=metadata.get('description', ''),
            parent_version=metadata.get('parent_version_id'),
            status=AlgorithmStatus(metadata.get('status', 'development')),
            partition_config=domain_analysis.get('partition', {}),
            evidence_weights=domain_analysis.get('evidence_weights', {}),
            coverage_thresholds=domain_analysis.get('coverage_thresholds', {}),
            behavioral_flags=domain_analysis.get('behavioral_flags', {}),
            created_by=metadata.get('created_by', ''),
            notes=metadata.get('notes', '')
        )


class AlgorithmVersionManager:
    """Manages algorithm versions with database integration"""
    
    def __init__(self, context: ApplicationContext):
        self.context = context
        self.db = context.db
        self.logger = logging.getLogger(__name__)
        self._ensure_tables()
    
    def _ensure_tables(self):
        """Ensure required database tables exist"""
        try:
            # Test if table exists by running a simple query
            self.db.execute_query("SELECT COUNT(*) FROM ecod_schema.algorithm_version LIMIT 1")
        except Exception as e:
            self.logger.warning(f"Algorithm version table may not exist: {e}")
    
    def register_version(self, algorithm: AlgorithmVersion) -> int:
        """Register new algorithm version"""
        # Check if version already exists
        existing = self.get_version(algorithm.version_id)
        if existing:
            raise ValidationError(f"Algorithm version {algorithm.version_id} already exists")
        
        # Resolve parent version ID if specified
        parent_id = None
        if algorithm.parent_version:
            parent = self.get_version(algorithm.parent_version)
            if not parent:
                raise ValidationError(f"Parent version {algorithm.parent_version} not found")
            parent_id = parent.database_id
        
        # Prepare data for insertion
        db_data = algorithm.to_database_dict()
        
        # Insert into database
        query = """
        INSERT INTO ecod_schema.algorithm_version 
        (version_id, name, description, parent_version_id, status, config_data, created_by)
        VALUES (%s, %s, %s, %s, %s, %s, %s)
        RETURNING id
        """
        
        result = self.db.execute_query(query, (
            db_data['version_id'],
            db_data['name'],
            db_data['description'],
            parent_id,
            db_data['status'],
            db_data['config_data'],
            db_data['created_by']
        ))
        
        database_id = result[0][0]
        algorithm.database_id = database_id
        
        self.logger.info(f"Registered algorithm version {algorithm.version_id} with ID {database_id}")
        return database_id
    
    def get_version(self, version_id: str) -> Optional[AlgorithmVersion]:
        """Get algorithm version by ID"""
        query = """
        SELECT av.*, pv.version_id as parent_version_id_str
        FROM ecod_schema.algorithm_version av
        LEFT JOIN ecod_schema.algorithm_version pv ON av.parent_version_id = pv.id
        WHERE av.version_id = %s
        """
        
        results = self.db.execute_dict_query(query, (version_id,))
        
        if not results:
            return None
        
        row = results[0]
        # Replace parent_version_id (int) with parent_version_id_str (string)
        if row['parent_version_id_str']:
            row['parent_version_id'] = row['parent_version_id_str']
        
        return AlgorithmVersion.from_database_row(row)
    
    def list_versions(self, status: Optional[AlgorithmStatus] = None) -> List[AlgorithmVersion]:
        """List algorithm versions"""
        query = """
        SELECT av.*, pv.version_id as parent_version_id_str
        FROM ecod_schema.algorithm_version av
        LEFT JOIN ecod_schema.algorithm_version pv ON av.parent_version_id = pv.id
        """
        
        params = []
        if status:
            query += " WHERE av.status = %s"
            params.append(status.value)
        
        query += " ORDER BY av.created_at DESC"
        
        results = self.db.execute_dict_query(query, params)
        
        versions = []
        for row in results:
            if row['parent_version_id_str']:
                row['parent_version_id'] = row['parent_version_id_str']
            versions.append(AlgorithmVersion.from_database_row(row))
        
        return versions
    
    def promote_version(self, version_id: str, new_status: AlgorithmStatus) -> bool:
        """Promote algorithm version to new status"""
        algorithm = self.get_version(version_id)
        if not algorithm:
            raise ValidationError(f"Algorithm version {version_id} not found")
        
        # Validate promotion path
        if algorithm.status == AlgorithmStatus.DEPRECATED:
            raise ValidationError("Cannot promote deprecated algorithm")
        
        if new_status == AlgorithmStatus.PRODUCTION:
            # Additional validation for production promotion
            if algorithm.status not in [AlgorithmStatus.TESTING]:
                raise ValidationError("Can only promote to production from testing status")
        
        # Update status in database
        query = """
        UPDATE ecod_schema.algorithm_version 
        SET status = %s 
        WHERE version_id = %s
        """
        
        self.db.execute_query(query, (new_status.value, version_id))
        
        self.logger.info(f"Promoted algorithm {version_id} from {algorithm.status.value} to {new_status.value}")
        return True
    
    def get_production_version(self) -> Optional[AlgorithmVersion]:
        """Get the current production algorithm version"""
        versions = self.list_versions(status=AlgorithmStatus.PRODUCTION)
        
        if not versions:
            return None
        
        if len(versions) > 1:
            self.logger.warning(f"Multiple production algorithms found: {[v.version_id for v in versions]}")
            # Return the most recently created one
            return max(versions, key=lambda v: v.created_at)
        
        return versions[0]
    
    def start_algorithm_run(self, version_id: str, batch_id: Optional[int] = None,
                          run_type: str = "production") -> int:
        """Start tracking an algorithm run"""
        algorithm = self.get_version(version_id)
        if not algorithm:
            raise ValidationError(f"Algorithm version {version_id} not found")
        
        query = """
        INSERT INTO ecod_schema.algorithm_run 
        (version_id, batch_id, run_type)
        VALUES (%s, %s, %s)
        RETURNING id
        """
        
        result = self.db.execute_query(query, (
            algorithm.database_id, batch_id, run_type
        ))
        
        run_id = result[0][0]
        self.logger.info(f"Started algorithm run {run_id} for version {version_id}")
        return run_id
    
    def complete_algorithm_run(self, run_id: int, results_summary: Optional[Dict[str, Any]] = None) -> None:
        """Mark algorithm run as complete"""
        query = """
        UPDATE ecod_schema.algorithm_run 
        SET completed_at = CURRENT_TIMESTAMP,
            status = 'completed',
            results_summary = %s
        WHERE id = %s
        """
        
        summary_json = json.dumps(results_summary) if results_summary else None
        
        self.db.execute_query(query, (summary_json, run_id))
        
        self.logger.info(f"Completed algorithm run {run_id}")
    
    def export_version(self, version_id: str, output_path: str) -> None:
        """Export algorithm version to configuration file"""
        algorithm = self.get_version(version_id)
        if not algorithm:
            raise ValidationError(f"Algorithm version {version_id} not found")
        
        import yaml
        
        config_data = {
            'algorithm': {
                'version_id': algorithm.version_id,
                'name': algorithm.name,
                'description': algorithm.description,
                'parent_version_id': algorithm.parent_version,
                'status': algorithm.status.value,
                'created_by': algorithm.created_by,
                'notes': algorithm.notes
            },
            'domain_analysis': {
                'partition': algorithm.partition_config,
                'evidence_weights': algorithm.evidence_weights,
                'coverage_thresholds': algorithm.coverage_thresholds,
                'behavioral_flags': algorithm.behavioral_flags
            }
        }
        
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            if output_path.suffix.lower() == '.json':
                json.dump(config_data, f, indent=2, default=str)
            else:
                yaml.dump(config_data, f, default_flow_style=False)
        
        self.logger.info(f"Exported algorithm {version_id} to {output_path}")
    
    def import_version(self, config_path: str) -> int:
        """Import algorithm version from configuration file"""
        algorithm = AlgorithmVersion.from_config_file(config_path)
        return self.register_version(algorithm)
