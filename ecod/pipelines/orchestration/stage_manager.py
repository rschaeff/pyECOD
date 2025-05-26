#!/usr/bin/env python3
"""
Pipeline stage management
"""
import logging
from datetime import datetime
from typing import Dict, Any, Optional, List

from ecod.core.context import ApplicationContext
from .models import PipelineStage, StageResult, PipelineMode


class StageManager:
    """Manages pipeline stage execution and dependencies"""
    
    # Stage dependencies
    STAGE_DEPENDENCIES = {
        PipelineStage.BLAST_CHAIN: [PipelineStage.INITIALIZATION],
        PipelineStage.BLAST_DOMAIN: [PipelineStage.INITIALIZATION],
        PipelineStage.HHSEARCH_PROFILE: [PipelineStage.BLAST_CHAIN, PipelineStage.BLAST_DOMAIN],
        PipelineStage.HHSEARCH_SEARCH: [PipelineStage.HHSEARCH_PROFILE],
        PipelineStage.HHSEARCH_REGISTRATION: [PipelineStage.HHSEARCH_SEARCH],
        PipelineStage.DOMAIN_SUMMARY: [PipelineStage.BLAST_DOMAIN],  # Can run after BLAST
        PipelineStage.DOMAIN_PARTITION: [PipelineStage.DOMAIN_SUMMARY],
        PipelineStage.CLASSIFICATION: [PipelineStage.DOMAIN_PARTITION]
    }
    
    # Stages by mode
    MODE_STAGES = {
        PipelineMode.FULL: [
            PipelineStage.INITIALIZATION,
            PipelineStage.BLAST_CHAIN,
            PipelineStage.BLAST_DOMAIN,
            PipelineStage.HHSEARCH_PROFILE,
            PipelineStage.HHSEARCH_SEARCH,
            PipelineStage.HHSEARCH_REGISTRATION,
            PipelineStage.DOMAIN_SUMMARY,
            PipelineStage.DOMAIN_PARTITION,
            PipelineStage.CLASSIFICATION,
            PipelineStage.COMPLETE
        ],
        PipelineMode.BLAST_ONLY: [
            PipelineStage.INITIALIZATION,
            PipelineStage.BLAST_CHAIN,
            PipelineStage.BLAST_DOMAIN,
            PipelineStage.DOMAIN_SUMMARY,
            PipelineStage.DOMAIN_PARTITION,
            PipelineStage.CLASSIFICATION,
            PipelineStage.COMPLETE
        ],
        PipelineMode.ANALYSIS_ONLY: [
            PipelineStage.INITIALIZATION,
            PipelineStage.DOMAIN_SUMMARY,
            PipelineStage.DOMAIN_PARTITION,
            PipelineStage.CLASSIFICATION,
            PipelineStage.COMPLETE
        ]
    }
    
    def __init__(self, context: ApplicationContext):
        self.context = context
        self.logger = logging.getLogger(__name__)
    
    def get_stages_for_mode(self, mode: PipelineMode) -> List[PipelineStage]:
        """Get stages to execute for a given mode"""
        return self.MODE_STAGES.get(mode, [])
    
    def can_execute_stage(self, stage: PipelineStage, 
                         completed_stages: List[PipelineStage]) -> bool:
        """Check if a stage can be executed based on dependencies"""
        dependencies = self.STAGE_DEPENDENCIES.get(stage, [])
        return all(dep in completed_stages for dep in dependencies)
    
    def get_stage_status(self, batch_id: int, stage: PipelineStage) -> Dict[str, Any]:
        """Get current status of a stage for a batch"""
        status = {
            'stage': stage.value,
            'ready': False,
            'completed': False,
            'in_progress': False,
            'failed': False,
            'skipped': False,
            'stats': {}
        }
        
        # Check stage-specific status in database
        if stage == PipelineStage.BLAST_CHAIN:
            status['stats'] = self._get_blast_status(batch_id, 'chain')
        elif stage == PipelineStage.BLAST_DOMAIN:
            status['stats'] = self._get_blast_status(batch_id, 'domain')
        elif stage in [PipelineStage.HHSEARCH_PROFILE, PipelineStage.HHSEARCH_SEARCH]:
            status['stats'] = self._get_hhsearch_status(batch_id)
        elif stage == PipelineStage.DOMAIN_SUMMARY:
            status['stats'] = self._get_summary_status(batch_id)
        elif stage == PipelineStage.DOMAIN_PARTITION:
            status['stats'] = self._get_partition_status(batch_id)
        
        # Determine overall status
        stats = status['stats']
        if stats.get('total', 0) > 0:
            if stats.get('completed', 0) == stats['total']:
                status['completed'] = True
            elif stats.get('in_progress', 0) > 0:
                status['in_progress'] = True
            elif stats.get('failed', 0) == stats['total']:
                status['failed'] = True
        
        return status
    
    def _get_blast_status(self, batch_id: int, blast_type: str) -> Dict[str, int]:
        """Get BLAST processing status"""
        file_type = f"{blast_type}_blast_result"
        
        query = """
        SELECT 
            COUNT(*) as total,
            COUNT(CASE WHEN pf.id IS NOT NULL THEN 1 END) as completed
        FROM ecod_schema.process_status ps
        LEFT JOIN ecod_schema.process_file pf ON (
            ps.id = pf.process_id AND 
            pf.file_type = %s AND 
            pf.file_exists = TRUE
        )
        WHERE ps.batch_id = %s
        """
        
        try:
            result = self.context.db.execute_dict_query(query, (file_type, batch_id))
            if result:
                return {
                    'total': result[0]['total'],
                    'completed': result[0]['completed'],
                    'pending': result[0]['total'] - result[0]['completed']
                }
        except Exception as e:
            self.logger.error(f"Error getting BLAST status: {e}")
        
        return {'total': 0, 'completed': 0, 'pending': 0}
    
    def _get_hhsearch_status(self, batch_id: int) -> Dict[str, int]:
        """Get HHSearch processing status"""
        query = """
        SELECT 
            COUNT(*) as total,
            COUNT(CASE WHEN pf_a3m.id IS NOT NULL THEN 1 END) as profiles_completed,
            COUNT(CASE WHEN pf_hhr.id IS NOT NULL THEN 1 END) as searches_completed,
            COUNT(CASE WHEN pf_xml.id IS NOT NULL THEN 1 END) as registered
        FROM ecod_schema.process_status ps
        LEFT JOIN ecod_schema.process_file pf_a3m ON (
            ps.id = pf_a3m.process_id AND 
            pf_a3m.file_type = 'a3m' AND 
            pf_a3m.file_exists = TRUE
        )
        LEFT JOIN ecod_schema.process_file pf_hhr ON (
            ps.id = pf_hhr.process_id AND 
            pf_hhr.file_type = 'hhr' AND 
            pf_hhr.file_exists = TRUE
        )
        LEFT JOIN ecod_schema.process_file pf_xml ON (
            ps.id = pf_xml.process_id AND 
            pf_xml.file_type = 'hh_xml' AND 
            pf_xml.file_exists = TRUE
        )
        WHERE ps.batch_id = %s
        """
        
        try:
            result = self.context.db.execute_dict_query(query, (batch_id,))
            if result:
                return {
                    'total': result[0]['total'],
                    'profiles_completed': result[0]['profiles_completed'],
                    'searches_completed': result[0]['searches_completed'],
                    'registered': result[0]['registered']
                }
        except Exception as e:
            self.logger.error(f"Error getting HHSearch status: {e}")
        
        return {'total': 0, 'profiles_completed': 0, 'searches_completed': 0, 'registered': 0}
    
    def _get_summary_status(self, batch_id: int) -> Dict[str, int]:
        """Get domain summary status"""
        query = """
        SELECT 
            COUNT(*) as total,
            COUNT(CASE WHEN pf.id IS NOT NULL THEN 1 END) as completed
        FROM ecod_schema.process_status ps
        LEFT JOIN ecod_schema.process_file pf ON (
            ps.id = pf.process_id AND 
            pf.file_type = 'domain_summary' AND 
            pf.file_exists = TRUE
        )
        WHERE ps.batch_id = %s
        """
        
        try:
            result = self.context.db.execute_dict_query(query, (batch_id,))
            if result:
                return {
                    'total': result[0]['total'],
                    'completed': result[0]['completed'],
                    'pending': result[0]['total'] - result[0]['completed']
                }
        except Exception as e:
            self.logger.error(f"Error getting summary status: {e}")
        
        return {'total': 0, 'completed': 0, 'pending': 0}
    
    def _get_partition_status(self, batch_id: int) -> Dict[str, int]:
        """Get domain partition status"""
        query = """
        SELECT 
            COUNT(*) as total,
            COUNT(CASE WHEN pf.id IS NOT NULL THEN 1 END) as completed
        FROM ecod_schema.process_status ps
        LEFT JOIN ecod_schema.process_file pf ON (
            ps.id = pf.process_id AND 
            pf.file_type = 'domain_partition' AND 
            pf.file_exists = TRUE
        )
        WHERE ps.batch_id = %s
        """
        
        try:
            result = self.context.db.execute_dict_query(query, (batch_id,))
            if result:
                return {
                    'total': result[0]['total'],
                    'completed': result[0]['completed'],
                    'pending': result[0]['total'] - result[0]['completed']
                }
        except Exception as e:
            self.logger.error(f"Error getting partition status: {e}")
        
        return {'total': 0, 'completed': 0, 'pending': 0}
