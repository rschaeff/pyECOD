#!/usr/bin/env python3
"""
Pipeline Orchestration Service

Provides high-level coordination of the ECOD pipeline stages with a
deterministic, easy-to-understand workflow.
"""
import os
import logging
from datetime import datetime
from typing import Dict, Any, List, Optional
import uuid

from ecod.core.context import ApplicationContext
from ecod.exceptions import PipelineError
from ecod.pipelines.blast_pipeline import BlastPipeline
from ecod.pipelines.hhsearch_pipeline import HHSearchPipeline
from ecod.pipelines.hhsearch.service import HHSearchRegistrationService
from ecod.pipelines.domain_analysis.summary.service import DomainSummaryService
from ecod.pipelines.domain_analysis.partition.service import DomainPartitionService

from .models import (
    PipelineMode, PipelineStage, StageResult, PipelineRun,
    OrchestrationConfig
)
from .stage_manager import StageManager


class PipelineOrchestrationService:
    """
    Orchestrates the ECOD pipeline with a clear, deterministic workflow.
    
    The service coordinates:
    1. BLAST searches (chain and domain)
    2. HHSearch profile generation and searching
    3. HHSearch result registration
    4. Domain summary generation
    5. Domain partitioning
    6. Classification (if implemented)
    
    All proteins follow the same path - no complex routing decisions.
    """
    
    def __init__(self, context: ApplicationContext,
                 config: Optional[OrchestrationConfig] = None):
        """
        Initialize orchestration service.
        
        Args:
            context: Application context
            config: Optional orchestration configuration
        """
        self.context = context
        self.db = context.db
        self.config = config or OrchestrationConfig()
        self.logger = logging.getLogger(__name__)
        
        # Initialize pipeline components
        self.blast_pipeline = BlastPipeline(context)
        self.hhsearch_pipeline = HHSearchPipeline(context)
        self.hhsearch_registration = HHSearchRegistrationService(context)
        self.summary_service = DomainSummaryService(context)
        self.partition_service = DomainPartitionService(context)
        
        # Stage manager
        self.stage_manager = StageManager(context)
        
        # Track active runs
        self.active_runs: Dict[str, PipelineRun] = {}
        
        self.logger.info("PipelineOrchestrationService initialized")
    
    def run_pipeline(self, batch_id: int, mode: PipelineMode = PipelineMode.FULL,
                    force_restart: bool = False) -> PipelineRun:
        """
        Run the pipeline for a batch.
        
        Args:
            batch_id: Batch ID to process
            mode: Pipeline execution mode
            force_restart: Force restart from beginning
            
        Returns:
            PipelineRun with execution results
        """
        # Create run ID
        run_id = f"{batch_id}_{mode.value}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        
        self.logger.info(f"Starting pipeline run {run_id} for batch {batch_id} in {mode.value} mode")
        
        # Get batch info
        batch_info = self._get_batch_info(batch_id)
        if not batch_info:
            raise PipelineError(f"Batch {batch_id} not found")
        
        # Create pipeline run
        run = PipelineRun(
            run_id=run_id,
            batch_id=batch_id,
            mode=mode,
            total_proteins=batch_info['total_items']
        )
        
        self.active_runs[run_id] = run
        
        try:
            # Get stages to execute
            stages = self.stage_manager.get_stages_for_mode(mode)
            
            # Check for existing progress
            if not force_restart and self.config.checkpoint_enabled:
                completed_stages = self._get_completed_stages(batch_id, mode)
                self.logger.info(f"Found {len(completed_stages)} completed stages")
            else:
                completed_stages = []
            
            # Execute stages
            for stage in stages:
                # Skip if already completed
                if stage in completed_stages and not force_restart:
                    self.logger.info(f"Skipping completed stage: {stage.value}")
                    continue
                
                # Check dependencies
                if not self.stage_manager.can_execute_stage(stage, completed_stages):
                    self.logger.warning(f"Cannot execute {stage.value} - dependencies not met")
                    if not self.config.continue_on_error:
                        break
                    continue
                
                # Execute stage
                self.logger.info(f"Executing stage: {stage.value}")
                stage_result = self._execute_stage(batch_id, stage, batch_info)
                run.add_stage_result(stage_result)
                
                # Check for failure
                if not stage_result.success:
                    self.logger.error(f"Stage {stage.value} failed: {stage_result.error}")
                    if not self.config.continue_on_error:
                        break
                else:
                    completed_stages.append(stage)
            
            # Finalize run
            run.finalize()
            
            # Log summary
            self.logger.info(
                f"Pipeline run {run_id} completed in {run.duration:.2f}s. "
                f"Success: {run.success}, Stages: {len(completed_stages)}/{len(stages)}"
            )
            
            return run
            
        except Exception as e:
            self.logger.error(f"Pipeline error: {str(e)}", exc_info=True)
            run.finalize()
            raise
        finally:
            # Remove from active runs
            self.active_runs.pop(run_id, None)
    
    def _execute_stage(self, batch_id: int, stage: PipelineStage,
                      batch_info: Dict[str, Any]) -> StageResult:
        """Execute a single pipeline stage"""
        result = StageResult(
            stage=stage,
            success=False,
            start_time=datetime.now()
        )
        
        try:
            if stage == PipelineStage.INITIALIZATION:
                # Just validation
                result.success = True
                result.details = {'batch_info': batch_info}
                
            elif stage == PipelineStage.BLAST_CHAIN:
                job_ids = self.blast_pipeline.run_chain_blast(
                    batch_id, 
                    self.config.blast_batch_size
                )
                result.items_processed = len(job_ids)
                result.success = True
                result.details = {'job_ids': job_ids}
                
            elif stage == PipelineStage.BLAST_DOMAIN:
                job_ids = self.blast_pipeline.run_domain_blast(
                    batch_id,
                    self.config.blast_batch_size
                )
                result.items_processed = len(job_ids)
                result.success = True
                result.details = {'job_ids': job_ids}
                
            elif stage == PipelineStage.HHSEARCH_PROFILE:
                job_ids = self.hhsearch_pipeline.generate_profiles(
                    batch_id,
                    threads=self.config.hhsearch_threads,
                    memory=self.config.hhsearch_memory
                )
                result.items_processed = len(job_ids)
                result.success = True
                result.details = {'job_ids': job_ids}
                
            elif stage == PipelineStage.HHSEARCH_SEARCH:
                job_ids = self.hhsearch_pipeline.run_hhsearch(
                    batch_id,
                    threads=self.config.hhsearch_threads,
                    memory=self.config.hhsearch_memory
                )
                result.items_processed = len(job_ids)
                result.success = True
                result.details = {'job_ids': job_ids}
                
            elif stage == PipelineStage.HHSEARCH_REGISTRATION:
                reg_results = self.hhsearch_registration.register_batch(batch_id)
                result.items_processed = reg_results.registered
                result.items_failed = reg_results.failed
                result.success = reg_results.registered > 0
                result.details = {
                    'registered': reg_results.registered,
                    'skipped': reg_results.skipped,
                    'failed': reg_results.failed
                }
                
            elif stage == PipelineStage.DOMAIN_SUMMARY:
                # Get batch path
                batch_path = batch_info.get('base_path', '')
                
                # Process summaries
                summary_results = self.summary_service.process_batch(
                    batch_id,
                    batch_path,
                    blast_only=(self.config.mode == PipelineMode.BLAST_ONLY)
                )
                
                result.items_processed = summary_results.success_count
                result.items_failed = summary_results.error_count
                result.success = summary_results.success_count > 0
                result.details = summary_results.get_summary()
                
            elif stage == PipelineStage.DOMAIN_PARTITION:
                # Get batch path
                batch_path = batch_info.get('base_path', '')
                
                # Process partitions
                partition_results = self.partition_service.partition_batch(
                    batch_id,
                    batch_path
                )
                
                result.items_processed = partition_results.success_count
                result.items_failed = partition_results.error_count
                result.success = partition_results.success_count > 0
                result.details = partition_results.get_summary()
                
            elif stage == PipelineStage.CLASSIFICATION:
                # Placeholder for classification stage
                self.logger.info("Classification stage not yet implemented")
                result.success = True
                result.details = {'status': 'not_implemented'}
                
            elif stage == PipelineStage.COMPLETE:
                # Just marks completion
                result.success = True
                result.details = {'completed_at': datetime.now().isoformat()}
            
            else:
                result.error = f"Unknown stage: {stage.value}"
                
        except Exception as e:
            result.error = str(e)
            self.logger.error(f"Stage {stage.value} failed: {str(e)}", exc_info=True)
        
        finally:
            result.finalize(result.success, result.error)
        
        return result
    
    def get_batch_status(self, batch_id: int) -> Dict[str, Any]:
        """
        Get comprehensive status for a batch.
        
        Returns:
            Dictionary with batch and stage status
        """
        batch_info = self._get_batch_info(batch_id)
        if not batch_info:
            return {'error': f'Batch {batch_id} not found'}
        
        # Get status for all stages
        stage_status = {}
        for stage in PipelineStage:
            if stage in [PipelineStage.INITIALIZATION, PipelineStage.COMPLETE]:
                continue
            stage_status[stage.value] = self.stage_manager.get_stage_status(batch_id, stage)
        
        # Check for active runs
        active_run = None
        for run_id, run in self.active_runs.items():
            if run.batch_id == batch_id:
                active_run = run.get_summary()
                break
        
        return {
            'batch': batch_info,
            'stages': stage_status,
            'active_run': active_run
        }
    
    def resume_pipeline(self, batch_id: int, mode: Optional[PipelineMode] = None) -> PipelineRun:
        """
        Resume a pipeline from where it left off.
        
        Args:
            batch_id: Batch ID to resume
            mode: Override the mode (uses last mode if not specified)
            
        Returns:
            PipelineRun
        """
        # Determine mode
        if not mode:
            mode = self._get_last_mode(batch_id) or PipelineMode.FULL
        
        self.logger.info(f"Resuming pipeline for batch {batch_id} in {mode.value} mode")
        
        # Run with checkpoint enabled
        return self.run_pipeline(batch_id, mode, force_restart=False)
    
    def _get_batch_info(self, batch_id: int) -> Optional[Dict[str, Any]]:
        """Get batch information"""
        query = """
        SELECT id, batch_name, base_path, ref_version, 
               total_items, status, type
        FROM ecod_schema.batch
        WHERE id = %s
        """
        
        try:
            rows = self.db.execute_dict_query(query, (batch_id,))
            return rows[0] if rows else None
        except Exception as e:
            self.logger.error(f"Error getting batch info: {e}")
            return None
    
    def _get_completed_stages(self, batch_id: int, mode: PipelineMode) -> List[PipelineStage]:
        """Determine which stages have been completed for a batch"""
        completed = [PipelineStage.INITIALIZATION]  # Always include initialization
        
        # Check each stage's completion status
        for stage in self.stage_manager.get_stages_for_mode(mode):
            if stage == PipelineStage.INITIALIZATION:
                continue
                
            status = self.stage_manager.get_stage_status(batch_id, stage)
            if status['completed']:
                completed.append(stage)
        
        return completed
    
    def _get_last_mode(self, batch_id: int) -> Optional[PipelineMode]:
        """Get the last mode used for a batch (from metadata or guess from files)"""
        # Check if HHSearch files exist
        query = """
        SELECT COUNT(*) 
        FROM ecod_schema.process_file pf
        JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
        WHERE ps.batch_id = %s AND pf.file_type = 'hhr'
        """
        
        try:
            result = self.db.execute_query(query, (batch_id,))
            if result and result[0][0] > 0:
                return PipelineMode.FULL
            else:
                return PipelineMode.BLAST_ONLY
        except Exception:
            return None
    
    def get_run_history(self, batch_id: Optional[int] = None, 
                       limit: int = 10) -> List[Dict[str, Any]]:
        """Get pipeline run history (would need database table for persistence)"""
        # For now, just return active runs
        runs = []
        for run in self.active_runs.values():
            if batch_id is None or run.batch_id == batch_id:
                runs.append(run.get_summary())
        
        return sorted(runs, key=lambda x: x['start_time'], reverse=True)[:limit]


# Convenience functions

def create_orchestrator(config_path: Optional[str] = None,
                       orchestration_config: Optional[Dict[str, Any]] = None) -> PipelineOrchestrationService:
    """Create a pipeline orchestrator"""
    if config_path:
        context = ApplicationContext(config_path)
    else:
        config_path = os.environ.get('ECOD_CONFIG_PATH', 'config/config.yml')
        context = ApplicationContext(config_path)
    
    config = OrchestrationConfig(**orchestration_config) if orchestration_config else None
    
    return PipelineOrchestrationService(context, config)


def run_full_pipeline(batch_id: int, config_path: Optional[str] = None) -> Dict[str, Any]:
    """Run full pipeline for a batch"""
    orchestrator = create_orchestrator(config_path)
    run = orchestrator.run_pipeline(batch_id, PipelineMode.FULL)
    return run.get_summary()


def run_blast_only_pipeline(batch_id: int, config_path: Optional[str] = None) -> Dict[str, Any]:
    """Run BLAST-only pipeline for a batch"""
    orchestrator = create_orchestrator(config_path)
    run = orchestrator.run_pipeline(batch_id, PipelineMode.BLAST_ONLY)
    return run.get_summary()
