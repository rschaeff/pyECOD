# ecod/pipelines/orchestrator.py
import logging
from typing import Dict, Any, List, Optional
from pathlib import Path
from datetime import datetime

from ..core.db_manager import DBManager
from ..core.job_manager import JobManager
from ..core.config import ConfigManager
from ..core.logging_config import LoggingManager
from .blast_pipeline import BlastPipeline
from .hhsearch_pipeline import HHSearchPipeline

class PipelineOrchestrator:
    def __init__(self, config_path: Optional[str] = None):
        try:
            # Initialize configuration
            self.config_manager = ConfigManager(config_path)
            self.config = self.config_manager.config
            
            # Setup logging
            self.logger = LoggingManager.get_logger("ecod.orchestrator")
            
            # Initialize managers
            self.db = DBManager(self.config_manager.get_db_config())
            self.job_manager = JobManager(self.config)
            
            # Initialize pipeline components
            self.blast = BlastPipeline(self.db, self.job_manager, self.config)
            self.hhsearch = HHSearchPipeline(self.db, self.job_manager, self.config)
            
            self.logger.info("Pipeline orchestrator initialized successfully")
            
        except ConfigurationError as e:
            # Log the error and re-raise
            if hasattr(self, 'logger'):
                self.logger.error(f"Configuration error: {str(e)}")
            else:
                # Fallback logging if logger not yet initialized
                logging.error(f"Configuration error: {str(e)}")
            raise
        except Exception as e:
            # Log any other initialization errors
            error_msg = f"Error initializing pipeline orchestrator: {str(e)}"
            if hasattr(self, 'logger'):
                self.logger.error(error_msg, exc_info=True)
            else:
                logging.error(error_msg, exc_info=True)
            raise PipelineError(error_msg) from e
        
        def run_full_pipeline(self, batch_size: int = 500, max_chains: int = None) -> None:
        """Run the complete pipeline with enhanced error handling"""
        self.logger.info("Starting full pipeline")
        pipeline_status = {
            'started': True,
            'status': 'starting',
            'steps_completed': [],
            'steps_failed': [],
            'batch_id': None
        }
        
        try:
            # 1. Get unclassified chains
            self.logger.info("Getting unclassified chains")
            try:
                chains = self.blast.get_unclassified_chains(max_chains)
                if not chains:
                    self.logger.info("No unclassified chains found, pipeline completed")
                    pipeline_status['status'] = 'completed_no_chains'
                    return
                    
                self.logger.info(f"Found {len(chains)} unclassified chains")
                pipeline_status['chains_found'] = len(chains)
            except Exception as e:
                error_msg = f"Error getting unclassified chains: {str(e)}"
                self.logger.error(error_msg, exc_info=True)
                pipeline_status['status'] = 'failed'
                pipeline_status['error'] = error_msg
                pipeline_status['step_failed'] = 'get_chains'
                raise PipelineError(error_msg) from e
            
            # 2. Create batch for BLAST
            try:
                self.logger.info(f"Creating BLAST batch for {len(chains)} chains")
                blast_batch = self.blast.create_batch(chains, batch_size)
                pipeline_status['batch_id'] = blast_batch.id
                pipeline_status['steps_completed'].append('create_batch')
            except Exception as e:
                error_msg = f"Error creating BLAST batch: {str(e)}"
                self.logger.error(error_msg, exc_info=True)
                pipeline_status['status'] = 'failed'
                pipeline_status['error'] = error_msg
                pipeline_status['step_failed'] = 'create_batch'
                raise PipelineError(error_msg) from e
            
            # 3. Run BLAST searches
            try:
                self.logger.info("Running chain-wise BLAST")
                chain_job_ids = self.blast.run_chain_blast(blast_batch.id, batch_size)
                pipeline_status['steps_completed'].append('chain_blast')
                pipeline_status['chain_jobs'] = len(chain_job_ids)
            except Exception as e:
                error_msg = f"Error running chain BLAST: {str(e)}"
                self.logger.error(error_msg, exc_info=True)
                pipeline_status['status'] = 'failed'
                pipeline_status['error'] = error_msg
                pipeline_status['step_failed'] = 'chain_blast'
                
                # Continue to next step despite error
                self.logger.info("Continuing pipeline despite chain BLAST error")
                pipeline_status['steps_failed'].append('chain_blast')
                chain_job_ids = []
            
            try:
                self.logger.info("Running domain-wise BLAST")
                domain_job_ids = self.blast.run_domain_blast(blast_batch.id, batch_size)
                pipeline_status['steps_completed'].append('domain_blast')
                pipeline_status['domain_jobs'] = len(domain_job_ids)
            except Exception as e:
                error_msg = f"Error running domain BLAST: {str(e)}"
                self.logger.error(error_msg, exc_info=True)
                pipeline_status['status'] = 'failed'
                pipeline_status['error'] = error_msg
                pipeline_status['step_failed'] = 'domain_blast'
                
                # Continue to next step despite error
                self.logger.info("Continuing pipeline despite domain BLAST error")
                pipeline_status['steps_failed'].append('domain_blast')
                domain_job_ids = []
            
            # 4. Create batch for HHSearch
            try:
                self.logger.info("Creating HHSearch batch")
                hhsearch_batch = self.hhsearch.create_batch(chains, batch_size)
                pipeline_status['hhsearch_batch_id'] = hhsearch_batch.id
                pipeline_status['steps_completed'].append('create_hhsearch_batch')
            except Exception as e:
                error_msg = f"Error creating HHSearch batch: {str(e)}"
                self.logger.error(error_msg, exc_info=True)
                pipeline_status['status'] = 'failed'
                pipeline_status['error'] = error_msg
                pipeline_status['step_failed'] = 'create_hhsearch_batch'
                
                # Continue to next step despite error
                self.logger.info("Continuing pipeline despite HHSearch batch creation error")
                pipeline_status['steps_failed'].append('create_hhsearch_batch')
            
            # 5. Generate profiles
            try:
                self.logger.info("Generating HHblits profiles")
                profile_job_ids = self.hhsearch.generate_profiles(hhsearch_batch.id, 8, "16G")
                pipeline_status['steps_completed'].append('generate_profiles')
                pipeline_status['profile_jobs'] = len(profile_job_ids)
            except Exception as e:
                error_msg = f"Error generating HHblits profiles: {str(e)}"
                self.logger.error(error_msg, exc_info=True)
                pipeline_status['status'] = 'failed'
                pipeline_status['error'] = error_msg
                pipeline_status['step_failed'] = 'generate_profiles'
                
                # Add to failed steps and continue
                pipeline_status['steps_failed'].append('generate_profiles')
            
            # Set final status
            if not pipeline_status['steps_failed']:
                pipeline_status['status'] = 'completed'
                self.logger.info(f"Pipeline started successfully!")
            else:
                pipeline_status['status'] = 'completed_with_errors'
                self.logger.warning(f"Pipeline completed with errors in steps: {', '.join(pipeline_status['steps_failed'])}")
            
            self.logger.info(f"- BLAST batch ID: {blast_batch.id}, submitted {len(chain_job_ids) + len(domain_job_ids)} jobs")
            if 'hhsearch_batch_id' in pipeline_status:
                self.logger.info(f"- HHSearch batch ID: {hhsearch_batch.id}, submitted {len(profile_job_ids)} jobs")
            
            return pipeline_status
            
        except Exception as e:
            error_msg = f"Unhandled error in pipeline: {str(e)}"
            self.logger.error(error_msg, exc_info=True)
            pipeline_status['status'] = 'failed'
            pipeline_status['error'] = error_msg
            pipeline_status['step_failed'] = 'unhandled'
            
            # Log pipeline status for debugging
            self.logger.error(f"Pipeline status at failure: {pipeline_status}")
            
            raise PipelineError(error_msg) from e
        
    # Update the check_status method in the orchestrator
    def check_status(self, batch_id: Optional[int] = None) -> None:
        """Check status of all running jobs"""
        self.logger.info("Checking job status")
        
        # Check BLAST jobs
        self.blast.check_job_status(batch_id)
        
        # Check HHSearch jobs
        self.hhsearch.check_status(batch_id)
        
        # Update batch completion status
        self._update_batch_status(batch_id)
        
        # Parse HHSearch results for completed jobs
        self._parse_completed_hhsearch_results(batch_id)
        
    def _parse_completed_hhsearch_results(self, batch_id: Optional[int] = None) -> None:
        """Parse and generate summaries for completed HHSearch results"""
        query = """
            SELECT 
                ps.id, p.pdb_id, p.chain_id, ps.relative_path, b.base_path
            FROM 
                ecod_schema.process_status ps
            JOIN
                ecod_schema.protein p ON ps.protein_id = p.id
            JOIN
                ecod_schema.batch b ON ps.batch_id = b.id
            JOIN
                ecod_schema.process_file pf ON ps.id = pf.process_id
            WHERE 
                ps.current_stage = 'hhsearch_complete'
                AND ps.status = 'success'
                AND pf.file_type = 'hhr'
                AND pf.file_exists = TRUE
                AND NOT EXISTS (
                    SELECT 1 FROM ecod_schema.process_file pf2
                    WHERE pf2.process_id = ps.id AND pf2.file_type = 'hh_summ_xml'
                )
        """
        
        if batch_id:
            query += " AND ps.batch_id = %s"
            rows = self.db.execute_dict_query(query, (batch_id,))
        else:
            rows = self.db.execute_dict_query(query)
        
        if rows:
            self.logger.info(f"Parsing {len(rows)} completed HHSearch results")
        self.hhsearch.parse_results(rows)
        
    def _update_batch_status(self, batch_id: Optional[int] = None) -> None:
        """Update batch completion status"""
        query = """
            UPDATE ecod_schema.batch
            SET completed_items = (
                SELECT COUNT(*) 
                FROM ecod_schema.process_status
                WHERE batch_id = ecod_schema.batch.id 
                AND status IN ('success', 'skipped')
            ),
            status = CASE 
                WHEN (SELECT COUNT(*) FROM ecod_schema.process_status 
                      WHERE batch_id = ecod_schema.batch.id 
                      AND status NOT IN ('success', 'skipped')) = 0 
                THEN 'completed' 
                ELSE status 
            END,
            completed_at = CASE 
                WHEN (SELECT COUNT(*) FROM ecod_schema.process_status 
                      WHERE batch_id = ecod_schema.batch.id 
                      AND status NOT IN ('success', 'skipped')) = 0 
                THEN CURRENT_TIMESTAMP 
                ELSE completed_at 
            END
        """
        
        if batch_id:
            query += " WHERE id = %s"
            self.db.execute_query(query, (batch_id,))
        else:
            query += " WHERE status != 'completed'"
            self.db.execute_query(query)
            
        # Get summary
        if batch_id:
            batch_query = """
                SELECT 
                    b.id, b.batch_name, b.status, b.total_items, b.completed_items
                FROM 
                    ecod_schema.batch b
                WHERE 
                    b.id = %s
            """
            rows = self.db.execute_dict_query(batch_query, (batch_id,))
        else:
            batch_query = """
                SELECT 
                    b.id, b.batch_name, b.status, b.total_items, b.completed_items
                FROM 
                    ecod_schema.batch b
                WHERE 
                    b.status != 'completed'
                ORDER BY 
                    b.id
            """
            rows = self.db.execute_dict_query(batch_query)
            
        for row in rows:
            self.logger.info(f"Batch {row['id']} ({row['batch_name']}): {row['completed_items']}/{row['total_items']} completed, status: {row['status']}")