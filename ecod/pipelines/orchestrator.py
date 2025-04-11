# ecod/pipelines/orchestrator.py
import logging
from typing import Dict, Any, List, Optional
from pathlib import Path
from datetime import datetime

from ecod.db import DBManager
from ecod.jobs import JobManager
from ecod.core.context import ApplicationContext
from ecod.config import ConfigManager
from ..core.logging_config import LoggingManager
from ecod.exceptions import ConfigurationError, PipelineError
from .blast_pipeline import BlastPipeline
from .hhsearch_pipeline import HHSearchPipeline
from ecod.pipelines.domain_analysis.routing import ProcessingRouter

class PipelineOrchestrator:
    """Orchestrates different pipeline components for ECOD protein domain classification"""
    
    def __init__(self, config_path: Optional[str] = None):
        """Initialize orchestrator with configuration
        
        Args:
            config_path: Path to configuration file
        """
        try:
            # Initialize application context
            self.context = ApplicationContext(config_path)
            
            # Get references to resources from context
            self.db = self.context.db
            self.config = self.context.config_manager.config
            self.config_manager = self.context.config_manager
                
            # Setup logging
            self.logger = LoggingManager.get_logger("ecod.orchestrator")
            
            # Initialize managers
            self.db = DBManager(self.config_manager.get_db_config())
            from ecod.jobs.factory import create_job_manager
            self.job_manager = create_job_manager(self.config)
            
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
    
    def run_full_pipeline(self, batch_size: int = 500, max_chains: int = None) -> Dict[str, Any]:
        """Run the complete pipeline with enhanced error handling
        
        Args:
            batch_size: Maximum items per batch job
            max_chains: Maximum number of chains to process
            
        Returns:
            Pipeline status dictionary
        """
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
                    return pipeline_status
                    
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
                hhsearch_batch_id = self.hhsearch.create_batch(chains, batch_size)
                pipeline_status['hhsearch_batch_id'] = hhsearch_batch_id
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
                hhsearch_batch_id = None
            
            # 5. Generate profiles
            if hhsearch_batch_id:
                try:
                    self.logger.info("Generating HHblits profiles")
                    profile_job_ids = self.hhsearch.generate_profiles(hhsearch_batch_id, 8, "16G")
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
            if hhsearch_batch_id:
                self.logger.info(f"- HHSearch batch ID: {hhsearch_batch_id}, submitted {len(profile_job_ids)} jobs")
            
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
    
    def run_blast_pipeline(self, batch_id: int) -> bool:
        """Run only the BLAST pipeline for an existing batch"""
        self.logger.info(f"Running BLAST pipeline for batch {batch_id}")
        try:
            # Run chain BLAST
            chain_job_ids = self.blast.run_chain_blast(batch_id) or []
            self.logger.info(f"Submitted {len(chain_job_ids)} chain BLAST jobs")
            
            # Run domain BLAST
            domain_job_ids = self.blast.run_domain_blast(batch_id) or []
            self.logger.info(f"Submitted {len(domain_job_ids)} domain BLAST jobs")
            
            # If no jobs were submitted, check if this is because everything was already processed
            all_job_ids = chain_job_ids + domain_job_ids
            if not all_job_ids:
                # Check for existing results
                chain_results = self._count_existing_results(batch_id, "chain_blast_result")
                domain_results = self._count_existing_results(batch_id, "domain_blast_result")
                
                if chain_results > 0 or domain_results > 0:
                    self.logger.info(f"Found {chain_results} chain and {domain_results} domain BLAST results")
                    # Parse and store BLAST results
                    result_counts = self.blast.parse_and_store_blast_results(batch_id)
                    self.logger.info(f"BLAST pipeline stored {result_counts.get('domain_hits', 0)} domain hits and {result_counts.get('chain_hits', 0)} chain hits")
                    return True
                
                return False
                
            # Process completed immediately (local jobs)
            # Parse and store BLAST results
            result_counts = self.blast.parse_and_store_blast_results(batch_id)
            self.logger.info(f"BLAST pipeline stored {result_counts.get('domain_hits', 0)} domain hits and {result_counts.get('chain_hits', 0)} chain hits")
            
            return True
        
      except Exception as e:
         self.logger.error(f"Error running BLAST pipeline: {str(e)}", exc_info=True)
        return False
    
    def run_hhsearch_pipeline(self, batch_id: int) -> bool:
        """Run only the HHSearch pipeline for an existing batch
        
        Args:
            batch_id: Batch ID to process
            
        Returns:
            True if successful
        """
        self.logger.info(f"Running HHSearch pipeline for batch {batch_id}")
        try:
            # Generate profiles
            profile_job_ids = self.hhsearch.generate_profiles(batch_id)
            self.logger.info(f"Submitted {len(profile_job_ids)} profile generation jobs")
            
            return len(profile_job_ids) > 0
        except Exception as e:
            self.logger.error(f"Error running HHSearch pipeline: {str(e)}", exc_info=True)
            return False
    
    def run_domain_analysis(self, batch_id: int, blast_only: bool= False) -> bool:
        """Run domain analysis for an existing batch
        
        Args:
            batch_id: Batch ID to process
            blast_only: Whether to use only BLAST results (no HHsearch)
            
        Returns:
            True if successful
        """
        self.logger.info(f"Running domain analysis for batch {batch_id}")
        try:
            # Import domain analysis pipeline
            from .domain_analysis.pipeline import DomainAnalysisPipeline
            
            # Initialize domain analysis pipeline
            domain_pipeline = DomainAnalysisPipeline(self.config_manager.config_path)
            
            # Run pipeline
            return domain_pipeline.run_pipeline(batch_id, blast_only)
        except Exception as e:
            self.logger.error(f"Error running domain analysis: {str(e)}", exc_info=True)
            return False
    
    def run_full_pipeline_for_batch(self, batch_id: int) -> bool:
        """Run the full pipeline for an existing batch
        
        Args:
            batch_id: Batch ID to process
            
        Returns:
            True if successful
        """
        self.logger.info(f"Running full pipeline for batch {batch_id}")
        
        # Run BLAST
        blast_success = self.run_blast_pipeline(batch_id)
        if not blast_success:
            self.logger.warning(f"BLAST pipeline failed for batch {batch_id}")
        
        # Run HHSearch
        hhsearch_success = self.run_hhsearch_pipeline(batch_id)
        if not hhsearch_success:
            self.logger.warning(f"HHSearch pipeline failed for batch {batch_id}")
        
        # Run domain analysis
        domain_success = self.run_domain_analysis(batch_id)
        if not domain_success:
            self.logger.warning(f"Domain analysis failed for batch {batch_id}")
        
        # Return overall success
        return blast_success or hhsearch_success or domain_success

    def run_partitioned_pipeline(self, batch_id: int) -> Dict[str, Any]:
        """Run the pipeline using the adaptive routing strategy"""
        results = {
            "status": "starting",
            "paths": {}
        }
        
        if not blast_success:
            # Check if this is because everything was already processed
            existing_results = self.blast._count_existing_results(batch_id, "chain_blast_result")
            if existing_results > 0:
                self.logger.info(f"Found {existing_results} existing BLAST results, continuing with routing")
                blast_success = True
                
        if not blast_success:
            results["status"] = "blast_failed"
            return results
        
        # Parse and store BLAST results
        self.logger.info(f"Parsing BLAST results for batch {batch_id}")
        blast_results = self.blast.parse_and_store_blast_results(batch_id)
        results["blast_results"] = blast_results
        
        # Assign proteins to processing paths based on BLAST results
        paths = router.assign_processing_paths(batch_id)
        results["paths"] = {k: len(v) for k, v in paths.items()}
        
        # Rest of the method remains unchanged
        # ...

        return results
    
    def _count_existing_results(self, batch_id: int, file_type: str) -> int:
        """Count the number of existing result files for a batch"""
        query = """
            SELECT COUNT(DISTINCT ps.id)
            FROM ecod_schema.process_status ps
            JOIN ecod_schema.process_file pf ON ps.id = pf.process_id
            WHERE ps.batch_id = %s
            AND pf.file_type = %s
            AND pf.file_exists = TRUE
        """
        
        try:
            rows = self.db.execute_query(query, (batch_id, file_type))
            if rows and rows[0][0]:
                return rows[0][0]
        except Exception as e:
            self.logger.error(f"Error counting existing {file_type} files: {e}")
        
        return 0

    def _process_blast_only_proteins(self, batch_id: int, protein_ids: List[int]) -> bool:
        """Process proteins using only BLAST results
        
        Args:
            batch_id: Batch ID
            protein_ids: List of protein IDs to process
            
        Returns:
            True if successful
        """
        self.logger.info(f"Processing {len(protein_ids)} proteins with BLAST-only path")
        
        # Create domain summary but skip HHSearch
        from ecod.pipelines.domain_analysis.pipeline import DomainAnalysisPipeline
        domain_pipeline = DomainAnalysisPipeline(self.config_manager.config_path)
        
        # Process with blast_only=True to skip HHSearch
        success = domain_pipeline.process_proteins(batch_id, protein_ids, blast_only=True)
        
        return success

    def _process_full_pipeline_proteins(self, batch_id: int, protein_ids: List[int],
                                       threads: int = 8, memory: str = "16G") -> bool:
        """Process proteins using the full pipeline including HHSearch
        
        Args:
            batch_id: Batch ID
            protein_ids: List of protein IDs to process
            threads: Number of threads for HHSearch
            memory: Memory allocation for HHSearch
            
        Returns:
            True if successful
        """
        self.logger.info(f"Processing {len(protein_ids)} proteins with full pipeline")
        
        # 1. Run HHSearch for these proteins
        from ecod.pipelines.hhsearch_pipeline import HHSearchPipeline
        hhsearch = HHSearchPipeline(self.db, self.job_manager, self.config)
        
        # Submit HHSearch jobs with custom resource allocation
        hhsearch_job_ids = hhsearch.process_specific_proteins(
            batch_id, protein_ids, threads, memory
        )
        
        if not hhsearch_job_ids:
            self.logger.warning("No HHSearch jobs submitted")
            return False
        
        # 2. Wait for HHSearch jobs to complete
        completed = self.job_manager.wait_for_jobs(hhsearch_job_ids, timeout=3600)
        
        if not completed:
            self.logger.warning("Not all HHSearch jobs completed successfully")
        
        # 3. Run domain analysis with full pipeline
        from ecod.pipelines.domain_analysis.pipeline import DomainAnalysisPipeline
        domain_pipeline = DomainAnalysisPipeline(self.config_manager.config_path)
        
        # Process with blast_only=False to include HHSearch results
        success = domain_pipeline.process_proteins(batch_id, protein_ids, blast_only=False)
        
        return success

    def check_status(self, batch_id: Optional[int] = None) -> None:
        """Check status of all running jobs
        
        Args:
            batch_id: Optional batch ID to check
        """
        self.logger.info("Checking job status")
        
        # Check all job status
        completed, failed, running = self.job_manager.slurm_manager.check_all_jobs(batch_id)
        self.logger.info(f"Job status: {completed} completed, {failed} failed, {running} running")
        
        # Update batch completion status
        self._update_batch_status(batch_id)
    
    def _update_batch_status(self, batch_id: Optional[int] = None) -> None:
        """Update batch completion status
        
        Args:
            batch_id: Optional batch ID to update
        """
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