# ecod/pipelines/orchestrator.py
import logging
from typing import Dict, Any, List, Optional
from pathlib import Path
from datetime import datetime

from ecod.core.context import ApplicationContext
from ecod.exceptions import ConfigurationError, PipelineError
from ecod.pipelines.blast_pipeline import BlastPipeline
from ecod.pipelines.hhsearch_pipeline import HHSearchPipeline
from ecod.pipelines.domain_analysis.routing import ProcessingRouter

class PipelineOrchestrator:
    """Orchestrates different pipeline components for ECOD protein domain classification"""
    
    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize orchestrator with application context
        
        Args:
            config_path: Optional path to configuration file
        """
        try:
            # Initialize application context
            self.context = ApplicationContext(config_path)
            
            # Initialize logger
            self.logger = logging.getLogger("ecod.orchestrator")
            
            # Initialize pipeline components with shared context
            self.blast = BlastPipeline(self.context)
            self.hhsearch = HHSearchPipeline(self.context)
            
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
        
        # Initialize router
        from ecod.core.context import ApplicationContext
        router = ProcessingRouter(self.context)
        
        try:
            # Check if this batch has already been routed
            query = """
            SELECT COUNT(*) FROM ecod_schema.protein_processing_path
            WHERE batch_id = %s
            """
            count_rows = self.db.execute_query(query, (batch_id,))
            already_routed = count_rows and count_rows[0][0] > 0
            
            if already_routed:
                self.logger.info(f"Batch {batch_id} already has routing information - using existing paths")
                # Get existing path counts
                query = """
                SELECT path_type, COUNT(*) 
                FROM ecod_schema.protein_processing_path
                WHERE batch_id = %s
                GROUP BY path_type
                """
                path_rows = self.db.execute_dict_query(query, (batch_id,))
                paths = {}
                for row in path_rows:
                    path_type = row['path_type']
                    count = row['count']
                    self.logger.info(f"Found {count} proteins in {path_type} path")
                    
                    # Get protein IDs for this path
                    query = """
                    SELECT protein_id
                    FROM ecod_schema.protein_processing_path
                    WHERE batch_id = %s AND path_type = %s
                    """
                    protein_rows = self.db.execute_query(query, (batch_id, path_type))
                    paths[path_type] = [row[0] for row in protein_rows]
                
                results["paths"] = {k: len(v) for k, v in paths.items()}
            else:
                # Run BLAST for all proteins
                blast_success = self.run_blast_pipeline(batch_id)
                
                if not blast_success:
                    results["status"] = "blast_failed"
                    return results
                
                # DUPLICATE CALL
                # Parse and store BLAST results
                #self.logger.info(f"Parsing BLAST results for batch {batch_id}")
                #blast_results = self.blast.parse_and_store_blast_results(batch_id)
                #results["blast_results"] = blast_results
                
                # Assign proteins to processing paths based on BLAST results
                paths = router.assign_processing_paths(batch_id)
                results["paths"] = {k: len(v) for k, v in paths.items()}
            
            # Process each path (whether newly created or existing)
            for path_type, protein_ids in paths.items():
                if not protein_ids:
                    continue
                    
                self.logger.info(f"Processing {len(protein_ids)} proteins with {path_type} path")
                    
                if path_type == 'blast_only':
                    # Process proteins that can be classified with just BLAST results
                    self._process_blast_only_proteins(batch_id, protein_ids)
                elif path_type == 'full_pipeline':
                    # Process proteins that need HHSearch
                    self._process_full_pipeline_proteins(batch_id, protein_ids)
            
            # Set final status
            results["status"] = "completed"
            return results
            
        except Exception as e:
            self.logger.error(f"Error in partitioned pipeline: {str(e)}", exc_info=True)
            results["status"] = "error"
            results["error"] = str(e)
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

    def register_file(self, process_id, file_type, file_path, file_exists=True):
        """Register a file with proper duplicate handling
        
        Args:
            process_id: Process ID
            file_type: File type
            file_path: Path to file (relative to batch directory)
            file_exists: Whether file exists
        """
        try:
            # Calculate file size if file exists
            file_size = 0
            if file_exists and os.path.exists(file_path):
                file_size = os.path.getsize(file_path)
            
            # Check if record already exists
            query = """
            SELECT id FROM ecod_schema.process_file
            WHERE process_id = %s AND file_type = %s
            """
            
            existing = self.db.execute_query(query, (process_id, file_type))
            
            if existing:
                # Update existing record
                self.db.update(
                    "ecod_schema.process_file",
                    {
                        "file_path": file_path,
                        "file_exists": file_exists,
                        "file_size": file_size
                    },
                    "id = %s",
                    (existing[0][0],)
                )
                self.logger.debug(f"Updated existing {file_type} file record for process {process_id}")
            else:
                # Insert new record
                self.db.insert(
                    "ecod_schema.process_file",
                    {
                        "process_id": process_id,
                        "file_type": file_type,
                        "file_path": file_path,
                        "file_exists": file_exists,
                        "file_size": file_size
                    }
                )
                self.logger.debug(f"Created new {file_type} file record for process {process_id}")
        except Exception as e:
            self.logger.warning(f"Error registering {file_type} file for process {process_id}: {e}")

    def update_protein_processing_state(context, process_id, new_state, error_message=None):
        """Update processing state with proper tracking and history"""
        try:
            # Get current state
            current_state = context.db.execute_query(
                "SELECT current_stage, status FROM ecod_schema.process_status WHERE id = %s", 
                (process_id,)
            )[0]
            
            # Record state transition in history table
            context.db.insert(
                "ecod_schema.process_history", 
                {
                    "process_id": process_id,
                    "previous_state": current_state[0],
                    "new_state": new_state,
                    "transition_time": "NOW()",
                    "notes": error_message
                }
            )
            
            # Update current state
            context.db.update(
                "ecod_schema.process_status",
                {
                    "current_stage": new_state,
                    "status": "error" if error_message else "processing",
                    "error_message": error_message
                },
                "id = %s",
                (process_id,)
            )
        except Exception as e:
            logger.error(f"Failed to update processing state: {str(e)}")

    def identify_stalled_processes(context, batch_id, time_threshold_hours=24):
        """Identify processes that have been stuck in a processing state for too long"""
        query = """
        SELECT ps.id, p.pdb_id, p.chain_id, ps.current_stage, ps.status, 
               EXTRACT(EPOCH FROM (NOW() - ps.updated_at))/3600 as hours_since_update
        FROM ecod_schema.process_status ps
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        WHERE ps.batch_id = %s
        AND ps.status = 'processing'
        AND (NOW() - ps.updated_at) > interval '%s hours'
        """
        
        stalled = context.db.execute_dict_query(query, (batch_id, time_threshold_hours))
        return stalled

    def get_batch_progress_summary(context, batch_id):
        """Get a summary of batch processing progress"""
        query = """
        SELECT
            b.batch_name,
            b.type,
            b.total_items,
            COUNT(ps.id) as processed_items,
            SUM(CASE WHEN ps.status = 'success' THEN 1 ELSE 0 END) as success_count,
            SUM(CASE WHEN ps.status = 'error' THEN 1 ELSE 0 END) as error_count,
            SUM(CASE WHEN ps.status = 'processing' THEN 1 ELSE 0 END) as processing_count,
            COUNT(DISTINCT CASE WHEN pf.file_type = 'fasta' THEN ps.id END) as fasta_count,
            COUNT(DISTINCT CASE WHEN pf.file_type = 'a3m' THEN ps.id END) as profile_count,
            COUNT(DISTINCT CASE WHEN pf.file_type = 'hhr' THEN ps.id END) as hhsearch_count
        FROM
            ecod_schema.batch b
        LEFT JOIN
            ecod_schema.process_status ps ON b.id = ps.batch_id
        LEFT JOIN
            ecod_schema.process_file pf ON ps.id = pf.process_id
        WHERE
            b.id = %s
        GROUP BY
            b.id, b.batch_name, b.type, b.total_items
        """
        
        summary = context.db.execute_dict_query(query, (batch_id,))
        return summary[0] if summary else None