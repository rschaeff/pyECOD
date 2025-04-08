#!/usr/bin/env python3
"""
Database-aware job manager for the ECOD pipeline.
"""
import logging
from typing import Dict, Any, List, Tuple, Optional, Union

from ecod.config import ConfigManager
from ecod.exceptions import JobSubmissionError, JobExecutionError
from ecod.db.manager import DBManager
from .base import JobManager
from .factory import create_job_manager
from .job import Job, JobItem

logger = logging.getLogger("ecod.jobs.db_job_manager")

class DatabaseJobManager:
    """Job manager with database integration"""
    
    def __init__(self, config_or_path: Union[Dict[str, Any], str]):
        """Initialize with configuration
        
        Args:
            config_or_path: Configuration dictionary or path to config file
        """
        # Load configuration
        if isinstance(config_or_path, str):
            self.config_manager = ConfigManager(config_or_path)
            self.config = self.config_manager.config
            self.db = DBManager(self.config_manager.get_db_config())
        else:
            self.config = config_or_path
            self.db = DBManager(self.config.get('database', {}))
        
        # Create job manager
        manager_type = self.config.get('job_manager', {}).get('type', 'local')
        self.job_manager = create_job_manager(self.config, manager_type)
        
        logger.info(f"Initialized database job manager with {manager_type} job manager")
    
    def create_job_script(self, commands: List[str], job_name: str, 
                         output_dir: str, **options) -> str:
        """Create a job script
        
        Args:
            commands: List of commands to run
            job_name: Name of the job
            output_dir: Directory for output files
            **options: Additional job options
            
        Returns:
            Path to the created job script
        """
        return self.job_manager.create_job_script(commands, job_name, output_dir, **options)
    
    def submit_job(self, script_path: str, batch_id: Optional[int] = None, 
                  job_type: str = "generic", process_ids: Optional[List[int]] = None) -> Optional[str]:
        """Submit a job and track in database
        
        Args:
            script_path: Path to job script
            batch_id: Optional batch ID for tracking
            job_type: Type of job
            process_ids: Optional list of process IDs to associate with job
            
        Returns:
            Job ID if successful, None otherwise
        """
        # Submit job to execution system
        slurm_job_id = self.job_manager.submit_job(script_path)
        
        if not slurm_job_id:
            logger.error(f"Failed to submit job with script {script_path}")
            return None
        
        # Record in database if batch_id provided
        if batch_id:
            try:
                # Create job record
                job_db_id = self.db.insert(
                    "ecod_schema.job",
                    {
                        "batch_id": batch_id,
                        "job_type": job_type,
                        "slurm_job_id": slurm_job_id,
                        "status": "submitted",
                        "items_count": len(process_ids) if process_ids else 1
                    },
                    "id"
                )
                
                # Record job items if process_ids provided
                if process_ids and job_db_id:
                    for process_id in process_ids:
                        self.db.insert(
                            "ecod_schema.job_item",
                            {
                                "job_id": job_db_id,
                                "process_id": process_id,
                                "status": "pending"
                            }
                        )
                        
                        # Update process status
                        self.db.update(
                            "ecod_schema.process_status",
                            {
                                "status": "processing"
                            },
                            "id = %s",
                            (process_id,)
                        )
                
                logger.info(f"Recorded job {slurm_job_id} in database with ID {job_db_id}")
            
            except Exception as e:
                logger.error(f"Error recording job in database: {str(e)}")
        
        return slurm_job_id
    
    def check_job_status(self, job_id: str) -> str:
        """Check job status
        
        Args:
            job_id: Job ID
            
        Returns:
            Status string
        """
        return self.job_manager.check_job_status(job_id)
    
    def cancel_job(self, job_id: str) -> bool:
        """Cancel a job
        
        Args:
            job_id: Job ID
            
        Returns:
            True if successful
        """
        return self.job_manager.cancel_job(job_id)
    
    def check_all_jobs(self, batch_id: Optional[int] = None) -> Tuple[int, int, int]:
        """Check status of all jobs and update database
        
        Args:
            batch_id: Optional batch ID filter
            
        Returns:
            Tuple of (completed, failed, running) counts
        """
        # Get pending jobs from database
        query = """
        SELECT 
            j.id, j.slurm_job_id, j.job_type, j.batch_id, j.status
        FROM 
            ecod_schema.job j
        WHERE 
            j.status = 'submitted'
        """
        
        if batch_id:
            query += " AND j.batch_id = %s"
            rows = self.db.execute_dict_query(query, (batch_id,))
        else:
            rows = self.db.execute_dict_query(query)
        
        completed = 0
        failed = 0
        running = 0
        
        for row in rows:
            job_id = row['id']
            slurm_job_id = row['slurm_job_id']
            
            # Check status in execution system
            status = self.job_manager.check_job_status(slurm_job_id)
            
            # Map status to database status
            db_status = "submitted"  # Default (still running)
            
            if status in ["COMPLETED", "COMPLETING"]:
                db_status = "completed"
                completed += 1
            elif status in ["FAILED", "TIMEOUT", "CANCELLED", "NODE_FAIL"]:
                db_status = "failed"
                failed += 1
            else:
                running += 1
            
            # Only update database if status changed
            if db_status != "submitted":
                self.db.update(
                    "ecod_schema.job",
                    {
                        "status": db_status,
                    },
                    "id = %s",
                    (job_id,)
                )
                
                # Update job items
                self._update_job_items(job_id, db_status)
        
        logger.info(f"Checked {len(rows)} jobs: {completed} completed, {failed} failed, {running} running")
        return completed, failed, running
    
    def _update_job_items(self, job_id: int, status: str) -> None:
        """Update status of items in a job
        
        Args:
            job_id: Job ID
            status: Job status
        """
        # Map job status to item status
        item_status = "pending"
        process_status = "processing"
        
        if status == "completed":
            item_status = "completed"
            process_status = "success"
        elif status == "failed":
            item_status = "failed"
            process_status = "error"
        
        # Update job items
        self.db.update(
            "ecod_schema.job_item",
            {
                "status": item_status
            },
            "job_id = %s",
            (job_id,)
        )
        
        # Get process IDs
        query = """
        SELECT process_id FROM ecod_schema.job_item WHERE job_id = %s
        """
        rows = self.db.execute_query(query, (job_id,))
        
        # Update process status
        for row in rows:
            process_id = row[0]
            self.db.update(
                "ecod_schema.process_status",
                {
                    "status": process_status
                },
                "id = %s",
                (process_id,)
            )
    
    def create_batch_jobs(self, items: List[Tuple[Any, str]], batch_size: int, 
                         job_template: Dict[str, Any], batch_id: Optional[int] = None) -> List[str]:
        """Create batch jobs
        
        Args:
            items: List of (item_id, item_path) tuples
            batch_size: Number of items per batch
            job_template: Template with job configuration
            batch_id: Optional batch ID for tracking
            
        Returns:
            List of job IDs
        """
        # Create batch jobs
        jobs = self.job_manager.create_batch_jobs(items, batch_size, job_template)
        
        job_ids = []
        for job in jobs:
            # Extract process IDs from items
            process_ids = [item_id for item_id, _ in job['items']]
            
            # Submit job
            job_id = self.submit_job(
                job['script_path'],
                batch_id,
                job_template.get('job_type', 'batch'),
                process_ids
            )
            
            if job_id:
                job_ids.append(job_id)
        
        return job_ids