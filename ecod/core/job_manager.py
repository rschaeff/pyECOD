#!/usr/bin/env python3
"""
Job manager for the ECOD pipeline
Provides a simplified interface to the enhanced SlurmJobManager
"""

import os
import logging
from typing import Dict, Any, List, Tuple, Optional, Union

from .slurm_job_manager import SlurmJobManager
from .exceptions import JobSubmissionError, JobExecutionError

class JobManager:
    """Job manager for the ECOD pipeline
    
    This is a wrapper around SlurmJobManager that maintains backward compatibility
    with the existing codebase while providing enhanced Slurm functionality.
    """
    
    def __init__(self, config: Union[Dict[str, Any], str]):
        """Initialize job manager
        
        Args:
            config: Configuration dictionary or path to config file
        """
        self.logger = logging.getLogger("ecod.job_manager")
        self.slurm_manager = SlurmJobManager(config)
        self.config = self.slurm_manager.config
    
    def create_job_script(self, commands: List[str], job_name: str, output_dir: str, 
                         threads: int = 8, memory: str = "8G", time: str = "24:00:00") -> str:
        """Create a job script
        
        Args:
            commands: List of commands to include in script
            job_name: Name of the job
            output_dir: Directory for output files
            threads: Number of threads/CPUs
            memory: Memory allocation
            time: Time limit
            
        Returns:
            Path to created job script
        """
        return self.slurm_manager.create_job_script(
            commands, job_name, output_dir, threads, memory, time
        )
        
    def submit_job(self, script_path: str) -> Optional[str]:
        """Submit a job
        
        Args:
            script_path: Path to job script
            
        Returns:
            Job ID if successful, None otherwise
        """
        return self.slurm_manager.submit_job(script_path)
            
    def check_job_status(self, job_id: str) -> str:
        """Check job status (simplified interface)
        
        Args:
            job_id: Job ID
            
        Returns:
            Status string
        """
        return self.slurm_manager.get_job_simple_status(job_id)
    
    def create_batch_jobs(self, items: List[Tuple[Any, str]], batch_size: int, 
                         job_template: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Create batch jobs
        
        Args:
            items: List of (item_id, item_path) tuples
            batch_size: Number of items per batch
            job_template: Template with job configuration
            
        Returns:
            List of job info dictionaries
        """
        return self.slurm_manager.create_batch_jobs(items, batch_size, job_template)