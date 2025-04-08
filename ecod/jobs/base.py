#!/usr/bin/env python3
"""
Base job management interfaces for the ECOD pipeline.
"""
import abc
from typing import Dict, Any, List, Tuple, Optional, Union

class JobManager(abc.ABC):
    """Base interface for job managers"""
    
    @abc.abstractmethod
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
        pass
    
    @abc.abstractmethod
    def submit_job(self, script_path: str) -> Optional[str]:
        """Submit a job
        
        Args:
            script_path: Path to job script
            
        Returns:
            Job ID if successful, None otherwise
        """
        pass
    
    @abc.abstractmethod
    def check_job_status(self, job_id: str) -> str:
        """Check job status
        
        Args:
            job_id: Job ID
            
        Returns:
            Status string
        """
        pass
    
    @abc.abstractmethod
    def cancel_job(self, job_id: str) -> bool:
        """Cancel a job
        
        Args:
            job_id: Job ID
            
        Returns:
            True if successful
        """
        pass
    
    @abc.abstractmethod
    def get_job_output(self, job_id: str, output_dir: str, job_name: str) -> Dict[str, Optional[str]]:
        """Get job output
        
        Args:
            job_id: Job ID
            output_dir: Directory containing output files
            job_name: Job name used in output file names
            
        Returns:
            Dictionary with stdout and stderr content
        """
        pass
    
    @abc.abstractmethod
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
        pass
    
    @abc.abstractmethod
    def check_all_jobs(self, batch_id: Optional[int] = None) -> Tuple[int, int, int]:
        """Check status of all jobs
        
        Args:
            batch_id: Optional batch ID filter
            
        Returns:
            Tuple of (completed, failed, running) counts
        """
        pass