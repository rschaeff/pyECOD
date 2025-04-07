#!/usr/bin/env python3
"""
Slurm Job Manager - Enhanced job submission and management for ECOD pipeline
"""

import os
import subprocess
import logging
import re
from typing import Dict, Any, List, Tuple, Optional, Union, Callable

from .config import ConfigManager
from .db_manager import DBManager
from .exceptions import JobSubmissionError, JobExecutionError


class SlurmJobManager:
    """Enhanced job manager for Slurm cluster integration"""
    
    def __init__(self, config: Union[Dict[str, Any], str]):
        """Initialize with configuration
        
        Args:
            config: Configuration dictionary or path to config file
        """
        if isinstance(config, str):
            self.config_manager = ConfigManager(config)
            self.config = self.config_manager.config
            self.db = DBManager(self.config_manager.get_db_config())
        else:
            self.config = config
            self.db = None
        
        self.logger = logging.getLogger("ecod.slurm_job_manager")
        
        # Check for slurm commands
        self._check_slurm_commands()
    
    def _check_slurm_commands(self) -> bool:
        """Check if slurm commands are available
        
        Returns:
            True if Slurm commands are available
        """
        try:
            # Check for sbatch and squeue
            for cmd in ['sbatch', 'squeue', 'sacct']:
                result = subprocess.run([cmd, '--version'], 
                                     stdout=subprocess.PIPE, 
                                     stderr=subprocess.PIPE)
                if result.returncode != 0:
                    self.logger.warning(f"Command {cmd} returned non-zero exit code: {result.returncode}")
                    self.logger.warning(f"Stderr: {result.stderr.decode('utf-8')}")
                    return False
                    
            self.logger.info("Slurm commands are available")
            return True
        except (subprocess.SubprocessError, FileNotFoundError) as e:
            self.logger.warning(f"Slurm commands not found: {str(e)}")
            self.logger.warning("Make sure Slurm is installed and in PATH.")
            return False
    
    def create_job_script(self, commands: List[str], job_name: str, output_dir: str, 
                        threads: int = 8, memory: str = "8G", time: str = "24:00:00") -> str:
        """Create a Slurm job script with error handling
        
        Args:
            commands: List of commands to include in the script
            job_name: Name of the job
            output_dir: Directory for output files
            threads: Number of threads/CPUs
            memory: Memory allocation
            time: Time limit
            
        Returns:
            Path to the created job script
            
        Raises:
            JobSubmissionError: If script creation fails
        """
        try:
            # Ensure output directory exists
            os.makedirs(output_dir, exist_ok=True)
            
            script_path = os.path.join(output_dir, f"{job_name}.sh")
            
            self.logger.debug(f"Creating job script: {script_path}")
            
            with open(script_path, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write(f"#SBATCH --job-name={job_name}\n")
                f.write(f"#SBATCH --output={output_dir}/{job_name}-%j.out\n")
                f.write(f"#SBATCH --error={output_dir}/{job_name}-%j.err\n")
                f.write(f"#SBATCH --mem={memory}\n")
                f.write(f"#SBATCH --cpus-per-task={threads}\n")
                f.write(f"#SBATCH --time={time}\n\n")
                
                # Add error handling to the job script
                f.write("# Error handling\n")
                f.write("set -e\n")
                f.write("trap 'echo \"Error occurred at line $LINENO\"' ERR\n\n")
                
                # Record job start
                f.write("echo \"Job started at $(date)\"\n")
                f.write(f"echo \"Job name: {job_name}\"\n\n")
                
                # Add module loading if configured
                modules = self.config.get('modules', [])
                if modules:
                    f.write("# Load required modules\n")
                    f.write("module purge\n")
                    for module in modules:
                        f.write(f"module load {module}\n")
                    f.write("\n")
                
                # Add commands
                f.write("# Job commands\n")
                for i, cmd in enumerate(commands):
                    f.write(f"echo \"Executing command {i+1}/{len(commands)}: {cmd}\"\n")
                    f.write(f"{cmd}\n")
                    f.write(f"echo \"Command {i+1} completed with exit code $?\"\n\n")
                
                # Record job end
                f.write("echo \"Job completed at $(date)\"\n")
                
            # Make executable
            os.chmod(script_path, 0o755)
            self.logger.info(f"Created job script: {script_path}")
            return script_path
            
        except (OSError, IOError) as e:
            error_msg = f"Error creating job script: {str(e)}"
            self.logger.error(error_msg)
            raise JobSubmissionError(error_msg) from e
        
    def submit_job(self, script_path: str) -> Optional[str]:
        """Submit a job to Slurm with enhanced error handling
        
        Args:
            script_path: Path to the job script
            
        Returns:
            Job ID if successful, None otherwise
            
        Raises:
            JobSubmissionError: If submission fails
        """
        if not os.path.exists(script_path):
            error_msg = f"Job script not found: {script_path}"
            self.logger.error(error_msg)
            raise JobSubmissionError(error_msg)
            
        try:
            self.logger.debug(f"Submitting job with script: {script_path}")
            result = subprocess.run(['sbatch', script_path], 
                                  capture_output=True, text=True, check=True)
            
            # Parse job ID from output (expected: "Submitted batch job 12345")
            output = result.stdout.strip()
            self.logger.debug(f"sbatch output: {output}")
            
            # More robust job ID extraction
            job_id_match = re.search(r'Submitted batch job (\d+)', output)
            
            if job_id_match:
                job_id = job_id_match.group(1)
                self.logger.info(f"Submitted job {job_id} with script {script_path}")
                return job_id
            else:
                error_msg = f"Failed to parse job ID from sbatch output: {output}"
                self.logger.error(error_msg)
                raise JobSubmissionError(error_msg)
                
        except subprocess.CalledProcessError as e:
            error_msg = f"Error submitting job: {e}\nSTDOUT: {e.stdout}\nSTDERR: {e.stderr}"
            self.logger.error(error_msg)
            raise JobSubmissionError(error_msg) from e
        except Exception as e:
            error_msg = f"Unexpected error submitting job: {str(e)}"
            self.logger.error(error_msg)
            raise JobSubmissionError(error_msg) from e

    def check_job_status(self, job_id: str) -> Dict[str, Any]:
        """Check status of a Slurm job with enhanced error handling and information
        
        Args:
            job_id: Slurm job ID
            
        Returns:
            Dictionary with job status information
        """
        status_info = {
            'state': 'UNKNOWN',
            'exit_code': None,
            'start_time': None,
            'end_time': None,
            'elapsed': None,
            'node': None,
            'error': None
        }
        
        try:
            # First try sacct for detailed information
            sacct_cmd = [
                'sacct', 
                '-j', job_id, 
                '--format=JobID,State,ExitCode,Start,End,Elapsed,NodeList', 
                '--parsable2', 
                '--noheader'
            ]
            
            self.logger.debug(f"Checking job status with sacct: {job_id}")
            result = subprocess.run(sacct_cmd, capture_output=True, text=True)
            
            if result.returncode == 0 and result.stdout.strip():
                lines = result.stdout.strip().split('\n')
                fields = lines[0].split('|')
                
                if len(fields) >= 7:
                    status_info['state'] = fields[1]
                    status_info['exit_code'] = fields[2]
                    status_info['start_time'] = fields[3]
                    status_info['end_time'] = fields[4]
                    status_info['elapsed'] = fields[5]
                    status_info['node'] = fields[6]
                    
                    self.logger.debug(f"Job {job_id} status: {status_info['state']}")
                    return status_info
            
            # If sacct doesn't return useful information, try squeue
            squeue_cmd = ['squeue', '-j', job_id, '--noheader', '-o', '%T|%N']
            
            self.logger.debug(f"Checking job status with squeue: {job_id}")
            result = subprocess.run(squeue_cmd, capture_output=True, text=True)
            
            if result.returncode == 0 and result.stdout.strip():
                fields = result.stdout.strip().split('|')
                status_info['state'] = fields[0]
                if len(fields) > 1:
                    status_info['node'] = fields[1]
                    
                self.logger.debug(f"Job {job_id} status from squeue: {status_info['state']}")
                return status_info
            
            # If not found in squeue, job is probably done or doesn't exist
            # Check scontrol for job information
            scontrol_cmd = ['scontrol', 'show', 'job', job_id]
            result = subprocess.run(scontrol_cmd, capture_output=True, text=True)
            
            if "Invalid job id specified" in result.stderr:
                # Job doesn't exist or has been purged from Slurm's memory
                # We assume it completed if we can't find any record
                status_info['state'] = "COMPLETED"
                status_info['error'] = "Job not found in Slurm database"
                self.logger.warning(f"Job {job_id} not found in Slurm database, assuming COMPLETED")
            else:
                # Job exists but has an unclear state
                status_info['state'] = "UNKNOWN"
                status_info['error'] = "Unable to determine job state from Slurm commands"
                self.logger.warning(f"Unable to determine state for job {job_id}")
            
            return status_info
                
        except subprocess.SubprocessError as e:
            error_msg = f"Error checking job status: {str(e)}"
            self.logger.error(error_msg)
            status_info['error'] = error_msg
            return status_info

    def get_job_simple_status(self, job_id: str) -> str:
        """Get simplified job status string
        
        Args:
            job_id: Slurm job ID
            
        Returns:
            Status string (COMPLETED, RUNNING, FAILED, etc.)
        """
        status_info = self.check_job_status(job_id)
        return status_info['state']

    def cancel_job(self, job_id: str) -> bool:
        """Cancel a running Slurm job with enhanced error handling
        
        Args:
            job_id: Slurm job ID
            
        Returns:
            True if cancellation was successful
            
        Raises:
            JobExecutionError: If cancellation fails
        """
        try:
            self.logger.debug(f"Cancelling job: {job_id}")
            result = subprocess.run(['scancel', job_id], 
                                  capture_output=True, text=True, check=True)
            self.logger.info(f"Cancelled job {job_id}")
            return True
        except subprocess.CalledProcessError as e:
            error_msg = f"Error cancelling job {job_id}: {str(e)}\nSTDOUT: {e.stdout}\nSTDERR: {e.stderr}"
            self.logger.error(error_msg)
            # Check if job doesn't exist or is already completed
            if "Invalid job id specified" in e.stderr:
                self.logger.info(f"Job {job_id} already completed or doesn't exist")
                return True
            else:
                raise JobExecutionError(error_msg) from e
        except Exception as e:
            error_msg = f"Unexpected error cancelling job {job_id}: {str(e)}"
            self.logger.error(error_msg)
            raise JobExecutionError(error_msg) from e

    def get_job_output(self, job_id: str, output_dir: str, job_name: str) -> Dict[str, Optional[str]]:
        """Get output and error logs from a job
        
        Args:
            job_id: Slurm job ID
            output_dir: Directory containing output files
            job_name: Job name used in output file names
            
        Returns:
            Dictionary with 'stdout' and 'stderr' content
        """
        result = {
            'stdout': None,
            'stderr': None
        }
        
        stdout_path = os.path.join(output_dir, f"{job_name}-{job_id}.out")
        stderr_path = os.path.join(output_dir, f"{job_name}-{job_id}.err")
        
        try:
            if os.path.exists(stdout_path):
                with open(stdout_path, 'r') as f:
                    result['stdout'] = f.read()
                self.logger.debug(f"Read stdout file: {stdout_path}")
            else:
                self.logger.warning(f"Stdout file not found: {stdout_path}")
                
            if os.path.exists(stderr_path):
                with open(stderr_path, 'r') as f:
                    result['stderr'] = f.read()
                self.logger.debug(f"Read stderr file: {stderr_path}")
            else:
                self.logger.warning(f"Stderr file not found: {stderr_path}")
                
        except (OSError, IOError) as e:
            self.logger.error(f"Error reading job output files: {str(e)}")
            # Don't raise exception, just return what we have
            
        return result
    
    def get_pending_jobs(self, batch_id: Optional[int] = None) -> List[Dict[str, Any]]:
        """Get pending jobs from database
        
        Args:
            batch_id: Optional batch ID filter
            
        Returns:
            List of pending job dictionaries
        """
        if not self.db:
            self.logger.error("Database connection not available")
            return []
            
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
            return self.db.execute_dict_query(query, (batch_id,))
        else:
            return self.db.execute_dict_query(query)
    
    def update_job_status(self, job_id: int, slurm_job_id: str) -> str:
        """Update job status in database
        
        Args:
            job_id: Database job ID
            slurm_job_id: Slurm job ID
            
        Returns:
            Updated status string
        """
        if not self.db:
            self.logger.error("Database connection not available")
            return "unknown"
            
        # Get current status from Slurm
        status_info = self.check_job_status(slurm_job_id)
        slurm_status = status_info['state']
        
        # Map Slurm status to our status
        if slurm_status in ["COMPLETED", "COMPLETING"]:
            db_status = "completed"
        elif slurm_status in ["FAILED", "TIMEOUT", "CANCELLED", "NODE_FAIL"]:
            db_status = "failed"
        elif slurm_status in ["PENDING", "RUNNING", "CONFIGURING"]:
            db_status = "submitted"  # Still in progress
        else:
            db_status = "unknown"
            self.logger.warning(f"Unknown job status: {slurm_status}")
        
        # Only update if status changed
        if db_status != "submitted":
            # Update job status
            self.db.update(
                "ecod_schema.job",
                {
                    "status": db_status,
                    "completion_time": "CURRENT_TIMESTAMP" if db_status in ["completed", "failed"] else None
                },
                "id = %s",
                (job_id,)
            )
            self.logger.info(f"Updated job {job_id} (Slurm job {slurm_job_id}) status to {db_status}")
        
        return db_status
    
    def check_all_jobs(self, batch_id: Optional[int] = None) -> Tuple[int, int, int]:
        """Check status of all pending jobs
        
        Args:
            batch_id: Optional batch ID filter
            
        Returns:
            Tuple of (completed, failed, still_running) counts
        """
        if not self.db:
            self.logger.error("Database connection not available")
            return (0, 0, 0)
            
        jobs = self.get_pending_jobs(batch_id)
        
        completed = 0
        failed = 0
        still_running = 0
        
        for job in jobs:
            job_id = job['id']
            slurm_job_id = job['slurm_job_id']
            
            # Update status
            status = self.update_job_status(job_id, slurm_job_id)
            
            if status == "completed":
                completed += 1
            elif status == "failed":
                failed += 1
            else:
                still_running += 1
        
        self.logger.info(f"Checked {len(jobs)} jobs: {completed} completed, {failed} failed, {still_running} still running")
        return completed, failed, still_running
    
    def create_batch_jobs(self, items: List[Tuple[Any, str]], batch_size: int, 
                        job_template: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Create batch jobs for a list of items
        
        Args:
            items: List of (item_id, item_path) tuples
            batch_size: Number of items per batch
            job_template: Template with keys: 'name', 'output_dir', 'command_template'
                         command_template should be a function that takes item_id, item_path
        
        Returns:
            List of job info dictionaries with 'script_path', 'items'
        """
        jobs = []
        
        # Process in batches
        for i in range(0, len(items), batch_size):
            batch_items = items[i:i+batch_size]
            batch_num = i // batch_size
            job_name = f"{job_template['name']}_{batch_num}"
            output_dir = job_template['output_dir']
            
            # Create commands for this batch
            commands = []
            for item_id, item_path in batch_items:
                cmd = job_template['command_template'](item_id, item_path)
                commands.append(cmd)
            
            # Create job script
            script_path = self.create_job_script(
                commands, 
                job_name, 
                output_dir,
                threads=job_template.get('threads', 8),
                memory=job_template.get('memory', '8G'),
                time=job_template.get('time', '24:00:00')
            )
            
            # Add job info
            jobs.append({
                'script_path': script_path,
                'items': batch_items,
                'name': job_name,
                'batch_num': batch_num
            })
            
        return jobs