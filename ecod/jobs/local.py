#!/usr/bin/env python3
"""
Local job manager for the ECOD pipeline.
Executes jobs on the local machine.
"""
import os
import subprocess
import logging
import time
import shlex
from typing import Dict, Any, List, Tuple, Optional, Union
import threading
import uuid

from ecod.exceptions import JobSubmissionError, JobExecutionError
from ecod.utils.file import ensure_dir, write_text_file
from .base import JobManager

logger = logging.getLogger("ecod.jobs.local")

class LocalJobManager(JobManager):
    """Job manager for local execution"""
    
    def __init__(self, config: Dict[str, Any]):
        """Initialize with configuration"""
        self.config = config
        self.running_jobs: Dict[str, Dict[str, Any]] = {}
    
    def create_job_script(self, commands: List[str], job_name: str, 
                         output_dir: str, **options) -> str:
        """Create a job script for local execution
        
        Args:
            commands: List of commands to run
            job_name: Name of the job
            output_dir: Directory for output files
            **options: Additional job options (ignored for local execution)
            
        Returns:
            Path to the created job script
        """
        ensure_dir(output_dir)
        script_path = os.path.join(output_dir, f"{job_name}.sh")
        
        script_content = "#!/bin/bash\n\n"
        script_content += f"# Job: {job_name}\n"
        script_content += f"# Created: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n"
        
        # Add error handling
        script_content += "# Error handling\n"
        script_content += "set -e\n"
        script_content += "trap 'echo \"Error occurred at line $LINENO\"' ERR\n\n"
        
        # Add commands
        script_content += "# Commands\n"
        for i, cmd in enumerate(commands):
            script_content += f"echo \"Running command {i+1}/{len(commands)}: {cmd}\"\n"
            script_content += f"{cmd}\n"
            script_content += f"echo \"Command {i+1} completed with exit code $?\"\n\n"
        
        # Write script file
        write_text_file(script_path, script_content)
        os.chmod(script_path, 0o755)
        
        logger.info(f"Created job script: {script_path}")
        return script_path
    
    def submit_job(self, script_path: str) -> Optional[str]:
        """Submit a job for local execution
        
        Args:
            script_path: Path to job script
            
        Returns:
            Job ID if successful, None otherwise
        """
        if not os.path.exists(script_path):
            raise JobSubmissionError(f"Job script not found: {script_path}")
        
        # Generate a unique job ID
        job_id = str(uuid.uuid4())
        
        # Prepare output files
        output_dir = os.path.dirname(script_path)
        job_name = os.path.basename(script_path).replace('.sh', '')
        stdout_path = os.path.join(output_dir, f"{job_name}.out")
        stderr_path = os.path.join(output_dir, f"{job_name}.err")
        
        # Start job in a separate thread
        thread = threading.Thread(
            target=self._run_job_thread,
            args=(job_id, script_path, stdout_path, stderr_path)
        )
        thread.daemon = True
        
        # Store job information
        self.running_jobs[job_id] = {
            'script_path': script_path,
            'stdout_path': stdout_path,
            'stderr_path': stderr_path,
            'start_time': time.time(),
            'thread': thread,
            'status': 'PENDING',
            'process': None
        }
        
        # Start the thread
        thread.start()
        logger.info(f"Submitted job {job_id} with script {script_path}")
        
        return job_id
    
    def _run_job_thread(self, job_id: str, script_path: str, 
                       stdout_path: str, stderr_path: str) -> None:
        """Run a job in a separate thread
        
        Args:
            job_id: Job ID
            script_path: Path to job script
            stdout_path: Path to stdout file
            stderr_path: Path to stderr file
        """
        try:
            # Open output files
            with open(stdout_path, 'w') as stdout_file, open(stderr_path, 'w') as stderr_file:
                # Update status
                self.running_jobs[job_id]['status'] = 'RUNNING'
                
                # Start process
                process = subprocess.Popen(
                    [script_path],
                    stdout=stdout_file,
                    stderr=stderr_file,
                    shell=False
                )
                
                # Store process
                self.running_jobs[job_id]['process'] = process
                
                # Wait for completion
                returncode = process.wait()
                
                # Update status
                if returncode == 0:
                    self.running_jobs[job_id]['status'] = 'COMPLETED'
                else:
                    self.running_jobs[job_id]['status'] = 'FAILED'
                
                self.running_jobs[job_id]['return_code'] = returncode
                self.running_jobs[job_id]['end_time'] = time.time()
                
                logger.info(f"Job {job_id} finished with status {self.running_jobs[job_id]['status']}")
                
        except Exception as e:
            # Handle any exceptions
            logger.error(f"Error running job {job_id}: {str(e)}")
            self.running_jobs[job_id]['status'] = 'FAILED'
            self.running_jobs[job_id]['error'] = str(e)
            self.running_jobs[job_id]['end_time'] = time.time()
    
    def check_job_status(self, job_id: str) -> str:
        """Check job status
        
        Args:
            job_id: Job ID
            
        Returns:
            Status string
        """
        if job_id not in self.running_jobs:
            return 'UNKNOWN'
        
        return self.running_jobs[job_id]['status']
    
    def cancel_job(self, job_id: str) -> bool:
        """Cancel a job
        
        Args:
            job_id: Job ID
            
        Returns:
            True if successful
        """
        if job_id not in self.running_jobs:
            logger.warning(f"Job {job_id} not found")
            return False
        
        job = self.running_jobs[job_id]
        
        if job['status'] not in ('PENDING', 'RUNNING'):
            logger.info(f"Job {job_id} is already in state {job['status']}, cannot cancel")
            return False
        
        # Try to terminate the process
        process = job.get('process')
        if process and process.poll() is None:
            process.terminate()
            
            # Give it a moment to terminate
            for _ in range(5):
                if process.poll() is not None:
                    break
                time.sleep(0.1)
                
            # Force kill if still running
            if process.poll() is None:
                process.kill()
        
        # Update status
        job['status'] = 'CANCELLED'
        job['end_time'] = time.time()
        
        logger.info(f"Cancelled job {job_id}")
        return True
    
    def get_job_output(self, job_id: str, output_dir: Optional[str] = None, 
                      job_name: Optional[str] = None) -> Dict[str, Optional[str]]:
        """Get job output
        
        Args:
            job_id: Job ID
            output_dir: Directory containing output files (optional if job is known)
            job_name: Job name used in output file names (optional if job is known)
            
        Returns:
            Dictionary with stdout and stderr content
        """
        result = {
            'stdout': None,
            'stderr': None
        }
        
        # Get output paths
        if job_id in self.running_jobs:
            stdout_path = self.running_jobs[job_id].get('stdout_path')
            stderr_path = self.running_jobs[job_id].get('stderr_path')
        elif output_dir and job_name:
            stdout_path = os.path.join(output_dir, f"{job_name}.out")
            stderr_path = os.path.join(output_dir, f"{job_name}.err")
        else:
            logger.warning(f"Job {job_id} not found and no output paths provided")
            return result
        
        # Read output files
        try:
            if stdout_path and os.path.exists(stdout_path):
                with open(stdout_path, 'r') as f:
                    result['stdout'] = f.read()
            
            if stderr_path and os.path.exists(stderr_path):
                with open(stderr_path, 'r') as f:
                    result['stderr'] = f.read()
        except Exception as e:
            logger.warning(f"Error reading job output files: {str(e)}")
        
        return result
    
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
        jobs = []
        
        # Get template values
        output_dir = job_template.get('output_dir', './output')
        name_template = job_template.get('name', 'batch_job')
        command_template = job_template.get('command_template')
        
        if not command_template:
            raise ValueError("command_template is required in job_template")
        
        # Create batches
        for i in range(0, len(items), batch_size):
            batch_items = items[i:i+batch_size]
            batch_num = i // batch_size
            job_name = f"{name_template}_{batch_num}"
            
            # Generate commands for this batch
            commands = []
            for item_id, item_path in batch_items:
                cmd = command_template(item_id, item_path)
                commands.append(cmd)
            
            # Create job script
            script_path = self.create_job_script(
                commands, 
                job_name, 
                output_dir,
                threads=job_template.get('threads', 1)
            )
            
            # Add job info
            jobs.append({
                'script_path': script_path,
                'items': batch_items,
                'name': job_name,
                'batch_num': batch_num
            })
        
        return jobs
    
    def check_all_jobs(self, batch_id: Optional[int] = None) -> Tuple[int, int, int]:
        """Check status of all jobs
        
        Args:
            batch_id: Optional batch ID filter (ignored for local jobs)
            
        Returns:
            Tuple of (completed, failed, running) counts
        """
        completed = 0
        failed = 0
        running = 0
        
        for job_id, job in list(self.running_jobs.items()):
            status = job['status']
            
            if status == 'COMPLETED':
                completed += 1
            elif status in ('FAILED', 'CANCELLED'):
                failed += 1
            elif status in ('PENDING', 'RUNNING'):
                running += 1
                
                # Check if thread is still alive
                thread = job.get('thread')
                if thread and not thread.is_alive():
                    # Thread ended but status wasn't updated
                    process = job.get('process')
                    if process and process.poll() is not None:
                        # Process ended
                        if process.returncode == 0:
                            job['status'] = 'COMPLETED'
                            completed += 1
                            running -= 1
                        else:
                            job['status'] = 'FAILED'
                            failed += 1
                            running -= 1
        
        return completed, failed, running