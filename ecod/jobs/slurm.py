#!/usr/bin/env python3
"""
SLURM job manager for the ECOD pipeline.
"""
import os
import subprocess
import logging
import re
import time
from typing import Dict, Any, List, Tuple, Optional, Union

from ecod.exceptions import JobSubmissionError, JobExecutionError
from ecod.utils.file import ensure_dir, write_text_file
from .base import JobManager

logger = logging.getLogger("ecod.jobs.slurm")

class SlurmJobManager(JobManager):
    """Job manager for SLURM cluster execution"""
    
    def __init__(self, config: Dict[str, Any]):
        """Initialize with configuration"""
        self.config = config
        
        # Check for SLURM commands
        self._check_slurm_commands()
    
    def _check_slurm_commands(self) -> bool:
        """Check if SLURM commands are available
        
        Returns:
            True if commands are available
        """
        try:
            # Check for basic SLURM commands
            for cmd in ['sbatch', 'squeue', 'sacct']:
                result = subprocess.run(
                    [cmd, '--version'],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    check=False
                )
                
                if result.returncode != 0:
                    logger.warning(f"Command {cmd} returned non-zero exit code: {result.returncode}")
                    logger.warning(f"Stderr: {result.stderr.decode('utf-8')}")
                    return False
            
            logger.info("SLURM commands are available")
            return True
        except (subprocess.SubprocessError, FileNotFoundError) as e:
            logger.warning(f"SLURM commands not found: {str(e)}")
            logger.warning("Make sure SLURM is installed and in PATH")
            return False
    
    def create_job_script(self, commands: List[str], job_name: str, 
                         output_dir: str, **options
    ) -> str:
        """Create a SLURM job script
        
        Args:
            commands: List of commands to run
            job_name: Name of the job
            output_dir: Directory for output files
            **options: Additional SLURM options
                - threads: Number of CPUs (default: 1)
                - memory: Memory allocation (default: '4G')
                - time: Time limit (default: '01:00:00')
            
        Returns:
            Path to the created job script
        """
        # Get SLURM options with defaults
        threads = options.get('threads', 1)
        memory = options.get('memory', '4G')
        time_limit = options.get('time', '01:00:00')
        
        # Ensure output directory exists
        ensure_dir(output_dir)
        script_path = os.path.join(output_dir, f"{job_name}.sh")
        
        # Create script content
        script_content = "#!/bin/bash\n\n"
        
        # Add SLURM directives
        script_content += f"#SBATCH --job-name={job_name}\n"
        script_content += f"#SBATCH --output={output_dir}/{job_name}-%j.out\n"
        script_content += f"#SBATCH --error={output_dir}/{job_name}-%j.err\n"
        script_content += f"#SBATCH --cpus-per-task={threads}\n"
        script_content += f"#SBATCH --mem={memory}\n"
        script_content += f"#SBATCH --time={time_limit}\n\n"
        
        # Add error handling
        script_content += "# Error handling\n"
        script_content += "set -e\n"
        script_content += "trap 'echo \"Error occurred at line $LINENO\"' ERR\n\n"
        
        # Record job start
        script_content += "echo \"Job started at $(date)\"\n"
        script_content += f"echo \"Job name: {job_name}\"\n\n"
        
        # Add module loading if configured
        modules = self.config.get('modules', [])
        if modules:
            script_content += "# Load required modules\n"
            script_content += "module purge\n"
            for module in modules:
                script_content += f"module load {module}\n"
            script_content += "\n"
        
        # Add commands
        script_content += "# Commands\n"
        for i, cmd in enumerate(commands):
            script_content += f"echo \"Running command {i+1}/{len(commands)}: {cmd}\"\n"
            script_content += f"{cmd}\n"
            script_content += f"echo \"Command {i+1} completed with exit code $?\"\n\n"
        
        # Record job end
        script_content += "echo \"Job completed at $(date)\"\n"
        
        # Write script file
        write_text_file(script_path, script_content)
        os.chmod(script_path, 0o755)
        
        logger.info(f"Created SLURM job script: {script_path}")
        return script_path
    
    def submit_job(self, script_path: str) -> Optional[str]:
        """Submit a job to SLURM
        
        Args:
            script_path: Path to job script
            
        Returns:
            SLURM job ID if successful, None otherwise
        """
        if not os.path.exists(script_path):
            raise JobSubmissionError(f"Job script not found: {script_path}")
        
        try:
            logger.debug(f"Submitting job with script: {script_path}")
            result = subprocess.run(
                ['sbatch', script_path],
                capture_output=True,
                text=True,
                check=True
            )
            
            # Parse job ID from output (expected: "Submitted batch job 12345")
            output = result.stdout.strip()
            logger.debug(f"sbatch output: {output}")
            
            job_id_match = re.search(r'Submitted batch job (\d+)', output)
            
            if job_id_match:
                job_id = job_id_match.group(1)
                logger.info(f"Submitted job {job_id} with script {script_path}")
                return job_id
            else:
                error_msg = f"Failed to parse job ID from sbatch output: {output}"
                logger.error(error_msg)
                raise JobSubmissionError(error_msg)
                
        except subprocess.CalledProcessError as e:
            error_msg = f"Error submitting job: {e}\nSTDOUT: {e.stdout}\nSTDERR: {e.stderr}"
            logger.error(error_msg)
            raise JobSubmissionError(error_msg) from e
        except Exception as e:
            error_msg = f"Unexpected error submitting job: {str(e)}"
            logger.error(error_msg)
            raise JobSubmissionError(error_msg) from e
    
    def check_job_status(self, job_id: str) -> str:
        """Check status of a SLURM job
        
        Args:
            job_id: SLURM job ID
            
        Returns:
            Status string (PENDING, RUNNING, COMPLETED, FAILED, etc.)
        """
        status_info = self._get_job_status_info(job_id)
        return status_info['state']
    
    def _get_job_status_info(self, job_id: str) -> Dict[str, Any]:
        """Get detailed status information for a SLURM job
        
        Args:
            job_id: SLURM job ID
            
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
            
            logger.debug(f"Checking job status with sacct: {job_id}")
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
                    
                    logger.debug(f"Job {job_id} status: {status_info['state']}")
                    return status_info
            
            # If sacct doesn't return useful information, try squeue
            squeue_cmd = ['squeue', '-j', job_id, '--noheader', '-o', '%T|%N']
            
            logger.debug(f"Checking job status with squeue: {job_id}")
            result = subprocess.run(squeue_cmd, capture_output=True, text=True)
            
            if result.returncode == 0 and result.stdout.strip():
                fields = result.stdout.strip().split('|')
                status_info['state'] = fields[0]
                if len(fields) > 1:
                    status_info['node'] = fields[1]
                    
                logger.debug(f"Job {job_id} status from squeue: {status_info['state']}")
                return status_info
            
            # If not found in squeue, job is probably done or doesn't exist
            # Check scontrol for job information
            scontrol_cmd = ['scontrol', 'show', 'job', job_id]
            result = subprocess.run(scontrol_cmd, capture_output=True, text=True)
            
            if "Invalid job id specified" in result.stderr:
                # Job doesn't exist or has been purged from SLURM's memory
                # We assume it completed if we can't find any record
                status_info['state'] = "COMPLETED"
                status_info['error'] = "Job not found in SLURM database"
                logger.warning(f"Job {job_id} not found in SLURM database, assuming COMPLETED")
            else:
                # Job exists but has an unclear state
                status_info['state'] = "UNKNOWN"
                status_info['error'] = "Unable to determine job state"
                logger.warning(f"Unable to determine state for job {job_id}")
            
            return status_info
                
        except subprocess.SubprocessError as e:
            error_msg = f"Error checking job status: {str(e)}"
            logger.error(error_msg)
            status_info['error'] = error_msg
            return status_info
    
    def cancel_job(self, job_id: str) -> bool:
        """Cancel a running SLURM job
        
        Args:
            job_id: SLURM job ID
            
        Returns:
            True if cancellation was successful
        """
        try:
            logger.debug(f"Cancelling job: {job_id}")
            result = subprocess.run(
                ['scancel', job_id], 
                capture_output=True, 
                text=True, 
                check=True
            )
            logger.info(f"Cancelled job {job_id}")
            return True
        except subprocess.CalledProcessError as e:
            error_msg = f"Error cancelling job {job_id}: {str(e)}\nSTDOUT: {e.stdout}\nSTDERR: {e.stderr}"
            logger.error(error_msg)
            
            # Check if job doesn't exist or is already completed
            if "Invalid job id specified" in e.stderr:
                logger.info(f"Job {job_id} already completed or doesn't exist")
                return True
            
            raise JobExecutionError(error_msg) from e
        except Exception as e:
            error_msg = f"Unexpected error cancelling job {job_id}: {str(e)}"
            logger.error(error_msg)
            raise JobExecutionError(error_msg) from e
    
    def get_job_output(self, job_id: str, output_dir: str, job_name: str) -> Dict[str, Optional[str]]:
        """Get output and error logs from a SLURM job
        
        Args:
            job_id: SLURM job ID
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
                logger.debug(f"Read stdout file: {stdout_path}")
            else:
                logger.warning(f"Stdout file not found: {stdout_path}")
                
            if os.path.exists(stderr_path):
                with open(stderr_path, 'r') as f:
                    result['stderr'] = f.read()
                logger.debug(f"Read stderr file: {stderr_path}")
            else:
                logger.warning(f"Stderr file not found: {stderr_path}")
                
        except (OSError, IOError) as e:
            logger.error(f"Error reading job output files: {str(e)}")
            
        return result
    
    def create_batch_jobs(self, items: List[Tuple[Any, str]], batch_size: int, 
                         job_template: Dict[str, Any]
    ) -> List[Dict[str, Any]]:
        """Create batch jobs for a list of items
        
        Args:
            items: List of (item_id, item_path) tuples
            batch_size: Number of items per batch
            job_template: Template with job configuration
            
        Returns:
            List of job info dictionaries with 'script_path', 'items'
        """
        jobs = []
        
        # Process in batches
        for i in range(0, len(items), batch_size):
            batch_items = items[i:i+batch_size]
            batch_num = i // batch_size
            job_name = f"{job_template.get('name', 'batch')}_{batch_num}"
            output_dir = job_template.get('output_dir', './output')
            
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
                threads=job_template.get('threads', 1),
                memory=job_template.get('memory', '4G'),
                time=job_template.get('time', '01:00:00')
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
        """Check status of all jobs based on SLURM information
        
        Args:
            batch_id: Optional batch ID filter (not used in this direct SLURM implementation)
                but kept for API compatibility with DatabaseJobManager
            
        Returns:
            Tuple of (completed, failed, running) counts
        """
        completed = 0
        failed = 0
        running = 0
        
        try:
            # Get all running and pending jobs using squeue
            logger.debug("Checking running jobs with squeue")
            squeue_cmd = ['squeue', '--noheader', '-o', '%i|%T']
            
            result = subprocess.run(squeue_cmd, capture_output=True, text=True)
            
            # Track job IDs we've seen
            running_job_ids = set()
            
            if result.returncode == 0 and result.stdout.strip():
                for line in result.stdout.strip().split('\n'):
                    if not line:
                        continue
                        
                    parts = line.strip().split('|')
                    if len(parts) < 2:
                        continue
                        
                    job_id = parts[0]
                    state = parts[1]
                    
                    running_job_ids.add(job_id)
                    
                    if state in ['RUNNING', 'PENDING', 'CONFIGURING', 'COMPLETING']:
                        running += 1
            
            # Get recently completed jobs using sacct
            # Look at jobs completed in the last 24 hours
            logger.debug("Checking completed jobs with sacct")
            sacct_cmd = [
                'sacct', 
                '--starttime', 'now-24hours',
                '--format=JobID,State', 
                '--parsable2', 
                '--noheader'
            ]
            
            result = subprocess.run(sacct_cmd, capture_output=True, text=True)
            
            if result.returncode == 0 and result.stdout.strip():
                for line in result.stdout.strip().split('\n'):
                    if not line:
                        continue
                        
                    parts = line.strip().split('|')
                    if len(parts) < 2:
                        continue
                    
                    # Skip array job summaries (they have JobID.batch)
                    job_id_part = parts[0]
                    if '.' in job_id_part:
                        continue
                        
                    job_id = job_id_part
                    state = parts[1]
                    
                    # Skip jobs we already counted in squeue
                    if job_id in running_job_ids:
                        continue
                    
                    if state in ['COMPLETED', 'COMPLETING']:
                        completed += 1
                    elif state in ['FAILED', 'TIMEOUT', 'CANCELLED', 'NODE_FAIL', 'OUT_OF_MEMORY']:
                        failed += 1
            
            logger.info(f"SLURM job status: {completed} completed, {failed} failed, {running} running")
            return completed, failed, running
                    
        except subprocess.SubprocessError as e:
            logger.error(f"Error checking SLURM job status: {str(e)}")
            return (0, 0, 0)

    def wait_for_jobs(job_manager, job_ids, check_interval=60, timeout=None, progress_callback=None):
        """Wait for SLURM jobs to complete with timeout and progress updates"""
        start_time = time.time()
        all_completed = False
        
        while not all_completed:
            status_counts = {"COMPLETED": 0, "FAILED": 0, "RUNNING": 0, "PENDING": 0}
            
            for job_id in job_ids:
                status = job_manager.check_job_status(job_id)
                if status in status_counts:
                    status_counts[status] += 1
            
            # Check if all jobs completed
            all_completed = (status_counts["COMPLETED"] + status_counts["FAILED"]) == len(job_ids)
            
            # Call progress callback if provided
            if progress_callback:
                progress_callback(status_counts)
            
            # Check timeout
            if timeout and (time.time() - start_time) > timeout:
                logger.warning(f"Timeout reached after {timeout} seconds while waiting for jobs")
                return False
            
            if not all_completed:
                time.sleep(check_interval)
        
        return all_completed