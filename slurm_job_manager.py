#!/usr/bin/env python3
"""
slurm_job_manager.py - Manage job submissions to Slurm cluster
"""

import os
import sys
import argparse
import logging
import subprocess
import time
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional

# Add parent directory to path if needed
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from ecod.core.config import ConfigManager
from ecod.core.db_manager import DBManager

def setup_logging(verbose: bool = False):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    )

class SlurmJobManager:
    """Class to manage Slurm cluster job submissions"""
    
    def __init__(self, config_path: str = None):
        """Initialize with configuration"""
        self.config_manager = ConfigManager(config_path)
        self.config = self.config_manager.config
        self.db_config = self.config_manager.get_db_config()
        self.db = DBManager(self.db_config)
        self.logger = logging.getLogger("ecod.slurm_job_manager")
        
        # Check for slurm commands
        self._check_slurm_commands()
    
    def _check_slurm_commands(self):
        """Check if slurm commands are available"""
        try:
            subprocess.run(['sbatch', '--version'], 
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            subprocess.run(['squeue', '--version'], 
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            self.logger.info("Slurm commands are available")
        except (subprocess.SubprocessError, FileNotFoundError):
            self.logger.warning("Slurm commands not found. Make sure Slurm is installed and in PATH.")
    
    def create_job_script(self, commands: List[str], job_name: str, output_dir: str, 
                        threads: int = 8, memory: str = "8G", time: str = "24:00:00") -> str:
        """Create a Slurm job script"""
        os.makedirs(output_dir, exist_ok=True)
        
        script_path = os.path.join(output_dir, f"{job_name}.sh")
        
        with open(script_path, 'w') as f:
            f.write(f"#SBATCH --time={time}\n\n")
            
            # Add module loading if configured
            modules = self.config.get('modules', [])
            if modules:
                f.write("module purge\n")
                for module in modules:
                    f.write(f"module load {module}\n")
                f.write("\n")
            
            # Add commands
            for cmd in commands:
                f.write(f"{cmd}\n")
                
        # Make executable
        os.chmod(script_path, 0o755)
        self.logger.info(f"Created job script: {script_path}")
        return script_path
        
    def submit_job(self, script_path: str) -> Optional[str]:
        """Submit a job to Slurm and return the job ID"""
        try:
            result = subprocess.run(['sbatch', script_path], 
                                  capture_output=True, text=True, check=True)
            # Expected output: "Submitted batch job 12345"
            job_id = result.stdout.strip().split()[-1]
            self.logger.info(f"Submitted job {job_id} with script {script_path}")
            return job_id
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Error submitting job: {e}")
            self.logger.error(f"STDOUT: {e.stdout}")
            self.logger.error(f"STDERR: {e.stderr}")
            return None
            
    def check_job_status(self, job_id: str) -> str:
        """Check status of a Slurm job"""
        try:
            result = subprocess.run(['sacct', '-j', job_id, '--format=State', 
                                   '--noheader', '--parsable2'],
                                  capture_output=True, text=True, check=True)
            # Split by newline and take first status (there might be multiple lines for job steps)
            status = result.stdout.strip().split('\n')[0]
            self.logger.debug(f"Job {job_id} status: {status}")
            return status
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Error checking job status: {e}")
            # If the job is not found, it's likely it hasn't started yet or has been completed
            # In this case, we check squeue to see if it's pending or running
            try:
                result = subprocess.run(['squeue', '-j', job_id, '--noheader'],
                                      capture_output=True, text=True, check=True)
                if result.stdout.strip():
                    # Job is in queue
                    return "PENDING"
                else:
                    # Job not in queue, assume completed
                    return "COMPLETED"
            except subprocess.CalledProcessError:
                # If both commands fail, return unknown
                return "UNKNOWN"
    
    def cancel_job(self, job_id: str) -> bool:
        """Cancel a running Slurm job"""
        try:
            subprocess.run(['scancel', job_id], check=True)
            self.logger.info(f"Cancelled job {job_id}")
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Error cancelling job {job_id}: {e}")
            return False
    
    def get_pending_jobs(self, batch_id: Optional[int] = None):
        """Get pending jobs from database"""
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
    
    def update_job_status(self, job_id: int, slurm_job_id: str):
        """Update job status in database"""
        # Get current status from Slurm
        status = self.check_job_status(slurm_job_id)
        
        # Map Slurm status to our status
        if status in ["COMPLETED", "COMPLETING"]:
            db_status = "completed"
        elif status in ["FAILED", "TIMEOUT", "CANCELLED", "NODE_FAIL"]:
            db_status = "failed"
        elif status in ["PENDING", "RUNNING", "CONFIGURING"]:
            db_status = "submitted"  # Still in progress
        else:
            db_status = "unknown"
            self.logger.warning(f"Unknown job status: {status}")
        
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
    
    def check_all_jobs(self, batch_id: Optional[int] = None):
        """Check status of all pending jobs"""
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

def main():
    parser = argparse.ArgumentParser(description='Slurm Job Manager for PyECOD')
    parser.add_argument('--config', type=str, default='config/config.yml', 
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int,
                      help='Specific batch ID to check')
    parser.add_argument('--check', action='store_true',
                      help='Check status of running jobs')
    parser.add_argument('--submit', type=str,
                      help='Submit a job script')
    parser.add_argument('--cancel', type=str,
                      help='Cancel a job by Slurm job ID')
    parser.add_argument('--interval', type=int, default=0,
                      help='Interval in seconds to continuously check jobs (0 = once)')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose)
    
    job_manager = SlurmJobManager(args.config)
    
    if args.submit:
        # Submit a job
        job_id = job_manager.submit_job(args.submit)
        if job_id:
            print(f"Submitted job with ID: {job_id}")
        else:
            print("Failed to submit job")
    
    elif args.cancel:
        # Cancel a job
        success = job_manager.cancel_job(args.cancel)
        if success:
            print(f"Cancelled job {args.cancel}")
        else:
            print(f"Failed to cancel job {args.cancel}")
    
    elif args.check:
        # Check job status
        if args.interval > 0:
            print(f"Checking jobs every {args.interval} seconds (press Ctrl+C to stop)")
            try:
                while True:
                    completed, failed, running = job_manager.check_all_jobs(args.batch_id)
                    print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - Completed: {completed}, Failed: {failed}, Running: {running}")
                    if completed + failed + running == 0:
                        print("No more jobs to check")
                        break
                    time.sleep(args.interval)
            except KeyboardInterrupt:
                print("\nStopped checking jobs")
        else:
            completed, failed, running = job_manager.check_all_jobs(args.batch_id)
            print(f"Completed: {completed}, Failed: {failed}, Running: {running}")
    
    else:
        parser.print_help()

if __name__ == "__main__":
    main()"#!/bin/bash\n")
            f.write(f"#SBATCH --job-name={job_name}\n")
            f.write(f"#SBATCH --output={output_dir}/{job_name}-%j.out\n")
            f.write(f"#SBATCH --error={output_dir}/{job_name}-%j.err\n")
            f.write(f"#SBATCH --mem={memory}\n")
            f.write(f"#SBATCH --cpus-per-task={threads}\n")
            f.write(