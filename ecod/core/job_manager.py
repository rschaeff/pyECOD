# ecod/core/job_manager.py
import os
import subprocess
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple

class JobManager:
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.logger = logging.getLogger("ecod.jobs")
        
    def create_job_script(self, commands: List[str], job_name: str, output_dir: str, 
                         threads: int = 8, memory: str = "8G", time: str = "24:00:00") -> str:
        """Create a Slurm job script"""
        script_path = os.path.join(output_dir, f"{job_name}.sh")
        
        with open(script_path, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write(f"#SBATCH --job-name={job_name}\n")
            f.write(f"#SBATCH --output={output_dir}/{job_name}-%j.out\n")
            f.write(f"#SBATCH --error={output_dir}/{job_name}-%j.err\n")
            f.write(f"#SBATCH --mem={memory}\n")
            f.write(f"#SBATCH --cpus-per-task={threads}\n")
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
        return script_path
        
    def submit_job(self, script_path: str) -> Optional[str]:
        """Submit a job to Slurm and return the job ID"""
        try:
            result = subprocess.run(['sbatch', script_path], 
                                  capture_output=True, text=True, check=True)
            job_id = result.stdout.strip().split()[-1]
            self.logger.info(f"Submitted job {job_id} with script {script_path}")
            return job_id
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Error submitting job: {e}")
            return None
            
    def check_job_status(self, job_id: str) -> str:
        """Check status of a Slurm job"""
        try:
            result = subprocess.run(['sacct', '-j', job_id, '--format=State', 
                                   '--noheader', '--parsable2'],
                                  capture_output=True, text=True, check=True)
            status = result.stdout.strip()
            return status
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Error checking job status: {e}")
            return "UNKNOWN"
            
    def create_batch_jobs(self, items: List[Tuple[str, str]], batch_size: int, 
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