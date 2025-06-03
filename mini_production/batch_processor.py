#!/usr/bin/env python3
"""
Production Batch Processor for Mini PyECOD
Emergency wrapper to scale mini to production workloads

Usage:
    python batch_processor.py --batch-size 100 --max-jobs 50
    python batch_processor.py --protein-list proteins.txt
    python batch_processor.py --status  # Monitor progress
"""

import os
import sys
import subprocess
import time
import json
import argparse
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass, asdict
import sqlite3
from datetime import datetime
import concurrent.futures

@dataclass
class JobStatus:
    """Track status of individual protein processing jobs"""
    protein_id: str
    slurm_job_id: Optional[str] = None
    status: str = "pending"  # pending, submitted, running, completed, failed
    submitted_at: Optional[str] = None
    completed_at: Optional[str] = None
    runtime_seconds: Optional[float] = None
    error_message: Optional[str] = None
    domains_found: Optional[int] = None
    output_file: Optional[str] = None

class ProductionBatchProcessor:
    """Scale mini_pyecod to production workloads"""
    
    def __init__(self, config_file: str = "config/production.yml"):
        self.config = self._load_config(config_file)
        self.status_db = self._init_status_database()
        self.mini_executable = self._find_mini_executable()
        
    def _load_config(self, config_file: str) -> Dict:
        """Load production configuration"""
        # Default configuration
        default_config = {
            "slurm": {
                "partition": "batch",
                "time": "1:00:00",
                "memory": "4G",
                "cpus": 1,
                "max_concurrent_jobs": 100
            },
            "paths": {
                "output_dir": "/tmp/mini_production_results",
                "log_dir": "/tmp/mini_production_logs",
                "job_scripts_dir": "/tmp/mini_production_jobs"
            },
            "processing": {
                "batch_size": 50,
                "retry_failed": True,
                "max_retries": 2,
                "timeout_hours": 2
            }
        }
        
        # TODO: Load from YAML file if it exists
        return default_config
    
    def _init_status_database(self) -> sqlite3.Connection:
        """Initialize SQLite database for job tracking"""
        db_path = "/tmp/mini_production_status.db"
        conn = sqlite3.connect(db_path, check_same_thread=False)
        
        conn.execute("""
            CREATE TABLE IF NOT EXISTS job_status (
                protein_id TEXT PRIMARY KEY,
                slurm_job_id TEXT,
                status TEXT,
                submitted_at TEXT,
                completed_at TEXT,
                runtime_seconds REAL,
                error_message TEXT,
                domains_found INTEGER,
                output_file TEXT
            )
        """)
        conn.commit()
        return conn
    
    def _find_mini_executable(self) -> str:
        """Find the mini_pyecod executable"""
        candidates = [
            "./pyecod_mini",
            "mini/pyecod_mini", 
            "../mini/pyecod_mini",
            "pyecod_mini"
        ]
        
        for candidate in candidates:
            if os.path.exists(candidate) and os.access(candidate, os.X_OK):
                return os.path.abspath(candidate)
        
        raise FileNotFoundError("Could not find pyecod_mini executable")
    
    def get_pending_proteins(self, source: str = "database", limit: int = None) -> List[str]:
        """Get list of proteins that need processing"""
        
        if source == "test":
            # Test with known working proteins
            test_proteins = [
                "8ovp_A", "8oni_L", "8p6i_L", "8p8o_H", "8oz3_B", "8p12_L",
                "8p2e_B", "8olg_A", "8p49_A"
            ]
            return test_proteins[:limit] if limit else test_proteins
        
        elif source == "database":
            # TODO: Query production database for unprocessed proteins
            # For now, return a test batch
            print("WARNING: Database integration not yet implemented")
            print("Using test proteins for initial validation")
            return self.get_pending_proteins("test", limit)
        
        else:
            raise ValueError(f"Unknown protein source: {source}")
    
    def create_job_script(self, protein_id: str) -> str:
        """Create SLURM job script for a protein"""
        
        # Ensure directories exist
        job_dir = Path(self.config["paths"]["job_scripts_dir"])
        log_dir = Path(self.config["paths"]["log_dir"])
        output_dir = Path(self.config["paths"]["output_dir"])
        
        for directory in [job_dir, log_dir, output_dir]:
            directory.mkdir(parents=True, exist_ok=True)
        
        # Job script content
        script_content = f"""#!/bin/bash
#SBATCH --partition={self.config["slurm"]["partition"]}
#SBATCH --time={self.config["slurm"]["time"]}
#SBATCH --memory={self.config["slurm"]["memory"]}
#SBATCH --cpus-per-task={self.config["slurm"]["cpus"]}
#SBATCH --job-name=mini_{protein_id}
#SBATCH --output={log_dir}/mini_{protein_id}_%j.out
#SBATCH --error={log_dir}/mini_{protein_id}_%j.err

# Record job start
echo "Job started: $(date)"
echo "Protein: {protein_id}"
echo "Node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"

# Change to mini directory (adjust path as needed)
cd {Path(self.mini_executable).parent}

# Run mini_pyecod
echo "Running: {self.mini_executable} {protein_id}"
{self.mini_executable} {protein_id}

# Check if output was created
OUTPUT_FILE="/tmp/{protein_id}_mini.domains.xml"
if [ -f "$OUTPUT_FILE" ]; then
    echo "Success: Output file created"
    # Copy to production results directory
    cp "$OUTPUT_FILE" {output_dir}/
    echo "Domain processing completed successfully"
else
    echo "Error: No output file created"
    exit 1
fi

echo "Job completed: $(date)"
"""
        
        # Write job script
        script_path = job_dir / f"mini_{protein_id}.sh"
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        # Make executable
        os.chmod(script_path, 0o755)
        
        return str(script_path)
    
    def submit_job(self, protein_id: str) -> str:
        """Submit SLURM job for a protein"""
        
        # Create job script
        script_path = self.create_job_script(protein_id)
        
        # Submit to SLURM
        try:
            result = subprocess.run(
                ["sbatch", script_path],
                capture_output=True,
                text=True,
                check=True
            )
            
            # Parse job ID from output (e.g., "Submitted batch job 12345")
            slurm_job_id = result.stdout.strip().split()[-1]
            
            # Update status database
            job_status = JobStatus(
                protein_id=protein_id,
                slurm_job_id=slurm_job_id,
                status="submitted",
                submitted_at=datetime.now().isoformat()
            )
            self._update_job_status(job_status)
            
            print(f"✓ Submitted {protein_id}: SLURM job {slurm_job_id}")
            return slurm_job_id
            
        except subprocess.CalledProcessError as e:
            error_msg = f"SLURM submission failed: {e.stderr}"
            job_status = JobStatus(
                protein_id=protein_id,
                status="failed",
                error_message=error_msg
            )
            self._update_job_status(job_status)
            print(f"✗ Failed to submit {protein_id}: {error_msg}")
            raise
    
    def _update_job_status(self, job_status: JobStatus):
        """Update job status in database"""
        with self.status_db:
            self.status_db.execute("""
                INSERT OR REPLACE INTO job_status 
                (protein_id, slurm_job_id, status, submitted_at, completed_at, 
                 runtime_seconds, error_message, domains_found, output_file)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                job_status.protein_id,
                job_status.slurm_job_id,
                job_status.status,
                job_status.submitted_at,
                job_status.completed_at,
                job_status.runtime_seconds,
                job_status.error_message,
                job_status.domains_found,
                job_status.output_file
            ))
    
    def process_batch(self, protein_list: List[str], max_concurrent: int = None) -> Dict[str, str]:
        """Process a batch of proteins"""
        
        if max_concurrent is None:
            max_concurrent = self.config["slurm"]["max_concurrent_jobs"]
        
        print(f"Processing batch of {len(protein_list)} proteins")
        print(f"Max concurrent jobs: {max_concurrent}")
        
        # Submit jobs with concurrency control
        submitted_jobs = {}
        active_jobs = 0
        
        for protein_id in protein_list:
            # Wait if we're at max capacity
            while active_jobs >= max_concurrent:
                print(f"At max capacity ({active_jobs}/{max_concurrent}), waiting...")
                time.sleep(10)
                active_jobs = self._count_active_jobs()
            
            try:
                slurm_job_id = self.submit_job(protein_id)
                submitted_jobs[protein_id] = slurm_job_id
                active_jobs += 1
                
                # Brief delay to avoid overwhelming SLURM
                time.sleep(1)
                
            except Exception as e:
                print(f"Failed to submit {protein_id}: {e}")
        
        print(f"✓ Submitted {len(submitted_jobs)} jobs")
        return submitted_jobs
    
    def _count_active_jobs(self) -> int:
        """Count currently active SLURM jobs for this user"""
        try:
            result = subprocess.run(
                ["squeue", "-u", os.getenv("USER", "unknown"), "-h"],
                capture_output=True,
                text=True
            )
            return len(result.stdout.strip().split('\n')) if result.stdout.strip() else 0
        except:
            return 0
    
    def monitor_progress(self, check_interval: int = 30) -> Dict:
        """Monitor progress of submitted jobs"""
        
        print("Monitoring job progress...")
        print("Press Ctrl+C to stop monitoring")
        
        try:
            while True:
                stats = self.get_batch_statistics()
                
                print(f"\n{'='*60}")
                print(f"Batch Progress Report - {datetime.now().strftime('%H:%M:%S')}")
                print(f"{'='*60}")
                print(f"Total jobs: {stats['total']}")
                print(f"Completed: {stats['completed']} ({stats['completed']/stats['total']*100:.1f}%)")
                print(f"Running: {stats['running']}")
                print(f"Failed: {stats['failed']}")
                print(f"Pending: {stats['pending']}")
                
                if stats['completed'] > 0:
                    print(f"Avg runtime: {stats['avg_runtime']:.1f}s")
                    print(f"Success rate: {(stats['completed']/(stats['completed']+stats['failed'])*100):.1f}%")
                
                # Check if all jobs are done
                if stats['running'] == 0 and stats['pending'] == 0:
                    print("\n✓ All jobs completed!")
                    break
                
                time.sleep(check_interval)
                
        except KeyboardInterrupt:
            print("\nMonitoring stopped by user")
            return self.get_batch_statistics()
    
    def get_batch_statistics(self) -> Dict:
        """Get current batch processing statistics"""
        
        cursor = self.status_db.execute("""
            SELECT status, COUNT(*) as count, AVG(runtime_seconds) as avg_runtime
            FROM job_status 
            GROUP BY status
        """)
        
        status_counts = {row[0]: row[1] for row in cursor.fetchall()}
        
        # Update job statuses by checking SLURM
        self._update_job_statuses_from_slurm()
        
        total = sum(status_counts.values())
        
        return {
            'total': total,
            'completed': status_counts.get('completed', 0),
            'running': status_counts.get('running', 0),
            'failed': status_counts.get('failed', 0),
            'pending': status_counts.get('pending', 0),
            'avg_runtime': self._get_average_runtime()
        }
    
    def _update_job_statuses_from_slurm(self):
        """Update job statuses by querying SLURM"""
        # TODO: Query squeue and sacct to update running/completed statuses
        pass
    
    def _get_average_runtime(self) -> float:
        """Get average runtime of completed jobs"""
        cursor = self.status_db.execute("""
            SELECT AVG(runtime_seconds) FROM job_status 
            WHERE status='completed' AND runtime_seconds IS NOT NULL
        """)
        result = cursor.fetchone()[0]
        return result if result else 0.0


def main():
    """Command-line interface for batch processing"""
    parser = argparse.ArgumentParser(
        description='Production Batch Processor for Mini PyECOD',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python batch_processor.py --test-batch 10        # Process 10 test proteins
  python batch_processor.py --batch-size 50        # Process 50 proteins from database  
  python batch_processor.py --status               # Show current status
  python batch_processor.py --monitor              # Monitor progress
        """
    )
    
    parser.add_argument('--test-batch', type=int,
                       help='Process N test proteins for validation')
    parser.add_argument('--batch-size', type=int, default=50,
                       help='Number of proteins to process from database')
    parser.add_argument('--max-jobs', type=int, default=20,
                       help='Maximum concurrent SLURM jobs')
    parser.add_argument('--status', action='store_true',
                       help='Show current batch status and exit')
    parser.add_argument('--monitor', action='store_true',
                       help='Monitor job progress continuously')
    parser.add_argument('--protein-list', type=str,
                       help='File containing list of proteins to process')
    
    args = parser.parse_args()
    
    # Initialize processor
    processor = ProductionBatchProcessor()
    
    if args.status:
        stats = processor.get_batch_statistics()
        print("Current Batch Status:")
        print(f"  Total: {stats['total']}")
        print(f"  Completed: {stats['completed']}")
        print(f"  Running: {stats['running']}")
        print(f"  Failed: {stats['failed']}")
        print(f"  Pending: {stats['pending']}")
        return
    
    if args.monitor:
        processor.monitor_progress()
        return
    
    # Get protein list
    if args.protein_list:
        with open(args.protein_list, 'r') as f:
            proteins = [line.strip() for line in f if line.strip()]
    elif args.test_batch:
        proteins = processor.get_pending_proteins("test", args.test_batch)
    else:
        proteins = processor.get_pending_proteins("database", args.batch_size)
    
    if not proteins:
        print("No proteins to process")
        return
    
    print(f"Processing {len(proteins)} proteins...")
    
    # Process batch
    submitted_jobs = processor.process_batch(proteins, args.max_jobs)
    
    print(f"\n✓ Batch submission completed")
    print(f"Submitted {len(submitted_jobs)} jobs")
    print("\nTo monitor progress:")
    print("  python batch_processor.py --monitor")
    print("  python batch_processor.py --status")


if __name__ == "__main__":
    main()
