#!/usr/bin/env python3
"""
Filesystem-First Batch Processor for Mini PyECOD
Scans batch directories, checks for missing mini results, submits SLURM jobs

Usage:
    python filesystem_batch_processor.py --scan-all                    # Scan all batches
    python filesystem_batch_processor.py --batch-name batch_036       # Specific batch
    python filesystem_batch_processor.py --test-proteins 10           # Test with 10 proteins
    python filesystem_batch_processor.py --status                     # Check status
    python filesystem_batch_processor.py --monitor                    # Monitor progress
"""

import os
import sys
import subprocess
import time
import json
import argparse
import sqlite3
import yaml
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Set
from dataclasses import dataclass, asdict
from datetime import datetime
import glob
import psycopg2
import psycopg2.extras

@dataclass
class ProteinJob:
    """Information about a protein processing job"""
    protein_id: str           # e.g., "8ovp_A"
    batch_name: str          # e.g., "batch_036_20250406_1424"  
    domain_summary_path: str # path to input .domain_summary.xml
    mini_output_path: str    # expected path to .mini.domains.xml
    slurm_job_id: Optional[str] = None
    status: str = "pending"   # pending, submitted, running, completed, failed
    submitted_at: Optional[str] = None
    completed_at: Optional[str] = None
    error_message: Optional[str] = None
    is_representative: bool = False

class FilesystemBatchProcessor:
    """Filesystem-first approach to scaling mini_pyecod"""
    
    def __init__(self, config_path: str = "config/config.local.yml"):
        self.config = self._load_config(config_path)
        self.status_db = self._init_tracking_database()
        self.mini_executable = self._find_mini_executable()
        self.db_conn = None

        # Initialize database connection if configured
        if self.config.get('representatives', {}).get('use_database'):
            self._init_db_connection()

    def _load_config(self, config_path: str) -> Dict:
        """Load configuration from YAML file"""
        try:
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)
                print(f"✓ Loaded config from {config_path}")
                return config
        except FileNotFoundError:
            print(f"⚠️  Config file {config_path} not found, using defaults")
            return self._default_config()

    def _default_config(self) -> Dict:
        """Default configuration"""
        return {
            "slurm": {
                "partition": "All",
                "time": "1:00:00",
                "memory": "4G",
                "cpus": 1,
                "max_concurrent_jobs": 50
            },
            "paths": {
                "batch_base_dir": "/data/ecod/pdb_updates/batches",
                "mini_executable": "./mini/pyecod_mini",
                "output_dir": "/tmp/mini_production_results",
                "log_dir": "/tmp/mini_production_logs",
                "job_scripts_dir": "/tmp/mini_production_jobs",
                "tracking_db": "/tmp/mini_production_status.db"
            },
            "processing": {
                "batch_size": 50,
                "scan_pattern": "*.develop291.domain_summary.xml",
                "output_pattern": "{protein_id}.mini.domains.xml"
            }
        }

    def _init_tracking_database(self) -> sqlite3.Connection:
        """Initialize SQLite tracking database"""
        db_path = self.config["paths"]["tracking_db"]
        conn = sqlite3.connect(db_path, check_same_thread=False)

        conn.execute("""
            CREATE TABLE IF NOT EXISTS protein_jobs (
                protein_id TEXT,
                batch_name TEXT,
                domain_summary_path TEXT,
                mini_output_path TEXT,
                slurm_job_id TEXT,
                status TEXT,
                submitted_at TEXT,
                completed_at TEXT,
                error_message TEXT,
                is_representative BOOLEAN,
                PRIMARY KEY (protein_id, batch_name)
            )
        """)
        conn.commit()
        return conn

    def _init_db_connection(self):
        """Initialize PostgreSQL connection for representative filtering"""
        try:
            db_config = self.config['database']
            self.db_conn = psycopg2.connect(**db_config)
            print("✓ Connected to PostgreSQL database")
        except Exception as e:
            print(f"⚠️  Database connection failed: {e}")
            print("Will use filesystem-only mode")

    def _find_mini_executable(self) -> str:
        """Find the mini_pyecod executable"""
        candidates = [
            self.config["paths"]["mini_executable"],
            "./pyecod_mini",
            "mini/pyecod_mini",
            "../mini/pyecod_mini"
        ]

        for candidate in candidates:
            if os.path.exists(candidate) and os.access(candidate, os.X_OK):
                abs_path = os.path.abspath(candidate)
                print(f"✓ Found mini executable: {abs_path}")
                return abs_path

        raise FileNotFoundError("Could not find pyecod_mini executable")

    def scan_batch_directories(self, batch_name_filter: Optional[str] = None) -> List[ProteinJob]:
        """Scan batch directories for proteins needing processing"""

        batch_base = Path(self.config["paths"]["batch_base_dir"])
        scan_pattern = self.config["processing"]["scan_pattern"]

        if not batch_base.exists():
            raise FileNotFoundError(f"Batch base directory not found: {batch_base}")

        print(f"🔍 Scanning batch directories in {batch_base}")

        # Find batch directories
        if batch_name_filter:
            batch_dirs = [batch_base / batch_name_filter]
        else:
            batch_dirs = [d for d in batch_base.iterdir() if d.is_dir()]

        protein_jobs = []

        for batch_dir in batch_dirs:
            batch_name = batch_dir.name
            domains_dir = batch_dir / "domains"
            mini_domains_dir = batch_dir / "mini_domains"

            if not domains_dir.exists():
                print(f"⚠️  No domains directory in {batch_name}")
                continue

            print(f"📁 Scanning batch {batch_name}")

            # Find domain summary files
            summary_files = list(domains_dir.glob(scan_pattern))
            print(f"   Found {len(summary_files)} domain summary files")

            for summary_file in summary_files:
                # Extract protein ID from filename
                # e.g., "8ovp_A.develop291.domain_summary.xml" -> "8ovp_A"
                protein_id = summary_file.name.split('.')[0]

                # Check if mini output already exists
                mini_output_file = mini_domains_dir / f"{protein_id}.mini.domains.xml"

                if mini_output_file.exists():
                    continue  # Already processed

                # Create job record
                protein_job = ProteinJob(
                    protein_id=protein_id,
                    batch_name=batch_name,
                    domain_summary_path=str(summary_file),
                    mini_output_path=str(mini_output_file)
                )
                protein_jobs.append(protein_job)

        print(f"✓ Found {len(protein_jobs)} proteins needing processing")
        return protein_jobs

    def filter_representatives(self, protein_jobs: List[ProteinJob]) -> List[ProteinJob]:
        """Filter to representative proteins using database"""

        if not self.db_conn:
            print("⚠️  No database connection, returning all proteins")
            return protein_jobs

        try:
            # Get representative protein IDs
            with self.db_conn.cursor() as cursor:
                cursor.execute("""
                    SELECT DISTINCT p.source_id
                    FROM ecod_schema.process_status ps
                    JOIN ecod_schema.protein p ON ps.protein_id = p.id
                    WHERE ps.is_representative = TRUE
                """)

                rep_source_ids = {row[0] for row in cursor.fetchall()}
                print(f"📊 Found {len(rep_source_ids)} representative proteins in database")

        except Exception as e:
            print(f"⚠️  Error querying representatives: {e}")
            return protein_jobs

        # Filter protein jobs
        filtered_jobs = []
        for job in protein_jobs:
            if job.protein_id in rep_source_ids:
                job.is_representative = True
                filtered_jobs.append(job)

        print(f"✓ Filtered to {len(filtered_jobs)} representative proteins")
        return filtered_jobs

    def create_job_script(self, protein_job: ProteinJob) -> str:
        """Create SLURM job script for a protein"""

        # Ensure directories exist
        job_dir = Path(self.config["paths"]["job_scripts_dir"])
        log_dir = Path(self.config["paths"]["log_dir"])
        job_dir.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)

        # Ensure mini_domains output directory exists
        mini_output_dir = Path(protein_job.mini_output_path).parent
        mini_output_dir.mkdir(parents=True, exist_ok=True)

        # SLURM job script with correct syntax for your cluster
        script_content = f"""#!/bin/bash
#SBATCH --partition={self.config["slurm"]["partition"]}
#SBATCH --time={self.config["slurm"]["time"]}
#SBATCH --mem={self.config["slurm"]["memory"]}
#SBATCH --cpus-per-task={self.config["slurm"]["cpus"]}
#SBATCH --job-name=mini_{protein_job.protein_id}
#SBATCH --output={log_dir}/mini_{protein_job.protein_id}_%j.out
#SBATCH --error={log_dir}/mini_{protein_job.protein_id}_%j.err

echo "=== Mini PyECOD Job ==="
echo "Job started: $(date)"
echo "Protein: {protein_job.protein_id}"
echo "Batch: {protein_job.batch_name}"
echo "Node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Domain summary: {protein_job.domain_summary_path}"
echo "Output target: {protein_job.mini_output_path}"

# Verify input file exists
if [ ! -f "{protein_job.domain_summary_path}" ]; then
    echo "ERROR: Domain summary file not found"
    exit 1
fi

# Change to mini directory 
cd {Path(self.mini_executable).parent}
echo "Working directory: $(pwd)"

# Run mini_pyecod with the specific input file
echo "Running: {self.mini_executable} {protein_job.protein_id}"
echo "Input file: {protein_job.domain_summary_path}"

# Execute mini
{self.mini_executable} {protein_job.protein_id}
EXIT_CODE=$?

echo "Mini exit code: $EXIT_CODE"

# Check if output was created (mini might create it in a temp location)
TEMP_OUTPUT="/tmp/{protein_job.protein_id}_mini.domains.xml"
if [ -f "$TEMP_OUTPUT" ]; then
    echo "✓ Found output file: $TEMP_OUTPUT"
    # Move to final location
    mv "$TEMP_OUTPUT" "{protein_job.mini_output_path}"
    echo "✓ Moved to: {protein_job.mini_output_path}"
elif [ -f "{protein_job.mini_output_path}" ]; then
    echo "✓ Output file already at target location"
else
    echo "✗ ERROR: No output file created"
    exit 1
fi

echo "Job completed: $(date)"
exit $EXIT_CODE
"""
        
        # Write script
        script_path = job_dir / f"mini_{protein_job.protein_id}_{protein_job.batch_name}.sh"
        with open(script_path, 'w') as f:
            f.write(script_content)
        os.chmod(script_path, 0o755)
        
        return str(script_path)
    
    def submit_job(self, protein_job: ProteinJob) -> str:
        """Submit SLURM job for a protein"""
        
        # Create job script
        script_path = self.create_job_script(protein_job)
        
        try:
            # Submit to SLURM
            result = subprocess.run(
                ["sbatch", script_path],
                capture_output=True,
                text=True,
                check=True
            )
            
            # Parse job ID
            slurm_job_id = result.stdout.strip().split()[-1]
            
            # Update protein job
            protein_job.slurm_job_id = slurm_job_id
            protein_job.status = "submitted"
            protein_job.submitted_at = datetime.now().isoformat()
            
            # Save to database
            self._save_protein_job(protein_job)
            
            print(f"✓ {protein_job.protein_id} -> SLURM job {slurm_job_id}")
            return slurm_job_id
            
        except subprocess.CalledProcessError as e:
            error_msg = f"SLURM submission failed: {e.stderr}"
            protein_job.status = "failed"
            protein_job.error_message = error_msg
            self._save_protein_job(protein_job)
            print(f"✗ {protein_job.protein_id}: {error_msg}")
            raise
    
    def _save_protein_job(self, protein_job: ProteinJob):
        """Save to PostgreSQL instead of SQLite"""
        with self.db_conn.cursor() as cursor:
            cursor.execute("""
                INSERT INTO mini_production.tracking_status
                (protein_id, slurm_job_id, status, submitted_at)
                VALUES (%s, %s, %s, NOW())
                ON CONFLICT (protein_id) DO UPDATE SET
                    slurm_job_id = EXCLUDED.slurm_job_id,
                    status = EXCLUDED.status,
                    submitted_at = EXCLUDED.submitted_at
            """, (protein_job.protein_id, protein_job.slurm_job_id, 'submitted'))
            self.db_conn.commit()
    
    def process_proteins(self, protein_jobs: List[ProteinJob], 
                        max_concurrent: int = None) -> Dict[str, str]:
        """Process a list of proteins with concurrency control"""
        
        if max_concurrent is None:
            max_concurrent = self.config["slurm"]["max_concurrent_jobs"]
        
        print(f"🚀 Processing {len(protein_jobs)} proteins")
        print(f"   Max concurrent: {max_concurrent}")
        
        submitted = {}
        for i, protein_job in enumerate(protein_jobs):
            # Check concurrency
            while self._count_active_jobs() >= max_concurrent:
                print(f"   At capacity, waiting... ({i+1}/{len(protein_jobs)})")
                time.sleep(10)
            
            try:
                slurm_job_id = self.submit_job(protein_job)
                submitted[protein_job.protein_id] = slurm_job_id
                time.sleep(1)  # Brief delay
                
            except Exception as e:
                print(f"✗ Failed {protein_job.protein_id}: {e}")
        
        print(f"✓ Submitted {len(submitted)}/{len(protein_jobs)} jobs")
        return submitted
    
    def _count_active_jobs(self) -> int:
        """Count active SLURM jobs for current user"""
        try:
            result = subprocess.run(
                ["squeue", "-u", os.getenv("USER", "unknown"), "-h"],
                capture_output=True, text=True
            )
            if result.stdout.strip():
                return len(result.stdout.strip().split('\n'))
            return 0
        except:
            return 0
    
    def get_statistics(self) -> Dict:
        """Get current processing statistics"""
        cursor = self.status_db.execute("""
            SELECT status, COUNT(*) as count 
            FROM protein_jobs 
            GROUP BY status
        """)
        
        status_counts = {row[0]: row[1] for row in cursor.fetchall()}
        total = sum(status_counts.values())
        
        return {
            'total': total,
            'completed': status_counts.get('completed', 0),
            'submitted': status_counts.get('submitted', 0),
            'failed': status_counts.get('failed', 0),
            'pending': status_counts.get('pending', 0)
        }
    
    def monitor_progress(self, interval: int = 30):
        """Monitor job progress"""
        print("📊 Monitoring progress (Ctrl+C to stop)")
        
        try:
            while True:
                stats = self.get_statistics()
                print(f"\n{datetime.now().strftime('%H:%M:%S')} - Progress:")
                print(f"  Total: {stats['total']}")
                print(f"  Completed: {stats['completed']}")
                print(f"  Running: {stats['submitted']}")
                print(f"  Failed: {stats['failed']}")
                
                if stats['total'] > 0:
                    pct = (stats['completed'] + stats['failed']) / stats['total'] * 100
                    print(f"  Progress: {pct:.1f}%")
                
                if stats['submitted'] == 0:
                    print("✓ All jobs completed!")
                    break
                
                time.sleep(interval)
                
        except KeyboardInterrupt:
            print("\n⏹️  Monitoring stopped")


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Filesystem-First Batch Processor for Mini PyECOD'
    )
    
    parser.add_argument('--scan-all', action='store_true',
                       help='Scan all batch directories')
    parser.add_argument('--batch-name', type=str,
                       help='Process specific batch')
    parser.add_argument('--test-proteins', type=int,
                       help='Process N test proteins (override scanning)')
    parser.add_argument('--reps-only', action='store_true',
                       help='Process only representative proteins')
    parser.add_argument('--max-jobs', type=int, default=20,
                       help='Maximum concurrent jobs')
    parser.add_argument('--status', action='store_true',
                       help='Show current status')
    parser.add_argument('--monitor', action='store_true',
                       help='Monitor progress')
    
    args = parser.parse_args()
    
    # Initialize processor
    processor = FilesystemBatchProcessor()
    
    if args.status:
        stats = processor.get_statistics()
        print("📊 Current Status:")
        for key, value in stats.items():
            print(f"   {key}: {value}")
        return
    
    if args.monitor:
        processor.monitor_progress()
        return
    
    # Get proteins to process
    if args.test_proteins:
        # Create test protein jobs
        test_proteins = ["8ovp_A", "8oni_L", "8p6i_L"][:args.test_proteins]
        protein_jobs = []
        for pid in test_proteins:
            protein_jobs.append(ProteinJob(
                protein_id=pid,
                batch_name="test_batch",
                domain_summary_path=f"/data/ecod/pdb_updates/batches/batch_036/domains/{pid}.develop291.domain_summary.xml",
                mini_output_path=f"/tmp/mini_test/{pid}.mini.domains.xml"
            ))
        print(f"🧪 Created {len(protein_jobs)} test jobs")
    else:
        # Scan filesystem
        protein_jobs = processor.scan_batch_directories(args.batch_name)
        
        # Filter to representatives if requested
        if args.reps_only:
            protein_jobs = processor.filter_representatives(protein_jobs)
    
    if not protein_jobs:
        print("No proteins to process")
        return
    
    # Process proteins
    submitted = processor.process_proteins(protein_jobs, args.max_jobs)
    
    print(f"\n✅ Batch processing completed")
    print(f"   Jobs submitted: {len(submitted)}")
    print("\nTo monitor:")
    print("   python filesystem_batch_processor.py --monitor")


if __name__ == "__main__":
    main()
