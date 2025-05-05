#!/usr/bin/env python3
"""
Run domain partitioning on all batches using SLURM job distribution
"""

import os
import sys
import logging
import argparse
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any, Optional

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.config import ConfigManager
from ecod.core.context import ApplicationContext
from ecod.core.slurm_job_manager import SlurmJobManager
from ecod.db import DBManager


class DomainPartitionSlurmRunner:
    """Distribute domain partitioning tasks across SLURM cluster"""
    
    def __init__(self, config_path: str = None):
        self.logger = logging.getLogger("domain_partition_slurm")
        
        # Create application context and slurm manager
        self.context = ApplicationContext(config_path)
        self.slurm = SlurmJobManager(config_path)
        self.db = self.context.db
        
        # Configuration from config file
        self.config = self.context.config
        
        # Default job resources
        self.default_resources = {
            "threads": 4,
            "memory": "8G", 
            "time": "4:00:00"
        }
        
    def get_all_batches(self, reference: str = None, status: str = None) -> List[Dict[str, Any]]:
        """Get all batches from database"""
        query = """
        SELECT id, batch_name, base_path, ref_version, created_at
        FROM ecod_schema.batch
        """
        
        where_clauses = []
        params = []
        
        if reference:
            where_clauses.append("ref_version = %s")
            params.append(reference)
            
        if status:
            where_clauses.append("status = %s")
            params.append(status)
            
        if where_clauses:
            query += " WHERE " + " AND ".join(where_clauses)
            
        query += " ORDER BY id"
        
        batches = self.db.execute_dict_query(query, tuple(params))
        return batches
    
    def check_batch_readiness(self, batch_id: int, blast_only: bool) -> Dict[str, Any]:
        """Check if batch has the required domain summary files"""
        summary_type = "blast_only_summary" if blast_only else "domain_summary"
        
        query = """
        SELECT
            COUNT(*) as total_proteins,
            SUM(CASE WHEN EXISTS (
                SELECT 1 FROM ecod_schema.process_file pf
                WHERE pf.process_id = ps.id
                AND pf.file_type = %s
                AND pf.file_exists = TRUE
            ) THEN 1 ELSE 0 END) as ready_proteins
        FROM
            ecod_schema.process_status ps
        WHERE
            ps.batch_id = %s
        """
        
        rows = self.db.execute_dict_query(query, (summary_type, batch_id))
        
        if rows:
            total = rows[0]['total_proteins']
            ready = rows[0]['ready_proteins']
            readiness = (ready / total * 100) if total > 0 else 0
            
            return {
                "batch_id": batch_id,
                "total_proteins": total,
                "ready_proteins": ready,
                "readiness_percent": readiness,
                "is_ready": readiness >= 90  # Consider ready if 90%+ have domain summaries
            }
        
        return {"batch_id": batch_id, "is_ready": False, "error": "No proteins found"}
    
    def create_batch_partition_job(self, batch_id: int, output_dir: str, 
                                 blast_only: bool = False, limit: int = None,
                                 resources: Dict[str, Any] = None) -> Dict[str, Any]:
        """Create a SLURM job to run domain partitioning for a batch
        
        Note: Non-blast-only runs automatically process only representative proteins
        """
        
        if resources is None:
            resources = self.default_resources.copy()
        
        # Create command for domain partitioning
        commands = [
            f"cd {self.context.config['paths']['project_root']}",
            f"python scripts/run_domain_partition.py "
            f"--config {self.context.config_path} "
            f"--batch-id {batch_id}"
        ]
        
        if blast_only:
            commands[-1] += " --blast-only"
        else:
            # Non-blast-only should only process representative proteins
            commands[-1] += " --reps-only"
            
        if limit:
            commands[-1] += f" --limit {limit}"
            
        # Add logging
        log_prefix = "blast_only" if blast_only else "full_pipeline"
        log_file = f"{output_dir}/logs/domain_partition_{log_prefix}_batch_{batch_id}.log"
        commands[-1] += f" --log-file {log_file}"
        
        # Create job script
        job_name = f"dp_batch_{batch_id}_{'blast' if blast_only else 'full'}"
        script_path = self.slurm.create_job_script(
            commands=commands,
            job_name=job_name,
            output_dir=output_dir,
            threads=resources['threads'],
            memory=resources['memory'],
            time=resources['time']
        )
        
        return {
            "script_path": script_path,
            "job_name": job_name,
            "commands": commands
        }
    
    def submit_all_batches(self, reference: str = None, blast_only: bool = False,
                          dry_run: bool = False, limit: int = None,
                          check_readiness: bool = True) -> List[Dict[str, Any]]:
        """Submit domain partitioning jobs for all batches
        
        Note: Non-blast-only jobs automatically process only representative proteins
        """
        
        # Create output directory for job scripts and logs
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = f"/tmp/domain_partition_jobs_{timestamp}"
        os.makedirs(f"{output_dir}/logs", exist_ok=True)
        
        self.logger.info(f"Job output directory: {output_dir}")
        
        # Get all batches
        batches = self.get_all_batches(reference=reference)
        self.logger.info(f"Found {len(batches)} total batches")
        
        results = []
        batches_ready = 0
        batches_submitted = 0
        
        for batch in batches:
            batch_id = batch['id']
            
            # Check if batch is ready
            if check_readiness:
                readiness = self.check_batch_readiness(batch_id, blast_only)
                if not readiness['is_ready']:
                    readiness_pct = readiness.get('readiness_percent', 0)
                    self.logger.info(f"Batch {batch_id} not ready "
                                   f"({readiness_pct:.1f}% - skipping)")
                    results.append({
                        "batch_id": batch_id,
                        "status": "not_ready",
                        "readiness": readiness
                    })
                    continue
                else:
                    self.logger.info(f"Batch {batch_id} is ready "
                                   f"({readiness['readiness_percent']:.1f}%)")
                    batches_ready += 1
            
            # Create job script
            job_info = self.create_batch_partition_job(
                batch_id=batch_id,
                output_dir=output_dir,
                blast_only=blast_only,
                limit=limit
            )
            
            # Submit job if not dry run
            if not dry_run:
                try:
                    job_id = self.slurm.submit_job(job_info['script_path'])
                    batches_submitted += 1
                    
                    results.append({
                        "batch_id": batch_id,
                        "status": "submitted",
                        "job_id": job_id,
                        "job_name": job_info['job_name']
                    })
                    
                    self.logger.info(f"Submitted job {job_id} for batch {batch_id}")
                    
                except Exception as e:
                    self.logger.error(f"Failed to submit batch {batch_id}: {e}")
                    results.append({
                        "batch_id": batch_id,
                        "status": "failed",
                        "error": str(e)
                    })
            else:
                results.append({
                    "batch_id": batch_id,
                    "status": "dry_run",
                    "job_script": job_info['script_path']
                })
        
        # Summary
        self.logger.info(f"Summary: {batches_ready} batches ready, "
                       f"{batches_submitted} jobs submitted")
        
        return results
    
    def monitor_jobs(self, job_ids: List[str], interval: int = 60) -> None:
        """Monitor status of submitted jobs"""
        import time
        
        self.logger.info(f"Monitoring {len(job_ids)} jobs...")
        
        completed = 0
        failed = 0
        
        while job_ids:
            # Check job status
            completed_jobs = []
            
            for job_id in job_ids:
                try:
                    status_info = self.slurm.check_job_status(job_id)
                    state = status_info.get('state', 'UNKNOWN')
                    
                    if state in ['COMPLETED', 'CANCELLED', 'FAILED', 'TIMEOUT']:
                        completed_jobs.append(job_id)
                        
                        if state == 'COMPLETED':
                            completed += 1
                            self.logger.info(f"Job {job_id} completed successfully")
                        else:
                            failed += 1
                            self.logger.error(f"Job {job_id} {state}")
                    
                except Exception as e:
                    self.logger.error(f"Error checking job {job_id}: {e}")
                    completed_jobs.append(job_id)
                    failed += 1
            
            # Remove completed jobs from monitoring
            for job_id in completed_jobs:
                job_ids.remove(job_id)
            
            # Status update
            if job_ids:
                self.logger.info(f"Jobs remaining: {len(job_ids)}, "
                               f"Completed: {completed}, Failed: {failed}")
                time.sleep(interval)
        
        # Final summary
        self.logger.info(f"All jobs finished - Completed: {completed}, Failed: {failed}")
    
    def run_in_batches(self, batch_size: int = 10, reference: str = None,
                      blast_only: bool = False, dry_run: bool = False,
                      limit: int = None, monitor: bool = True) -> None:
        """Run domain partitioning in batches to avoid overwhelming the scheduler
        
        Note: Non-blast-only processing automatically restricts to representative proteins
        """
        
        all_batches = self.get_all_batches(reference=reference)
        total_batches = len(all_batches)
        
        for i in range(0, total_batches, batch_size):
            # Get slice of batches to process
            batch_slice = all_batches[i:i + batch_size]
            batch_ids = [b['id'] for b in batch_slice]
            
            self.logger.info(f"Processing batch group {i//batch_size + 1}: "
                           f"batches {batch_ids[0]}-{batch_ids[-1]}")
            
            # Create filter function to process only these batches
            def batch_filter(batch):
                return batch['id'] in batch_ids
            
            # Submit jobs for this group
            results = self.submit_all_batches(
                reference=reference,
                blast_only=blast_only,
                dry_run=dry_run,
                limit=limit
            )
            
            # Extract job IDs from successful submissions
            job_ids = [r['job_id'] for r in results if r.get('status') == 'submitted']
            
            # Monitor this batch group if requested
            if monitor and job_ids and not dry_run:
                self.logger.info(f"Waiting for batch group to complete...")
                self.monitor_jobs(job_ids)
                self.logger.info(f"Batch group {i//batch_size + 1} complete")
                
                # Brief pause between batch groups
                if i + batch_size < total_batches:
                    import time
                    time.sleep(30)  # 30 second pause between groups


def setup_logging(log_file: str = None, verbose: bool = False) -> None:
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    
    handlers = [logging.StreamHandler()]
    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Run domain partitioning on all batches using SLURM')
    
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    
    # Operation modes
    parser.add_argument('--submit', action='store_true',
                      help='Submit jobs for all ready batches')
    parser.add_argument('--monitor', type=str, nargs='+',
                      help='Monitor specific job IDs')
    parser.add_argument('--check-status', action='store_true',
                      help='Check readiness status of all batches')
    
    # Options
    parser.add_argument('--reference', type=str,
                      help='Filter batches by reference version')
    parser.add_argument('--blast-only', action='store_true',
                      help='Use only BLAST results (no HHSearch). Default: process only representative proteins with full pipeline.')
    parser.add_argument('--limit', type=int,
                      help='Limit proteins per batch (for testing)')
    parser.add_argument('--dry-run', action='store_true',
                      help='Create job scripts but dont submit them')
    parser.add_argument('--no-monitor', action='store_true',
                      help='Submit jobs without monitoring')
    parser.add_argument('--batch-size', type=int, default=10,
                      help='Number of batches to process simultaneously')
    
    # Logging options
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.log_file, args.verbose)
    logger = logging.getLogger("domain_partition_slurm")
    
    # Initialize runner
    runner = DomainPartitionSlurmRunner(args.config)
    
    # Handle different operations
    if args.monitor:
        # Monitor specific jobs
        runner.monitor_jobs(args.monitor)
        
    elif args.check_status:
        # Check status of all batches
        batches = runner.get_all_batches(reference=args.reference)
        
        ready_count = 0
        for batch in batches:
            readiness = runner.check_batch_readiness(batch['id'], args.blast_only)
            
            if readiness['is_ready']:
                ready_count += 1
                status = "READY"
            else:
                pct = readiness.get('readiness_percent', 0)
                status = f"NOT READY ({pct:.1f}%)"
            
            logger.info(f"Batch {batch['id']:3d}: {status}")
        
        logger.info(f"Summary: {ready_count}/{len(batches)} batches ready")
        
    elif args.submit:
        # Submit jobs for all ready batches
        runner.run_in_batches(
            batch_size=args.batch_size,
            reference=args.reference,
            blast_only=args.blast_only,
            dry_run=args.dry_run,
            limit=args.limit,
            monitor=not args.no_monitor
        )
        
    else:
        parser.print_help()
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
