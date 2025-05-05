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
from ecod.jobs import DatabaseJobManager
from ecod.db import DBManager


class DomainPartitionSlurmRunner:
    """Distribute domain partitioning tasks across SLURM cluster using DatabaseJobManager"""

    def __init__(self, config_path: str = None):
        self.logger = logging.getLogger("domain_partition_slurm")

        # Make config path absolute before creating contexts
        if config_path:
            config_path = os.path.abspath(config_path)

        # Create application context
        self.context = ApplicationContext(config_path)

        # Create database-aware job manager
        self.job_manager = DatabaseJobManager(config_path)
        self.db = self.context.db

        # Configuration from config file
        self.config = self.context.config

        # Store absolute path to config for use in job scripts
        self.config_path = config_path

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

    def submit_batch_partition_job(self, batch_id: int, output_dir: str,
                                 blast_only: bool = False, limit: int = None,
                                 resources: Dict[str, Any] = None) -> Optional[str]:
        """Submit a SLURM job to run domain partitioning for a batch"""

        if resources is None:
            resources = self.default_resources.copy()

        # Get project root directory (using current working directory if not configured)
        project_root = self.context.config.get('paths', {}).get('project_root')
        if not project_root:
            # Use current working directory
            project_root = os.getcwd()

        # Create command for domain partitioning
        commands = [
            f"cd {project_root}",
            f"python scripts/run_domain_partition.py "
            f"--config {self.config_path} "
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
        script_path = self.job_manager.create_job_script(
            commands=commands,
            job_name=job_name,
            output_dir=output_dir,
            threads=resources['threads'],
            memory=resources['memory'],
            time=resources['time']
        )

        # Get process IDs for this batch
        query = """
        SELECT id FROM ecod_schema.process_status
        WHERE batch_id = %s
        """
        rows = self.db.execute_query(query, (batch_id,))
        process_ids = [row[0] for row in rows]

        # Submit job with database tracking
        job_id = self.job_manager.submit_job(
            script_path=script_path,
            batch_id=batch_id,
            job_type="domain_partition",
            process_ids=process_ids
        )

        return job_id

    def submit_all_batches(self, reference: str = None, blast_only: bool = False,
                          dry_run: bool = False, limit: int = None,
                          check_readiness: bool = True) -> List[Dict[str, Any]]:
        """Submit domain partitioning jobs for all batches"""

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

            # Submit job if not dry run
            if not dry_run:
                try:
                    job_id = self.submit_batch_partition_job(
                        batch_id=batch_id,
                        output_dir=output_dir,
                        blast_only=blast_only,
                        limit=limit
                    )

                    if job_id:
                        batches_submitted += 1
                        results.append({
                            "batch_id": batch_id,
                            "status": "submitted",
                            "job_id": job_id
                        })
                        self.logger.info(f"Submitted job {job_id} for batch {batch_id}")
                    else:
                        results.append({
                            "batch_id": batch_id,
                            "status": "failed",
                            "error": "Job submission failed"
                        })

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
                    "status": "dry_run"
                })

        # Summary
        self.logger.info(f"Summary: {batches_ready} batches ready, "
                       f"{batches_submitted} jobs submitted")

        return results

    def monitor_jobs(self, batch_id: Optional[int] = None, interval: int = 60) -> None:
        """Monitor status of jobs in database"""
        import time

        while True:
            try:
                # Check all jobs using database job manager
                completed, failed, running = self.job_manager.check_all_jobs(batch_id)

                if running == 0:
                    self.logger.info(f"All jobs completed: {completed} successful, {failed} failed")
                    break
                else:
                    self.logger.info(f"Jobs: {completed} completed, {failed} failed, {running} running")
                    time.sleep(interval)
            except Exception as e:
                self.logger.error(f"Error monitoring jobs: {e}")
                time.sleep(interval)

    def run_in_batches(self, batch_size: int = 10, reference: str = None,
                      blast_only: bool = False, dry_run: bool = False,
                      limit: int = None, monitor: bool = True) -> None:
        """Run domain partitioning in batches to avoid overwhelming the scheduler"""

        all_batches = self.get_all_batches(reference=reference)
        total_batches = len(all_batches)

        for i in range(0, total_batches, batch_size):
            # Get slice of batches to process
            batch_slice = all_batches[i:i + batch_size]
            batch_ids = [b['id'] for b in batch_slice]

            self.logger.info(f"Processing batch group {i//batch_size + 1}: "
                           f"batches {batch_ids[0]}-{batch_ids[-1]}")

            # Submit jobs for this group
            results = []
            for batch in batch_slice:
                if not dry_run:
                    result = self.submit_batch_partition_job(
                        batch_id=batch['id'],
                        output_dir=f"/tmp/domain_partition_jobs_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
                        blast_only=blast_only,
                        limit=limit
                    )
                    results.append(result)
                else:
                    results.append(None)

            # Monitor this batch group if requested
            if monitor and any(results) and not dry_run:
                self.logger.info(f"Waiting for batch group to complete...")
                # Monitor only jobs from this group
                submitted_batch_ids = [batch_ids[j] for j, r in enumerate(results) if r]
                if submitted_batch_ids:
                    self.monitor_jobs(batch_id=None, interval=60)
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
    parser.add_argument('--monitor', type=int, nargs='*',
                      help='Monitor jobs (optionally filter by batch ID)')
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

    # Use absolute path for config if provided
    if args.config and not os.path.isabs(args.config):
        config_path = os.path.abspath(args.config)
    else:
        config_path = args.config

    # Setup logging
    setup_logging(args.log_file, args.verbose)
    logger = logging.getLogger("domain_partition_slurm")

    # Initialize runner
    runner = DomainPartitionSlurmRunner(config_path)

    # Handle different operations
    if args.monitor is not None:
        # Monitor jobs
        batch_id = args.monitor[0] if args.monitor else None
        runner.monitor_jobs(batch_id)
        
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
