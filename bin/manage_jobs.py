#!/usr/bin/env python3
"""
Job Management Utility for pyECOD
"""

import os
import sys
import argparse
import logging
import time
from typing import Dict, Any, Optional

# Add parent directory to path to allow imports from ecod package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

from ecod.core.slurm_job_manager import SlurmJobManager
from ecod.config import ConfigManager

def setup_logging(verbose: bool = False):
    """Configure logging
    
    Args:
        verbose: Enable verbose logging
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    )

def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='ECOD Job Management Utility')
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
    parser.add_argument('--list', action='store_true',
                      help='List pending jobs')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose)
    
    job_manager = SlurmJobManager(args.config)
    
    if args.submit:
        # Submit a job
        try:
            job_id = job_manager.submit_job(args.submit)
            if job_id:
                print(f"Submitted job with ID: {job_id}")
            else:
                print("Failed to submit job")
                return 1
        except Exception as e:
            print(f"Error submitting job: {str(e)}")
            return 1
    
    elif args.cancel:
        # Cancel a job
        try:
            success = job_manager.cancel_job(args.cancel)
            if success:
                print(f"Cancelled job {args.cancel}")
            else:
                print(f"Failed to cancel job {args.cancel}")
                return 1
        except Exception as e:
            print(f"Error cancelling job: {str(e)}")
            return 1
    
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
            try:
                completed, failed, running = job_manager.check_all_jobs(args.batch_id)
                print(f"Completed: {completed}, Failed: {failed}, Running: {running}")
            except Exception as e:
                print(f"Error checking jobs: {str(e)}")
                return 1
    
    elif args.list:
        # List pending jobs
        try:
            jobs = job_manager.get_pending_jobs(args.batch_id)
            if not jobs:
                print("No pending jobs found")
            else:
                print(f"Found {len(jobs)} pending jobs:")
                for job in jobs:
                    print(f"  Job ID: {job['id']}, Slurm Job ID: {job['slurm_job_id']}, Type: {job['job_type']}, Batch: {job['batch_id']}")
        except Exception as e:
            print(f"Error listing jobs: {str(e)}")
            return 1
    
    else:
        parser.print_help()
    
    return 0

if __name__ == "__main__":
    sys.exit(main())