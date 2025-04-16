#!/usr/bin/env python3
"""
Script to submit multiple repopulate_ecod_schema.py jobs as a SLURM array job
"""

import os
import sys
import argparse
import subprocess
from ecod.jobs import SlurmJobManager

def main():
    parser = argparse.ArgumentParser(description="Submit ECOD batch jobs to SLURM as an array job")
    parser.add_argument("--config", default="config/config.yml", help="Path to configuration file")
    parser.add_argument("--output-dir", default="./slurm_jobs", help="Directory for job outputs")
    parser.add_argument("--threads", default=4, type=int, help="Number of threads per job")
    parser.add_argument("--memory", default="8G", help="Memory per job")
    parser.add_argument("--time", default="24:00:00", help="Time limit per job")
    parser.add_argument("--start-batch", default=1, type=int, help="Starting batch number")
    parser.add_argument("--end-batch", default=49, type=int, help="Ending batch number")
    parser.add_argument("--date-code", default="20250406_1424", help="Date code for batch names")
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Create array script
    array_script_path = os.path.join(args.output_dir, "ecod_array_job.sh")
    
    array_script_content = """#!/bin/bash
#SBATCH --job-name=ecod_batch_array
#SBATCH --output={output_dir}/ecod_batch_%A_%a.out
#SBATCH --error={output_dir}/ecod_batch_%A_%a.err
#SBATCH --cpus-per-task={threads}
#SBATCH --mem={memory}
#SBATCH --time={time}
#SBATCH --array={array_range}

# Error handling
set -e
trap 'echo "Error occurred at line $LINENO"' ERR

# Format the batch number with leading zeros
BATCH_NUM=$(printf "%03d" $SLURM_ARRAY_TASK_ID)

# Set the batch name
BATCH_NAME="ecod_batch_${{BATCH_NUM}}_{date_code}"

# Run the script
echo "Processing batch: $BATCH_NAME"
python scripts/repopulate_ecod_schema.py --batch-name $BATCH_NAME

# Exit with the script's exit code
exit $?
""".format(
        output_dir=args.output_dir,
        threads=args.threads,
        memory=args.memory,
        time=args.time,
        array_range=f"{args.start_batch}-{args.end_batch}",
        date_code=args.date_code
    )
    
    # Write array script
    with open(array_script_path, 'w') as f:
        f.write(array_script_content)
    os.chmod(array_script_path, 0o755)
    
    # Submit the array job
    try:
        result = subprocess.run(
            ['sbatch', array_script_path],
            capture_output=True,
            text=True,
            check=True
        )
        job_id = result.stdout.strip().split()[-1]
        print(f"Successfully submitted array job {job_id} for batches {args.start_batch} to {args.end_batch}")
    except subprocess.CalledProcessError as e:
        print(f"Error submitting array job: {e.stderr}")
        sys.exit(1)

if __name__ == "__main__":
    main()