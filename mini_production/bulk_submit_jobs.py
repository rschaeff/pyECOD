#!/usr/bin/env python3
"""
Bulk Job Submission Script

Submit ALL jobs at once, let SLURM handle concurrency and dependencies.
This script runs quickly and exits - no long-running process needed.

Usage:
    python bulk_submit_jobs.py --scan-all --reps-only --max-concurrent 40
"""

import argparse
import subprocess
from pathlib import Path
from typing import List
import time

def submit_all_jobs_bulk(protein_jobs: List, max_concurrent: int = 40):
    """Submit all jobs at once with SLURM concurrency control"""
    
    print(f"🚀 Bulk submitting {len(protein_jobs)} jobs")
    print(f"   SLURM will manage concurrency (max {max_concurrent})")
    
    # Create job arrays or submit with dependencies
    job_ids = []
    
    for i, protein_job in enumerate(protein_jobs):
        # Create job script (same as before)
        script_path = create_job_script(protein_job)
        
        # Submit with dependency on previous jobs to limit concurrency
        submit_cmd = ["sbatch"]
        
        # Limit concurrency by creating dependencies
        if len(job_ids) >= max_concurrent:
            # Wait for some jobs to finish before submitting more
            dependency_jobs = job_ids[-max_concurrent:]
            depend_str = f"afterany:{':'.join(dependency_jobs)}"
            submit_cmd.extend(["--dependency", depend_str])
        
        submit_cmd.append(script_path)
        
        try:
            result = subprocess.run(submit_cmd, capture_output=True, text=True, check=True)
            job_id = result.stdout.strip().split()[-1]
            job_ids.append(job_id)
            
            # Update tracking
            save_protein_job_to_tracking(protein_job, job_id)
            
            if i % 100 == 0:
                print(f"   Submitted {i+1}/{len(protein_jobs)} jobs")
                
        except subprocess.CalledProcessError as e:
            print(f"✗ Failed to submit {protein_job.protein_id}: {e}")
    
    print(f"\n✅ Bulk submission complete!")
    print(f"   Total jobs submitted: {len(job_ids)}")
    print(f"   SLURM job IDs: {job_ids[0]} to {job_ids[-1]}")
    print(f"\nMonitor with:")
    print(f"   squeue -u $USER")
    print(f"   python filesystem_batch_processor.py --monitor")

def submit_with_job_arrays(protein_jobs: List, array_size: int = 50):
    """Alternative: Use SLURM job arrays for better efficiency"""
    
    print(f"🚀 Submitting {len(protein_jobs)} jobs as arrays (size {array_size})")
    
    # Split jobs into arrays
    arrays = [protein_jobs[i:i+array_size] for i in range(0, len(protein_jobs), array_size)]
    
    for array_idx, job_array in enumerate(arrays):
        # Create array job script
        array_script = create_array_job_script(job_array, array_idx)
        
        # Submit array
        submit_cmd = [
            "sbatch", 
            f"--array=0-{len(job_array)-1}",
            array_script
        ]
        
        try:
            result = subprocess.run(submit_cmd, capture_output=True, text=True, check=True)
            array_job_id = result.stdout.strip().split()[-1]
            
            print(f"   Array {array_idx}: {len(job_array)} jobs -> {array_job_id}")
            
            # Update tracking for all jobs in array
            for job_idx, protein_job in enumerate(job_array):
                individual_job_id = f"{array_job_id}_{job_idx}"
                save_protein_job_to_tracking(protein_job, individual_job_id)
                
        except subprocess.CalledProcessError as e:
            print(f"✗ Failed to submit array {array_idx}: {e}")
    
    print(f"\n✅ All arrays submitted!")

def create_array_job_script(protein_jobs: List, array_idx: int) -> str:
    """Create a job array script that processes multiple proteins"""
    
    script_content = f"""#!/bin/bash
#SBATCH --partition=All
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=mini_array_{array_idx}
#SBATCH --output=/tmp/mini_production_logs/array_{array_idx}_%A_%a.out
#SBATCH --error=/tmp/mini_production_logs/array_{array_idx}_%A_%a.err

# Array of protein IDs
PROTEINS=(
{chr(10).join(f'    "{job.protein_id}"' for job in protein_jobs)}
)

# Array of input paths  
INPUTS=(
{chr(10).join(f'    "{job.domain_summary_path}"' for job in protein_jobs)}
)

# Array of output paths
OUTPUTS=(
{chr(10).join(f'    "{job.mini_output_path}"' for job in protein_jobs)}
)

# Get this job's assignments
PROTEIN_ID=${{PROTEINS[$SLURM_ARRAY_TASK_ID]}}
INPUT_FILE=${{INPUTS[$SLURM_ARRAY_TASK_ID]}}
OUTPUT_FILE=${{OUTPUTS[$SLURM_ARRAY_TASK_ID]}}

echo "=== Mini PyECOD Array Job ==="
echo "Array Job ID: $SLURM_ARRAY_JOB_ID"
echo "Task ID: $SLURM_ARRAY_TASK_ID" 
echo "Protein: $PROTEIN_ID"
echo "Input: $INPUT_FILE"
echo "Output: $OUTPUT_FILE"

# Verify input exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    exit 1
fi

# Run mini
cd {Path("./mini").absolute()}
./pyecod_mini $PROTEIN_ID

# Move output
TEMP_OUTPUT="/tmp/${{PROTEIN_ID}}_mini.domains.xml"
if [ -f "$TEMP_OUTPUT" ]; then
    mkdir -p $(dirname "$OUTPUT_FILE")
    mv "$TEMP_OUTPUT" "$OUTPUT_FILE"
    echo "✓ Success: $OUTPUT_FILE"
else
    echo "✗ Error: No output created"
    exit 1
fi
"""
    
    script_path = f"/tmp/mini_production_jobs/array_{array_idx}.sh"
    Path(script_path).parent.mkdir(parents=True, exist_ok=True)
    
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    Path(script_path).chmod(0o755)
    return script_path

def main():
    """Bulk submission main"""
    parser = argparse.ArgumentParser(description='Bulk Submit Mini PyECOD Jobs')
    
    parser.add_argument('--scan-all', action='store_true')
    parser.add_argument('--reps-only', action='store_true') 
    parser.add_argument('--max-concurrent', type=int, default=40)
    parser.add_argument('--use-arrays', action='store_true',
                       help='Use job arrays instead of dependencies')
    parser.add_argument('--array-size', type=int, default=50,
                       help='Size of each job array')
    
    args = parser.parse_args()
    
    # Use existing scanner from filesystem_batch_processor
    from filesystem_batch_processor import FilesystemBatchProcessor
    
    processor = FilesystemBatchProcessor()
    protein_jobs = processor.scan_batch_directories()
    
    if args.reps_only:
        protein_jobs = processor.filter_representatives(protein_jobs)
    
    if not protein_jobs:
        print("No jobs to submit")
        return
    
    # Submit jobs
    if args.use_arrays:
        submit_with_job_arrays(protein_jobs, args.array_size)
    else:
        submit_all_jobs_bulk(protein_jobs, args.max_concurrent)

if __name__ == "__main__":
    main()
