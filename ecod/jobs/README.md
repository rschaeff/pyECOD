# ECOD Jobs Module

A comprehensive job management system for the ECOD (Evolutionary Classification of Protein Domains) pipeline, supporting both local execution and SLURM cluster environments with database integration.

## Overview

The `ecod.jobs` module provides a unified interface for submitting, monitoring, and managing computational jobs in the ECOD protein domain analysis pipeline. It supports multiple execution backends and includes database integration for tracking job progress across large-scale batch processing operations.

## Features

- **Multiple Execution Backends**: Local threading-based execution and SLURM cluster integration
- **Database Integration**: Track job status, items, and progress in PostgreSQL
- **Batch Processing**: Efficient handling of large sets of computational tasks
- **Job Monitoring**: Real-time status checking and progress tracking
- **Error Handling**: Robust error handling and logging throughout the pipeline
- **Flexible Configuration**: YAML-based configuration with sensible defaults

## Architecture

```
ecod.jobs/
├── base.py           # Abstract JobManager interface
├── local.py          # Local execution using threading
├── slurm.py          # SLURM cluster execution
├── db_job_manager.py # Database-aware job manager wrapper
├── factory.py        # Job manager factory
├── job.py            # Job and JobItem data models
└── __init__.py       # Module exports
```

## Quick Start

### Basic Usage

```python
from ecod.jobs import create_job_manager

# Create a job manager from configuration
config = {
    'job_manager': {'type': 'local'},  # or 'slurm'
    'database': {...}  # database config
}

job_manager = create_job_manager(config)

# Create and submit a simple job
script_path = job_manager.create_job_script(
    commands=["echo 'Hello World'", "sleep 10"],
    job_name="hello_world",
    output_dir="/tmp/jobs"
)

job_id = job_manager.submit_job(script_path)
print(f"Submitted job: {job_id}")

# Check job status
status = job_manager.check_job_status(job_id)
print(f"Job status: {status}")
```

### Database-Integrated Usage

```python
from ecod.jobs import DatabaseJobManager

# Initialize with database integration
db_job_manager = DatabaseJobManager("config/config.yml")

# Submit job with database tracking
job_id = db_job_manager.submit_job(
    script_path,
    batch_id=123,
    job_type="domain_analysis",
    process_ids=[1, 2, 3, 4, 5]  # Track specific process items
)

# Check all jobs for a batch
completed, failed, running = db_job_manager.check_all_jobs(batch_id=123)
print(f"Batch 123: {completed} completed, {failed} failed, {running} running")
```

## Job Manager Types

### LocalJobManager

Executes jobs on the local machine using Python threading. Suitable for development, testing, and small-scale processing.

**Features:**
- Thread-based parallel execution
- Real-time job monitoring
- Local file-based output capture
- No external dependencies

**Configuration:**
```yaml
job_manager:
  type: local
```

### SlurmJobManager

Integrates with SLURM workload manager for cluster execution. Ideal for large-scale computational workloads.

**Features:**
- SLURM job submission with `sbatch`
- Status monitoring with `squeue` and `sacct`
- Resource allocation (CPU, memory, time limits)
- Module loading support
- Job cancellation with `scancel`

**Configuration:**
```yaml
job_manager:
  type: slurm

# Optional: modules to load
modules:
  - bioinfo
  - blast/2.13.0
  - python/3.9
```

**Requirements:**
- SLURM commands (`sbatch`, `squeue`, `sacct`, `scancel`) in PATH
- Appropriate cluster permissions

## Batch Processing

Process large datasets efficiently by grouping items into batches:

```python
# Define items to process
items = [
    (1, "/data/protein1.fasta"),
    (2, "/data/protein2.fasta"),
    (3, "/data/protein3.fasta"),
    # ... more items
]

# Create batch job template
job_template = {
    "name": "protein_analysis",
    "output_dir": "/results/batch_jobs",
    "command_template": lambda item_id, path: f"analyze_protein.py --id {item_id} --input {path}",
    "threads": 4,
    "memory": "8G",
    "time": "02:00:00"
}

# Create batch jobs (2 items per job)
jobs = job_manager.create_batch_jobs(items, batch_size=2, job_template=job_template)

# Submit all jobs
job_ids = []
for job in jobs:
    job_id = job_manager.submit_job(job["script_path"])
    job_ids.append(job_id)
    print(f"Submitted batch {job['batch_num']} with job ID {job_id}")
```

## Database Schema

The module integrates with these database tables:

```sql
-- Job tracking
CREATE TABLE ecod_schema.job (
    id SERIAL PRIMARY KEY,
    batch_id INTEGER,
    job_type VARCHAR,
    slurm_job_id VARCHAR,
    status VARCHAR DEFAULT 'submitted',
    items_count INTEGER,
    submission_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    completion_time TIMESTAMP
);

-- Individual job items
CREATE TABLE ecod_schema.job_item (
    id SERIAL PRIMARY KEY,
    job_id INTEGER REFERENCES ecod_schema.job(id),
    process_id INTEGER,
    status VARCHAR DEFAULT 'pending'
);

-- Process status tracking
CREATE TABLE ecod_schema.process_status (
    id SERIAL PRIMARY KEY,
    batch_id INTEGER,
    status VARCHAR,
    -- other process-specific fields
);
```

## Configuration

### Complete Configuration Example

```yaml
# Database configuration
database:
  host: localhost
  port: 5432
  database: ecod
  user: ecod_user
  password: ecod_pass

# Job manager configuration
job_manager:
  type: slurm  # or 'local'

# SLURM-specific configuration
modules:
  - bioinfo
  - blast/2.13.0
  - python/3.9

# Pipeline configuration
pipeline:
  batch_size: 100
  default_threads: 8
  default_memory: "16G"
  default_time: "04:00:00"
```

### Environment Variables

The module respects these environment variables:

- `SLURM_JOB_ID`: Current SLURM job ID (when running within a job)
- `SLURM_SUBMIT_DIR`: SLURM submission directory
- `SLURM_CLUSTER_NAME`: Cluster name for logging

## API Reference

### JobManager (Abstract Base Class)

The base interface implemented by all job managers.

#### Methods

- `create_job_script(commands, job_name, output_dir, **options) -> str`
- `submit_job(script_path) -> Optional[str]`
- `check_job_status(job_id) -> str`
- `cancel_job(job_id) -> bool`
- `get_job_output(job_id, output_dir, job_name) -> Dict[str, Optional[str]]`
- `create_batch_jobs(items, batch_size, job_template) -> List[Dict[str, Any]]`
- `check_all_jobs(batch_id=None) -> Tuple[int, int, int]`

### DatabaseJobManager

Enhanced job manager with database integration.

#### Additional Methods

- `submit_job(script_path, batch_id=None, job_type="generic", process_ids=None) -> Optional[str]`

### Job Status Values

**SLURM Status Values:**
- `PENDING`: Job is waiting in queue
- `RUNNING`: Job is currently executing
- `COMPLETED`: Job finished successfully
- `FAILED`: Job terminated with error
- `TIMEOUT`: Job exceeded time limit
- `CANCELLED`: Job was cancelled by user
- `NODE_FAIL`: Job failed due to node failure

**Local Status Values:**
- `PENDING`: Job queued for execution
- `RUNNING`: Job is executing
- `COMPLETED`: Job finished successfully
- `FAILED`: Job terminated with error
- `CANCELLED`: Job was cancelled

## Error Handling

The module includes comprehensive error handling:

```python
from ecod.exceptions import JobSubmissionError, JobExecutionError

try:
    job_id = job_manager.submit_job(script_path)
except JobSubmissionError as e:
    logger.error(f"Failed to submit job: {e}")
except JobExecutionError as e:
    logger.error(f"Job execution error: {e}")
```

## Logging

Configure logging to monitor job operations:

```python
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("ecod.jobs")

# The module will log:
# - Job submissions and completions
# - Status changes
# - Error conditions
# - Resource usage (for SLURM)
```

## Best Practices

1. **Resource Allocation**: Specify appropriate CPU, memory, and time limits for SLURM jobs
2. **Batch Size Optimization**: Choose batch sizes that balance efficiency and resource usage
3. **Error Handling**: Always handle exceptions from job operations
4. **Status Monitoring**: Regularly check job status for long-running operations
5. **Output Management**: Ensure adequate disk space for job output files
6. **Database Transactions**: Use database transactions for consistency when tracking multiple jobs

## Integration with ECOD Pipeline

This job management system is designed to work seamlessly with the broader ECOD pipeline:

```python
# Example: Domain analysis pipeline integration
from ecod.pipelines.domain_analysis import DomainAnalysisPipeline
from ecod.jobs import DatabaseJobManager

pipeline = DomainAnalysisPipeline(config_path="config/config.yml")
job_manager = DatabaseJobManager(config_path="config/config.yml")

# Process a batch of proteins
batch_id = 123
proteins = pipeline.get_pending_proteins(batch_id)

# Create analysis jobs
job_template = {
    "name": "domain_analysis",
    "command_template": lambda protein_id, _: f"analyze_domains.py --protein-id {protein_id}",
    "threads": 8,
    "memory": "16G",
    "time": "02:00:00"
}

items = [(p.id, p.sequence_path) for p in proteins]
job_ids = job_manager.create_batch_jobs(items, batch_size=50, job_template=job_template)

# Submit and monitor
for job_id in job_ids:
    job_manager.submit_job(job_id, batch_id=batch_id, job_type="domain_analysis")
```

## Contributing

When extending this module:

1. Implement the `JobManager` abstract interface for new backends
2. Add comprehensive error handling and logging
3. Include unit tests for new functionality
4. Update configuration documentation
5. Follow the existing code style and patterns

## Dependencies

- Python 3.7+
- PostgreSQL (for database integration)
- SLURM (for cluster execution)
- PyYAML (for configuration)
- Standard library: `subprocess`, `threading`, `uuid`, `logging`
