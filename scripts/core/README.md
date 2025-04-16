# pyECOD Core Scripts

This directory contains the core scripts responsible for executing the primary functions of the pyECOD pipeline.

## Primary Scripts

### generate_domain_summary_v2.py

The primary script for generating domain summaries with comprehensive XML creation:
- Creates properly structured XML domain summary files
- Handles both chain and domain BLAST results
- Processes special cases like peptide sequences
- Includes robust file path resolution
- Integrates with database for file tracking

**Usage:**
```bash
python generate_domain_summary_v2.py --config config/config.yml --batch-id 123 --protein-id 456 --blast-only
```

### generate_batch_domain_summaries.py

Batch processing script for generating domain summaries for all proteins in a batch:
- Processes multiple proteins in parallel (optional)
- Tracks completion and updates batch status
- Handles error aggregation
- Provides batch-level statistics

**Usage:**
```bash
python generate_batch_domain_summaries.py --config config/config.yml --batch-id 123 --threads 4 --blast-only
```

### run_domain_analysis.py

Executes the domain analysis pipeline on a specified batch:
- Can run with only BLAST results (using --blast-only flag)
- Provides batch-level status tracking
- Validates database records before processing

**Usage:**
```bash
python run_domain_analysis.py --config config/config.yml --batch-id 123 --blast-only --log-file logs/batch_123.log
```

### run_domain_analysis_v2.py

Enhanced domain analysis script with additional options:
- Supports targeted testing of the partition model
- Provides partition-only mode for batches with completed summaries
- Offers validation-only mode to check database and file structures
- Includes detailed error reporting

**Usage:**
```bash
python run_domain_analysis_v2.py --config config/config.yml --batch-id 123 --partition-only --log-file logs/batch_123.log
```

### run_hhsearch.py

Executes the HHSearch pipeline on proteins with BLAST results:
- Generates HHblits profiles for more sensitive search
- Runs HHSearch against the reference database
- Handles job submission and tracking
- Suitable for proteins without sufficient BLAST hits

**Usage:**
```bash
python run_hhsearch.py --config config/config.yml --batch-id 123 --threads 8 --memory 16G
```

### run_indexed_batches.py

Processes indexed batches for domain summary generation:
- Identifies indexed batches in the database
- Checks existing summary status before processing
- Can skip completed batches with --skip-complete flag
- Provides batch filtering options

**Usage:**
```bash
python run_indexed_batches.py --config config/config.yml --blast-only --threads 8 --skip-complete
```

### run_pipeline.py

Orchestrates the full pyECOD processing pipeline:
- Manages the complete workflow from FASTA to domains
- Handles batch creation and registration
- Coordinates BLAST and HHSearch pipelines
- Provides job tracking and status updates

**Usage:**
```bash
python run_pipeline.py --config config/config.yml --limit 100 --batch-size 10 --threads 8
```

## Common Parameters

- `--config` - Path to configuration file (default: config/config.yml)
- `--batch-id` - ID of the batch to process
- `--blast-only` - Use only BLAST results (no HHSearch)
- `--verbose` / `-v` - Enable verbose output
- `--log-file` - Specify log file location
- `--threads` - Number of threads for parallel processing