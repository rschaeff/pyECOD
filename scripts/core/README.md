# pyECOD Core Domain Summary Scripts

This directory contains the core scripts responsible for generating domain summaries in the pyECOD pipeline.

## Scripts

### generate_domain_summary_v2.py

The primary script for generating domain summaries with comprehensive XML creation. This script:

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

Batch processing script for generating domain summaries for all proteins in a batch. This script:

- Processes multiple proteins in parallel (optional)
- Tracks completion and updates batch status
- Handles error aggregation
- Provides batch-level statistics

**Usage:**
```bash
python generate_batch_domain_summaries.py --config config/config.yml --batch-id 123 --threads 4 --blast-only
```

## Common Parameters

- `--config` - Path to configuration file
- `--batch-id` - ID of the batch to process
- `--blast-only` - Generate blast-only summaries (skip HHSearch)
- `--output-dir` - Override output directory
- `--verbose` / `-v` - Enable verbose output
- `--log-file` - Specify log file location

## Integration

These scripts integrate with the broader pyECOD architecture through:

1. The ApplicationContext pattern for configuration and database access
2. The standard database schema for tracking processing status
3. The XML-based domain summary format