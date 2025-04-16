# pyECOD Utilities Scripts

This directory contains utility scripts for maintaining, fixing, and reorganizing domain summaries in the pyECOD pipeline.

## Scripts

### fix_domain_summaries.py

Repairs incorrectly formatted domain summary files by:
- Extracting data from raw BLAST XML files
- Creating properly structured domain summary XML
- Updating database records to point to new files
- Maintaining proper file organization

**Usage:**
```bash
python fix_domain_summaries.py --config config/config.yml --batch-id 123 [--apply]
```

### cleanup_domain_summaries.py

Relocates domain summaries to standardized locations:
- Finds domain summaries across the batch directory
- Moves them to a standard `domains/` subdirectory
- Updates database records with new file paths
- Handles various naming patterns

**Usage:**
```bash
python cleanup_domain_summaries.py --config config/config.yml --batch-id 123 [--apply]
```

### migrate_domain_summaries.py

Migrates domain summaries to a new directory structure:
- Processes files across multiple batches (optional)
- Preserves the original files until migration is confirmed
- Updates database records to point to the new location
- Provides progress tracking with tqdm

**Usage:**
```bash
python migrate_domain_summaries.py --config config/config.yml [--batch-id 123] [--dry-run]
```

### remove_invalid_summaries.py

Removes domain summaries that fail validation:
- Identifies XML files that cannot be properly parsed
- Checks for required elements like blast_summ
- Updates database records to mark files as non-existent
- Can target specific proteins with the `--filter` option

**Usage:**
```bash
python remove_invalid_summaries.py --config config/config.yml --batch-id 123 [--apply] [--filter 8gh6_P]
```

### debug_domain_blast_lookup.py

Debug tool for BLAST file path lookup issues:
- Diagnoses problems with domain_blast_result file paths
- Tests file path resolution methods
- Checks database records against filesystem
- Helps identify mismatches in file types

**Usage:**
```bash
python debug_domain_blast_lookup.py --config config/config.yml --batch-id 123 --protein-id 456 [--apply-fix]
```

## Common Parameters

- `--config` - Path to configuration file
- `--batch-id` - ID of the batch to process
- `--dry-run` / `--apply` - Preview changes or actually apply them
- `--verbose` / `-v` - Enable verbose output
- `--log-file` - Specify log file location

## Features

These utility scripts provide several key capabilities:

1. **File Repair and Maintenance**
   - Fix incorrectly formatted XML files
   - Standardize file locations
   - Clean up invalid files

2. **Database Alignment**
   - Update database records to match actual file locations
   - Ensure consistency between filesystem and database

3. **Debugging Support**
   - Diagnose file path resolution issues
   - Identify missing files and broken links
   - Test resolution algorithms