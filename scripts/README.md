# pyECOD Scripts

This repository contains scripts for the pyECOD (Evolutionary Classification of Domains) pipeline, specifically focused on domain summary generation and analysis.

## Directory Structure

The scripts are organized into four main directories:

### [core/](core/)
Contains the core scripts for generating domain summaries:
- `generate_domain_summary_v2.py` - Primary script for domain summary generation
- `generate_batch_domain_summaries.py` - Batch processing for domain summary generation

### [utilities/](utilities/)
Utility scripts for maintaining, fixing, and reorganizing domain summaries:
- `fix_domain_summaries.py` - Repairs incorrectly formatted files
- `cleanup_domain_summaries.py` - Relocates files to standardized locations
- `migrate_domain_summaries.py` - Migrates to new directory structure
- `remove_invalid_summaries.py` - Removes files that fail validation
- `debug_domain_blast_lookup.py` - Debug tool for file path issues

### [analysis/](analysis/)
Scripts for analyzing domains, batches, and quality metrics:
- `analyze_domain_summaries.py` - Analyzes domain summary content
- `analyze_batch_domain_summaries.py` - Analyzes success/failure patterns
- `analyze_blast_summaries.py` - Analyzes XML structure
- `analyze_batch_dates.py` - Analyzes PDB structure deposition dates
- `analyze_batch_pdb_ids.py` - Analyzes PDB IDs for structure age

### [inspection/](inspection/)
Scripts for examining individual files and debugging:
- `inspect_domain_summary_xml.py` - Inspects individual XML files
- `debug_domain_summary.py` - Comprehensive domain analysis tool
- `check_domain_summaries_xml.py` - Checks XML files across batches

## Common Parameters

Most scripts support these common parameters:
- `--config` - Path to configuration file (default: config/config.yml)
- `--batch-id` - ID of the batch to process
- `--verbose` / `-v` - Enable verbose output
- `--log-file` - Specify log file location
- `--dry-run` / `--apply` - Preview changes or actually apply them

## Domain Summary XML Format

The domain summary files use a consistent XML format with:
- Root element: `<blast_summ_doc>`
- Main container: `<blast_summ pdb="XXXX" chain="X">`
- BLAST results sections:
  - `<chain_blast_run>` - Chain-level BLAST results
  - `<blast_run>` - Domain-level BLAST results
- Hit information with query regions in `<query_reg>` elements
- Domain definitions (when available)

## Database Integration

Scripts interact with a PostgreSQL database with these key tables:
- `ecod_schema.batch` - Batch information
- `ecod_schema.protein` - Protein information
- `ecod_schema.process_status` - Processing status tracking
- `ecod_schema.process_file` - File location tracking

## Getting Started

1. Set up the proper environment:
   ```bash
   conda env create -f environment.yml
   conda activate pyecod
   ```

2. Configure your database connection in `config/config.yml`

3. Run the appropriate scripts based on your needs:
   ```bash
   # Generate domain summaries for a batch
   python core/generate_batch_domain_summaries.py --config config/config.yml --batch-id 123
   
   # Analyze domain summaries
   python analysis/analyze_domain_summaries.py --config config/config.yml --batch-id 123
   
   # Fix domain summaries
   python utilities/fix_domain_summaries.py --config config/config.yml --batch-id 123 --apply
   ```

## Contributing

When adding new scripts, please follow these guidelines:
1. Use the ApplicationContext pattern for configuration
2. Follow the established command-line interface pattern
3. Implement proper logging
4. Add appropriate unit tests
5. Document the script in the relevant README.md