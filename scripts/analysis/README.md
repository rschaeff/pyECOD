# pyECOD Analysis Scripts

This directory contains scripts for analyzing domains, batches, and quality metrics in the pyECOD pipeline.

## Scripts

### analyze_domain_summaries.py

Analyzes domain summary content across batches:
- Examines XML structure and content
- Provides statistics on BLAST hits and domains
- Samples representative files for detailed inspection
- Creates comprehensive reports on domain distribution

**Usage:**
```bash
python analyze_domain_summaries.py --config config/config.yml [--batch-id 123] [--sample-size 10] [--detailed]
```

### analyze_batch_domain_summaries.py

Analyzes success and failure patterns in domain summaries:
- Identifies potential issues and failure patterns
- Examines domain classifications (T-groups, H-groups)
- Checks file type discrepancies
- Provides detailed statistics on domain summary quality

**Usage:**
```bash
python analyze_batch_domain_summaries.py --config config/config.yml --batch-id 123 [--sample-size 50]
```

### analyze_blast_summaries.py

Analyzes XML structure of BLAST summary files:
- Validates the expected `blast_summ_doc` structure
- Identifies incorrect or malformed XML structures
- Counts occurrences of different root elements
- Lists examples of invalid files for investigation

**Usage:**
```bash
python analyze_blast_summaries.py --config config/config.yml --batch-id 123 [--max-examples 10]
```

### analyze_batch_dates.py

Analyzes PDB structure deposition dates within batches:
- Queries the `pdb_analysis` schema for dates
- Categorizes structures by date ranges
- Generates recommendations for testing based on recency
- Creates CSV reports of date distributions

**Usage:**
```bash
python analyze_batch_dates.py --config config/config.yml [--batch-ids 123 456] [--min-date 2023-01-01]
```

### analyze_batch_pdb_ids.py

Analyzes PDB IDs to estimate structure age:
- Uses PDB ID prefixes as proxies for structure age
- Identifies batches with recent structures
- Recommends candidate proteins for testing
- Generates CSV reports of PDB ID distributions

**Usage:**
```bash
python analyze_batch_pdb_ids.py --config config/config.yml [--batch-ids 123 456] [--min-recent 50]
```

## Common Parameters

- `--config` - Path to configuration file
- `--batch-id` / `--batch-ids` - ID(s) of the batch(es) to analyze
- `--sample-size` - Number of files to sample for detailed analysis
- `--output` - Output file for detailed results (usually JSON)
- `--verbose` / `-v` - Enable verbose output

## Integration

These analysis scripts provide critical insights for:
1. Quality assurance of domain summaries
2. Identification of processing bottlenecks
3. Testing and validation of the pipeline
4. Selection of representative test cases