# pyECOD Analysis Scripts

This directory contains scripts for analyzing domains, batches, and quality metrics in the pyECOD pipeline.

## Batch Analysis Scripts

### batch_blast_analysis.py

Analyzes batch proteins for BLAST-only partition suitability:
- Examines domain summary XML files for BLAST hit patterns
- Generates statistics on peptides and chains with/without hits
- Creates detailed reports and visualizations
- Recommends processing strategies based on hit coverage

**Usage:**
```bash
python batch_blast_analysis.py --config config/config.yml --batch-id 123 --output-dir results/batch123
```

### batch_domain_analysis.py

Analyzes domain generation results for specific batches:
- Examines quality of domain generation
- Identifies proteins with no domains and their BLAST hit status
- Breaks down chain vs. domain BLAST usage
- Identifies candidates for full HHsearch pipelines

**Usage:**
```bash
python batch_domain_analysis.py --config config/config.yml --batch-id 123 --output-dir analysis_results
```

### batch_domain_analyzer.py

Advanced analyzer for determining partition strategy:
- Performs detailed analysis of BLAST hit confidence
- Evaluates HHsearch candidacy with sophisticated criteria
- Generates actionable recommendations for processing
- Creates processing scripts for execution

**Usage:**
```bash
python batch_domain_analyzer.py --config config/config.yml --batch-id 123 --output-dir batch_123_analysis
```

### batch_domain_partition_analyzer.py

Focused analyzer for domain partition suitability:
- Examines hit quality and distribution
- Generates comprehensive partition strategy reports
- Creates visualizations of partition categories
- Provides recommendations for optimal processing

**Usage:**
```bash
python batch_domain_partition_analyzer.py --config config/config.yml --batch-id 123 --peptide-threshold 30
```

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

### analyze_routing_confidence_scores.py

Analyzes confidence s