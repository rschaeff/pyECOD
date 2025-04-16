# pyECOD Inspection Scripts

This directory contains scripts for examining individual domain summary files and debugging specific issues in the pyECOD pipeline.

## Scripts

### inspect_domain_summary_xml.py

Inspects individual XML domain summary files:
- Displays basic file structure and content
- Shows BLAST hit information
- Lists domains and their boundaries
- Validates XML structure

**Usage:**
```bash
python inspect_domain_summary_xml.py /path/to/domain_summary.xml [--verbose]
```

### debug_domain_summary.py

Comprehensive debug tool for domain summary processing issues:
- Analyzes domain boundary determination
- Provides detailed hit quality assessment
- Tests multiple threshold combinations
- Suggests optimal domain partitioning
- Exports detailed JSON reports for further analysis

**Usage:**
```bash
python debug_domain_summary.py --file /path/to/domain_summary.xml [--output results.json] [--test-thresholds]
```

### check_domain_summaries_xml.py

Checks domain summary XML files across batches:
- Validates XML structure for multiple files
- Provides batch-level statistics on XML quality
- Shows common XML elements and their frequencies
- Identifies patterns in XML structure issues

**Usage:**
```bash
python check_domain_summaries_xml.py --config config/config.yml --batch-id 123 [--sample-size 5] [--all-batches]
```

## Features

These inspection tools provide several key capabilities:

1. **Individual File Debugging**
   - Examine specific problematic files in detail
   - Validate XML structure and content
   - Test domain boundary algorithms

2. **Domain Boundary Analysis**
   - Debug domain boundary determination logic
   - Test different threshold parameters
   - Examine evidence for domain boundaries

3. **BLAST Hit Exploration**
   - Inspect raw BLAST hit data
   - Verify query region parsing
   - Examine E-values and coverage

4. **Data Visualization**
   - Output formatted reports
   - Generate JSON for further analysis
   - Provide human-readable summaries

## Integration

These inspection tools are particularly useful for:
1. Debugging specific problematic proteins
2. Validating changes to the domain summary algorithm
3. Understanding how BLAST hits influence domain boundaries
4. Training new developers on the domain summary format