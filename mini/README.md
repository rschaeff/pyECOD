# mini_pyecod - Minimal Domain Partitioning Implementation

A clean-room implementation of protein domain partitioning that focuses on the core algorithm without infrastructure dependencies.

## Quick Start

```bash
# Test a single protein
python quick_test.py 8ovp_A

# Find good test cases from your data
python find_real_test_cases.py

# Run batch tests with report
python run_batch_tests.py --auto-find 10
```

## Core Algorithm

The partitioner uses **residue blocking** inspired by the Perl implementation:

1. **Domains are exclusive** - each residue belongs to ONE domain
2. **Evidence processed by quality** - best evidence first (confidence, e-value)
3. **Coverage thresholds**:
   - `NEW_COVERAGE > 70%` - Hit must cover >70% unused residues
   - `OLD_COVERAGE < 10%` - Hit must overlap <10% with existing domains
4. **Residues are blocked** - once assigned, they can't be reassigned

## Key Files

- `partitioner.py` - Core domain partitioning algorithm
- `parser.py` - Parse domain summary XML files
- `models.py` - Minimal data models (Evidence, Domain)
- `writer.py` - Write domain partition XML files
- `quick_test.py` - Test single proteins
- `run_batch_tests.py` - Test multiple proteins with reporting
- `find_real_test_cases.py` - Find test cases from batch data

## Testing Workflow

### 1. Find Test Cases
```bash
python find_real_test_cases.py
```

This analyzes your batch data and recommends proteins in different categories:
- Single domain (baseline)
- Two domains (decomposition test)
- Three domains (medium complexity)
- Many domains (high complexity)
- Discontinuous domains

### 2. Test Individual Proteins
```bash
python quick_test.py 8opd_Aa
```

Shows detailed output including:
- Evidence statistics
- Domain partitioning decisions
- Final domains with ranges

### 3. Run Batch Tests
```bash
# Test recommended proteins
./run_test_suite.sh

# Or use batch runner
python run_batch_tests.py --auto-find 20
```

Generates:
- Human-readable report (`test_report_*.txt`)
- JSON results (`test_results_*.json`)
- Individual domain XML files

## Understanding Output

### Partitioner Output
```
✓ DOMAIN 1: 2ia4 @ 5-256
  Source: hhsearch, Confidence: 0.99
  Coverage: 95.0% new, 0.0% overlap
  Residues assigned: 252, Remaining: 266

Rejection summary (47 evidence items rejected):
  Insufficient New Coverage: 35
  Too Much Overlap: 10
  Too Small: 2

Final coverage: 485/518 residues (93.6%)
```

### Domain XML Output
```xml
<domain_partition pdb_id="8ovp" chain_id="A" reference="mini_pyecod" is_classified="true">
  <domains>
    <domain id="d1" range="5-256" family="2ia4" source="hhsearch" evidence_count="1" is_discontinuous="false" size="252"/>
    <domain id="d2" range="262-368" family="1ytf" source="domain_blast" evidence_count="1" is_discontinuous="false" size="107"/>
    <domain id="d3" range="376-514" family="3mfh" source="hhsearch" evidence_count="1" is_discontinuous="false" size="139"/>
  </domains>
</domain_partition>
```

## Configuration

Key thresholds in `partitioner.py`:

```python
NEW_COVERAGE_THRESHOLD = 0.7   # Hit must be >70% new residues
OLD_COVERAGE_THRESHOLD = 0.1   # Hit must have <10% overlap
MIN_DOMAIN_SIZE = 20          # Minimum domain size
```

## Architecture Rules

To maintain clean architecture:

1. **NO imports from pipelines/** - This is standalone
2. **NO database access** - Works with files only  
3. **NO external dependencies** beyond stdlib + lxml
4. **ONLY import from ecod.core** - Shared utilities like SequenceRange
5. **Must work standalone** - No infrastructure required

## Troubleshooting

**No domains found:**
- Check evidence quality thresholds
- Verify sequence length estimation
- Look at rejection statistics

**Too many small domains:**
- Increase `MIN_DOMAIN_SIZE`
- Check evidence quality filtering

**Poor coverage:**
- Lower confidence thresholds for evidence filtering
- Check if evidence covers the full sequence

**Test failures:**
- Verify domain summary XML exists
- Check for valid evidence in XML
- Ensure proper batch directory path
