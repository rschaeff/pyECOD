# Mini PyECOD Test Suite

This directory contains the formal test suite for mini_pyecod domain partitioning.

## Structure

```
tests/
├── __init__.py              # Package initialization
├── conftest.py              # Pytest configuration and fixtures
├── test_cases.py            # Formal test cases (8ovp_A, etc.)
├── test_core.py             # Core algorithm unit tests
├── test_decomposition.py    # Chain BLAST decomposition tests
├── test_visualization.py    # PyMOL visualization tests
└── README.md               # This file
```

## Test Categories

### 1. Formal Test Cases (`test_cases.py`)
These are the "gold standard" tests that validate the algorithm works correctly on known proteins:

- **8ovp_A**: GFP-PBP fusion with chain BLAST decomposition (3 domains expected)
- Additional validated cases as they are added

### 2. Core Algorithm Tests (`test_core.py`)
Unit tests for the partitioning algorithm:

- Evidence parsing and filtering
- Residue blocking logic
- Domain assignment and merging
- Coverage calculations

### 3. Decomposition Tests (`test_decomposition.py`)
Tests for chain BLAST decomposition:

- Alignment parsing from BLAST XML
- Domain definition loading
- Mapping query positions to reference domains
- Discontinuous domain handling

### 4. Visualization Tests (`test_visualization.py`)
Tests for PyMOL comparison generation:

- Domain file parsing (old vs new format)
- PyMOL script generation
- Structure file location

## Running Tests

### All Tests
```bash
cd mini
python -m pytest tests/
```

### Specific Categories
```bash
# Just formal test cases
python -m pytest tests/test_cases.py

# Core algorithm only
python -m pytest tests/test_core.py

# Skip slow tests
python -m pytest tests/ -m "not slow"

# Integration tests only
python -m pytest tests/ -m integration
```

### Specific Test Cases
```bash
# Run just the 8ovp_A test
python tests/test_cases.py --protein 8ovp_A

# Run with verbose output
python tests/test_cases.py --protein 8ovp_A --verbose
```

## Requirements

### Required Files
Tests require reference data to be set up first:

```bash
# Generate reference files from ECOD range cache
python range_cache_parser.py --cache-file /path/to/cache --output-dir test_data
```

Required files:
- `test_data/domain_lengths.csv`
- `test_data/protein_lengths.csv`
- `test_data/domain_definitions.csv`

### Required Data
Tests need access to:
- Batch directory with domain summary XML files
- BLAST XML files (for decomposition tests)
- PDB structure files (for visualization tests)

## Adding New Test Cases

### 1. Validate the Protein
First, manually validate that mini_pyecod works correctly:

```bash
./pyecod_mini NEW_PROTEIN_ID --verbose
```

Check the results and determine expected domain count and architecture.

### 2. Add to Formal Test Cases
Edit `test_cases.py` and add a new `TestCase`:

```python
TestCase(
    protein_id="NEW_PROTEIN_ID",
    description="Description of protein architecture",
    expected_domain_count=N,
    expected_domains=[
        ExpectedDomain(
            family="family_name",
            approximate_range="start-end",
            min_size=min_residues,
            max_size=max_residues,
            discontinuous=False,
            notes="Optional notes"
        ),
        # Add more expected domains...
    ],
    requires_decomposition=False,  # True if chain BLAST decomposition needed
    requires_blast_alignments=False,  # True if BLAST XML needed
    notes="Optional notes about the test case"
)
```

### 3. Test the New Case
```bash
python tests/test_cases.py --protein NEW_PROTEIN_ID --verbose
```

### 4. Document the Case
Add notes about why this test case is important and what it validates.

## Test Fixtures

The test suite provides several fixtures (defined in `conftest.py`):

- `batch_dir`: Default batch directory path
- `test_data_dir`: Directory containing reference CSV files
- `reference_data`: Pre-loaded domain/protein lengths and definitions
- `sample_evidence`: Parsed evidence for 8ovp_A
- `temp_output_dir`: Temporary directory for test outputs

## Continuous Integration

For CI/CD systems, tests can be run with specific configurations:

```bash
# Quick tests only (skip slow decomposition tests)
python -m pytest tests/ -m "not slow"

# Core functionality only (no visualization)
python -m pytest tests/ -m "not visualization"

# Generate JUnit XML for CI
python -m pytest tests/ --junitxml=test-results.xml
```

## Test Data Management

### Updating Reference Data
When ECOD releases new versions:

1. Update the range cache file path in `range_cache_parser.py`
2. Regenerate reference files:
   ```bash
   python range_cache_parser.py --cache-file NEW_CACHE --output-dir test_data
   ```
3. Run tests to ensure compatibility:
   ```bash
   python -m pytest tests/test_cases.py
   ```

### Adding New Reference Proteins
To support tests for new proteins:

1. Ensure the protein exists in the ECOD range cache
2. Regenerate reference files to include the new protein
3. The test framework will automatically use the updated references

## Debugging Failed Tests

### 1. Check Reference Data
```bash
# Verify reference files exist and have expected proteins
head test_data/domain_lengths.csv
grep PROTEIN_ID test_data/protein_lengths.csv
```

### 2. Run with Verbose Output
```bash
python tests/test_cases.py --protein FAILING_PROTEIN --verbose
```

### 3. Check Domain Summary XML
Ensure the protein's domain summary file exists and is parseable:
```bash
ls /path/to/batch/domains/PROTEIN_ID.develop291.domain_summary.xml
python -c "import xml.etree.ElementTree as ET; ET.parse('path/to/file')"
```

### 4. Validate Evidence Parsing
Check that evidence is being parsed correctly:
```bash
# Run the parser directly to see what evidence is found
python -c "
from mini.parser import parse_domain_summary
evidence = parse_domain_summary('path/to/xml', {}, {}, {}, verbose=True)
print(f'Found {len(evidence)} evidence items')
"
```

## Performance Notes

- **Formal test cases** (8ovp_A, etc.): ~10-30 seconds each
- **Core unit tests**: ~1-5 seconds total
- **Decomposition tests**: ~5-15 seconds (marked as slow)
- **Visualization tests**: ~5-10 seconds

Use markers to skip slow tests during development:
```bash
python -m pytest tests/ -m "not slow"
```
