# pyecod_mini Setup Instructions

## Installation

1. **Save the main script** as `mini/pyecod_mini.py`
2. **Save the wrapper script** as `mini/pyecod_mini` (no extension)
3. **Make it executable**:
   ```bash
   chmod +x mini/pyecod_mini
   ```

4. **Add to PATH** (optional):
   ```bash
   # Add this to your ~/.bashrc or ~/.zshrc
   export PATH="$PATH:/path/to/your/pyecod/mini"
   ```

## File Structure

Your mini directory should look like:
```
mini/
├── pyecod_mini              # Executable wrapper
├── pyecod_mini.py           # Main entry point
├── __init__.py              # Module init
├── models.py                # Data models
├── parser.py                # XML parsing
├── partitioner.py           # Core algorithm
├── decomposer.py            # Chain BLAST decomposition
├── blast_parser.py          # BLAST XML parsing
├── writer.py                # Output writing
├── test_data/               # Reference data
│   ├── domain_definitions.csv
│   ├── domain_lengths.csv
│   └── protein_lengths.csv
└── scripts/                 # Legacy scripts
    └── quick_test.py
```

## Basic Usage

```bash
# Simple usage - uses all intelligent defaults
./pyecod_mini 8ovp_A

# With verbose output
./pyecod_mini 8ovp_A --verbose

# Use specific batch
./pyecod_mini 8ovp_A --batch-id 036

# Show available batches
./pyecod_mini --list-batches

# Validate configuration
./pyecod_mini --validate
```

## Expected Output

```
Found 180 evidence items:
  chain_blast: 6
  domain_blast: 8
  hhsearch: 166

Decomposition status:
  Chain BLAST evidence: 6
  With alignment data: 6
  Domain definitions: ✓

Partitioning domains...

==================================================
RESULTS: 3 domains found
==================================================

1. Domain d1:
   Family: e2ia4A1
   Range: 2-248
   Size: 247 residues
   Source: chain_blast_decomposed

2. Domain d2:
   Family: e2ia4A2
   Range: 491-517
   Size: 27 residues
   Source: chain_blast_decomposed

3. Domain d3:
   Family: 6dgv
   Range: 253-499
   Size: 247 residues
   Source: chain_blast

Total coverage: 521/569 residues (91.6%)

✓ Output written to: /tmp/8ovp_A_mini.domains.xml
```

## Configuration

The tool automatically finds:

- **Batch directories**: `/data/ecod/pdb_updates/batches/ecod_batch_*`
- **Domain summaries**: `{batch}/domains/{protein_id}.develop291.domain_summary.xml`
- **BLAST XML**: `{batch}/blast/chain/{protein_id}.develop291.xml`
- **Reference data**: `mini/test_data/*.csv`

## Troubleshooting

### "No batch directories found"
Check that `/data/ecod/pdb_updates/batches/` exists and contains `ecod_batch_*` directories.

### "Domain summary file not found"
Verify the protein ID format (e.g., `8ovp_A`) and that the batch contains domain files.

### "No evidence with reference lengths found"
Update your `test_data/*.csv` files with more complete reference data.

### "Domain definitions not found"
Chain BLAST decomposition requires domain definitions. Check `test_data/domain_definitions.csv`.

## Migration from quick_test.py

Replace:
```bash
python scripts/quick_test.py 8ovp_A \
  --domain-lengths test_data/domain_lengths.csv \
  --protein-lengths test_data/protein_lengths.csv \
  --domain-definitions test_data/domain_definitions.csv \
  --blast-dir /data/ecod/.../blast/chain \
  --verbose
```

With:
```bash
./pyecod_mini 8ovp_A --verbose
```

All the file paths are discovered automatically!
