# Mini Production - Emergency Scale-Up for pyECOD

This directory contains the production wrapper to scale the proven mini_pyecod algorithm (6/6 regression tests) to process ~40k representative proteins within 2-3 weeks.

## Quick Start

### 1. Install Dependencies
```bash
pip install -r mini_production/requirements.txt
```

### 2. Configure Database Access
```bash
# Copy template and edit with your credentials
cp mini_production/config.template.yml config.local.yml
# Edit config.local.yml with your actual database credentials

# NEVER commit config.local.yml - it contains passwords!
```

Example config.local.yml:
```yaml
database:
  host: dione
  port: 45000
  database: ecod_protein
  user: ecod
  password: "your_password_here"
```

### 3. Validate Setup
```bash
python mini_production/validate_setup.py
```
This checks:
- Config file validity
- Mini executable location
- Batch directory access
- SLURM availability
- Database connectivity
- Test protein discovery

### 4. Test Run (Start Here!)
```bash
# Test with 3 proteins
python mini_production/filesystem_batch_processor.py --test-proteins 3

# Monitor progress
python mini_production/filesystem_batch_processor.py --monitor
```

### 5. Production Scale-Up
```bash
# Process all representative proteins
python mini_production/filesystem_batch_processor.py --scan-all --reps-only --max-jobs 50

# Process specific batch
python mini_production/filesystem_batch_processor.py --batch-name batch_036_20250406_1424 --reps-only

# Check status anytime
python mini_production/filesystem_batch_processor.py --status
```

## How It Works

### Filesystem-First Approach
1. **Scan**: `/data/ecod/pdb_updates/batches/*/domains/` for `*.develop291.domain_summary.xml`
2. **Check**: If corresponding `mini_domains/*.mini.domains.xml` exists
3. **Filter**: Optionally filter to representative proteins via database query
4. **Submit**: Create SLURM jobs for missing results
5. **Track**: SQLite database tracks job progress independently

### File Structure Created
```
/data/ecod/pdb_updates/batches/{batch_name}/
├── domains/                    # Input (existing)
│   └── 8ovp_A.develop291.domain_summary.xml
├── mini_domains/              # Output (created)
│   └── 8ovp_A.mini.domains.xml
└── ...

/tmp/mini_production_*          # Temporary files
├── jobs/                      # SLURM job scripts
├── logs/                      # Job output logs
├── results/                   # Backup results
└── status.db                  # SQLite tracking
```

### Database Integration
- **Read**: Query `ecod_schema.process_status` for representative proteins
- **Track**: Independent SQLite database for mini-specific progress
- **Write**: Later import results to `pdb_analysis.partition_*` tables

## Key Features

### Representative Protein Filtering
```bash
# Database-driven filtering (recommended)
--reps-only

# Queries:
SELECT DISTINCT p.source_id
FROM ecod_schema.process_status ps
JOIN ecod_schema.protein p ON ps.protein_id = p.id
WHERE ps.is_representative = TRUE
```

### Duplicate Handling
- Processes proteins in each batch separately
- Same protein in multiple batches = multiple mini results
- Deduplication happens later during final database import

### SLURM Integration
- Automatic job script generation
- Concurrency control (`--max-jobs`)
- Resource allocation (memory, time, CPUs)
- Progress monitoring via `squeue`/`sacct`

### Error Handling
- Robust file existence checking
- SLURM submission error handling
- Job status tracking and retry capability
- Graceful degradation if database unavailable

## Commands Reference

### Scanning and Processing
```bash
# Scan all batches, show what would be processed
python filesystem_batch_processor.py --scan-all --reps-only

# Process specific batch
python filesystem_batch_processor.py --batch-name batch_036_20250406_1424

# Test mode with known proteins
python filesystem_batch_processor.py --test-proteins 5
```

## Monitoring and Import

### Real-Time Progress Monitoring
```bash
# Check current status
python mini_production/monitor_progress.py

# Continuous monitoring
python mini_production/monitor_progress.py --watch

# Show results ready for import
python mini_production/monitor_progress.py --ready-import

# Detailed batch statistics
python mini_production/monitor_progress.py --stats
```

### Import Results to Production Database

**⚠️ IMPORTANT: Handle Data Collisions Safely**

The production database may already contain results from the main pyECOD system. Always check for collisions first:

```bash
# 1. Detect potential collisions
python mini_production/detect_collisions.py

# 2. Check specific batch
python mini_production/detect_collisions.py --batch-name ecod_batch_031_20250406_1424

# 3. Summary only
python mini_production/detect_collisions.py --summary
```

**Collision Strategies:**
- **`separate`** ✅ (Recommended): Create separate records with `process_version='mini_pyecod_1.0'`
- **`skip`**: Skip proteins that already have records (preserves existing data)
- **`update`**: Overwrite existing records (⚠️ dangerous!)
- **`check`**: Check for collisions without importing

**Safe Import Commands:**
```bash
# Test collision check without importing
python mini_production/import_results.py --check-collisions --limit 10

# Safe import with separate records (recommended)
python mini_production/import_results.py --import-all --collision-strategy separate --limit 100

# Skip existing records
python mini_production/import_results.py --import-all --collision-strategy skip --limit 100

# Import specific batch safely
python mini_production/import_results.py --batch-name ecod_batch_031_20250406_1424 --collision-strategy separate

# Verify imported data quality
python mini_production/import_results.py --verify
```

**Query Mini Results in Database:**
```sql
-- Find mini results
SELECT * FROM pdb_analysis.partition_proteins
WHERE process_version = 'mini_pyecod_1.0';

-- Compare mini vs main results
SELECT
    pdb_id, chain_id,
    process_version,
    is_classified,
    domains_with_evidence,
    timestamp
FROM pdb_analysis.partition_proteins
WHERE pdb_id = '8ovp' AND chain_id = 'A'
ORDER BY timestamp DESC;
```

### Complete Safe Workflow
```bash
# 1. Monitor progress
python mini_production/monitor_progress.py

# 2. Check for collisions BEFORE importing
python mini_production/detect_collisions.py --summary

# 3. Check what's ready for import
python mini_production/monitor_progress.py --ready-import

# 4. Test collision detection without importing
python mini_production/import_results.py --check-collisions --limit 10

# 5. Safe test import with small batch
python mini_production/import_results.py --import-all --collision-strategy separate --limit 50

# 6. Verify quality
python mini_production/import_results.py --verify

# 7. Scale up import safely
python mini_production/import_results.py --import-all --collision-strategy separate --limit 1000
```

### Configuration
```bash
# Validate setup
python validate_setup.py

# Check config
cat mini_production/config.local.yml
```

## Troubleshooting

### Common Issues

**1. "Mini executable not found"**
```bash
# Check location
ls -la mini/pyecod_mini
ls -la ./pyecod_mini

# Update config.local.yml paths.mini_executable
```

**2. "Batch directories not found"** 
```bash
# Check permissions
ls -la /data/ecod/pdb_updates/batches/

# Verify batch structure
ls /data/ecod/pdb_updates/batches/*/domains/ | head
```

**3. "Database connection failed"**
```bash
# Test manually
psql -h dione -p 5432 -U ecod -d ecod_protein

# Check config.local.yml database section
```

**4. "SLURM submission failed"**
```bash
# Check SLURM
sinfo
sbatch --version

# Check partition access
sinfo -p All
```

### Recovery

**Reset tracking database:**
```bash
rm /tmp/mini_production_status.db
# Re-run to recreate
```

**Clean temporary files:**
```bash
rm -rf /tmp/mini_production_*
```

**Check individual job:**
```bash
# Find log files
ls /tmp/mini_production_logs/mini_8ovp_A_*.{out,err}

# Check job status  
squeue -j <job_id>
sacct -j <job_id>
```

## Expected Timeline

### Week 1: Validation (Current)
- ✅ Setup validation and test runs
- ✅ Process 100-500 test proteins
- ✅ Verify filesystem approach works
- ✅ Database integration functional

### Week 2: Scale-Up
- 🎯 Process 1,000+ representative proteins
- 🎯 Achieve <5% failure rate
- 🎯 Stakeholder demo with progress dashboard
- 🎯 Quality validation of results

### Week 3: Research Value
- 🎯 Process 10,000+ representative proteins
- 🎯 Research analysis and comparison
- 🎯 Scientific insights documentation
- 🎯 Project continuation secured

## Architecture Benefits

1. **Independence**: Works without modifying main pyECOD system
2. **Proven Algorithm**: 6/6 regression tests, 80% domain accuracy
3. **Filesystem-First**: Minimal database dependencies
4. **Scalable**: SLURM cluster processing
5. **Trackable**: Real-time progress monitoring
6. **Recoverable**: Robust error handling and retry
7. **Fast**: Immediate results vs months of main system debugging

This approach saves the project by generating research value immediately while maintaining the proven mini algorithm quality.
