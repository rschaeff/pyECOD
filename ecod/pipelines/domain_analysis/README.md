```markdown
# Domain Analysis Pipeline

## Overview

The Domain Analysis Pipeline is a comprehensive system for identifying, classifying, and analyzing protein domains in the ECOD (Evolutionary Classification of Protein Domains) database. It processes evidence from multiple sources (BLAST, HHSearch, self-comparison) to determine domain boundaries and assign ECOD classifications.

## Architecture

```
domain_analysis/
├── pipeline.py              # Main orchestrator
├── partition.py             # Domain boundary identification
├── hhresult_registrar.py    # HHSearch result registration
├── routing.py               # Batch processing router
└── summary/                 # Evidence collection module
    ├── service.py           # High-level API
    ├── generator.py         # Core orchestration
    ├── file_locator.py      # File discovery
    ├── filters.py           # Quality filtering
    ├── models.py            # Data structures
    └── processors/          # Evidence extractors
        ├── base.py          # Abstract interface
        ├── blast.py         # BLAST processor
        ├── hhsearch.py      # HHSearch processor
        └── self_comparison.py # Self-comparison processor
```

## Key Components

### 1. Domain Analysis Pipeline (`pipeline.py`)

The main orchestrator that coordinates the entire domain analysis workflow:

```python
from ecod.pipelines.domain_analysis import DomainAnalysisPipeline
from ecod.core.context import ApplicationContext

# Initialize pipeline
context = ApplicationContext()
pipeline = DomainAnalysisPipeline(context)

# Run complete pipeline for a batch
result = pipeline.run_pipeline(
    batch_id=123,
    blast_only=False,      # Include HHSearch
    limit=None,            # Process all proteins
    partition_only=False,  # Run both summary and partition
    reps_only=False,       # Process all, not just representatives
    reset_failed=True      # Reset previously failed proteins
)
```

### 2. Domain Summary Service (`summary/`)

Modern, modular service for collecting and processing evidence:

```python
from ecod.pipelines.domain_analysis.summary import DomainSummaryService

# Create service
service = DomainSummaryService(context)

# Process single protein
result = service.process_protein(
    pdb_id="1abc",
    chain_id="A",
    job_dump_dir="/data/jobs/123",
    blast_only=False,
    force_overwrite=False
)

# Process batch
proteins = [("1abc", "A"), ("2def", "B"), ("3ghi", "C")]
batch_results = service.process_batch(
    proteins=proteins,
    job_dump_dir="/data/jobs/123"
)
```

### 3. Domain Partition (`partition.py`)

Identifies domain boundaries and assigns classifications using Evidence models:

```python
from ecod.pipelines.domain_analysis.partition import DomainPartition

# Initialize partition module
partition = DomainPartition(context)

# Process domains for a protein
result = partition.process_protein_domains(
    pdb_id="1abc",
    chain_id="A",
    domain_summary_path="/path/to/summary.xml",
    output_dir="/data/output",
    reference="develop291"
)

# Returns DomainPartitionResult with:
# - domains: List[DomainModel] with boundaries and classifications
# - success: bool
# - coverage: float
# - is_classified/is_unclassified/is_peptide: bool
```

### 4. HHSearch Result Registrar (`hhresult_registrar.py`)

Handles HHSearch result file discovery, conversion, and registration:

```python
from ecod.pipelines.domain_analysis.hhresult_registrar import HHResultRegistrar

registrar = HHResultRegistrar(context)

# Register all HHSearch results for a batch
count = registrar.register_batch_results(
    batch_id=123,
    force_regenerate=False  # Use existing XML if available
)

# Register specific proteins
count = registrar.register_specific_chains(
    batch_id=123,
    chain_ids=["1abc_A", "2def_B"],
    force_regenerate=True
)
```

### 5. Processing Router (`routing.py`)

Intelligent routing for batch processing based on BLAST confidence:

```python
from ecod.pipelines.domain_analysis.routing import ProcessingRouter

router = ProcessingRouter(context)

# Assign proteins to processing paths
paths = router.assign_processing_paths(batch_id=123)
# Returns: {
#     "blast_only": [...],     # High-confidence BLAST hits
#     "full_pipeline": [...]   # Need HHSearch analysis
# }

# Prioritize full pipeline proteins
priorities = router.prioritize_full_pipeline_proteins(batch_id=123)
# Returns: {
#     "high_priority": [...],    # Complex/novel domains
#     "standard_priority": [...] # Typical cases
# }
```

## Data Flow

```
1. Input: Protein sequences from batch
   ↓
2. Evidence Collection (summary/)
   - BLAST results (chain & domain)
   - HHSearch results
   - Self-comparison (DALI/HHrepid)
   ↓
3. Evidence Processing
   - Quality filtering
   - Discontinuous domain stitching
   - Confidence calculation
   ↓
4. Domain Partition (partition.py)
   - Boundary identification
   - Overlap resolution
   - Classification assignment
   ↓
5. Output: DomainPartitionResult
   - Domain boundaries
   - ECOD classifications
   - Coverage statistics
```

## Configuration

### Application Configuration (`config.yml`)

```yaml
domain_analysis:
  # Evidence collection
  evidence:
    min_confidence: 0.3
    hsp_evalue_threshold: 0.005
    hit_coverage_threshold: 0.7
  
  # Domain partitioning
  partition:
    high_confidence_threshold: 0.95
    medium_confidence_threshold: 0.7
    overlap_threshold: 0.3
    gap_tolerance: 20
  
  # Processing options
  summary_service:
    max_batch_workers: 4
    use_multiprocessing: false
    save_summaries: true
    
  # Routing
  routing:
    confidence_threshold: 0.9
```

### Processing Options

```python
from ecod.pipelines.domain_analysis.summary import SummaryOptions

options = SummaryOptions(
    blast_only=False,              # Include HHSearch
    force_overwrite=True,          # Overwrite existing files
    skip_filtering=False,          # Apply quality filters
    stitch_discontinuous=True,     # Stitch discontinuous domains
    max_gap_for_stitching=30,      # Max gap between segments
    min_confidence=0.3,            # Minimum evidence confidence
    min_coverage=0.0,              # Minimum coverage requirement
    max_evalue=10.0               # Maximum acceptable e-value
)
```

## Usage Examples

### Example 1: Complete Pipeline for a Batch

```python
from ecod.pipelines.domain_analysis import DomainAnalysisPipeline
from ecod.core.context import ApplicationContext

# Initialize
context = ApplicationContext()
pipeline = DomainAnalysisPipeline(context)

# Process batch
result = pipeline.run_pipeline(
    batch_id=123,
    blast_only=False,
    limit=100,  # Process first 100 proteins
    reps_only=True  # Only representative proteins
)

# Check results
print(f"Success: {result.success}")
print(f"Proteins processed: {result.processing_stats['total_proteins']}")
print(f"Domains identified: {result.partition_stats['domains_created']}")
```

### Example 2: Process Specific Proteins

```python
# Process specific proteins by ID
success = pipeline.process_proteins(
    batch_id=123,
    protein_ids=[1, 2, 3],  # Database protein IDs
    blast_only=False,
    partition_only=False
)
```

### Example 3: Single Protein Analysis

```python
# Analyze a single protein
result = pipeline.analyze_domain(
    pdb_id="1abc",
    chain_id="A",
    output_dir="/data/output",
    reference="develop291",
    blast_only=False
)

# Access results
print(f"Status: {result['status']}")
print(f"Domains found: {len(result['models']['domains'])}")
for domain in result['models']['domains']:
    print(f"  Domain: {domain.range}, Classification: {domain.get_classification_level()}")
```

### Example 4: Partition-Only Mode

```python
# Run only partition step (assumes summaries exist)
result = pipeline.run_pipeline(
    batch_id=123,
    partition_only=True  # Skip summary generation
)
```

## Database Schema

The pipeline interacts with several database tables:

```sql
-- Process tracking
ecod_schema.batch               -- Batch information
ecod_schema.process_status      -- Per-protein processing status
ecod_schema.process_file        -- File registration
ecod_schema.protein            -- Protein information

-- Reference data
pdb_analysis.domain            -- ECOD domain classifications
pdb_analysis.protein           -- Reference protein data
```

## File Naming Conventions

The pipeline expects and generates files following these patterns:

```
# Input files
{pdb_id}_{chain_id}.fasta                          # Sequence
{pdb_id}_{chain_id}.chain_blast.xml                # Chain BLAST
{pdb_id}_{chain_id}.domain_blast.xml               # Domain BLAST
{pdb_id}_{chain_id}.{reference}.hhsearch.xml       # HHSearch
{pdb_id}_{chain_id}.self_comp.xml                  # Self-comparison

# Output files
{pdb_id}_{chain_id}.{reference}.domain_summary.xml     # Evidence summary
{pdb_id}_{chain_id}.{reference}.domains.xml            # Domain partition
```

## Error Handling

The pipeline implements comprehensive error handling:

```python
try:
    result = pipeline.run_pipeline(batch_id=123)
except PipelineError as e:
    # Pipeline-specific errors
    logger.error(f"Pipeline error: {e}")
except ValidationError as e:
    # Input validation errors
    logger.error(f"Validation error: {e}")
except Exception as e:
    # Unexpected errors
    logger.error(f"Unexpected error: {e}")
```

## Migration Guide

### From Old Dictionary-Based Approach

```python
# Old approach (DEPRECATED)
from ecod.pipelines.domain_analysis.summary import DomainSummary
summary = DomainSummary(context)
result = summary.create_summary(pdb_id, chain_id, reference, job_dir)

# New approach
from ecod.pipelines.domain_analysis.summary import DomainSummaryService
service = DomainSummaryService(context)
result = service.process_protein(pdb_id, chain_id, job_dir)
```

### Key Changes

1. **Evidence Model**: All evidence now uses the unified `Evidence` class instead of separate `BlastHit`/`HHSearchHit` models
2. **Service Architecture**: Modular service-based design instead of monolithic classes
3. **Better Error Handling**: Results include success status and detailed error messages
4. **Enhanced Statistics**: Comprehensive tracking of processing metrics
5. **Flexible Configuration**: Fine-grained control over processing options

## Performance Considerations

- **Batch Processing**: Use `process_batch()` for multiple proteins to benefit from connection pooling
- **Multiprocessing**: Enable for CPU-intensive batch processing:
  ```python
  service_config = {
      'use_multiprocessing': True,
      'max_batch_workers': 4
  }
  ```
- **Caching**: File parsing results are cached during batch processing
- **Database Queries**: Classification lookups are cached to reduce database load

## Monitoring and Logging

The pipeline provides detailed logging at multiple levels:

```python
import logging

# Enable debug logging
logging.getLogger("ecod.pipelines.domain_analysis").setLevel(logging.DEBUG)

# Monitor specific components
logging.getLogger("ecod.pipelines.domain_analysis.partition").setLevel(logging.INFO)
logging.getLogger("ecod.pipelines.domain_analysis.summary").setLevel(logging.DEBUG)
```

## Troubleshooting

### Common Issues

1. **Missing Files**: Check file naming conventions and paths
2. **No Domains Found**: Verify evidence quality thresholds
3. **Classification Failures**: Ensure reference database is populated
4. **Memory Issues**: Reduce batch size or disable multiprocessing

### Debug Mode

```python
# Enable detailed statistics
result = pipeline.run_pipeline(
    batch_id=123,
    collect_detailed_stats=True
)

# Access detailed metrics
print(result.processing_stats)
print(result.summary_stats)
print(result.partition_stats)
```

## API Reference

See individual module documentation:
- [DomainSummaryService API](summary/README.md)
- [DomainPartition API](partition.py)
- [ProcessingRouter API](routing.py)

## Contributing

When adding new evidence sources:

1. Create a processor in `summary/processors/`
2. Implement the `EvidenceProcessor` interface
3. Register in `generator.py`
4. Update file discovery in `file_locator.py`

## License

Part of the ECOD project. See main LICENSE file.
```

This README provides a comprehensive overview of the domain analysis pipeline, including:

1. **Architecture overview** with component descriptions
2. **Usage examples** for common scenarios
3. **Configuration options** and file conventions
4. **Migration guide** from the old approach
5. **Performance tips** and troubleshooting
6. **API references** for deeper dives

The documentation emphasizes the new Evidence-based architecture while providing clear migration paths from the legacy dictionary-based approach.
