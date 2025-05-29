# Domain Analysis Pipeline

## Overview

The Domain Analysis Pipeline is a comprehensive system for identifying, classifying, and analyzing protein domains in the ECOD (Evolutionary Classification of Protein Domains) database. It processes evidence from multiple sources (BLAST, HHSearch, self-comparison) to determine domain boundaries and assign ECOD classifications using a modern service-based architecture.

## Architecture

```
domain_analysis/
├── pipeline.py              # Main orchestrator
├── hhresult_registrar.py    # HHSearch result registration
├── routing.py               # Batch processing router
├── partition_legacy.py      # DEPRECATED - Legacy partition code
├── summary/                 # Evidence collection service
│   ├── service.py           # High-level summary API
│   ├── generator.py         # Core orchestration
│   ├── file_locator.py      # File discovery
│   ├── filters.py           # Quality filtering
│   ├── models.py            # Data structures
│   └── processors/          # Evidence extractors
│       ├── base.py          # Abstract interface
│       ├── blast.py         # BLAST processor
│       ├── hhsearch.py      # HHSearch processor
│       └── self_comparison.py # Self-comparison processor
└── partition/               # Domain boundary identification service
    ├── service.py           # High-level partition API
    ├── analyzer.py          # Evidence analysis
    ├── processor.py         # Boundary processing
    ├── tracker.py           # Status tracking
    └── models.py            # Data structures
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

### 2. Domain Summary Service (`summary/service.py`)

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

### 3. Domain Partition Service (`partition/service.py`)

High-level service for identifying domain boundaries and assigning classifications:

```python
from ecod.pipelines.domain_analysis.partition import DomainPartitionService

# Create service
service = DomainPartitionService(context)

# Process single protein
result = service.partition_protein(
    pdb_id="1abc",
    chain_id="A",
    summary_path="/path/to/summary.xml",
    output_dir="/data/output",
    process_id=123
)

# Process batch
batch_results = service.partition_batch(
    batch_id=123,
    batch_path="/data/batch/123",
    limit=100,
    representatives_only=False
)
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
2. HHSearch Registration (hhresult_registrar.py)
   - Convert HHR files to XML
   - Register in database
   ↓
3. Evidence Collection (summary/service.py)
   - BLAST results (chain & domain)
   - HHSearch results
   - Self-comparison (DALI/HHrepid)
   ↓
4. Evidence Processing (summary/)
   - Quality filtering
   - Evidence validation
   - Confidence calculation
   ↓
5. Domain Partition (partition/service.py)
   - Boundary identification
   - Overlap resolution
   - Classification assignment
   ↓
6. Output: DomainPartitionResult
   - Domain boundaries (DomainModel objects)
   - ECOD classifications
   - Coverage statistics
   - Processing metadata
```

## Data Models

The pipeline uses modern, unified data models:

### Evidence Model
```python
from ecod.models.pipeline.evidence import Evidence

# Unified evidence from any source
evidence = Evidence(
    type="hhsearch",           # blast, hhsearch, self_comparison
    source_id="e4hluA1",       # Source identifier
    domain_id="e4hluA1",       # Domain identifier
    query_range="5-256",       # Query coordinates
    hit_range="1-252",         # Hit coordinates
    confidence=0.999,          # Confidence score
    evalue=1e-50,             # E-value
    t_group="2004.1.1",       # ECOD classification
    h_group="2004.1",
    x_group="2004",
    a_group="a.17"
)
```

### Domain Model
```python
from ecod.models.pipeline.domain import DomainModel

# Complete domain representation
domain = DomainModel(
    id="1abc_A_d1",
    start=5,
    end=256,
    range="5-256",
    source="hhsearch",
    confidence=0.999,
    t_group="2004.1.1",
    h_group="2004.1",
    x_group="2004",
    a_group="a.17",
    evidence=[evidence1, evidence2]  # Supporting evidence
)
```

### Partition Result
```python
from ecod.models.pipeline.partition import DomainPartitionResult

# Complete processing result
result = DomainPartitionResult(
    pdb_id="1abc",
    chain_id="A",
    reference="develop291",
    domains=[domain1, domain2],      # List of DomainModel objects
    success=True,
    is_classified=True,
    coverage=0.95,
    processing_time=2.34,
    domain_file="/path/to/output.xml"
)
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

  # Summary service
  summary_service:
    max_batch_workers: 4
    use_multiprocessing: false
    save_summaries: true

  # Partition service
  partition_service:
    max_workers: 4
    use_multiprocessing: false
    track_status: true

  # HHSearch registration
  hhsearch:
    force_regenerate: false
    validate_files: true

  # Routing
  routing:
    confidence_threshold: 0.9
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

### Example 2: Service-Based Processing

```python
from ecod.pipelines.domain_analysis.summary import DomainSummaryService
from ecod.pipelines.domain_analysis.partition import DomainPartitionService

# Create services
summary_service = DomainSummaryService(context)
partition_service = DomainPartitionService(context)

# Process summary
summary_result = summary_service.process_protein(
    pdb_id="1abc",
    chain_id="A",
    job_dump_dir="/data/jobs/123"
)

# Process partition using summary
partition_result = partition_service.partition_protein(
    pdb_id="1abc",
    chain_id="A",
    summary_path="/data/jobs/123/domains/1abc_A.develop291.evidence_summary.xml",
    output_dir="/data/jobs/123"
)

# Access results
print(f"Found {len(partition_result.domains)} domains")
for domain in partition_result.domains:
    print(f"  Domain: {domain.range}, Classification: {domain.get_classification_level()}")
```

### Example 3: HHSearch Registration

```python
from ecod.pipelines.domain_analysis import HHResultRegistrar

# Register HHSearch results before processing
registrar = HHResultRegistrar(context)

# Register batch results
registered_count = registrar.register_batch_results(
    batch_id=123,
    force_regenerate=False
)

print(f"Registered {registered_count} HHSearch result files")
```

### Example 4: Single Protein Analysis

```python
# Complete analysis for single protein
result = pipeline.analyze_domain(
    pdb_id="1abc",
    chain_id="A",
    output_dir="/data/output",
    reference="develop291",
    blast_only=False
)

# Access results
print(f"Status: {result['status']}")
print(f"Summary file: {result['files'].get('summary')}")
print(f"Partition file: {result['files'].get('partition')}")

# Access models
if 'summary' in result['models']:
    summary_model = result['models']['summary']
    print(f"Evidence sources: {len(summary_model.evidence_sources)}")

if 'domains' in result['models']:
    for domain in result['models']['domains']:
        print(f"Domain {domain.id}: {domain.range}")
```

## Database Schema

The pipeline interacts with several database tables:

```sql
-- Process tracking
ecod_schema.batch               -- Batch information
ecod_schema.process_status      -- Per-protein processing status
ecod_schema.process_file        -- File registration
ecod_schema.protein            -- Protein information

-- Processing paths
ecod_schema.protein_processing_path  -- Routing decisions
ecod_schema.blast_confidence_metrics -- BLAST confidence analysis

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
{pdb_id}_{chain_id}.{reference}.hhsearch.xml       # HHSearch XML
{pdb_id}_{chain_id}.{reference}.hhr                # HHSearch HHR

# Output files
{pdb_id}_{chain_id}.{reference}.evidence_summary.xml    # Evidence summary
{pdb_id}_{chain_id}.{reference}.domains.xml             # Domain partition
```

## Error Handling

The pipeline implements comprehensive error handling:

```python
from ecod.exceptions import PipelineError, ValidationError

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

### From Legacy Code

**Old approach (DEPRECATED):**
```python
# DEPRECATED - Don't use
from ecod.pipelines.domain_analysis.partition_legacy import DomainPartition
partition = DomainPartition(context)
result = partition.process_protein_domains(...)
```

**New approach:**
```python
# Use service-based architecture
from ecod.pipelines.domain_analysis.partition import DomainPartitionService
service = DomainPartitionService(context)
result = service.partition_protein(...)
```

### Key Changes

1. **Service Architecture**: Replaced monolithic classes with focused services
2. **Unified Models**: All evidence uses the `Evidence` class, domains use `DomainModel`
3. **Better Error Handling**: Results include success status and detailed error messages
4. **Enhanced Configuration**: Fine-grained control over processing options
5. **Status Tracking**: Comprehensive tracking of processing status in database
6. **HHSearch Integration**: Dedicated registrar for HHSearch result processing

## Performance Considerations

- **Service Caching**: Evidence and classification data is cached between requests
- **Batch Processing**: Use `process_batch()` methods for multiple proteins
- **Multiprocessing**: Enable for CPU-intensive batch processing:
  ```python
  service_config = {
      'use_multiprocessing': True,
      'max_batch_workers': 4
  }
  ```
- **Database Optimization**: Classification lookups are cached to reduce database load
- **File Validation**: Files are validated before processing to avoid unnecessary work

## Monitoring and Logging

The pipeline provides detailed logging at multiple levels:

```python
import logging

# Enable debug logging
logging.getLogger("ecod.pipelines.domain_analysis").setLevel(logging.DEBUG)

# Monitor specific services
logging.getLogger("ecod.pipelines.domain_analysis.summary").setLevel(logging.INFO)
logging.getLogger("ecod.pipelines.domain_analysis.partition").setLevel(logging.DEBUG)
```

## API Reference

- [DomainSummaryService API](summary/README.md)
- [DomainPartitionService API](partition/service.py)
- [HHResultRegistrar API](hhresult_registrar.py)
- [ProcessingRouter API](routing.py)

## Contributing

When extending the pipeline:

1. **Use Services**: Extend existing services rather than creating new monolithic classes
2. **Follow Models**: Use `Evidence` and `DomainModel` consistently
3. **Add Tests**: Include unit tests for new components
4. **Update Documentation**: Keep service APIs documented
5. **Handle Errors**: Use appropriate exception types and result objects

## Troubleshooting

### Common Issues

1. **Missing HHSearch Files**: Run `HHResultRegistrar.register_batch_results()` first
2. **Service Initialization**: Ensure context is properly configured
3. **Memory Issues**: Reduce batch size or disable multiprocessing
4. **Classification Failures**: Verify reference database is populated

### Debug Mode

```python
# Enable detailed statistics and logging
summary_service = DomainSummaryService(context, {
    'generator': {'collect_detailed_stats': True}
})

partition_service = DomainPartitionService(context, {
    'track_status': True,
    'save_intermediate': True
})

# Get detailed statistics
stats = summary_service.get_service_statistics()
print(stats)
```

## License

Part of the ECOD project. See main LICENSE file.
