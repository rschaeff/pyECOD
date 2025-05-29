# Domain Partition Module

The Domain Partition Module is a core component of the ECOD (Evolutionary Classification of Protein Domains) pipeline, responsible for identifying domain boundaries from evidence and assigning ECOD classifications. This module has been redesigned with a modern, service-oriented architecture that provides better maintainability, testability, and performance.

## Overview

The partition module takes evidence collected by the [summary module](../summary/) and processes it to:

1. **Identify Domain Boundaries**: Analyze evidence to determine where domains start and end
2. **Validate Reference Coverage**: Ensure evidence adequately covers reference domains
3. **Handle Discontinuous Domains**: Process domains that span multiple sequence segments
4. **Resolve Overlaps**: Manage conflicts between overlapping domain predictions
5. **Assign Classifications**: Apply ECOD hierarchical classifications (T, H, X, A groups)
6. **Track Processing**: Monitor status and performance throughout the pipeline

## Architecture

```
partition/
├── service.py              # High-level API and batch processing
├── processor.py            # Core domain identification algorithms
├── analyzer.py             # Evidence analysis and validation
├── tracker.py              # Status tracking and monitoring
├── reference_analyzer.py   # Reference coverage validation
├── models.py              # Data models and configuration
├── migration.md           # Migration guide from legacy system
└── utils.py               # Utility functions
```

## Core Components

### 1. Service Layer (`service.py`)

**Purpose**: Provides the main API for domain partitioning with both single protein and batch processing capabilities.

**Key Classes**:
- `DomainPartitionService`: Main service interface

**Common Methods**:
```python
# Single protein processing
result = service.partition_protein(
    pdb_id="1abc", 
    chain_id="A",
    summary_path="/path/to/summary.xml",
    output_dir="/path/to/output"
)

# Batch processing
batch_results = service.partition_batch(
    batch_id=123,
    batch_path="/path/to/batch",
    limit=100,
    representatives_only=True
)

# Reprocess failed proteins
retry_results = service.reprocess_failed(batch_id=123, batch_path="/path/to/batch")
```

**Usage Concerns**:
- Requires database connection for classification lookups
- Handles file I/O and can process large batches
- Memory usage scales with batch size when using parallel processing
- Automatic status tracking can generate significant database activity

### 2. Processing Core (`processor.py`)

**Purpose**: Implements the core algorithms for domain boundary identification, overlap resolution, and classification assignment.

**Key Classes**:
- `PartitionProcessor`: Main processing engine
- `DiscontinuousDomainCandidate`: Handles multi-segment domains

**Common Methods**:
```python
# Process evidence into domains
result = processor.process_evidence(evidence_list, context)

# Resolve overlapping domains
resolved = processor.resolve_domain_overlaps(candidates, sequence_length)

# Assign classifications
processor.assign_domain_classifications(domain_models)
```

**Usage Concerns**:
- **Reference Coverage**: Requires database access for optimal performance
- **Memory Usage**: Classification cache can grow large with extensive processing
- **Discontinuous Domains**: Special handling for domains with gaps
- **Protected Domains**: High-confidence domains are preserved during overlap resolution

### 3. Evidence Analysis (`analyzer.py`)

**Purpose**: Comprehensive evidence validation, quality control, and parsing of domain summary files.

**Key Classes**:
- `EvidenceAnalyzer`: Main analysis engine
- `ClassificationCache`: Performance optimization for database lookups
- `EvidenceQualityMetrics`: Detailed quality assessment

**Common Methods**:
```python
# Parse domain summary files
summary_data = analyzer.parse_domain_summary(file_path)

# Extract and validate evidence
evidence_list = analyzer.extract_evidence_with_classification(summary_data)

# Validate individual evidence
validation = analyzer.validate_evidence(evidence, context)

# Generate quality metrics
metrics = analyzer.calculate_quality_metrics(evidence_list, sequence_length, processing_time)
```

**Usage Concerns**:
- **XML Parsing**: Robust error handling for malformed files, but can be memory-intensive
- **Validation Levels**: STRICT mode may reject useful evidence; LENIENT mode more forgiving
- **Caching**: Classification cache improves performance but consumes memory
- **Parallel Processing**: Available but increases complexity and resource usage

### 4. Reference Coverage Analysis (`reference_analyzer.py`)

**Purpose**: Validates that evidence adequately covers reference domains and suggests boundary extensions.

**Key Classes**:
- `ReferenceCoverageAnalyzer`: Coverage validation engine
- `EvidenceWithCoverage`: Enhanced evidence with coverage information

**Common Methods**:
```python
# Analyze evidence coverage
enhanced_evidence = analyzer.analyze_evidence_coverage(evidence)

# Get reference domain information
ref_info = analyzer.get_reference_info(domain_id)

# Suggest extended boundaries
new_start, new_end = analyzer.suggest_extended_boundaries(
    evidence, current_start, current_end, sequence_length
)
```

**Usage Concerns**:
- **Database Dependency**: Requires access to reference domain database
- **Coverage Thresholds**: Conservative settings may reject valid evidence
- **Discontinuous References**: Complex coverage calculations for multi-segment domains
- **Boundary Extension**: Can significantly alter domain boundaries

### 5. Status Tracking (`tracker.py`)

**Purpose**: Comprehensive monitoring of processing status, error handling, and performance metrics.

**Key Classes**:
- `StatusTracker`: Main tracking engine
- `ProcessInfo`: Individual process information
- `PerformanceMetrics`: System performance monitoring

**Common Methods**:
```python
# Track process lifecycle
tracker.start_process(process_id, pdb_id, chain_id, batch_id)
tracker.update_process_stage(process_id, stage, status, progress)
tracker.add_process_error(process_id, error_message, retry=True)

# Monitor performance
metrics = tracker.get_performance_metrics()
error_summary = tracker.get_error_summary(hours=24)

# Batch management
batch_status = tracker.get_batch_status(batch_id)
tracker.update_batch_completion_status(batch_id, representatives_only=False)
```

**Usage Concerns**:
- **Database Resilience**: Continues operation when database is unavailable
- **Memory Growth**: In-memory tracking can accumulate over time
- **Cleanup**: Regular cleanup of old processes needed
- **Performance Impact**: Extensive tracking can slow processing

### 6. Configuration (`models.py`)

**Purpose**: Unified data models and configuration management for all partition components.

**Key Classes**:
- `PartitionOptions`: Comprehensive configuration management
- `ValidationResult`: Validation outcome tracking
- `BatchPartitionResults`: Batch processing results
- `PartitionContext`: Processing context information

**Common Configuration**:
```python
# Create custom options
options = PartitionOptions(
    validation_level=ValidationLevel.STRICT,
    min_reference_coverage=0.8,
    overlap_threshold=0.3,
    merge_gap_tolerance=20,
    use_cache=True,
    parallel_processing=True,
    max_workers=4
)

# Validate configuration
options.validate()  # Throws ValueError if invalid

# Factory methods for common configurations
blast_options = PartitionOptions.create_blast_only()
strict_options = PartitionOptions.create_strict()
```

**Usage Concerns**:
- **Option Validation**: Always call `validate()` after creating custom options
- **Coverage Settings**: Reference coverage thresholds significantly impact results
- **Performance Trade-offs**: Caching and parallel processing use more resources
- **Compatibility**: Enum values must match string representations exactly

## Integration with Summary Module

The partition module is tightly integrated with the summary module:

```python
# Typical workflow
from ecod.pipelines.domain_analysis.summary import DomainSummaryService
from ecod.pipelines.domain_analysis.partition import DomainPartitionService

# 1. Generate evidence summary
summary_service = DomainSummaryService(context)
summary_result = summary_service.process_protein(pdb_id, chain_id, job_dir)

# 2. Partition domains from summary
partition_service = DomainPartitionService(context)
partition_result = partition_service.partition_protein(
    pdb_id, chain_id, summary_result.output_file, output_dir
)
```

**Data Flow**:
1. Summary module collects raw evidence (BLAST, HHSearch, self-comparison)
2. Summary module creates domain_summary.xml file
3. Partition module parses domain_summary.xml
4. Partition module identifies domain boundaries
5. Partition module creates domains.xml output

## Pipeline Integration

The partition module fits into the broader domain analysis pipeline:

```
Input Proteins
      ↓
Evidence Collection (summary module)
  - BLAST searches
  - HHSearch analysis  
  - Self-comparison
      ↓
Domain Summary Creation (summary module)
      ↓
Domain Partitioning (partition module) ← YOU ARE HERE
  - Evidence validation
  - Boundary identification
  - Overlap resolution
  - Classification assignment
      ↓
Domain Classification (downstream)
      ↓
Final Domain Models
```

## Quick Start Examples

### Basic Single Protein Processing

```python
from ecod.pipelines.domain_analysis.partition import partition_single_protein

# Simple case - use defaults
result = partition_single_protein(
    pdb_id="1abc",
    chain_id="A", 
    summary_path="/data/summaries/1abc_A.domain_summary.xml",
    output_dir="/data/output"
)

if result.success:
    print(f"Found {len(result.domains)} domains")
    for domain in result.domains:
        print(f"  {domain.id}: {domain.range} ({domain.get_classification_level()})")
else:
    print(f"Processing failed: {result.error}")
```

### Batch Processing with Custom Options

```python
from ecod.pipelines.domain_analysis.partition import create_service, PartitionOptions

# Create service with custom configuration
service = create_service()
service.default_options = PartitionOptions(
    min_reference_coverage=0.8,  # Require 80% coverage
    overlap_threshold=0.2,       # Allow 20% overlap
    use_cache=True,             # Enable caching
    parallel_processing=True,    # Use multiple cores
    max_workers=8               # Use 8 worker threads
)

# Process batch
batch_results = service.partition_batch(
    batch_id=123,
    batch_path="/data/batch_123",
    limit=1000,                 # Process first 1000 proteins
    representatives_only=True   # Only representative proteins
)

print(f"Batch Results: {batch_results.success_count}/{batch_results.total} succeeded")
print(f"Found {batch_results.total_domains_found} total domains")
```

### Advanced Evidence Analysis

```python
from ecod.pipelines.domain_analysis.partition import EvidenceAnalyzer, PartitionOptions

# Create analyzer with strict validation
options = PartitionOptions(validation_level=ValidationLevel.STRICT)
analyzer = EvidenceAnalyzer(options)

# Analyze domain summary
analysis_result = analyzer.analyze_domain_summary(
    file_path="/data/summaries/complex_protein.xml",
    protein_id="complex_protein",
    sequence_length=500
)

if analysis_result['success']:
    # Access detailed metrics
    quality = analysis_result['quality_metrics']
    print(f"Evidence coverage: {quality['sequence_coverage']:.1%}")
    print(f"Classification consistency: {quality['classification_consistency']:.1%}")
    
    # Review domain suggestions
    for suggestion in analysis_result['domain_suggestions']:
        print(f"Suggested domain: {suggestion['range']} (confidence: {suggestion['confidence']:.2f})")
```

## Performance Considerations

### Memory Usage
- **Evidence Lists**: Scale with number of BLAST/HHSearch hits
- **Classification Cache**: Grows with number of unique domains processed
- **Batch Processing**: Memory usage increases with batch size
- **Cleanup**: Regular cleanup of old processes prevents memory leaks

### Database Load
- **Classification Lookups**: Can generate many database queries
- **Status Updates**: Frequent updates during batch processing
- **Connection Pooling**: Important for batch processing performance
- **Caching**: Significantly reduces database load

### File I/O
- **XML Parsing**: Can be slow for very large summary files
- **Parallel Processing**: Improves throughput but increases file handle usage
- **Network Storage**: Performance depends on storage system characteristics

## Error Handling and Recovery

### Common Error Scenarios

1. **Malformed XML Files**:
   ```python
   # The analyzer includes robust XML repair capabilities
   result = analyzer.parse_domain_summary(problematic_file)
   if 'error' in result:
       print(f"XML parsing failed: {result['error']}")
   ```

2. **Database Connection Issues**:
   ```python
   # Service continues operation without database
   service = DomainPartitionService(context)
   if not service.validate_setup()['database']:
       print("Warning: Database unavailable, using reduced functionality")
   ```

3. **Missing Reference Domains**:
   ```python
   # Coverage analyzer handles missing references gracefully
   enhanced_evidence = coverage_analyzer.analyze_evidence_coverage(evidence)
   if enhanced_evidence.coverage_warning:
       print(f"Coverage warning: {enhanced_evidence.coverage_warning}")
   ```

### Recovery Strategies

- **Retry Failed Processes**: Use `reprocess_failed()` to retry failed proteins
- **Adjust Validation**: Use LENIENT validation for problematic datasets
- **Reduce Coverage Requirements**: Lower thresholds for difficult cases
- **Batch Subdivision**: Process smaller batches to reduce memory pressure

## Migration from Legacy System

If migrating from the old monolithic `DomainPartition` class, see [`migration.md`](migration.md) for detailed guidance.

**Key Changes**:
- Service-oriented architecture replaces monolithic class
- Evidence-based models replace dictionary-based data
- Comprehensive validation with configurable levels
- Reference coverage validation is now standard
- Better error handling and recovery

## Configuration Reference

### Key Configuration Options

```python
PartitionOptions(
    # Processing mode
    validation_level=ValidationLevel.NORMAL,    # LENIENT, NORMAL, STRICT
    blast_only=False,                          # Use only BLAST evidence
    representatives_only=False,                # Process only representative proteins
    
    # Domain constraints  
    min_domain_size=20,                        # Minimum domain size (residues)
    max_domain_size=2000,                      # Maximum domain size (residues)
    
    # Overlap handling
    overlap_threshold=0.3,                     # Maximum allowed overlap (0-1)
    merge_gap_tolerance=20,                    # Merge domains within N residues
    
    # Reference coverage (NEW)
    min_reference_coverage=0.7,                # Minimum coverage required
    strict_reference_coverage=0.9,             # High-confidence threshold
    extend_to_reference_size=True,             # Extend domains to match reference
    
    # Performance
    use_cache=True,                           # Enable classification caching
    parallel_processing=False,                # Use multiple threads/processes
    max_workers=4,                           # Number of parallel workers
    
    # Output
    save_intermediate=False,                  # Save intermediate files
    include_evidence_in_output=True          # Include evidence in domain XML
)
```

## Troubleshooting

### Common Issues

1. **"No domains found"**:
   - Check evidence quality with `analyzer.calculate_quality_metrics()`
   - Lower `min_reference_coverage` or `min_evidence_confidence`
   - Use `ValidationLevel.LENIENT` mode

2. **"Database connection failed"**:
   - Service will continue with reduced functionality
   - Classification assignments may be incomplete
   - Check database credentials and network connectivity

3. **"XML parsing failed"**:
   - Check file permissions and format
   - Analyzer includes automatic repair for common issues
   - Review file creation process if repair fails

4. **High memory usage**:
   - Reduce batch size for batch processing
   - Disable caching if memory is limited
   - Call `service.clear_all_caches()` periodically

5. **Slow performance**:
   - Enable parallel processing for large batches
   - Ensure database connection pooling is configured
   - Use SSD storage for file I/O intensive operations

### Debug Mode

```python
import logging

# Enable detailed logging
logging.getLogger("ecod.pipelines.domain_analysis.partition").setLevel(logging.DEBUG)

# Monitor specific components
logging.getLogger("ecod.pipelines.domain_analysis.partition.processor").setLevel(logging.INFO)

# Get detailed statistics
service_stats = service.get_service_statistics()
print(json.dumps(service_stats, indent=2))
```

## Best Practices

1. **Configuration Management**: Always validate options before use
2. **Error Handling**: Check result.success before accessing results
3. **Resource Management**: Use context managers for services
4. **Batch Processing**: Monitor memory usage and adjust batch sizes
5. **Database Operations**: Handle connection failures gracefully
6. **Performance Monitoring**: Use built-in statistics for optimization
7. **Testing**: Validate with known good datasets before production use

## Dependencies

- **Core ECOD Models**: `ecod.models.pipeline.evidence`, `ecod.models.pipeline.domain`
- **Database**: `ecod.db.DBManager` for classification lookups
- **Configuration**: `ecod.core.context.ApplicationContext`
- **Summary Module**: For evidence input (domain_summary.xml files)
- **Pipeline Models**: For result structures and validation

## License

Part of the ECOD project. See main LICENSE file for details.
