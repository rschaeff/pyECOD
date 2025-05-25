# Domain Partition Module - Migration Guide

This guide helps you migrate from the old monolithic `partition.py` to the new service-oriented architecture.

## Overview of Changes

The monolithic `DomainPartition` class has been refactored into several focused components:

- **`DomainPartitionService`**: High-level service interface (replaces main class)
- **`PartitionProcessor`**: Core domain identification algorithms
- **`EvidenceAnalyzer`**: Evidence validation and analysis
- **`StatusTracker`**: Database status tracking
- **`PartitionOptions`**: Configuration management

## Quick Migration Examples

### 1. Basic Single Protein Processing

**Old Code:**
```python
from ecod.pipelines.domain_analysis.partition import DomainPartition

partition = DomainPartition(context)
result = partition.process_protein_domains(
    pdb_id="1abc",
    chain_id="A",
    domain_summary_path="/path/to/summary.xml",
    output_dir="/path/to/output",
    reference="develop291"
)
```

**New Code:**
```python
from ecod.pipelines.domain_analysis.partition import DomainPartitionService

service = DomainPartitionService(context)
result = service.partition_protein(
    pdb_id="1abc",
    chain_id="A",
    summary_path="/path/to/summary.xml",
    output_dir="/path/to/output"
)
```

**Even Simpler:**
```python
from ecod.pipelines.domain_analysis.partition import partition_single_protein

result = partition_single_protein(
    pdb_id="1abc",
    chain_id="A",
    summary_path="/path/to/summary.xml",
    output_dir="/path/to/output"
)
```

### 2. Batch Processing

**Old Code:**
```python
partition = DomainPartition(context)
results = partition.process_batch(
    batch_id=123,
    batch_path="/path/to/batch",
    reference="develop291",
    blast_only=False,
    limit=100,
    reps_only=True
)
```

**New Code:**
```python
service = DomainPartitionService(context)
batch_results = service.partition_batch(
    batch_id=123,
    batch_path="/path/to/batch",
    limit=100,
    representatives_only=True
)

# Access individual results
for result in batch_results.results:
    print(f"{result.pdb_id}_{result.chain_id}: {len(result.domains)} domains")
```

### 3. Custom Configuration

**Old Code:**
```python
partition = DomainPartition(context)
partition.high_confidence_threshold = 0.98
partition.overlap_threshold = 0.2
partition.gap_tolerance = 30
```

**New Code:**
```python
from ecod.pipelines.domain_analysis.partition import (
    DomainPartitionService, PartitionOptions, ValidationLevel
)

# Create custom options
options = PartitionOptions(
    validation_level=ValidationLevel.STRICT,
    overlap_threshold=0.2,
    merge_gap_tolerance=30,
    min_evidence_confidence=0.98,
    use_cache=True,
    save_intermediate=True
)

# Create service with custom options
service = DomainPartitionService(context)
service.default_options = options

# Or pass options per-call
result = service.partition_protein(
    pdb_id="1abc",
    chain_id="A",
    summary_path="/path/to/summary.xml",
    output_dir="/path/to/output",
    overlap_threshold=0.2,
    merge_gap_tolerance=30
)
```

### 4. Direct Evidence Analysis

**Old Code:**
```python
partition = DomainPartition(context)
# Hidden internal methods
evidence = partition._extract_evidence_from_summary(summary_data)
validated = partition._validate_evidence(evidence, "context")
```

**New Code:**
```python
from ecod.pipelines.domain_analysis.partition import EvidenceAnalyzer

analyzer = EvidenceAnalyzer(options)
summary_data = analyzer.parse_domain_summary(summary_path)
evidence_list = analyzer.extract_evidence_with_classification(summary_data)

# Validate individual evidence
for evidence in evidence_list:
    validation = analyzer.validate_evidence(evidence, "my_context")
    if not validation.is_valid:
        print(f"Invalid: {validation.get_summary()}")
```

### 5. Status Tracking

**Old Code:**
```python
# Status tracking was embedded in main class
partition._update_process_status(process_id, "stage", "status")
```

**New Code:**
```python
from ecod.pipelines.domain_analysis.partition import StatusTracker

tracker = StatusTracker(db_manager)
tracker.update_process_status(
    process_id=123,
    stage="domain_partition_processing",
    status="processing"
)

# Get batch progress
progress = tracker.get_batch_progress(batch_id)
print(f"Batch progress: {progress['complete']}/{progress['total']}")
```

## Key Differences

### 1. Configuration Management
- Old: Direct attribute modification
- New: `PartitionOptions` class with validation

### 2. Error Handling
- Old: Mixed return types and exceptions
- New: Consistent `DomainPartitionResult` objects with success/error fields

### 3. Validation
- Old: Scattered validation logic
- New: Centralized in `EvidenceAnalyzer` with configurable levels

### 4. Database Operations
- Old: Mixed with processing logic
- New: Isolated in `StatusTracker`

### 5. Caching
- Old: Class-level dictionaries
- New: Structured `ClassificationCache` with statistics

## Advanced Usage

### Using Components Individually

```python
from ecod.pipelines.domain_analysis.partition import (
    EvidenceAnalyzer, PartitionProcessor, PartitionOptions
)

# Create components
options = PartitionOptions(validation_level=ValidationLevel.NORMAL)
analyzer = EvidenceAnalyzer(options)
processor = PartitionProcessor(options, analyzer)

# Use analyzer standalone
summary_data = analyzer.parse_domain_summary(summary_path)
evidence = analyzer.extract_evidence_with_classification(summary_data)

# Use processor standalone
from ecod.pipelines.domain_analysis.partition import PartitionContext

context = PartitionContext(
    pdb_id="1abc",
    chain_id="A",
    reference="develop291",
    sequence_length=250
)

result = processor.process_evidence(evidence, context)
```

### Service Configuration

```python
# Configure service behavior
service_config = {
    'max_workers': 8,
    'use_multiprocessing': True,
    'batch_size': 50,
    'save_intermediate': True,
    'track_status': True
}

service = DomainPartitionService(context, service_config)
```

### Monitoring and Statistics

```python
# Get service statistics
stats = service.get_service_statistics()
print(f"Proteins processed: {stats['service']['proteins_processed']}")
print(f"Cache hit rate: {stats['analyzer']['hit_rate']:.1f}%")

# Clear caches
service.clear_all_caches()

# Validate setup
validation = service.validate_setup()
if not all(validation.values()):
    print("Setup issues:", [k for k, v in validation.items() if not v])
```

## Backward Compatibility

The old `DomainPartition` class is still available but deprecated:

```python
# This still works but shows deprecation warning
from ecod.pipelines.domain_analysis.partition import DomainPartition

partition = DomainPartition(context)
result = partition.process_protein_domains(...)
```

## Benefits of the New Architecture

1. **Testability**: Each component can be unit tested independently
2. **Flexibility**: Mix and match components as needed
3. **Performance**: Better caching and optional parallel processing
4. **Maintainability**: Clear separation of concerns
5. **Extensibility**: Easy to add new validation rules or processing strategies

## Migration Checklist

- [ ] Update imports to use new service classes
- [ ] Replace direct attribute access with `PartitionOptions`
- [ ] Update method calls to new names
- [ ] Handle `BatchPartitionResults` instead of raw lists
- [ ] Use proper validation levels instead of boolean flags
- [ ] Update error handling to check `result.success`
- [ ] Consider using convenience functions for simple cases
- [ ] Update any custom validation logic to use `EvidenceAnalyzer`
- [ ] Move any database operations to use `StatusTracker`
- [ ] Test thoroughly with your existing data

## Getting Help

For questions or issues during migration:
1. Check the docstrings in the new modules
2. Review the test cases for usage examples
3. Use the backward compatibility wrapper temporarily
4. Report any missing functionality as an issue
