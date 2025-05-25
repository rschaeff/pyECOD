# Domain Summary Module

## Overview

The Domain Summary module is responsible for collecting, processing, and aggregating evidence from multiple sources (BLAST, HHSearch, self-comparison) to support domain identification and classification in protein structures.

## Architecture

This module follows a clean architecture pattern with clear separation of concerns:

```
summary/
├── service.py          # High-level API
├── generator.py        # Core orchestration logic
├── file_locator.py     # File discovery and validation
├── processors/         # Evidence extraction
├── filters.py          # Quality control
└── models.py          # Data structures
```

### Design Principles

- **Single Responsibility**: Each component has one clear purpose
- **Dependency Injection**: Components are loosely coupled
- **Strategy Pattern**: Processors are pluggable and extensible
- **Value Objects**: Immutable data structures prevent invalid states
- **Service Layer**: Clean public API hides implementation details

## Key Components

### DomainSummaryService

The main entry point for the module. Provides a simple, high-level API for processing proteins.

```python
from ecod.pipelines.domain_analysis.summary import DomainSummaryService

service = DomainSummaryService(context)
result = service.process_protein(pdb_id="1abc", chain_id="A", job_dump_dir="/data/jobs/123")
```

### Evidence Processors

Modular processors that convert different file formats into standardized Evidence objects:

- **BlastEvidenceProcessor**: Processes BLAST XML results
- **HHSearchEvidenceProcessor**: Processes HHSearch XML results  
- **SelfComparisonProcessor**: Processes self-comparison (DALI/HHrepid) results

Each processor implements the `EvidenceProcessor` interface:

```python
class EvidenceProcessor(ABC):
    @abstractmethod
    def process(self, file_path: Path) -> List[Evidence]:
        """Process file and return Evidence objects"""
        
    @abstractmethod
    def validate_file(self, file_path: Path) -> bool:
        """Validate file before processing"""
```

### Evidence Quality Filter

Applies quality control to evidence based on configurable thresholds:

- Confidence score filtering
- Source-specific rules (e.g., chain BLAST gets slightly lower threshold)
- Coverage requirements
- E-value thresholds

### File Locator

Handles file discovery with fallback strategies:

1. Query database for registered file paths
2. Check standard naming conventions
3. Search alternative locations
4. Validate file existence and format

## Usage

### Basic Usage

```python
from ecod.pipelines.domain_analysis.summary import DomainSummaryService
from ecod.core.context import ApplicationContext

# Initialize service
context = ApplicationContext(config_path="config.yml")
service = DomainSummaryService(context)

# Process single protein
result = service.process_protein(
    pdb_id="1abc",
    chain_id="A", 
    job_dump_dir="/data/jobs/123",
    blast_only=False  # Include HHSearch
)

# Check results
if result.success:
    print(f"Found {len(result.domains)} domains")
    print(f"Coverage: {result.coverage:.2%}")
```

### Batch Processing

```python
# Process multiple proteins
proteins = [("1abc", "A"), ("2def", "B"), ("3ghi", "C")]

batch_results = service.process_batch(
    proteins=proteins,
    job_dump_dir="/data/jobs/123",
    blast_only=True,  # BLAST-only mode for speed
    min_confidence=0.5  # Higher confidence threshold
)

print(f"Processed {batch_results.success_count}/{batch_results.total} proteins")
```

### Advanced Configuration

```python
from ecod.pipelines.domain_analysis.summary import SummaryOptions

# Custom options
options = SummaryOptions(
    blast_only=False,
    force_overwrite=True,
    include_evidence_details=True,
    min_confidence=0.3,
    max_gap_for_stitching=30
)

result = service.process_protein(
    pdb_id="1abc",
    chain_id="A",
    job_dump_dir="/data/jobs/123",
    options=options
)
```

## Data Model

### Input

- **ProteinIdentifier**: Immutable identifier for a protein chain
- **SummaryOptions**: Configuration options for processing
- **SequenceInfo**: Protein sequence and metadata

### Output  

- **DomainPartitionResult**: Complete results with domains and evidence
- **Evidence**: Standardized evidence from any source
- **EvidenceSummary**: Aggregated evidence statistics

## File Naming Conventions

The module expects files following these naming patterns:

- FASTA: `{pdb_id}_{chain_id}.fasta`
- Chain BLAST: `{pdb_id}_{chain_id}.chain_blast.xml` 
- Domain BLAST: `{pdb_id}_{chain_id}.domain_blast.xml`
- HHSearch: `{pdb_id}_{chain_id}.{reference}.hhsearch.xml`
- Self-comparison: `{pdb_id}_{chain_id}.self_comp.xml`

## Migration from Legacy Code

### Old Way (Deprecated)

```python
from ecod.pipelines.domain_analysis import DomainSummary

summary = DomainSummary(context)
result = summary.create_summary(pdb_id, chain_id, reference, job_dump_dir)
```

### New Way

```python
from ecod.pipelines.domain_analysis.summary import DomainSummaryService

service = DomainSummaryService(context)
result = service.process_protein(pdb_id, chain_id, job_dump_dir)
```

### Key Differences

1. **Evidence-based**: All results are Evidence objects, not mixed dictionaries/XML
2. **Cleaner API**: Single service entry point instead of multiple classes
3. **Better Error Handling**: Exceptions are caught and reported in results
4. **Configurable**: Fine-grained control over processing options
5. **Testable**: Components can be mocked for unit testing

## Extending the Module

### Adding a New Evidence Source

1. Create a new processor in `processors/`:

```python
from .base import EvidenceProcessor

class MyCustomProcessor(EvidenceProcessor):
    def validate_file(self, file_path: Path) -> bool:
        # Validation logic
        
    def process(self, file_path: Path) -> List[Evidence]:
        # Processing logic
```

2. Register in the generator:

```python
self.processors['my_custom'] = MyCustomProcessor(self.logger)
```

3. Add to evidence collection in `_collect_all_evidence()`

### Custom Filtering Rules

Extend `EvidenceQualityFilter` to add source-specific filtering:

```python
def _filter_my_custom(self, evidence_list: List[Evidence]) -> List[Evidence]:
    # Custom filtering logic
    return filtered_list
```

## Performance Considerations

- File validation is performed before parsing to avoid processing invalid files
- Evidence objects auto-calculate confidence scores on creation
- Batch processing reuses initialized components
- Database queries are minimized through caching

## Error Handling

The module uses a fail-safe approach:

- Individual evidence source failures don't stop processing
- Errors are logged and reported in results
- Missing files are noted but don't cause exceptions
- Invalid evidence is filtered out rather than causing failures

## Testing

The modular design facilitates testing:

```python
# Mock processor for testing
mock_processor = MockEvidenceProcessor(evidence_list)
generator.processors['blast'] = mock_processor

# Test file locator with mock database
mock_db = Mock(DBManager)
locator = EvidenceFileLocator(mock_db, "/test/dir")
```

## Future Enhancements

- Parallel processing of evidence sources
- Caching of processed evidence
- Plugin system for custom processors
- REST API endpoint for remote processing
- Machine learning-based evidence filtering
