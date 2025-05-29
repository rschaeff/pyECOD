# Evidence Processors Package

The processors package provides a flexible, extensible framework for parsing and extracting evidence from various bioinformatics analysis results. It's designed to convert heterogeneous file formats into standardized `Evidence` objects for downstream analysis in the domain analysis pipeline.

## Overview

The package follows a plugin-based architecture where each processor handles specific file formats and analysis types. All processors inherit from common base classes that provide validation, caching, error handling, and standardized interfaces.

### Core Components

- **Base Classes**: Abstract interfaces and common utilities
- **BLAST Processors**: Handle various BLAST result formats (XML)
- **HHSearch Processors**: Process HHSearch results (XML and HHR formats)
- **Composite Processors**: Combine multiple processors for multi-format support

## Main Processing Modes

### 1. BLAST Evidence Processing

The BLAST processors handle different types of BLAST analyses:

#### Domain BLAST (`BlastEvidenceProcessor`)
Processes domain-vs-domain database searches with support for:
- HSP filtering and e-value thresholds
- Discontinuous domain detection through HSP stitching
- Multi-segment evidence creation for complex domains

```python
from ecod.pipelines.domain_analysis.summary.processors.blast import BlastEvidenceProcessor

# Create processor with custom thresholds
processor = BlastEvidenceProcessor(
    blast_type="domain_blast",
    hsp_evalue_threshold=0.001,
    hit_coverage_threshold=0.7
)

# Process a BLAST XML file
result = processor.process(Path("blast_results.xml"))
print(f"Found {result.evidence_count} evidence items")
```

#### Chain BLAST (`ChainBlastProcessor`)
Specialized for full-chain comparisons with:
- Residue usage tracking to prevent overlaps
- Coverage validation for chain-level matches
- Length difference tolerance checking

```python
from ecod.pipelines.domain_analysis.summary.processors.blast import ChainBlastProcessor

processor = ChainBlastProcessor(
    hsp_evalue_threshold=0.005,
    hit_coverage_threshold=0.8,
    hit_diff_tolerance=50
)

results = processor.process_batch([
    Path("chain1.xml"),
    Path("chain2.xml")
])
```

#### PSI-BLAST (`PSIBlastProcessor`)
Handles iterative PSI-BLAST results:
- Iteration-specific processing
- Configurable iteration ranges
- Evidence tagging with iteration numbers

### 2. HHSearch Evidence Processing

The HHSearch processors support multiple output formats:

#### XML Format Processing
```python
from ecod.pipelines.domain_analysis.summary.processors.hhsearch import HHSearchEvidenceProcessor

processor = HHSearchEvidenceProcessor(
    probability_threshold=50.0,  # Minimum probability (0-100)
    max_hits=100
)

result = processor.process(Path("hhsearch_results.xml"))
```

#### HHR Format Processing
```python
# Same processor handles both XML and HHR formats automatically
result = processor.process(Path("hhsearch_results.hhr"))

# Or use format-specific processor
from ecod.pipelines.domain_analysis.summary.processors.hhsearch import HHRProcessor

hhr_processor = HHRProcessor(probability_threshold=30.0)
result = hhr_processor.process(Path("results.hhr"))
```

## Architecture and Extensibility

### Base Class Hierarchy

```
EvidenceProcessor (ABC)
├── XMLProcessor
│   ├── BlastEvidenceProcessor
│   │   ├── ChainBlastProcessor
│   │   └── PSIBlastProcessor
│   └── HHSearchEvidenceProcessor
│       ├── HHSearchXMLProcessor
│       └── HHRProcessor
└── CompositeProcessor
```

### Key Abstract Methods

Every processor must implement:

```python
@property
@abstractmethod
def supported_formats(self) -> Set[str]:
    """Return supported file extensions"""
    pass

@property
@abstractmethod
def evidence_type(self) -> str:
    """Return evidence type identifier"""
    pass

@abstractmethod
def validate_content(self, content: Any) -> ValidationResult:
    """Validate parsed file content"""
    pass

@abstractmethod
def extract_evidence(self, content: Any, file_path: Path) -> List[Evidence]:
    """Extract Evidence objects from content"""
    pass
```

### Creating Custom Processors

To add support for a new analysis tool, extend the appropriate base class:

```python
from .base import XMLProcessor, ValidationResult
from ecod.models.pipeline import Evidence

class CustomToolProcessor(XMLProcessor):
    """Process CustomTool XML results"""
    
    def __init__(self, score_threshold: float = 10.0, **kwargs):
        super().__init__(**kwargs)
        self.score_threshold = score_threshold
        self.expected_root_tag = "CustomToolOutput"
    
    @property
    def supported_formats(self) -> Set[str]:
        return {'.xml', '.ctool'}
    
    @property
    def evidence_type(self) -> str:
        return "custom_tool"
    
    def validate_content(self, content: ET.Element) -> ValidationResult:
        result = super().validate_content(content)
        if not result.valid:
            return result
            
        # Add custom validation logic
        hits = content.findall(".//Hit")
        result.has_content = len(hits) > 0
        result.metadata['hit_count'] = len(hits)
        
        return result
    
    def extract_evidence(self, content: ET.Element, file_path: Path) -> List[Evidence]:
        evidence_list = []
        
        for hit in content.findall(".//Hit"):
            score = float(hit.get('score', '0'))
            if score < self.score_threshold:
                continue
                
            evidence = Evidence(
                type=self.evidence_type,
                source_id=hit.get('id', ''),
                score=score,
                # ... other attributes
            )
            evidence_list.append(evidence)
        
        return evidence_list
```

### Non-XML Processors

For non-XML formats, inherit directly from `EvidenceProcessor`:

```python
from .base import EvidenceProcessor

class CustomTextProcessor(EvidenceProcessor):
    """Process custom text format"""
    
    @property
    def supported_formats(self) -> Set[str]:
        return {'.txt', '.custom'}
    
    def _parse_file(self, file_path: Path) -> Any:
        """Override to implement custom parsing"""
        with open(file_path, 'r') as f:
            # Custom parsing logic
            return self._parse_custom_format(f.read())
    
    def validate_content(self, content: Any) -> ValidationResult:
        # Custom validation
        pass
    
    def extract_evidence(self, content: Any, file_path: Path) -> List[Evidence]:
        # Custom evidence extraction
        pass
```

## Processing Pipeline Integration

### Batch Processing

All processors support batch operations:

```python
# Process multiple files at once
file_paths = [
    Path("result1.xml"),
    Path("result2.xml"),
    Path("result3.xml")
]

results = processor.process_batch(file_paths)

# Access individual results
for file_path, result in results.items():
    if result.success:
        print(f"{file_path}: {result.evidence_count} evidence items")
    else:
        print(f"{file_path}: Error - {result.error}")
```

### Factory Functions

Use factory functions for dynamic processor creation:

```python
from ecod.pipelines.domain_analysis.summary.processors.blast import create_blast_processor

# Create processor based on analysis type
processor = create_blast_processor(
    blast_type="chain_blast",
    hsp_evalue_threshold=0.001
)
```

### Composite Processing

Combine multiple processors for files with mixed formats:

```python
from .base import CompositeProcessor
from .blast import BlastEvidenceProcessor
from .hhsearch import HHSearchEvidenceProcessor

composite = CompositeProcessor([
    BlastEvidenceProcessor(),
    HHSearchEvidenceProcessor()
])

# Automatically detects format and uses appropriate processor
result = composite.process(file_path)
```

## Configuration and Validation

### Input Validation

All processors perform comprehensive validation:

- **File existence and format checking**
- **Content structure validation**
- **Format-specific requirements**
- **Threshold-based filtering**

```python
# Validate before processing
validation = processor.validate_file(Path("results.xml"))

if validation.valid:
    print("File is valid for processing")
    print(f"File size: {validation.file_size} bytes")
    print(f"Format valid: {validation.format_valid}")
else:
    print(f"Validation failed: {validation.error}")
    for warning in validation.warnings:
        print(f"Warning: {warning}")
```

### Caching and Performance

Processors include built-in caching for improved performance:

```python
# Check cache status
cache_info = processor.get_cache_info()
print(f"Cached files: {cache_info['file_cache_size']}")
print(f"Memory usage estimate: {cache_info['memory_usage_estimate']} chars")

# Clear cache when needed
processor.clear_cache()
```

## Error Handling and Logging

### Processing Results

Every processing operation returns a detailed `ProcessingResult`:

```python
result = processor.process(file_path)

print(f"Status: {result.status}")
print(f"Evidence count: {result.evidence_count}")
print(f"Success: {result.success}")

if result.error:
    print(f"Error: {result.error}")

for warning in result.warnings:
    print(f"Warning: {warning}")

# Access metadata
print(f"File size: {result.metadata.get('file_size', 'unknown')}")
print(f"Evidence type: {result.metadata.get('evidence_type', 'unknown')}")
```

### Logging Configuration

Processors use Python's standard logging:

```python
import logging

# Configure logging for detailed debugging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger('domain_analysis.processors')

processor = BlastEvidenceProcessor(logger=logger)
```

## Best Practices

### Processor Configuration

1. **Set appropriate thresholds** based on your analysis requirements
2. **Use batch processing** for multiple files to improve performance
3. **Validate inputs** before processing in production pipelines
4. **Monitor cache usage** for long-running processes

### Extending the Framework

1. **Follow the abstract interface** defined by base classes
2. **Implement comprehensive validation** for new formats
3. **Add appropriate logging** for debugging and monitoring
4. **Include unit tests** for new processors
5. **Document configuration options** and usage patterns

### Performance Considerations

1. **Use caching wisely** - clear when memory becomes a concern
2. **Process in batches** rather than individual files when possible
3. **Set reasonable limits** (max_hits, file size thresholds)
4. **Profile memory usage** for large-scale processing

## Integration with Pipeline

The processors integrate seamlessly with the broader domain analysis pipeline:

```python
# Typical pipeline integration
from ecod.pipelines.domain_analysis.summary.processors import (
    create_blast_processor,
    HHSearchEvidenceProcessor
)

def process_analysis_results(result_dir: Path):
    """Process all analysis results in a directory"""
    blast_processor = create_blast_processor("domain_blast")
    hhsearch_processor = HHSearchEvidenceProcessor()
    
    all_evidence = []
    
    # Process BLAST results
    for blast_file in result_dir.glob("*.blast.xml"):
        result = blast_processor.process(blast_file)
        if result.success:
            all_evidence.extend(result.evidence)
    
    # Process HHSearch results  
    for hh_file in result_dir.glob("*.hhr"):
        result = hhsearch_processor.process(hh_file)
        if result.success:
            all_evidence.extend(result.evidence)
    
    return all_evidence
```

This modular, extensible design allows the evidence processing system to grow with new analysis tools and formats while maintaining consistent interfaces and robust error handling.
