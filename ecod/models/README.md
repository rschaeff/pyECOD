# ECOD Models Module

This module provides data models for the ECOD (Evolutionary Classification of protein Domains) pipeline, with a focus on domain analysis, evidence processing, and result management.

## 🚀 Quick Start - Use These Models

**For new development, use the gold standard pipeline models:**

```python
# ✅ RECOMMENDED - Gold Standard Models
from ecod.models import Evidence, DomainModel, DomainPartitionResult

# Create evidence
evidence = Evidence(
    type="hhsearch",
    domain_id="e1abc12",
    probability=0.95,
    query_range="10-150"
)

# Create domain
domain = DomainModel(
    id="domain_1",
    start=10,
    end=150,
    range="10-150",
    evidence=[evidence]
)

# Create partition result
result = DomainPartitionResult(
    pdb_id="1abc",
    chain_id="A",
    reference="develop291",
    domains=[domain]
)
```

## 📋 Model Priority and Status

### 🟢 Active Models (Use These)

| Model Category | Models | Purpose |
|---------------|--------|---------|
| **Pipeline Models** | `Evidence`, `DomainModel`, `DomainPartitionResult` | Domain analysis and processing |
| **Infrastructure** | `XmlSerializable` | Base classes and interfaces |
| **Job Management** | `Batch`, `ProcessStatus`, `Job`, `JobItem` | Workflow and processing management |

### 🟡 Legacy Models (Deprecated)

| Model Category | Models | Status | Migration Path |
|---------------|--------|---------|---------------|
| **Old Pipeline** | `BlastHit`, `HHSearchHit`, `DomainSummaryModel` | Deprecated with warnings | → Use `Evidence` |
| **Protein Models** | `Protein`, `ProteinSequence`, `ProteinStructure` | Being deprecated | → Use database models directly |

## 🏗️ Architecture Overview

```
ecod/models/
├── __init__.py                 # Main entry point
├── base.py                     # Base classes (XmlSerializable)
├── job.py                      # Job management models (active)
├── pipeline.py                 # Legacy models (deprecated)
├── protein.py                  # Protein models (deprecated)
└── pipeline/                   # Gold standard models
    ├── __init__.py
    ├── evidence.py             # Unified evidence model
    ├── domain.py               # Domain model
    └── partition.py            # Partition results
```

## 🎯 Gold Standard Pipeline Models

### Evidence Model

The `Evidence` model consolidates all evidence types (BLAST, HHSearch, etc.) into a unified interface:

```python
from ecod.models import Evidence

# HHSearch evidence
hhsearch_evidence = Evidence(
    type="hhsearch",
    domain_id="e1abc12",
    probability=0.95,
    query_range="10-150",
    t_group="1.1.1",
    h_group="1.1"
)

# BLAST evidence  
blast_evidence = Evidence(
    type="domain_blast",
    domain_id="e2def34",
    evalue=1e-20,
    query_range="50-200"
)

# Confidence is automatically calculated
print(f"HHSearch confidence: {hhsearch_evidence.confidence:.3f}")
print(f"BLAST confidence: {blast_evidence.confidence:.3f}")
```

### Domain Model

The `DomainModel` provides comprehensive domain representation:

```python
from ecod.models import DomainModel

domain = DomainModel(
    id="domain_1",
    start=10,
    end=150,
    range="10-150",
    source="hhsearch",
    t_group="1.1.1",
    h_group="1.1",
    evidence=[hhsearch_evidence, blast_evidence]
)

# Rich functionality
print(f"Domain size: {domain.size}")
print(f"Classification level: {domain.get_classification_level()}")
print(f"Confidence: {domain.confidence:.3f}")
print(f"Is fully classified: {domain.is_fully_classified()}")
```

### Partition Result Model

The `DomainPartitionResult` manages complete analysis results:

```python
from ecod.models import DomainPartitionResult

result = DomainPartitionResult(
    pdb_id="1abc",
    chain_id="A", 
    reference="develop291",
    domains=[domain],
    sequence_length=200
)

# Automatic analysis
print(f"Coverage: {result.coverage:.1%}")
print(f"Classification status: {result.is_classified}")
print(f"Domain count: {len(result.domains)}")

# Save result
result.save("/path/to/output")
```

## 🔄 Migration Guide

### From BlastHit/HHSearchHit → Evidence

```python
# ❌ OLD (deprecated)
from ecod.models import BlastHit, HHSearchHit

blast_hit = BlastHit.from_xml(xml_element)
hhsearch_hit = HHSearchHit.from_xml(xml_element)

# ✅ NEW (recommended)
from ecod.models import Evidence

blast_evidence = Evidence.from_blast_xml(xml_element, "domain_blast")
hhsearch_evidence = Evidence.from_hhsearch_xml(xml_element)
```

### From Multiple Domain Models → DomainModel

```python
# ❌ OLD (deprecated)
from ecod.models.domain import DomainModel as OldDomainModel
from ecod.models.domain_analysis.domain_model import DomainModel as AnalysisDomainModel

# ✅ NEW (consolidated)
from ecod.models import DomainModel

domain = DomainModel(
    id="domain_1",
    start=10,
    end=150,
    range="10-150"
)
```

### From Protein Models → Database Direct

```python
# ❌ OLD (being deprecated)
from ecod.models import Protein, ProteinSequence

# ✅ NEW (recommended)
# Access protein data directly through database layer
# Use database models or raw queries instead
```

## 📚 Key Features

### Automatic Confidence Calculation

Evidence and domains automatically calculate confidence scores:

```python
evidence = Evidence(
    type="hhsearch",
    probability=0.95  # High probability
)
# Confidence automatically calculated: ~0.95

evidence = Evidence(
    type="domain_blast", 
    evalue=1e-20  # Very good E-value
)
# Confidence automatically calculated: ~0.85
```

### XML Serialization

All pipeline models support XML serialization:

```python
# To XML
xml_element = evidence.to_xml()
xml_string = ET.tostring(xml_element, encoding='unicode')

# From XML
evidence = Evidence.from_xml(xml_element)
```

### Dictionary Conversion

Models provide dictionary conversion for APIs and serialization:

```python
evidence_dict = evidence.to_dict()
evidence_restored = Evidence.from_dict(evidence_dict)

domain_dict = domain.to_dict()
domain_restored = DomainModel.from_dict(domain_dict)
```

## 🛠️ Advanced Usage

### Custom Evidence Types

```python
custom_evidence = Evidence(
    type="custom_method",
    source_id="custom_123",
    confidence=0.8,  # Explicit confidence
    extra_attributes={
        "method_specific_score": 42,
        "custom_parameter": "value"
    }
)
```

### Domain Analysis

```python
# Check overlaps between domains
if domain1.overlaps(domain2):
    overlap_size = domain1.overlap_size(domain2)
    print(f"Domains overlap by {overlap_size} residues")

# Merge compatible domains
merged_domain = domain1.merge_with(domain2)

# Split domain at position
part1, part2 = domain.split_at(100)
```

### Partition Analysis

```python
result = DomainPartitionResult(...)

# Get analysis statistics
stats = result.get_summary_stats()
print(f"Average confidence: {stats['average_confidence']:.3f}")

# Find overlapping domains
overlaps = result.get_overlapping_domains()

# Get domains by source
hhsearch_domains = result.get_domains_by_source("hhsearch")
```

## ⚠️ Important Notes

### Deprecation Timeline

- **Legacy pipeline.py models**: Deprecated, emit warnings
- **Protein models**: Being phased out  
- **Migration deadline**: Target removal in next major version

### Thread Safety

Models are generally thread-safe for reading, but not for concurrent modification. Use appropriate locking for multi-threaded access.

### Performance Considerations

- Evidence confidence calculation has computational cost
- Use explicit confidence setting for performance-critical code
- Batch processing is more efficient than individual model operations

## 🧪 Testing

```python
# Test evidence creation
evidence = Evidence(type="test", confidence=0.5)
assert evidence.confidence == 0.5

# Test domain functionality  
domain = DomainModel(id="test", start=1, end=100, range="1-100")
assert domain.size == 100
assert not domain.is_classified()  # No classification set

# Test partition results
result = DomainPartitionResult(
    pdb_id="test", 
    chain_id="A", 
    reference="test",
    domains=[domain]
)
assert len(result.domains) == 1
```

## 📖 Examples and Patterns

For complete usage examples, refer to the actual implementation in the pipeline modules and their test files.

## 🤝 Contributing

When adding new models:

1. Inherit from `XmlSerializable` for serialization support
2. Provide `to_dict()` and `from_dict()` methods
3. Add comprehensive type hints
4. Include validation in `__post_init__`
5. Write tests for all functionality
6. Update this README

## 📋 Version Information

- **Model Version**: 2.0.0 (Consolidated)
- **Status**: Gold Standard
- **Python Requirements**: 3.8+
- **Dependencies**: xml.etree.ElementTree, dataclasses, typing

---

**Quick Reference**: Use `Evidence`, `DomainModel`, and `DomainPartitionResult` for all new development. Avoid deprecated models in `pipeline.py` and `protein.py`.
