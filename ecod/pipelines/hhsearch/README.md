# ECOD Pipeline Services Documentation

## Table of Contents
1. [HHSearch Pipeline Service](#hhsearch-pipeline-service)
2. [Pipeline Orchestration Service](#pipeline-orchestration-service)

---

## HHSearch Pipeline Service

### Overview

The HHSearch Pipeline Service handles the processing, validation, and registration of HHSearch results within the ECOD system. It converts HHR (HHSearch result) files to standardized XML format and manages file discovery across standard and legacy locations.

### Architecture

```
ecod/pipelines/hhsearch/
├── __init__.py          # Package exports
├── file_manager.py      # File discovery and migration
├── models.py            # Data models
├── processor.py         # HHR parsing and XML conversion
└── service.py           # High-level registration service
```

### Core Components

#### 1. HHSearchRegistrationService
The main service class that coordinates HHSearch result processing.

**Key Features:**
- Batch and individual chain processing
- Parallel processing support
- File validation and migration
- Progress tracking and statistics

**Usage:**
```python
from ecod.pipelines.hhsearch.service import create_service

# Create service instance
service = create_service(config_path='config/config.yml')

# Register results for a batch
results = service.register_batch(
    batch_id=123,
    limit=100,  # Optional: process only 100 chains
    chain_ids=['1abc_A', '2def_B']  # Optional: specific chains
)

# Check results
print(f"Registered: {results.registered}")
print(f"Failed: {results.failed}")
print(f"Success rate: {results.success_rate}%")
```

#### 2. HHRToXMLConverter
Converts parsed HHR data to XML format.

**Key Features:**
- Configurable probability threshold filtering
- Preserves alignment details and consensus sequences
- Generates well-formatted XML with metadata

**XML Structure:**
```xml
<hh_summ_doc>
    <metadata>
        <pdb_id>1abc</pdb_id>
        <chain_id>A</chain_id>
        <reference>ecod70</reference>
        <creation_date>2024-01-15 10:30:00</creation_date>
        <min_probability>20.0</min_probability>
    </metadata>
    <hh_hit_list>
        <hh_hit hit_num="1" hit_id="e1xyzA1" probability="99.5" 
                e_value="1.2e-10" score="125.3">
            <query_range start="10" end="150">10-150</query_range>
            <template_seqid_range start="5" end="145" 
                                  identity="35.2" 
                                  coverage="0.95">5-145</template_seqid_range>
            <alignment>
                <query_ali>MKLTPEHVFI...</query_ali>
                <template_ali>MKLSPDHVFI...</template_ali>
                <match_sequence>MKL+P+HVFI...</match_sequence>
            </alignment>
        </hh_hit>
    </hh_hit_list>
</hh_summ_doc>
```

#### 3. HHSearchFileManager
Manages file discovery across standard and legacy locations.

**Key Features:**
- Searches multiple file location patterns
- Validates HHSearch vs HHblits output
- Supports file migration to standard locations

#### 4. HHSearchProcessor
Integrates HHSearch results with other evidence to create domain summaries.

**Key Methods:**
- `process_batch()` - Process all chains in a batch
- `process_protein()` - Process a specific protein
- `_process_chain()` - Core processing logic for a single chain

### Data Models

#### RegistrationStatus
```python
class RegistrationStatus(Enum):
    PENDING = "pending"
    FOUND = "found"
    CONVERTING = "converting"
    REGISTERED = "registered"
    FAILED = "failed"
    SKIPPED = "skipped"
```

#### ServiceConfig
```python
@dataclass
class ServiceConfig:
    force_regenerate: bool = False      # Force XML regeneration
    min_probability: float = 20.0       # Min HHSearch probability
    validate_hhsearch: bool = True      # Validate HHR files
    migrate_legacy_files: bool = True   # Migrate from legacy paths
    parallel_processing: bool = True    # Enable parallel processing
    max_workers: int = 4               # Max parallel workers
    process_timeout: int = 300         # Timeout in seconds
```

### Workflow

1. **File Discovery**
   - Search standard paths first
   - Check legacy locations if needed
   - Validate file content

2. **Validation**
   - Ensure file is HHSearch output (not HHblits)
   - Check for ECOD domain hits
   - Verify file integrity

3. **Conversion**
   - Parse HHR file
   - Filter hits by probability threshold
   - Generate structured XML

4. **Registration**
   - Save XML to standard location
   - Register files in database
   - Update process status

### Configuration

The service uses configuration from `config.yml`:

```yaml
hhsearch:
  min_probability: 20.0
  validate_files: true
  parallel:
    enabled: true
    max_workers: 4
```

---

