## Pipeline Orchestration Service

### Overview

The Pipeline Orchestration Service provides high-level coordination of the entire ECOD pipeline with a deterministic, easy-to-understand workflow. It manages the execution of pipeline stages, tracks progress, and handles error recovery.

### Architecture

```
ecod/pipelines/orchestration/
├── models.py         # Data models for orchestration
├── service.py        # Main orchestration service
└── stage_manager.py  # Stage dependency management
```

### Pipeline Modes

#### 1. FULL Mode
Complete pipeline execution including all stages:
- BLAST searches (chain and domain)
- HHSearch profile generation and searching
- Domain analysis and partitioning

#### 2. BLAST_ONLY Mode
Simplified pipeline without HHSearch:
- BLAST searches only
- Domain analysis based on BLAST evidence
- Faster execution for preliminary analysis

#### 3. ANALYSIS_ONLY Mode
Analysis using existing evidence:
- No new searches performed
- Uses existing BLAST/HHSearch results
- Domain summary and partitioning only

### Core Components

#### 1. PipelineOrchestrationService
Main service class that coordinates pipeline execution.

**Usage:**
```python
from ecod.pipelines.orchestration.service import create_orchestrator

# Create orchestrator
orchestrator = create_orchestrator(config_path='config/config.yml')

# Run full pipeline
run = orchestrator.run_pipeline(
    batch_id=123,
    mode=PipelineMode.FULL,
    force_restart=False  # Resume from checkpoint
)

# Check status
status = orchestrator.get_batch_status(batch_id=123)
print(f"Current stage: {status['active_run']['current_stage']}")

# Resume interrupted pipeline
run = orchestrator.resume_pipeline(batch_id=123)
```

#### 2. StageManager
Manages stage dependencies and execution order.

**Stage Dependencies:**
```
INITIALIZATION
├── BLAST_CHAIN
├── BLAST_DOMAIN
│   └── DOMAIN_SUMMARY
│       └── DOMAIN_PARTITION
│           └── CLASSIFICATION
└── HHSEARCH_PROFILE (requires BLAST_CHAIN + BLAST_DOMAIN)
    └── HHSEARCH_SEARCH
        └── HHSEARCH_REGISTRATION
            └── DOMAIN_SUMMARY
```

### Pipeline Stages

#### PipelineStage Enum
```python
class PipelineStage(Enum):
    INITIALIZATION = "initialization"
    BLAST_CHAIN = "blast_chain"
    BLAST_DOMAIN = "blast_domain"
    HHSEARCH_PROFILE = "hhsearch_profile"
    HHSEARCH_SEARCH = "hhsearch_search"
    HHSEARCH_REGISTRATION = "hhsearch_registration"
    DOMAIN_SUMMARY = "domain_summary"
    DOMAIN_PARTITION = "domain_partition"
    CLASSIFICATION = "classification"
    COMPLETE = "complete"
```

### Data Models

#### PipelineRun
Tracks a complete pipeline execution:
```python
@dataclass
class PipelineRun:
    run_id: str
    batch_id: int
    mode: PipelineMode
    total_proteins: int
    start_time: datetime
    end_time: Optional[datetime]
    stages: List[StageResult]
    
    # Properties
    is_complete: bool
    duration: float
    success: bool
    current_stage: Optional[PipelineStage]
```

#### StageResult
Results from executing a pipeline stage:
```python
@dataclass
class StageResult:
    stage: PipelineStage
    success: bool
    start_time: datetime
    end_time: Optional[datetime]
    items_processed: int
    items_failed: int
    error: Optional[str]
    details: Dict[str, Any]
```

#### OrchestrationConfig
```python
@dataclass
class OrchestrationConfig:
    # Execution settings
    checkpoint_enabled: bool = True     # Save progress
    continue_on_error: bool = True      # Continue if stage fails
    parallel_stages: bool = False       # Parallel stage execution
    
    # Stage settings
    blast_batch_size: int = 100
    hhsearch_threads: int = 8
    hhsearch_memory: str = "16G"
    domain_analysis_batch_size: int = 100
    
    # Timeouts
    stage_timeout: int = 3600          # 1 hour
    total_timeout: int = 86400         # 24 hours
```

### Workflow Example

```python
# 1. Create orchestrator with custom config
config = {
    'checkpoint_enabled': True,
    'continue_on_error': True,
    'hhsearch_threads': 16,
    'hhsearch_memory': '32G'
}
orchestrator = create_orchestrator(orchestration_config=config)

# 2. Run pipeline
run = orchestrator.run_pipeline(
    batch_id=123,
    mode=PipelineMode.FULL
)

# 3. Monitor progress
while not run.is_complete:
    status = orchestrator.get_batch_status(123)
    print(f"Stage: {status['active_run']['current_stage']}")
    print(f"Progress: {status['stages']}")
    time.sleep(60)

# 4. Get final results
summary = run.get_summary()
print(f"Pipeline completed in {summary['duration']}s")
print(f"Success: {summary['success']}")
```

### Error Handling and Recovery

The orchestration service provides robust error handling:

1. **Checkpoint System**
   - Progress saved after each stage
   - Can resume from last successful stage
   - No reprocessing of completed work

2. **Error Recovery**
   - `continue_on_error` allows pipeline to continue despite failures
   - Failed items tracked separately
   - Detailed error messages preserved

3. **Stage Status Tracking**
   ```python
   status = orchestrator.get_batch_status(batch_id)
   # Returns:
   {
       'batch': {...},
       'stages': {
           'blast_chain': {
               'completed': True,
               'stats': {'total': 100, 'completed': 100}
           },
           'hhsearch_search': {
               'in_progress': True,
               'stats': {'total': 100, 'completed': 45}
           }
       },
       'active_run': {...}
   }
   ```

### Convenience Functions

```python
# Run full pipeline
from ecod.pipelines.orchestration.service import run_full_pipeline
result = run_full_pipeline(batch_id=123)

# Run BLAST-only pipeline
from ecod.pipelines.orchestration.service import run_blast_only_pipeline
result = run_blast_only_pipeline(batch_id=123)
```

### Best Practices

1. **Use Checkpoints**: Enable `checkpoint_enabled` for long-running pipelines
2. **Monitor Progress**: Check `get_batch_status()` regularly
3. **Handle Failures**: Set `continue_on_error=True` for resilient execution
4. **Resource Management**: Adjust `max_workers` and memory settings based on system
5. **Logging**: Monitor logs for detailed execution information

### Integration with Other Services

The orchestration service integrates with:
- **BlastPipeline**: For BLAST searches
- **HHSearchPipeline**: For HMM profile searches
- **HHSearchRegistrationService**: For result processing
- **DomainSummaryService**: For evidence integration
- **DomainPartitionService**: For domain boundary determination
