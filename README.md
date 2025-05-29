# pyECOD: Evolutionary Classification of Protein Domains Pipeline

**Status: Active Development** | **Not Released**

A comprehensive computational pipeline for automated protein domain identification and classification using the ECOD (Evolutionary Classification of Protein Domains) hierarchical system. pyECOD processes protein sequences through multiple evidence sources to identify domain boundaries and assign ECOD classifications.

## Overview

pyECOD transforms protein sequences into classified domain models through a multi-stage pipeline that combines:
- **BLAST searches** against ECOD reference databases
- **HHSearch profile-based** searches for remote homology detection
- **Evidence integration** from multiple computational approaches
- **Domain boundary identification** using sophisticated algorithms
- **ECOD classification assignment** following the hierarchical system (A > X > H > T groups)

### What pyECOD Does

```
Input: Protein Sequences (FASTA)
    ↓
Evidence Collection:
├── BLAST searches (chain + domain databases)
├── HHSearch profile generation & searching
└── Self-comparison analysis
    ↓
Evidence Processing:
├── Quality filtering & validation
├── Confidence scoring
└── Coverage analysis
    ↓
Domain Analysis:
├── Boundary identification
├── Overlap resolution
└── Classification assignment
    ↓
Output: Classified Domain Models (XML/Database)
```

## Architecture Overview

pyECOD follows a modern, service-oriented architecture with clear separation of concerns:

```
pyECOD/
├── ecod/
│   ├── config/          # YAML configuration system
│   ├── db/              # PostgreSQL database layer
│   ├── jobs/            # Job management (local/SLURM)
│   ├── core/            # Application context & orchestration
│   ├── models/          # Data models (Evidence, Domain, etc.)
│   ├── pipelines/       # Processing pipelines
│   │   ├── blast_pipeline.py
│   │   ├── hhsearch_pipeline.py
│   │   ├── domain_analysis/    # Evidence → Domains
│   │   └── orchestration/      # Pipeline coordination
│   └── cli/             # Command-line interface (in development)
└── scripts/             # Legacy scripts (being migrated to CLI)
```

### Core Components

#### 1. Configuration System (`ecod.config`)
- **Hierarchical YAML configuration** with environment variable overrides
- **Secure secrets management** via local config files
- **Multi-environment support** (development, production, testing)

#### 2. Database Layer (`ecod.db`)
- **PostgreSQL-based** with repository pattern
- **Schema migrations** and data migration utilities
- **Connection pooling** and transaction management
- **Type-safe queries** with comprehensive error handling

#### 3. Job Management (`ecod.jobs`)
- **Multiple backends**: Local threading and SLURM cluster
- **Database integration** for progress tracking
- **Batch processing** with configurable parallelization
- **Error recovery** and job monitoring

#### 4. Data Models (`ecod.models`)
- **Unified Evidence model** consolidating BLAST, HHSearch, and other sources
- **Rich Domain model** with classification hierarchy
- **Partition Results** with coverage statistics and metadata
- **XML serialization** and database persistence

#### 5. Processing Pipelines (`ecod.pipelines`)

**BLAST Pipeline**:
- Chain-level and domain-level database searches
- Batch job submission and monitoring
- Standardized result parsing and validation

**HHSearch Pipeline**:
- HMM profile generation via HHblits
- Profile-based database searching
- Result conversion and registration

**Domain Analysis Pipeline**:
- **Summary Service**: Evidence collection, validation, and integration
- **Partition Service**: Domain boundary identification and classification
- **Evidence Processors**: Pluggable parsers for different analysis tools
- **Quality Control**: Multi-level validation and confidence scoring
- **Reference Coverage**: Validation against known domain structures
- **Overlap Resolution**: Sophisticated conflict resolution algorithms

**Orchestration Service**:
- **Multi-mode execution**: FULL, BLAST_ONLY, ANALYSIS_ONLY
- **Checkpoint system** for resumable processing
- **Stage management** with dependency tracking
- **Error recovery** and progress monitoring

## Data Flow

### Stage 1: Evidence Collection
```
Protein Sequence (FASTA)
    ↓
┌─── BLAST Chain Search ───────┐
│    └→ Chain similarities      │
├─── BLAST Domain Search ──────┤
│    └→ Domain matches         │ → Evidence Integration
├─── HHSearch Profile Gen ─────┤
│    └→ HMM Profile (.a3m)     │
└─── HHSearch Database Search ─┘
     └→ Remote homologs (.hhr)
```

### Stage 2: Evidence Processing
```
Raw Evidence Files
    ↓
┌─── Validation & Quality Control ───┐
│    ├→ E-value thresholds           │
│    ├→ Coverage requirements        │
│    └→ Confidence scoring           │
├─── Evidence Unification ──────────┤
│    └→ Standardized Evidence objects│
└─── Classification Lookup ─────────┘
     └→ ECOD hierarchy integration
```

### Stage 3: Domain Analysis
```
Validated Evidence
    ↓
┌─── Boundary Identification ────┐
│    ├→ Coverage analysis        │
│    ├→ Overlap resolution       │
│    └→ Discontinuous domains    │
├─── Classification Assignment ──┤
│    ├→ T-group (topology)       │
│    ├→ H-group (homology)       │
│    ├→ X-group (possible)       │
│    └→ A-group (architecture)   │
└─── Quality Assessment ────────┘
     └→ Coverage, confidence stats
```

## Domain Analysis: The Core Logic

The domain analysis pipeline represents the most sophisticated component of pyECOD, where raw computational evidence is transformed into biologically meaningful domain models. This two-stage process combines evidence integration with advanced boundary detection algorithms.

### Two-Stage Architecture

#### Stage 1: Evidence Summary (`ecod.pipelines.domain_analysis.summary`)

**Purpose**: Transform heterogeneous analysis results into standardized, validated evidence objects.

```python
# Evidence processing workflow
from ecod.pipelines.domain_analysis.summary import DomainSummaryService

service = DomainSummaryService(context)
result = service.process_protein(
    pdb_id="1abc",
    chain_id="A",
    job_dump_dir="/data/batch_123"
)

# Output: Unified evidence summary XML file
# Contains: All evidence objects with confidence scores
```

**Key Components**:

**Evidence Processors** - Pluggable architecture for different analysis tools:
```python
# BLAST Evidence Processor
blast_processor = BlastEvidenceProcessor(
    blast_type="domain_blast",
    hsp_evalue_threshold=0.001,    # E-value filtering
    hit_coverage_threshold=0.7,    # Coverage requirements
    max_gap_for_stitching=30      # HSP stitching for discontinuous domains
)

# HHSearch Evidence Processor
hhsearch_processor = HHSearchEvidenceProcessor(
    probability_threshold=50.0,    # Minimum HHSearch probability
    max_hits=100,                 # Limit number of hits processed
    validate_alignments=True      # Alignment quality validation
)
```

**Evidence Quality Control**:
- **Multi-level validation**: File format, content structure, biological plausibility
- **Confidence scoring**: Automatic calculation based on statistical significance
- **Source-specific filtering**: Different thresholds for BLAST vs HHSearch vs self-comparison
- **Coverage validation**: Ensures evidence adequately spans predicted domains

**File Discovery and Validation**:
```python
# Intelligent file location with fallback strategies
file_locator = EvidenceFileLocator(context.db, job_dump_dir)

# 1. Database-registered file paths
# 2. Standard naming conventions
# 3. Alternative locations (legacy paths)
# 4. Format validation before processing
```

#### Stage 2: Domain Partition (`ecod.pipelines.domain_analysis.partition`)

**Purpose**: Convert validated evidence into domain boundaries with ECOD classifications.

```python
# Domain partitioning workflow
from ecod.pipelines.domain_analysis.partition import DomainPartitionService

service = DomainPartitionService(context)
result = service.partition_protein(
    pdb_id="1abc",
    chain_id="A",
    summary_path="/data/batch_123/1abc_A.evidence_summary.xml",
    output_dir="/data/output"
)

# Output: Domain models with boundaries and classifications
# Format: XML file with DomainModel objects
```

### Advanced Domain Boundary Algorithms

#### Evidence Analysis and Validation
```python
# Multi-criteria evidence assessment
evidence_analyzer = EvidenceAnalyzer(options)

validation_result = evidence_analyzer.validate_evidence(evidence, context)
# Checks:
# - Statistical significance (E-values, probabilities)
# - Biological plausibility (length, coverage)
# - Consistency across evidence sources
# - Reference domain compatibility
```

#### Reference Coverage Validation
```python
# Ensures evidence adequately covers known domain structures
coverage_analyzer = ReferenceCoverageAnalyzer(context.db)

enhanced_evidence = coverage_analyzer.analyze_evidence_coverage(evidence)
# Analysis includes:
# - Coverage percentage against reference domain
# - Boundary extension suggestions
# - Discontinuous domain detection
# - Multi-segment domain handling
```

#### Boundary Identification Logic
```python
# Core domain boundary detection
processor = PartitionProcessor(context, options)

# 1. Evidence clustering by sequence position
position_clusters = processor.cluster_evidence_by_position(evidence_list)

# 2. Boundary candidate generation
boundary_candidates = processor.generate_boundary_candidates(
    clusters, sequence_length
)

# 3. Optimal boundary selection
optimal_boundaries = processor.select_optimal_boundaries(
    candidates, evidence_quality_metrics
)
```

#### Overlap Resolution Algorithms
```python
# Sophisticated conflict resolution when domains overlap
overlap_resolver = OverlapResolver(options)

resolved_domains = overlap_resolver.resolve_overlaps(
    domain_candidates,
    resolution_strategy="confidence_based"  # or "coverage_based", "length_based"
)

# Resolution strategies:
# - High-confidence domains are protected
# - Lower-confidence domains adjusted or merged
# - Biological constraints enforced (minimum domain size)
# - Gap tolerance for discontinuous domains
```

#### Discontinuous Domain Handling
```python
# Special processing for domains with gaps/insertions
discontinuous_processor = DiscontinuousDomainProcessor()

# Detect multi-segment domains
segments = discontinuous_processor.identify_segments(evidence, gaps_threshold=50)

# Create discontinuous domain candidates
discontinuous_domains = discontinuous_processor.create_candidates(
    segments, min_segment_size=20
)
```

### Classification Assignment Logic

#### ECOD Hierarchy Integration
```python
# Database-driven classification lookup with caching
classification_cache = ClassificationCache(context.db)

for domain in domains:
    # Retrieve ECOD classifications from evidence
    classifications = classification_cache.get_classifications(
        domain.evidence_domain_ids
    )

    # Assign hierarchical classifications
    domain.t_group = select_consensus_classification(classifications, "t_group")
    domain.h_group = select_consensus_classification(classifications, "h_group")
    domain.x_group = select_consensus_classification(classifications, "x_group")
    domain.a_group = select_consensus_classification(classifications, "a_group")
```

#### Confidence-Based Classification
```python
# Classification confidence based on evidence agreement
def calculate_classification_confidence(domain):
    evidence_groups = group_evidence_by_classification(domain.evidence)

    # Metrics:
    # - Evidence consensus (how many sources agree)
    # - Source diversity (BLAST + HHSearch + self-comparison)
    # - Statistical significance (combined E-values/probabilities)
    # - Reference coverage (how well evidence spans domain)

    return weighted_confidence_score(evidence_groups)
```

### Quality Control and Validation

#### Multi-Level Validation
```python
# Validation levels: LENIENT, NORMAL, STRICT
validation_options = PartitionOptions(
    validation_level=ValidationLevel.STRICT,
    min_reference_coverage=0.8,      # Require 80% reference coverage
    min_evidence_confidence=0.5,     # Minimum evidence confidence
    overlap_threshold=0.3,           # Maximum allowed overlap
    min_domain_size=20,              # Minimum domain length
    max_domain_size=2000             # Maximum domain length
)
```

#### Quality Metrics Calculation
```python
# Comprehensive quality assessment
quality_metrics = analyzer.calculate_quality_metrics(
    evidence_list, sequence_length, processing_time
)

# Metrics include:
# - sequence_coverage: Fraction of sequence covered by evidence
# - evidence_diversity: Number of different evidence sources
# - classification_consistency: Agreement between evidence classifications
# - statistical_significance: Combined significance scores
# - boundary_confidence: Reliability of predicted boundaries
```

### Processing Modes and Optimization

#### Evidence Source Routing
```python
# Intelligent processing path assignment based on evidence quality
router = ProcessingRouter(context)

processing_paths = router.assign_processing_paths(batch_id=123)
# Returns:
{
    "blast_only": [protein_ids],      # High-confidence BLAST evidence
    "full_pipeline": [protein_ids],   # Need HHSearch for classification
    "manual_review": [protein_ids]    # Complex cases requiring inspection
}
```

#### Batch Processing Optimization
```python
# Parallel processing with resource management
batch_results = service.partition_batch(
    batch_id=123,
    batch_path="/data/batch_123",
    limit=1000,                      # Process first 1000 proteins
    representatives_only=True,        # Focus on representative structures
    options=PartitionOptions(
        parallel_processing=True,     # Enable multiprocessing
        max_workers=8,               # Parallel worker count
        use_cache=True,              # Classification caching
        track_status=True            # Database status tracking
    )
)
```

### Error Handling and Recovery

#### Robust Error Management
```python
# Processing continues despite individual failures
try:
    result = service.partition_protein(pdb_id, chain_id, summary_path, output_dir)
    if result.success:
        domains_found += len(result.domains)
    else:
        failed_proteins.append((pdb_id, chain_id, result.error))
except Exception as e:
    # Log error but continue batch processing
    logger.error(f"Unexpected error for {pdb_id}_{chain_id}: {e}")
    continue
```

#### Status Tracking and Recovery
```python
# Database-driven progress tracking
status_tracker = StatusTracker(context.db)

# Track processing lifecycle
status_tracker.start_process(process_id, pdb_id, chain_id, batch_id)
status_tracker.update_process_stage(process_id, "evidence_analysis", "in_progress")
status_tracker.update_process_stage(process_id, "boundary_identification", "completed")

# Recovery capabilities
failed_processes = status_tracker.get_failed_processes(batch_id)
retry_results = service.reprocess_failed(batch_id, batch_path)
```

### Advanced Features

#### Configurable Processing Options
```python
# Flexible configuration for different use cases
research_options = PartitionOptions(
    validation_level=ValidationLevel.LENIENT,     # More permissive
    min_reference_coverage=0.6,                   # Lower coverage requirement
    extend_to_reference_size=True,                # Extend boundaries to match reference
    merge_gap_tolerance=30,                       # Merge nearby domains
    include_evidence_in_output=True               # Detailed evidence in output
)

production_options = PartitionOptions(
    validation_level=ValidationLevel.STRICT,      # Higher quality standards
    min_reference_coverage=0.9,                   # Strict coverage requirement
    parallel_processing=True,                     # High throughput
    use_cache=True,                              # Performance optimization
    track_status=True                            # Monitoring
)
```

#### Integration with Broader Pipeline
```python
# Seamless integration with orchestration service
from ecod.pipelines.orchestration.service import create_orchestrator

orchestrator = create_orchestrator()
run = orchestrator.run_pipeline(
    batch_id=123,
    mode=PipelineMode.FULL  # Includes complete domain analysis
)

# Domain analysis automatically receives:
# - BLAST evidence from blast pipeline
# - HHSearch evidence from hhsearch pipeline
# - Integrated evidence from summary service
# - Produces classified domain models for downstream use
```

This sophisticated domain analysis system represents the culmination of computational evidence into biologically meaningful domain classifications, with extensive quality control, error handling, and optimization for large-scale processing.

## Key Features

### Evidence Integration
- **Multi-source evidence** from BLAST, HHSearch, and self-comparison
- **Confidence-based weighting** of different evidence types
- **Quality filtering** with configurable thresholds
- **Coverage validation** against reference domains

### Domain Boundary Detection
- **Sophisticated algorithms** for identifying domain boundaries
- **Overlap resolution** between competing domain predictions
- **Discontinuous domain support** for complex multi-segment domains
- **Reference coverage validation** ensures biological relevance

### ECOD Classification
- **Hierarchical assignment** following ECOD taxonomy
- **Classification confidence** scoring
- **Database integration** with reference ECOD domains
- **Automatic validation** against known classifications

### Batch Processing
- **SLURM cluster integration** for high-throughput processing
- **Parallel execution** with configurable worker counts
- **Progress tracking** and status monitoring
- **Error recovery** and failed job reprocessing

## Pipeline Modes

### FULL Mode (Complete Analysis)
- BLAST searches + HHSearch analysis
- Complete evidence integration
- Highest accuracy domain identification
- Recommended for final classifications

### BLAST_ONLY Mode (Fast Processing)
- BLAST evidence only
- Faster execution (no HHSearch)
- Good for preliminary analysis
- Suitable for high-confidence cases

### ANALYSIS_ONLY Mode (Re-analysis)
- Uses existing evidence files
- No new searches performed
- Domain analysis and classification only
- Useful for parameter tuning

## Dependencies

### Core Requirements
- **Python 3.8+**
- **PostgreSQL 12+** (database backend)
- **PyYAML** (configuration management)
- **psycopg2** (PostgreSQL connectivity)

### Bioinformatics Tools
- **BLAST+** (blastp, makeblastdb)
- **HH-suite** (hhblits, hhsearch, hhmake)
- **NCBI databases** (for profile generation)

### Optional Dependencies
- **SLURM** (for cluster job submission)
- **Module system** (environment management on clusters)

### Reference Databases
- **ECOD domain database** (for BLAST searches)
- **ECOD chain database** (for chain-level searches)
- **ECOD HMM database** (for HHSearch)
- **UniClust database** (for profile generation)

## Current Status

### ✅ Completed Components
- Core configuration and database systems
- Job management framework
- Data models and serialization
- BLAST and HHSearch pipelines
- Domain analysis services
- Pipeline orchestration

### 🚧 In Development
- **CLI migration**: Moving from standalone scripts to unified CLI
- **Performance optimization**: Caching and parallel processing improvements
- **Documentation**: API documentation and user guides
- **Testing**: Comprehensive test suite

### 📋 Planned Features
- **Interactive CLI modes** for complex workflows
- **Performance monitoring** and optimization tools
- **REST API** for remote processing
- **Machine learning integration** for evidence scoring

## Module Interactions

### Application Context
The `ApplicationContext` serves as the central orchestration point:
```python
from ecod.core.context import ApplicationContext

# Initialize shared resources
context = ApplicationContext("config.yml")
# Provides: database, configuration, job_manager

# All services use the same context
blast_pipeline = BlastPipeline(context)
domain_service = DomainAnalysisService(context)
orchestrator = PipelineOrchestrator(context)
```

### Service Integration
```python
# Typical workflow integration
orchestrator = PipelineOrchestrator(context)
run = orchestrator.run_pipeline(
    batch_id=123,
    mode=PipelineMode.FULL
)

# Individual service usage
summary_service = DomainSummaryService(context)
partition_service = DomainPartitionService(context)

# Services share database, configuration, job management
```

## File Organization

### Input Files
```
{pdb_id}_{chain_id}.fasta                    # Protein sequence
```

### Intermediate Files
```
{pdb_id}_{chain_id}.chain_blast.xml          # Chain BLAST results
{pdb_id}_{chain_id}.domain_blast.xml         # Domain BLAST results
{pdb_id}_{chain_id}.{ref}.a3m                # HHblits profile
{pdb_id}_{chain_id}.{ref}.hhr                # HHsearch results
{pdb_id}_{chain_id}.{ref}.hhsearch.xml       # Processed HHsearch
```

### Output Files
```
{pdb_id}_{chain_id}.{ref}.evidence_summary.xml   # Integrated evidence
{pdb_id}_{chain_id}.{ref}.domains.xml            # Final domain models
```

## Configuration Example

```yaml
# Database connection
database:
  host: "localhost"
  port: 5432
  database: "ecod_pipeline"
  user: "ecod_user"

# Tool paths
tools:
  blast_path: "blastp"
  hhblits_path: "hhblits"
  hhsearch_path: "hhsearch"

# Reference databases
reference:
  current_version: "develop291"
  chain_db: "/data/ecod/ecod_chains.fasta"
  domain_db: "/data/ecod/ecod_domains.fasta"
  ecod_hh_db: "/data/ecod/ecod_domains_hhm"

# Pipeline settings
orchestration:
  checkpoint_enabled: true
  continue_on_error: true
  hhsearch_threads: 8
  hhsearch_memory: "16G"
```

## Development Philosophy

pyECOD emphasizes:
- **Service-oriented architecture** over monolithic scripts
- **Type-safe data models** with rich functionality
- **Comprehensive error handling** and recovery
- **Database-driven** progress tracking and results storage
- **Configurable processing** for different use cases
- **Maintainable, testable code** with clear APIs

---

**Note**: This repository is under active development and not yet released. The codebase represents a modern rewrite of the ECOD pipeline with improved architecture, error handling, and maintainability.
