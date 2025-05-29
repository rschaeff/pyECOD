# Script-to-CLI Migration: Targeted Entry Points

**Status: In Progress** - These large, feature-rich scripts represent the primary targets for CLI integration, replacing standalone script-based entry points with unified CLI commands.

## Overview

The pyECOD project includes several comprehensive standalone scripts that provide essential functionality for protein domain analysis. These scripts are being systematically integrated into the unified CLI system to provide better organization, consistency, and user experience.

## Target Scripts for CLI Migration

### 1. Domain Partition Processing (`scripts/domain_partition_run.py`)

**Size: ~2,400 lines** - The largest and most comprehensive script in the pipeline.

#### Current Script Structure
```bash
python scripts/domain_partition_run.py <mode> <action> [options]
```

**Modes and Actions:**
- `process batch` - Process entire batches
- `process all` - Process multiple batches with SLURM support
- `process specific` - Process specific proteins by process ID
- `process single` - Process individual proteins
- `analyze status` - Check batch processing status
- `analyze protein` - Analyze specific protein status  
- `analyze counts` - Generate domain count statistics
- `repair failed` - Reset and retry failed processes
- `repair missing` - Regenerate missing domain files
- `monitor batch` - Monitor batch processing progress
- `cleanup non-representative` - Clean up status tracking

#### CLI Integration Mapping
This script maps directly to the `ecod domain` command group:

| Script Command | CLI Command |
|----------------|-------------|
| `python scripts/domain_partition_run.py process batch --batch-id 123` | `ecod domain partition batch --batch-id 123` |
| `python scripts/domain_partition_run.py process all --use-slurm` | `ecod domain partition all --use-slurm` |
| `python scripts/domain_partition_run.py analyze status --batch-ids 123 124` | `ecod domain analyze status --batch-ids 123 124` |
| `python scripts/domain_partition_run.py repair failed --batch-id 123 --rerun` | `ecod domain repair failed --batch-id 123 --rerun` |
| `python scripts/domain_partition_run.py monitor batch --batch-id 123` | `ecod domain monitor --batch-id 123` |

#### Key Features Being Migrated
- **Service-Based Architecture**: Uses `DomainPartitionService` for high-level operations
- **SLURM Integration**: Full support for cluster job submission and monitoring
- **Comprehensive Analysis**: Detailed status tracking and statistics generation
- **Repair Capabilities**: Automated fixing of failed processes and missing files
- **Real-time Monitoring**: Progress tracking with configurable intervals

### 2. Domain Diagnostic Tools (`scripts/domain_partition_tools.py`)

**Size: ~800 lines** - Specialized diagnostic utilities for troubleshooting domain analysis issues.

#### Current Script Structure
```bash
python scripts/domain_partition_tools.py <mode> [options]
```

**Modes:**
- `process` - Diagnose specific process issues
- `batch` - Diagnose entire batch problems
- `failures` - Analyze failure patterns
- `readiness` - Check batch readiness for processing

#### CLI Integration Mapping
This script integrates into the `ecod domain diagnose` commands:

| Script Command | CLI Command |
|----------------|-------------|
| `python scripts/domain_partition_tools.py process --process-id 456` | `ecod domain diagnose process --process-id 456` |
| `python scripts/domain_partition_tools.py batch --batch-id 123` | `ecod domain diagnose batch --batch-id 123` |
| `python scripts/domain_partition_tools.py failures --limit 20` | `ecod domain analyze failures --limit 20` |
| `python scripts/domain_partition_tools.py readiness --batch-id 123` | `ecod domain diagnose readiness --batch-id 123` |

#### Key Features Being Migrated
- **Deep Diagnostics**: File system vs database consistency checks
- **Path Resolution**: Automatic detection and fixing of legacy file paths
- **Evidence Analysis**: Detailed examination of domain evidence availability
- **Recommendation Engine**: Automated suggestions for fixing identified issues

### 3. HHSearch Processing Tools (`scripts/hhsearch_tools.py`)

**Size: ~2,100 lines** - Comprehensive HHSearch pipeline management and analysis.

#### Current Script Structure
```bash
python scripts/hhsearch_tools.py <mode> [action] [options]
```

**Modes and Actions:**
- `process db` - Process using database backend
- `process fs` - Process using filesystem backend  
- `process register` - Register HHR files using HHResultRegistrar
- `collate` - Combine HHSearch results with BLAST evidence
- `collate-all` - Process all batches
- `parallel-collate` - SLURM-based parallel processing
- `analyze content` - Examine file content and consistency
- `analyze missing` - Find missing files
- `repair missing` - Fix missing HHR/XML files
- `repair fix_summaries` - Update domain summaries with HHSearch data
- `repair db_paths` - Fix database file path inconsistencies
- `batches` - Process multiple batches in parallel

#### CLI Integration Mapping
This script maps to the `ecod hhsearch` command group:

| Script Command | CLI Command |
|----------------|-------------|
| `python scripts/hhsearch_tools.py process db --batch-id 123` | `ecod hhsearch process db --batch-id 123` |
| `python scripts/hhsearch_tools.py collate --batch-id 123` | `ecod hhsearch collate batch --batch-id 123` |
| `python scripts/hhsearch_tools.py parallel-collate --wait` | `ecod hhsearch collate parallel --wait` |
| `python scripts/hhsearch_tools.py analyze content --batch-path /data/batch` | `ecod hhsearch analyze content --batch-path /data/batch` |
| `python scripts/hhsearch_tools.py repair missing --batch-path /data/batch` | `ecod hhsearch repair missing --batch-path /data/batch` |

#### Key Features Being Migrated
- **Multiple Backends**: Database, filesystem, and registration processing modes
- **Parallel Processing**: ThreadPoolExecutor and SLURM-based parallelization
- **Content Analysis**: Deep inspection of HHR, XML, and summary files
- **Comprehensive Repair**: Automated fixing of missing or corrupted files
- **Evidence Collation**: Integration of HHSearch results with BLAST evidence

## Migration Benefits

### 1. Unified Interface
**Before (Multiple Scripts):**
```bash
python scripts/domain_partition_run.py process batch --batch-id 123
python scripts/domain_partition_tools.py readiness --batch-id 123  
python scripts/hhsearch_tools.py collate --batch-id 123
```

**After (Unified CLI):**
```bash
ecod domain partition batch --batch-id 123
ecod domain diagnose readiness --batch-id 123
ecod hhsearch collate batch --batch-id 123
```

### 2. Consistent Option Handling
- Standardized `--config`, `--verbose`, `--log-file` options across all commands
- Uniform help system with `ecod <group> <command> --help`
- Consistent error handling and logging throughout

### 3. Improved Discoverability
- Hierarchical command structure makes functionality easier to find
- Built-in help system shows available commands and options
- Logical grouping by functional area

### 4. Better Integration
- Shared configuration management
- Unified error handling and logging
- Common database connection handling
- Standardized output formatting

## Migration Strategy

### Phase 1: Core Integration (Current)
- [x] CLI framework established with command groups
- [x] Base command infrastructure implemented
- [x] Domain and HHSearch command groups defined
- [ ] Script functionality ported to CLI commands

### Phase 2: Feature Parity
- [ ] All script modes and actions available as CLI commands
- [ ] Complete option coverage for all existing functionality
- [ ] Backward compatibility maintained during transition

### Phase 3: Enhanced Features  
- [ ] Improved help and documentation
- [ ] Better error messages and user guidance
- [ ] Enhanced progress reporting and status display
- [ ] Integration testing with existing workflows

### Phase 4: Script Deprecation
- [ ] Scripts marked as deprecated with migration guidance
- [ ] Documentation updated to use CLI commands
- [ ] Migration guide for existing users
- [ ] Scripts removed after sufficient transition period

## Usage Examples: Before and After

### Complex Domain Analysis Workflow

**Before (Multiple Scripts):**
```bash
# Check batch readiness
python scripts/domain_partition_tools.py readiness --batch-id 123 --detailed

# Run domain partition
python scripts/domain_partition_run.py process batch --batch-id 123 --workers 8

# Check for failures  
python scripts/domain_partition_tools.py failures --batch-id 123 --limit 10

# Repair failed processes
python scripts/domain_partition_run.py repair failed --batch-id 123 --rerun

# Monitor progress
python scripts/domain_partition_run.py monitor batch --batch-id 123 --interval 60
```

**After (Unified CLI):**
```bash
# Check batch readiness
ecod domain diagnose readiness --batch-id 123 --detailed

# Run domain partition  
ecod domain partition batch --batch-id 123 --workers 8

# Check for failures
ecod domain analyze failures --batch-id 123 --limit 10

# Repair failed processes
ecod domain repair failed --batch-id 123 --rerun

# Monitor progress
ecod domain monitor --batch-id 123 --interval 60
```

### HHSearch Processing Pipeline

**Before (Script):**
```bash
# Process HHR files to XML
python scripts/hhsearch_tools.py process register --batch-id 123 --force

# Analyze content for issues
python scripts/hhsearch_tools.py analyze content --batch-path /data/batch_123

# Repair missing files
python scripts/hhsearch_tools.py repair missing --batch-path /data/batch_123

# Collate with BLAST evidence
python scripts/hhsearch_tools.py collate --batch-id 123

# Process multiple batches in parallel
python scripts/hhsearch_tools.py parallel-collate --batch-ids 120 121 122 --wait
```

**After (CLI):**
```bash
# Process HHR files to XML
ecod hhsearch process register --batch-id 123 --force

# Analyze content for issues  
ecod hhsearch analyze content --batch-path /data/batch_123

# Repair missing files
ecod hhsearch repair missing --batch-path /data/batch_123

# Collate with BLAST evidence
ecod hhsearch collate batch --batch-id 123

# Process multiple batches in parallel
ecod hhsearch collate parallel --batch-ids 120 121 122 --wait
```

## Implementation Status

### Currently Available (CLI)
- ✅ Basic command structure and parsing
- ✅ Configuration management integration
- ✅ Database connection handling
- ✅ Error handling framework
- ✅ Logging integration

### In Development
- 🚧 Complete domain partition command integration
- 🚧 Full HHSearch command implementation  
- 🚧 Diagnostic tool integration
- 🚧 SLURM job management integration
- 🚧 Progress monitoring and status reporting

### Planned
- 📋 Advanced help system with examples
- 📋 Shell completion scripts
- 📋 Interactive mode for complex workflows
- 📋 Configuration validation and guidance
- 📋 Performance monitoring and optimization

## Technical Considerations

### Code Reuse Strategy
The migration prioritizes code reuse rather than rewriting:
- Existing service classes (`DomainPartitionService`, `HHSearchProcessor`) are used directly
- Core business logic remains unchanged
- CLI commands act as thin wrappers around existing functionality

### Backward Compatibility
During the transition period:
- Original scripts remain functional
- Both interfaces accept the same configuration files
- Database schema and file formats are unchanged
- Existing workflows continue to work

### Testing Strategy
- CLI commands are tested against the same test cases as original scripts
- Integration tests verify end-to-end functionality
- Performance benchmarks ensure no regression
- User acceptance testing with existing workflows

## Conclusion

The migration of these large, complex scripts to the unified CLI represents a significant improvement in the pyECOD user experience. By maintaining full feature parity while providing better organization and consistency, the CLI will make the powerful domain analysis capabilities more accessible and easier to use.

The hierarchical command structure, unified configuration management, and comprehensive help system will significantly reduce the learning curve for new users while providing existing users with a more intuitive and powerful interface to the same underlying functionality.
