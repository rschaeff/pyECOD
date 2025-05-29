# pyECOD Core Modules

This directory contains the core infrastructure for the pyECOD pipeline. **The primary and actively maintained module is `context.py`** - all other modules are marked for deprecation or incorporation into other packages.

## 🎯 **PRIMARY MODULE: ApplicationContext**

### `context.py` ⭐ **[ACTIVELY USED]**

The **ApplicationContext** is the central orchestration point for the pyECOD application, providing initialized access to all shared resources.

**Key Features:**
- Initializes database connections, configuration, and job managers
- Provides centralized access to shared resources
- Handles force overwrite settings for pipeline operations
- Thread-safe resource management

**Basic Usage:**
```python
from ecod.core.context import ApplicationContext

# Initialize with default config
context = ApplicationContext()

# Initialize with specific config file
context = ApplicationContext("path/to/config.yml")

# Access shared resources
db = context.db
config = context.config
job_manager = context.job_manager

# Control file overwrite behavior
context.set_force_overwrite(True)
if context.is_force_overwrite():
    # Proceed with overwriting existing files
    pass
```

**Typical Integration Pattern:**
```python
class MyPipeline:
    def __init__(self, context=None):
        self.context = context or ApplicationContext()
        self.db = self.context.db
        self.config = self.context.config
        self.logger = logging.getLogger("my.pipeline")
```

---

## 🚨 **DEPRECATED / TRANSITIONAL MODULES**

The following modules are **marked for deprecation or incorporation into other packages**. Use ApplicationContext instead for new development.

### `base_pipeline.py` ➜ *Moving to `ecod.pipelines`*
- Pipeline base classes being moved to the pipelines package
- Use ApplicationContext directly in new pipeline components

### `batch_logger.py` ➜ *Moving to `ecod.logging`*  
- Structured logging functionality being consolidated
- Use standard Python logging with ApplicationContext

### `logging_config.py` ➜ *Moving to `ecod.logging`*
- Logging configuration being centralized elsewhere
- ApplicationContext handles logging setup

### `cli_utils.py` ➜ *Moving to `ecod.cli`*
- CLI utilities being moved to dedicated CLI package
- Use with ApplicationContext for CLI applications

### `command_utils.py` ➜ *Moving to `ecod.utils`*
- Command execution utilities being moved to utils
- Will integrate with ApplicationContext job management

### `validation.py` ➜ *Moving to `ecod.utils`*
- Data validation functions being reorganized
- Will be accessible through ApplicationContext

---

## 💡 **RECOMMENDED USAGE PATTERNS**

### Standard Pipeline Component Pattern
```python
from ecod.core.context import ApplicationContext
import logging

class DomainAnalysisPipeline:
    def __init__(self, context=None):
        self.context = context or ApplicationContext()
        self.db = self.context.db
        self.config = self.context.config
        self.job_manager = self.context.job_manager
        self.logger = logging.getLogger("ecod.domain_analysis")
        
    def process_batch(self, batch_id: int):
        # Use context resources directly
        chains = self.db.execute_query(
            "SELECT * FROM chains WHERE batch_id = %s", 
            (batch_id,)
        )
        # Process chains...
```

### CLI Application Pattern  
```python
from ecod.core.context import ApplicationContext

def main():
    # Initialize application context - this is all you need
    context = ApplicationContext("config.yml")
    
    # Everything else comes from context
    db = context.db
    config = context.config
    job_manager = context.job_manager
    
    # Control file overwrite behavior
    context.set_force_overwrite(True)
    
    # Process data using context resources
    results = process_data(context)
    return results

if __name__ == "__main__":
    main()
```

### Simple Service Pattern
```python
from ecod.core.context import ApplicationContext

class ProteinService:
    def __init__(self, context=None):
        self.context = context or ApplicationContext()
        
    def get_protein_domains(self, pdb_id: str, chain_id: str):
        return self.context.db.execute_query("""
            SELECT * FROM domains 
            WHERE pdb_id = %s AND chain_id = %s
        """, (pdb_id, chain_id))
```

---

## 🔗 **Dependencies**

ApplicationContext integrates with:
- **Database**: PostgreSQL via `ecod.db.DBManager`
- **Configuration**: YAML configuration via `ecod.config.ConfigManager`  
- **Job Management**: SLURM integration via `ecod.jobs`
- **Exceptions**: Custom exceptions via `ecod.exceptions`

---

## 📋 **Migration Notes**

- **For new development**: Use `ApplicationContext` directly
- **For existing code**: Gradually migrate away from deprecated modules
- **Database access**: Always go through `context.db`
- **Configuration**: Always use `context.config`
- **Job management**: Always use `context.job_manager`

---

## ⚠️ **Important**

**DO NOT** use the deprecated modules for new development. The ApplicationContext provides all necessary functionality and will remain the stable API for pyECOD core services. Other functionality is being reorganized into more appropriate packages (`ecod.pipelines`, `ecod.cli`, `ecod.utils`, `ecod.logging`).catch errors before processing
