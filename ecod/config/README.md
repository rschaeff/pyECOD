# pyECOD Configuration System Documentation

## Overview

The pyECOD pipeline uses a flexible, hierarchical configuration system that supports multiple configuration sources with automatic merging and validation. The system is designed to handle both development and production environments while keeping sensitive information secure.

## Configuration Sources (Priority Order)

The configuration system loads settings from multiple sources in the following priority order (later sources override earlier ones):

1. **Default Configuration** - Built-in defaults
2. **Main Configuration File** - `config.yml` or specified file
3. **Local Configuration File** - `config.local.yml` (for secrets and local overrides)
4. **Environment Variables** - Runtime overrides with `ECOD_` prefix

## Configuration File Structure

### Main Configuration File (`config.yml`)

```yaml
# Database connection settings
database:
  host: "dione"
  port: 45000
  database: "ecod_protein"
  user: "ecod"
  # password: DO NOT PUT PASSWORDS HERE - use config.local.yml

# File paths and directories
paths:
  output_dir: "./output"
  blast_db: "/data/blast/nr"

# Tool executable paths
tools:
  blast_path: "blastp"
  hhblits_path: "hhblits"
  hhsearch_path: "hhsearch"
  hhmake_path: "hhmake"

# Reference database settings
reference:
  current_version: "develop291"
  uniclust_db: "/data/uniclust/UniRef30_2020_06"
  chain_db: "/data/ecod/chain_blast_db"
  domain_db: "/data/ecod/domain_blast_db"

# Logging configuration
logging:
  level: "INFO"
  format: "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
  log_dir: "./logs"

# Pipeline execution settings
pipeline:
  batch_size: 100
  threads: 8
  memory: "8G"
  timeout: 3600
  check_interval: 300
  force_overwrite: false
```

### Local Configuration File (`config.local.yml`)

**⚠️ IMPORTANT: This file contains secrets and should NEVER be committed to version control!**

```yaml
# Local configuration file - contains secrets and local overrides
# This file should be in .gitignore and not committed to version control

database:
  password: "your_secure_database_password"
  host: "localhost"  # Override for local development

paths:
  output_dir: "/home/username/ecod_output"  # Local output directory

# Local development overrides
pipeline:
  threads: 4  # Fewer threads for local development
  batch_size: 10  # Smaller batches for testing
```

## Environment Variable Overrides

You can override any configuration value using environment variables with the `ECOD_` prefix:

### Simple Values
```bash
export ECOD_DATABASE="my_database"
export ECOD_THREADS=16
```

### Nested Values (use double underscore)
```bash
export ECOD_DATABASE__HOST="production-db-server"
export ECOD_DATABASE__PASSWORD="secure_password"
export ECOD_PIPELINE__BATCH_SIZE=500
export ECOD_TOOLS__BLAST_PATH="/usr/local/bin/blastp"
```

### Type Conversion
Environment variables are automatically converted to appropriate types:
- `"true"`, `"yes"`, `"1"` → `True`
- `"false"`, `"no"`, `"0"` → `False`
- Numeric strings → `int` or `float`
- Everything else → `str`

## API Usage

### Basic Usage

```python
from ecod.config import ConfigManager

# Load configuration (searches for config.yml in current directory)
config = ConfigManager()

# Load from specific file
config = ConfigManager("path/to/config.yml")

# Access configuration values
db_host = config.get("database.host")
output_dir = config.get("paths.output_dir")
threads = config.get("pipeline.threads", default=8)
```

### Specialized Getters

```python
# Get database configuration as dictionary
db_config = config.get_db_config()
# Returns: {'host': 'dione', 'port': 45000, 'database': 'ecod_protein', ...}

# Get file paths
output_dir = config.get_path("output_dir")
blast_db = config.get_path("blast_db", default="/default/blast/db")

# Get tool paths
blast_path = config.get_tool_path("blast")  # Gets tools.blast_path
hhblits_path = config.get_tool_path("hhblits")

# Get reference paths
version = config.get_reference_path("current_version")
uniclust = config.get_reference_path("uniclust_db")

# Check feature flags
if config.is_enabled("experimental_feature"):
    # Feature is enabled
    pass
```

### Database Connection Example

```python
import psycopg2
from ecod.config import ConfigManager

config = ConfigManager("config/config.yml")
db_config = config.get_db_config()

# Connect to database using configuration
conn = psycopg2.connect(**db_config)
```

## Configuration Schema Validation

The configuration system includes schema validation to catch configuration errors early:

### Validation Rules

- **Required Fields**: Some fields are mandatory (e.g., database connection info)
- **Type Checking**: Values must match expected types (int, str, bool)
- **Section Validation**: Required sections must be present

### Handling Validation Errors

```python
config = ConfigManager("config.yml")
# Validation happens automatically during initialization
# Errors are logged, but configuration still loads with warnings
```

Common validation errors:
- Missing required database connection parameters
- Invalid data types (e.g., string where integer expected)
- Missing required configuration sections

## Security Best Practices

### 1. Use Local Configuration for Secrets

**DO:**
```yaml
# config.local.yml (not in version control)
database:
  password: "actual_password"
```

**DON'T:**
```yaml
# config.yml (in version control)
database:
  password: "secret123"  # ❌ NEVER DO THIS
```

### 2. Environment Variables for Production

For production deployments, use environment variables for sensitive data:

```bash
# Production environment
export ECOD_DATABASE__PASSWORD="$DB_PASSWORD"
export ECOD_DATABASE__HOST="$DB_HOST"
```

### 3. .gitignore Configuration

Ensure your `.gitignore` includes:
```gitignore
# Local configuration files (contain secrets)
config.local.yml
config.local.yaml
*.local.yml
*.local.yaml

# Environment files
.env
.env.local
```

## File Structure and Location

```
project_root/
├── config/
│   ├── config.yml          # Main configuration (versioned)
│   ├── config.local.yml    # Local secrets (NOT versioned)
│   ├── production.yml      # Production overrides
│   └── development.yml     # Development settings
├── ecod/
│   └── config/
│       ├── __init__.py
│       ├── manager.py      # ConfigManager class
│       ├── schema.py       # Validation schema
│       └── defaults.py     # Default values
└── .gitignore              # Must exclude *.local.yml
```

## Configuration Examples

### Development Environment

```yaml
# config.local.yml (development)
database:
  host: "localhost"
  password: "dev_password"

paths:
  output_dir: "/tmp/ecod_dev"

pipeline:
  batch_size: 10
  threads: 2
```

### Production Environment

```bash
# Production environment variables
export ECOD_DATABASE__HOST="prod-db.example.com"
export ECOD_DATABASE__PASSWORD="secure_prod_password"
export ECOD_PATHS__OUTPUT_DIR="/data/ecod/production"
export ECOD_PIPELINE__BATCH_SIZE=1000
export ECOD_PIPELINE__THREADS=32
```

### Testing Environment

```yaml
# config.test.yml
database:
  database: "ecod_test"

paths:
  output_dir: "/tmp/ecod_test"

pipeline:
  batch_size: 5
  timeout: 60
```

## Common Patterns

### Pipeline Scripts

```python
#!/usr/bin/env python3
import sys
from ecod.config import ConfigManager

def main():
    config_path = sys.argv[1] if len(sys.argv) > 1 else "config/config.yml"
    config = ConfigManager(config_path)
    
    # Use configuration
    batch_size = config.get("pipeline.batch_size")
    threads = config.get("pipeline.threads")
    
    # Process data...

if __name__ == "__main__":
    main()
```

### Database Utilities

```python
from ecod.config import ConfigManager
import psycopg2

class DatabaseManager:
    def __init__(self, config_path=None):
        self.config = ConfigManager(config_path)
        self.db_config = self.config.get_db_config()
    
    def get_connection(self):
        return psycopg2.connect(**self.db_config)
```

## Troubleshooting

### Common Issues

1. **Missing Local Config**: If secrets aren't loading, check that `config.local.yml` exists and has correct format

2. **Environment Variables Not Working**: 
   - Ensure proper `ECOD_` prefix
   - Use double underscore `__` for nested values
   - Check variable names match configuration structure

3. **Validation Errors**: 
   - Check required fields are present
   - Verify data types match schema expectations
   - Review error messages in logs

4. **File Not Found**: 
   - Verify configuration file path is correct
   - Check file permissions
   - Ensure working directory is correct

### Debug Configuration Loading

```python
import logging
logging.basicConfig(level=logging.DEBUG)

from ecod.config import ConfigManager
config = ConfigManager("config.yml")  # Will show debug info
```

## Default Configuration Reference

The system includes these built-in defaults:

```python
DEFAULT_CONFIG = {
    'database': {
        'database': 'ecod_protein',
        'host': 'dione',
        'port': 45000,
        'user': 'ecod',
    },
    'paths': {
        'output_dir': './output',
    },
    'tools': {
        'blast_path': 'blastp',
        'hhblits_path': 'hhblits',
        'hhsearch_path': 'hhsearch',
        'hhmake_path': 'hhmake',
    },
    'reference': {
        'current_version': 'develop291',
    },
    'logging': {
        'level': 'INFO',
        'format': '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    },
    'pipeline': {
        'batch_size': 100,
        'threads': 8,
        'memory': '8G',
        'timeout': 3600,
        'check_interval': 300,
        'force_overwrite': False
    }
}
```

These defaults are always loaded first and can be overridden by any other configuration source.
