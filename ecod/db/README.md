# ECOD Database Module

The ECOD Database Module provides a comprehensive database abstraction layer for the pyECOD pipeline, handling protein domain classification data storage, retrieval, and management.

## Overview

This module implements a repository pattern architecture for managing:
- **Protein data**: Structure, sequence, and metadata
- **Domain classifications**: ECOD hierarchical classifications
- **Processing workflows**: Batch processing and job tracking
- **Reference data**: ECOD versions and resources

The module is built on PostgreSQL and provides type-safe, error-handled database operations with automatic connection management and transaction support.

## Architecture

```
ecod/db/
├── __init__.py                 # Module exports
├── manager.py                  # Core database manager
├── migration_manager.py        # Schema migration system
├── data_migration.py          # Data migration utilities
├── utils.py                   # Database utilities
├── migrations/                # SQL migration files
│   └── 001_create_initial_schema.sql
└── repositories/              # Data access repositories
    ├── domain_repository.py   # Domain data operations
    └── protein_repository.py  # Protein data operations
```

### Core Components

1. **DBManager**: Central database connection and query management
2. **Repository Classes**: Type-safe data access objects
3. **Migration System**: Schema versioning and updates
4. **Error Handling**: Comprehensive exception management

## Quick Start

### 1. Basic Setup

```python
from ecod.db import DBManager
from ecod.db.repositories.protein_repository import ProteinRepository
from ecod.db.repositories.domain_repository import DomainRepository

# Database configuration
config = {
    'host': 'localhost',
    'port': 5432,
    'database': 'ecod_pipeline',
    'user': 'ecod_user',
    'password': 'your_password'
}

# Initialize database manager
db_manager = DBManager(config)

# Initialize repositories
protein_repo = ProteinRepository(db_manager)
domain_repo = DomainRepository(db_manager)
```

### 2. Schema Migrations

```python
from ecod.db.migration_manager import MigrationManager

# Run schema migrations
migration_manager = MigrationManager(config, "ecod/db/migrations")
migration_manager.apply_migrations()
```

### 3. Basic Operations

```python
# Get a protein by PDB and chain ID
protein = protein_repo.get_by_pdb_chain("1abc", "A")

# Search for domains by classification
domains = domain_repo.get_by_classification("t_group", "2004.1.1")

# Create a new protein
from ecod.models.protein import Protein, ProteinSequence

protein = Protein(
    pdb_id="1xyz",
    chain_id="A",
    source_id="1xyz_A",
    length=250
)
protein_id = protein_repo.create(protein)
```

## Database Schema

### Core Tables

#### Protein Data
- `ecod_schema.protein`: Protein chain metadata
- `ecod_schema.protein_sequence`: Amino acid sequences
- `ecod_schema.protein_structure`: Structure metadata (resolution, method)

#### Domain Classification
- `ecod_schema.domain`: Domain boundaries and classifications
- `ecod_schema.domain_sequence`: Domain-specific sequences
- `ecod_schema.domain_dssp_detail`: Secondary structure information

#### Processing Workflow
- `ecod_schema.batch`: Processing batch metadata
- `ecod_schema.process_status`: Per-protein processing status
- `ecod_schema.process_file`: File tracking and validation
- `ecod_schema.job`: SLURM job tracking
- `ecod_schema.job_item`: Individual job items

#### Reference Data
- `ecod_schema.ecod_version`: ECOD database versions
- `ecod_schema.reference_resource`: Reference file locations

### Key Relationships

```sql
protein 1:1 protein_sequence
protein 1:n domain
domain 1:1 domain_sequence
domain 1:1 domain_dssp_detail
batch 1:n process_status
process_status 1:n process_file
```

## Repository Pattern

### Protein Repository

```python
# Retrieve operations
protein = protein_repo.get_by_id(123)
protein = protein_repo.get_by_source_id("1abc_A")
proteins = protein_repo.search({"length_min": 100, "length_max": 500})

# CRUD operations
protein_id = protein_repo.create(protein)
success = protein_repo.update(protein)
success = protein_repo.delete(protein_id)

# Specialized queries
unprocessed = protein_repo.get_unprocessed_proteins(limit=100)
count = protein_repo.count_proteins({"pdb_id": "1abc"})
```

### Domain Repository

```python
# Classification-based queries
domains = domain_repo.get_by_classification("h_group", "2004.1")
domains = domain_repo.get_by_protein_id(123)

# Search with criteria
domains = domain_repo.search({
    "length_min": 50,
    "representative": True,
    "has_structure": True
})

# Evidence data
sequence = domain_repo.get_sequence(domain_id)
dssp = domain_repo.get_dssp_details(domain_id)
```

## Error Handling

The module provides structured exception handling:

```python
from ecod.exceptions import ConnectionError, QueryError, DatabaseError, ValidationError

try:
    protein = protein_repo.get_by_id(123)
except ConnectionError as e:
    # Database connection failed
    logger.error(f"Connection failed: {e.message}")
    
except QueryError as e:
    # SQL query error
    logger.error(f"Query failed: {e.message}")
    logger.debug(f"Query: {e.context['query']}")
    
except ValidationError as e:
    # Data validation error
    logger.error(f"Validation failed: {e.message}")
```

## Advanced Usage

### Transaction Management

```python
def transfer_domain_classification(domain_id, new_classification):
    def transaction_callback(cursor):
        # Update domain classification
        cursor.execute(
            "UPDATE ecod_schema.domain SET t_group = %s WHERE id = %s",
            (new_classification, domain_id)
        )
        
        # Log the change
        cursor.execute(
            "INSERT INTO ecod_schema.classification_history (domain_id, old_value, new_value) VALUES (%s, %s, %s)",
            (domain_id, old_classification, new_classification)
        )
        
        return cursor.rowcount
    
    return db_manager.execute_transaction(transaction_callback)
```

### Custom Queries

```python
# Complex analytical queries
query = """
SELECT 
    t_group,
    COUNT(*) as domain_count,
    AVG(length) as avg_length
FROM ecod_schema.domain 
WHERE is_manual_rep = TRUE
GROUP BY t_group
ORDER BY domain_count DESC
"""

results = db_manager.execute_dict_query(query)
```

### Batch Operations

```python
# Process large datasets efficiently
def process_protein_batch(batch_size=1000):
    offset = 0
    while True:
        proteins = protein_repo.get_all(limit=batch_size, offset=offset)
        if not proteins:
            break
            
        # Process batch
        for protein in proteins:
            # Perform analysis
            pass
            
        offset += batch_size
```

## Data Migration

### Schema Migrations

Create migration files in `ecod/db/migrations/`:

```sql
-- 002_add_new_classification_level.sql
ALTER TABLE ecod_schema.domain 
ADD COLUMN f_group VARCHAR(50);

CREATE INDEX idx_domain_f_group 
ON ecod_schema.domain (f_group);
```

Apply migrations:

```python
migration_manager = MigrationManager(config, "ecod/db/migrations")
migration_manager.apply_migrations()
```

### Data Migration

```python
from ecod.db.data_migration import DataMigration

# Migrate from old schema to new schema
migration = DataMigration(old_config, new_config)
migration.migrate_protein_data()
migration.migrate_reference_data()
```

## Configuration

### Database Configuration

```python
config = {
    'host': 'localhost',
    'port': 5432,
    'database': 'ecod_pipeline',
    'user': 'ecod_user',
    'password': 'secure_password',
    'connect_timeout': 10,
    'application_name': 'pyECOD'
}
```

### Connection Pooling

For high-throughput applications, consider connection pooling:

```python
import psycopg2.pool

pool = psycopg2.pool.ThreadedConnectionPool(
    minconn=1,
    maxconn=20,
    **config
)

# Custom DBManager with pooling
class PooledDBManager(DBManager):
    def __init__(self, pool):
        self.pool = pool
        super().__init__(config)
    
    @contextmanager
    def get_connection(self):
        conn = self.pool.getconn()
        try:
            yield conn
            conn.commit()
        except Exception:
            conn.rollback()
            raise
        finally:
            self.pool.putconn(conn)
```

## Performance Optimization

### Query Optimization

1. **Use appropriate indexes**: The schema includes indexes on frequently queried columns
2. **Batch operations**: Use batch inserts for large datasets
3. **Limit result sets**: Always use appropriate LIMIT clauses
4. **Use prepared statements**: The module automatically uses parameterized queries

### Database Tuning

```sql
-- Recommended PostgreSQL settings for ECOD pipeline
-- In postgresql.conf:

shared_buffers = 256MB                # 25% of RAM
effective_cache_size = 1GB           # 75% of RAM
work_mem = 64MB                      # Per-connection work memory
maintenance_work_mem = 256MB         # For maintenance operations
checkpoint_completion_target = 0.9   # Smooth checkpoints
wal_buffers = 16MB                   # WAL buffer size
```

## Monitoring and Logging

### Database Logging

```python
import logging

# Configure database logging
logging.getLogger('ecod.db').setLevel(logging.INFO)
logging.getLogger('ecod.db.manager').setLevel(logging.DEBUG)
```

### Query Performance Monitoring

```python
# Enable query timing in development
import time

class TimedDBManager(DBManager):
    def execute_query(self, query, params=None):
        start_time = time.time()
        result = super().execute_query(query, params)
        elapsed = time.time() - start_time
        
        if elapsed > 1.0:  # Log slow queries
            self.logger.warning(f"Slow query ({elapsed:.2f}s): {query[:100]}...")
            
        return result
```

## Testing

### Unit Tests

```python
import pytest
from ecod.db import DBManager
from ecod.db.repositories.protein_repository import ProteinRepository

@pytest.fixture
def db_manager():
    test_config = {
        'host': 'localhost',
        'database': 'ecod_test',
        'user': 'test_user',
        'password': 'test_password'
    }
    return DBManager(test_config)

def test_protein_creation(db_manager):
    repo = ProteinRepository(db_manager)
    
    protein = Protein(
        pdb_id="TEST",
        chain_id="A",
        source_id="TEST_A",
        length=100
    )
    
    protein_id = repo.create(protein)
    assert protein_id is not None
    
    retrieved = repo.get_by_id(protein_id)
    assert retrieved.source_id == "TEST_A"
```

### Integration Tests

```python
def test_domain_classification_workflow():
    # Test complete workflow from protein creation to domain classification
    pass
```

## Best Practices

### Error Handling
- Always catch specific exceptions
- Log errors with appropriate context
- Use transactions for multi-table operations

### Performance
- Use appropriate query limits
- Batch large operations
- Monitor slow queries

### Security
- Use parameterized queries (automatically handled)
- Implement proper access controls
- Regularly update database credentials

### Maintainability
- Use repository pattern for data access
- Write comprehensive tests
- Document complex queries
- Use meaningful variable names

## Contributing

When adding new functionality:

1. **Follow the repository pattern**: Create repository classes for new entity types
2. **Add migrations**: Create SQL migration files for schema changes
3. **Write tests**: Include unit and integration tests
4. **Update documentation**: Document new methods and usage patterns
5. **Handle errors**: Implement appropriate exception handling

## Troubleshooting

### Common Issues

1. **Connection Errors**
   ```python
   # Check database connectivity
   if not db_manager.schema_exists('ecod_schema'):
       print("Schema not found - run migrations")
   ```

2. **Performance Issues**
   ```sql
   -- Check for missing indexes
   SELECT schemaname, tablename, attname, n_distinct, correlation
   FROM pg_stats
   WHERE schemaname = 'ecod_schema'
   ORDER BY n_distinct DESC;
   ```

3. **Migration Issues**
   ```python
   # Check migration status
   migration_manager = MigrationManager(config, "ecod/db/migrations")
   applied = migration_manager._get_applied_migrations(conn)
   print(f"Applied migrations: {applied}")
   ```

For additional support, check the logs or contact the development team.
