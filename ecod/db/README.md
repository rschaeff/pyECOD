# ECOD Database Module

The ECOD Database Module provides database connectivity and query execution for the pyECOD pipeline, handling protein domain classification data storage and retrieval using raw SQL queries.

## Overview

This module provides a PostgreSQL database abstraction layer built around:
- **DBManager**: Core database connection and query execution
- **Raw SQL Queries**: Direct SQL for complex operations and performance
- **Schema Awareness**: Validation and safety checking for queries
- **Migration System**: Schema versioning and updates
- **Error Handling**: Comprehensive exception management

## Architecture

```
ecod/db/
├── __init__.py                 # Module exports (DBManager)
├── manager.py                  # Core database manager with schema awareness
├── migration_manager.py        # Schema migration system
├── data_migration.py          # Data migration utilities
├── utils.py                   # Database utilities
└── migrations/                # SQL migration files
    └── 001_create_initial_schema.sql
```

## Quick Start

### 1. Basic Setup

```python
from ecod.db import DBManager

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
```

### 2. Basic Query Operations

```python
# Execute queries that return data
proteins = db_manager.execute_dict_query("""
    SELECT id, pdb_id, chain_id, source_id, length
    FROM ecod_schema.protein
    WHERE length > %s
""", (100,))

# Execute queries that don't return data (INSERT, UPDATE, DELETE)
db_manager.execute_query("""
    UPDATE ecod_schema.process_status
    SET status = 'completed'
    WHERE batch_id = %s
""", (batch_id,))

# Insert with returning value
protein_id = db_manager.insert(
    "ecod_schema.protein",
    {
        "pdb_id": "1abc",
        "chain_id": "A",
        "source_id": "1abc_A",
        "length": 250
    },
    returning="id"
)
```

### 3. Schema Migrations

```python
from ecod.db.migration_manager import MigrationManager

# Run schema migrations
migration_manager = MigrationManager(config, "ecod/db/migrations")
migration_manager.apply_migrations()
```

## Database Schema

### Core Tables

#### Protein Data
- `ecod_schema.protein`: Protein chain metadata
- `ecod_schema.protein_sequence`: Amino acid sequences

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
batch 1:n process_status
process_status 1:n process_file
job 1:n job_item
```

## Raw SQL Query Patterns

The pyECOD codebase uses raw SQL queries for flexibility and performance. Here are common patterns:

### Data Retrieval

```python
# Get proteins for a batch
def get_batch_proteins(db, batch_id):
    return db.execute_dict_query("""
        SELECT p.id, p.pdb_id, p.chain_id, p.source_id, ps.status
        FROM ecod_schema.protein p
        JOIN ecod_schema.process_status ps ON p.id = ps.protein_id
        WHERE ps.batch_id = %s
        ORDER BY p.pdb_id, p.chain_id
    """, (batch_id,))

# Get processing status summary
def get_batch_status_summary(db, batch_id):
    return db.execute_dict_query("""
        SELECT
            current_stage,
            status,
            COUNT(*) as count
        FROM ecod_schema.process_status
        WHERE batch_id = %s
        GROUP BY current_stage, status
        ORDER BY current_stage, status
    """, (batch_id,))
```

### Data Updates

```python
# Update processing status
def update_process_status(db, process_id, stage, status, error_msg=None):
    return db.update(
        "ecod_schema.process_status",
        {
            "current_stage": stage,
            "status": status,
            "error_message": error_msg,
            "updated_at": "CURRENT_TIMESTAMP"
        },
        "id = %s",
        (process_id,)
    )

# Bulk status update
def mark_batch_complete(db, batch_id):
    db.execute_query("""
        UPDATE ecod_schema.batch
        SET status = 'completed', completed_at = CURRENT_TIMESTAMP
        WHERE id = %s
    """, (batch_id,))
```

### Complex Analytics

```python
# Domain classification statistics
def get_classification_stats(db):
    return db.execute_dict_query("""
        SELECT
            t_group,
            COUNT(*) as domain_count,
            COUNT(DISTINCT protein_id) as protein_count,
            AVG(length) as avg_length
        FROM ecod_schema.domain
        WHERE t_group IS NOT NULL
        GROUP BY t_group
        ORDER BY domain_count DESC
    """)
```

## Schema Awareness Features

The enhanced DBManager includes schema awareness to prevent runtime errors:

### Column Validation

```python
# Check which columns exist
validation = db_manager.validate_columns_exist('ecod_schema', 'protein',
                                              ['id', 'pdb_id', 'fake_column'])
# Returns: {'id': True, 'pdb_id': True, 'fake_column': False}

# Get only existing columns
safe_columns = db_manager.get_safe_columns('ecod_schema', 'protein',
                                          ['id', 'pdb_id', 'fake_column'])
# Returns: ['id', 'pdb_id']
```

### Safe Query Building

```python
# Build queries with validated columns
query = db_manager.build_safe_select(
    'ecod_schema', 'protein',
    ['id', 'pdb_id', 'source_id'],
    where_clause='length > %s',
    limit=100
)

if query:
    results = db_manager.execute_dict_query(query, (200,))
```

### Query Safety Validation

```python
# Validate query before execution
safety = db_manager.validate_query_safety("""
    SELECT id, pdb_id FROM ecod_schema.protein
    JOIN ecod_schema.nonexistent_table n ON p.id = n.protein_id
""")

if not safety['safe']:
    print(f"Query issues: {safety['issues']}")
    # ['Table ecod_schema.nonexistent_table does not exist']
```

### Cache Management

```python
# Refresh schema cache
db_manager.refresh_column_cache('ecod_schema')          # Refresh schema
db_manager.refresh_column_cache('ecod_schema', 'protein')  # Refresh table
db_manager.refresh_column_cache()                       # Refresh all common schemas
```

## Error Handling

The module provides structured exception handling:

```python
from ecod.exceptions import ConnectionError, QueryError, DatabaseError

try:
    results = db_manager.execute_dict_query(query, params)
except ConnectionError as e:
    logger.error(f"Database connection failed: {e.message}")

except QueryError as e:
    logger.error(f"Query failed: {e.message}")
    logger.debug(f"Query: {e.context['query']}")

except DatabaseError as e:
    logger.error(f"Database error: {e.message}")
```

## Advanced Usage

### Transaction Management

```python
def transfer_classification(db, domain_id, old_group, new_group):
    def transaction_callback(cursor):
        # Update domain classification
        cursor.execute("""
            UPDATE ecod_schema.domain
            SET t_group = %s
            WHERE id = %s AND t_group = %s
        """, (new_group, domain_id, old_group))

        if cursor.rowcount != 1:
            raise ValueError("Domain not found or classification mismatch")

        # Log the change
        cursor.execute("""
            INSERT INTO ecod_schema.classification_history
            (domain_id, old_value, new_value, changed_at)
            VALUES (%s, %s, %s, CURRENT_TIMESTAMP)
        """, (domain_id, old_group, new_group))

        return cursor.rowcount

    return db.execute_transaction(transaction_callback)
```

### Batch Processing

```python
def process_large_dataset(db, batch_size=1000):
    """Process large datasets efficiently"""
    offset = 0

    while True:
        proteins = db.execute_dict_query("""
            SELECT id, pdb_id, chain_id, source_id
            FROM ecod_schema.protein
            ORDER BY id
            LIMIT %s OFFSET %s
        """, (batch_size, offset))

        if not proteins:
            break

        # Process batch
        for protein in proteins:
            process_protein(protein)

        offset += batch_size
```

### Schema Validation in Functions

```python
def get_proteins_safe(db, requested_columns, min_length=100):
    """Get proteins with column validation"""

    # Validate columns exist
    safe_columns = db.get_safe_columns('ecod_schema', 'protein', requested_columns)
    if not safe_columns:
        return []

    # Build and execute query
    query = f"""
        SELECT {', '.join(safe_columns)}
        FROM ecod_schema.protein
        WHERE length >= %s
    """

    return db.execute_dict_query(query, (min_length,))
```

## Data Migration

### Schema Migrations

Create migration files in `ecod/db/migrations/`:

```sql
-- 002_add_domain_classification.sql
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

### Performance Configuration

```python
# For high-throughput applications
import psycopg2.pool

pool = psycopg2.pool.ThreadedConnectionPool(
    minconn=1,
    maxconn=20,
    **config
)

# Use connection pooling for concurrent processing
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

## CLI Tools

### Schema Checker

Create `scripts/check_schema.py`:

```bash
# Check table columns
python scripts/check_schema.py columns ecod_schema protein

# Validate specific columns
python scripts/check_schema.py columns ecod_schema protein --check id pdb_id fake_column

# List tables in schema
python scripts/check_schema.py tables ecod_schema

# Refresh schema cache
python scripts/check_schema.py refresh --schema ecod_schema
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
# Monitor slow queries
import time

def execute_with_timing(db, query, params=None):
    start_time = time.time()
    result = db.execute_dict_query(query, params)
    elapsed = time.time() - start_time

    if elapsed > 1.0:
        logger.warning(f"Slow query ({elapsed:.2f}s): {query[:100]}...")

    return result
```

## Testing

### Unit Tests

```python
import pytest
from ecod.db import DBManager

@pytest.fixture
def db_manager():
    test_config = {
        'host': 'localhost',
        'database': 'ecod_test',
        'user': 'test_user',
        'password': 'test_password'
    }
    return DBManager(test_config)

def test_protein_query(db_manager):
    # Test basic query
    proteins = db_manager.execute_dict_query("""
        SELECT id, pdb_id, chain_id
        FROM ecod_schema.protein
        LIMIT 5
    """)

    assert isinstance(proteins, list)

def test_schema_validation(db_manager):
    # Test column validation
    validation = db_manager.validate_columns_exist(
        'ecod_schema', 'protein',
        ['id', 'pdb_id', 'fake_column']
    )

    assert validation['id'] == True
    assert validation['pdb_id'] == True
    assert validation['fake_column'] == False
```

## Best Practices

### Query Safety
- Use parameterized queries (automatically handled by DBManager)
- Validate table/column existence before executing complex queries
- Use schema awareness features for dynamic query building

### Performance
- Use appropriate query limits
- Batch large operations
- Monitor slow queries
- Use indexes effectively

### Error Handling
- Always catch specific exceptions
- Log errors with appropriate context
- Use transactions for multi-table operations

### Maintainability
- Use meaningful variable names in queries
- Comment complex SQL
- Write tests for critical database operations
- Document query patterns

## Troubleshooting

### Common Issues

1. **Connection Errors**
   ```python
   # Check database connectivity
   if not db_manager.schema_exists('ecod_schema'):
       print("Schema not found - run migrations")
   ```

2. **Missing Tables/Columns**
   ```python
   # Validate before executing
   safety = db_manager.validate_query_safety(query)
   if not safety['safe']:
       print(f"Issues: {safety['issues']}")
   ```

3. **Performance Issues**
   ```sql
   -- Check for missing indexes
   SELECT schemaname, tablename, attname, n_distinct, correlation
   FROM pg_stats
   WHERE schemaname = 'ecod_schema'
   ORDER BY n_distinct DESC;
   ```

4. **Migration Issues**
   ```python
   # Check migration status
   migration_manager = MigrationManager(config, "ecod/db/migrations")
   # Check logs for specific migration errors
   ```

For additional support, check the logs or use the schema validation tools to diagnose issues.
