#!/bin/bash
# Database setup script for domain partition tests

set -e

# PostgreSQL paths - adjust if needed
PG_BIN_PATH="/sw/apps/postgresql-17.4/bin"
if [ -d "$PG_BIN_PATH" ]; then
    export PATH="$PG_BIN_PATH:$PATH"
fi

# Database configuration - adjust as needed
DB_NAME="ecod_test"
DB_USER="test_user"
DB_PASS="test_pass"
DB_HOST="${PGHOST:-localhost}"     # Use PGHOST env var or default to lotta
DB_PORT="${PGPORT:-5432}"      # Use PGPORT env var or default to 5432
DB_SUPERUSER="${PGUSER:-rschaeff}"  # Use PGUSER env var or default to rschaeff

echo "Setting up test database on $DB_HOST:$DB_PORT as user $DB_SUPERUSER..."

# Check if PostgreSQL utilities are available
if ! command -v createdb &> /dev/null; then
    echo "Error: PostgreSQL utilities not found in PATH"
    echo "Please add PostgreSQL bin directory to PATH:"
    echo "export PATH="/sw/apps/postgresql-17.4/bin:\$PATH""
    exit 1
fi

# Test connection first
echo "Testing PostgreSQL connection..."
if ! pg_isready -p $DB_PORT -h $DB_HOST; then
    echo "Error: Cannot connect to PostgreSQL on $DB_HOST:$DB_PORT"
    echo "Please check:"
    echo "1. PostgreSQL is running"
    echo "2. Correct host (set PGHOST=<host> if not lotta)"
    echo "3. Correct port (set PGPORT=<port> if not 5432)"
    echo "4. PostgreSQL is accepting connections"
    echo ""
    echo "Try testing connection manually:"
    echo "  psql -U $DB_SUPERUSER -h $DB_HOST -p $DB_PORT -d postgres -c '\l'"
    exit 1
fi

# Check if we can connect as superuser
echo "Testing superuser connection..."
if ! psql -U $DB_SUPERUSER -h $DB_HOST -p $DB_PORT -d postgres -c "SELECT version();" > /dev/null 2>&1; then
    echo "Error: Cannot connect as user $DB_SUPERUSER"
    echo "You may need to:"
    echo "1. Use a different superuser (set PGUSER=<username>)"
    echo "2. Set up password authentication"
    echo "3. Check pg_hba.conf settings"
    echo ""
    echo "Try manual connection:"
    echo "  psql -U $DB_SUPERUSER -h $DB_HOST -p $DB_PORT -d postgres"
    exit 1
fi

# Create test database
echo "Creating test database..."
createdb -U $DB_SUPERUSER -h $DB_HOST -p $DB_PORT $DB_NAME || echo "Database '$DB_NAME' may already exist"

# Check if test database was created or already exists
if ! psql -U $DB_SUPERUSER -h $DB_HOST -p $DB_PORT -d $DB_NAME -c "SELECT 1;" > /dev/null 2>&1; then
    echo "Error: Cannot access test database $DB_NAME"
    exit 1
fi

# Create test user
echo "Creating test user..."
psql -U $DB_SUPERUSER -h $DB_HOST -p $DB_PORT -d $DB_NAME -c "
    DO \$\$
    BEGIN
        IF NOT EXISTS (SELECT FROM pg_catalog.pg_roles WHERE rolname = '$DB_USER') THEN
            CREATE USER $DB_USER WITH PASSWORD '$DB_PASS';
            RAISE NOTICE 'User $DB_USER created';
        ELSE
            RAISE NOTICE 'User $DB_USER already exists';
        END IF;
    END
    \$\$;
"

# Grant permissions
echo "Granting permissions..."
psql -U $DB_SUPERUSER -h $DB_HOST -p $DB_PORT -d $DB_NAME -c "
    GRANT ALL PRIVILEGES ON DATABASE $DB_NAME TO $DB_USER;
    GRANT ALL PRIVILEGES ON SCHEMA public TO $DB_USER;
    GRANT CREATE ON SCHEMA public TO $DB_USER;
"

# Create test schema (minimal version of production schema)
psql -U rschaeff -d $DB_NAME -c "
CREATE SCHEMA IF NOT EXISTS ecod_schema;
CREATE SCHEMA IF NOT EXISTS pdb_analysis;

-- Basic tables for testing
CREATE TABLE IF NOT EXISTS ecod_schema.batch (
    id SERIAL PRIMARY KEY,
    batch_name VARCHAR(255),
    base_path TEXT,
    ref_version VARCHAR(50),
    status VARCHAR(50),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS ecod_schema.protein (
    id SERIAL PRIMARY KEY,
    pdb_id VARCHAR(4),
    chain_id VARCHAR(10),
    sequence_length INTEGER,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS ecod_schema.process_status (
    id SERIAL PRIMARY KEY,
    batch_id INTEGER REFERENCES ecod_schema.batch(id),
    protein_id INTEGER REFERENCES ecod_schema.protein(id),
    current_stage VARCHAR(100),
    status VARCHAR(50),
    error_message TEXT,
    is_representative BOOLEAN DEFAULT FALSE,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS ecod_schema.process_file (
    id SERIAL PRIMARY KEY,
    process_id INTEGER REFERENCES ecod_schema.process_status(id),
    file_type VARCHAR(100),
    file_path TEXT,
    file_exists BOOLEAN DEFAULT FALSE,
    file_size BIGINT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Grant permissions on all tables
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA ecod_schema TO $DB_USER;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA pdb_analysis TO $DB_USER;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA ecod_schema TO $DB_USER;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA pdb_analysis TO $DB_USER;
"

echo "Test database setup complete"
