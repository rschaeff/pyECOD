#!/bin/bash
# setup_database.sh - Setup the PyECOD database schema

set -e  # Exit on any error

# Configuration
DB_NAME="ecod_protein"
DB_USER="ecod"
DB_PASSWORD="your_secure_password"
DB_HOST="localhost"
DB_PORT="5432"

# Create the database user and database
echo "Creating database user and database..."
sudo -u postgres psql -c "CREATE USER $DB_USER WITH PASSWORD '$DB_PASSWORD';"
sudo -u postgres psql -c "CREATE DATABASE $DB_NAME OWNER $DB_USER;"

# Create the required schemas
echo "Creating schemas..."
PGPASSWORD=$DB_PASSWORD psql -h $DB_HOST -p $DB_PORT -U $DB_USER -d $DB_NAME -c "CREATE SCHEMA IF NOT EXISTS ecod_schema;"
PGPASSWORD=$DB_PASSWORD psql -h $DB_HOST -p $DB_PORT -U $DB_USER -d $DB_NAME -c "CREATE SCHEMA IF NOT EXISTS pdb_analysis;"

# Run the migrations
echo "Running migrations..."
python migrate.py --config config/config.yml --migrations-dir ecod/db/migrations --verbose

# Install the schema from dump (for pdb_analysis schema)
echo "Installing schema from dump file..."
PGPASSWORD=$DB_PASSWORD psql -h $DB_HOST -p $DB_PORT -U $DB_USER -d $DB_NAME -f pdb_analysis.schema

echo "Database setup complete!"