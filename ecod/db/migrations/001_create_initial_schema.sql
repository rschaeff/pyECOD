# ecod/db/migrations/001_create_initial_schema.sql
-- Create schema if not exists
CREATE SCHEMA IF NOT EXISTS ecod_schema;

-- Core data tables
CREATE TABLE ecod_schema.protein (
    id SERIAL PRIMARY KEY,
    pdb_id VARCHAR(10) NOT NULL,
    chain_id VARCHAR(10) NOT NULL,
    source_id VARCHAR(50) UNIQUE NOT NULL,
    length INTEGER NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (pdb_id, chain_id)
);

CREATE TABLE ecod_schema.protein_sequence (
    id SERIAL PRIMARY KEY,
    protein_id INTEGER REFERENCES ecod_schema.protein(id),
    sequence TEXT NOT NULL,
    md5_hash VARCHAR(32) UNIQUE NOT NULL
);

-- Processing workflow tables
CREATE TABLE ecod_schema.batch (
    id SERIAL PRIMARY KEY,
    batch_name VARCHAR(100) NOT NULL,
    base_path TEXT NOT NULL,
    type VARCHAR(20) NOT NULL, -- 'blast', 'hhsearch', 'classification'
    ref_version VARCHAR(50) NOT NULL,
    total_items INTEGER NOT NULL,
    completed_items INTEGER DEFAULT 0,
    status VARCHAR(20) DEFAULT 'created',
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    completed_at TIMESTAMP
);

CREATE TABLE ecod_schema.process_status (
    id SERIAL PRIMARY KEY,
    protein_id INTEGER REFERENCES ecod_schema.protein(id),
    batch_id INTEGER REFERENCES ecod_schema.batch(id),
    current_stage VARCHAR(50) NOT NULL,
    status VARCHAR(20) DEFAULT 'pending',
    is_representative BOOLEAN DEFAULT FALSE,
    relative_path TEXT,
    error_message TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE ecod_schema.process_file (
    id SERIAL PRIMARY KEY,
    process_id INTEGER REFERENCES ecod_schema.process_status(id),
    file_type VARCHAR(50) NOT NULL,
    file_path TEXT NOT NULL,
    file_exists BOOLEAN DEFAULT FALSE,
    file_size BIGINT,
    last_checked TIMESTAMP
);

-- Job tracking tables
CREATE TABLE ecod_schema.job (
    id SERIAL PRIMARY KEY,
    batch_id INTEGER REFERENCES ecod_schema.batch(id),
    job_type VARCHAR(50) NOT NULL,
    slurm_job_id VARCHAR(50),
    status VARCHAR(20) DEFAULT 'submitted',
    items_count INTEGER,
    submission_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    completion_time TIMESTAMP
);

CREATE TABLE ecod_schema.job_item (
    id SERIAL PRIMARY KEY,
    job_id INTEGER REFERENCES ecod_schema.job(id),
    process_id INTEGER REFERENCES ecod_schema.process_status(id),
    status VARCHAR(20) DEFAULT 'pending'
);

-- Reference data tables
CREATE TABLE ecod_schema.ecod_version (
    id SERIAL PRIMARY KEY,
    version_name VARCHAR(50) UNIQUE NOT NULL,
    release_date DATE,
    is_current BOOLEAN DEFAULT FALSE
);

CREATE TABLE ecod_schema.reference_resource (
    id SERIAL PRIMARY KEY,
    version_id INTEGER REFERENCES ecod_schema.ecod_version(id),
    resource_type VARCHAR(50) NOT NULL,
    resource_path TEXT NOT NULL,
    UNIQUE (version_id, resource_type)
);