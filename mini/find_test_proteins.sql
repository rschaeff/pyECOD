-- SQL queries to find diverse test proteins from ecod_schema
-- Run these against your PostgreSQL database to find actual proteins in your batches

-- 1. Find proteins with domain summaries in batch 36
-- This gives you the universe of testable proteins
SELECT 
    p.pdb_id || '_' || p.chain_id as protein_id,
    p.length,
    ps.current_stage,
    ps.status,
    pf.file_path
FROM ecod_schema.protein p
JOIN ecod_schema.process_status ps ON p.id = ps.protein_id
JOIN ecod_schema.process_file pf ON ps.id = pf.process_id
WHERE ps.batch_id = 36 
    AND pf.file_type = 'domain_summary'
    AND pf.file_exists = true
ORDER BY p.length
LIMIT 20;

-- 2. Find small proteins (potential single domain)
-- Good for simple test cases
SELECT 
    p.pdb_id || '_' || p.chain_id as protein_id,
    p.length,
    COUNT(DISTINCT pf.file_type) as file_types
FROM ecod_schema.protein p
JOIN ecod_schema.process_status ps ON p.id = ps.protein_id
JOIN ecod_schema.process_file pf ON ps.id = pf.process_id
WHERE ps.batch_id = 36 
    AND p.length BETWEEN 50 AND 150
    AND pf.file_exists = true
GROUP BY p.id, p.pdb_id, p.chain_id, p.length
HAVING COUNT(DISTINCT CASE WHEN pf.file_type = 'domain_summary' THEN 1 END) > 0
ORDER BY p.length
LIMIT 10;

-- 3. Find large proteins (potential multi-domain)  
-- Good for complex test cases
SELECT 
    p.pdb_id || '_' || p.chain_id as protein_id,
    p.length,
    COUNT(DISTINCT pf.file_type) as file_types
FROM ecod_schema.protein p
JOIN ecod_schema.process_status ps ON p.id = ps.protein_id
JOIN ecod_schema.process_file pf ON ps.id = pf.process_id
WHERE ps.batch_id = 36 
    AND p.length > 500
    AND pf.file_exists = true
GROUP BY p.id, p.pdb_id, p.chain_id, p.length
HAVING COUNT(DISTINCT CASE WHEN pf.file_type = 'domain_summary' THEN 1 END) > 0
ORDER BY p.length DESC
LIMIT 10;

-- 4. Find proteins with all evidence types
-- These are well-characterized and good for testing
SELECT 
    p.pdb_id || '_' || p.chain_id as protein_id,
    p.length,
    COUNT(DISTINCT CASE WHEN pf.file_type = 'fasta' THEN 1 END) as has_fasta,
    COUNT(DISTINCT CASE WHEN pf.file_type = 'chain_blast' THEN 1 END) as has_chain_blast,
    COUNT(DISTINCT CASE WHEN pf.file_type = 'domain_blast' THEN 1 END) as has_domain_blast,
    COUNT(DISTINCT CASE WHEN pf.file_type = 'hhblits_a3m' THEN 1 END) as has_hhblits,
    COUNT(DISTINCT CASE WHEN pf.file_type = 'hhsearch' THEN 1 END) as has_hhsearch,
    COUNT(DISTINCT CASE WHEN pf.file_type = 'domain_summary' THEN 1 END) as has_summary
FROM ecod_schema.protein p
JOIN ecod_schema.process_status ps ON p.id = ps.protein_id
JOIN ecod_schema.process_file pf ON ps.id = pf.process_id
WHERE ps.batch_id = 36 
    AND pf.file_exists = true
GROUP BY p.id, p.pdb_id, p.chain_id, p.length
HAVING COUNT(DISTINCT CASE WHEN pf.file_type = 'domain_summary' THEN 1 END) > 0
    AND COUNT(DISTINCT pf.file_type) >= 5  -- Has most evidence types
ORDER BY p.length
LIMIT 20;

-- 5. Find proteins at specific stages
-- Good for understanding the pipeline state
SELECT 
    ps.current_stage,
    COUNT(*) as protein_count,
    AVG(p.length) as avg_length,
    MIN(p.length) as min_length,
    MAX(p.length) as max_length
FROM ecod_schema.protein p
JOIN ecod_schema.process_status ps ON p.id = ps.protein_id
WHERE ps.batch_id = 36
GROUP BY ps.current_stage
ORDER BY protein_count DESC;

-- 6. Find proteins by size distribution
-- Helps ensure test coverage across sizes
SELECT 
    CASE 
        WHEN p.length < 100 THEN 'very_small'
        WHEN p.length < 200 THEN 'small'
        WHEN p.length < 400 THEN 'medium'
        WHEN p.length < 800 THEN 'large'
        ELSE 'very_large'
    END as size_category,
    COUNT(*) as protein_count,
    array_agg(DISTINCT p.pdb_id || '_' || p.chain_id ORDER BY p.length LIMIT 5) as example_proteins
FROM ecod_schema.protein p
JOIN ecod_schema.process_status ps ON p.id = ps.protein_id
JOIN ecod_schema.process_file pf ON ps.id = pf.process_id
WHERE ps.batch_id = 36 
    AND pf.file_type = 'domain_summary'
    AND pf.file_exists = true
GROUP BY size_category
ORDER BY 
    CASE size_category
        WHEN 'very_small' THEN 1
        WHEN 'small' THEN 2
        WHEN 'medium' THEN 3
        WHEN 'large' THEN 4
        WHEN 'very_large' THEN 5
    END;

-- 7. Find potential multi-chain complexes
-- Good for future multi-chain testing
SELECT 
    p.pdb_id,
    COUNT(DISTINCT p.chain_id) as chain_count,
    string_agg(DISTINCT p.chain_id, ', ' ORDER BY p.chain_id) as chains,
    array_agg(DISTINCT p.length ORDER BY p.length) as chain_lengths
FROM ecod_schema.protein p
JOIN ecod_schema.process_status ps ON p.id = ps.protein_id
WHERE ps.batch_id = 36
GROUP BY p.pdb_id
HAVING COUNT(DISTINCT p.chain_id) > 1
ORDER BY chain_count DESC
LIMIT 10;

-- 8. Get batch overview
-- Understand what's in your batch
SELECT 
    b.id,
    b.batch_name,
    b.total_items,
    b.completed_items,
    b.status,
    b.created_at,
    COUNT(DISTINCT ps.protein_id) as unique_proteins
FROM ecod_schema.batch b
LEFT JOIN ecod_schema.process_status ps ON b.id = ps.batch_id
WHERE b.id = 36
GROUP BY b.id, b.batch_name, b.total_items, b.completed_items, b.status, b.created_at;

-- 9. Find proteins with specific file patterns
-- Useful for chain BLAST decomposition candidates
WITH protein_files AS (
    SELECT 
        p.pdb_id || '_' || p.chain_id as protein_id,
        p.length,
        bool_or(pf.file_type = 'chain_blast' AND pf.file_exists) as has_chain_blast,
        bool_or(pf.file_type = 'domain_blast' AND pf.file_exists) as has_domain_blast,
        bool_or(pf.file_type = 'domain_summary' AND pf.file_exists) as has_summary
    FROM ecod_schema.protein p
    JOIN ecod_schema.process_status ps ON p.id = ps.protein_id
    JOIN ecod_schema.process_file pf ON ps.id = pf.process_id
    WHERE ps.batch_id = 36
    GROUP BY p.id, p.pdb_id, p.chain_id, p.length
)
SELECT 
    protein_id,
    length,
    CASE 
        WHEN has_chain_blast AND length > 400 THEN 'chain_blast_decomposition_candidate'
        WHEN has_domain_blast AND NOT has_chain_blast THEN 'domain_blast_only'
        WHEN has_chain_blast AND has_domain_blast THEN 'full_blast_coverage'
        ELSE 'other'
    END as test_category
FROM protein_files
WHERE has_summary = true
ORDER BY length DESC
LIMIT 20;

-- 10. Random sample for diversity
-- Get a random sample across different sizes
SELECT 
    p.pdb_id || '_' || p.chain_id as protein_id,
    p.length,
    ps.current_stage
FROM ecod_schema.protein p
JOIN ecod_schema.process_status ps ON p.id = ps.protein_id
JOIN ecod_schema.process_file pf ON ps.id = pf.process_id
WHERE ps.batch_id = 36 
    AND pf.file_type = 'domain_summary'
    AND pf.file_exists = true
ORDER BY RANDOM()
LIMIT 20;
