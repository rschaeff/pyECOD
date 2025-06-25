-- Fixed Analysis function to resolve type mismatch
CREATE OR REPLACE FUNCTION pdb_analysis.analyze_propagation_potential(
    p_source_version text DEFAULT 'mini_pyecod_v2.0'
) RETURNS TABLE(
    sequence_md5 text,
    donor_count integer,
    target_count integer,
    total_domains integer,
    avg_confidence numeric,
    propagation_factor numeric
) AS $$
BEGIN
    RETURN QUERY
    WITH donor_sequences AS (
        SELECT DISTINCT
            ps.sequence_md5,
            COUNT(*) as donors,
            SUM(pp.domains_with_evidence) as total_domains,
            AVG(COALESCE(pd.confidence, 0))::numeric as avg_confidence  -- Fixed: cast to numeric
        FROM pdb_analysis.partition_proteins pp
        JOIN pdb_analysis.protein p ON pp.pdb_id = p.pdb_id AND pp.chain_id = p.chain_id
        JOIN pdb_analysis.protein_sequence ps ON p.id = ps.protein_id
        LEFT JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
        WHERE pp.process_version = p_source_version
          AND pp.is_classified = true
          AND pp.domains_with_evidence > 0
          AND ps.sequence_md5 IS NOT NULL
        GROUP BY ps.sequence_md5
    ),
    target_sequences AS (
        SELECT 
            ps.sequence_md5,
            COUNT(*) as targets
        FROM pdb_analysis.protein p
        JOIN pdb_analysis.protein_sequence ps ON p.id = ps.protein_id
        LEFT JOIN pdb_analysis.partition_proteins pp ON p.pdb_id = pp.pdb_id 
            AND p.chain_id = pp.chain_id 
            AND pp.process_version LIKE '%mini_pyecod%'
        WHERE ps.sequence_md5 IS NOT NULL
          AND pp.id IS NULL  -- No existing mini classification
        GROUP BY ps.sequence_md5
    )
    SELECT 
        d.sequence_md5,
        d.donors::integer,
        COALESCE(t.targets, 0)::integer,
        d.total_domains::integer,
        d.avg_confidence,
        CASE WHEN d.donors > 0 THEN COALESCE(t.targets, 0)::numeric / d.donors ELSE 0::numeric END
    FROM donor_sequences d
    LEFT JOIN target_sequences t ON d.sequence_md5 = t.sequence_md5
    WHERE COALESCE(t.targets, 0) > 0  -- Only sequences with propagation potential
    ORDER BY COALESCE(t.targets, 0) DESC, d.total_domains DESC;
END;
$$ LANGUAGE plpgsql;

-- Quick diagnostic function to troubleshoot the original error
CREATE OR REPLACE FUNCTION pdb_analysis.diagnose_propagation_types() 
RETURNS TABLE(
    check_name text,
    data_type text,
    sample_value text,
    row_count bigint
) AS $$
BEGIN
    -- Check confidence data types
    RETURN QUERY
    SELECT 
        'confidence_column_type'::text,
        pg_typeof(confidence)::text,
        confidence::text,
        1::bigint
    FROM pdb_analysis.partition_domains 
    LIMIT 1;
    
    -- Check AVG result type
    RETURN QUERY
    SELECT 
        'avg_confidence_type'::text,
        pg_typeof(AVG(confidence))::text,
        AVG(confidence)::text,
        COUNT(*)
    FROM pdb_analysis.partition_domains pd
    JOIN pdb_analysis.partition_proteins pp ON pd.protein_id = pp.id
    WHERE pp.process_version = 'mini_pyecod_v2.0';
    
    -- Check sequence data availability
    RETURN QUERY
    SELECT 
        'sequence_data_check'::text,
        'availability'::text,
        CASE WHEN COUNT(*) > 0 THEN 'available' ELSE 'missing' END::text,
        COUNT(*)
    FROM pdb_analysis.protein_sequence ps
    JOIN pdb_analysis.protein p ON ps.protein_id = p.id
    JOIN pdb_analysis.partition_proteins pp ON p.pdb_id = pp.pdb_id AND p.chain_id = pp.chain_id
    WHERE pp.process_version = 'mini_pyecod_v2.0';
END;
$$ LANGUAGE plpgsql;

-- Simplified analysis query for immediate testing
CREATE OR REPLACE FUNCTION pdb_analysis.simple_propagation_analysis()
RETURNS TABLE(
    total_v2_proteins bigint,
    sequences_with_md5 bigint,
    potential_targets bigint,
    propagatable_sequences bigint
) AS $$
BEGIN
    RETURN QUERY
    SELECT 
        (SELECT COUNT(*) FROM pdb_analysis.partition_proteins WHERE process_version = 'mini_pyecod_v2.0'),
        (SELECT COUNT(DISTINCT ps.sequence_md5) 
         FROM pdb_analysis.protein_sequence ps
         JOIN pdb_analysis.protein p ON ps.protein_id = p.id
         JOIN pdb_analysis.partition_proteins pp ON p.pdb_id = pp.pdb_id AND p.chain_id = pp.chain_id
         WHERE pp.process_version = 'mini_pyecod_v2.0' AND ps.sequence_md5 IS NOT NULL),
        (SELECT COUNT(*) 
         FROM pdb_analysis.protein p
         JOIN pdb_analysis.protein_sequence ps ON p.id = ps.protein_id
         LEFT JOIN pdb_analysis.partition_proteins pp ON p.pdb_id = pp.pdb_id 
             AND p.chain_id = pp.chain_id 
             AND pp.process_version LIKE '%mini_pyecod%'
         WHERE ps.sequence_md5 IS NOT NULL AND pp.id IS NULL),
        (SELECT COUNT(DISTINCT ps.sequence_md5)
         FROM pdb_analysis.protein p
         JOIN pdb_analysis.protein_sequence ps ON p.id = ps.protein_id
         LEFT JOIN pdb_analysis.partition_proteins pp ON p.pdb_id = pp.pdb_id 
             AND p.chain_id = pp.chain_id 
             AND pp.process_version LIKE '%mini_pyecod%'
         WHERE ps.sequence_md5 IS NOT NULL AND pp.id IS NULL
         AND EXISTS (
             SELECT 1 FROM pdb_analysis.protein p2
             JOIN pdb_analysis.protein_sequence ps2 ON p2.id = ps2.protein_id
             JOIN pdb_analysis.partition_proteins pp2 ON p2.pdb_id = pp2.pdb_id AND p2.chain_id = pp2.chain_id
             WHERE pp2.process_version = 'mini_pyecod_v2.0' AND ps2.sequence_md5 = ps.sequence_md5
         ));
END;
$$ LANGUAGE plpgsql;
