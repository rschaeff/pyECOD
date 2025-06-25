-- SQL-Only Sequence Propagation Procedure
-- Propagates classifications from high-quality mini results to sequence-identical proteins
-- This is a cleaner approach than mixing propagation with import logic

-- Function to propagate mini classifications to sequence-identical proteins
CREATE OR REPLACE FUNCTION pdb_analysis.propagate_mini_classifications(
    p_source_version text DEFAULT 'mini_pyecod_v2.0',
    p_target_version text DEFAULT 'mini_pyecod_v2.0_propagated',
    p_max_targets_per_sequence integer DEFAULT 10,
    p_dry_run boolean DEFAULT true
) RETURNS TABLE(
    sequence_md5 text,
    donor_protein_id text,
    target_count integer,
    domains_propagated integer,
    success boolean,
    message text
) AS $$
DECLARE
    propagation_record RECORD;
    target_record RECORD;
    domain_record RECORD;
    new_partition_id integer;
    domains_created integer;
    total_propagated integer := 0;
    total_targets integer := 0;
BEGIN
    RAISE NOTICE 'Starting sequence propagation: % -> % (max %/seq) [%]', 
        p_source_version, p_target_version, p_max_targets_per_sequence,
        CASE WHEN p_dry_run THEN 'DRY RUN' ELSE 'EXECUTING' END;

    -- Find all high-quality mini donors with their sequence MD5s
    FOR propagation_record IN
        SELECT DISTINCT
            ps.sequence_md5,
            pp.pdb_id as donor_pdb_id,
            pp.chain_id as donor_chain_id,
            pp.id as donor_partition_id,
            pp.domains_with_evidence,
            pp.is_classified,
            pp.batch_id,
            pp.reference_version,
            pp.sequence_length,
            pp.coverage,
            pp.residues_assigned,
            pp.fully_classified_domains,
            pp.algorithm_version,
            pp.git_commit_hash,
            pp.process_parameters
        FROM pdb_analysis.partition_proteins pp
        JOIN pdb_analysis.protein p ON pp.pdb_id = p.pdb_id AND pp.chain_id = p.chain_id
        JOIN pdb_analysis.protein_sequence ps ON p.id = ps.protein_id
        WHERE pp.process_version = p_source_version
          AND pp.is_classified = true
          AND pp.domains_with_evidence > 0
          AND ps.sequence_md5 IS NOT NULL
        ORDER BY pp.domains_with_evidence DESC, pp.coverage DESC
    LOOP
        domains_created := 0;
        
        -- Find target proteins with same sequence but no classification
        FOR target_record IN
            SELECT 
                p.id as target_protein_db_id,
                p.pdb_id as target_pdb_id,
                p.chain_id as target_chain_id,
                p.source_id as target_source_id,
                p.length as target_length
            FROM pdb_analysis.protein p
            JOIN pdb_analysis.protein_sequence ps ON p.id = ps.protein_id
            LEFT JOIN pdb_analysis.partition_proteins pp ON p.pdb_id = pp.pdb_id 
                AND p.chain_id = pp.chain_id 
                AND pp.process_version LIKE '%mini_pyecod%'
            WHERE ps.sequence_md5 = propagation_record.sequence_md5
              AND p.source_id != (propagation_record.donor_pdb_id || '_' || propagation_record.donor_chain_id)
              AND pp.id IS NULL  -- No existing mini classification
            ORDER BY p.id
            LIMIT p_max_targets_per_sequence
        LOOP
            IF NOT p_dry_run THEN
                -- Create partition_proteins record for target
                INSERT INTO pdb_analysis.partition_proteins (
                    pdb_id, chain_id, batch_id, reference_version, is_classified,
                    sequence_length, coverage, residues_assigned, domains_with_evidence,
                    fully_classified_domains, process_version, algorithm_version,
                    git_commit_hash, process_parameters
                ) VALUES (
                    target_record.target_pdb_id,
                    target_record.target_chain_id,
                    propagation_record.batch_id,  -- This should be integer from donor
                    propagation_record.reference_version,
                    propagation_record.is_classified,
                    target_record.target_length,
                    propagation_record.coverage,
                    propagation_record.residues_assigned,
                    propagation_record.domains_with_evidence,
                    propagation_record.fully_classified_domains,
                    p_target_version,  -- Mark as propagated
                    propagation_record.algorithm_version,
                    propagation_record.git_commit_hash,
                    jsonb_build_object(
                        'propagated_from', propagation_record.donor_pdb_id || '_' || propagation_record.donor_chain_id,
                        'propagation_method', 'sequence_md5_match',
                        'source_sequence_md5', propagation_record.sequence_md5,
                        'propagation_timestamp', CURRENT_TIMESTAMP
                    ) || COALESCE(propagation_record.process_parameters, '{}'::jsonb)
                ) RETURNING id INTO new_partition_id;

                -- Copy all domains from donor to target
                FOR domain_record IN
                    SELECT * FROM pdb_analysis.partition_domains
                    WHERE protein_id = propagation_record.donor_partition_id
                    ORDER BY domain_number
                LOOP
                    INSERT INTO pdb_analysis.partition_domains (
                        protein_id, domain_number, domain_id, start_pos, end_pos, range,
                        source, source_id, confidence, t_group, h_group, x_group, a_group,
                        is_manual_rep, is_f70, is_f40, is_f99, created_at,
                        primary_evidence_type, primary_evidence_id, evidence_evalue,
                        evidence_query_range, evidence_hit_range, original_range,
                        optimization_actions, pdb_range, pdb_start, pdb_end, length,
                        hit_ecod_domain_id
                    ) VALUES (
                        new_partition_id,
                        domain_record.domain_number,
                        target_record.target_source_id || '_d' || domain_record.domain_number || '_propagated',
                        domain_record.start_pos,
                        domain_record.end_pos,
                        domain_record.range,
                        'mini_pyecod_propagated',  -- Mark source as propagated
                        domain_record.source_id,
                        domain_record.confidence,
                        domain_record.t_group,
                        domain_record.h_group,
                        domain_record.x_group,
                        domain_record.a_group,
                        domain_record.is_manual_rep,
                        domain_record.is_f70,
                        domain_record.is_f40,
                        domain_record.is_f99,
                        CURRENT_TIMESTAMP,
                        domain_record.primary_evidence_type,
                        domain_record.primary_evidence_id,
                        domain_record.evidence_evalue,
                        domain_record.evidence_query_range,
                        domain_record.evidence_hit_range,
                        domain_record.original_range,
                        domain_record.optimization_actions,
                        domain_record.pdb_range,
                        domain_record.pdb_start,
                        domain_record.pdb_end,
                        domain_record.length,
                        domain_record.hit_ecod_domain_id
                    );

                    domains_created := domains_created + 1;
                END LOOP;

                -- Copy domain evidence (optional - links to original evidence)
                INSERT INTO pdb_analysis.domain_evidence (
                    domain_id, evidence_type, source_id, domain_ref_id, hit_id,
                    pdb_id, chain_id, confidence, probability, evalue, score,
                    hsp_count, is_discontinuous, t_group, h_group, x_group, a_group,
                    query_range, hit_range, created_at
                )
                SELECT 
                    pd_new.id as domain_id,
                    de.evidence_type,
                    de.source_id,
                    de.domain_ref_id,
                    de.hit_id,
                    de.pdb_id,
                    de.chain_id,
                    de.confidence,
                    de.probability,
                    de.evalue,
                    de.score,
                    de.hsp_count,
                    de.is_discontinuous,
                    de.t_group,
                    de.h_group,
                    de.x_group,
                    de.a_group,
                    de.query_range,
                    de.hit_range,
                    CURRENT_TIMESTAMP
                FROM pdb_analysis.partition_domains pd_old
                JOIN pdb_analysis.domain_evidence de ON pd_old.id = de.domain_id
                JOIN pdb_analysis.partition_domains pd_new ON pd_new.protein_id = new_partition_id 
                    AND pd_new.domain_number = pd_old.domain_number
                WHERE pd_old.protein_id = propagation_record.donor_partition_id;

            ELSE
                -- Dry run: just count domains
                SELECT COUNT(*) INTO domains_created
                FROM pdb_analysis.partition_domains
                WHERE protein_id = propagation_record.donor_partition_id;
            END IF;

            total_targets := total_targets + 1;
            total_propagated := total_propagated + domains_created;
        END LOOP;

        -- Return result for this sequence
        RETURN QUERY SELECT 
            propagation_record.sequence_md5,
            propagation_record.donor_pdb_id || '_' || propagation_record.donor_chain_id,
            (SELECT COUNT(*) FROM (
                SELECT 1 FROM pdb_analysis.protein p
                JOIN pdb_analysis.protein_sequence ps ON p.id = ps.protein_id
                LEFT JOIN pdb_analysis.partition_proteins pp ON p.pdb_id = pp.pdb_id 
                    AND p.chain_id = pp.chain_id 
                    AND pp.process_version LIKE '%mini_pyecod%'
                WHERE ps.sequence_md5 = propagation_record.sequence_md5
                  AND p.source_id != (propagation_record.donor_pdb_id || '_' || propagation_record.donor_chain_id)
                  AND pp.id IS NULL
                LIMIT p_max_targets_per_sequence
            ) AS targets)::integer,
            domains_created,
            true,
            CASE WHEN p_dry_run THEN 'DRY RUN' ELSE 'SUCCESS' END;
            
    END LOOP;

    RAISE NOTICE 'Propagation complete: % targets, % domains propagated', 
        total_targets, total_propagated;
        
    -- Return summary
    RETURN QUERY SELECT 
        'SUMMARY'::text,
        'ALL_SEQUENCES'::text,
        total_targets,
        total_propagated,
        true,
        format('Total: %s targets, %s domains %s', 
            total_targets, total_propagated,
            CASE WHEN p_dry_run THEN '(DRY RUN)' ELSE 'propagated' END
        );
END;
$$ LANGUAGE plpgsql;

-- Convenience function for most common use case
CREATE OR REPLACE FUNCTION pdb_analysis.propagate_mini_v2_classifications(
    p_dry_run boolean DEFAULT true
) RETURNS TABLE(
    sequence_md5 text,
    donor_protein_id text,
    target_count integer,
    domains_propagated integer,
    success boolean,
    message text
) AS $$
BEGIN
    RETURN QUERY SELECT * FROM pdb_analysis.propagate_mini_classifications(
        'mini_pyecod_v2.0',
        'mini_pyecod_v2.0_propagated',
        10,
        p_dry_run
    );
END;
$$ LANGUAGE plpgsql;

-- Analysis function to assess propagation potential
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
            AVG(COALESCE(pd.confidence, 0)) as avg_confidence
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
        CASE WHEN d.donors > 0 THEN COALESCE(t.targets, 0)::numeric / d.donors ELSE 0 END
    FROM donor_sequences d
    LEFT JOIN target_sequences t ON d.sequence_md5 = t.sequence_md5
    WHERE COALESCE(t.targets, 0) > 0  -- Only sequences with propagation potential
    ORDER BY COALESCE(t.targets, 0) DESC, d.total_domains DESC;
END;
$$ LANGUAGE plpgsql;

-- Example usage:

-- 1. Analyze propagation potential first:
-- SELECT * FROM pdb_analysis.analyze_propagation_potential('mini_pyecod_v2.0');

-- 2. Run dry run to see what would be propagated:
-- SELECT * FROM pdb_analysis.propagate_mini_v2_classifications(true);

-- 3. Execute actual propagation:
-- SELECT * FROM pdb_analysis.propagate_mini_v2_classifications(false);

-- 4. Check results:
-- SELECT process_version, COUNT(*) as protein_count, SUM(domains_with_evidence) as total_domains
-- FROM pdb_analysis.partition_proteins 
-- WHERE process_version LIKE '%mini_pyecod%'
-- GROUP BY process_version
-- ORDER BY process_version;
