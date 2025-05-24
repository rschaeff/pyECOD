#!/usr/bin/env python3
"""
Investigate Domain Summary Generation Failures

This script specifically investigates why 84.6% of proteins fail at the domain summary 
generation stage (fasta→blast but no summary).

Usage:
    python scripts/investigate_domain_summary_failures.py --config config.yml [options]
"""

import os
import sys
import logging
import argparse
import json
import psycopg2
from psycopg2.extras import RealDictCursor
import yaml
import glob
from datetime import datetime
from typing import List, Dict, Any, Optional
from collections import Counter, defaultdict


def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    format_str = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=level, format=format_str)
    return logging.getLogger(__name__)


def parse_config(config_path):
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def get_db_connection(config):
    db_config = config.get('database', {})
    try:
        conn = psycopg2.connect(
            host=db_config.get('host', 'dione'),
            port=db_config.get('port', 45000),
            dbname=db_config.get('name', 'ecod_protein'),
            user=db_config.get('user', 'ecod'),
            password=db_config.get('password', '')
        )
        return conn
    except psycopg2.Error as e:
        logging.error(f"Database connection error: {e}")
        raise


def get_domain_summary_failures(conn, batch_ids, limit=None):
    """Get proteins that have BLAST results but no domain summary."""
    
    batch_filter = f"AND ps.batch_id IN ({','.join(map(str, batch_ids))})"
    limit_clause = f"LIMIT {limit}" if limit else ""
    
    query = f"""
    WITH file_status AS (
        SELECT 
            ps.id as process_id,
            ps.protein_id,
            ep.source_id,
            ep.pdb_id,
            ep.chain_id,
            ep.length,
            ps.batch_id,
            ps.current_stage,
            ps.status,
            ps.error_message,
            ps.updated_at,
            b.batch_name,
            b.base_path,
            
            -- File existence flags
            bool_or(pf.file_type = 'fasta' AND pf.file_exists) as has_fasta,
            bool_or(pf.file_type LIKE '%blast%' AND pf.file_exists) as has_blast,
            bool_or(pf.file_type = 'domain_summary' AND pf.file_exists) as has_summary,
            
            -- File paths for investigation
            string_agg(CASE WHEN pf.file_type = 'fasta' AND pf.file_exists THEN pf.file_path END, ', ') as fasta_paths,
            string_agg(CASE WHEN pf.file_type LIKE '%blast%' AND pf.file_exists THEN pf.file_path END, ', ') as blast_paths,
            string_agg(CASE WHEN pf.file_type = 'domain_summary' THEN pf.file_path END, ', ') as summary_paths,
            
            -- File sizes
            max(CASE WHEN pf.file_type LIKE '%blast%' AND pf.file_exists THEN pf.file_size END) as blast_file_size
            
        FROM ecod_schema.process_status ps
        JOIN ecod_schema.protein ep ON ps.protein_id = ep.id
        JOIN ecod_schema.batch b ON ps.batch_id = b.id
        LEFT JOIN ecod_schema.process_file pf ON ps.id = pf.process_id
        WHERE 1=1 {batch_filter}
        GROUP BY ps.id, ps.protein_id, ep.source_id, ep.pdb_id, ep.chain_id, 
                 ep.length, ps.batch_id, ps.current_stage, ps.status, 
                 ps.error_message, ps.updated_at, b.batch_name, b.base_path
    )
    SELECT *
    FROM file_status
    WHERE has_fasta = true 
      AND has_blast = true 
      AND has_summary = false
    ORDER BY length, batch_id, source_id
    {limit_clause}
    """
    
    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(query)
        return cur.fetchall()


def analyze_blast_file_quality(blast_file_path):
    """Analyze a BLAST file to understand if it's valid and has results."""
    
    if not os.path.exists(blast_file_path):
        return {"error": "File does not exist", "path": blast_file_path}
    
    try:
        file_size = os.path.getsize(blast_file_path)
        
        if file_size == 0:
            return {"error": "Empty file", "size": 0, "path": blast_file_path}
        
        # Read the file and analyze content
        with open(blast_file_path, 'r') as f:
            content = f.read(10000)  # First 10KB should be enough to assess
        
        analysis = {
            "valid_file": True,
            "size": file_size,
            "path": blast_file_path
        }
        
        # Check for BLAST XML format
        if content.startswith('<?xml'):
            analysis["format"] = "xml"
            analysis["has_hits"] = "<Hit>" in content
            analysis["hit_count"] = content.count("<Hit>")
        # Check for tabular format
        elif '\t' in content:
            analysis["format"] = "tabular"
            lines = content.split('\n')
            analysis["line_count"] = len([l for l in lines if l.strip()])
            analysis["has_hits"] = len([l for l in lines if l.strip() and not l.startswith('#')]) > 0
        # Check for text format
        else:
            analysis["format"] = "text"
            analysis["has_hits"] = "Score =" in content or "Expect =" in content
        
        # Look for error indicators
        if "error" in content.lower() or "failed" in content.lower():
            analysis["has_errors"] = True
            analysis["error_snippets"] = [line.strip() for line in content.split('\n') 
                                        if 'error' in line.lower() or 'failed' in line.lower()]
        
        return analysis
        
    except Exception as e:
        return {"error": f"File read error: {str(e)}", "path": blast_file_path}


def check_expected_domain_summary_path(protein_data, base_data_dir="/data/ecod/pdb_updates"):
    """Check if domain summary should exist at expected path."""
    
    batch_path = os.path.join(base_data_dir, "batches", str(protein_data['batch_id']))
    expected_summary_path = os.path.join(
        batch_path, 
        "domains", 
        f"{protein_data['pdb_id']}_{protein_data['chain_id']}.develop291.domains_v14.xml"
    )
    
    return {
        "expected_path": expected_summary_path,
        "exists": os.path.exists(expected_summary_path),
        "directory_exists": os.path.exists(os.path.dirname(expected_summary_path))
    }


def analyze_length_30_42_pattern(conn, batch_ids):
    """Specifically analyze the 30-42 residue length range with high error rates."""
    
    query = f"""
    SELECT 
        ep.length,
        ep.source_id,
        ep.pdb_id,
        ep.chain_id,
        ps.status,
        ps.error_message,
        ps.current_stage,
        
        -- Check what files exist
        bool_or(pf.file_type = 'fasta' AND pf.file_exists) as has_fasta,
        bool_or(pf.file_type LIKE '%blast%' AND pf.file_exists) as has_blast,
        bool_or(pf.file_type = 'domain_summary' AND pf.file_exists) as has_summary,
        
        string_agg(CASE WHEN pf.file_type LIKE '%blast%' AND pf.file_exists THEN pf.file_path END, ', ') as blast_paths
        
    FROM ecod_schema.process_status ps
    JOIN ecod_schema.protein ep ON ps.protein_id = ep.id
    LEFT JOIN ecod_schema.process_file pf ON ps.id = pf.process_id
    WHERE ps.batch_id IN ({','.join(map(str, batch_ids))})
      AND ep.length BETWEEN 30 AND 42
      AND ps.status = 'error'
    GROUP BY ep.length, ep.source_id, ep.pdb_id, ep.chain_id, 
             ps.status, ps.error_message, ps.current_stage
    ORDER BY ep.length, ep.source_id
    LIMIT 50
    """
    
    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(query)
        return cur.fetchall()


def investigate_domain_summary_generation_process(protein_data, base_data_dir="/data/ecod/pdb_updates"):
    """Try to understand why domain summary generation failed for a specific protein."""
    
    print(f"\n{'='*80}")
    print(f"INVESTIGATING DOMAIN SUMMARY FAILURE: {protein_data['source_id']}")
    print(f"{'='*80}")
    
    print(f"📋 Basic Info:")
    print(f"   PDB ID: {protein_data['pdb_id']}")
    print(f"   Chain: {protein_data['chain_id']}")
    print(f"   Length: {protein_data['length']}")
    print(f"   Batch: {protein_data['batch_id']} ({protein_data['batch_name']})")
    print(f"   Status: {protein_data['status']}")
    print(f"   Stage: {protein_data['current_stage']}")
    print(f"   Error: {protein_data['error_message'] or 'None recorded'}")
    
    # Check files
    print(f"\n📁 File Analysis:")
    print(f"   ✅ FASTA: {protein_data['has_fasta']}")
    print(f"   ✅ BLAST: {protein_data['has_blast']}")
    print(f"   ❌ Summary: {protein_data['has_summary']}")
    
    # Analyze BLAST files
    if protein_data['blast_paths']:
        blast_paths = [p.strip() for p in protein_data['blast_paths'].split(',') if p.strip()]
        for blast_path in blast_paths:
            print(f"\n🔍 BLAST File Analysis: {os.path.basename(blast_path)}")
            blast_analysis = analyze_blast_file_quality(blast_path)
            
            if 'error' in blast_analysis:
                print(f"   ❌ {blast_analysis['error']}")
            else:
                print(f"   ✅ Format: {blast_analysis.get('format', 'unknown')}")
                print(f"   📊 Size: {blast_analysis['size']} bytes")
                print(f"   🎯 Has hits: {blast_analysis.get('has_hits', 'unknown')}")
                if 'hit_count' in blast_analysis:
                    print(f"   📈 Hit count: {blast_analysis['hit_count']}")
                if blast_analysis.get('has_errors'):
                    print(f"   ⚠️  Contains errors: {blast_analysis['error_snippets']}")
    
    # Check expected domain summary location
    summary_check = check_expected_domain_summary_path(protein_data, base_data_dir)
    print(f"\n📂 Expected Domain Summary:")
    print(f"   Path: {summary_check['expected_path']}")
    print(f"   Exists: {summary_check['exists']}")
    print(f"   Directory exists: {summary_check['directory_exists']}")
    
    if not summary_check['directory_exists']:
        print(f"   ❌ Domain summary directory missing!")
    
    # Check for any domain summary files in the directory
    if summary_check['directory_exists']:
        domain_dir = os.path.dirname(summary_check['expected_path'])
        pattern = f"{protein_data['pdb_id']}_{protein_data['chain_id']}*.xml"
        matching_files = glob.glob(os.path.join(domain_dir, pattern))
        
        if matching_files:
            print(f"   🔍 Found related files:")
            for file_path in matching_files:
                file_size = os.path.getsize(file_path) if os.path.exists(file_path) else 0
                print(f"     - {os.path.basename(file_path)} ({file_size} bytes)")
        else:
            print(f"   ❌ No domain summary files found for this protein")
    
    return {
        'source_id': protein_data['source_id'],
        'blast_analysis': blast_analysis if protein_data['blast_paths'] else None,
        'summary_check': summary_check,
        'failure_hypothesis': determine_failure_hypothesis(protein_data, blast_analysis if protein_data['blast_paths'] else None, summary_check)
    }


def determine_failure_hypothesis(protein_data, blast_analysis, summary_check):
    """Determine the most likely reason for domain summary generation failure."""
    
    hypotheses = []
    
    # Length-based hypothesis
    if protein_data['length'] <= 30:
        hypotheses.append("Peptide: Too short for domain analysis")
    elif 30 <= protein_data['length'] <= 42:
        hypotheses.append("Transition zone: Length in problematic 30-42 range")
    
    # File-based hypotheses
    if not summary_check['directory_exists']:
        hypotheses.append("Infrastructure: Domain summary directory missing")
    
    if blast_analysis:
        if blast_analysis.get('error'):
            hypotheses.append(f"BLAST file issue: {blast_analysis['error']}")
        elif blast_analysis.get('size', 0) < 100:
            hypotheses.append("BLAST results: File too small, likely no hits")
        elif not blast_analysis.get('has_hits', True):
            hypotheses.append("BLAST results: No significant hits found")
        elif blast_analysis.get('has_errors'):
            hypotheses.append("BLAST processing: Contains error messages")
    
    # Process status hypotheses
    if protein_data['error_message']:
        if 'timeout' in protein_data['error_message'].lower():
            hypotheses.append("Resource: Processing timeout")
        elif 'memory' in protein_data['error_message'].lower():
            hypotheses.append("Resource: Memory exhaustion")
        elif 'permission' in protein_data['error_message'].lower():
            hypotheses.append("Infrastructure: File permission issue")
    
    if not hypotheses:
        hypotheses.append("Unknown: Domain summary generation failed for unclear reason")
    
    return hypotheses


def main():
    parser = argparse.ArgumentParser(description='Investigate domain summary generation failures')
    parser.add_argument('--config', required=True, help='Path to configuration file')
    parser.add_argument('--batch-ids', help='Comma-separated batch IDs (default: problematic batches)')
    parser.add_argument('--sample-size', type=int, default=20, help='Number of failures to investigate in detail')
    parser.add_argument('--focus-length-30-42', action='store_true', help='Focus on the problematic 30-42 length range')
    parser.add_argument('--output', help='Output JSON file path')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')

    args = parser.parse_args()

    logger = setup_logging(args.verbose)

    try:
        config = parse_config(args.config)
    except Exception as e:
        logger.error(f"Error parsing config file: {str(e)}")
        sys.exit(1)

    # Parse batch IDs
    if args.batch_ids:
        try:
            batch_ids = [int(bid.strip()) for bid in args.batch_ids.split(',')]
        except ValueError:
            logger.error("Invalid batch IDs format")
            sys.exit(1)
    else:
        batch_ids = [18, 19, 20, 21, 22, 23, 25, 26, 30, 33, 43, 50, 51]

    try:
        conn = get_db_connection(config)
    except Exception as e:
        logger.error(f"Error connecting to database: {str(e)}")
        sys.exit(1)

    try:
        if args.focus_length_30_42:
            logger.info("Analyzing the problematic 30-42 length range...")
            length_30_42_failures = analyze_length_30_42_pattern(conn, batch_ids)
            
            print(f"\n{'='*80}")
            print(f"LENGTH 30-42 ANALYSIS - {len(length_30_42_failures)} failures found")
            print(f"{'='*80}")
            
            for protein in length_30_42_failures[:10]:  # Investigate first 10
                investigate_domain_summary_generation_process(protein)
        
        else:
            logger.info(f"Retrieving domain summary failures from batches {batch_ids}...")
            failures = get_domain_summary_failures(conn, batch_ids, args.sample_size)
            logger.info(f"Found {len(failures)} proteins with domain summary generation failures")

            if not failures:
                logger.info("No domain summary failures found")
                return

            # Investigate each failure
            investigation_results = []
            for i, protein in enumerate(failures):
                logger.info(f"Investigating failure {i+1}/{len(failures)}: {protein['source_id']}")
                result = investigate_domain_summary_generation_process(protein)
                investigation_results.append(result)

            # Summarize findings
            print(f"\n{'='*80}")
            print(f"SUMMARY OF DOMAIN SUMMARY FAILURES")
            print(f"{'='*80}")
            
            # Count hypotheses
            hypothesis_counter = Counter()
            for result in investigation_results:
                for hypothesis in result['failure_hypothesis']:
                    hypothesis_counter[hypothesis] += 1
            
            print(f"Top failure hypotheses:")
            for hypothesis, count in hypothesis_counter.most_common():
                percentage = (count / len(investigation_results)) * 100
                print(f"  {hypothesis}: {count} cases ({percentage:.1f}%)")

            # Output results
            if args.output:
                output_data = {
                    'timestamp': datetime.now().isoformat(),
                    'batch_ids': batch_ids,
                    'total_investigated': len(investigation_results),
                    'failure_hypotheses': dict(hypothesis_counter),
                    'detailed_results': investigation_results
                }
                
                with open(args.output, 'w') as f:
                    json.dump(output_data, f, indent=2, default=str)
                logger.info(f"Results written to {args.output}")

    except Exception as e:
        logger.error(f"Error during investigation: {str(e)}")
        if args.verbose:
            import traceback
            logger.error(traceback.format_exc())
        sys.exit(1)

    finally:
        conn.close()


if __name__ == "__main__":
    main()
