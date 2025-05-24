#!/usr/bin/env python3
"""
Characterize Domain Partitioning Failures

This script analyzes proteins that consistently fail domain partitioning to identify
systematic patterns and potential algorithm improvements.

Usage:
    python scripts/characterize_failure_modes.py --config config.yml [options]
"""

import os
import sys
import logging
import argparse
import json
import psycopg2
from psycopg2.extras import RealDictCursor
import yaml
from datetime import datetime
from typing import List, Dict, Any, Optional
from collections import defaultdict, Counter
import statistics


def setup_logging(verbose=False):
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    format_str = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=level, format=format_str)
    return logging.getLogger(__name__)


def parse_config(config_path):
    """Parse configuration file."""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def get_db_connection(config):
    """Create database connection from config."""
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


def get_failed_proteins(conn, batch_ids=None):
    """Get proteins that consistently fail domain partitioning."""
    
    batch_filter = ""
    if batch_ids:
        batch_filter = f"AND ps.batch_id IN ({','.join(map(str, batch_ids))})"
    
    query = f"""
    SELECT 
        ep.id,
        ep.source_id,
        ep.pdb_id,
        ep.chain_id,
        ep.name,
        ep.type,
        ep.tax_id,
        ep.length,
        ps.batch_id,
        ps.current_stage,
        ps.status,
        ps.error_message,
        ps.is_representative,
        ps.updated_at,
        
        -- Check for various file types
        COUNT(CASE WHEN pf.file_type = 'fasta' AND pf.file_exists = true THEN 1 END) as has_fasta,
        COUNT(CASE WHEN pf.file_type LIKE '%blast%' AND pf.file_exists = true THEN 1 END) as has_blast,
        COUNT(CASE WHEN pf.file_type LIKE '%hhsearch%' AND pf.file_exists = true THEN 1 END) as has_hhsearch,
        COUNT(CASE WHEN pf.file_type = 'domain_summary' AND pf.file_exists = true THEN 1 END) as has_summary,
        COUNT(CASE WHEN pf.file_type LIKE '%partition%' AND pf.file_exists = true THEN 1 END) as has_partition,
        
        -- Get sample file paths for inspection
        array_agg(DISTINCT pf.file_path) FILTER (WHERE pf.file_type = 'domain_summary') as summary_paths,
        array_agg(DISTINCT pf.file_path) FILTER (WHERE pf.file_type LIKE '%partition%') as partition_paths
        
    FROM ecod_schema.process_status ps
    JOIN ecod_schema.protein ep ON ps.protein_id = ep.id
    LEFT JOIN ecod_schema.process_file pf ON ps.id = pf.process_id
    LEFT JOIN pdb_analysis.partition_proteins pp ON (
        ep.source_id = (pp.pdb_id || '_' || pp.chain_id)
        AND pp.batch_id = ps.batch_id
    )
    WHERE ps.status = 'error' 
       OR (ps.current_stage = 'domain_partition_failed')
       OR (pp.id IS NULL AND ps.current_stage != 'initialized')
       {batch_filter}
    GROUP BY ep.id, ep.source_id, ep.pdb_id, ep.chain_id, ep.name, ep.type, 
             ep.tax_id, ep.length, ps.batch_id, ps.current_stage, ps.status, 
             ps.error_message, ps.is_representative, ps.updated_at
    ORDER BY ep.length, ep.pdb_id, ep.chain_id
    """
    
    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(query)
        return cur.fetchall()


def analyze_length_patterns(failed_proteins):
    """Analyze length patterns in failed proteins."""
    lengths = [p['length'] for p in failed_proteins if p['length']]
    
    if not lengths:
        return {}
    
    return {
        'total_failed': len(lengths),
        'length_stats': {
            'min': min(lengths),
            'max': max(lengths),
            'mean': statistics.mean(lengths),
            'median': statistics.median(lengths),
            'std_dev': statistics.stdev(lengths) if len(lengths) > 1 else 0
        },
        'length_distribution': {
            'peptide_size (<30)': len([l for l in lengths if l < 30]),
            'very_short (30-50)': len([l for l in lengths if 30 <= l < 50]),
            'short (50-100)': len([l for l in lengths if 50 <= l < 100]),
            'medium (100-300)': len([l for l in lengths if 100 <= l < 300]),
            'long (300-500)': len([l for l in lengths if 300 <= l < 500]),
            'very_long (>500)': len([l for l in lengths if l >= 500])
        }
    }


def analyze_error_patterns(failed_proteins):
    """Analyze error message patterns."""
    error_counter = Counter()
    stage_counter = Counter()
    
    for protein in failed_proteins:
        if protein['error_message']:
            # Categorize error messages
            error = protein['error_message'].lower()
            
            if 'timeout' in error:
                error_counter['timeout'] += 1
            elif 'memory' in error or 'oom' in error:
                error_counter['memory_error'] += 1
            elif 'parse' in error or 'xml' in error:
                error_counter['parsing_error'] += 1
            elif 'file' in error and 'not found' in error:
                error_counter['file_not_found'] += 1
            elif 'permission' in error:
                error_counter['permission_error'] += 1
            elif 'blast' in error:
                error_counter['blast_error'] += 1
            elif 'hhsearch' in error or 'hhblits' in error:
                error_counter['hhsearch_error'] += 1
            else:
                error_counter['other_error'] += 1
        
        stage_counter[protein['current_stage']] += 1
    
    return {
        'error_categories': dict(error_counter),
        'failure_stages': dict(stage_counter)
    }


def analyze_file_availability(failed_proteins):
    """Analyze which processing steps completed before failure."""
    file_patterns = defaultdict(int)
    
    for protein in failed_proteins:
        pattern = []
        if protein['has_fasta'] > 0:
            pattern.append('fasta')
        if protein['has_blast'] > 0:
            pattern.append('blast')
        if protein['has_hhsearch'] > 0:
            pattern.append('hhsearch')
        if protein['has_summary'] > 0:
            pattern.append('summary')
        if protein['has_partition'] > 0:
            pattern.append('partition')
        
        file_patterns['→'.join(pattern) if pattern else 'no_files'] += 1
    
    return dict(file_patterns)


def analyze_taxonomy_patterns(failed_proteins):
    """Analyze taxonomic patterns in failures."""
    tax_counter = Counter()
    
    for protein in failed_proteins:
        if protein['tax_id']:
            tax_counter[protein['tax_id']] += 1
        else:
            tax_counter['unknown'] += 1
    
    return {
        'top_tax_ids': dict(tax_counter.most_common(20)),
        'total_unique_tax_ids': len(tax_counter)
    }


def analyze_representative_bias(failed_proteins):
    """Check if failures are biased toward representative/non-representative proteins."""
    rep_counter = Counter()
    
    for protein in failed_proteins:
        rep_counter[protein['is_representative']] += 1
    
    return dict(rep_counter)


def get_successful_comparison(conn, batch_ids=None):
    """Get comparable data from successful proteins for context."""
    
    batch_filter = ""
    if batch_ids:
        batch_filter = f"AND ps.batch_id IN ({','.join(map(str, batch_ids))})"
    
    query = f"""
    SELECT 
        ep.length,
        ep.tax_id,
        ps.is_representative
    FROM ecod_schema.process_status ps
    JOIN ecod_schema.protein ep ON ps.protein_id = ep.id
    JOIN pdb_analysis.partition_proteins pp ON (
        ep.source_id = (pp.pdb_id || '_' || pp.chain_id)
        AND pp.batch_id = ps.batch_id
    )
    WHERE ps.status != 'error' 
       AND ps.current_stage NOT IN ('domain_partition_failed', 'error')
       AND pp.is_classified = true
       {batch_filter}
    LIMIT 5000
    """
    
    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(query)
        return cur.fetchall()


def generate_sample_investigation_list(failed_proteins, sample_size=20):
    """Generate a list of specific proteins for manual investigation."""
    
    # Categorize failures for sampling
    categories = {
        'peptide_length': [p for p in failed_proteins if p['length'] and p['length'] < 30],
        'short_proteins': [p for p in failed_proteins if p['length'] and 30 <= p['length'] < 100],
        'medium_proteins': [p for p in failed_proteins if p['length'] and 100 <= p['length'] < 300],
        'long_proteins': [p for p in failed_proteins if p['length'] and p['length'] >= 300],
        'no_blast': [p for p in failed_proteins if p['has_blast'] == 0],
        'no_summary': [p for p in failed_proteins if p['has_summary'] == 0 and p['has_blast'] > 0]
    }
    
    samples = {}
    for category, proteins in categories.items():
        if proteins:
            sample_count = min(sample_size // len(categories), len(proteins))
            samples[category] = [
                {
                    'source_id': p['source_id'],
                    'pdb_id': p['pdb_id'],
                    'chain_id': p['chain_id'],
                    'length': p['length'],
                    'batch_id': p['batch_id'],
                    'error_message': p['error_message'],
                    'summary_paths': p['summary_paths']
                }
                for p in proteins[:sample_count]
            ]
    
    return samples


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Characterize domain partitioning failures')
    parser.add_argument('--config', required=True, help='Path to configuration file')
    parser.add_argument('--batch-ids', help='Comma-separated batch IDs to analyze (default: all problematic batches)')
    parser.add_argument('--output', help='Output JSON file path')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')

    args = parser.parse_args()

    # Set up logging
    logger = setup_logging(args.verbose)

    # Parse config
    try:
        config = parse_config(args.config)
    except Exception as e:
        logger.error(f"Error parsing config file: {str(e)}")
        sys.exit(1)

    # Parse batch IDs
    batch_ids = None
    if args.batch_ids:
        try:
            batch_ids = [int(bid.strip()) for bid in args.batch_ids.split(',')]
            logger.info(f"Analyzing batch IDs: {batch_ids}")
        except ValueError:
            logger.error("Invalid batch IDs format. Use comma-separated integers.")
            sys.exit(1)
    else:
        # Default problematic batches from the analysis
        batch_ids = [18, 19, 20, 21, 22, 23, 25, 26, 30, 33, 43, 50, 51]
        logger.info(f"Using default problematic batch IDs: {batch_ids}")

    # Get database connection
    try:
        conn = get_db_connection(config)
    except Exception as e:
        logger.error(f"Error connecting to database: {str(e)}")
        sys.exit(1)

    try:
        # Get failed proteins
        logger.info("Retrieving failed proteins...")
        failed_proteins = get_failed_proteins(conn, batch_ids)
        logger.info(f"Found {len(failed_proteins)} failed proteins")

        if not failed_proteins:
            logger.info("No failed proteins found")
            return

        # Get successful proteins for comparison
        logger.info("Retrieving successful proteins for comparison...")
        successful_proteins = get_successful_comparison(conn, batch_ids)
        logger.info(f"Retrieved {len(successful_proteins)} successful proteins for comparison")

        # Perform analyses
        logger.info("Analyzing failure patterns...")
        
        analysis_results = {
            'metadata': {
                'timestamp': datetime.now().isoformat(),
                'total_failed_proteins': len(failed_proteins),
                'total_successful_proteins': len(successful_proteins),
                'analyzed_batches': batch_ids
            },
            'length_analysis': analyze_length_patterns(failed_proteins),
            'error_analysis': analyze_error_patterns(failed_proteins),
            'file_availability': analyze_file_availability(failed_proteins),
            'taxonomy_patterns': analyze_taxonomy_patterns(failed_proteins),
            'representative_bias': analyze_representative_bias(failed_proteins),
            'sample_proteins': generate_sample_investigation_list(failed_proteins)
        }

        # Add comparison with successful proteins
        if successful_proteins:
            successful_lengths = [p['length'] for p in successful_proteins if p['length']]
            analysis_results['comparison'] = {
                'successful_length_stats': {
                    'min': min(successful_lengths),
                    'max': max(successful_lengths),
                    'mean': statistics.mean(successful_lengths),
                    'median': statistics.median(successful_lengths)
                } if successful_lengths else {},
                'successful_representative_distribution': dict(Counter(p['is_representative'] for p in successful_proteins))
            }

        # Output results
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(analysis_results, f, indent=2, default=str)
            logger.info(f"Analysis results written to {args.output}")
        else:
            print(json.dumps(analysis_results, indent=2, default=str))

        # Print summary
        print("\n" + "="*80)
        print("FAILURE CHARACTERIZATION SUMMARY")
        print("="*80)
        
        print(f"\nTotal Failed Proteins: {len(failed_proteins)}")
        
        if analysis_results['length_analysis']:
            length_dist = analysis_results['length_analysis']['length_distribution']
            print(f"\nLength Distribution of Failures:")
            for category, count in length_dist.items():
                percentage = (count / len(failed_proteins)) * 100
                print(f"  {category}: {count} ({percentage:.1f}%)")
        
        print(f"\nFile Completion Patterns:")
        for pattern, count in analysis_results['file_availability'].items():
            percentage = (count / len(failed_proteins)) * 100
            print(f"  {pattern}: {count} ({percentage:.1f}%)")
        
        print(f"\nError Categories:")
        for category, count in analysis_results['error_analysis']['error_categories'].items():
            percentage = (count / len(failed_proteins)) * 100
            print(f"  {category}: {count} ({percentage:.1f}%)")

    except Exception as e:
        logger.error(f"Error during analysis: {str(e)}")
        if args.verbose:
            import traceback
            logger.error(traceback.format_exc())
        sys.exit(1)

    finally:
        conn.close()


if __name__ == "__main__":
    main()
