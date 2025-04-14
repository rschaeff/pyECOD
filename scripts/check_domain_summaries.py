#!/usr/bin/env python3
"""
check_domain_summaries.py - Evaluate quality of generated domain summaries
"""

import os
import sys
import json
import logging
import argparse
from typing import Dict, List, Any, Optional, Tuple

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext

def setup_logging(verbose: bool = False, log_file: Optional[str] = None):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    
    handlers = [logging.StreamHandler()]
    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )

def get_batch_summary_stats(batch_id: int, context: ApplicationContext) -> Dict[str, Any]:
    """Get statistics about domain summaries for a batch"""
    logger = logging.getLogger("ecod.summary_check")
    
    # Query to get basic summary stats
    query = """
    SELECT 
        b.id as batch_id,
        b.batch_name,
        b.total_items,
        COUNT(DISTINCT ps.id) as total_processes,
        COUNT(DISTINCT CASE WHEN pf.file_type = 'domain_summary' AND pf.file_exists = TRUE THEN ps.id END) as summary_count,
        COUNT(DISTINCT CASE WHEN ps.status = 'success' THEN ps.id END) as success_count,
        COUNT(DISTINCT CASE WHEN ps.status = 'error' THEN ps.id END) as error_count,
        COUNT(DISTINCT CASE WHEN ps.current_stage = 'domain_summary' THEN ps.id END) as summary_stage_count
    FROM 
        ecod_schema.batch b
    JOIN 
        ecod_schema.process_status ps ON b.id = ps.batch_id
    LEFT JOIN 
        ecod_schema.process_file pf ON ps.id = pf.process_id
    WHERE 
        b.id = %s
    GROUP BY 
        b.id, b.batch_name, b.total_items
    """
    
    result = context.db.execute_dict_query(query, (batch_id,))
    
    if not result:
        logger.warning(f"No data found for batch {batch_id}")
        return {}
    
    stats = result[0]
    logger.info(f"Batch {batch_id} ({stats['batch_name']}): {stats['summary_count']}/{stats['total_items']} summaries generated")
    
    # Get sample of successful summaries for inspection
    sample_query = """
    SELECT 
        p.pdb_id, 
        p.chain_id, 
        ps.id as process_id,
        pf.file_path
    FROM 
        ecod_schema.process_status ps
    JOIN 
        ecod_schema.protein p ON ps.protein_id = p.id
    JOIN 
        ecod_schema.process_file pf ON ps.id = pf.process_id
    WHERE 
        ps.batch_id = %s
        AND pf.file_type = 'domain_summary'
        AND pf.file_exists = TRUE
    ORDER BY 
        RANDOM()
    LIMIT 5
    """
    
    sample_results = context.db.execute_dict_query(sample_query, (batch_id,))
    
    # Add sample to stats (initialize the key first)
    stats['samples'] = []
    if sample_results:
        stats['samples'] = sample_results
    
    # Analyze summary file sizes
    size_query = """
    SELECT 
        MIN(pf.file_size) as min_size,
        MAX(pf.file_size) as max_size,
        AVG(pf.file_size) as avg_size,
        PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY pf.file_size) as median_size
    FROM 
        ecod_schema.process_file pf
    JOIN 
        ecod_schema.process_status ps ON pf.process_id = ps.id
    WHERE 
        ps.batch_id = %s
        AND pf.file_type = 'domain_summary'
        AND pf.file_exists = TRUE
    """
    
    size_results = context.db.execute_dict_query(size_query, (batch_id,))
    
    if size_results:
        stats['file_sizes'] = size_results[0]
    
    return stats

def check_summary_content(file_path: str, base_path: str) -> Dict[str, Any]:
    """Check the content of a domain summary file"""
    full_path = os.path.join(base_path, file_path)
    
    if not os.path.exists(full_path):
        return {
            'exists': False,
            'error': f"File not found: {full_path}"
        }
    
    try:
        with open(full_path, 'r') as f:
            content = json.load(f)
        
        # Basic content checks
        return {
            'exists': True,
            'valid_json': True,
            'has_domains': 'domains' in content and len(content.get('domains', [])) > 0,
            'domain_count': len(content.get('domains', [])),
            'has_blast_results': 'blast_results' in content and len(content.get('blast_results', [])) > 0,
            'blast_count': len(content.get('blast_results', [])),
            'has_errors': content.get('errors', []) != []
        }
    except json.JSONDecodeError:
        return {
            'exists': True,
            'valid_json': False,
            'error': "Invalid JSON content"
        }
    except Exception as e:
        return {
            'exists': True,
            'valid_json': False,
            'error': str(e)
        }

def analyze_batch_summaries(batch_id: int, context: ApplicationContext) -> Dict[str, Any]:
    """Analyze domain summaries for a batch"""
    logger = logging.getLogger("ecod.summary_check")
    
    # Get batch statistics
    stats = get_batch_summary_stats(batch_id, context)
    
    if not stats:
        return {}
    
    # Get batch base path
    path_query = """
    SELECT base_path FROM ecod_schema.batch WHERE id = %s
    """
    path_result = context.db.execute_query(path_query, (batch_id,))
    base_path = path_result[0][0] if path_result else None
    
    if not base_path:
        logger.warning(f"Could not determine base path for batch {batch_id}")
        return stats
    
    # Initialize content_analysis list
    stats['content_analysis'] = []
    
    # Check sample summary contents
    for sample in stats.get('samples', []):
        file_path = sample['file_path']
        content_info = check_summary_content(file_path, base_path)
        
        stats['content_analysis'].append({
            'pdb_id': sample['pdb_id'],
            'chain_id': sample['chain_id'],
            'file_path': file_path,
            'content': content_info
        })
    
    # Check for any errors in the sample
    error_count = sum(1 for s in stats['content_analysis'] if not s['content'].get('valid_json', False) or s['content'].get('error', False))
    
    if error_count > 0:
        logger.warning(f"Found {error_count} errors in sample summaries")
    
    # Count domains found in sample
    domain_counts = [s['content'].get('domain_count', 0) for s in stats['content_analysis'] if s['content'].get('has_domains', False)]
    if domain_counts:
        avg_domains = sum(domain_counts) / len(domain_counts)
        logger.info(f"Average domains per chain in sample: {avg_domains:.2f}")
        stats['avg_domains_per_chain'] = avg_domains
    
    return stats

def print_summary_details(file_path: str, base_path: str):
    """Print details of a domain summary file"""
    full_path = os.path.join(base_path, file_path)
    
    if not os.path.exists(full_path):
        print(f"File not found: {full_path}")
        return
    
    try:
        with open(full_path, 'r') as f:
            content = json.load(f)
        
        print(f"Summary file: {full_path}")
        print(f"PDB ID: {content.get('pdb_id')}")
        print(f"Chain ID: {content.get('chain_id')}")
        print(f"Sequence length: {len(content.get('sequence', ''))}")
        print(f"Number of domains: {len(content.get('domains', []))}")
        print(f"Number of BLAST results: {len(content.get('blast_results', []))}")
        
        print("\nDomains:")
        for i, domain in enumerate(content.get('domains', [])):
            print(f"  Domain {i+1}:")
            print(f"    Range: {domain.get('range')}")
            print(f"    Score: {domain.get('score')}")
            print(f"    Evidence: {domain.get('evidence_type')}")
            
        if content.get('errors', []):
            print("\nErrors:")
            for error in content.get('errors', []):
                print(f"  {error}")
    
    except Exception as e:
        print(f"Error reading file: {e}")

def main():
    """Main function to check domain summary quality"""
    parser = argparse.ArgumentParser(description='Check quality of domain summaries')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int,
                      help='Check specific batch ID')
    parser.add_argument('--all-batches', action='store_true',
                      help='Check all indexed batches')
    parser.add_argument('--output', type=str,
                      help='Output file for detailed results (JSON)')
    parser.add_argument('--examine-file', type=str,
                      help='Path to specific summary file to examine in detail')
    parser.add_argument('--verbose', '-v', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose)
    logger = logging.getLogger("ecod.summary_check")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # If examining a specific file
    if args.examine_file:
        # Get base path for the specific batch
        if args.batch_id:
            path_query = """
            SELECT base_path FROM ecod_schema.batch WHERE id = %s
            """
            path_result = context.db.execute_query(path_query, (args.batch_id,))
            base_path = path_result[0][0] if path_result else '.'
            
            print_summary_details(args.examine_file, base_path)
            return 0
        else:
            logger.error("--examine-file requires --batch-id to determine base path")
            return 1
    
    # Determine which batches to check
    batch_ids = []
    
    if args.batch_id:
        batch_ids = [args.batch_id]
    elif args.all_batches:
        # Get all indexed batches
        query = """
        SELECT id FROM ecod_schema.batch WHERE status = 'indexed' AND total_items = 5000 ORDER BY id
        """
        results = context.db.execute_query(query)
        batch_ids = [row[0] for row in results]
    else:
        # Default: get batches with domain summaries
        query = """
        SELECT DISTINCT b.id
        FROM ecod_schema.batch b
        JOIN ecod_schema.process_status ps ON b.id = ps.batch_id
        JOIN ecod_schema.process_file pf ON ps.id = pf.process_id
        WHERE pf.file_type = 'domain_summary' AND pf.file_exists = TRUE
        ORDER BY b.id
        """
        results = context.db.execute_query(query)
        batch_ids = [row[0] for row in results]
    
    if not batch_ids:
        logger.warning("No batches found to check")
        return 1
    
    logger.info(f"Checking domain summaries for {len(batch_ids)} batches")
    
    # Analyze each batch
    all_results = {}
    for batch_id in batch_ids:
        logger.info(f"Analyzing batch {batch_id}")
        batch_stats = analyze_batch_summaries(batch_id, context)
        all_results[str(batch_id)] = batch_stats
    
    # Write detailed results to output file if requested
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(all_results, f, indent=2)
        logger.info(f"Detailed results written to {args.output}")
    
    # Summary of all batches
    logger.info("\nSummary of all batches:")
    logger.info("-----------------------")
    
    total_files = 0
    total_valid = 0
    total_with_domains = 0
    total_chains = 0
    
    for batch_id, stats in all_results.items():
        if not stats:
            logger.info(f"Batch {batch_id}: No data")
            continue
            
        summary_count = stats.get('summary_count', 0)
        total_items = stats.get('total_items', 0)
        percentage = (summary_count / total_items * 100) if total_items > 0 else 0
        
        logger.info(f"Batch {batch_id} ({stats.get('batch_name', 'unknown')}): {summary_count}/{total_items} summaries ({percentage:.1f}%)")
        
        # Update totals
        total_files += summary_count
        total_chains += total_items
        
        # Report on samples if available
        if 'content_analysis' in stats:
            valid_samples = sum(1 for s in stats['content_analysis'] if s['content'].get('valid_json', False))
            samples_with_domains = sum(1 for s in stats['content_analysis'] if s['content'].get('has_domains', False))
            
            total_valid += valid_samples
            total_with_domains += samples_with_domains
            
            logger.info(f"  Sample quality: {valid_samples}/{len(stats['content_analysis'])} valid, {samples_with_domains}/{len(stats['content_analysis'])} with domains")
            
            if 'avg_domains_per_chain' in stats:
                logger.info(f"  Avg domains per chain: {stats['avg_domains_per_chain']:.2f}")
    
    # Overall summary
    if total_chains > 0:
        logger.info("\nOverall Summary:")
        logger.info(f"Total domain summaries: {total_files}/{total_chains} ({total_files/total_chains*100:.1f}%)")
        
        sample_count = len(all_results) * 5  # 5 samples per batch
        if sample_count > 0:
            logger.info(f"Valid JSON: {total_valid}/{sample_count} ({total_valid/sample_count*100:.1f}%)")
            logger.info(f"With domains: {total_with_domains}/{sample_count} ({total_with_domains/sample_count*100:.1f}%)")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())