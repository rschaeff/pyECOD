#!/usr/bin/env python3
"""
simple_check_summaries.py - Simple check of domain summary completion
"""

import os
import sys
import json
import logging
import argparse
from typing import Dict, List, Any, Optional

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

def main():
    """Main function to check domain summary quality"""
    parser = argparse.ArgumentParser(description='Check domain summary completion')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int,
                      help='Check specific batch ID')
    parser.add_argument('--check-sample', action='store_true',
                      help='Check sample of domain files')
    parser.add_argument('--sample-size', type=int, default=3,
                      help='Number of sample files to check per batch')
    parser.add_argument('--output', type=str,
                      help='Write results to output file')
    parser.add_argument('--verbose', '-v', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose)
    logger = logging.getLogger("ecod.summary_check")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get batches to check
    batch_query = """
    SELECT 
        b.id, b.batch_name, b.total_items, b.status, b.completed_items,
        b.base_path
    FROM 
        ecod_schema.batch b
    """
    
    if args.batch_id:
        batch_query += " WHERE b.id = %s"
        batches = context.db.execute_query(batch_query, (args.batch_id,))
    else:
        batch_query += " ORDER BY b.id"
        batches = context.db.execute_query(batch_query)
    
    if not batches:
        logger.warning("No batches found")
        return 1
    
    logger.info(f"Found {len(batches)} batches")
    
    # Output file
    output_data = []
    
    # Check each batch
    completed_batches = 0
    indexed_batches = 0
    created_batches = 0
    
    # Batch summary totals
    total_batches = len(batches)
    total_chains = 0
    total_summaries = 0
    
    for batch in batches:
        batch_id = batch[0]
        batch_name = batch[1]
        total_items = batch[2]
        status = batch[3]
        completed_items = batch[4]
        base_path = batch[5]
        
        # Get summary count for this batch
        summary_query = """
        SELECT 
            COUNT(DISTINCT pf.id) as summary_count
        FROM 
            ecod_schema.process_status ps
        JOIN 
            ecod_schema.batch b ON ps.batch_id = b.id
        LEFT JOIN 
            ecod_schema.process_file pf ON ps.id = pf.process_id
        WHERE 
            b.id = %s
            AND pf.file_type = 'domain_summary'
            AND pf.file_exists = TRUE
        """
        
        summary_result = context.db.execute_query(summary_query, (batch_id,))
        summary_count = summary_result[0][0] if summary_result else 0
        
        # Update totals
        total_chains += total_items
        total_summaries += summary_count
        
        # Calculate completion percentage
        completion_pct = (summary_count / total_items) * 100 if total_items > 0 else 0
        
        # Count by status
        if status == 'completed':
            completed_batches += 1
        elif status == 'indexed':
            indexed_batches += 1
        elif status == 'created':
            created_batches += 1
            
        # Batch data for output
        batch_data = {
            'id': batch_id,
            'name': batch_name,
            'status': status,
            'total_items': total_items,
            'completed_items': completed_items,
            'summary_count': summary_count,
            'completion_pct': completion_pct
        }
        
        # Check sample files if requested
        if args.check_sample and summary_count > 0:
            # Get sample domain summary files
            sample_query = """
            SELECT 
                p.pdb_id, 
                p.chain_id,
                pf.file_path
            FROM 
                ecod_schema.process_file pf
            JOIN 
                ecod_schema.process_status ps ON pf.process_id = ps.id
            JOIN 
                ecod_schema.protein p ON ps.protein_id = p.id
            WHERE 
                ps.batch_id = %s
                AND pf.file_type = 'domain_summary'
                AND pf.file_exists = TRUE
            ORDER BY 
                RANDOM()
            LIMIT %s
            """
            
            sample_files = context.db.execute_query(sample_query, (batch_id, args.sample_size))
            
            # Check each sample file
            sample_data = []
            for sample in sample_files:
                pdb_id = sample[0]
                chain_id = sample[1]
                file_path = sample[2]
                
                full_path = os.path.join(base_path, file_path)
                
                # Try to analyze the file
                try:
                    if os.path.exists(full_path):
                        with open(full_path, 'r') as f:
                            content = json.load(f)
                            
                        # Basic file stats
                        file_size = os.path.getsize(full_path)
                        domains = content.get('domains', [])
                        blast_results = content.get('blast_results', [])
                        errors = content.get('errors', [])
                        
                        sample_data.append({
                            'pdb_id': pdb_id,
                            'chain_id': chain_id,
                            'file_size': file_size,
                            'domain_count': len(domains),
                            'blast_count': len(blast_results),
                            'has_errors': len(errors) > 0,
                            'valid_json': True
                        })
                        
                        if args.verbose:
                            logger.info(f"  File {pdb_id}_{chain_id}: {len(domains)} domains, {len(blast_results)} BLAST hits")
                    else:
                        sample_data.append({
                            'pdb_id': pdb_id,
                            'chain_id': chain_id,
                            'error': 'File not found'
                        })
                        
                        if args.verbose:
                            logger.warning(f"  File {pdb_id}_{chain_id}: Not found at {full_path}")
                            
                except json.JSONDecodeError:
                    sample_data.append({
                        'pdb_id': pdb_id,
                        'chain_id': chain_id,
                        'error': 'Invalid JSON'
                    })
                    
                    if args.verbose:
                        logger.warning(f"  File {pdb_id}_{chain_id}: Invalid JSON")
                        
                except Exception as e:
                    sample_data.append({
                        'pdb_id': pdb_id,
                        'chain_id': chain_id,
                        'error': str(e)
                    })
                    
                    if args.verbose:
                        logger.warning(f"  File {pdb_id}_{chain_id}: Error - {str(e)}")
            
            # Add sample data to batch data
            batch_data['samples'] = sample_data
            
            # Calculate sample statistics
            if sample_data:
                valid_files = sum(1 for s in sample_data if s.get('valid_json', False))
                files_with_domains = sum(1 for s in sample_data if s.get('domain_count', 0) > 0)
                domain_counts = [s.get('domain_count', 0) for s in sample_data if s.get('valid_json', False)]
                
                if domain_counts:
                    avg_domains = sum(domain_counts) / len(domain_counts)
                    batch_data['avg_domains_per_chain'] = avg_domains
                
                batch_data['valid_files'] = valid_files
                batch_data['files_with_domains'] = files_with_domains
        
        # Report results
        status_symbol = "✅" if status == "completed" else "⏳" if status == "indexed" else "🆕"
        summary_symbol = "✅" if summary_count >= total_items else "⚠️" if summary_count > 0 else "❌"
        
        logger.info(f"{status_symbol} Batch {batch_id} ({batch_name}): {status}")
        logger.info(f"  {summary_symbol} Domain summaries: {summary_count}/{total_items} ({completion_pct:.1f}%)")
        
        if 'avg_domains_per_chain' in batch_data:
            logger.info(f"  Avg domains per chain: {batch_data['avg_domains_per_chain']:.2f}")
        
        output_data.append(batch_data)
    
    # Overall summary
    logger.info("\nOverall Summary:")
    logger.info(f"Total batches: {total_batches}")
    logger.info(f"  Completed: {completed_batches}")
    logger.info(f"  Indexed: {indexed_batches}")
    logger.info(f"  Created: {created_batches}")
    
    # Domain summary stats
    if total_chains > 0:
        overall_pct = (total_summaries / total_chains) * 100
        logger.info(f"Total domain summaries: {total_summaries}/{total_chains} ({overall_pct:.1f}%)")
    
    # Write output file if requested
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(output_data, f, indent=2)
        logger.info(f"Results written to {args.output}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())