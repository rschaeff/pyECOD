#!/usr/bin/env python3
"""
analyze_batch_pdb_ids.py - Analyze PDB IDs in batches

This script examines PDB IDs in each batch to estimate their age distribution 
and processing status.
"""

import os
import sys
import logging
import argparse
import re
from typing import Dict, Any, Optional, List, Tuple
from collections import defaultdict

# Add parent directory to path to allow imports
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

def analyze_pdb_ids(context: Any, batch_ids: Optional[List[int]] = None) -> List[Dict[str, Any]]:
    """
    Analyze PDB IDs for proteins in batches
    
    Args:
        context: Application context
        batch_ids: Optional list of batch IDs to analyze (None for all)
        
    Returns:
        List of batch information dictionaries with PDB ID statistics
    """
    logger = logging.getLogger("ecod.analyze_batch_pdb_ids")
    
    # Query to get all batches if none specified
    if not batch_ids:
        batch_query = """
        SELECT id, batch_name, type, total_items
        FROM ecod_schema.batch
        ORDER BY id
        """
        batches = context.db.execute_query(batch_query)
    else:
        # Query for specific batches
        batch_query = """
        SELECT id, batch_name, type, total_items
        FROM ecod_schema.batch
        WHERE id IN %s
        ORDER BY id
        """
        batches = context.db.execute_query(batch_query, (tuple(batch_ids),))
    
    batch_results = []
    
    for batch in batches:
        batch_id = batch[0]
        batch_name = batch[1]
        batch_type = batch[2]
        total_items = batch[3]
        
        logger.info(f"Analyzing batch: {batch_name} (ID: {batch_id})")
        
        # Query PDB IDs in this batch
        id_query = """
        SELECT 
            p.pdb_id
        FROM 
            ecod_schema.process_status ps
        JOIN 
            ecod_schema.protein p ON ps.protein_id = p.id
        WHERE 
            ps.batch_id = %s
        """
        
        id_results = context.db.execute_query(id_query, (batch_id,))
        
        if not id_results:
            logger.warning(f"No PDB IDs found for batch {batch_id}")
            continue
        
        # Extract PDB IDs
        pdb_ids = [result[0] for result in id_results]
        
        # Analyze PDB ID distributions
        id_stats = analyze_pdb_id_distribution(pdb_ids)
        
        # Check for domain summaries
        summary_query = """
        SELECT COUNT(DISTINCT ps.protein_id)
        FROM ecod_schema.process_status ps
        JOIN ecod_schema.process_file pf ON ps.id = pf.process_id
        WHERE ps.batch_id = %s AND pf.file_type = 'domain_summary'
        """
        
        summary_result = context.db.execute_query(summary_query, (batch_id,))
        summary_count = summary_result[0][0] if summary_result else 0
        
        # Calculate completion percentage
        completion_pct = (summary_count / total_items * 100) if total_items > 0 else 0
        
        # Find some proteins without domain summaries for testing
        test_candidates_query = """
        SELECT 
            p.id, p.pdb_id, p.chain_id, p.source_id
        FROM 
            ecod_schema.protein p
        JOIN 
            ecod_schema.process_status ps ON p.id = ps.protein_id
        LEFT JOIN (
            SELECT DISTINCT process_id
            FROM ecod_schema.process_file
            WHERE file_type = 'domain_summary'
        ) ds ON ps.id = ds.process_id
        WHERE 
            ps.batch_id = %s
            AND ds.process_id IS NULL
        LIMIT 5
        """
        
        test_candidates = context.db.execute_query(test_candidates_query, (batch_id,))
        
        # Add to results
        batch_results.append({
            "id": batch_id,
            "name": batch_name,
            "type": batch_type,
            "total_items": total_items,
            "id_stats": id_stats,
            "summary_count": summary_count,
            "completion_pct": completion_pct,
            "test_candidates": test_candidates
        })
    
    return batch_results

def analyze_pdb_id_distribution(pdb_ids: List[str]) -> Dict[str, Any]:
    """
    Analyze distribution of PDB IDs
    
    Args:
        pdb_ids: List of PDB IDs
        
    Returns:
        Dictionary with distribution statistics
    """
    # Initialize counters
    id_ranges = {
        "1xxx": 0,  # Very old (1990s-early 2000s)
        "2xxx": 0,  # 2000s
        "3xxx": 0,  # Early 2010s
        "4xxx": 0,  # Mid 2010s
        "5xxx": 0,  # Late 2010s
        "6xxx": 0,  # Early 2020s
        "7xxx": 0,  # 2020s
        "8xxx": 0,  # Newest
        "other": 0  # Non-standard format
    }
    
    # Count PDB IDs by prefix
    for pdb_id in pdb_ids:
        if not pdb_id or len(pdb_id) < 4:
            id_ranges["other"] += 1
            continue
            
        prefix = pdb_id[0].lower()
        
        if prefix == '1':
            id_ranges["1xxx"] += 1
        elif prefix == '2':
            id_ranges["2xxx"] += 1
        elif prefix == '3':
            id_ranges["3xxx"] += 1
        elif prefix == '4':
            id_ranges["4xxx"] += 1
        elif prefix == '5':
            id_ranges["5xxx"] += 1
        elif prefix == '6':
            id_ranges["6xxx"] += 1
        elif prefix == '7':
            id_ranges["7xxx"] += 1
        elif prefix == '8':
            id_ranges["8xxx"] += 1
        else:
            id_ranges["other"] += 1
    
    # Calculate percentages
    total = len(pdb_ids)
    percentages = {}
    
    for range_name, count in id_ranges.items():
        percentages[range_name] = (count / total * 100) if total > 0 else 0
    
    return {
        "counts": id_ranges,
        "percentages": percentages,
        "total": total,
        "recent_count": id_ranges["7xxx"] + id_ranges["8xxx"],
        "recent_percent": (id_ranges["7xxx"] + id_ranges["8xxx"]) / total * 100 if total > 0 else 0
    }

def main():
    """Main function to analyze PDB IDs in batches"""
    parser = argparse.ArgumentParser(description='Analyze PDB IDs in batches')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-ids', type=int, nargs='+',
                      help='Specific batch IDs to analyze (omit for all)')
    parser.add_argument('--min-recent', type=int, default=0,
                      help='Minimum percentage of recent PDB IDs (7xxx-8xxx)')
    parser.add_argument('--max-count', type=int, default=20,
                      help='Maximum number of batches to display')
    parser.add_argument('--output', type=str,
                      help='Output file path for CSV report')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger("ecod.analyze_batch_pdb_ids")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Analyze batches
    batch_results = analyze_pdb_ids(context, args.batch_ids)
    
    # Filter by minimum recent percentage
    if args.min_recent > 0:
        batch_results = [b for b in batch_results if b["id_stats"]["recent_percent"] >= args.min_recent]
    
    # Sort by recency (percentage of 7xxx and 8xxx IDs)
    batch_results.sort(key=lambda b: b["id_stats"]["recent_percent"], reverse=True)
    
    # Limit to maximum count
    if args.max_count and len(batch_results) > args.max_count:
        batch_results = batch_results[:args.max_count]
    
    # Print results
    print("\nBatch PDB ID Analysis Report")
    print("===========================\n")
    print(f"{'ID':<5} {'Batch Name':<30} {'Type':<10} {'Total':<8} {'Recent %':<8} {'Completion':<10}")
    print(f"{'-'*5} {'-'*30} {'-'*10} {'-'*8} {'-'*8} {'-'*10}")
    
    for batch in batch_results:
        recent_pct = batch["id_stats"]["recent_percent"]
        
        print(f"{batch['id']:<5} {batch['name'][:30]:<30} {batch['type']:<10} {batch['total_items']:<8} "
              f"{recent_pct:>7.1f}% {batch['completion_pct']:>9.1f}%")
    
    print("\nPDB ID Distribution:")
    print(f"{'ID':<5} {'Batch Name':<25} {'1xxx':<7} {'2xxx':<7} {'3xxx':<7} {'4xxx':<7} {'5xxx':<7} {'6xxx':<7} {'7xxx':<7} {'8xxx':<7}")
    print(f"{'-'*5} {'-'*25} {'-'*7} {'-'*7} {'-'*7} {'-'*7} {'-'*7} {'-'*7} {'-'*7} {'-'*7}")
    
    for batch in batch_results:
        counts = batch["id_stats"]["counts"]
        
        print(f"{batch['id']:<5} {batch['name'][:25]:<25} "
              f"{counts['1xxx']:>6} {counts['2xxx']:>7} {counts['3xxx']:>7} {counts['4xxx']:>7} "
              f"{counts['5xxx']:>7} {counts['6xxx']:>7} {counts['7xxx']:>7} {counts['8xxx']:>7}")
    
    # Write CSV if output path provided
    if args.output:
        try:
            with open(args.output, 'w') as f:
                # Write header
                f.write("batch_id,batch_name,type,total_items,summary_count,completion_pct,")
                f.write("recent_count,recent_percent,")
                f.write("1xxx,2xxx,3xxx,4xxx,5xxx,6xxx,7xxx,8xxx,other\n")
                
                # Write data
                for batch in batch_results:
                    counts = batch["id_stats"]["counts"]
                    
                    f.write(f"{batch['id']},{batch['name']},{batch['type']},{batch['total_items']},")
                    f.write(f"{batch['summary_count']},{batch['completion_pct']:.1f},")
                    f.write(f"{batch['id_stats']['recent_count']},{batch['id_stats']['recent_percent']:.1f},")
                    f.write(f"{counts['1xxx']},{counts['2xxx']},{counts['3xxx']},{counts['4xxx']},")
                    f.write(f"{counts['5xxx']},{counts['6xxx']},{counts['7xxx']},{counts['8xxx']},{counts['other']}\n")
                
                logger.info(f"Report written to {args.output}")
        except Exception as e:
            logger.error(f"Error writing CSV report: {str(e)}")
    
    # Recommend batches
    print("\nRecommended Batches for Testing:")
    recommended = [b for b in batch_results if b["id_stats"]["recent_count"] > 0 and b["completion_pct"] < 100]
    recommended.sort(key=lambda b: (b["id_stats"]["recent_count"], -b["completion_pct"]), reverse=True)
    
    for i, batch in enumerate(recommended[:5]):
        recent_count = batch["id_stats"]["recent_count"]
        recent_pct = batch["id_stats"]["recent_percent"]
        
        print(f"{i+1}. Batch {batch['id']} ({batch['name']}): {recent_count} recent structures ({recent_pct:.1f}%), {batch['completion_pct']:.1f}% complete")
        
        # Show example test candidates
        if batch["test_candidates"]:
            print("   Example proteins for testing:")
            for j, candidate in enumerate(batch["test_candidates"]):
                protein_id = candidate[0]
                pdb_id = candidate[1]
                chain_id = candidate[2]
                source_id = candidate[3]
                
                print(f"     {j+1}. Protein ID: {protein_id}, PDB: {pdb_id}_{chain_id}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())