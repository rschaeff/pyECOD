#!/usr/bin/env python3
"""
analyze_batch_dates.py - Analyze deposition dates in batches

This script examines proteins in each batch to determine the earliest and latest
PDB deposition dates, helping identify batches with newer structures.
"""

import os
import sys
import logging
import argparse
from typing import Dict, Any, Optional, List
from datetime import datetime

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

def analyze_batch_dates(context: Any, batch_ids: Optional[List[int]] = None) -> List[Dict[str, Any]]:
    """
    Analyze deposition dates for proteins in batches
    
    Args:
        context: Application context
        batch_ids: Optional list of batch IDs to analyze (None for all)
        
    Returns:
        List of batch information dictionaries with date statistics
    """
    logger = logging.getLogger("ecod.analyze_batch_dates")
    
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
        
        # Query PDB deposition dates through pdb_analysis schema
        date_query = """
        SELECT 
            ps.deposition_date,
            COUNT(*) as structure_count
        FROM 
            ecod_schema.process_status es
        JOIN 
            ecod_schema.protein ep ON es.protein_id = ep.id
        LEFT JOIN 
            pdb_analysis.protein p ON ep.source_id = p.source_id
        LEFT JOIN 
            pdb_analysis.protein_structure ps ON p.id = ps.protein_id
        WHERE 
            es.batch_id = %s
            AND ps.deposition_date IS NOT NULL
        GROUP BY 
            ps.deposition_date
        ORDER BY 
            ps.deposition_date
        """
        
        date_results = context.db.execute_query(date_query, (batch_id,))
        
        if not date_results:
            logger.warning(f"No deposition date information found for batch {batch_id}")
            continue
        
        # Extract dates and counts
        dates = [result[0] for result in date_results]
        counts = [result[1] for result in date_results]
        
        # Calculate statistics
        earliest_date = min(dates) if dates else None
        latest_date = max(dates) if dates else None
        
        # Count proteins by date ranges
        date_ranges = {
            "pre_2010": 0,
            "2010_2015": 0,
            "2016_2020": 0,
            "2021_present": 0
        }
        
        for date, count in zip(dates, counts):
            year = date.year
            if year < 2010:
                date_ranges["pre_2010"] += count
            elif year < 2016:
                date_ranges["2010_2015"] += count
            elif year < 2021:
                date_ranges["2016_2020"] += count
            else:
                date_ranges["2021_present"] += count
        
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
        
        # Add to results
        batch_results.append({
            "id": batch_id,
            "name": batch_name,
            "type": batch_type,
            "total_items": total_items,
            "earliest_date": earliest_date,
            "latest_date": latest_date,
            "date_ranges": date_ranges,
            "summary_count": summary_count,
            "completion_pct": completion_pct
        })
    
    return batch_results

def main():
    """Main function to analyze batch dates"""
    parser = argparse.ArgumentParser(description='Analyze deposition dates in batches')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-ids', type=int, nargs='+',
                      help='Specific batch IDs to analyze (omit for all)')
    parser.add_argument('--min-date', type=str,
                      help='Minimum deposition date to consider (YYYY-MM-DD)')
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
    logger = logging.getLogger("ecod.analyze_batch_dates")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Set minimum date if provided
    min_date = None
    if args.min_date:
        try:
            min_date = datetime.strptime(args.min_date, '%Y-%m-%d').date()
        except ValueError:
            logger.error(f"Invalid date format: {args.min_date}. Please use YYYY-MM-DD.")
            return 1
    
    # Analyze batches
    batch_results = analyze_batch_dates(context, args.batch_ids)
    
    # Filter by minimum date if provided
    if min_date:
        batch_results = [b for b in batch_results if b["latest_date"] and b["latest_date"] >= min_date]
    
    # Sort by latest deposition date (newest first)
    batch_results.sort(key=lambda b: b["latest_date"] if b["latest_date"] else datetime.min.date(), reverse=True)
    
    # Limit to maximum count
    if args.max_count and len(batch_results) > args.max_count:
        batch_results = batch_results[:args.max_count]
    
    # Print results
    print("\nBatch Date Analysis Report")
    print("=========================\n")
    print(f"{'ID':<5} {'Batch Name':<30} {'Type':<10} {'Total':<8} {'Earliest':<12} {'Latest':<12} {'Completion':<10} {'Recent':<8}")
    print(f"{'-'*5} {'-'*30} {'-'*10} {'-'*8} {'-'*12} {'-'*12} {'-'*10} {'-'*8}")
    
    for batch in batch_results:
        earliest = batch["earliest_date"].strftime('%Y-%m-%d') if batch["earliest_date"] else 'N/A'
        latest = batch["latest_date"].strftime('%Y-%m-%d') if batch["latest_date"] else 'N/A'
        recent_count = batch["date_ranges"]["2021_present"]
        
        print(f"{batch['id']:<5} {batch['name'][:30]:<30} {batch['type']:<10} {batch['total_items']:<8} "
              f"{earliest:<12} {latest:<12} {batch['completion_pct']:>8.1f}% {recent_count:>8}")
    
    print("\nDate Range Distribution:")
    print(f"{'ID':<5} {'Batch Name':<30} {'Pre-2010':<10} {'2010-2015':<10} {'2016-2020':<10} {'2021+':<10}")
    print(f"{'-'*5} {'-'*30} {'-'*10} {'-'*10} {'-'*10} {'-'*10}")
    
    for batch in batch_results:
        ranges = batch["date_ranges"]
        print(f"{batch['id']:<5} {batch['name'][:30]:<30} {ranges['pre_2010']:>10} {ranges['2010_2015']:>10} "
              f"{ranges['2016_2020']:>10} {ranges['2021_present']:>10}")
    
    # Write CSV if output path provided
    if args.output:
        try:
            with open(args.output, 'w') as f:
                # Write header
                f.write("batch_id,batch_name,type,total_items,earliest_date,latest_date,summary_count,completion_pct,")
                f.write("pre_2010,2010_2015,2016_2020,2021_present\n")
                
                # Write data
                for batch in batch_results:
                    earliest = batch["earliest_date"].strftime('%Y-%m-%d') if batch["earliest_date"] else 'N/A'
                    latest = batch["latest_date"].strftime('%Y-%m-%d') if batch["latest_date"] else 'N/A'
                    ranges = batch["date_ranges"]
                    
                    f.write(f"{batch['id']},{batch['name']},{batch['type']},{batch['total_items']},{earliest},{latest},")
                    f.write(f"{batch['summary_count']},{batch['completion_pct']:.1f},")
                    f.write(f"{ranges['pre_2010']},{ranges['2010_2015']},{ranges['2016_2020']},{ranges['2021_present']}\n")
                
                logger.info(f"Report written to {args.output}")
        except Exception as e:
            logger.error(f"Error writing CSV report: {str(e)}")
    
    # Recommend batches
    print("\nRecommended Batches for Testing:")
    recommended = [b for b in batch_results if b["date_ranges"]["2021_present"] > 0 and b["completion_pct"] < 100]
    recommended.sort(key=lambda b: b["date_ranges"]["2021_present"], reverse=True)
    
    for i, batch in enumerate(recommended[:5]):
        print(f"{i+1}. Batch {batch['id']} ({batch['name']}): {batch['date_ranges']['2021_present']} structures from 2021+ ({batch['completion_pct']:.1f}% complete)")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())