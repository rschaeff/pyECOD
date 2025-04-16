#!/usr/bin/env python3
"""
check_batch.py - Query the status of existing ECOD batches
"""

import os
import sys
import argparse
import logging
from typing import Dict, Any, Optional

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))

from ecod.core.context import ApplicationContext
from ecod.exceptions import ECODError
from ecod.error_handlers import handle_exceptions

def setup_logging(verbose: bool = False):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

def get_batch_status(context: ApplicationContext, batch_id: Optional[int] = None) -> Dict[str, Any]:
    """Get status of batches
    
    Args:
        context: Application context
        batch_id: Optional batch ID to check specific batch
        
    Returns:
        Dictionary with batch status information
    """
    db = context.db
    
    if batch_id:
        # Query for specific batch
        query = """
        SELECT 
            b.id, b.batch_name, b.base_path, b.type, b.ref_version,
            b.total_items, b.completed_items, b.status, 
            b.created_at, b.completed_at
        FROM 
            ecod_schema.batch b
        WHERE 
            b.id = %s
        """
        batch_rows = db.execute_dict_query(query, (batch_id,))
        
        if not batch_rows:
            return {"error": f"Batch {batch_id} not found"}
        
        # Get process status breakdown
        query = """
        SELECT 
            status, current_stage, COUNT(*) as count
        FROM 
            ecod_schema.process_status
        WHERE 
            batch_id = %s
        GROUP BY 
            status, current_stage
        ORDER BY 
            status, current_stage
        """
        status_rows = db.execute_dict_query(query, (batch_id,))
        
        # Get file type counts
        query = """
        SELECT 
            pf.file_type, COUNT(*) as count
        FROM 
            ecod_schema.process_file pf
        JOIN
            ecod_schema.process_status ps ON pf.process_id = ps.id
        WHERE 
            ps.batch_id = %s
            AND pf.file_exists = TRUE
        GROUP BY 
            pf.file_type
        ORDER BY 
            pf.file_type
        """
        file_rows = db.execute_dict_query(query, (batch_id,))
        
        # Format results
        batch = batch_rows[0]
        result = {
            "batch": batch,
            "status_breakdown": status_rows,
            "file_types": file_rows
        }
        return result
    else:
        # Query all batches
        query = """
        SELECT 
            b.id, b.batch_name, b.base_path, b.type, b.ref_version,
            b.total_items, b.completed_items, b.status, 
            b.created_at, b.completed_at
        FROM 
            ecod_schema.batch b
        ORDER BY 
            b.id DESC
        """
        batch_rows = db.execute_dict_query(query)
        return {"batches": batch_rows}

@handle_exceptions(exit_on_error=True)
def main():
    parser = argparse.ArgumentParser(description='Check ECOD Batch Status')
    parser.add_argument('--config', type=str, default='config/config.yml', 
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int,
                      help='Specific batch ID to check')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose)
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get batch status
    status = get_batch_status(context, args.batch_id)
    
    # Display results
    if "error" in status:
        print(f"Error: {status['error']}")
        return 1
    
    if "batches" in status:
        # Display all batches
        print("Available batches:")
        print("="*80)
        for batch in status["batches"]:
            progress = (batch['completed_items'] / batch['total_items']) * 100 if batch['total_items'] > 0 else 0
            print(f"Batch {batch['id']}: {batch['batch_name']}")
            print(f"  Type: {batch['type']}")
            print(f"  Status: {batch['status']}")
            print(f"  Progress: {batch['completed_items']}/{batch['total_items']} ({progress:.2f}%)")
            print(f"  Created: {batch['created_at']}")
            if batch['completed_at']:
                print(f"  Completed: {batch['completed_at']}")
            print("-"*80)
    else:
        # Display specific batch details
        batch = status["batch"]
        progress = (batch['completed_items'] / batch['total_items']) * 100 if batch['total_items'] > 0 else 0
        
        print(f"Batch {batch['id']}: {batch['batch_name']}")
        print(f"  Type: {batch['type']}")
        print(f"  Reference: {batch['ref_version']}")
        print(f"  Status: {batch['status']}")
        print(f"  Progress: {batch['completed_items']}/{batch['total_items']} ({progress:.2f}%)")
        print(f"  Created: {batch['created_at']}")
        if batch['completed_at']:
            print(f"  Completed: {batch['completed_at']}")
        print(f"  Path: {batch['base_path']}")
        
        print("\nProcess Status Breakdown:")
        print("-"*80)
        for row in status["status_breakdown"]:
            print(f"  {row['status']} - {row['current_stage']}: {row['count']}")
        
        print("\nFile Types:")
        print("-"*80)
        for row in status["file_types"]:
            print(f"  {row['file_type']}: {row['count']}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())