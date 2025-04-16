#!/usr/bin/env python3
"""
examine_invalid_xml.py - Examine and display invalid domain summary XML files

This script finds and displays the content of invalid XML files in a batch
to help diagnose what's making them invalid.
"""

import os
import sys
import argparse
import logging
import xml.etree.ElementTree as ET
from typing import Dict, List, Any, Optional

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

# Import application context
from ecod.core.context import ApplicationContext

def setup_logging(verbose: bool = False):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

def is_valid_xml(file_path: str) -> bool:
    """Check if a file contains valid XML"""
    try:
        # Try to parse the XML
        ET.parse(file_path)
        return True
    except ET.ParseError:
        return False
    except Exception:
        return False

def get_file_size(file_path: str) -> int:
    """Get file size in bytes"""
    try:
        return os.path.getsize(file_path)
    except:
        return -1

def examine_file(file_path: str) -> Dict[str, Any]:
    """Examine a file and return its characteristics"""
    result = {
        'file_path': file_path,
        'exists': os.path.exists(file_path),
        'size': get_file_size(file_path),
        'is_valid_xml': False,
        'content': None,
        'error': None
    }
    
    if result['exists'] and result['size'] > 0:
        # Check if it's valid XML
        result['is_valid_xml'] = is_valid_xml(file_path)
        
        # Get file content
        try:
            with open(file_path, 'r') as f:
                content = f.read()
                # Limit content length for display
                if len(content) > 2000:
                    result['content'] = content[:2000] + "... [truncated]"
                else:
                    result['content'] = content
        except Exception as e:
            result['error'] = f"Error reading file: {str(e)}"
    
    return result

def find_invalid_xml_files(context, batch_id: int, limit: int = 5) -> List[Dict[str, Any]]:
    """Find invalid XML files in a batch"""
    logger = logging.getLogger('ecod.examine')
    
    # Get batch info
    batch_query = """
    SELECT base_path FROM ecod_schema.batch WHERE id = %s
    """
    batch_result = context.db.execute_query(batch_query, (batch_id,))
    
    if not batch_result:
        logger.error(f"Batch {batch_id} not found")
        return []
    
    base_path = batch_result[0][0]
    
    # Get domain summary files
    summary_query = """
    SELECT 
        p.pdb_id, 
        p.chain_id, 
        pf.file_path,
        pf.file_exists,
        pf.file_size
    FROM 
        ecod_schema.process_file pf
    JOIN 
        ecod_schema.process_status ps ON pf.process_id = ps.id
    JOIN 
        ecod_schema.protein p ON ps.protein_id = p.id
    WHERE 
        ps.batch_id = %s
        AND pf.file_type = 'domain_summary'
    ORDER BY 
        pf.file_size ASC  -- Sort by size to find potentially empty or minimal files
    LIMIT 100  -- Get a good sample to check
    """
    
    summary_files = context.db.execute_dict_query(summary_query, (batch_id,))
    
    # Examine files
    results = []
    count_valid = 0
    count_invalid = 0
    
    for summary in summary_files:
        file_path = os.path.join(base_path, summary['file_path'])
        
        # Examine file
        file_result = examine_file(file_path)
        
        # Add protein info
        file_result['pdb_id'] = summary['pdb_id']
        file_result['chain_id'] = summary['chain_id']
        file_result['db_file_exists'] = summary['file_exists']
        file_result['db_file_size'] = summary['file_size']
        
        # Track validity
        if file_result['is_valid_xml']:
            count_valid += 1
        else:
            count_invalid += 1
            results.append(file_result)
            
            # Stop if we have enough invalid files
            if count_invalid >= limit:
                break
    
    logger.info(f"Found {count_valid} valid and {count_invalid} invalid XML files")
    return results

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Examine invalid domain summary XML files')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to analyze')
    parser.add_argument('--limit', type=int, default=5,
                      help='Maximum number of invalid files to examine')
    parser.add_argument('--verbose', '-v', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    
    setup_logging(args.verbose)
    logger = logging.getLogger('ecod.examine')
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Find invalid XML files
    logger.info(f"Examining domain summary files for batch {args.batch_id}")
    invalid_files = find_invalid_xml_files(context, args.batch_id, args.limit)
    
    if not invalid_files:
        print(f"No invalid XML files found in batch {args.batch_id}")
        return 0
    
    # Display results
    print(f"\nFound {len(invalid_files)} invalid domain summary XML files:")
    for i, file_info in enumerate(invalid_files):
        print(f"\nFile {i+1}: {file_info['pdb_id']}_{file_info['chain_id']}")
        print(f"Path: {file_info['file_path']}")
        print(f"Size: {file_info['size']} bytes")
        print(f"DB file exists: {file_info['db_file_exists']}")
        print(f"DB file size: {file_info['db_file_size']}")
        
        if file_info['error']:
            print(f"Error: {file_info['error']}")
        
        print("\nFile Content:")
        print("-" * 70)
        print(file_info['content'])
        print("-" * 70)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())