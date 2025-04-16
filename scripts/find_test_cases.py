#!/usr/bin/env python3
"""
find_test_cases.py - Find test proteins for the blast-only domain partition pipeline

This script identifies proteins that have BLAST files but no domain summary or partition files,
suitable for testing the domain summary generation process.
"""

import os
import sys
import logging
import argparse
import xml.etree.ElementTree as ET
from typing import Dict, List, Any, Optional
from pathlib import Path

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
    
    return logging.getLogger("ecod.find_test_cases")

def check_xml_validity(file_path: str) -> Dict[str, Any]:
    """
    Check if an XML file is valid and extract basic information
    
    Args:
        file_path: Path to XML file
        
    Returns:
        Dictionary with validation results
    """
    result = {
        "file_path": file_path,
        "exists": os.path.exists(file_path),
        "valid_xml": False,
        "file_size": 0,
        "error": None
    }
    
    if not result["exists"]:
        result["error"] = "File does not exist"
        return result
    
    try:
        result["file_size"] = os.path.getsize(file_path)
        
        # Check if file is empty
        if result["file_size"] == 0:
            result["error"] = "File is empty"
            return result
        
        # Try to parse XML
        tree = ET.parse(file_path)
        root = tree.getroot()
        
        # Mark as valid XML
        result["valid_xml"] = True
        result["root_tag"] = root.tag
        
        # Extract information if it's a BLAST result
        if root.tag == 'BlastOutput':
            hits = root.findall(".//Hit")
            result["hit_count"] = len(hits)
        
        # Extract information if it's a domain summary
        if root.tag == 'blast_summ_doc':
            blast_summ = root.find('blast_summ')
            if blast_summ is not None:
                result["pdb_id"] = blast_summ.get('pdb')
                result["chain_id"] = blast_summ.get('chain')
            
            hits_elem = root.find('.//hits')
            if hits_elem is not None:
                hit_elems = hits_elem.findall('hit')
                result["hit_count"] = len(hit_elems)
        
        return result
    
    except ET.ParseError as e:
        result["error"] = f"XML parsing error: {str(e)}"
        return result
    except Exception as e:
        result["error"] = f"Error checking file: {str(e)}"
        return result

def find_test_cases(context, batch_id: int, limit: int = 5) -> List[Dict[str, Any]]:
    """
    Find test cases that have BLAST files but no domain summary files
    
    Args:
        context: Application context
        batch_id: Batch ID to search in
        limit: Maximum number of test cases to return
        
    Returns:
        List of protein dictionaries with metadata
    """
    logger = logging.getLogger("ecod.find_test_cases")
    
    # First get batch information
    batch_query = """
    SELECT 
        id, batch_name, base_path, ref_version
    FROM 
        ecod_schema.batch
    WHERE 
        id = %s
    """
    
    batch_result = context.db.execute_query(batch_query, (batch_id,))
    
    if not batch_result:
        logger.error(f"Batch {batch_id} not found")
        return []
    
    batch_info = {
        'id': batch_result[0][0],
        'name': batch_result[0][1],
        'base_path': batch_result[0][2],
        'reference': batch_result[0][3]
    }
    
    logger.info(f"Analyzing batch {batch_id} ({batch_info['name']})")
    
    # Query for proteins that have BLAST files but no domain summary
    query = """
    SELECT 
        p.id, p.pdb_id, p.chain_id, p.length,
        ps.id as process_id, ps.current_stage, ps.status
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    JOIN 
        ecod_schema.batch b ON ps.batch_id = b.id
    JOIN 
        ecod_schema.process_file pf_blast ON ps.id = pf_blast.process_id
    LEFT JOIN 
        ecod_schema.process_file pf_summary ON ps.id = pf_summary.process_id 
        AND pf_summary.file_type = 'domain_summary'
    WHERE 
        b.id = %s
        AND pf_blast.file_type = 'blast_result'
        AND pf_blast.file_exists = TRUE
        AND (pf_summary.id IS NULL OR pf_summary.file_exists = FALSE)
    GROUP BY
        p.id, p.pdb_id, p.chain_id, p.length, ps.id, ps.current_stage, ps.status
    ORDER BY 
        p.id
    LIMIT 100
    """
    
    candidates = context.db.execute_query(query, (batch_id,))
    
    if not candidates:
        logger.info(f"No proteins found with BLAST files but no domain summary in batch {batch_id}")
        return []
    
    logger.info(f"Found {len(candidates)} candidate proteins")
    
    # Process candidates
    results = []
    processed = 0
    
    for candidate in candidates:
        protein_id = candidate[0]
        pdb_id = candidate[1]
        chain_id = candidate[2]
        length = candidate[3]
        process_id = candidate[4]
        
        # Skip processing if we've found enough test cases
        if len(results) >= limit:
            break
        
        # Get file paths
        file_query = """
        SELECT 
            file_type, file_path, file_exists
        FROM 
            ecod_schema.process_file
        WHERE 
            process_id = %s
        """
        
        file_results = context.db.execute_query(file_query, (process_id,))
        
        file_paths = {}
        for file_row in file_results:
            file_type = file_row[0]
            file_path = file_row[1]
            file_exists = file_row[2]
            
            file_paths[file_type] = {
                'path': file_path,
                'exists': file_exists
            }
        
        # Check if chain and domain BLAST files exist
        chain_blast_path = None
        if 'blast_result' in file_paths:
            rel_path = file_paths['blast_result']['path']
            full_path = os.path.join(batch_info['base_path'], rel_path)
            full_path = os.path.normpath(full_path)
            
            if os.path.exists(full_path):
                chain_blast_path = full_path
        
        domain_blast_path = None
        if 'domain_blast_result' in file_paths:
            rel_path = file_paths['domain_blast_result']['path']
            full_path = os.path.join(batch_info['base_path'], rel_path)
            full_path = os.path.normpath(full_path)
            
            if os.path.exists(full_path):
                domain_blast_path = full_path
        
        # Skip if neither BLAST file exists
        if not chain_blast_path and not domain_blast_path:
            logger.debug(f"Skipping {pdb_id}_{chain_id} - no BLAST files found")
            continue
        
        # Verify XML validity
        chain_blast_valid = False
        domain_blast_valid = False
        
        if chain_blast_path:
            chain_check = check_xml_validity(chain_blast_path)
            chain_blast_valid = chain_check.get("valid_xml", False)
        
        if domain_blast_path:
            domain_check = check_xml_validity(domain_blast_path)
            domain_blast_valid = domain_check.get("valid_xml", False)
        
        # Skip if neither BLAST file is valid
        if not chain_blast_valid and not domain_blast_valid:
            logger.debug(f"Skipping {pdb_id}_{chain_id} - no valid BLAST files")
            continue
        
        # Add to results
        result = {
            'protein_id': protein_id,
            'pdb_id': pdb_id,
            'chain_id': chain_id,
            'length': length,
            'process_id': process_id,
            'chain_blast_path': chain_blast_path,
            'domain_blast_path': domain_blast_path,
            'chain_blast_valid': chain_blast_valid,
            'domain_blast_valid': domain_blast_valid
        }
        
        results.append(result)
        processed += 1
        
        # Log progress
        if processed % 10 == 0:
            logger.info(f"Processed {processed} candidates, found {len(results)} test cases")
    
    # Sort results to get diverse test cases
    # Prioritize proteins with both valid BLAST files
    results.sort(key=lambda x: (not (x['chain_blast_valid'] and x['domain_blast_valid']), not x['domain_blast_valid']))
    
    # Take only the requested number of test cases
    final_results = results[:limit]
    
    return final_results

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Find test cases for domain summary generation')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, default=31,
                      help='Batch ID to analyze (default: 31)')
    parser.add_argument('--limit', type=int, default=5,
                      help='Number of test cases to find (default: 5)')
    parser.add_argument('--output', type=str,
                      help='Output file for test cases')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('--verbose', '-v', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    logger = setup_logging(args.verbose, args.log_file)
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Find test cases
    test_cases = find_test_cases(context, args.batch_id, args.limit)
    
    if not test_cases:
        logger.error("No suitable test cases found")
        return 1
    
    # Display test cases
    print("\nFound test cases for domain summary generation:\n")
    print(f"{'Protein ID':<10} {'PDB Chain':<10} {'Length':<8} {'Chain BLAST':<12} {'Domain BLAST':<12}")
    print("-" * 60)
    
    for case in test_cases:
        chain_blast = "Valid" if case['chain_blast_valid'] else "Invalid" if case['chain_blast_path'] else "Missing"
        domain_blast = "Valid" if case['domain_blast_valid'] else "Invalid" if case['domain_blast_path'] else "Missing"
        
        print(f"{case['protein_id']:<10} {case['pdb_id']}_{case['chain_id']:<5} {case['length']:<8} {chain_blast:<12} {domain_blast:<12}")
    
    # Generate test commands
    print("\nTest commands for domain summary generation:\n")
    
    for case in test_cases:
        cmd = f"python scripts/generate_domain_summary_v2.py --config config/config.yml --batch-id {args.batch_id} --protein-id {case['protein_id']} --blast-only -v"
        print(cmd)
    
    # Write output file if requested
    if args.output:
        import json
        with open(args.output, 'w') as f:
            json.dump(test_cases, f, indent=2)
        logger.info(f"Test cases written to {args.output}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())