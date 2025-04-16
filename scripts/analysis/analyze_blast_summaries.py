#!/usr/bin/env python3
"""
analyze_blast_summaries.py - Analyze the XML structure of blast_summ files in a batch
"""

import os
import sys
import logging
import argparse
import xml.etree.ElementTree as ET
from collections import Counter
from typing import Dict, List, Tuple, Optional, Set

# Add parent directory to path
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

def validate_blast_summ_structure(file_path: str) -> Tuple[bool, str, Set[str]]:
    """
    Validate if a file has the correct blast_summ structure
    
    Returns:
        Tuple of (is_valid, structure_type, root_elements)
    """
    try:
        tree = ET.parse(file_path)
        root = tree.getroot()
        
        # Get root tag and all top-level elements
        root_tag = root.tag
        root_elements = {child.tag for child in root}
        
        # Check for expected structure
        if root_tag == "blast_summ_doc" and "blast_summ" in root_elements:
            # This is the correct structure
            return True, "blast_summ_doc", root_elements
        elif root_tag == "BlastOutput":
            # This is a raw BLAST output
            return False, "blast_output", root_elements
        else:
            # Some other structure
            return False, root_tag, root_elements
    
    except Exception as e:
        return False, f"error: {str(e)}", set()

def analyze_batch_summaries(batch_id: int, context: ApplicationContext) -> Dict:
    """Analyze blast summary files for a batch"""
    logger = logging.getLogger("ecod.analyze_summaries")
    
    # Get batch path
    batch_query = "SELECT base_path FROM ecod_schema.batch WHERE id = %s"
    batch_result = context.db.execute_query(batch_query, (batch_id,))
    
    if not batch_result:
        logger.error(f"Batch {batch_id} not found")
        return {}
        
    batch_path = batch_result[0][0]
    
    # Query for all domain summary files in the batch
    query = """
    SELECT 
        p.id, p.pdb_id, p.chain_id, pf.file_path
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    JOIN 
        ecod_schema.process_file pf ON ps.id = pf.process_id
    WHERE 
        ps.batch_id = %s 
        AND pf.file_type = 'domain_summary'
        AND pf.file_exists = TRUE
    ORDER BY p.pdb_id, p.chain_id
    """
    
    summaries = context.db.execute_query(query, (batch_id,))
    
    if not summaries:
        logger.error(f"No domain summaries found for batch {batch_id}")
        return {}
    
    logger.info(f"Found {len(summaries)} domain summaries to analyze")
    
    # Analyze each summary file
    results = {
        "total": len(summaries),
        "valid": 0,
        "invalid": 0,
        "errors": 0,
        "structure_types": Counter(),
        "root_elements": Counter(),
        "invalid_files": []
    }
    
    for i, (protein_id, pdb_id, chain_id, file_path) in enumerate(summaries):
        # Prepare absolute path if needed
        abs_path = os.path.join(batch_path, file_path) if not os.path.isabs(file_path) else file_path
        
        if not os.path.exists(abs_path):
            logger.warning(f"File does not exist: {abs_path}")
            results["errors"] += 1
            results["invalid_files"].append({
                "pdb_id": pdb_id,
                "chain_id": chain_id,
                "file_path": file_path,
                "issue": "File not found"
            })
            continue
        
        # Validate file structure
        is_valid, structure_type, root_elements = validate_blast_summ_structure(abs_path)
        
        # Update statistics
        results["structure_types"][structure_type] += 1
        
        for elem in root_elements:
            results["root_elements"][elem] += 1
        
        if is_valid:
            results["valid"] += 1
        else:
            results["invalid"] += 1
            results["invalid_files"].append({
                "pdb_id": pdb_id,
                "chain_id": chain_id,
                "file_path": file_path,
                "structure": structure_type
            })
        
        # Progress update
        if (i + 1) % 100 == 0 or (i + 1) == len(summaries):
            logger.info(f"Analyzed {i+1}/{len(summaries)} files. Valid: {results['valid']}, Invalid: {results['invalid']}")
    
    return results

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Analyze XML structure of blast_summ files')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to analyze')
    parser.add_argument('--output-file', type=str,
                      help='Write results to this file')
    parser.add_argument('--max-examples', type=int, default=10,
                      help='Maximum number of invalid file examples to list')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    
    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger("ecod.analyze_summaries")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Analyze batch summaries
    results = analyze_batch_summaries(args.batch_id, context)
    
    if not results:
        return 1
    
    # Display results
    logger.info("\n--- Analysis Results ---")
    logger.info(f"Total files analyzed: {results['total']}")
    logger.info(f"Valid structure (blast_summ_doc): {results['valid']} ({results['valid']/results['total']*100:.1f}%)")
    logger.info(f"Invalid structure: {results['invalid']} ({results['invalid']/results['total']*100:.1f}%)")
    logger.info(f"Files with errors: {results['errors']}")
    
    logger.info("\nStructure types found:")
    for struct_type, count in results['structure_types'].most_common():
        logger.info(f"  {struct_type}: {count} ({count/results['total']*100:.1f}%)")
    
    logger.info("\nRoot elements found:")
    for elem, count in results['root_elements'].most_common():
        logger.info(f"  {elem}: {count}")
    
    # Show examples of invalid files
    if results['invalid_files']:
        examples_count = min(len(results['invalid_files']), args.max_examples)
        logger.info(f"\nExamples of invalid files (showing {examples_count} of {len(results['invalid_files'])}):")
        
        for i, invalid_file in enumerate(results['invalid_files'][:examples_count]):
            logger.info(f"  {i+1}. {invalid_file['pdb_id']}_{invalid_file['chain_id']}: {invalid_file.get('structure', 'unknown')}")
            logger.info(f"     Path: {invalid_file['file_path']}")
    
    # Write results to file if requested
    if args.output_file:
        try:
            import json
            os.makedirs(os.path.dirname(args.output_file), exist_ok=True)
            
            # Convert Counter objects to dictionaries for JSON serialization
            results['structure_types'] = dict(results['structure_types'])
            results['root_elements'] = dict(results['root_elements'])
            
            with open(args.output_file, 'w') as f:
                json.dump(results, f, indent=2)
            
            logger.info(f"\nResults written to {args.output_file}")
        except Exception as e:
            logger.error(f"Error writing results file: {str(e)}")
    
    return 0 if results['valid'] == results['total'] else 1

if __name__ == "__main__":
    sys.exit(main())