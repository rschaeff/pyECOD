#!/usr/bin/env python3
"""
check_blast_hits.py - Check if BLAST files contain hits
"""

import os
import sys
import logging
import argparse
import xml.etree.ElementTree as ET
from typing import Dict, Any, Optional

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

def check_blast_file(file_path: str) -> Dict[str, Any]:
    """Check a BLAST XML file for hits"""
    results = {
        "file_path": file_path,
        "exists": os.path.exists(file_path),
        "hit_count": 0,
        "hits": []
    }
    
    if not results["exists"]:
        return results
    
    try:
        tree = ET.parse(file_path)
        root = tree.getroot()
        
        # Check for iterations (there's usually one per query)
        for iteration in root.findall(".//Iteration"):
            # Look for hits in this iteration
            hits = iteration.findall(".//Hit")
            
            results["hit_count"] += len(hits)
            
            # Get details for each hit
            for hit in hits:
                hit_info = {
                    "num": hit.findtext("Hit_num", ""),
                    "id": hit.findtext("Hit_id", ""),
                    "def": hit.findtext("Hit_def", ""),
                    "hsps": []
                }
                
                # Get HSPs for this hit
                for hsp in hit.findall(".//Hsp"):
                    hsp_info = {
                        "evalue": hsp.findtext("Hsp_evalue", ""),
                        "query_from": hsp.findtext("Hsp_query-from", ""),
                        "query_to": hsp.findtext("Hsp_query-to", ""),
                        "hit_from": hsp.findtext("Hsp_hit-from", ""),
                        "hit_to": hsp.findtext("Hsp_hit-to", "")
                    }
                    hit_info["hsps"].append(hsp_info)
                
                results["hits"].append(hit_info)
        
    except Exception as e:
        results["error"] = str(e)
    
    return results

def main():
    """Check BLAST files for hits"""
    parser = argparse.ArgumentParser(description='Check BLAST files for hits')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID')
    parser.add_argument('--protein-id', type=int, required=True,
                      help='Protein ID')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger("ecod.check_blast_hits")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get protein and batch information
    query = """
    SELECT 
        p.id, p.source_id, p.pdb_id, p.chain_id,
        ps.id as process_id
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    WHERE 
        p.id = %s AND ps.batch_id = %s
    """
    
    result = context.db.execute_query(query, (args.protein_id, args.batch_id))
    
    if not result:
        logger.error(f"Protein {args.protein_id} not found in batch {args.batch_id}")
        return 1
    
    protein_info = result[0]
    pdb_id = protein_info[2]
    chain_id = protein_info[3]
    process_id = protein_info[4]
    pdb_chain = f"{pdb_id}_{chain_id}"
    
    logger.info(f"Checking BLAST hits for protein: {pdb_chain}")
    
    # Get file paths
    file_query = """
    SELECT 
        file_type, file_path
    FROM 
        ecod_schema.process_file
    WHERE 
        process_id = %s AND file_type IN ('chain_blast_result', 'domain_blast_result')
    """
    
    files = context.db.execute_query(file_query, (process_id,))
    
    if not files:
        logger.error(f"No BLAST files found for protein {pdb_chain}")
        return 1
    
    for file_info in files:
        file_type = file_info[0]
        rel_path = file_info[1]
        
        # Get batch path
        batch_query = """
        SELECT base_path
        FROM ecod_schema.batch b
        JOIN ecod_schema.process_status ps ON b.id = ps.batch_id
        WHERE ps.id = %s
        """
        
        batch_result = context.db.execute_query(batch_query, (process_id,))
        
        if not batch_result:
            logger.error(f"Could not find batch path for process {process_id}")
            continue
        
        batch_path = batch_result[0][0]
        
        # Construct full path
        if os.path.isabs(rel_path):
            full_path = rel_path
        else:
            full_path = os.path.join(batch_path, rel_path)
        
        logger.info(f"Checking {file_type}: {full_path}")
        
        # Check the BLAST file
        results = check_blast_file(full_path)
        
        if "error" in results:
            logger.error(f"Error checking {file_type}: {results['error']}")
            continue
        
        logger.info(f"File exists: {results['exists']}")
        logger.info(f"Hit count: {results['hit_count']}")
        
        if results['hit_count'] == 0:
            logger.warning(f"No hits found in {file_type}")
        else:
            # Log hit details
            logger.info(f"BLAST hits in {file_type}:")
            for i, hit in enumerate(results['hits'][:5]):  # Show first 5 hits
                logger.info(f"  Hit {i+1}: {hit['def']}")
                for j, hsp in enumerate(hit['hsps'][:3]):  # Show first 3 HSPs per hit
                    logger.info(f"    HSP {j+1}: E-value={hsp['evalue']}, Query={hsp['query_from']}-{hsp['query_to']}, Hit={hsp['hit_from']}-{hsp['hit_to']}")
            
            if len(results['hits']) > 5:
                logger.info(f"  ... and {len(results['hits']) - 5} more hits")
    
    # Also check the domain summary
    summary_query = """
    SELECT file_path
    FROM ecod_schema.process_file
    WHERE process_id = %s AND file_type = 'domain_summary'
    """
    
    summary_result = context.db.execute_query(summary_query, (process_id,))
    
    if summary_result:
        summary_path = summary_result[0][0]
        
        if os.path.isabs(summary_path):
            full_summary_path = summary_path
        else:
            full_summary_path = os.path.join(batch_path, summary_path)
        
        logger.info(f"Checking domain summary: {full_summary_path}")
        
        try:
            tree = ET.parse(full_summary_path)
            root = tree.getroot()
            
            # Check for chain blast hits
            chain_hits = root.findall(".//chain_blast_run/hits/hit")
            logger.info(f"Chain BLAST hits in summary: {len(chain_hits)}")
            
            # Check for domain blast hits
            domain_hits = root.findall(".//blast_run/hits/hit")
            logger.info(f"Domain BLAST hits in summary: {len(domain_hits)}")
            
            if len(chain_hits) == 0 and len(domain_hits) == 0:
                logger.warning("Domain summary contains no hits")
                
                # Check if hit filtering is too strict
                logger.info("Looking at DomainSummary processing parameters...")
                logger.info("Possible issues:")
                logger.info("1. hsp_evalue_threshold might be too stringent (default: 0.005)")
                logger.info("2. hit_coverage_threshold might be too high (default: 0.7)")
                logger.info("3. Query length is very small (25 residues), which might lead to no significant hits")
                
        except Exception as e:
            logger.error(f"Error checking domain summary: {str(e)}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())