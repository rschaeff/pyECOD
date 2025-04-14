#!/usr/bin/env python3
"""
check_domain_summaries_xml.py - Check domain summary XML files across batches
"""

import os
import sys
import logging
import argparse
import xml.etree.ElementTree as ET
from typing import Dict, List, Any, Optional
from collections import Counter

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

def check_xml_file(file_path: str) -> Dict[str, Any]:
    """Check a domain summary XML file"""
    result = {
        "file_path": file_path,
        "valid_xml": False,
        "blast_hits": 0,
        "domains": 0,
        "error": None
    }
    
    if not os.path.exists(file_path):
        result["error"] = "File not found"
        return result
    
    try:
        # Parse XML
        tree = ET.parse(file_path)
        root = tree.getroot()
        
        # Mark as valid XML
        result["valid_xml"] = True
        
        # Get basic information
        if root.tag == 'blast_summ_doc':
            blast_summ = root.find('blast_summ')
            if blast_summ is not None:
                result["pdb_id"] = blast_summ.get('pdb')
                result["chain_id"] = blast_summ.get('chain')
        
        # Count BLAST hits
        hits_elem = root.find('.//hits')
        if hits_elem is not None:
            hit_elems = hits_elem.findall('hit')
            result["blast_hits"] = len(hit_elems)
        
        # Count domains if present
        domains_elem = root.find('.//domains')
        if domains_elem is not None:
            domain_elems = domains_elem.findall('domain')
            result["domains"] = len(domain_elems)
        
        return result
    
    except ET.ParseError as e:
        result["error"] = f"XML parsing error: {str(e)}"
        return result
    except Exception as e:
        result["error"] = f"Error: {str(e)}"
        return result

def analyze_batch(batch_id: int, context: ApplicationContext, sample_size: int = 5) -> Dict[str, Any]:
    """Analyze XML files in a batch"""
    logger = logging.getLogger("ecod.xml_check")
    
    # Get batch info
    batch_query = """
    SELECT 
        b.id, b.batch_name, b.base_path, b.status, b.total_items
    FROM 
        ecod_schema.batch b
    WHERE 
        b.id = %s
    """
    
    batch_result = context.db.execute_query(batch_query, (batch_id,))
    if not batch_result:
        logger.error(f"Batch {batch_id} not found")
        return {}
    
    batch_info = {
        "id": batch_result[0][0],
        "name": batch_result[0][1],
        "base_path": batch_result[0][2],
        "status": batch_result[0][3],
        "total_items": batch_result[0][4]
    }
    
    # Get files count
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