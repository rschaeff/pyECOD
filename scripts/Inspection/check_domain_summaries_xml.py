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
        AND pf.file_exists = TRUE
    """
    
    summary_result = context.db.execute_query(summary_query, (batch_id,))
    summary_count = summary_result[0][0] if summary_result else 0
    
    batch_info["summary_count"] = summary_count
    completion_pct = (summary_count / batch_info["total_items"]) * 100 if batch_info["total_items"] > 0 else 0
    batch_info["completion_pct"] = completion_pct
    
    # Get sample files
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
    
    sample_files = context.db.execute_query(sample_query, (batch_id, sample_size))
    
    # Analyze sample files
    sample_results = []
    valid_count = 0
    blast_hit_counts = []
    domain_counts = []
    
    for sample in sample_files:
        pdb_id = sample[0]
        chain_id = sample[1]
        file_path = sample[2]
        
        full_path = os.path.join(batch_info["base_path"], file_path)
        
        # Check XML file
        file_result = check_xml_file(full_path)
        file_result["pdb_id"] = pdb_id
        file_result["chain_id"] = chain_id
        
        sample_results.append(file_result)
        
        if file_result["valid_xml"]:
            valid_count += 1
            blast_hit_counts.append(file_result["blast_hits"])
            domain_counts.append(file_result["domains"])
    
    # Calculate statistics
    batch_info["sample_size"] = len(sample_results)
    batch_info["valid_xml_count"] = valid_count
    batch_info["valid_xml_pct"] = (valid_count / len(sample_results)) * 100 if sample_results else 0
    
    if blast_hit_counts:
        batch_info["avg_blast_hits"] = sum(blast_hit_counts) / len(blast_hit_counts)
        batch_info["min_blast_hits"] = min(blast_hit_counts)
        batch_info["max_blast_hits"] = max(blast_hit_counts)
    
    if domain_counts:
        batch_info["avg_domains"] = sum(domain_counts) / len(domain_counts)
        batch_info["min_domains"] = min(domain_counts)
        batch_info["max_domains"] = max(domain_counts)
    
    batch_info["samples"] = sample_results
    
    return batch_info

def main():
    """Main function to check domain summary XML files"""
    parser = argparse.ArgumentParser(description='Check domain summary XML files')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int,
                      help='Check specific batch ID')
    parser.add_argument('--sample-size', type=int, default=5,
                      help='Number of files to sample per batch')
    parser.add_argument('--output', type=str,
                      help='Output file for detailed results')
    parser.add_argument('--all-batches', action='store_true',
                      help='Check all batches')
    parser.add_argument('--verbose', '-v', action='store_true',
                      help='Enable verbose output')
    parser.add_argument('--examine-file', type=str,
                      help='Path to a specific file to examine')
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose)
    logger = logging.getLogger("ecod.xml_check")
    
    # If examining a specific file
    if args.examine_file:
        result = check_xml_file(args.examine_file)
        print(f"\nFile: {os.path.basename(args.examine_file)}")
        print(f"Valid XML: {result['valid_xml']}")
        print(f"BLAST hits: {result['blast_hits']}")
        print(f"Domains: {result['domains']}")
        if result["error"]:
            print(f"Error: {result['error']}")
        return 0
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get batches to check
    if args.batch_id:
        batch_ids = [args.batch_id]
    elif args.all_batches:
        # Get all batches
        query = """
        SELECT id FROM ecod_schema.batch ORDER BY id
        """
        batch_results = context.db.execute_query(query)
        batch_ids = [row[0] for row in batch_results]
    else:
        # Get completed and indexed batches
        query = """
        SELECT id FROM ecod_schema.batch 
        WHERE status IN ('completed', 'indexed') 
        ORDER BY id
        """
        batch_results = context.db.execute_query(query)
        batch_ids = [row[0] for row in batch_results]
    
    # Check batches
    batch_results = []
    
    for batch_id in batch_ids:
        logger.info(f"Analyzing batch {batch_id}")
        batch_info = analyze_batch(batch_id, context, args.sample_size)
        
        if batch_info:
            batch_results.append(batch_info)
            
            # Print batch summary
            status_symbol = "✅" if batch_info["status"] == "completed" else "⏳" if batch_info["status"] == "indexed" else "🆕"
            summary_symbol = "✅" if batch_info["summary_count"] >= batch_info["total_items"] else "⚠️" if batch_info["summary_count"] > 0 else "❌"
            
            logger.info(f"{status_symbol} Batch {batch_info['id']} ({batch_info['name']}): {batch_info['status']}")
            logger.info(f"  {summary_symbol} Domain summaries: {batch_info['summary_count']}/{batch_info['total_items']} ({batch_info['completion_pct']:.1f}%)")
            
            # Print sample results
            if "valid_xml_count" in batch_info:
                logger.info(f"  Sample analysis: {batch_info['valid_xml_count']}/{batch_info['sample_size']} valid XML files ({batch_info['valid_xml_pct']:.1f}%)")
                
                if "avg_blast_hits" in batch_info:
                    logger.info(f"  BLAST hits: avg={batch_info['avg_blast_hits']:.1f}, min={batch_info['min_blast_hits']}, max={batch_info['max_blast_hits']}")
                
                if "avg_domains" in batch_info:
                    logger.info(f"  Domains: avg={batch_info['avg_domains']:.1f}, min={batch_info['min_domains']}, max={batch_info['max_domains']}")
            
            # Print errors if any
            errors = [s for s in batch_info["samples"] if not s["valid_xml"]]
            if errors:
                logger.warning(f"  Invalid XML files found in sample: {len(errors)}/{batch_info['sample_size']}")
                for e in errors[:3]:  # Show first 3 errors
                    logger.warning(f"    {e['pdb_id']}_{e['chain_id']}: {e['error']}")
                
                if len(errors) > 3:
                    logger.warning(f"    ... and {len(errors) - 3} more")
    
    # Calculate overall statistics
    if batch_results:
        total_batches = len(batch_results)
        completed_batches = sum(1 for b in batch_results if b["status"] == "completed")
        indexed_batches = sum(1 for b in batch_results if b["status"] == "indexed")
        created_batches = sum(1 for b in batch_results if b["status"] == "created")
        
        total_items = sum(b["total_items"] for b in batch_results)
        total_summaries = sum(b["summary_count"] for b in batch_results)
        
        # Valid XML statistics
        total_samples = sum(b["sample_size"] for b in batch_results)
        valid_samples = sum(b.get("valid_xml_count", 0) for b in batch_results)
        
        logger.info("\nOverall Summary:")
        logger.info(f"Total batches: {total_batches}")
        logger.info(f"  Completed: {completed_batches}")
        logger.info(f"  Indexed: {indexed_batches}")
        logger.info(f"  Created: {created_batches}")
        
        if total_items > 0:
            overall_pct = (total_summaries / total_items) * 100
            logger.info(f"Total domain summaries: {total_summaries}/{total_items} ({overall_pct:.1f}%)")
        
        if total_samples > 0:
            valid_pct = (valid_samples / total_samples) * 100
            logger.info(f"Valid XML files in samples: {valid_samples}/{total_samples} ({valid_pct:.1f}%)")
    
    # Write output file if requested
    if args.output and batch_results:
        import json
        
        # Convert to serializable format
        for batch in batch_results:
            # Convert sample objects to dicts
            samples = []
            for sample in batch.get("samples", []):
                sample_dict = {}
                for k, v in sample.items():
                    sample_dict[k] = v
                samples.append(sample_dict)
            batch["samples"] = samples
        
        with open(args.output, 'w') as f:
            json.dump(batch_results, f, indent=2)
        
        logger.info(f"Results written to {args.output}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())