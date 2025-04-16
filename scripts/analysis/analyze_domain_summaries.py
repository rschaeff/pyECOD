#!/usr/bin/env python3
"""
analyze_domain_summaries.py - Analyze domain summary content across batches
"""

import os
import sys
import json
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

def analyze_xml_file(file_path: str, detailed: bool = False) -> Dict[str, Any]:
    """Analyze a domain summary XML file"""
    result = {
        "file_path": file_path,
        "valid_xml": False,
        "has_blast_hits": False,
        "has_domains": False,
        "blast_hit_count": 0,
        "domain_count": 0,
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
                
                # Get chain BLAST run info
                chain_blast = blast_summ.find('chain_blast_run')
                if chain_blast is not None:
                    result["blast_program"] = chain_blast.get('program')
                    result["blast_version"] = chain_blast.get('version')
                    
                    # Get blast DB
                    blast_db = chain_blast.find('blast_db')
                    if blast_db is not None and blast_db.text:
                        result["blast_db"] = blast_db.text
                    
                    # Get query length
                    query_len = chain_blast.find('query_len')
                    if query_len is not None and query_len.text:
                        result["query_length"] = int(query_len.text)
        
        # Analyze BLAST hits
        hits_elem = root.find('.//hits')
        if hits_elem is not None:
            hit_elems = hits_elem.findall('hit')
            result["blast_hit_count"] = len(hit_elems)
            result["has_blast_hits"] = len(hit_elems) > 0
            
            # Analyze top hits if detailed
            if detailed and hit_elems:
                top_hits = []
                for i, hit in enumerate(hit_elems[:5]):  # Get top 5 hits
                    hit_info = {
                        "pdb_id": hit.get('pdb_id'),
                        "chain_id": hit.get('chain_id'),
                        "evalue": hit.get('evalues'),
                        "hsp_count": hit.get('hsp_count')
                    }
                    
                    # Get regions
                    query_reg = hit.find('query_reg')
                    hit_reg = hit.find('hit_reg')
                    
                    if query_reg is not None and query_reg.text:
                        hit_info["query_region"] = query_reg.text
                    
                    if hit_reg is not None and hit_reg.text:
                        hit_info["hit_region"] = hit_reg.text
                    
                    top_hits.append(hit_info)
                
                result["top_hits"] = top_hits
        
        # Analyze domains if present
        domains_elem = root.find('.//domains')
        if domains_elem is not None:
            domain_elems = domains_elem.findall('domain')
            result["domain_count"] = len(domain_elems)
            result["has_domains"] = len(domain_elems) > 0
            
            # Analyze domains if detailed
            if detailed and domain_elems:
                domains = []
                for domain in domain_elems:
                    domain_info = {
                        "id": domain.get('id'),
                        "range": domain.get('range')
                    }
                    
                    # Get classification if available
                    classification = domain.find('classification')
                    if classification is not None:
                        domain_info["t_group"] = classification.get('t_group')
                        domain_info["h_group"] = classification.get('h_group')
                        domain_info["x_group"] = classification.get('x_group')
                        domain_info["a_group"] = classification.get('a_group')
                    
                    domains.append(domain_info)
                
                result["domains"] = domains
                
                # Count classifications
                t_groups = Counter()
                h_groups = Counter()
                
                for domain in domains:
                    if "t_group" in domain and domain["t_group"]:
                        t_groups[domain["t_group"]] += 1
                    if "h_group" in domain and domain["h_group"]:
                        h_groups[domain["h_group"]] += 1
                
                result["t_group_counts"] = dict(t_groups)
                result["h_group_counts"] = dict(h_groups)
        
        return result
    
    except ET.ParseError as e:
        result["error"] = f"XML parsing error: {str(e)}"
        return result
    except Exception as e:
        result["error"] = f"Error: {str(e)}"
        return result

def analyze_batch_files(batch_id: int, context: ApplicationContext, 
                       sample_size: int = 10, detailed: bool = False) -> Dict[str, Any]:
    """Analyze domain summary files in a batch"""
    logger = logging.getLogger("ecod.analyze")
    
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
    
    # Get summary count
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
    batch_info["summary_count"] = summary_result[0][0] if summary_result else 0
    
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
    
    if not sample_files:
        logger.warning(f"No domain summary files found for batch {batch_id}")
        return batch_info
    
    # Analyze sample files
    valid_xml_count = 0
    with_blast_hits = 0
    with_domains = 0
    blast_hit_counts = []
    domain_counts = []
    samples = []
    
    for sample in sample_files:
        pdb_id = sample[0]
        chain_id = sample[1]
        file_path = sample[2]
        
        full_path = os.path.join(batch_info["base_path"], file_path)
        
        # Analyze file
        file_result = analyze_xml_file(full_path, detailed)
        file_result["pdb_id"] = pdb_id
        file_result["chain_id"] = chain_id
        
        samples.append(file_result)
        
        if file_result["valid_xml"]:
            valid_xml_count += 1
            
            if file_result["has_blast_hits"]:
                with_blast_hits += 1
                blast_hit_counts.append(file_result["blast_hit_count"])
            
            if file_result["has_domains"]:
                with_domains += 1
                domain_counts.append(file_result["domain_count"])
    
    # Calculate statistics
    batch_info["sample_size"] = len(samples)
    batch_info["valid_xml_count"] = valid_xml_count
    batch_info["with_blast_hits"] = with_blast_hits
    batch_info["with_domains"] = with_domains
    
    if valid_xml_count > 0:
        batch_info["valid_xml_pct"] = (valid_xml_count / len(samples)) * 100
        batch_info["with_blast_hits_pct"] = (with_blast_hits / valid_xml_count) * 100
        batch_info["with_domains_pct"] = (with_domains / valid_xml_count) * 100
    
    if blast_hit_counts:
        batch_info["avg_blast_hits"] = sum(blast_hit_counts) / len(blast_hit_counts)
        batch_info["min_blast_hits"] = min(blast_hit_counts)
        batch_info["max_blast_hits"] = max(blast_hit_counts)
    
    if domain_counts:
        batch_info["avg_domains"] = sum(domain_counts) / len(domain_counts)
        batch_info["min_domains"] = min(domain_counts)
        batch_info["max_domains"] = max(domain_counts)
    
    # Add detailed sample results if requested
    if detailed:
        batch_info["samples"] = samples
    
    return batch_info

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Analyze domain summary content')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int,
                      help='Analyze specific batch ID')
    parser.add_argument('--sample-size', type=int, default=10,
                      help='Number of files to sample per batch')
    parser.add_argument('--detailed', action='store_true',
                      help='Include detailed analysis of samples')
    parser.add_argument('--output', type=str,
                      help='Output file for results')
    parser.add_argument('--file', type=str,
                      help='Analyze a specific file')
    parser.add_argument('--verbose', '-v', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    
    setup_logging(args.verbose)
    logger = logging.getLogger("ecod.analyze")
    
    # Analyze a specific file if requested
    if args.file:
        result = analyze_xml_file(args.file, detailed=True)
        print(f"\nFile: {os.path.basename(args.file)}")
        print(f"Valid XML: {result['valid_xml']}")
        
        if result["valid_xml"]:
            print(f"PDB ID: {result.get('pdb_id', 'Unknown')}")
            print(f"Chain ID: {result.get('chain_id', 'Unknown')}")
            print(f"Query length: {result.get('query_length', 'Unknown')}")
            print(f"BLAST hits: {result['blast_hit_count']}")
            print(f"Domains: {result['domain_count']}")
            
            if result.get('top_hits'):
                print("\nTop BLAST hits:")
                for i, hit in enumerate(result['top_hits']):
                    print(f"  Hit {i+1}: {hit['pdb_id']}_{hit['chain_id']}")
                    print(f"    E-value: {hit['evalue']}")
                    print(f"    Query region: {hit.get('query_region', 'Unknown')}")
                    print(f"    Hit region: {hit.get('hit_region', 'Unknown')}")
            
            if result.get('domains'):
                print("\nDomains:")
                for i, domain in enumerate(result['domains']):
                    print(f"  Domain {i+1}: {domain.get('id', 'Unknown')}")
                    print(f"    Range: {domain.get('range', 'Unknown')}")
                    if 't_group' in domain:
                        print(f"    Classification: {domain.get('t_group', '')}, {domain.get('h_group', '')}")
        
        if result["error"]:
            print(f"Error: {result['error']}")
        
        return 0
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get batches to analyze
    if args.batch_id:
        batch_ids = [args.batch_id]
    else:
        # Get completed batches
        query = """
        SELECT id FROM ecod_schema.batch 
        WHERE status = 'completed' 
        ORDER BY id
        """
        batch_results = context.db.execute_query(query)
        batch_ids = [row[0] for row in batch_results]
    
    if not batch_ids:
        logger.error("No batches found to analyze")
        return 1
    
    # Analyze batches
    batch_results = []
    
    for batch_id in batch_ids:
        logger.info(f"Analyzing batch {batch_id}")
        result = analyze_batch_files(
            batch_id, 
            context, 
            sample_size=args.sample_size,
            detailed=args.detailed
        )
        
        if result:
            batch_results.append(result)
            
            # Display summary
            logger.info(f"Batch {batch_id} ({result['name']}):")
            logger.info(f"  Summary count: {result['summary_count']}/{result['total_items']} ({result['summary_count']/result['total_items']*100:.1f}%)")
            
            if "valid_xml_count" in result:
                logger.info(f"  Valid XML: {result['valid_xml_count']}/{result['sample_size']} ({result['valid_xml_pct']:.1f}%)")
            
            if "with_blast_hits" in result:
                logger.info(f"  With BLAST hits: {result['with_blast_hits']}/{result['valid_xml_count']} ({result['with_blast_hits_pct']:.1f}%)")
                if "avg_blast_hits" in result:
                    logger.info(f"  BLAST hits: avg={result['avg_blast_hits']:.1f}, min={result['min_blast_hits']}, max={result['max_blast_hits']}")
            
            if "with_domains" in result:
                logger.info(f"  With domains: {result['with_domains']}/{result['valid_xml_count']} ({result['with_domains_pct']:.1f}%)")
                if "avg_domains" in result:
                    logger.info(f"  Domains: avg={result['avg_domains']:.1f}, min={result['min_domains']}, max={result['max_domains']}")
    
    # Calculate overall statistics
    if batch_results:
        logger.info("\nOverall Statistics:")
        logger.info(f"Batches analyzed: {len(batch_results)}")
        
        # Calculate totals
        total_items = sum(batch["total_items"] for batch in batch_results)
        total_summaries = sum(batch["summary_count"] for batch in batch_results)
        
        # Sample stats
        total_samples = sum(batch["sample_size"] for batch in batch_results)
        valid_samples = sum(batch["valid_xml_count"] for batch in batch_results)
        samples_with_blast = sum(batch.get("with_blast_hits", 0) for batch in batch_results)
        samples_with_domains = sum(batch.get("with_domains", 0) for batch in batch_results)
        
        logger.info(f"Total summaries: {total_summaries}/{total_items} ({total_summaries/total_items*100:.1f}%)")
        logger.info(f"Valid XML files: {valid_samples}/{total_samples} ({valid_samples/total_samples*100:.1f}%)")
        
        if valid_samples > 0:
            logger.info(f"Files with BLAST hits: {samples_with_blast}/{valid_samples} ({samples_with_blast/valid_samples*100:.1f}%)")
            logger.info(f"Files with domains: {samples_with_domains}/{valid_samples} ({samples_with_domains/valid_samples*100:.1f}%)")
    
    # Write output file if requested
    if args.output and batch_results:
        with open(args.output, 'w') as f:
            json.dump(batch_results, f, indent=2)
        
        logger.info(f"Results written to {args.output}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())