#!/usr/bin/env python3
"""
analyze_batch_domain_summaries.py - Analyze success and failure patterns in domain summaries

This script analyzes a specific batch to identify patterns in domain summary generation,
including successful summaries, those with no hits, and potential failure modes.
"""

import os
import sys
import argparse
import logging
import xml.etree.ElementTree as ET
from typing import Dict, List, Any, Optional
from collections import Counter
import json

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

# Import application context (assuming standard pyECOD structure)
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

def get_batch_info(context, batch_id: int) -> Dict[str, Any]:
    """Get basic information about the batch"""
    query = """
    SELECT 
        b.id, b.batch_name, b.base_path, b.type, b.ref_version, 
        b.total_items, b.completed_items, b.status
    FROM 
        ecod_schema.batch b
    WHERE 
        b.id = %s
    """
    
    result = context.db.execute_dict_query(query, (batch_id,))
    if not result:
        return {}
    
    return result[0]

def get_process_status_summary(context, batch_id: int) -> Dict[str, int]:
    """Get summary of process status for the batch"""
    query = """
    SELECT 
        current_stage, status, COUNT(*) as count
    FROM 
        ecod_schema.process_status
    WHERE 
        batch_id = %s
    GROUP BY 
        current_stage, status
    ORDER BY 
        current_stage, status
    """
    
    results = context.db.execute_dict_query(query, (batch_id,))
    
    # Convert to a more structured format
    summary = {}
    for row in results:
        stage = row['current_stage']
        status = row['status']
        count = row['count']
        
        if stage not in summary:
            summary[stage] = {}
        
        summary[stage][status] = count
    
    return summary

def get_file_type_counts(context, batch_id: int) -> Dict[str, Dict[str, int]]:
    """Get counts of different file types for the batch"""
    query = """
    SELECT 
        file_type, 
        COUNT(*) as total_count,
        SUM(CASE WHEN file_exists = TRUE THEN 1 ELSE 0 END) as exists_count
    FROM 
        ecod_schema.process_file pf
    JOIN 
        ecod_schema.process_status ps ON pf.process_id = ps.id
    WHERE 
        ps.batch_id = %s
    GROUP BY 
        file_type
    ORDER BY 
        file_type
    """
    
    results = context.db.execute_dict_query(query, (batch_id,))
    
    counts = {}
    for row in results:
        counts[row['file_type']] = {
            'total': row['total_count'],
            'exists': row['exists_count']
        }
    
    return counts

def get_failure_patterns(context, batch_id: int) -> Dict[str, Any]:
    """Identify patterns in failures based on process_status"""
    # Get error messages with counts
    error_query = """
    SELECT 
        error_message, COUNT(*) as count
    FROM 
        ecod_schema.process_status
    WHERE 
        batch_id = %s
        AND status = 'error'
        AND error_message IS NOT NULL
    GROUP BY 
        error_message
    ORDER BY 
        count DESC
    """
    
    error_results = context.db.execute_dict_query(error_query, (batch_id,))
    
    # Get proteins with blast results but no domain summary
    missing_summary_query = """
    SELECT 
        COUNT(*) as count
    FROM 
        ecod_schema.process_status ps
    WHERE 
        ps.batch_id = %s
        AND EXISTS (
            SELECT 1 FROM ecod_schema.process_file pf1
            WHERE pf1.process_id = ps.id 
            AND pf1.file_type = 'chain_blast_result' 
            AND pf1.file_exists = TRUE
        )
        AND NOT EXISTS (
            SELECT 1 FROM ecod_schema.process_file pf2
            WHERE pf2.process_id = ps.id 
            AND pf2.file_type = 'domain_summary'
        )
    """
    
    missing_summary_result = context.db.execute_dict_query(missing_summary_query, (batch_id,))
    
    # Get sample of proteins with blast results but no domain summary
    sample_query = """
    SELECT 
        p.pdb_id, p.chain_id, ps.current_stage, ps.status, ps.error_message
    FROM 
        ecod_schema.process_status ps
    JOIN 
        ecod_schema.protein p ON ps.protein_id = p.id
    WHERE 
        ps.batch_id = %s
        AND EXISTS (
            SELECT 1 FROM ecod_schema.process_file pf1
            WHERE pf1.process_id = ps.id 
            AND pf1.file_type = 'chain_blast_result' 
            AND pf1.file_exists = TRUE
        )
        AND NOT EXISTS (
            SELECT 1 FROM ecod_schema.process_file pf2
            WHERE pf2.process_id = ps.id 
            AND pf2.file_type = 'domain_summary'
        )
    LIMIT 10
    """
    
    sample_results = context.db.execute_dict_query(sample_query, (batch_id,))
    
    return {
        'error_patterns': error_results,
        'missing_summary_count': missing_summary_result[0]['count'] if missing_summary_result else 0,
        'sample_missing_summaries': sample_results
    }

def analyze_domain_summary_xml(file_path: str) -> Dict[str, Any]:
    """Analyze a domain summary XML file"""
    result = {
        'valid_xml': False,
        'has_blast_hits': False,
        'blast_hit_count': 0,
        'has_domains': False,
        'domain_count': 0,
        'error': None
    }
    
    if not os.path.exists(file_path):
        result['error'] = "File not found"
        return result
    
    try:
        # Parse XML
        tree = ET.parse(file_path)
        root = tree.getroot()
        
        # Mark as valid XML
        result['valid_xml'] = True
        
        # Check for blast_summ element
        blast_summ = root.find('blast_summ')
        if blast_summ is not None:
            result['pdb_id'] = blast_summ.get('pdb')
            result['chain_id'] = blast_summ.get('chain')
        
        # Check for chain blast hits
        chain_blast = root.find('.//chain_blast_run')
        if chain_blast is not None:
            chain_hits = chain_blast.find('hits')
            if chain_hits is not None:
                chain_hit_elems = chain_hits.findall('hit')
                result['chain_blast_hit_count'] = len(chain_hit_elems)
                result['has_chain_blast_hits'] = len(chain_hit_elems) > 0
        
        # Check for domain blast hits
        domain_blast = root.find('.//blast_run')
        if domain_blast is not None:
            domain_hits = domain_blast.find('hits')
            if domain_hits is not None:
                domain_hit_elems = domain_hits.findall('hit')
                result['domain_blast_hit_count'] = len(domain_hit_elems)
                result['has_domain_blast_hits'] = len(domain_hit_elems) > 0
        
        # Combine blast hit results
        result['blast_hit_count'] = result.get('chain_blast_hit_count', 0) + result.get('domain_blast_hit_count', 0)
        result['has_blast_hits'] = result.get('has_chain_blast_hits', False) or result.get('has_domain_blast_hits', False)
        
        # Check for domains
        domains_elem = root.find('.//domains')
        if domains_elem is not None:
            domain_elems = domains_elem.findall('domain')
            result['domain_count'] = len(domain_elems)
            result['has_domains'] = len(domain_elems) > 0
            
            # Extract domain info
            if domain_elems:
                domains = []
                for domain in domain_elems:
                    domain_info = {
                        'id': domain.get('id'),
                        'range': domain.get('range')
                    }
                    
                    # Get classification if available
                    classification = domain.find('classification')
                    if classification is not None:
                        domain_info['t_group'] = classification.get('t_group')
                        domain_info['h_group'] = classification.get('h_group')
                        domain_info['x_group'] = classification.get('x_group')
                        domain_info['a_group'] = classification.get('a_group')
                    
                    domains.append(domain_info)
                
                result['domains'] = domains
                
                # Get domain classifications
                t_groups = Counter()
                h_groups = Counter()
                
                for domain in domains:
                    if 't_group' in domain and domain['t_group']:
                        t_groups[domain['t_group']] += 1
                    if 'h_group' in domain and domain['h_group']:
                        h_groups[domain['h_group']] += 1
                
                if t_groups:
                    result['t_group_counts'] = dict(t_groups)
                if h_groups:
                    result['h_group_counts'] = dict(h_groups)
        
        return result
    
    except ET.ParseError as e:
        result['error'] = f"XML parsing error: {str(e)}"
        return result
    except Exception as e:
        result['error'] = f"Error: {str(e)}"
        return result

def analyze_batch_domain_summaries(context, batch_id: int, sample_size: int = 50) -> Dict[str, Any]:
    """Analyze domain summaries across a batch"""
    logger = logging.getLogger('ecod.domain_analysis')
    
    # Get batch info
    batch_info = get_batch_info(context, batch_id)
    if not batch_info:
        logger.error(f"Batch {batch_id} not found")
        return {}
    
    base_path = batch_info['base_path']
    
    # Get a sample of proteins with domain summaries
    summary_query = """
    SELECT 
        p.pdb_id, p.chain_id, pf.file_path
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
    
    summaries = context.db.execute_dict_query(summary_query, (batch_id, sample_size))
    
    # Analyze each summary
    summary_results = []
    valid_xml_count = 0
    with_blast_hits = 0
    with_domains = 0
    no_hits_but_valid = 0
    error_count = 0
    
    for summary in summaries:
        file_path = os.path.join(base_path, summary['file_path'])
        result = analyze_domain_summary_xml(file_path)
        
        # Add protein info
        result['pdb_id'] = summary['pdb_id']
        result['chain_id'] = summary['chain_id']
        result['file_path'] = summary['file_path']
        
        summary_results.append(result)
        
        # Update counts
        if result['valid_xml']:
            valid_xml_count += 1
            
            if result['has_blast_hits']:
                with_blast_hits += 1
            else:
                no_hits_but_valid += 1
                
            if result['has_domains']:
                with_domains += 1
        
        if result['error']:
            error_count += 1
    
    # Get patterns for proteins missing domain summaries
    missing_patterns = get_failure_patterns(context, batch_id)
    
    # Compile results
    analysis_results = {
        'batch_id': batch_id,
        'batch_name': batch_info['batch_name'],
        'summary_sample_size': len(summaries),
        'valid_xml_count': valid_xml_count,
        'with_blast_hits': with_blast_hits,
        'with_domains': with_domains,
        'no_hits_but_valid': no_hits_but_valid,
        'error_count': error_count,
        'file_type_counts': get_file_type_counts(context, batch_id),
        'process_status_summary': get_process_status_summary(context, batch_id),
        'missing_summaries': missing_patterns,
        'sample_results': summary_results
    }
    
    return analysis_results

def identify_potential_issues(analysis_results: Dict[str, Any]) -> List[str]:
    """Identify potential issues based on analysis results"""
    issues = []
    
    # Check file type discrepancies
    file_counts = analysis_results.get('file_type_counts', {})
    
    chain_blast = file_counts.get('chain_blast_result', {}).get('exists', 0)
    domain_blast = file_counts.get('domain_blast_result', {}).get('exists', 0)
    domain_summary = file_counts.get('domain_summary', {}).get('exists', 0)
    domain_file = file_counts.get('domain_file', {}).get('exists', 0)
    
    if chain_blast > domain_summary:
        issues.append(f"Missing domain summaries: {chain_blast - domain_summary} proteins have chain BLAST results but no domain summary")
    
    if domain_blast > domain_summary:
        issues.append(f"Missing domain summaries: {domain_blast - domain_summary} proteins have domain BLAST results but no domain summary")
    
    if domain_summary > domain_file:
        issues.append(f"Missing domain files: {domain_summary - domain_file} proteins have domain summaries but no domain partition file")
    
    # Check for XML parsing errors
    if analysis_results.get('valid_xml_count', 0) < analysis_results.get('summary_sample_size', 0):
        invalid_count = analysis_results.get('summary_sample_size', 0) - analysis_results.get('valid_xml_count', 0)
        issues.append(f"XML parsing issues: {invalid_count} domain summary files could not be parsed correctly")
    
    # Check for summaries with no blast hits
    no_hits = analysis_results.get('no_hits_but_valid', 0)
    if no_hits > 0:
        issues.append(f"No BLAST hits: {no_hits} domain summaries contain no BLAST hits")
    
    # Check for error messages in process_status
    error_patterns = analysis_results.get('missing_summaries', {}).get('error_patterns', [])
    if error_patterns:
        for pattern in error_patterns[:3]:  # Top 3 error patterns
            issues.append(f"Error pattern: {pattern.get('error_message', 'Unknown error')} (affects {pattern.get('count', 0)} proteins)")
    
    return issues

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Analyze domain summary success and failure patterns')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to analyze')
    parser.add_argument('--sample-size', type=int, default=50,
                      help='Number of domain summaries to analyze')
    parser.add_argument('--output', type=str,
                      help='Output file for JSON results')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('--verbose', '-v', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger('ecod.domain_analysis')
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Analyze batch
    logger.info(f"Analyzing domain summaries for batch {args.batch_id}")
    analysis_results = analyze_batch_domain_summaries(context, args.batch_id, args.sample_size)
    
    if not analysis_results:
        logger.error(f"Failed to analyze batch {args.batch_id}")
        return 1
    
    # Identify potential issues
    issues = identify_potential_issues(analysis_results)
    analysis_results['identified_issues'] = issues
    
    # Display summary
    batch_name = analysis_results.get('batch_name', f'Batch {args.batch_id}')
    print(f"\nAnalysis Results for {batch_name} (ID: {args.batch_id})")
    print("-" * 70)
    
    # File counts
    file_counts = analysis_results.get('file_type_counts', {})
    print("File Counts:")
    for file_type, counts in file_counts.items():
        print(f"  {file_type}: {counts.get('exists', 0)}/{counts.get('total', 0)} files exist")
    
    # Domain summary sample analysis
    print("\nDomain Summary Analysis:")
    print(f"  Sample size: {analysis_results.get('summary_sample_size', 0)} summaries")
    print(f"  Valid XML: {analysis_results.get('valid_xml_count', 0)} summaries")
    print(f"  With BLAST hits: {analysis_results.get('with_blast_hits', 0)} summaries")
    print(f"  With domains: {analysis_results.get('with_domains', 0)} summaries")
    print(f"  No hits but valid XML: {analysis_results.get('no_hits_but_valid', 0)} summaries")
    
    # Missing summaries
    missing_count = analysis_results.get('missing_summaries', {}).get('missing_summary_count', 0)
    print(f"\nMissing Summaries: {missing_count} proteins have BLAST results but no domain summary")
    
    # Top error patterns
    error_patterns = analysis_results.get('missing_summaries', {}).get('error_patterns', [])
    if error_patterns:
        print("\nTop Error Patterns:")
        for i, pattern in enumerate(error_patterns[:5]):  # Top 5 error patterns
            print(f"  {i+1}. {pattern.get('error_message', 'Unknown error')}: {pattern.get('count', 0)} occurrences")
    
    # Identified issues
    if issues:
        print("\nIdentified Issues:")
        for i, issue in enumerate(issues):
            print(f"  {i+1}. {issue}")
    
    # Write output file if requested
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(analysis_results, f, indent=2)
        
        logger.info(f"Results written to {args.output}")
        print(f"\nDetailed results written to {args.output}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())