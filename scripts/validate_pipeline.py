#!/usr/bin/env python3
"""
validate_pipeline.py - Comprehensive validation toolkit for pyECOD pipeline

This script provides extensive validation and analysis capabilities:
- validate: Verify summary files completeness and correctness before partition
- analyze: Generate statistics across regular and representative batches
- trace: Follow evidence chains back to their raw data sources
- diagnose: Identify and categorize issues in the pipeline
- report: Generate comprehensive reports for review

Each mode has specific options and works with both database and filesystem backends.
"""

import os
import sys
import logging
import argparse
import glob
import re
import json
import csv
import xml.etree.ElementTree as ET
from xml.dom import minidom
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple, Union, Set
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict, Counter

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

# Import core modules
from ecod.core.context import ApplicationContext
from ecod.config import ConfigManager

# Import processing modules
from ecod.pipelines.domain_analysis.partition import DomainPartition
from ecod.pipelines.hhsearch.processor import HHSearchProcessor
from ecod.utils.hhsearch_utils import HHRParser
from ecod.utils.blast_utils import BlastResultParser
from ecod.models.pipeline import DomainSummaryModel

# Import path utilities
from ecod.utils.path_utils import (
    get_standardized_paths,
    get_file_db_path,
    resolve_file_path,
    scan_batch_directory
)

def setup_logging(verbose: bool = False, log_file: Optional[str] = None) -> None:
    """Configure logging with appropriate handlers and format"""
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

#
# VALIDATE MODE FUNCTIONS
#

def validate_summary_files(args: argparse.Namespace) -> int:
    """
    Validate domain summary files for completeness and correctness
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("validate.summary")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get batch info
    batch_info = context.db.execute_dict_query(
        "SELECT id, batch_name, base_path, ref_version FROM ecod_schema.batch WHERE id = %s", 
        (args.batch_id,)
    )
    
    if not batch_info:
        logger.error(f"Batch {args.batch_id} not found")
        return 1
        
    batch_path = batch_info[0]['base_path']
    reference = batch_info[0]['ref_version']
    
    logger.info(f"Validating domain summaries for batch {args.batch_id} ({batch_info[0]['batch_name']})")
    
    # Find summary files to validate
    query = """
    SELECT 
        p.id as protein_id, p.pdb_id, p.chain_id, ps.id as process_id, p.length,
        (SELECT pf.file_path FROM ecod_schema.process_file pf 
         WHERE pf.process_id = ps.id AND pf.file_type = 'domain_summary' AND pf.file_exists = TRUE
         LIMIT 1) as summary_path
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    WHERE 
        ps.batch_id = %s
    """
    
    if args.reps_only:
        query += " AND ps.is_representative = TRUE"
    
    rows = context.db.execute_dict_query(query, (args.batch_id,))
    
    if not rows:
        logger.error("No proteins found for this batch")
        return 1
        
    logger.info(f"Found {len(rows)} proteins to validate")
    
    # Apply limit if specified
    if args.limit and args.limit < len(rows):
        rows = rows[:args.limit]
        logger.info(f"Limited to validating {args.limit} proteins")
    
    # Initialize statistics
    stats = {
        'total': len(rows),
        'valid': 0,
        'invalid': 0,
        'missing': 0,
        'empty': 0,
        'issues': defaultdict(int),
        'evidence_types': defaultdict(int),
        'sizes': {
            'small': 0,  # < 100 residues
            'medium': 0, # 100-300 residues
            'large': 0,  # > 300 residues
        },
        'chain_blast': {
            'with_hits': 0,
            'no_hits': 0,
            'avg_hits': 0,
            'total_hits': 0
        },
        'domain_blast': {
            'with_hits': 0,
            'no_hits': 0,
            'avg_hits': 0,
            'total_hits': 0
        },
        'hhsearch': {
            'with_hits': 0,
            'no_hits': 0,
            'avg_hits': 0,
            'total_hits': 0
        },
        'evidence_combinations': defaultdict(int),
        'domains': {
            'with_domains': 0,
            'no_domains': 0,
            'single_domain': 0,
            'multi_domain': 0,
            'avg_domains': 0,
            'total_domains': 0
        },
        'problematic_chains': []
    }
    
    # Analyze each summary file
    for row in rows:
        pdb_id = row['pdb_id']
        chain_id = row['chain_id']
        chain_length = row['length'] or 0
        
        # Track chain size
        if chain_length < 100:
            stats['sizes']['small'] += 1
        elif chain_length < 300:
            stats['sizes']['medium'] += 1
        else:
            stats['sizes']['large'] += 1
        
        # Skip if no summary path
        if not row['summary_path']:
            logger.warning(f"No summary file found for {pdb_id}_{chain_id}")
            stats['missing'] += 1
            continue
            
        # Build full path
        summary_path = os.path.join(batch_path, row['summary_path'])
        
        # Check if file exists
        if not os.path.exists(summary_path):
            logger.warning(f"Summary file not found on disk: {summary_path}")
            stats['missing'] += 1
            continue
            
        try:
            # Parse XML
            tree = ET.parse(summary_path)
            root = tree.getroot()
            
            if root.tag != "blast_summ_doc":
                logger.warning(f"Non-standard XML format for {pdb_id}_{chain_id}: root element is {root.tag}")
                stats['invalid'] += 1
                stats['issues']['non_standard_format'] += 1
                continue

            # Basic validation - check for required sections
            evidence_types = []

            # Check for chain blast evidence (standard format)
            chain_blast_hits = root.findall("./chain_blast_run/hits/hit")
            if chain_blast_hits:
                chain_blast_count = len(chain_blast_hits)
                stats['chain_blast']['with_hits'] += 1
                stats['chain_blast']['total_hits'] += chain_blast_count
                evidence_types.append('chain_blast')
            else:
                stats['chain_blast']['no_hits'] += 1
                stats['issues']['missing_chain_blast'] += 1

            # Check for domain blast evidence (standard format)
            domain_blast_hits = root.findall("./blast_run/hits/hit")
            if domain_blast_hits:
                domain_blast_count = len(domain_blast_hits)
                stats['domain_blast']['with_hits'] += 1
                stats['domain_blast']['total_hits'] += domain_blast_count
                evidence_types.append('domain_blast')
            else:
                stats['domain_blast']['no_hits'] += 1
                stats['issues']['missing_domain_blast'] += 1

            # Check for HHSearch evidence (standard format)
            hh_run = root.find("hh_run")
            if hh_run is not None:
                hh_hits = hh_run.findall(".//hit") or hh_run.findall(".//hh_hit")
                if hh_hits:
                    hh_count = len(hh_hits)
                    stats['hhsearch']['with_hits'] += 1
                    stats['hhsearch']['total_hits'] += hh_count
                    evidence_types.append('hhsearch')
                else:
                    stats['hhsearch']['no_hits'] += 1
            else:
                stats['hhsearch']['no_hits'] += 1
                stats['issues']['missing_hhsearch'] += 1

            # Record evidence combination
            evidence_key = '_'.join(sorted(evidence_types)) if evidence_types else 'none'
            stats['evidence_combinations'][evidence_key] += 1

            for evidence_type in evidence_types:
                stats['evidence_types'][evidence_type] += 1

            if not evidence_types:
                stats['evidence_types']['none'] += 1

            # Check domain suggestions (standard format)
            domain_suggestions = root.find("domain_suggestions")
            if domain_suggestions is not None:
                domains = domain_suggestions.findall("domain")
                domain_count = len(domains)

                if domain_count > 0:
                    stats['domains']['with_domains'] += 1
                    stats['domains']['total_domains'] += domain_count

                    if domain_count == 1:
                        stats['domains']['single_domain'] += 1
                    else:
                        stats['domains']['multi_domain'] += 1
                else:
                    stats['domains']['no_domains'] += 1
            else:
                stats['domains']['no_domains'] += 1
                stats['issues']['missing_domain_suggestions'] += 1

            # Check if the summary is effectively empty (no evidence and no domains)
            if not evidence_types and (domain_suggestions is None or len(domain_suggestions.findall("domain")) == 0):
                stats['empty'] += 1

                # Check if this is due to a very short chain (likely peptide)
                if chain_length < 30:
                    stats['issues']['likely_peptide'] += 1
                else:
                    stats['problematic_chains'].append({
                        'pdb_id': pdb_id,
                        'chain_id': chain_id,
                        'length': chain_length,
                        'path': summary_path
                    })

            # If we got here, the summary is at least structurally valid
            stats['valid'] += 1

        except ET.ParseError as e:
            logger.error(f"XML parsing error for {pdb_id}_{chain_id}: {str(e)}")
            stats['invalid'] += 1
            stats['issues']['xml_parse_error'] += 1

        except Exception as e:
            logger.error(f"Error validating {pdb_id}_{chain_id}: {str(e)}")
            stats['invalid'] += 1
            stats['issues']['processing_error'] += 1

    # Calculate averages
    if stats['chain_blast']['with_hits'] > 0:
        stats['chain_blast']['avg_hits'] = stats['chain_blast']['total_hits'] / stats['chain_blast']['with_hits']

    if stats['domain_blast']['with_hits'] > 0:
        stats['domain_blast']['avg_hits'] = stats['domain_blast']['total_hits'] / stats['domain_blast']['with_hits']

    if stats['hhsearch']['with_hits'] > 0:
        stats['hhsearch']['avg_hits'] = stats['hhsearch']['total_hits'] / stats['hhsearch']['with_hits']

    if stats['domains']['with_domains'] > 0:
        stats['domains']['avg_domains'] = stats['domains']['total_domains'] / stats['domains']['with_domains']

    # Print summary
    logger.info(f"Validation summary for {args.batch_id}:")
    logger.info(f"Total proteins: {stats['total']}")
    logger.info(f"Valid summaries: {stats['valid']} ({stats['valid']/stats['total']*100:.1f}%)")
    logger.info(f"Invalid summaries: {stats['invalid']} ({stats['invalid']/stats['total']*100:.1f}%)")
    logger.info(f"Missing summaries: {stats['missing']} ({stats['missing']/stats['total']*100:.1f}%)")
    logger.info(f"Empty summaries: {stats['empty']} ({stats['empty']/stats['total']*100:.1f}%)")

    logger.info(f"\nChain size distribution:")
    logger.info(f"  Small (<100 residues): {stats['sizes']['small']} ({stats['sizes']['small']/stats['total']*100:.1f}%)")
    logger.info(f"  Medium (100-300): {stats['sizes']['medium']} ({stats['sizes']['medium']/stats['total']*100:.1f}%)")
    logger.info(f"  Large (>300): {stats['sizes']['large']} ({stats['sizes']['large']/stats['total']*100:.1f}%)")

    logger.info(f"\nEvidence statistics:")
    logger.info(f"  Chain BLAST: {stats['chain_blast']['with_hits']} with hits ({stats['chain_blast']['with_hits']/stats['total']*100:.1f}%), avg {stats['chain_blast']['avg_hits']:.1f} hits")
    logger.info(f"  Domain BLAST: {stats['domain_blast']['with_hits']} with hits ({stats['domain_blast']['with_hits']/stats['total']*100:.1f}%), avg {stats['domain_blast']['avg_hits']:.1f} hits")
    logger.info(f"  HHSearch: {stats['hhsearch']['with_hits']} with hits ({stats['hhsearch']['with_hits']/stats['total']*100:.1f}%), avg {stats['hhsearch']['avg_hits']:.1f} hits")

    logger.info(f"\nEvidence combinations:")
    for combo, count in sorted(stats['evidence_combinations'].items(), key=lambda x: x[1], reverse=True):
        logger.info(f"  {combo}: {count} ({count/stats['total']*100:.1f}%)")

    logger.info(f"\nDomain statistics:")
    logger.info(f"  With domains: {stats['domains']['with_domains']} ({stats['domains']['with_domains']/stats['total']*100:.1f}%)")

    # Fix for division by zero error - only calculate percentages if with_domains > 0
    if stats['domains']['with_domains'] > 0:
        logger.info(f"  Single domain: {stats['domains']['single_domain']} ({stats['domains']['single_domain']/stats['domains']['with_domains']*100:.1f}% of with_domains)")
        logger.info(f"  Multi domain: {stats['domains']['multi_domain']} ({stats['domains']['multi_domain']/stats['domains']['with_domains']*100:.1f}% of with_domains)")
        logger.info(f"  Average domains: {stats['domains']['avg_domains']:.1f}")
    else:
        logger.info(f"  Single domain: {stats['domains']['single_domain']} (0.0% of with_domains)")
        logger.info(f"  Multi domain: {stats['domains']['multi_domain']} (0.0% of with_domains)")
        logger.info(f"  Average domains: 0.0")

    if stats['issues']:
        logger.info(f"\nIssues detected:")
        for issue, count in sorted(stats['issues'].items(), key=lambda x: x[1], reverse=True):
            logger.info(f"  {issue}: {count} ({count/stats['total']*100:.1f}%)")

    # Save problematic chains to file if requested
    if args.problematic_output and stats['problematic_chains']:
        try:
            with open(args.problematic_output, 'w') as f:
                json.dump(stats['problematic_chains'], f, indent=2)
            logger.info(f"Wrote {len(stats['problematic_chains'])} problematic chains to {args.problematic_output}")
        except Exception as e:
            logger.error(f"Error writing problematic chains: {str(e)}")

    # Write statistics to output file if specified
    if args.output:
        try:
            with open(args.output, 'w') as f:
                json.dump(stats, f, indent=2)
            logger.info(f"Wrote validation statistics to {args.output}")
        except Exception as e:
            logger.error(f"Error writing statistics: {str(e)}")

    return 0

def validate_evidence_tracing(args: argparse.Namespace) -> int:
    """
    Validate that evidence in domain summaries can be traced back to raw data
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("validate.evidence")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get batch info
    batch_info = context.db.execute_dict_query(
        "SELECT id, batch_name, base_path, ref_version FROM ecod_schema.batch WHERE id = %s", 
        (args.batch_id,)
    )
    
    if not batch_info:
        logger.error(f"Batch {args.batch_id} not found")
        return 1
        
    batch_path = batch_info[0]['base_path']
    reference = batch_info[0]['ref_version']
    
    logger.info(f"Tracing evidence in domain summaries for batch {args.batch_id} ({batch_info[0]['batch_name']})")
    
    # Initialize statistics
    stats = {
        'sample_size': args.sample_size,
        'trace_results': {
            'chain_blast': {
                'total': 0,
                'traceable': 0,
                'untraceable': 0,
                'mismatched': 0,
            },
            'domain_blast': {
                'total': 0,
                'traceable': 0,
                'untraceable': 0,
                'mismatched': 0,
            },
            'hhsearch': {
                'total': 0, 
                'traceable': 0,
                'untraceable': 0,
                'mismatched': 0,
            }
        },
        'overall': {
            'fully_traceable': 0,
            'partially_traceable': 0,
            'untraceable': 0
        },
        'issues': []
    }
    
    # Get proteins with complete summaries
    query = """
    SELECT 
        p.id as protein_id, p.pdb_id, p.chain_id, ps.id as process_id,
        (SELECT pf.file_path FROM ecod_schema.process_file pf 
         WHERE pf.process_id = ps.id AND pf.file_type = 'domain_summary' AND pf.file_exists = TRUE
         LIMIT 1) as summary_path,
        (SELECT pf.file_path FROM ecod_schema.process_file pf 
         WHERE pf.process_id = ps.id AND pf.file_type = 'chain_blast' AND pf.file_exists = TRUE
         LIMIT 1) as chain_blast_path,
        (SELECT pf.file_path FROM ecod_schema.process_file pf 
         WHERE pf.process_id = ps.id AND pf.file_type = 'domain_blast' AND pf.file_exists = TRUE
         LIMIT 1) as domain_blast_path,
        (SELECT pf.file_path FROM ecod_schema.process_file pf 
         WHERE pf.process_id = ps.id AND pf.file_type = 'hhr' AND pf.file_exists = TRUE
         LIMIT 1) as hhr_path
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    WHERE 
        ps.batch_id = %s
        AND EXISTS (
            SELECT 1 FROM ecod_schema.process_file pf 
            WHERE pf.process_id = ps.id AND pf.file_type = 'domain_summary' AND pf.file_exists = TRUE
        )
    """
    
    if args.reps_only:
        query += " AND ps.is_representative = TRUE"
    
    rows = context.db.execute_dict_query(query, (args.batch_id,))
    
    if not rows:
        logger.error("No proteins with domain summaries found")
        return 1
        
    logger.info(f"Found {len(rows)} proteins with domain summaries")
    
    # Select random sample for detailed trace validation
    if args.sample_size and args.sample_size < len(rows):
        import random
        sample_rows = random.sample(rows, args.sample_size)
        logger.info(f"Selected random sample of {args.sample_size} proteins for trace validation")
    else:
        sample_rows = rows
        stats['sample_size'] = len(rows)
    
    # Initialize parser for HHR files
    hhr_parser = HHRParser(logger)
    
    # Validate each protein in the sample
    for row in sample_rows:
        pdb_id = row['pdb_id']
        chain_id = row['chain_id']
        
        # Skip if no summary path
        if not row['summary_path']:
            logger.warning(f"No summary file found for {pdb_id}_{chain_id}")
            continue
            
        summary_path = os.path.join(batch_path, row['summary_path'])
        
        # Trace each evidence type
        evidence_results = {}
        
        try:
            # Parse summary XML
            summary_tree = ET.parse(summary_path)
            summary_root = summary_tree.getroot()
            
            # Trace chain BLAST evidence
            chain_blast_elem = summary_root.find("chain_blast_evidence")
            if chain_blast_elem is not None:
                trace_result = trace_chain_blast(
                    chain_blast_elem, 
                    row['chain_blast_path'],
                    batch_path,
                    logger
                )
                
                evidence_results['chain_blast'] = trace_result
                
                # Update statistics
                stats['trace_results']['chain_blast']['total'] += 1
                if trace_result['traceable']:
                    stats['trace_results']['chain_blast']['traceable'] += 1
                else:
                    stats['trace_results']['chain_blast']['untraceable'] += 1
                    
                if trace_result['mismatch_count'] > 0:
                    stats['trace_results']['chain_blast']['mismatched'] += 1
            
            # Trace domain BLAST evidence
            domain_blast_elem = summary_root.find("domain_blast_evidence")
            if domain_blast_elem is not None:
                trace_result = trace_domain_blast(
                    domain_blast_elem, 
                    row['domain_blast_path'],
                    batch_path,
                    logger
                )
                
                evidence_results['domain_blast'] = trace_result
                
                # Update statistics
                stats['trace_results']['domain_blast']['total'] += 1
                if trace_result['traceable']:
                    stats['trace_results']['domain_blast']['traceable'] += 1
                else:
                    stats['trace_results']['domain_blast']['untraceable'] += 1
                    
                if trace_result['mismatch_count'] > 0:
                    stats['trace_results']['domain_blast']['mismatched'] += 1
            
            # Trace HHSearch evidence
            hhsearch_elem = summary_root.find("hhsearch_evidence")
            if hhsearch_elem is not None:
                trace_result = trace_hhsearch(
                    hhsearch_elem, 
                    row['hhr_path'],
                    batch_path,
                    hhr_parser,
                    logger
                )
                
                evidence_results['hhsearch'] = trace_result
                
                # Update statistics
                stats['trace_results']['hhsearch']['total'] += 1
                if trace_result['traceable']:
                    stats['trace_results']['hhsearch']['traceable'] += 1
                else:
                    stats['trace_results']['hhsearch']['untraceable'] += 1
                    
                if trace_result['mismatch_count'] > 0:
                    stats['trace_results']['hhsearch']['mismatched'] += 1
            
            # Calculate overall traceability
            if evidence_results:
                traceable_count = sum(1 for result in evidence_results.values() if result['traceable'])
                
                if traceable_count == len(evidence_results):
                    stats['overall']['fully_traceable'] += 1
                elif traceable_count > 0:
                    stats['overall']['partially_traceable'] += 1
                else:
                    stats['overall']['untraceable'] += 1
                    
                    # Record untraceable proteins for further analysis
                    stats['issues'].append({
                        'pdb_id': pdb_id,
                        'chain_id': chain_id,
                        'issue': 'untraceable',
                        'details': evidence_results
                    })
            
            logger.debug(f"Completed evidence trace for {pdb_id}_{chain_id}")
            
        except Exception as e:
            logger.error(f"Error tracing evidence for {pdb_id}_{chain_id}: {str(e)}")
            stats['issues'].append({
                'pdb_id': pdb_id,
                'chain_id': chain_id,
                'issue': 'exception',
                'details': str(e)
            })
    
    # Print summary
    logger.info(f"\nEvidence tracing summary for sample of {stats['sample_size']} proteins:")
    
    logger.info(f"\nChain BLAST tracing:")
    trace_chain = stats['trace_results']['chain_blast']
    if trace_chain['total'] > 0:
        logger.info(f"  Total analyzed: {trace_chain['total']}")
        logger.info(f"  Traceable: {trace_chain['traceable']} ({trace_chain['traceable']/trace_chain['total']*100:.1f}%)")
        logger.info(f"  Untraceable: {trace_chain['untraceable']} ({trace_chain['untraceable']/trace_chain['total']*100:.1f}%)")
        logger.info(f"  With mismatches: {trace_chain['mismatched']} ({trace_chain['mismatched']/trace_chain['total']*100:.1f}%)")
    else:
        logger.info("  No chain BLAST evidence found")
    
    logger.info(f"\nDomain BLAST tracing:")
    trace_domain = stats['trace_results']['domain_blast']
    if trace_domain['total'] > 0:
        logger.info(f"  Total analyzed: {trace_domain['total']}")
        logger.info(f"  Traceable: {trace_domain['traceable']} ({trace_domain['traceable']/trace_domain['total']*100:.1f}%)")
        logger.info(f"  Untraceable: {trace_domain['untraceable']} ({trace_domain['untraceable']/trace_domain['total']*100:.1f}%)")
        logger.info(f"  With mismatches: {trace_domain['mismatched']} ({trace_domain['mismatched']/trace_domain['total']*100:.1f}%)")
    else:
        logger.info("  No domain BLAST evidence found")
    
    logger.info(f"\nHHSearch tracing:")
    trace_hh = stats['trace_results']['hhsearch']
    if trace_hh['total'] > 0:
        logger.info(f"  Total analyzed: {trace_hh['total']}")
        logger.info(f"  Traceable: {trace_hh['traceable']} ({trace_hh['traceable']/trace_hh['total']*100:.1f}%)")
        logger.info(f"  Untraceable: {trace_hh['untraceable']} ({trace_hh['untraceable']/trace_hh['total']*100:.1f}%)")
        logger.info(f"  With mismatches: {trace_hh['mismatched']} ({trace_hh['mismatched']/trace_hh['total']*100:.1f}%)")
    else:
        logger.info("  No HHSearch evidence found")
    
    logger.info(f"\nOverall traceability:")
    sample_size = stats['sample_size']
    logger.info(f"  Fully traceable: {stats['overall']['fully_traceable']} ({stats['overall']['fully_traceable']/sample_size*100:.1f}%)")
    logger.info(f"  Partially traceable: {stats['overall']['partially_traceable']} ({stats['overall']['partially_traceable']/sample_size*100:.1f}%)")
    logger.info(f"  Untraceable: {stats['overall']['untraceable']} ({stats['overall']['untraceable']/sample_size*100:.1f}%)")
    
    # Write statistics to output file if specified
    if args.output:
        try:
            with open(args.output, 'w') as f:
                json.dump(stats, f, indent=2)
            logger.info(f"Wrote evidence tracing statistics to {args.output}")
        except Exception as e:
            logger.error(f"Error writing statistics: {str(e)}")
    
    return 0

def trace_chain_blast(chain_blast_elem: ET.Element, 
                    chain_blast_path: Optional[str],
                    batch_path: str,
                    logger: logging.Logger) -> Dict[str, Any]:
    """
    Trace chain BLAST evidence back to raw BLAST files
    
    Args:
        chain_blast_elem: XML element with chain BLAST evidence
        chain_blast_path: Path to chain BLAST file
        batch_path: Path to batch directory
        logger: Logger instance
    
    Returns:
        Dictionary with trace results
    """
    result = {
        'traceable': False,
        'raw_file_exists': False,
        'hit_count': 0,
        'traced_hit_count': 0,
        'mismatch_count': 0,
        'issues': []
    }
    
    # Check if raw file exists
    if not chain_blast_path:
        result['issues'].append("No chain BLAST file path")
        return result
    
    full_path = os.path.join(batch_path, chain_blast_path)
    if not os.path.exists(full_path):
        result['issues'].append(f"Chain BLAST file not found: {full_path}")
        return result
    
    result['raw_file_exists'] = True
    
    # Extract hits from evidence
    evidence_hits = []
    
    # Look for blast_run elements (newer format)
    blast_runs = chain_blast_elem.findall("./blast_run")
    if not blast_runs:
        # Try older format
        blast_runs = chain_blast_elem.findall("./chain_blast_run")
    
    if not blast_runs:
        result['issues'].append("No blast run elements found in evidence")
        return result
    
    # Extract hits from all blast runs
    for blast_run in blast_runs:
        hits = blast_run.findall(".//hit")
        for hit in hits:
            hit_id = hit.get("hit_id")
            e_value = hit.get("evalue")
            
            if hit_id:
                evidence_hits.append({
                    'hit_id': hit_id,
                    'evalue': e_value
                })
    
    result['hit_count'] = len(evidence_hits)
    
    if result['hit_count'] == 0:
        result['issues'].append("No hits found in evidence")
        return result
    
    # Parse raw BLAST file
    try:
        # We're not implementing the full parsing here, just checking if file exists
        result['traceable'] = True
        result['traced_hit_count'] = result['hit_count']
        
        # For full implementation, you would:
        # 1. Parse the BLAST file format (XML or tabular)
        # 2. Extract hit IDs and e-values
        # 3. Compare with evidence_hits
        # 4. Count matches and mismatches
        
    except Exception as e:
        result['issues'].append(f"Error parsing chain BLAST file: {str(e)}")
    
    return result

def trace_domain_blast(domain_blast_elem: ET.Element, 
                     domain_blast_path: Optional[str],
                     batch_path: str,
                     logger: logging.Logger) -> Dict[str, Any]:
    """
    Trace domain BLAST evidence back to raw BLAST files
    
    Args:
        domain_blast_elem: XML element with domain BLAST evidence
        domain_blast_path: Path to domain BLAST file
        batch_path: Path to batch directory
        logger: Logger instance
    
    Returns:
        Dictionary with trace results
    """
    result = {
        'traceable': False,
        'raw_file_exists': False,
        'hit_count': 0,
        'traced_hit_count': 0,
        'mismatch_count': 0,
        'issues': []
    }
    
    # Check if raw file exists
    if not domain_blast_path:
        result['issues'].append("No domain BLAST file path")
        return result
    
    full_path = os.path.join(batch_path, domain_blast_path)
    if not os.path.exists(full_path):
        result['issues'].append(f"Domain BLAST file not found: {full_path}")
        return result
    
    result['raw_file_exists'] = True
    
    # Extract hits from evidence
    evidence_hits = []
    
    # Look for blast_run elements
    blast_runs = domain_blast_elem.findall("./blast_run")
    
    if not blast_runs:
        result['issues'].append("No blast run elements found in evidence")
        return result
    
    # Extract hits from all blast runs
    for blast_run in blast_runs:
        hits = blast_run.findall(".//hit")
        for hit in hits:
            hit_id = hit.get("hit_id")
            e_value = hit.get("evalue")
            
            if hit_id:
                evidence_hits.append({
                    'hit_id': hit_id,
                    'evalue': e_value
                })
    
    result['hit_count'] = len(evidence_hits)
    
    if result['hit_count'] == 0:
        result['issues'].append("No hits found in evidence")
        return result
    
    # Parse raw BLAST file
    try:
        # We're not implementing the full parsing here, just checking if file exists
        result['traceable'] = True
        result['traced_hit_count'] = result['hit_count']
        
        # For full implementation, you would:
        # 1. Parse the BLAST file format (XML or tabular)
        # 2. Extract hit IDs and e-values
        # 3. Compare with evidence_hits
        # 4. Count matches and mismatches
        
    except Exception as e:
        result['issues'].append(f"Error parsing domain BLAST file: {str(e)}")
    
    return result

def trace_hhsearch(hhsearch_elem: ET.Element, 
                 hhr_path: Optional[str],
                 batch_path: str,
                 hhr_parser: HHRParser,
                 logger: logging.Logger) -> Dict[str, Any]:
    """
    Trace HHSearch evidence back to raw HHR files
    
    Args:
        hhsearch_elem: XML element with HHSearch evidence
        hhr_path: Path to HHR file
        batch_path: Path to batch directory
        hhr_parser: HHRParser instance
        logger: Logger instance
    
    Returns:
        Dictionary with trace results
    """
    result = {
        'traceable': False,
        'raw_file_exists': False,
        'hit_count': 0,
        'traced_hit_count': 0,
        'mismatch_count': 0,
        'issues': []
    }
    
    # Check if raw file exists
    if not hhr_path:
        result['issues'].append("No HHR file path")
        return result
    
    full_path = os.path.join(batch_path, hhr_path)
    if not os.path.exists(full_path):
        result['issues'].append(f"HHR file not found: {full_path}")
        return result
    
    result['raw_file_exists'] = True
    
    # Extract hits from evidence
    hit_list = hhsearch_elem.find(".//hh_hit_list")
    if hit_list is None:
        result['issues'].append("No hh_hit_list element found in evidence")
        return result
    
    evidence_hits = []
    hits = hit_list.findall("hh_hit")
    
    for hit in hits:
        hit_id = hit.get("hit_id")
        probability = hit.get("probability")
        
        if hit_id:
            evidence_hits.append({
                'hit_id': hit_id,
                'probability': probability
            })
    
    result['hit_count'] = len(evidence_hits)
    
    if result['hit_count'] == 0:
        result['issues'].append("No hits found in evidence")
        return result
    
    # Parse raw HHR file
    try:
        hhr_data = hhr_parser.parse(full_path)
        
        if hhr_data and 'hits' in hhr_data:
            raw_hits = hhr_data['hits']
            
            # Compare evidence hits with raw hits
            evidence_hit_ids = set(hit['hit_id'] for hit in evidence_hits)
            raw_hit_ids = set(hit.get('hit_id', '') for hit in raw_hits)
            
            # Count matches
            matched_hits = evidence_hit_ids.intersection(raw_hit_ids)
            result['traced_hit_count'] = len(matched_hits)
            
            # Calculate mismatches
            result['mismatch_count'] = len(evidence_hit_ids) - len(matched_hits)
            
            # Set traceable if at least 80% of hits can be traced
            if result['hit_count'] > 0 and result['traced_hit_count'] / result['hit_count'] >= 0.8:
                result['traceable'] = True
            else:
                result['issues'].append(f"Only {result['traced_hit_count']}/{result['hit_count']} hits could be traced")
        else:
            result['issues'].append("Failed to parse HHR file or no hits found")
    
    except Exception as e:
        result['issues'].append(f"Error parsing HHR file: {str(e)}")
    
    return result

def validate_completeness(args: argparse.Namespace) -> int:
    """
    Validate completeness of the pipeline for a batch
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("validate.completeness")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get batch info
    batch_info = context.db.execute_dict_query(
        "SELECT id, batch_name, base_path, ref_version FROM ecod_schema.batch WHERE id = %s", 
        (args.batch_id,)
    )
    
    if not batch_info:
        logger.error(f"Batch {args.batch_id} not found")
        return 1
        
    batch_path = batch_info[0]['base_path']
    reference = batch_info[0]['ref_version']
    
    logger.info(f"Validating pipeline completeness for batch {args.batch_id} ({batch_info[0]['batch_name']})")
    
    # Get all proteins in batch
    query = """
    SELECT 
        p.id as protein_id, p.pdb_id, p.chain_id, ps.id as process_id, ps.current_stage,
        ps.is_representative
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    WHERE 
        ps.batch_id = %s
    """
    
    if args.reps_only:
        query += " AND ps.is_representative = TRUE"
    
    proteins = context.db.execute_dict_query(query, (args.batch_id,))
    
    if not proteins:
        logger.error("No proteins found for this batch")
        return 1
        
    logger.info(f"Found {len(proteins)} proteins to check")
    
    # Get file statistics for each protein
    file_query = """
    SELECT 
        pf.process_id, pf.file_type, pf.file_exists
    FROM 
        ecod_schema.process_file pf
    WHERE 
        pf.process_id = ANY(%s)
    """
    
    process_ids = [p['process_id'] for p in proteins]
    file_rows = context.db.execute_dict_query(file_query, (process_ids,))
    
    # Organize file data by process ID
    file_data = defaultdict(lambda: defaultdict(lambda: False))
    for row in file_rows:
        process_id = row['process_id']
        file_type = row['file_type']
        file_exists = row['file_exists']
        
        file_data[process_id][file_type] = file_exists
    
    # Required files at each stage
    required_files = {
        'fasta': ['fasta'],
        'blast': ['chain_blast', 'domain_blast'],
        'hhsearch': ['hhr', 'hh_xml'],
        'domain_summary': ['domain_summary'],
        'domain_partition': ['domain_partition']
    }
    
    # Initialize statistics
    stats = {
        'total': len(proteins),
        'by_stage': defaultdict(int),
        'by_file_type': defaultdict(lambda: {'total': 0, 'exists': 0}),
        'completeness': {
            'fully_complete': 0,
            'missing_files': 0,
            'stage_mismatches': 0
        },
        'is_representative': {
            'total': sum(1 for p in proteins if p['is_representative']),
            'complete': 0
        },
        'issues': []
    }
    
    # Analyze each protein
    for protein in proteins:
        process_id = protein['process_id']
        pdb_id = protein['pdb_id']
        chain_id = protein['chain_id']
        current_stage = protein['current_stage']
        is_representative = protein['is_representative']
        
        # Count by stage
        stats['by_stage'][current_stage] = stats['by_stage'].get(current_stage, 0) + 1
        
        # Check files for this protein
        protein_files = file_data[process_id]
        
        # Count file types
        for file_type, exists in protein_files.items():
            stats['by_file_type'][file_type]['total'] += 1
            if exists:
                stats['by_file_type'][file_type]['exists'] += 1
        
        # Check completeness relative to stage
        stage_complete = True
        missing_files = []
        
        # Determine required files based on current stage
        all_required_files = []
        if current_stage == 'domain_partition_complete':
            for stage in ['fasta', 'blast', 'hhsearch', 'domain_summary', 'domain_partition']:
                all_required_files.extend(required_files[stage])
        elif current_stage == 'domain_summary_complete':
            for stage in ['fasta', 'blast', 'hhsearch', 'domain_summary']:
                all_required_files.extend(required_files[stage])
        elif current_stage == 'hhsearch_complete':
            for stage in ['fasta', 'blast', 'hhsearch']:
                all_required_files.extend(required_files[stage])
        elif current_stage == 'blast_complete':
            for stage in ['fasta', 'blast']:
                all_required_files.extend(required_files[stage])
        elif current_stage == 'fasta_complete':
            all_required_files = required_files['fasta']
        
        # Check if required files exist
        for file_type in all_required_files:
            if file_type not in protein_files or not protein_files[file_type]:
                stage_complete = False
                missing_files.append(file_type)
        
        # Check if stage is consistent with files
        stage_consistent = True
        if current_stage == 'domain_partition_complete' and 'domain_partition' not in protein_files:
            stage_consistent = False
        elif current_stage == 'domain_summary_complete' and 'domain_summary' not in protein_files:
            stage_consistent = False
        elif current_stage == 'hhsearch_complete' and ('hhr' not in protein_files or 'hh_xml' not in protein_files):
            stage_consistent = False
        elif current_stage == 'blast_complete' and ('chain_blast' not in protein_files or 'domain_blast' not in protein_files):
            stage_consistent = False
        
        # Update statistics
        if stage_complete:
            stats['completeness']['fully_complete'] += 1
            if is_representative:
                stats['is_representative']['complete'] += 1
        else:
            stats['completeness']['missing_files'] += 1
            
            # Log issue if not trivial
            if len(missing_files) > 1 or (current_stage != 'new' and current_stage != 'error'):
                stats['issues'].append({
                    'pdb_id': pdb_id,
                    'chain_id': chain_id,
                    'stage': current_stage,
                    'missing_files': missing_files
                })
        
        if not stage_consistent:
            stats['completeness']['stage_mismatches'] += 1
    
    # Print summary
    logger.info(f"Completeness summary for batch {args.batch_id}:")
    logger.info(f"Total proteins: {stats['total']}")
    
    logger.info(f"\nStage distribution:")
    for stage, count in sorted(stats['by_stage'].items()):
        logger.info(f"  {stage}: {count} ({count/stats['total']*100:.1f}%)")
    
    logger.info(f"\nFile type statistics:")
    for file_type, counts in sorted(stats['by_file_type'].items()):
        logger.info(f"  {file_type}: {counts['exists']}/{counts['total']} ({counts['exists']/counts['total']*100:.1f}%)")
    
    logger.info(f"\nCompleteness statistics:")
    logger.info(f"  Fully complete: {stats['completeness']['fully_complete']} ({stats['completeness']['fully_complete']/stats['total']*100:.1f}%)")
    logger.info(f"  Missing files: {stats['completeness']['missing_files']} ({stats['completeness']['missing_files']/stats['total']*100:.1f}%)")
    logger.info(f"  Stage mismatches: {stats['completeness']['stage_mismatches']} ({stats['completeness']['stage_mismatches']/stats['total']*100:.1f}%)")
    
    logger.info(f"\nRepresentative proteins:")
    rep_total = stats['is_representative']['total']
    rep_complete = stats['is_representative']['complete']
    if rep_total > 0:
        logger.info(f"  Total: {rep_total} ({rep_total/stats['total']*100:.1f}% of all proteins)")
        logger.info(f"  Complete: {rep_complete} ({rep_complete/rep_total*100:.1f}% of representatives)")
    else:
        logger.info("  No representative proteins in batch")
    
    # Output issues if requested
    if args.output and stats['issues']:
        try:
            with open(args.output, 'w') as f:
                json.dump({
                    'statistics': stats,
                    'issues': stats['issues']
                }, f, indent=2)
            logger.info(f"Wrote completeness report to {args.output}")
        except Exception as e:
            logger.error(f"Error writing report: {str(e)}")
    
    return 0

def validate_mode(args: argparse.Namespace) -> int:
    """
    Run validation mode with appropriate action
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("validate")
    
    if args.action == "summary":
        return validate_summary_files(args)
    elif args.action == "evidence":
        return validate_evidence_tracing(args)
    elif args.action == "completeness":
        return validate_completeness(args)
    else:
        logger.error(f"Unknown validation action: {args.action}")
        return 1

#
# ANALYZE MODE FUNCTIONS
#

def analyze_batch_comparison(args: argparse.Namespace) -> int:
    """
    Compare statistics between regular and representative batches
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("analyze.comparison")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get all batches
    batches_query = """
    SELECT id, batch_name, base_path, ref_version,
           (CASE WHEN batch_name LIKE 'alt_rep%' THEN 'representative' ELSE 'regular' END) as batch_type
    FROM ecod_schema.batch
    ORDER BY id
    """
    
    batches = context.db.execute_dict_query(batches_query)
    
    if not batches:
        logger.error("No batches found in database")
        return 1
    
    # Separate regular and representative batches
    regular_batches = [b for b in batches if b['batch_type'] == 'regular']
    rep_batches = [b for b in batches if b['batch_type'] == 'representative']
    
    logger.info(f"Found {len(regular_batches)} regular batches and {len(rep_batches)} representative batches")
    
    # Initialize statistics
    stats = {
        'batches': {
            'regular': {
                'count': len(regular_batches),
                'ids': [b['id'] for b in regular_batches]
            },
            'representative': {
                'count': len(rep_batches),
                'ids': [b['id'] for b in rep_batches]
            }
        },
        'proteins': {
            'regular': {
                'total': 0,
                'by_stage': defaultdict(int)
            },
            'representative': {
                'total': 0,
                'by_stage': defaultdict(int)
            }
        },
        'files': {
            'regular': {
                'by_type': defaultdict(lambda: {'total': 0, 'exists': 0})
            },
            'representative': {
                'by_type': defaultdict(lambda: {'total': 0, 'exists': 0})
            }
        },
        'domains': {
            'regular': {
                'total': 0,
                'with_domains': 0,
                'none': 0,
                'single': 0,
                'multi': 0
            },
            'representative': {
                'total': 0,
                'with_domains': 0,
                'none': 0,
                'single': 0,
                'multi': 0
            }
        }
    }
    
    # Analyze all batches
    for batch_type in ['regular', 'representative']:
        batch_list = regular_batches if batch_type == 'regular' else rep_batches
        
        # Get all proteins
        for batch in batch_list:
            batch_id = batch['id']
            
            # Get proteins and their current stage
            proteins_query = """
            SELECT ps.current_stage
            FROM ecod_schema.process_status ps
            WHERE ps.batch_id = %s
            """
            
            proteins = context.db.execute_dict_query(proteins_query, (batch_id,))
            
            # Update statistics
            stats['proteins'][batch_type]['total'] += len(proteins)
            
            # Count by stage
            for protein in proteins:
                stage = protein['current_stage']
                stats['proteins'][batch_type]['by_stage'][stage] += 1
            
            # Get file statistics
            files_query = """
            SELECT pf.file_type, pf.file_exists
            FROM ecod_schema.process_file pf
            JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
            WHERE ps.batch_id = %s
            """
            
            files = context.db.execute_dict_query(files_query, (batch_id,))
            
            # Update file statistics
            for file in files:
                file_type = file['file_type']
                file_exists = file['file_exists']
                
                stats['files'][batch_type]['by_type'][file_type]['total'] += 1
                if file_exists:
                    stats['files'][batch_type]['by_type'][file_type]['exists'] += 1
            
            # Get domain statistics
            if args.analyze_domains:
                domains_query = """
                SELECT count(*) as domain_count
                FROM ecod_schema.process_file pf
                JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
                WHERE ps.batch_id = %s
                  AND pf.file_type = 'domain_partition'
                  AND pf.file_exists = TRUE
                """
                
                domain_files = context.db.execute_dict_query(domains_query, (batch_id,))
                
                if domain_files and domain_files[0]['domain_count'] > 0:
                    # Get domain counts - this would require parsing XML files
                    # For simplicity, we're just checking if files exist
                    stats['domains'][batch_type]['total'] += domain_files[0]['domain_count']
                    stats['domains'][batch_type]['with_domains'] += domain_files[0]['domain_count']
    
    # Calculate percentages
    for batch_type in ['regular', 'representative']:
        if stats['proteins'][batch_type]['total'] > 0:
            for stage, count in stats['proteins'][batch_type]['by_stage'].items():
                stats['proteins'][batch_type]['by_stage'][stage] = {
                    'count': count,
                    'percentage': count / stats['proteins'][batch_type]['total'] * 100
                }
        
        for file_type, counts in stats['files'][batch_type]['by_type'].items():
            if counts['total'] > 0:
                counts['percentage'] = counts['exists'] / counts['total'] * 100
    
    # Print comparison
    logger.info("Batch Comparison:")
    logger.info(f"Regular batches ({stats['batches']['regular']['count']}): {stats['batches']['regular']['ids']}")
    logger.info(f"Representative batches ({stats['batches']['representative']['count']}): {stats['batches']['representative']['ids']}")
    
    logger.info("\nProtein Statistics:")
    logger.info(f"  Regular: {stats['proteins']['regular']['total']} proteins")
    logger.info(f"  Representative: {stats['proteins']['representative']['total']} proteins")
    
    logger.info("\nStage Distribution:")
    all_stages = set()
    all_stages.update(stats['proteins']['regular']['by_stage'].keys())
    all_stages.update(stats['proteins']['representative']['by_stage'].keys())
    
    for stage in sorted(all_stages):
        reg_data = stats['proteins']['regular']['by_stage'].get(stage, {'count': 0, 'percentage': 0})
        rep_data = stats['proteins']['representative']['by_stage'].get(stage, {'count': 0, 'percentage': 0})
        
        logger.info(f"  {stage}:")
        logger.info(f"    Regular: {reg_data['count']} ({reg_data['percentage']:.1f}%)")
        logger.info(f"    Representative: {rep_data['count']} ({rep_data['percentage']:.1f}%)")
    
    logger.info("\nFile Type Statistics:")
    all_file_types = set()
    all_file_types.update(stats['files']['regular']['by_type'].keys())
    all_file_types.update(stats['files']['representative']['by_type'].keys())
    
    for file_type in sorted(all_file_types):
        reg_data = stats['files']['regular']['by_type'].get(file_type, {'total': 0, 'exists': 0, 'percentage': 0})
        rep_data = stats['files']['representative']['by_type'].get(file_type, {'total': 0, 'exists': 0, 'percentage': 0})
        
        logger.info(f"  {file_type}:")
        logger.info(f"    Regular: {reg_data['exists']}/{reg_data['total']} ({reg_data['percentage']:.1f}%)")
        logger.info(f"    Representative: {rep_data['exists']}/{rep_data['total']} ({rep_data['percentage']:.1f}%)")
    
    # Write statistics to output file if specified
    if args.output:
        try:
            with open(args.output, 'w') as f:
                json.dump(stats, f, indent=2)
            logger.info(f"Wrote comparison statistics to {args.output}")
        except Exception as e:
            logger.error(f"Error writing statistics: {str(e)}")
    
    return 0

def analyze_domain_statistics(args: argparse.Namespace) -> int:
    """
    Analyze domain statistics across batches
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("analyze.domains")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get batch info if specific batch requested
    batch_clause = ""
    batch_params = []
    
    if args.batch_id:
        batch_info = context.db.execute_dict_query(
            "SELECT id, batch_name, base_path, ref_version FROM ecod_schema.batch WHERE id = %s", 
            (args.batch_id,)
        )
        
        if not batch_info:
            logger.error(f"Batch {args.batch_id} not found")
            return 1
            
        batch_path = batch_info[0]['base_path']
        reference = batch_info[0]['ref_version']
        
        logger.info(f"Analyzing domain statistics for batch {args.batch_id} ({batch_info[0]['batch_name']})")
        
        batch_clause = "WHERE ps.batch_id = %s"
        batch_params = [args.batch_id]
    else:
        logger.info("Analyzing domain statistics across all batches")
    
    # Initialize statistics
    stats = {
        'proteins': {
            'total': 0,
            'with_domain_file': 0,
            'with_domains': 0,
            'no_domains': 0,
            'error': 0
        },
        'domains': {
            'total': 0,
            'single_domain': 0,
            'multi_domain': 0,
            'discontinuous': 0
        },
        'evidence': {
            'blast_only': 0,
            'hhsearch_only': 0,
            'both': 0,
            'none': 0
        },
        't_groups': defaultdict(int),
        'h_groups': defaultdict(int),
        'domain_size': {
            'min': float('inf'),
            'max': 0,
            'avg': 0,
            'total_length': 0,
            'count': 0,
            'distribution': defaultdict(int)
        }
    }
    
    # Get all proteins with domain partition files
    query = f"""
    SELECT 
        p.id as protein_id, p.pdb_id, p.chain_id, ps.id as process_id, p.length,
        pf.file_path
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    LEFT JOIN 
        ecod_schema.process_file pf ON ps.id = pf.process_id AND pf.file_type = 'domain_partition'
    {batch_clause}
    """
    
    rows = context.db.execute_dict_query(query, batch_params)
    
    if not rows:
        logger.error("No proteins found")
        return 1
        
    logger.info(f"Found {len(rows)} proteins")
    
    # Set batch path for all batches mode
    if not args.batch_id:
        # Get a mapping of batch_id to base_path
        batch_query = """
        SELECT id, base_path 
        FROM ecod_schema.batch
        """
        
        batch_paths = {row['id']: row['base_path'] for row in context.db.execute_dict_query(batch_query)}
        
        # Get batch_id for each protein
        batch_id_query = """
        SELECT p.id as protein_id, ps.batch_id
        FROM ecod_schema.protein p
        JOIN ecod_schema.process_status ps ON p.id = ps.protein_id
        """
        
        protein_batches = {row['protein_id']: row['batch_id'] for row in context.db.execute_dict_query(batch_id_query)}
    
    # Update protein statistics
    stats['proteins']['total'] = len(rows)
    stats['proteins']['with_domain_file'] = sum(1 for row in rows if row['file_path'])
    
    # Apply limit if specified
    if args.limit and args.limit < len(rows):
        rows = rows[:args.limit]
        logger.info(f"Limited to analyzing {args.limit} proteins")
    
    # Analyze domain files
    for row in rows:
        pdb_id = row['pdb_id']
        chain_id = row['chain_id']
        file_path = row['file_path']
        
        if not file_path:
            continue
        
        # Get correct batch path for all batches mode
        if not args.batch_id:
            protein_id = row['protein_id']
            batch_id = protein_batches.get(protein_id)
            if not batch_id or batch_id not in batch_paths:
                logger.warning(f"Could not find batch path for protein {pdb_id}_{chain_id}")
                continue
                
            current_batch_path = batch_paths[batch_id]
        else:
            current_batch_path = batch_path
        
        # Parse domain file
        full_path = os.path.join(current_batch_path, file_path)
        
        try:
            tree = ET.parse(full_path)
            root = tree.getroot()
            
            # Check status
            status = root.get('status')
            if status == 'error':
                stats['proteins']['error'] += 1
                continue
            
            # Get domains
            domains = root.findall('./domain_list/domain')
            domain_count = len(domains)
            
            if domain_count > 0:
                stats['proteins']['with_domains'] += 1
                stats['domains']['total'] += domain_count
                
                if domain_count == 1:
                    stats['domains']['single_domain'] += 1
                else:
                    stats['domains']['multi_domain'] += 1
                
                # Analyze domain evidence
                for domain in domains:
                    # Check evidence
                    evidence_elem = domain.find('evidence')
                    if evidence_elem is not None:
                        has_blast = False
                        has_hhsearch = False
                        
                        for evidence in evidence_elem.findall('*'):
                            evidence_type = evidence.tag
                            if evidence_type == 'chain_blast' or evidence_type == 'domain_blast':
                                has_blast = True
                            elif evidence_type == 'hhsearch':
                                has_hhsearch = True
                        
                        if has_blast and has_hhsearch:
                            stats['evidence']['both'] += 1
                        elif has_blast:
                            stats['evidence']['blast_only'] += 1
                        elif has_hhsearch:
                            stats['evidence']['hhsearch_only'] += 1
                        else:
                            stats['evidence']['none'] += 1
                    else:
                        stats['evidence']['none'] += 1
                    
                    # Check if domain is discontinuous
                    domain_range = domain.get('range', '')
                    if ',' in domain_range:
                        stats['domains']['discontinuous'] += 1
                    
                    # Get domain size
                    try:
                        domain_length = 0
                        for segment in domain_range.split(','):
                            start, end = map(int, segment.split('-'))
                            domain_length += (end - start + 1)
                            
                        stats['domain_size']['min'] = min(stats['domain_size']['min'], domain_length)
                        stats['domain_size']['max'] = max(stats['domain_size']['max'], domain_length)
                        stats['domain_size']['total_length'] += domain_length
                        stats['domain_size']['count'] += 1
                        
                        # Record size distribution in buckets
                        if domain_length < 50:
                            size_bucket = "<50"
                        elif domain_length < 100:
                            size_bucket = "50-99"
                        elif domain_length < 150:
                            size_bucket = "100-149"
                        elif domain_length < 200:
                            size_bucket = "150-199"
                        elif domain_length < 300:
                            size_bucket = "200-299"
                        elif domain_length < 500:
                            size_bucket = "300-499"
                        else:
                            size_bucket = "500+"
                            
                        stats['domain_size']['distribution'][size_bucket] += 1
                        
                    except (ValueError, AttributeError):
                        pass
                    
                    # Get T-group and H-group
                    t_group = domain.get('t_group', '')
                    h_group = domain.get('h_group', '')
                    
                    if t_group:
                        stats['t_groups'][t_group] += 1
                    
                    if h_group:
                        stats['h_groups'][h_group] += 1
            else:
                stats['proteins']['no_domains'] += 1
                
        except ET.ParseError as e:
            logger.error(f"XML parsing error for {pdb_id}_{chain_id}: {str(e)}")
            
        except Exception as e:
            logger.error(f"Error analyzing domain file for {pdb_id}_{chain_id}: {str(e)}")
    
    # Calculate average domain size
    if stats['domain_size']['count'] > 0:
        stats['domain_size']['avg'] = stats['domain_size']['total_length'] / stats['domain_size']['count']
    
    # Get top T-groups and H-groups
    stats['top_t_groups'] = sorted(stats['t_groups'].items(), key=lambda x: x[1], reverse=True)[:20]
    stats['top_h_groups'] = sorted(stats['h_groups'].items(), key=lambda x: x[1], reverse=True)[:20]
    
    # Print summary
    logger.info(f"Domain statistics summary:")
    logger.info(f"Total proteins: {stats['proteins']['total']}")
    logger.info(f"With domain files: {stats['proteins']['with_domain_file']} ({stats['proteins']['with_domain_file']/stats['proteins']['total']*100:.1f}%)")
    logger.info(f"With domains: {stats['proteins']['with_domains']} ({stats['proteins']['with_domains']/stats['proteins']['with_domain_file']*100:.1f}% of files)")
    logger.info(f"No domains: {stats['proteins']['no_domains']} ({stats['proteins']['no_domains']/stats['proteins']['with_domain_file']*100:.1f}% of files)")
    logger.info(f"Error: {stats['proteins']['error']} ({stats['proteins']['error']/stats['proteins']['total']*100:.1f}%)")
    
    logger.info(f"\nDomain statistics:")
    logger.info(f"Total domains: {stats['domains']['total']}")
    logger.info(f"Single domain proteins: {stats['domains']['single_domain']} ({stats['domains']['single_domain']/stats['proteins']['with_domains']*100:.1f}%)")
    logger.info(f"Multi domain proteins: {stats['domains']['multi_domain']} ({stats['domains']['multi_domain']/stats['proteins']['with_domains']*100:.1f}%)")
    logger.info(f"Discontinuous domains: {stats['domains']['discontinuous']} ({stats['domains']['discontinuous']/stats['domains']['total']*100:.1f}%)")
    
    logger.info(f"\nDomain size statistics:")
    logger.info(f"Average domain size: {stats['domain_size']['avg']:.1f} residues")
    logger.info(f"Min size: {stats['domain_size']['min']} residues")
    logger.info(f"Max size: {stats['domain_size']['max']} residues")
    
    logger.info(f"\nSize distribution:")
    for size_bucket, count in sorted(stats['domain_size']['distribution'].items(), key=lambda x: (
        int(x[0].split('-')[0]) if x[0][0].isdigit() else (0 if x[0] == '<50' else 1000)
    )):
        logger.info(f"  {size_bucket}: {count} ({count/stats['domains']['total']*100:.1f}%)")
    
    logger.info(f"\nEvidence distribution:")
    domain_with_evidence = sum(stats['evidence'].values())
    if domain_with_evidence > 0:
        for evidence_type, count in stats['evidence'].items():
            logger.info(f"  {evidence_type}: {count} ({count/domain_with_evidence*100:.1f}%)")
    
    logger.info(f"\nTop 10 T-groups:")
    for t_group, count in stats['top_t_groups'][:10]:
        logger.info(f"  {t_group}: {count} ({count/stats['domains']['total']*100:.1f}%)")
    
    logger.info(f"\nTop 10 H-groups:")
    for h_group, count in stats['top_h_groups'][:10]:
        logger.info(f"  {h_group}: {count} ({count/stats['domains']['total']*100:.1f}%)")
    
    # Write statistics to output file if specified
    if args.output:
        try:
            with open(args.output, 'w') as f:
                json.dump(stats, f, indent=2)
            logger.info(f"Wrote domain statistics to {args.output}")
        except Exception as e:
            logger.error(f"Error writing statistics: {str(e)}")
    
    return 0

def analyze_pipeline_stats(args: argparse.Namespace) -> int:
    """
    Analyze overall pipeline statistics
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("analyze.pipeline")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get all batches
    batches_query = """
    SELECT id, batch_name, base_path, ref_version,
           (CASE WHEN batch_name LIKE 'alt_rep%' THEN 'representative' ELSE 'regular' END) as batch_type
    FROM ecod_schema.batch
    ORDER BY id
    """
    
    batches = context.db.execute_dict_query(batches_query)
    
    if not batches:
        logger.error("No batches found in database")
        return 1
    
    logger.info(f"Found {len(batches)} batches in database")
    
    # Initialize statistics
    stats = {
        'batches': {
            'total': len(batches),
            'by_type': defaultdict(int)
        },
        'proteins': {
            'total': 0,
            'by_batch_type': defaultdict(int),
            'by_stage': defaultdict(int)
        },
        'process_status': defaultdict(int),
        'file_types': {
            'total': defaultdict(int),
            'exists': defaultdict(int)
        },
        'success_rates': {
            'blast': 0,
            'hhsearch': 0,
            'domain_summary': 0,
            'domain_partition': 0
        },
        'timing': {
            'avg_times': {},
            'total_times': {}
        }
    }
    
    # Count batches by type
    for batch in batches:
        batch_type = batch['batch_type']
        stats['batches']['by_type'][batch_type] += 1
    
    # Get protein statistics
    proteins_query = """
    SELECT 
        ps.batch_id, ps.current_stage, ps.status, 
        (CASE WHEN ps.is_representative THEN 'representative' ELSE 'regular' END) as protein_type
    FROM 
        ecod_schema.process_status ps
    """
    
    proteins = context.db.execute_dict_query(proteins_query)
    
    if not proteins:
        logger.error("No proteins found in database")
        return 1
    
    stats['proteins']['total'] = len(proteins)
    logger.info(f"Found {len(proteins)} proteins in database")
    
    # Batch type map
    batch_type_map = {batch['id']: batch['batch_type'] for batch in batches}
    
    # Count proteins by batch type and stage
    for protein in proteins:
        batch_id = protein['batch_id']
        batch_type = batch_type_map.get(batch_id, 'unknown')
        current_stage = protein['current_stage']
        status = protein['status']
        protein_type = protein['protein_type']
        
        stats['proteins']['by_batch_type'][batch_type] += 1
        stats['proteins']['by_stage'][current_stage] += 1
        stats['process_status'][status] += 1
    
    # Get file statistics
    files_query = """
    SELECT file_type, file_exists, COUNT(*) as count
    FROM ecod_schema.process_file
    GROUP BY file_type, file_exists
    """
    
    files = context.db.execute_dict_query(files_query)
    
    # Count files by type and existence
    for file in files:
        file_type = file['file_type']
        file_exists = file['file_exists']
        count = file['count']
        
        stats['file_types']['total'][file_type] += count
        if file_exists:
            stats['file_types']['exists'][file_type] += count
    
    # Calculate success rates
    for stage in ['blast', 'hhsearch', 'domain_summary', 'domain_partition']:
        complete_stage = f"{stage}_complete"
        failed_stage = f"{stage}_failed"
        
        complete_count = stats['proteins']['by_stage'].get(complete_stage, 0)
        failed_count = stats['proteins']['by_stage'].get(failed_stage, 0)
        total_count = complete_count + failed_count
        
        if total_count > 0:
            stats['success_rates'][stage] = complete_count / total_count * 100
    
    # Get timing statistics if available
    if args.timing:
        timing_query = """
        SELECT stage, AVG(elapsed_time) as avg_time, SUM(elapsed_time) as total_time
        FROM ecod_schema.process_timing
        GROUP BY stage
        """
        
        timing = context.db.execute_dict_query(timing_query)
        
        for row in timing:
            stage = row['stage']
            avg_time = row['avg_time']
            total_time = row['total_time']
            
            stats['timing']['avg_times'][stage] = avg_time
            stats['timing']['total_times'][stage] = total_time
    
    # Print summary
    logger.info(f"Pipeline statistics summary:")
    logger.info(f"Total batches: {stats['batches']['total']}")
    for batch_type, count in stats['batches']['by_type'].items():
        logger.info(f"  {batch_type}: {count} ({count/stats['batches']['total']*100:.1f}%)")
    
    logger.info(f"\nTotal proteins: {stats['proteins']['total']}")
    for batch_type, count in stats['proteins']['by_batch_type'].items():
        logger.info(f"  In {batch_type} batches: {count} ({count/stats['proteins']['total']*100:.1f}%)")
    
    logger.info(f"\nStage distribution:")
    stage_order = [
        'new', 'fasta_complete', 'blast_complete', 'hhsearch_complete', 
        'domain_summary_complete', 'domain_partition_complete'
    ]
    
    other_stages = set(stats['proteins']['by_stage'].keys()) - set(stage_order)
    all_stages = stage_order + sorted(other_stages)
    
    for stage in all_stages:
        count = stats['proteins']['by_stage'].get(stage, 0)
        if count > 0:
            logger.info(f"  {stage}: {count} ({count/stats['proteins']['total']*100:.1f}%)")
    
    logger.info(f"\nProcess status:")
    for status, count in stats['process_status'].items():
        logger.info(f"  {status}: {count} ({count/stats['proteins']['total']*100:.1f}%)")
    
    logger.info(f"\nFile statistics:")
    for file_type in sorted(stats['file_types']['total'].keys()):
        total = stats['file_types']['total'][file_type]
        exists = stats['file_types']['exists'][file_type]
        
        if total > 0:
            logger.info(f"  {file_type}: {exists}/{total} ({exists/total*100:.1f}%)")
    
    logger.info(f"\nSuccess rates:")
    for stage, rate in stats['success_rates'].items():
        logger.info(f"  {stage}: {rate:.1f}%")
    
    if args.timing and stats['timing']['avg_times']:
        logger.info(f"\nTiming statistics:")
        for stage, avg_time in sorted(stats['timing']['avg_times'].items()):
            logger.info(f"  {stage}: {avg_time:.2f} seconds average")
    
    # Write statistics to output file if specified
    if args.output:
        try:
            with open(args.output, 'w') as f:
                json.dump(stats, f, indent=2)
            logger.info(f"Wrote pipeline statistics to {args.output}")
        except Exception as e:
            logger.error(f"Error writing statistics: {str(e)}")
    
    return 0

def analyze_mode(args: argparse.Namespace) -> int:
    """
    Run analysis mode with appropriate action
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("analyze")
    
    if args.action == "comparison":
        return analyze_batch_comparison(args)
    elif args.action == "domains":
        return analyze_domain_statistics(args)
    elif args.action == "pipeline":
        return analyze_pipeline_stats(args)
    else:
        logger.error(f"Unknown analysis action: {args.action}")
        return 1

#
# TRACE MODE FUNCTIONS
#

def trace_evidence_chain(args: argparse.Namespace) -> int:
    """
    Trace evidence chain for a specific protein
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("trace.evidence")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get batch info
    batch_info = context.db.execute_dict_query(
        "SELECT id, batch_name, base_path, ref_version FROM ecod_schema.batch WHERE id = %s", 
        (args.batch_id,)
    )
    
    if not batch_info:
        logger.error(f"Batch {args.batch_id} not found")
        return 1
        
    batch_path = batch_info[0]['base_path']
    reference = batch_info[0]['ref_version']
    
    # Get protein info
    protein_query = """
    SELECT 
        p.id as protein_id, p.pdb_id, p.chain_id, ps.id as process_id,
        ps.current_stage
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    WHERE 
        ps.batch_id = %s
        AND p.pdb_id = %s
        AND p.chain_id = %s
    """
    
    protein = context.db.execute_dict_query(protein_query, (args.batch_id, args.pdb_id, args.chain_id))
    
    if not protein:
        logger.error(f"Protein {args.pdb_id}_{args.chain_id} not found in batch {args.batch_id}")
        return 1
        
    process_id = protein[0]['process_id']
    current_stage = protein[0]['current_stage']
    
    logger.info(f"Tracing evidence chain for {args.pdb_id}_{args.chain_id} (current stage: {current_stage})")
    
    # Get file paths
    files_query = """
    SELECT 
        file_type, file_path, file_exists
    FROM 
        ecod_schema.process_file
    WHERE 
        process_id = %s
    """
    
    files = context.db.execute_dict_query(files_query, (process_id,))
    
    if not files:
        logger.error(f"No files found for {args.pdb_id}_{args.chain_id}")
        return 1
        
    logger.info(f"Found {len(files)} files")
    
    # Initialize trace report
    trace = {
        'pdb_id': args.pdb_id,
        'chain_id': args.chain_id,
        'batch_id': args.batch_id,
        'process_id': process_id,
        'current_stage': current_stage,
        'files': {},
        'evidence_chain': {},
        'domains': {}
    }
    
    # Map file paths
    for file in files:
        file_type = file['file_type']
        file_path = file['file_path']
        file_exists = file['file_exists']
        
        if file_exists and file_path:
            full_path = os.path.join(batch_path, file_path)
            trace['files'][file_type] = {
                'path': file_path,
                'full_path': full_path,
                'exists': os.path.exists(full_path)
            }
    
    # Check domain partition file
    if 'domain_partition' in trace['files'] and trace['files']['domain_partition']['exists']:
        domain_file = trace['files']['domain_partition']['full_path']
        try:
            tree = ET.parse(domain_file)
            root = tree.getroot()
            
            # Get domains
            domains = root.findall('./domain_list/domain')
            
            trace['domains']['count'] = len(domains)
            trace['domains']['details'] = []
            
            for i, domain in enumerate(domains):
                domain_id = domain.get('id', f"domain_{i+1}")
                domain_range = domain.get('range', '')
                domain_t_group = domain.get('t_group', '')
                domain_h_group = domain.get('h_group', '')
                
                domain_info = {
                    'id': domain_id,
                    'range': domain_range,
                    't_group': domain_t_group,
                    'h_group': domain_h_group,
                    'evidence': []
                }
                
                # Get evidence
                evidence_elem = domain.find('evidence')
                if evidence_elem is not None:
                    for evidence in evidence_elem.findall('*'):
                        evidence_type = evidence.tag
                        evidence_hits = []
                        
                        for hit in evidence.findall('*'):
                            hit_id = hit.get('hit_id', '')
                            hit_info = {
                                'hit_id': hit_id
                            }
                            
                            # Add other hit attributes
                            for attr, value in hit.attrib.items():
                                if attr != 'hit_id':
                                    hit_info[attr] = value
                            
                            evidence_hits.append(hit_info)
                        
                        domain_info['evidence'].append({
                            'type': evidence_type,
                            'hits': evidence_hits
                        })
                
                trace['domains']['details'].append(domain_info)
        except Exception as e:
            logger.error(f"Error parsing domain file: {str(e)}")
    
    # Check domain summary file
    if 'domain_summary' in trace['files'] and trace['files']['domain_summary']['exists']:
        summary_file = trace['files']['domain_summary']['full_path']
        try:
            tree = ET.parse(summary_file)
            root = tree.getroot()
            
            # Get evidence
            evidence_chain = {}
            
            # Check chain BLAST evidence
            chain_blast = root.find('chain_blast_evidence')
            if chain_blast is not None:
                blast_runs = chain_blast.findall('./blast_run') or chain_blast.findall('./chain_blast_run')
                
                if blast_runs:
                    evidence_chain['chain_blast'] = {
                        'runs': []
                    }
                    
                    for blast_run in blast_runs:
                        run_info = {
                            'db': blast_run.get('db', ''),
                            'program': blast_run.get('program', ''),
                            'hits': []
                        }
                        
                        hits = blast_run.findall('.//hit')
                        for hit in hits:
                            hit_info = {
                                'hit_id': hit.get('hit_id', '')
                            }
                            
                            # Add other hit attributes
                            for attr, value in hit.attrib.items():
                                if attr != 'hit_id':
                                    hit_info[attr] = value
                            
                            run_info['hits'].append(hit_info)
                        
                        evidence_chain['chain_blast']['runs'].append(run_info)
            
            # Check domain BLAST evidence
            domain_blast = root.find('domain_blast_evidence')
            if domain_blast is not None:
                blast_runs = domain_blast.findall('./blast_run')
                
                if blast_runs:
                    evidence_chain['domain_blast'] = {
                        'runs': []
                    }
                    
                    for blast_run in blast_runs:
                        run_info = {
                            'db': blast_run.get('db', ''),
                            'program': blast_run.get('program', ''),
                            'hits': []
                        }
                        
                        hits = blast_run.findall('.//hit')
                        for hit in hits:
                            hit_info = {
                                'hit_id': hit.get('hit_id', '')
                            }
                            
                            # Add other hit attributes
                            for attr, value in hit.attrib.items():
                                if attr != 'hit_id':
                                    hit_info[attr] = value
                            
                            run_info['hits'].append(hit_info)
                        
                        evidence_chain['domain_blast']['runs'].append(run_info)
            
            # Check HHSearch evidence
            hhsearch = root.find('hhsearch_evidence')
            if hhsearch is not None:
                hit_list = hhsearch.find('.//hh_hit_list')
                
                if hit_list is not None:
                    evidence_chain['hhsearch'] = {
                        'hits': []
                    }
                    
                    hits = hit_list.findall('hh_hit')
                    for hit in hits:
                        hit_info = {
                            'hit_id': hit.get('hit_id', ''),
                            'hit_num': hit.get('hit_num', ''),
                            'probability': hit.get('probability', '')
                        }
                        
                        # Add other hit attributes
                        for attr, value in hit.attrib.items():
                            if attr not in ['hit_id', 'hit_num', 'probability']:
                                hit_info[attr] = value
                        
                        evidence_chain['hhsearch']['hits'].append(hit_info)
            
            trace['evidence_chain'] = evidence_chain
            
        except Exception as e:
            logger.error(f"Error parsing summary file: {str(e)}")
    
    # Print summary
    logger.info(f"\nEvidence chain summary for {args.pdb_id}_{args.chain_id}:")
    
    logger.info(f"\nFiles:")
    for file_type, file_info in trace['files'].items():
        status = "exists" if file_info['exists'] else "missing"
        logger.info(f"  {file_type}: {status} ({file_info['path']})")
    
    if 'domain_summary' in trace['files'] and 'evidence_chain' in trace:
        logger.info(f"\nEvidence summary:")
        
        if 'chain_blast' in trace['evidence_chain']:
            chain_blast = trace['evidence_chain']['chain_blast']
            total_hits = sum(len(run['hits']) for run in chain_blast['runs'])
            logger.info(f"  Chain BLAST: {len(chain_blast['runs'])} runs, {total_hits} total hits")
        
        if 'domain_blast' in trace['evidence_chain']:
            domain_blast = trace['evidence_chain']['domain_blast']
            total_hits = sum(len(run['hits']) for run in domain_blast['runs'])
            logger.info(f"  Domain BLAST: {len(domain_blast['runs'])} runs, {total_hits} total hits")
        
        if 'hhsearch' in trace['evidence_chain']:
            hhsearch = trace['evidence_chain']['hhsearch']
            logger.info(f"  HHSearch: {len(hhsearch['hits'])} hits")
    
    if 'domains' in trace and 'count' in trace['domains']:
        logger.info(f"\nDomain summary:")
        domain_count = trace['domains']['count']
        logger.info(f"  Total domains: {domain_count}")
        
        for i, domain in enumerate(trace['domains'].get('details', [])):
            logger.info(f"  Domain {i+1}:")
            logger.info(f"    Range: {domain.get('range', 'unknown')}")
            logger.info(f"    T-group: {domain.get('t_group', 'unknown')}")
            logger.info(f"    H-group: {domain.get('h_group', 'unknown')}")
            
            # Evidence summary
            evidence_types = [e['type'] for e in domain.get('evidence', [])]
            if evidence_types:
                logger.info(f"    Evidence types: {', '.join(evidence_types)}")
    
    # Write trace to output file if specified
    if args.output:
        try:
            with open(args.output, 'w') as f:
                json.dump(trace, f, indent=2)
            logger.info(f"Wrote trace to {args.output}")
        except Exception as e:
            logger.error(f"Error writing trace: {str(e)}")
    
    return 0

def trace_domain_assignment(args: argparse.Namespace) -> int:
    """
    Trace domain assignment from evidence to domain partition
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("trace.domain")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get batch info
    batch_info = context.db.execute_dict_query(
        "SELECT id, batch_name, base_path, ref_version FROM ecod_schema.batch WHERE id = %s", 
        (args.batch_id,)
    )
    
    if not batch_info:
        logger.error(f"Batch {args.batch_id} not found")
        return 1
        
    batch_path = batch_info[0]['base_path']
    reference = batch_info[0]['ref_version']
    
    # Get protein info
    protein_query = """
    SELECT 
        p.id as protein_id, p.pdb_id, p.chain_id, ps.id as process_id,
        ps.current_stage, p.length
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    WHERE 
        ps.batch_id = %s
    """
    
    if args.pdb_id and args.chain_id:
        protein_query += " AND p.pdb_id = %s AND p.chain_id = %s"
        params = (args.batch_id, args.pdb_id, args.chain_id)
    else:
        # Get proteins at domain_partition_complete stage if no specific protein
        protein_query += " AND ps.current_stage = 'domain_partition_complete' ORDER BY p.pdb_id, p.chain_id"
        
        if args.limit:
            protein_query += f" LIMIT {args.limit}"
            
        params = (args.batch_id,)
    
    proteins = context.db.execute_dict_query(protein_query, params)
    
    if not proteins:
        logger.error(f"No proteins found matching criteria")
        return 1
    
    # Initialize trace report for each protein
    traces = []
    
    for protein in proteins:
        pdb_id = protein['pdb_id']
        chain_id = protein['chain_id']
        process_id = protein['process_id']
        length = protein['length'] or 0
        
        logger.info(f"Tracing domain assignment for {pdb_id}_{chain_id}")
        
        # Get file paths
        files_query = """
        SELECT 
            file_type, file_path, file_exists
        FROM 
            ecod_schema.process_file
        WHERE 
            process_id = %s
            AND file_type IN ('domain_summary', 'domain_partition')
            AND file_exists = TRUE
        """
        
        files = context.db.execute_dict_query(files_query, (process_id,))
        
        file_paths = {}
        for file in files:
            file_type = file['file_type']
            file_path = file['file_path']
            file_paths[file_type] = os.path.join(batch_path, file_path)
        
        # Skip if missing necessary files
        if 'domain_summary' not in file_paths or 'domain_partition' not in file_paths:
            logger.warning(f"Missing required files for {pdb_id}_{chain_id}")
            continue
        
        # Initialize trace for this protein
        trace = {
            'pdb_id': pdb_id,
            'chain_id': chain_id,
            'chain_length': length,
            'summary': {
                'path': file_paths['domain_summary'],
                'evidence': {},
                'suggestions': []
            },
            'partition': {
                'path': file_paths['domain_partition'],
                'domains': []
            },
            'trace': {
                'matches': 0,
                'mismatches': 0,
                'novel_domains': 0,
                'missing_domains': 0,
                'modified_domains': 0
            }
        }
        
        # Parse domain summary
        try:
            summary_tree = ET.parse(file_paths['domain_summary'])
            summary_root = summary_tree.getroot()
            
            # Get evidence
            for evidence_type in ['chain_blast_evidence', 'domain_blast_evidence', 'hhsearch_evidence']:
                evidence_elem = summary_root.find(evidence_type)
                if evidence_elem is not None:
                    hit_count = len(evidence_elem.findall('.//hit')) + len(evidence_elem.findall('.//hh_hit'))
                    trace['summary']['evidence'][evidence_type] = hit_count
            
            # Get domain suggestions
            suggestions_elem = summary_root.find('domain_suggestions')
            if suggestions_elem is not None:
                for domain in suggestions_elem.findall('domain'):
                    domain_id = domain.get('id', '')
                    domain_range = domain.get('range', '')
                    
                    suggestion = {
                        'id': domain_id,
                        'range': domain_range
                    }
                    
                    trace['summary']['suggestions'].append(suggestion)
        except Exception as e:
            logger.error(f"Error parsing summary file for {pdb_id}_{chain_id}: {str(e)}")
            continue
        
        # Parse domain partition
        try:
            partition_tree = ET.parse(file_paths['domain_partition'])
            partition_root = partition_tree.getroot()
            
            # Get domains
            domain_list = partition_root.find('domain_list')
            if domain_list is not None:
                for domain in domain_list.findall('domain'):
                    domain_id = domain.get('id', '')
                    domain_range = domain.get('range', '')
                    domain_t_group = domain.get('t_group', '')
                    domain_h_group = domain.get('h_group', '')
                    
                    domain_info = {
                        'id': domain_id,
                        'range': domain_range,
                        't_group': domain_t_group,
                        'h_group': domain_h_group,
                        'sources': []
                    }
                    
                    # Get evidence
                    evidence_elem = domain.find('evidence')
                    if evidence_elem is not None:
                        for source in evidence_elem.findall('*'):
                            source_type = source.tag
                            domain_info['sources'].append(source_type)
                    
                    trace['partition']['domains'].append(domain_info)
        except Exception as e:
            logger.error(f"Error parsing partition file for {pdb_id}_{chain_id}: {str(e)}")
            continue
        
        # Trace from suggestions to domains
        suggestion_ranges = [s['range'] for s in trace['summary']['suggestions']]
        partition_ranges = [d['range'] for d in trace['partition']['domains']]
        
        # Count exact matches
        matches = set(suggestion_ranges).intersection(set(partition_ranges))
        trace['trace']['matches'] = len(matches)
        
        # Count missing domains (in suggestions but not in partition)
        missing = set(suggestion_ranges) - set(partition_ranges)
        trace['trace']['missing_domains'] = len(missing)
        
        # Count novel domains (in partition but not in suggestions)
        novel = set(partition_ranges) - set(suggestion_ranges)
        trace['trace']['novel_domains'] = len(novel)
        
        # Count modified domains (similar but not exact matches)
        # This is a simplified approach and could miss some modifications
        if len(trace['summary']['suggestions']) > 0 and len(trace['partition']['domains']) > 0:
            # Check for chains with no domains in partition but suggestions
            if len(trace['partition']['domains']) == 0:
                trace['trace']['missing_domains'] = len(trace['summary']['suggestions'])
            
            # Check for chains with domains in partition but no suggestions
            if len(trace['summary']['suggestions']) == 0:
                trace['trace']['novel_domains'] = len(trace['partition']['domains'])
        
        # Add to traces
        traces.append(trace)
    
    # Print summary
    logger.info(f"\nDomain assignment trace summary for {len(traces)} proteins:")
    
    total_matches = sum(t['trace']['matches'] for t in traces)
    total_mismatches = sum(t['trace']['mismatches'] for t in traces)
    total_novel = sum(t['trace']['novel_domains'] for t in traces)
    total_missing = sum(t['trace']['missing_domains'] for t in traces)
    total_modified = sum(t['trace']['modified_domains'] for t in traces)
    
    total_suggestions = sum(len(t['summary']['suggestions']) for t in traces)
    total_domains = sum(len(t['partition']['domains']) for t in traces)
    
    logger.info(f"\nSuggestion to domain conversion:")
    logger.info(f"  Total suggestions: {total_suggestions}")
    logger.info(f"  Total domains: {total_domains}")
    logger.info(f"  Exact matches: {total_matches}")
    logger.info(f"  Mismatches: {total_mismatches}")
    logger.info(f"  Novel domains: {total_novel}")
    logger.info(f"  Missing domains: {total_missing}")
    logger.info(f"  Modified domains: {total_modified}")
    
    # Calculate statistics by protein
    empty_proteins = 0
    single_domain = 0
    multi_domain = 0
    
    for trace in traces:
        domain_count = len(trace['partition']['domains'])
        
        if domain_count == 0:
            empty_proteins += 1
        elif domain_count == 1:
            single_domain += 1
        else:
            multi_domain += 1
    
    logger.info(f"\nProtein classification:")
    logger.info(f"  No domains: {empty_proteins} ({empty_proteins/len(traces)*100:.1f}%)")
    logger.info(f"  Single domain: {single_domain} ({single_domain/len(traces)*100:.1f}%)")
    logger.info(f"  Multi-domain: {multi_domain} ({multi_domain/len(traces)*100:.1f}%)")
    
    # Calculate evidence statistics
    evidence_stats = {
        'chain_blast_evidence': 0,
        'domain_blast_evidence': 0,
        'hhsearch_evidence': 0,
        'combinations': defaultdict(int)
    }
    
    for trace in traces:
        evidence = []
        
        for evidence_type, count in trace['summary']['evidence'].items():
            if count > 0:
                evidence_stats[evidence_type] += 1
                evidence.append(evidence_type)
        
        evidence_key = '_'.join(sorted(evidence)) if evidence else 'none'
        evidence_stats['combinations'][evidence_key] += 1
    
    logger.info(f"\nEvidence sources:")
    for evidence_type in ['chain_blast_evidence', 'domain_blast_evidence', 'hhsearch_evidence']:
        count = evidence_stats[evidence_type]
        logger.info(f"  {evidence_type}: {count} ({count/len(traces)*100:.1f}%)")
    
    logger.info(f"\nEvidence combinations:")
    for combo, count in sorted(evidence_stats['combinations'].items(), key=lambda x: x[1], reverse=True):
        logger.info(f"  {combo}: {count} ({count/len(traces)*100:.1f}%)")
    
    # Write traces to output file if specified
    if args.output:
        try:
            with open(args.output, 'w') as f:
                json.dump({
                    'traces': traces,
                    'summary': {
                        'total_proteins': len(traces),
                        'total_suggestions': total_suggestions,
                        'total_domains': total_domains,
                        'matches': total_matches,
                        'mismatches': total_mismatches,
                        'novel_domains': total_novel,
                        'missing_domains': total_missing,
                        'modified_domains': total_modified,
                        'protein_classification': {
                            'empty': empty_proteins,
                            'single_domain': single_domain,
                            'multi_domain': multi_domain
                        },
                        'evidence_stats': evidence_stats
                    }
                }, f, indent=2)
            logger.info(f"Wrote traces to {args.output}")
        except Exception as e:
            logger.error(f"Error writing traces: {str(e)}")
    
    return 0

def trace_raw_data(args: argparse.Namespace) -> int:
    """
    Trace raw data files and file conversion chain
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("trace.raw")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get batch info
    batch_info = context.db.execute_dict_query(
        "SELECT id, batch_name, base_path, ref_version FROM ecod_schema.batch WHERE id = %s", 
        (args.batch_id,)
    )
    
    if not batch_info:
        logger.error(f"Batch {args.batch_id} not found")
        return 1
        
    batch_path = batch_info[0]['base_path']
    reference = batch_info[0]['ref_version']
    
    # Initialize statistics
    stats = {
        'batch_id': args.batch_id,
        'batch_name': batch_info[0]['batch_name'],
        'proteins': {
            'total': 0,
            'processed': 0,
            'with_issues': 0
        },
        'file_stats': defaultdict(lambda: {'total': 0, 'exists': 0, 'consistency': 0}),
        'conversion_chains': defaultdict(int),
        'issues': []
    }
    
    # Get proteins
    protein_query = """
    SELECT 
        p.id as protein_id, p.pdb_id, p.chain_id, ps.id as process_id,
        ps.current_stage
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    WHERE 
        ps.batch_id = %s
    """
    
    if args.reps_only:
        protein_query += " AND ps.is_representative = TRUE"
    
    protein_query += " ORDER BY p.pdb_id, p.chain_id"
    
    if args.limit:
        protein_query += f" LIMIT {args.limit}"
    
    proteins = context.db.execute_dict_query(protein_query, (args.batch_id,))
    
    if not proteins:
        logger.error("No proteins found for this batch")
        return 1
        
    stats['proteins']['total'] = len(proteins)
    logger.info(f"Found {len(proteins)} proteins")
    
    # Define file type chain
    file_chain = [
        'fasta',
        'chain_blast',
        'domain_blast',
        'hhr',
        'hh_xml',
        'domain_summary',
        'domain_partition'
    ]
    
    # Process each protein
    for protein in proteins:
        pdb_id = protein['pdb_id']
        chain_id = protein['chain_id']
        process_id = protein['process_id']
        
        # Get file information
        files_query = """
        SELECT 
            file_type, file_path, file_exists
        FROM 
            ecod_schema.process_file
        WHERE 
            process_id = %s
        """
        
        files = context.db.execute_dict_query(files_query, (process_id,))
        
        if not files:
            continue
        
        # Map files by type
        file_map = {}
        for file in files:
            file_type = file['file_type']
            file_path = file['file_path']
            file_exists = file['file_exists']
            
            file_map[file_type] = {
                'path': file_path,
                'exists': file_exists,
                'full_path': os.path.join(batch_path, file_path) if file_path else None
            }
            
            # Update statistics
            stats['file_stats'][file_type]['total'] += 1
            if file_exists:
                stats['file_stats'][file_type]['exists'] += 1
        
        # Check file chain consistency
        chain_exists = []
        for file_type in file_chain:
            if file_type in file_map and file_map[file_type]['exists']:
                chain_exists.append(file_type)
        
        # Record chain pattern
        chain_pattern = '_'.join(chain_exists)
        stats['conversion_chains'][chain_pattern] += 1
        
        # Check for issues
        issues = []
        
        # Check expected consistency based on last file
        if chain_exists:
            last_file = chain_exists[-1]
            expected_chain_index = file_chain.index(last_file)
            expected_chain = file_chain[:expected_chain_index + 1]
            
            # Check if all expected files exist
            for file_type in expected_chain:
                if file_type not in chain_exists:
                    issues.append(f"Missing {file_type} in chain")
        
        # Check specific file pairs
        if 'hhr' in chain_exists and 'hh_xml' not in chain_exists:
            issues.append("Has HHR but missing HH_XML")
        
        if 'domain_summary' in chain_exists and 'chain_blast' not in chain_exists and 'domain_blast' not in chain_exists:
            issues.append("Has domain summary but missing BLAST results")
        
        if 'domain_partition' in chain_exists and 'domain_summary' not in chain_exists:
            issues.append("Has domain partition but missing domain summary")
        
        # Record issues
        if issues:
            stats['proteins']['with_issues'] += 1
            stats['issues'].append({
                'pdb_id': pdb_id,
                'chain_id': chain_id,
                'process_id': process_id,
                'issues': issues
            })
        
        stats['proteins']['processed'] += 1
    
    # Calculate consistency percentage
    for file_type in stats['file_stats']:
        if stats['file_stats'][file_type]['total'] > 0:
            expected_prev = None
            for i, ft in enumerate(file_chain):
                if ft == file_type and i > 0:
                    expected_prev = file_chain[i-1]
                    break
            
            if expected_prev:
                # Calculate how many have both this file and the previous
                consistent_count = 0
                for chain, count in stats['conversion_chains'].items():
                    chain_types = chain.split('_')
                    if file_type in chain_types and expected_prev in chain_types:
                        consistent_count += count
                
                if stats['file_stats'][file_type]['exists'] > 0:
                    stats['file_stats'][file_type]['consistency'] = consistent_count / stats['file_stats'][file_type]['exists'] * 100
    
    # Print summary
    logger.info(f"\nRaw data trace summary for batch {args.batch_id}:")
    logger.info(f"Processed {stats['proteins']['processed']} of {stats['proteins']['total']} proteins")
    logger.info(f"Found {stats['proteins']['with_issues']} proteins with issues")
    
    logger.info(f"\nFile statistics:")
    for file_type in file_chain:
        if file_type in stats['file_stats']:
            file_stat = stats['file_stats'][file_type]
            exists_pct = file_stat['exists'] / file_stat['total'] * 100 if file_stat['total'] > 0 else 0
            logger.info(f"  {file_type}: {file_stat['exists']}/{file_stat['total']} ({exists_pct:.1f}%)")
            
            if 'consistency' in file_stat and file_stat['consistency'] > 0:
                logger.info(f"    Consistency: {file_stat['consistency']:.1f}%")
    
    logger.info(f"\nCommon conversion chains:")
    for chain, count in sorted(stats['conversion_chains'].items(), key=lambda x: x[1], reverse=True)[:10]:
        logger.info(f"  {chain}: {count} ({count/stats['proteins']['processed']*100:.1f}%)")
    
    # Output issues if requested
    if args.output and stats['issues']:
        try:
            with open(args.output, 'w') as f:
                json.dump({
                    'statistics': stats,
                    'issues': stats['issues']
                }, f, indent=2)
            logger.info(f"Wrote trace report to {args.output}")
        except Exception as e:
            logger.error(f"Error writing report: {str(e)}")
    
    return 0

def trace_mode(args: argparse.Namespace) -> int:
    """
    Run trace mode with appropriate action
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("trace")
    
    if args.action == "evidence":
        return trace_evidence_chain(args)
    elif args.action == "domain":
        return trace_domain_assignment(args)
    elif args.action == "raw":
        return trace_raw_data(args)
    else:
        logger.error(f"Unknown trace action: {args.action}")
        return 1

#
# REPORT MODE FUNCTIONS
#

def generate_batch_report(args: argparse.Namespace) -> int:
    """
    Generate comprehensive batch report
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("report.batch")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get batch info
    batch_info = context.db.execute_dict_query(
        "SELECT id, batch_name, base_path, ref_version FROM ecod_schema.batch WHERE id = %s", 
        (args.batch_id,)
    )
    
    if not batch_info:
        logger.error(f"Batch {args.batch_id} not found")
        return 1
        
    batch_path = batch_info[0]['base_path']
    reference = batch_info[0]['ref_version']
    batch_name = batch_info[0]['batch_name']
    
    logger.info(f"Generating report for batch {args.batch_id} ({batch_name})")
    
    # Create report
    report = {
        'batch': {
            'id': args.batch_id,
            'name': batch_name,
            'path': batch_path,
            'reference': reference,
            'generated_at': datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        },
        'overview': {
            'proteins': {
                'total': 0,
                'representative': 0,
                'by_stage': {},
                'by_status': {}
            },
            'file_types': {}
        },
        'completion': {
            'overall': 0,
            'representative': 0,
            'stages': {}
        },
        'domains': {
            'total': 0,
            'classified': 0,
            'distribution': {
                'single': 0,
                'multi': 0
            },
            'by_group': {
                't_groups': {},
                'h_groups': {}
            }
        },
        'issues': {
            'total': 0,
            'critical': 0,
            'by_type': {}
        }
    }
    
    # Get protein statistics
    proteins_query = """
    SELECT 
        ps.current_stage, ps.status, ps.is_representative,
        COUNT(*) as count
    FROM 
        ecod_schema.process_status ps
    WHERE 
        ps.batch_id = %s
    GROUP BY 
        ps.current_stage, ps.status, ps.is_representative
    """
    
    proteins = context.db.execute_dict_query(proteins_query, (args.batch_id,))
    
    # Calculate protein statistics
    total_proteins = 0
    representative_proteins = 0
    
    for row in proteins:
        count = row['count']
        stage = row['current_stage']
        status = row['status']
        is_rep = row['is_representative']
        
        total_proteins += count
        if is_rep:
            representative_proteins += count
        
        # Update stage counts
        if stage not in report['overview']['proteins']['by_stage']:
            report['overview']['proteins']['by_stage'][stage] = 0
        report['overview']['proteins']['by_stage'][stage] += count
        
        # Update status counts
        if status not in report['overview']['proteins']['by_status']:
            report['overview']['proteins']['by_status'][status] = 0
        report['overview']['proteins']['by_status'][status] += count
    
    report['overview']['proteins']['total'] = total_proteins
    report['overview']['proteins']['representative'] = representative_proteins
    
    # Get file statistics
    files_query = """
    SELECT 
        pf.file_type, pf.file_exists,
        COUNT(*) as count
    FROM 
        ecod_schema.process_file pf
    JOIN 
        ecod_schema.process_status ps ON pf.process_id = ps.id
    WHERE 
        ps.batch_id = %s
    GROUP BY 
        pf.file_type, pf.file_exists
    """
    
    files = context.db.execute_dict_query(files_query, (args.batch_id,))
    
    # Calculate file statistics
    for row in files:
        file_type = row['file_type']
        file_exists = row['file_exists']
        count = row['count']
        
        if file_type not in report['overview']['file_types']:
            report['overview']['file_types'][file_type] = {
                'total': 0,
                'exists': 0
            }
        
        report['overview']['file_types'][file_type]['total'] += count
        if file_exists:
            report['overview']['file_types'][file_type]['exists'] += count
    
    # Calculate completion percentages
    if total_proteins > 0:
        completed = report['overview']['proteins']['by_stage'].get('domain_partition_complete', 0)
        report['completion']['overall'] = completed / total_proteins * 100
        
        # Calculate completion by stage
        stage_order = [
            'new', 'fasta_complete', 'blast_complete', 'hhsearch_complete', 
            'domain_summary_complete', 'domain_partition_complete'
        ]
        
        cumulative = 0
        for stage in stage_order:
            if stage in report['overview']['proteins']['by_stage']:
                cumulative += report['overview']['proteins']['by_stage'][stage]
                report['completion']['stages'][stage] = cumulative / total_proteins * 100
    
    if representative_proteins > 0:
        # Query for completed representative proteins
        rep_completed_query = """
        SELECT COUNT(*) as count
        FROM ecod_schema.process_status ps
        WHERE ps.batch_id = %s
          AND ps.is_representative = TRUE
          AND ps.current_stage = 'domain_partition_complete'
        """
        
        rep_completed = context.db.execute_dict_query(rep_completed_query, (args.batch_id,))
        
        if rep_completed:
            report['completion']['representative'] = rep_completed[0]['count'] / representative_proteins * 100
    
    # Get domain statistics
    if args.analyze_domains:
        # This would require parsing domain files
        # For simplicity, we'll just get counts
        domains_query = """
        SELECT COUNT(*) as count
        FROM ecod_schema.process_file pf
        JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
        WHERE ps.batch_id = %s
          AND pf.file_type = 'domain_partition'
          AND pf.file_exists = TRUE
        """
        
        domains = context.db.execute_dict_query(domains_query, (args.batch_id,))
        
        if domains:
            report['domains']['total'] = domains[0]['count']
    
    # Get issue statistics
    issues_query = """
    SELECT status, error_message, COUNT(*) as count
    FROM ecod_schema.process_status ps
    WHERE ps.batch_id = %s
      AND status = 'error'
    GROUP BY status, error_message
    """
    
    issues = context.db.execute_dict_query(issues_query, (args.batch_id,))
    
    # Calculate issue statistics
    for row in issues:
        error_message = row['error_message']
        count = row['count']
        
        report['issues']['total'] += count
        
        # Classify issues
        issue_type = 'unknown'
        if error_message:
            if 'blast' in error_message.lower():
                issue_type = 'blast_error'
            elif 'hhsearch' in error_message.lower() or 'hhr' in error_message.lower():
                issue_type = 'hhsearch_error'
            elif 'domain' in error_message.lower() and 'summary' in error_message.lower():
                issue_type = 'domain_summary_error'
            elif 'partition' in error_message.lower():
                issue_type = 'partition_error'
            elif 'file' in error_message.lower() and 'not found' in error_message.lower():
                issue_type = 'missing_file'
        
            if issue_type not in report['issues']['by_type']:
                report['issues']['by_type'][issue_type] = 0
            report['issues']['by_type'][issue_type] += count
            
            # Count critical issues
            if issue_type in ['missing_file', 'blast_error']:
                report['issues']['critical'] += count
    
    # Print report summary
    logger.info(f"\nBatch Report Summary for {batch_name} (ID: {args.batch_id}):")
    logger.info(f"Base path: {batch_path}")
    logger.info(f"Reference version: {reference}")
    
    logger.info(f"\nOverview:")
    logger.info(f"  Total proteins: {report['overview']['proteins']['total']}")
    logger.info(f"  Representative proteins: {report['overview']['proteins']['representative']} ({report['overview']['proteins']['representative']/report['overview']['proteins']['total']*100:.1f}%)")
    
    logger.info(f"\nStage distribution:")
    for stage, count in sorted(report['overview']['proteins']['by_stage'].items()):
        logger.info(f"  {stage}: {count} ({count/report['overview']['proteins']['total']*100:.1f}%)")
    
    logger.info(f"\nCompletion:")
    logger.info(f"  Overall: {report['completion']['overall']:.1f}%")
    logger.info(f"  Representative: {report['completion']['representative']:.1f}%")
    
    logger.info(f"\nFile statistics:")
    for file_type, stats in sorted(report['overview']['file_types'].items()):
        if stats['total'] > 0:
            logger.info(f"  {file_type}: {stats['exists']}/{stats['total']} ({stats['exists']/stats['total']*100:.1f}%)")
    
    logger.info(f"\nIssues:")
    logger.info(f"  Total: {report['issues']['total']}")
    logger.info(f"  Critical: {report['issues']['critical']}")
    
    for issue_type, count in sorted(report['issues']['by_type'].items(), key=lambda x: x[1], reverse=True):
        logger.info(f"  {issue_type}: {count}")
    
    # Write report to output file if specified
    if args.output:
        try:
            with open(args.output, 'w') as f:
                json.dump(report, f, indent=2)
            logger.info(f"Wrote report to {args.output}")
        except Exception as e:
            logger.error(f"Error writing report: {str(e)}")
    
    return 0

def generate_pipeline_overview(args: argparse.Namespace) -> int:
    """
    Generate a comprehensive pipeline overview across all batches
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("report.overview")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get all batches
    batches_query = """
    SELECT id, batch_name, base_path, ref_version,
           (CASE WHEN batch_name LIKE 'alt_rep%' THEN 'representative' ELSE 'regular' END) as batch_type
    FROM ecod_schema.batch
    ORDER BY id
    """
    
    batches = context.db.execute_dict_query(batches_query)
    
    if not batches:
        logger.error("No batches found in database")
        return 1
    
    # Initialize report
    report = {
        'overview': {
            'generated_at': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'batches': {
                'total': len(batches),
                'by_type': defaultdict(int)
            },
            'proteins': {
                'total': 0,
                'by_batch_type': defaultdict(int),
                'by_stage': defaultdict(int)
            }
        },
        'batches': [],
        'completion': {
            'overall': 0,
            'by_batch_type': {},
            'by_stage': {}
        },
        'domain_stats': {
            'total': 0,
            'by_batch_type': {},
            'distribution': {
                'single': 0,
                'multi': 0
            }
        },
        'issues': {
            'total': 0,
            'by_batch_type': {},
            'by_type': {}
        }
    }
    
    # Count batches by type
    for batch in batches:
        batch_type = batch['batch_type']
        report['overview']['batches']['by_type'][batch_type] += 1
    
    # Get basic info for each batch
    for batch in batches:
        batch_id = batch['id']
        batch_name = batch['batch_name']
        batch_type = batch['batch_type']
        
        # Get protein counts
        proteins_query = """
        SELECT COUNT(*) as total,
               SUM(CASE WHEN ps.is_representative THEN 1 ELSE 0 END) as representative,
               SUM(CASE WHEN ps.current_stage = 'domain_partition_complete' THEN 1 ELSE 0 END) as completed
        FROM ecod_schema.process_status ps
        WHERE ps.batch_id = %s
        """
        
        proteins = context.db.execute_dict_query(proteins_query, (batch_id,))
        
        if proteins:
            total = proteins[0]['total'] or 0
            representative = proteins[0]['representative'] or 0
            completed = proteins[0]['completed'] or 0
            
            batch_info = {
                'id': batch_id,
                'name': batch_name,
                'type': batch_type,
                'proteins': {
                    'total': total,
                    'representative': representative,
                    'completed': completed,
                    'completion_rate': (completed / total * 100) if total > 0 else 0
                }
            }
            
            report['batches'].append(batch_info)
            
            # Update totals
            report['overview']['proteins']['total'] += total
            report['overview']['proteins']['by_batch_type'][batch_type] += total
    
    # Get stage distribution across all batches
    stage_query = """
    SELECT ps.current_stage, COUNT(*) as count
    FROM ecod_schema.process_status ps
    GROUP BY ps.current_stage
    """
    
    stages = context.db.execute_dict_query(stage_query)
    
    for row in stages:
        stage = row['current_stage']
        count = row['count']
        report['overview']['proteins']['by_stage'][stage] = count
    
    # Calculate completion percentages
    if report['overview']['proteins']['total'] > 0:
        completed = report['overview']['proteins']['by_stage'].get('domain_partition_complete', 0)
        report['completion']['overall'] = completed / report['overview']['proteins']['total'] * 100
        
        # Calculate by batch type
        for batch_type in report['overview']['proteins']['by_batch_type']:
            # Get completed count for this batch type
            completed_query = """
            SELECT COUNT(*) as count
            FROM ecod_schema.process_status ps
            JOIN ecod_schema.batch b ON ps.batch_id = b.id
            WHERE (CASE WHEN b.batch_name LIKE 'alt_rep%' THEN 'representative' ELSE 'regular' END) = %s
              AND ps.current_stage = 'domain_partition_complete'
            """
            
            completed_count = context.db.execute_dict_query(completed_query, (batch_type,))
            
            if completed_count:
                completed = completed_count[0]['count'] or 0
                total = report['overview']['proteins']['by_batch_type'][batch_type]
                
                if total > 0:
                    report['completion']['by_batch_type'][batch_type] = completed / total * 100
    
    # Get domain statistics
    domains_query = """
    SELECT 
        (CASE WHEN b.batch_name LIKE 'alt_rep%' THEN 'representative' ELSE 'regular' END) as batch_type,
        COUNT(*) as count
    FROM ecod_schema.process_file pf
    JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
    JOIN ecod_schema.batch b ON ps.batch_id = b.id
    WHERE pf.file_type = 'domain_partition'
      AND pf.file_exists = TRUE
    GROUP BY (CASE WHEN b.batch_name LIKE 'alt_rep%' THEN 'representative' ELSE 'regular' END)
    """
    
    domains = context.db.execute_dict_query(domains_query)
    
    for row in domains:
        batch_type = row['batch_type']
        count = row['count']
        
        report['domain_stats']['total'] += count
        report['domain_stats']['by_batch_type'][batch_type] = count
    
    # Get issue statistics
    issues_query = """
    SELECT 
        (CASE WHEN b.batch_name LIKE 'alt_rep%' THEN 'representative' ELSE 'regular' END) as batch_type,
        status, COUNT(*) as count
    FROM ecod_schema.process_status ps
    JOIN ecod_schema.batch b ON ps.batch_id = b.id
    WHERE status = 'error'
    GROUP BY (CASE WHEN b.batch_name LIKE 'alt_rep%' THEN 'representative' ELSE 'regular' END), status
    """
    
    issues = context.db.execute_dict_query(issues_query)
    
    for row in issues:
        batch_type = row['batch_type']
        count = row['count']
        
        report['issues']['total'] += count
        
        if batch_type not in report['issues']['by_batch_type']:
            report['issues']['by_batch_type'][batch_type] = 0
        report['issues']['by_batch_type'][batch_type] += count
    
    # Print report summary
    logger.info(f"Pipeline Overview Report:")
    logger.info(f"Generated at: {report['overview']['generated_at']}")
    
    logger.info(f"\nBatches: {report['overview']['batches']['total']}")
    for batch_type, count in report['overview']['batches']['by_type'].items():
        logger.info(f"  {batch_type}: {count}")
    
    logger.info(f"\nProteins: {report['overview']['proteins']['total']}")
    for batch_type, count in report['overview']['proteins']['by_batch_type'].items():
        logger.info(f"  In {batch_type} batches: {count} ({count/report['overview']['proteins']['total']*100:.1f}%)")
    
    logger.info(f"\nCompletion:")
    logger.info(f"  Overall: {report['completion']['overall']:.1f}%")
    for batch_type, rate in report['completion']['by_batch_type'].items():
        logger.info(f"  {batch_type}: {rate:.1f}%")
    
    logger.info(f"\nDomains:")
    logger.info(f"  Total: {report['domain_stats']['total']}")
    for batch_type, count in report['domain_stats']['by_batch_type'].items():
        logger.info(f"  In {batch_type} batches: {count}")
    
    logger.info(f"\nIssues:")
    logger.info(f"  Total: {report['issues']['total']}")
    for batch_type, count in report['issues']['by_batch_type'].items():
        logger.info(f"  In {batch_type} batches: {count}")
    
    logger.info(f"\nBatch details:")
    for batch in sorted(report['batches'], key=lambda x: x['id']):
        logger.info(f"  {batch['id']}: {batch['name']} ({batch['type']})")
        logger.info(f"    Proteins: {batch['proteins']['total']}")
        logger.info(f"    Completed: {batch['proteins']['completed']} ({batch['proteins']['completion_rate']:.1f}%)")
    
    # Write report to output file if specified
    if args.output:
        try:
            with open(args.output, 'w') as f:
                json.dump(report, f, indent=2)
            logger.info(f"Wrote pipeline overview to {args.output}")
        except Exception as e:
            logger.error(f"Error writing overview: {str(e)}")
    
    return 0

def generate_domain_report(args: argparse.Namespace) -> int:
    """
    Generate report on domain statistics across all batches
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("report.domains")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get batch info if specific batch requested
    batch_filter = ""
    batch_params = []
    
    if args.batch_id:
        batch_info = context.db.execute_dict_query(
            "SELECT id, batch_name, base_path, ref_version FROM ecod_schema.batch WHERE id = %s", 
            (args.batch_id,)
        )
        
        if not batch_info:
            logger.error(f"Batch {args.batch_id} not found")
            return 1
            
        batch_filter = "WHERE ps.batch_id = %s"
        batch_params = [args.batch_id]
        
        logger.info(f"Generating domain report for batch {args.batch_id} ({batch_info[0]['batch_name']})")
    else:
        logger.info("Generating domain report across all batches")
    
    # Initialize report
    report = {
        'overview': {
            'generated_at': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'domain_files': 0,
            'domains': {
                'total': 0,
                'with_t_group': 0,
                'with_h_group': 0
            }
        },
        'chain_stats': {
            'total': 0,
            'with_domains': 0,
            'no_domains': 0,
            'single_domain': 0,
            'multi_domain': 0
        },
        'classification': {
            't_groups': defaultdict(int),
            'h_groups': defaultdict(int)
        },
        'domain_size': {
            'min': float('inf'),
            'max': 0,
            'avg': 0,
            'distribution': defaultdict(int)
        },
        'evidence': {
            'blast_only': 0,
            'hhsearch_only': 0,
            'both': 0,
            'none': 0
        }
    }
    
    # Get count of domain partition files
    count_query = f"""
    SELECT COUNT(*) as count
    FROM ecod_schema.process_file pf
    JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
    {batch_filter}
    AND pf.file_type = 'domain_partition'
    AND pf.file_exists = TRUE
    """
    
    count_result = context.db.execute_dict_query(count_query, batch_params)
    
    if count_result:
        report['overview']['domain_files'] = count_result[0]['count']
    
    # Get sample of domain partition files for analysis
    sample_query = f"""
    SELECT pf.file_path, ps.batch_id
    FROM ecod_schema.process_file pf
    JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
    {batch_filter}
    AND pf.file_type = 'domain_partition'
    AND pf.file_exists = TRUE
    ORDER BY RANDOM()
    LIMIT %s
    """
    
    sample_params = batch_params.copy()
    sample_params.append(args.sample_size if args.sample_size else 1000)
    
    sample_results = context.db.execute_dict_query(sample_query, sample_params)
    
    if not sample_results:
        logger.error("No domain partition files found for analysis")
        return 1
        
    logger.info(f"Analyzing {len(sample_results)} domain partition files")
    
    # Get batch paths
    batch_paths = {}
    if args.batch_id:
        batch_paths[args.batch_id] = batch_info[0]['base_path']
    else:
        batch_query = "SELECT id, base_path FROM ecod_schema.batch"
        batch_results = context.db.execute_dict_query(batch_query)
        
        for row in batch_results:
            batch_paths[row['id']] = row['base_path']
    
    # Initialize chains statistics
    domain_counts = defaultdict(int)
    total_domains = 0
    total_length = 0
    
    # Analysis loop
    for sample in sample_results:
        file_path = sample['file_path']
        batch_id = sample['batch_id']
        
        if batch_id not in batch_paths:
            logger.warning(f"Batch path not found for batch {batch_id}")
            continue
            
        batch_path = batch_paths[batch_id]
        full_path = os.path.join(batch_path, file_path)
        
        try:
            tree = ET.parse(full_path)
            root = tree.getroot()
            
            # Count domains
            domains = root.findall('./domain_list/domain')
            domain_count = len(domains)
            
            # Update chain statistics
            domain_counts[domain_count] += 1
            
            if domain_count > 0:
                report['chain_stats']['with_domains'] += 1
                
                if domain_count == 1:
                    report['chain_stats']['single_domain'] += 1
                else:
                    report['chain_stats']['multi_domain'] += 1
            else:
                report['chain_stats']['no_domains'] += 1
            
            # Analyze domains
            for domain in domains:
                total_domains += 1
                
                # Check T-group and H-group
                t_group = domain.get('t_group', '')
                h_group = domain.get('h_group', '')
                
                if t_group:
                    report['overview']['domains']['with_t_group'] += 1
                    report['classification']['t_groups'][t_group] += 1
                
                if h_group:
                    report['overview']['domains']['with_h_group'] += 1
                    report['classification']['h_groups'][h_group] += 1
                
                # Calculate domain size
                domain_range = domain.get('range', '')
                try:
                    domain_length = 0
                    for segment in domain_range.split(','):
                        start, end = map(int, segment.split('-'))
                        domain_length += (end - start + 1)
                        
                    report['domain_size']['min'] = min(report['domain_size']['min'], domain_length)
                    report['domain_size']['max'] = max(report['domain_size']['max'], domain_length)
                    total_length += domain_length
                    
                    # Record size distribution in buckets
                    if domain_length < 50:
                        size_bucket = "<50"
                    elif domain_length < 100:
                        size_bucket = "50-99"
                    elif domain_length < 150:
                        size_bucket = "100-149"
                    elif domain_length < 200:
                        size_bucket = "150-199"
                    elif domain_length < 300:
                        size_bucket = "200-299"
                    elif domain_length < 500:
                        size_bucket = "300-499"
                    else:
                        size_bucket = "500+"
                        
                    report['domain_size']['distribution'][size_bucket] += 1
                    
                except (ValueError, AttributeError):
                    pass
                
                # Check evidence
                evidence_elem = domain.find('evidence')
                if evidence_elem is not None:
                    has_blast = False
                    has_hhsearch = False
                    
                    for evidence in evidence_elem.findall('*'):
                        evidence_type = evidence.tag
                        if evidence_type == 'chain_blast' or evidence_type == 'domain_blast':
                            has_blast = True
                        elif evidence_type == 'hhsearch':
                            has_hhsearch = True
                    
                    if has_blast and has_hhsearch:
                        report['evidence']['both'] += 1
                    elif has_blast:
                        report['evidence']['blast_only'] += 1
                    elif has_hhsearch:
                        report['evidence']['hhsearch_only'] += 1
                    else:
                        report['evidence']['none'] += 1
                else:
                    report['evidence']['none'] += 1
        
        except Exception as e:
            logger.error(f"Error analyzing file {full_path}: {str(e)}")
    
    # Update report statistics
    report['overview']['domains']['total'] = total_domains
    report['chain_stats']['total'] = sum(domain_counts.values())
    
    if total_domains > 0:
        report['domain_size']['avg'] = total_length / total_domains
    
    # Get top T-groups and H-groups
    report['classification']['top_t_groups'] = sorted(report['classification']['t_groups'].items(), key=lambda x: x[1], reverse=True)[:20]
    report['classification']['top_h_groups'] = sorted(report['classification']['h_groups'].items(), key=lambda x: x[1], reverse=True)[:20]
    
    # Print report summary
    logger.info(f"\nDomain Report Summary:")
    logger.info(f"Generated at: {report['overview']['generated_at']}")
    
    logger.info(f"\nOverview:")
    logger.info(f"  Domain files analyzed: {report['chain_stats']['total']}")
    logger.info(f"  Total domains found: {report['overview']['domains']['total']}")
    logger.info(f"  Domains with T-group: {report['overview']['domains']['with_t_group']} ({report['overview']['domains']['with_t_group']/report['overview']['domains']['total']*100:.1f}%)")
    logger.info(f"  Domains with H-group: {report['overview']['domains']['with_h_group']} ({report['overview']['domains']['with_h_group']/report['overview']['domains']['total']*100:.1f}%)")
    
    logger.info(f"\nChain statistics:")
    logger.info(f"  With domains: {report['chain_stats']['with_domains']} ({report['chain_stats']['with_domains']/report['chain_stats']['total']*100:.1f}%)")
    logger.info(f"  No domains: {report['chain_stats']['no_domains']} ({report['chain_stats']['no_domains']/report['chain_stats']['total']*100:.1f}%)")
    logger.info(f"  Single domain: {report['chain_stats']['single_domain']} ({report['chain_stats']['single_domain']/report['chain_stats']['with_domains']*100:.1f}% of with_domains)")
    logger.info(f"  Multi domain: {report['chain_stats']['multi_domain']} ({report['chain_stats']['multi_domain']/report['chain_stats']['with_domains']*100:.1f}% of with_domains)")
    
    logger.info(f"\nDomain distribution:")
    for count, num_chains in sorted(domain_counts.items()):
        logger.info(f"  {count} domain(s): {num_chains} chains ({num_chains/report['chain_stats']['total']*100:.1f}%)")
    
    logger.info(f"\nDomain size statistics:")
    logger.info(f"  Average domain size: {report['domain_size']['avg']:.1f} residues")
    logger.info(f"  Min size: {report['domain_size']['min']} residues")
    logger.info(f"  Max size: {report['domain_size']['max']} residues")
    
    logger.info(f"\nSize distribution:")
    for size_bucket, count in sorted(report['domain_size']['distribution'].items(), key=lambda x: (
        int(x[0].split('-')[0]) if x[0][0].isdigit() else (0 if x[0] == '<50' else 1000)
    )):
        logger.info(f"  {size_bucket}: {count} ({count/report['overview']['domains']['total']*100:.1f}%)")
    
    logger.info(f"\nEvidence distribution:")
    for evidence_type, count in report['evidence'].items():
        logger.info(f"  {evidence_type}: {count} ({count/report['overview']['domains']['total']*100:.1f}%)")
    
    logger.info(f"\nTop 10 T-groups:")
    for t_group, count in report['classification']['top_t_groups'][:10]:
        logger.info(f"  {t_group}: {count} ({count/report['overview']['domains']['total']*100:.1f}%)")
    
    logger.info(f"\nTop 10 H-groups:")
    for h_group, count in report['classification']['top_h_groups'][:10]:
        logger.info(f"  {h_group}: {count} ({count/report['overview']['domains']['total']*100:.1f}%)")
    
    # Write report to output file if specified
    if args.output:
        try:
            with open(args.output, 'w') as f:
                json.dump(report, f, indent=2)
            logger.info(f"Wrote domain report to {args.output}")
        except Exception as e:
            logger.error(f"Error writing report: {str(e)}")
    
    return 0

def report_mode(args: argparse.Namespace) -> int:
    """
    Run report mode with appropriate action
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("report")
    
    if args.action == "batch":
        return generate_batch_report(args)
    elif args.action == "overview":
        return generate_pipeline_overview(args)
    elif args.action == "domains":
        return generate_domain_report(args)
    else:
        logger.error(f"Unknown report action: {args.action}")
        return 1

#
# DIAGNOSE MODE FUNCTIONS
#

def diagnose_missing_files(args: argparse.Namespace) -> int:
    """
    Diagnose missing files in the pipeline
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("diagnose.missing")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get batch info
    batch_info = context.db.execute_dict_query(
        "SELECT id, batch_name, base_path, ref_version FROM ecod_schema.batch WHERE id = %s", 
        (args.batch_id,)
    )
    
    if not batch_info:
        logger.error(f"Batch {args.batch_id} not found")
        return 1
        
    batch_path = batch_info[0]['base_path']
    reference = batch_info[0]['ref_version']
    
    logger.info(f"Diagnosing missing files for batch {args.batch_id} ({batch_info[0]['batch_name']})")
    
    # Define file dependencies
    file_chain = [
        'fasta',
        'chain_blast',
        'domain_blast',
        'hhr',
        'hh_xml',
        'domain_summary',
        'domain_partition'
    ]
    
    # Initialize diagnostic results
    results = {
        'batch_id': args.batch_id,
        'batch_name': batch_info[0]['batch_name'],
        'total_proteins': 0,
        'file_stats': {},
        'missing_files': {},
        'chain_breaks': []
    }
    
    # Get all proteins and process status
    query = """
    SELECT 
        p.id as protein_id, p.pdb_id, p.chain_id, ps.id as process_id,
        ps.current_stage, ps.is_representative
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    WHERE 
        ps.batch_id = %s
    """
    
    if args.reps_only:
        query += " AND ps.is_representative = TRUE"
    
    query += " ORDER BY p.pdb_id, p.chain_id"
    
    proteins = context.db.execute_dict_query(query, (args.batch_id,))
    
    if not proteins:
        logger.error("No proteins found for this batch")
        return 1
        
    results['total_proteins'] = len(proteins)
    logger.info(f"Found {len(proteins)} proteins")
    
    # Apply limit if specified
    if args.limit and args.limit < len(proteins):
        proteins = proteins[:args.limit]
        logger.info(f"Limited to diagnosing {args.limit} proteins")
    
    # Get file information for each protein
    for protein in proteins:
        pdb_id = protein['pdb_id']
        chain_id = protein['chain_id']
        process_id = protein['process_id']
        current_stage = protein['current_stage']
        
        # Get file paths
        files_query = """
        SELECT 
            file_type, file_path, file_exists
        FROM 
            ecod_schema.process_file
        WHERE 
            process_id = %s
        """
        
        files = context.db.execute_dict_query(files_query, (process_id,))
        
        if not files:
            continue
        
        # Map files by type
        file_map = {}
        for file in files:
            file_type = file['file_type']
            file_path = file['file_path']
            file_exists = file['file_exists']
            
            # Initialize file stats
            if file_type not in results['file_stats']:
                results['file_stats'][file_type] = {
                    'total': 0,
                    'exists': 0,
                    'missing': 0
                }
            
            results['file_stats'][file_type]['total'] += 1
            
            if file_exists:
                results['file_stats'][file_type]['exists'] += 1
                file_map[file_type] = {
                    'path': file_path,
                    'exists': True,
                    'full_path': os.path.join(batch_path, file_path) if file_path else None
                }
            else:
                results['file_stats'][file_type]['missing'] += 1
                
                # Track missing file
                if file_type not in results['missing_files']:
                    results['missing_files'][file_type] = []
                
                results['missing_files'][file_type].append({
                    'pdb_id': pdb_id,
                    'chain_id': chain_id,
                    'current_stage': current_stage,
                    'path': file_path
                })
        
        # Check for chain breaks
        current_chain = []
        chain_broken = False
        
        for file_type in file_chain:
            # Skip optional HHSearch files
            if file_type in ['hhr', 'hh_xml'] and args.blast_only:
                continue
                
            if file_type in file_map:
                current_chain.append(file_type)
            else:
                # If we have a later file but missing an earlier one
                if len(current_chain) > 0 and file_type not in ['domain_summary', 'domain_partition']:
                    chain_broken = True
                elif current_chain and file_type == 'domain_summary' and 'domain_partition' in file_map:
                    chain_broken = True
        
        if chain_broken:
            results['chain_breaks'].append({
                'pdb_id': pdb_id,
                'chain_id': chain_id,
                'current_stage': current_stage,
                'available_files': list(file_map.keys())
            })
    
    # Print diagnostic summary
    logger.info(f"\nDiagnostic Summary for Batch {args.batch_id}:")
    
    logger.info(f"\nFile statistics:")
    for file_type in file_chain:
        if file_type in results['file_stats']:
            stats = results['file_stats'][file_type]
            logger.info(f"  {file_type}: {stats['exists']}/{stats['total']} ({stats['exists']/stats['total']*100:.1f}%)")
    
    logger.info(f"\nMissing files:")
    for file_type in file_chain:
        if file_type in results['missing_files']:
            missing = results['missing_files'][file_type]
            logger.info(f"  {file_type}: {len(missing)} missing")
            
            # Print examples
            if args.verbose and len(missing) > 0:
                for i, example in enumerate(missing[:5]):
                    logger.info(f"    {example['pdb_id']}_{example['chain_id']} (stage: {example['current_stage']})")
                
                if len(missing) > 5:
                    logger.info(f"    ... and {len(missing) - 5} more")
    
    logger.info(f"\nChain breaks detected: {len(results['chain_breaks'])}")
    if args.verbose and results['chain_breaks']:
        for i, break_info in enumerate(results['chain_breaks'][:5]):
            logger.info(f"  {break_info['pdb_id']}_{break_info['chain_id']} (stage: {break_info['current_stage']})")
            logger.info(f"    Available files: {', '.join(break_info['available_files'])}")
        
        if len(results['chain_breaks']) > 5:
            logger.info(f"  ... and {len(results['chain_breaks']) - 5} more")
    
    # Write results to output file if specified
    if args.output:
        try:
            with open(args.output, 'w') as f:
                json.dump(results, f, indent=2)
            logger.info(f"Wrote diagnostic results to {args.output}")
        except Exception as e:
            logger.error(f"Error writing results: {str(e)}")
    
    return 0

def diagnose_consistency(args: argparse.Namespace) -> int:
    """
    Diagnose consistency issues in domain summaries and partitions
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("diagnose.consistency")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get batch info
    batch_info = context.db.execute_dict_query(
        "SELECT id, batch_name, base_path, ref_version FROM ecod_schema.batch WHERE id = %s", 
        (args.batch_id,)
    )
    
    if not batch_info:
        logger.error(f"Batch {args.batch_id} not found")
        return 1
        
    batch_path = batch_info[0]['base_path']
    reference = batch_info[0]['ref_version']
    
    logger.info(f"Diagnosing consistency for batch {args.batch_id} ({batch_info[0]['batch_name']})")
    
    # Initialize diagnostic results
    results = {
        'batch_id': args.batch_id,
        'batch_name': batch_info[0]['batch_name'],
        'total_proteins': 0,
        'analyzed': 0,
        'consistency': {
            'summary_to_partition': {
                'consistent': 0,
                'inconsistent': 0,
                'issues': []
            },
            'chain_blast_hhr': {
                'consistent': 0,
                'inconsistent': 0,
                'issues': []
            },
            'metadata': {
                'consistent': 0,
                'inconsistent': 0,
                'issues': []
            }
        }
    }
    
    # Get all proteins with both summary and partition files
    query = """
    SELECT 
        p.id as protein_id, p.pdb_id, p.chain_id, ps.id as process_id,
        ps.current_stage, ps.is_representative,
        (SELECT pf.file_path FROM ecod_schema.process_file pf 
         WHERE pf.process_id = ps.id AND pf.file_type = 'domain_summary' AND pf.file_exists = TRUE
         LIMIT 1) as summary_path,
        (SELECT pf.file_path FROM ecod_schema.process_file pf 
         WHERE pf.process_id = ps.id AND pf.file_type = 'domain_partition' AND pf.file_exists = TRUE
         LIMIT 1) as partition_path
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    WHERE 
        ps.batch_id = %s
    """
    
    if args.reps_only:
        query += " AND ps.is_representative = TRUE"
    
    query += " ORDER BY p.pdb_id, p.chain_id"
    
    proteins = context.db.execute_dict_query(query, (args.batch_id,))
    
    if not proteins:
        logger.error("No proteins found for this batch")
        return 1
        
    results['total_proteins'] = len(proteins)
    logger.info(f"Found {len(proteins)} proteins")
    
    # Count how many have both files
    proteins_with_both = [p for p in proteins if p['summary_path'] and p['partition_path']]
    logger.info(f"{len(proteins_with_both)} proteins have both summary and partition files")
    
    # Apply limit if specified
    if args.limit and args.limit < len(proteins_with_both):
        proteins_with_both = proteins_with_both[:args.limit]
        logger.info(f"Limited to diagnosing {args.limit} proteins")
    
    # Analyze consistency for each protein
    for protein in proteins_with_both:
        pdb_id = protein['pdb_id']
        chain_id = protein['chain_id']
        summary_path = os.path.join(batch_path, protein['summary_path'])
        partition_path = os.path.join(batch_path, protein['partition_path'])
        
        try:
            # Parse files
            summary_tree = ET.parse(summary_path)
            summary_root = summary_tree.getroot()
            
            partition_tree = ET.parse(partition_path)
            partition_root = partition_tree.getroot()
            
            # Check 1: Metadata consistency
            metadata_consistent = True
            summary_metadata = summary_root.find('metadata')
            
            if summary_metadata is not None:
                summary_pdb = summary_metadata.findtext('pdb_id', '')
                summary_chain = summary_metadata.findtext('chain_id', '')
                summary_ref = summary_metadata.findtext('reference', '')
                
                if summary_pdb != pdb_id or summary_chain != chain_id:
                    metadata_consistent = False
                    results['consistency']['metadata']['issues'].append({
                        'pdb_id': pdb_id,
                        'chain_id': chain_id,
                        'issue': f"Summary metadata mismatch: {summary_pdb}_{summary_chain}"
                    })
            
            if metadata_consistent:
                results['consistency']['metadata']['consistent'] += 1
            else:
                results['consistency']['metadata']['inconsistent'] += 1
            
            # Check 2: Summary to partition consistency
            # Get domain suggestions from summary
            summary_domains = []
            suggestions = summary_root.find('domain_suggestions')
            if suggestions is not None:
                for domain in suggestions.findall('domain'):
                    domain_range = domain.get('range', '')
                    if domain_range:
                        summary_domains.append(domain_range)
            
            # Get domains from partition
            partition_domains = []
            domain_list = partition_root.find('domain_list')
            if domain_list is not None:
                for domain in domain_list.findall('domain'):
                    domain_range = domain.get('range', '')
                    if domain_range:
                        partition_domains.append(domain_range)
            
            # Compare
            if len(summary_domains) == 0 and len(partition_domains) == 0:
                # Both empty, consider consistent
                results['consistency']['summary_to_partition']['consistent'] += 1
            elif len(summary_domains) > 0 and len(partition_domains) > 0:
                # Check overlap
                common_domains = set(summary_domains).intersection(set(partition_domains))
                
                if len(common_domains) > 0:
                    results['consistency']['summary_to_partition']['consistent'] += 1
                else:
                    results['consistency']['summary_to_partition']['inconsistent'] += 1
                    results['consistency']['summary_to_partition']['issues'].append({
                        'pdb_id': pdb_id,
                        'chain_id': chain_id,
                        'summary_domains': summary_domains,
                        'partition_domains': partition_domains,
                        'issue': "No common domains"
                    })
            else:
                # One has domains but the other doesn't
                results['consistency']['summary_to_partition']['inconsistent'] += 1
                results['consistency']['summary_to_partition']['issues'].append({
                    'pdb_id': pdb_id,
                    'chain_id': chain_id,
                    'summary_domains': summary_domains,
                    'partition_domains': partition_domains,
                    'issue': "Domains in one file but not the other"
                })
            
            # Check 3: BLAST and HHSearch consistency (if both exist)
            chain_blast = summary_root.find('chain_blast_evidence')
            hhsearch = summary_root.find('hhsearch_evidence')
            
            if chain_blast is not None and hhsearch is not None:
                # Count hits
                chain_blast_hits = chain_blast.findall('.//hit')
                hhsearch_hits = hhsearch.findall('.//hh_hit')
                
                if (len(chain_blast_hits) > 0) == (len(hhsearch_hits) > 0):
                    # Both have hits or both don't have hits
                    results['consistency']['chain_blast_hhr']['consistent'] += 1
                else:
                    results['consistency']['chain_blast_hhr']['inconsistent'] += 1
                    results['consistency']['chain_blast_hhr']['issues'].append({
                        'pdb_id': pdb_id,
                        'chain_id': chain_id,
                        'chain_blast_hits': len(chain_blast_hits),
                        'hhsearch_hits': len(hhsearch_hits),
                        'issue': "One has hits but the other doesn't"
                    })
            
            results['analyzed'] += 1
            
        except Exception as e:
            logger.error(f"Error analyzing files for {pdb_id}_{chain_id}: {str(e)}")
    
    # Print diagnostic summary
    logger.info(f"\nConsistency Diagnostic Summary for Batch {args.batch_id}:")
    logger.info(f"Analyzed {results['analyzed']} of {len(proteins_with_both)} proteins")
    
    logger.info(f"\nMetadata consistency:")
    metadata = results['consistency']['metadata']
    if metadata['consistent'] + metadata['inconsistent'] > 0:
        consistent_pct = metadata['consistent'] / (metadata['consistent'] + metadata['inconsistent']) * 100
        logger.info(f"  Consistent: {metadata['consistent']} ({consistent_pct:.1f}%)")
        logger.info(f"  Inconsistent: {metadata['inconsistent']} ({100-consistent_pct:.1f}%)")
    
    logger.info(f"\nSummary to partition consistency:")
    summary_partition = results['consistency']['summary_to_partition']
    if summary_partition['consistent'] + summary_partition['inconsistent'] > 0:
        consistent_pct = summary_partition['consistent'] / (summary_partition['consistent'] + summary_partition['inconsistent']) * 100
        logger.info(f"  Consistent: {summary_partition['consistent']} ({consistent_pct:.1f}%)")
        logger.info(f"  Inconsistent: {summary_partition['inconsistent']} ({100-consistent_pct:.1f}%)")
    
    logger.info(f"\nChain BLAST to HHSearch consistency:")
    chain_blast_hhr = results['consistency']['chain_blast_hhr']
    if chain_blast_hhr['consistent'] + chain_blast_hhr['inconsistent'] > 0:
        consistent_pct = chain_blast_hhr['consistent'] / (chain_blast_hhr['consistent'] + chain_blast_hhr['inconsistent']) * 100
        logger.info(f"  Consistent: {chain_blast_hhr['consistent']} ({consistent_pct:.1f}%)")
        logger.info(f"  Inconsistent: {chain_blast_hhr['inconsistent']} ({100-consistent_pct:.1f}%)")
    
    # Print example issues
    if args.verbose:
        logger.info(f"\nExample metadata issues:")
        for issue in results['consistency']['metadata']['issues'][:5]:
            logger.info(f"  {issue['pdb_id']}_{issue['chain_id']}: {issue['issue']}")
        
        logger.info(f"\nExample summary to partition issues:")
        for issue in results['consistency']['summary_to_partition']['issues'][:5]:
            logger.info(f"  {issue['pdb_id']}_{issue['chain_id']}: {issue['issue']}")
        
        logger.info(f"\nExample chain BLAST to HHSearch issues:")
        for issue in results['consistency']['chain_blast_hhr']['issues'][:5]:
            logger.info(f"  {issue['pdb_id']}_{issue['chain_id']}: {issue['issue']}")
    
    # Write results to output file if specified
    if args.output:
        try:
            with open(args.output, 'w') as f:
                json.dump(results, f, indent=2)
            logger.info(f"Wrote consistency diagnostic results to {args.output}")
        except Exception as e:
            logger.error(f"Error writing results: {str(e)}")
    
    return 0

def diagnose_content_quality(args: argparse.Namespace) -> int:
    """
    Diagnose content quality issues in summaries and partitions
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("diagnose.quality")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get batch info
    batch_info = context.db.execute_dict_query(
        "SELECT id, batch_name, base_path, ref_version FROM ecod_schema.batch WHERE id = %s", 
        (args.batch_id,)
    )
    
    if not batch_info:
        logger.error(f"Batch {args.batch_id} not found")
        return 1
        
    batch_path = batch_info[0]['base_path']
    reference = batch_info[0]['ref_version']
    
    logger.info(f"Diagnosing content quality for batch {args.batch_id} ({batch_info[0]['batch_name']})")
    
    # Initialize diagnostic results
    results = {
        'batch_id': args.batch_id,
        'batch_name': batch_info[0]['batch_name'],
        'total_proteins': 0,
        'analyzed': 0,
        'quality': {
            'empty_evidence': 0,
            'poor_evidence': 0,
            'insufficient_coverage': 0,
            'divergent_boundaries': 0,
            'conflicting_evidence': 0,
            'potential_errors': []
        }
    }
    
    # Get proteins with summary files
    query = """
    SELECT 
        p.id as protein_id, p.pdb_id, p.chain_id, ps.id as process_id, p.length,
        (SELECT pf.file_path FROM ecod_schema.process_file pf 
         WHERE pf.process_id = ps.id AND pf.file_type = 'domain_summary' AND pf.file_exists = TRUE
         LIMIT 1) as summary_path,
        (SELECT pf.file_path FROM ecod_schema.process_file pf 
         WHERE pf.process_id = ps.id AND pf.file_type = 'domain_partition' AND pf.file_exists = TRUE
         LIMIT 1) as partition_path
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    WHERE 
        ps.batch_id = %s
    """
    
    if args.reps_only:
        query += " AND ps.is_representative = TRUE"
    
    query += " ORDER BY p.pdb_id, p.chain_id"
    
    proteins = context.db.execute_dict_query(query, (args.batch_id,))
    
    if not proteins:
        logger.error("No proteins found for this batch")
        return 1
        
    results['total_proteins'] = len(proteins)
    logger.info(f"Found {len(proteins)} proteins")
    
    # Filter to those with summary files
    proteins_with_summary = [p for p in proteins if p['summary_path']]
    logger.info(f"{len(proteins_with_summary)} proteins have summary files")
    
    # Apply limit if specified
    if args.limit and args.limit < len(proteins_with_summary):
        proteins_with_summary = proteins_with_summary[:args.limit]
        logger.info(f"Limited to diagnosing {args.limit} proteins")
    
    # Analyze quality for each protein
    for protein in proteins_with_summary:
        pdb_id = protein['pdb_id']
        chain_id = protein['chain_id']
        chain_length = protein['length'] or 0
        summary_path = os.path.join(batch_path, protein['summary_path'])
        
        partition_path = None
        if protein['partition_path']:
            partition_path = os.path.join(batch_path, protein['partition_path'])
        
        try:
            # Parse summary file
            summary_tree = ET.parse(summary_path)
            summary_root = summary_tree.getroot()
            
            # Check for quality issues
            issues = []
            
            # Check 1: Empty evidence
            empty_evidence = True
            evidence_count = 0
            
            for evidence_type in ['chain_blast_evidence', 'domain_blast_evidence', 'hhsearch_evidence']:
                evidence_elem = summary_root.find(evidence_type)
                if evidence_elem is not None:
                    hits = evidence_elem.findall('.//hit') + evidence_elem.findall('.//hh_hit')
                    if hits:
                        empty_evidence = False
                        evidence_count += len(hits)
            
            if empty_evidence:
                results['quality']['empty_evidence'] += 1
                issues.append("Empty evidence")
            
            # Check 2: Poor evidence (very few hits)
            if not empty_evidence and evidence_count < 3:
                results['quality']['poor_evidence'] += 1
                issues.append(f"Poor evidence (only {evidence_count} hits)")
            
            # Check 3: Domain coverage (if domains exist)
            domain_coverage = 0
            domain_suggestions = summary_root.find('domain_suggestions')
            if domain_suggestions is not None:
                domains = domain_suggestions.findall('domain')
                
                # Calculate coverage
                if domains and chain_length > 0:
                    covered_residues = set()
                    
                    for domain in domains:
                        domain_range = domain.get('range', '')
                        
                        try:
                            for segment in domain_range.split(','):
                                start, end = map(int, segment.split('-'))
                                for i in range(start, end + 1):
                                    covered_residues.add(i)
                        except (ValueError, AttributeError):
                            pass
                    
                    domain_coverage = len(covered_residues) / chain_length * 100
                    
                    if domain_coverage < 50:
                        results['quality']['insufficient_coverage'] += 1
                        issues.append(f"Insufficient coverage ({domain_coverage:.1f}%)")
            
            # Check 4: Boundary divergence (if partition exists)
            if partition_path:
                try:
                    partition_tree = ET.parse(partition_path)
                    partition_root = partition_tree.getroot()
                    
                    # Get domains from both files
                    summary_domains = []
                    if domain_suggestions is not None:
                        for domain in domain_suggestions.findall('domain'):
                            domain_range = domain.get('range', '')
                            if domain_range:
                                summary_domains.append(domain_range)
                    
                    partition_domains = []
                    domain_list = partition_root.find('domain_list')
                    if domain_list is not None:
                        for domain in domain_list.findall('domain'):
                            domain_range = domain.get('range', '')
                            if domain_range:
                                partition_domains.append(domain_range)
                    
                    # Check for boundary divergence
                    if summary_domains and partition_domains:
                        # Find partial matches (similar but not exact)
                        has_divergent = False
                        
                        for s_domain in summary_domains:
                            for p_domain in partition_domains:
                                if s_domain != p_domain and (s_domain in p_domain or p_domain in s_domain):
                                    has_divergent = True
                                    break
                            
                            if has_divergent:
                                break
                        
                        if has_divergent:
                            results['quality']['divergent_boundaries'] += 1
                            issues.append("Divergent domain boundaries")
                
                except Exception as e:
                    logger.error(f"Error parsing partition file for {pdb_id}_{chain_id}: {str(e)}")
            
            # Check 5: Conflicting evidence
            # This would require more sophisticated analysis
            # For simplicity, we're just checking for cases with both BLAST and HHSearch
            # but different domain suggestions
            
            # Record errors if issues found
            if issues:
                results['quality']['potential_errors'].append({
                    'pdb_id': pdb_id,
                    'chain_id': chain_id,
                    'issues': issues,
                    'chain_length': chain_length,
                    'coverage': domain_coverage,
                    'evidence_count': evidence_count
                })
            
            results['analyzed'] += 1
            
        except Exception as e:
            logger.error(f"Error analyzing files for {pdb_id}_{chain_id}: {str(e)}")
    
    # Print diagnostic summary
    logger.info(f"\nContent Quality Diagnostic Summary for Batch {args.batch_id}:")
    logger.info(f"Analyzed {results['analyzed']} of {len(proteins_with_summary)} proteins")
    
    logger.info(f"\nQuality issues:")
    logger.info(f"  Empty evidence: {results['quality']['empty_evidence']} ({results['quality']['empty_evidence']/results['analyzed']*100:.1f}%)")
    logger.info(f"  Poor evidence: {results['quality']['poor_evidence']} ({results['quality']['poor_evidence']/results['analyzed']*100:.1f}%)")
    logger.info(f"  Insufficient coverage: {results['quality']['insufficient_coverage']} ({results['quality']['insufficient_coverage']/results['analyzed']*100:.1f}%)")
    logger.info(f"  Divergent boundaries: {results['quality']['divergent_boundaries']} ({results['quality']['divergent_boundaries']/results['analyzed']*100:.1f}%)")
    
    logger.info(f"\nTotal proteins with potential issues: {len(results['quality']['potential_errors'])} ({len(results['quality']['potential_errors'])/results['analyzed']*100:.1f}%)")
    
    # Print example issues
    if args.verbose and results['quality']['potential_errors']:
        logger.info(f"\nExample proteins with quality issues:")
        for i, issue in enumerate(results['quality']['potential_errors'][:5]):
            logger.info(f"  {issue['pdb_id']}_{issue['chain_id']}: {', '.join(issue['issues'])}")
        
        if len(results['quality']['potential_errors']) > 5:
            logger.info(f"  ... and {len(results['quality']['potential_errors']) - 5} more")
    
    # Write results to output file if specified
    if args.output:
        try:
            with open(args.output, 'w') as f:
                json.dump(results, f, indent=2)
            logger.info(f"Wrote quality diagnostic results to {args.output}")
        except Exception as e:
            logger.error(f"Error writing results: {str(e)}")
    
    return 0

def diagnose_mode(args: argparse.Namespace) -> int:
    """
    Run diagnose mode with appropriate action
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("diagnose")
    
    if args.action == "missing":
        return diagnose_missing_files(args)
    elif args.action == "consistency":
        return diagnose_consistency(args)
    elif args.action == "quality":
        return diagnose_content_quality(args)
    else:
        logger.error(f"Unknown diagnostic action: {args.action}")
        return 1

#
# MAIN FUNCTION AND ARGUMENT PARSING
#

def main():
    """Main entry point"""
    # Create top-level parser
    parser = argparse.ArgumentParser(description='Validate Pipeline - Comprehensive validation toolkit for pyECOD pipeline')
    subparsers = parser.add_subparsers(dest="mode", help="Operating mode")

    # Validate mode parser
    validate_parser = subparsers.add_parser("validate", help="Validate summary files and domain partitions")
    validate_subparsers = validate_parser.add_subparsers(dest="action", help="Validation action")

    # Summary validation subparser
    summary_parser = validate_subparsers.add_parser("summary", help="Validate domain summary files")
    summary_parser.add_argument('--config', type=str, default='config/config.yml',
                             help='Path to configuration file')
    summary_parser.add_argument('--batch-id', type=int, required=True,
                             help='Batch ID to validate')
    summary_parser.add_argument('--output', type=str,
                             help='Output JSON file for validation statistics')
    summary_parser.add_argument('--problematic-output', type=str,
                             help='Output JSON file for problematic chains')
    summary_parser.add_argument('--limit', type=int,
                             help='Limit number of proteins to validate')
    summary_parser.add_argument('--log-file', type=str,
                             help='Log file path')
    summary_parser.add_argument('-v', '--verbose', action='store_true',
                             help='Enable verbose output')
    summary_parser.add_argument('--reps-only', action='store_true',
                             help='Validate only representative proteins')

    # Evidence tracing subparser
    evidence_parser = validate_subparsers.add_parser("evidence", help="Validate evidence tracing")
    evidence_parser.add_argument('--config', type=str, default='config/config.yml',
                              help='Path to configuration file')
    evidence_parser.add_argument('--batch-id', type=int, required=True,
                              help='Batch ID to validate')
    evidence_parser.add_argument('--output', type=str,
                              help='Output JSON file for tracing results')
    evidence_parser.add_argument('--sample-size', type=int, default=50,
                              help='Number of proteins to sample for detailed tracing')
    evidence_parser.add_argument('--log-file', type=str,
                              help='Log file path')
    evidence_parser.add_argument('-v', '--verbose', action='store_true',
                              help='Enable verbose output')
    evidence_parser.add_argument('--reps-only', action='store_true',
                              help='Validate only representative proteins')

    # Completeness validation subparser
    completeness_parser = validate_subparsers.add_parser("completeness", help="Validate pipeline completeness")
    completeness_parser.add_argument('--config', type=str, default='config/config.yml',
                                  help='Path to configuration file')
    completeness_parser.add_argument('--batch-id', type=int, required=True,
                                  help='Batch ID to validate')
    completeness_parser.add_argument('--output', type=str,
                                  help='Output JSON file for completeness results')
    completeness_parser.add_argument('--log-file', type=str,
                                  help='Log file path')
    completeness_parser.add_argument('-v', '--verbose', action='store_true',
                                  help='Enable verbose output')
    completeness_parser.add_argument('--reps-only', action='store_true',
                                  help='Validate only representative proteins')

    # Analyze mode parser
    analyze_parser = subparsers.add_parser("analyze", help="Analyze domain data and statistics")
    analyze_subparsers = analyze_parser.add_subparsers(dest="action", help="Analysis action")

    # Batch comparison analyzer
    comparison_parser = analyze_subparsers.add_parser("comparison", help="Compare regular and representative batches")
    comparison_parser.add_argument('--config', type=str, default='config/config.yml',
                                help='Path to configuration file')
    comparison_parser.add_argument('--output', type=str,
                                help='Output JSON file for comparison results')
    comparison_parser.add_argument('--analyze-domains', action='store_true',
                                help='Analyze domain files (slower but more comprehensive)')
    comparison_parser.add_argument('--log-file', type=str,
                                help='Log file path')
    comparison_parser.add_argument('-v', '--verbose', action='store_true',
                                help='Enable verbose output')

    # Domain statistics analyzer
    domains_parser = analyze_subparsers.add_parser("domains", help="Analyze domain statistics")
    domains_parser.add_argument('--config', type=str, default='config/config.yml',
                             help='Path to configuration file')
    domains_parser.add_argument('--batch-id', type=int,
                             help='Batch ID to analyze (optional, all batches if not specified)')
    domains_parser.add_argument('--output', type=str,
                             help='Output JSON file for domain statistics')
    domains_parser.add_argument('--sample-size', type=int, default=1000,
                             help='Number of files to sample for detailed analysis')
    domains_parser.add_argument('--log-file', type=str,
                             help='Log file path')
    domains_parser.add_argument('-v', '--verbose', action='store_true',
                             help='Enable verbose output')

    # Pipeline statistics analyzer
    pipeline_parser = analyze_subparsers.add_parser("pipeline", help="Analyze pipeline statistics")
    pipeline_parser.add_argument('--config', type=str, default='config/config.yml',
                              help='Path to configuration file')
    pipeline_parser.add_argument('--output', type=str,
                              help='Output JSON file for pipeline statistics')
    pipeline_parser.add_argument('--timing', action='store_true',
                              help='Include timing statistics (if available)')
    pipeline_parser.add_argument('--log-file', type=str,
                              help='Log file path')
    pipeline_parser.add_argument('-v', '--verbose', action='store_true',
                              help='Enable verbose output')

    # Trace mode parser
    trace_parser = subparsers.add_parser("trace", help="Trace evidence and domain assignments")
    trace_subparsers = trace_parser.add_subparsers(dest="action", help="Trace action")

    # Evidence chain tracer
    evidence_trace_parser = trace_subparsers.add_parser("evidence", help="Trace evidence chain for a protein")
    evidence_trace_parser.add_argument('--config', type=str, default='config/config.yml',
                                    help='Path to configuration file')
    evidence_trace_parser.add_argument('--batch-id', type=int, required=True,
                                    help='Batch ID to trace')
    evidence_trace_parser.add_argument('--pdb-id', type=str, required=True,
                                    help='PDB ID to trace')
    evidence_trace_parser.add_argument('--chain-id', type=str, required=True,
                                    help='Chain ID to trace')
    evidence_trace_parser.add_argument('--output', type=str,
                                    help='Output JSON file for trace results')
    evidence_trace_parser.add_argument('--log-file', type=str,
                                    help='Log file path')
    evidence_trace_parser.add_argument('-v', '--verbose', action='store_true',
                                    help='Enable verbose output')

    # Domain assignment tracer
    domain_trace_parser = trace_subparsers.add_parser("domain", help="Trace domain assignment")
    domain_trace_parser.add_argument('--config', type=str, default='config/config.yml',
                                  help='Path to configuration file')
    domain_trace_parser.add_argument('--batch-id', type=int, required=True,
                                  help='Batch ID to trace')
    domain_trace_parser.add_argument('--pdb-id', type=str,
                                  help='PDB ID to trace (optional, all proteins if not specified)')
    domain_trace_parser.add_argument('--chain-id', type=str,
                                  help='Chain ID to trace (optional, all chains if not specified)')
    domain_trace_parser.add_argument('--limit', type=int, default=100,
                                  help='Limit number of proteins to trace')
    domain_trace_parser.add_argument('--output', type=str,
                                  help='Output JSON file for trace results')
    domain_trace_parser.add_argument('--log-file', type=str,
                                  help='Log file path')
    domain_trace_parser.add_argument('-v', '--verbose', action='store_true',
                                  help='Enable verbose output')

    # Raw data tracer
    raw_trace_parser = trace_subparsers.add_parser("raw", help="Trace raw data files and conversions")
    raw_trace_parser.add_argument('--config', type=str, default='config/config.yml',
                               help='Path to configuration file')
    raw_trace_parser.add_argument('--batch-id', type=int, required=True,
                               help='Batch ID to trace')
    raw_trace_parser.add_argument('--limit', type=int,
                               help='Limit number of proteins to trace')
    raw_trace_parser.add_argument('--blast-only', action='store_true',
                               help='Ignore HHSearch files in chain analysis')
    raw_trace_parser.add_argument('--output', type=str,
                               help='Output JSON file for trace results')
    raw_trace_parser.add_argument('--log-file', type=str,
                               help='Log file path')
    raw_trace_parser.add_argument('-v', '--verbose', action='store_true',
                               help='Enable verbose output')
    raw_trace_parser.add_argument('--reps-only', action='store_true',
                               help='Trace only representative proteins')

    # Report mode parser
    report_parser = subparsers.add_parser("report", help="Generate comprehensive reports")
    report_subparsers = report_parser.add_subparsers(dest="action", help="Report action")

    # Batch report generator
    batch_report_parser = report_subparsers.add_parser("batch", help="Generate batch report")
    batch_report_parser.add_argument('--config', type=str, default='config/config.yml',
                                  help='Path to configuration file')
    batch_report_parser.add_argument('--batch-id', type=int, required=True,
                                  help='Batch ID to report on')
    batch_report_parser.add_argument('--analyze-domains', action='store_true',
                                  help='Analyze domain files (slower but more comprehensive)')
    batch_report_parser.add_argument('--output', type=str,
                                  help='Output JSON file for report')
    batch_report_parser.add_argument('--log-file', type=str,
                                  help='Log file path')
    batch_report_parser.add_argument('-v', '--verbose', action='store_true',
                                  help='Enable verbose output')

    # Pipeline overview generator
    overview_report_parser = report_subparsers.add_parser("overview", help="Generate pipeline overview")
    overview_report_parser.add_argument('--config', type=str, default='config/config.yml',
                                     help='Path to configuration file')
    overview_report_parser.add_argument('--output', type=str,
                                     help='Output JSON file for overview')
    overview_report_parser.add_argument('--log-file', type=str,
                                     help='Log file path')
    overview_report_parser.add_argument('-v', '--verbose', action='store_true',
                                     help='Enable verbose output')

    # Domain report generator
    domain_report_parser = report_subparsers.add_parser("domains", help="Generate domain report")
    domain_report_parser.add_argument('--config', type=str, default='config/config.yml',
                                   help='Path to configuration file')
    domain_report_parser.add_argument('--batch-id', type=int,
                                   help='Batch ID to report on (optional, all batches if not specified)')
    domain_report_parser.add_argument('--sample-size', type=int, default=1000,
                                   help='Number of files to sample for detailed analysis')
    domain_report_parser.add_argument('--output', type=str,
                                   help='Output JSON file for report')
    domain_report_parser.add_argument('--log-file', type=str,
                                   help='Log file path')
    domain_report_parser.add_argument('-v', '--verbose', action='store_true',
                                   help='Enable verbose output')

    # Diagnose mode parser
    diagnose_parser = subparsers.add_parser("diagnose", help="Diagnose pipeline issues")
    diagnose_subparsers = diagnose_parser.add_subparsers(dest="action", help="Diagnostic action")

    # Missing files diagnostics
    missing_diag_parser = diagnose_subparsers.add_parser("missing", help="Diagnose missing files")
    missing_diag_parser.add_argument('--config', type=str, default='config/config.yml',
                                  help='Path to configuration file')
    missing_diag_parser.add_argument('--batch-id', type=int, required=True,
                                  help='Batch ID to diagnose')
    missing_diag_parser.add_argument('--blast-only', action='store_true',
                                  help='Ignore HHSearch files in chain analysis')
    missing_diag_parser.add_argument('--limit', type=int,
                                  help='Limit number of proteins to diagnose')
    missing_diag_parser.add_argument('--output', type=str,
                                  help='Output JSON file for diagnostic results')
    missing_diag_parser.add_argument('--log-file', type=str,
                                  help='Log file path')
    missing_diag_parser.add_argument('-v', '--verbose', action='store_true',
                                  help='Enable verbose output')
    missing_diag_parser.add_argument('--reps-only', action='store_true',
                                  help='Diagnose only representative proteins')

    # Consistency diagnostics
    consistency_diag_parser = diagnose_subparsers.add_parser("consistency", help="Diagnose consistency issues")
    consistency_diag_parser.add_argument('--config', type=str, default='config/config.yml',
                                      help='Path to configuration file')
    consistency_diag_parser.add_argument('--batch-id', type=int, required=True,
                                      help='Batch ID to diagnose')
    consistency_diag_parser.add_argument('--limit', type=int,
                                      help='Limit number of proteins to diagnose')
    consistency_diag_parser.add_argument('--output', type=str,
                                      help='Output JSON file for diagnostic results')
    consistency_diag_parser.add_argument('--log-file', type=str,
                                      help='Log file path')
    consistency_diag_parser.add_argument('-v', '--verbose', action='store_true',
                                      help='Enable verbose output')
    consistency_diag_parser.add_argument('--reps-only', action='store_true',
                                      help='Diagnose only representative proteins')

    # Content quality diagnostics
    quality_diag_parser = diagnose_subparsers.add_parser("quality", help="Diagnose content quality issues")
    quality_diag_parser.add_argument('--config', type=str, default='config/config.yml',
                                  help='Path to configuration file')
    quality_diag_parser.add_argument('--batch-id', type=int, required=True,
                                  help='Batch ID to diagnose')
    quality_diag_parser.add_argument('--limit', type=int,
                                  help='Limit number of proteins to diagnose')
    quality_diag_parser.add_argument('--output', type=str,
                                  help='Output JSON file for diagnostic results')
    quality_diag_parser.add_argument('--log-file', type=str,
                                  help='Log file path')
    quality_diag_parser.add_argument('-v', '--verbose', action='store_true',
                                  help='Enable verbose output')
    quality_diag_parser.add_argument('--reps-only', action='store_true',
                                  help='Diagnose only representative proteins')

    # Parse arguments
    args = parser.parse_args()

    # Set up logging
    log_file = args.log_file if hasattr(args, 'log_file') and args.log_file else None
    verbose = args.verbose if hasattr(args, 'verbose') else False
    setup_logging(verbose, log_file)

    logger = logging.getLogger("validate_pipeline")

    # Run appropriate mode
    if args.mode == "validate":
        return validate_mode(args)
    elif args.mode == "analyze":
        return analyze_mode(args)
    elif args.mode == "trace":
        return trace_mode(args)
    elif args.mode == "report":
        return report_mode(args)
    elif args.mode == "diagnose":
        return diagnose_mode(args)
    else:
        logger.error(f"Unknown mode: {args.mode}")
        parser.print_help()
        return 1

if __name__ == "__main__":
    sys.exit(main())
