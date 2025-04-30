#!/usr/bin/env python3
"""
domain_partition_tools.py - Unified Domain Partition toolkit

This script provides a comprehensive set of tools for working with domain partition:
- partition: Convert domain summaries to domain partition files
- analyze: Examine and diagnose domain coverage and quality
- repair: Fix missing or problematic domain files
- batches: Process multiple batches in parallel
- validate: Validate summary inputs before partition

Each mode has specific options and can work with different backends (database or filesystem).
"""

import os
import sys
import logging
import argparse
import glob
import random
import re
import json
import xml.etree.ElementTree as ET
from xml.dom import minidom
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple, Union
from concurrent.futures import ThreadPoolExecutor, as_completed

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

# Import core modules
from ecod.core.context import ApplicationContext
from ecod.config import ConfigManager

# Import processing modules
from ecod.pipelines.domain_analysis.partition import DomainPartition
from ecod.models.pipeline import DomainSummaryModel
from ecod.models import BlastHit, HHSearchHit

# Import path utilities
from ecod.utils.path_utils import (
    get_standardized_paths,
    get_file_db_path,
    resolve_file_path,
    find_files_with_legacy_paths,
    migrate_file_to_standard_path,
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
# PARTITION MODE FUNCTIONS
#

def partition_all_summaries(args: argparse.Namespace) -> int:
    """
    Partition all domain summaries in a batch
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("ecod.partition.batch")
    
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
    
    logger.info(f"Processing summaries for batch {args.batch_id} ({batch_info[0]['batch_name']})")
    
    # Initialize partition component
    partition = DomainPartition(context)
    
    # Process for specific representative IDs if provided
    if hasattr(args, 'process_ids') and args.process_ids:
        logger.info(f"Processing specific process IDs: {len(args.process_ids)}")
        result = partition.process_specific_ids(
            args.batch_id, 
            args.process_ids, 
            batch_path, 
            reference, 
            args.blast_only
        )
        
        return 0 if result else 1
    
    # Otherwise process full batch
    result = partition.process_batch(
        args.batch_id, 
        batch_path, 
        reference, 
        args.blast_only,
        args.limit,
        args.reps_only
    )
    
    return 0 if result else 1

def partition_via_filesystem(args: argparse.Namespace) -> int:
    """
    Partition domain summaries using filesystem instead of database
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("ecod.partition.fs")
    
    batch_path = args.batch_path
    reference = args.ref_version
    
    logger.info(f"Processing domain summaries from filesystem at {batch_path}")
    
    # Find domain summary files
    summary_pattern = os.path.join(batch_path, "domains", f"*.{reference}.domain_summary*.xml")
    summary_files = glob.glob(summary_pattern)
    
    if not summary_files:
        logger.warning(f"No domain summary files found matching pattern: {summary_pattern}")
        
        # Try alternative patterns
        alt_patterns = [
            os.path.join(batch_path, "domains", f"*.domain_summary.xml"),
            os.path.join(batch_path, "domains", f"*.domain_summary.{reference}.xml"),
            os.path.join(batch_path, "domains", f"*.{reference}.domains.xml"),
        ]
        
        for pattern in alt_patterns:
            alt_files = glob.glob(pattern)
            if alt_files:
                logger.info(f"Found {len(alt_files)} summary files with alternative pattern: {pattern}")
                summary_files = alt_files
                break
                
        if not summary_files:
            logger.error("No domain summary files found with any pattern")
            return 1
    
    logger.info(f"Found {len(summary_files)} domain summary files")
    
    if args.limit:
        summary_files = summary_files[:args.limit]
        logger.info(f"Limited processing to {args.limit} files")
    
    # Initialize application context for configuration
    context = ApplicationContext(args.config) if args.config else ApplicationContext()
    
    # Process domain partition files
    success_count = 0
    error_count = 0
    
    # Initialize domain partition with context
    partition = DomainPartition(context)
    
    # Set force overwrite if needed
    if args.force:
        context.set_force_overwrite(True)
    
    # Process each summary file
    for summary_file in summary_files:
        # Extract PDB and chain ID from filename
        filename = os.path.basename(summary_file)
        match = re.match(r'([^_]+)_([^.]+)', filename)
        
        if not match:
            logger.warning(f"Could not parse PDB and chain ID from filename: {filename}")
            error_count += 1
            continue
            
        pdb_id, chain_id = match.groups()
        
        try:
            logger.info(f"Processing {pdb_id}_{chain_id}")
            
            # Run partition
            domain_file = partition.partition_domains(
                pdb_id,
                chain_id,
                batch_path,
                'struct_seqid',  # Default input mode
                reference,
                args.blast_only
            )
            
            if domain_file:
                success_count += 1
                logger.info(f"Created domain file: {domain_file}")
            else:
                error_count += 1
                logger.warning(f"Failed to create domain file for {pdb_id}_{chain_id}")
                
        except Exception as e:
            error_count += 1
            logger.error(f"Error processing {pdb_id}_{chain_id}: {str(e)}")
    
    logger.info(f"Processed {len(summary_files)} summary files: {success_count} succeeded, {error_count} failed")
    
    return 0 if success_count > 0 else 1

def partition_mode(args: argparse.Namespace) -> int:
    """
    Run partition mode with appropriate backend
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("ecod.partition")
    
    # Choose between database or filesystem backend
    if args.backend == "db":
        return partition_all_summaries(args)
    elif args.backend == "fs":
        return partition_via_filesystem(args)
    else:
        logger.error(f"Unknown backend: {args.backend}")
        return 1

#
# ANALYZE MODE FUNCTIONS
#

def analyze_domain_distribution(args: argparse.Namespace) -> int:
    """
    Analyze domain distribution and statistics for a batch
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("ecod.analyze.distribution")
    
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
    
    logger.info(f"Analyzing domain distribution for batch {args.batch_id} ({batch_info[0]['batch_name']})")
    
    # Find all domain partition files
    domain_files = []
    
    # Option 1: Use database to find files
    query = """
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
        AND pf.file_type = 'domain_partition'
        AND pf.file_exists = TRUE
    """
    
    rows = context.db.execute_dict_query(query, (args.batch_id,))
    
    for row in rows:
        file_path = os.path.join(batch_path, row['file_path'])
        if os.path.exists(file_path):
            domain_files.append({
                'pdb_id': row['pdb_id'],
                'chain_id': row['chain_id'],
                'path': file_path
            })
    
    # Option 2: If no files found in database, scan filesystem
    if not domain_files:
        logger.info("No domain files found in database, scanning filesystem")
        
        domains_dir = os.path.join(batch_path, "domains")
        if os.path.exists(domains_dir):
            pattern = os.path.join(domains_dir, f"*.{reference}.domains*.xml")
            for file_path in glob.glob(pattern):
                filename = os.path.basename(file_path)
                match = re.match(r'([^_]+)_([^.]+)', filename)
                
                if match:
                    pdb_id, chain_id = match.groups()
                    domain_files.append({
                        'pdb_id': pdb_id,
                        'chain_id': chain_id,
                        'path': file_path
                    })
    
    if not domain_files:
        logger.error("No domain files found")
        return 1
        
    logger.info(f"Found {len(domain_files)} domain files")
    
    # Initialize counters and statistics
    stats = {
        'total_chains': len(domain_files),
        'chains_with_domains': 0,
        'total_domains': 0,
        'domain_counts': {},  # Number of chains with N domains
        'single_domain': 0,
        'multi_domain': 0,
        'unclassified': 0,
        'discontinuous': 0,
        'domain_sizes': [],  # List of domain sizes for histogram
        't_groups': {},  # T-group distribution
        'h_groups': {},  # H-group distribution
        'coverage': []  # Chain coverage percentages
    }
    
    # Process each domain file
    for domain_file in domain_files:
        try:
            # Parse XML
            tree = ET.parse(domain_file['path'])
            root = tree.getroot()
            
            # Check if unclassified
            if root.get('status') == 'unclassified':
                stats['unclassified'] += 1
                continue
                
            # Extract domains
            domains = root.findall('./domain_list/domain')
            domain_count = len(domains)
            
            # Update statistics
            stats['total_domains'] += domain_count
            
            if domain_count > 0:
                stats['chains_with_domains'] += 1
                
                # Track domain count distribution
                if domain_count not in stats['domain_counts']:
                    stats['domain_counts'][domain_count] = 0
                stats['domain_counts'][domain_count] += 1
                
                # Single vs multi-domain
                if domain_count == 1:
                    stats['single_domain'] += 1
                else:
                    stats['multi_domain'] += 1
                
                # Process each domain
                for domain in domains:
                    # Check for discontinuous domains
                    domain_range = domain.get('range', '')
                    if ',' in domain_range:
                        stats['discontinuous'] += 1
                    
                    # Extract domain size
                    try:
                        # Parse range like "1-100" or "1-50,60-100"
                        size = 0
                        for segment in domain_range.split(','):
                            start, end = map(int, segment.split('-'))
                            size += (end - start + 1)
                        stats['domain_sizes'].append(size)
                    except (ValueError, AttributeError):
                        pass
                    
                    # Extract T-group and H-group
                    t_group = domain.get('t_group', '')
                    if t_group:
                        if t_group not in stats['t_groups']:
                            stats['t_groups'][t_group] = 0
                        stats['t_groups'][t_group] += 1
                    
                    h_group = domain.get('h_group', '')
                    if h_group:
                        if h_group not in stats['h_groups']:
                            stats['h_groups'][h_group] = 0
                        stats['h_groups'][h_group] += 1
            
            # Extract coverage
            coverage_elem = root.find('./statistics/coverage')
            if coverage_elem is not None:
                try:
                    coverage = float(coverage_elem.text)
                    stats['coverage'].append(coverage * 100)  # Convert to percentage
                except (ValueError, TypeError):
                    pass
            
        except Exception as e:
            logger.error(f"Error processing domain file {domain_file['path']}: {str(e)}")
    
    # Calculate summary statistics
    if stats['chains_with_domains'] > 0:
        stats['avg_domains_per_chain'] = stats['total_domains'] / stats['chains_with_domains']
    else:
        stats['avg_domains_per_chain'] = 0
    
    if stats['domain_sizes']:
        stats['avg_domain_size'] = sum(stats['domain_sizes']) / len(stats['domain_sizes'])
        stats['min_domain_size'] = min(stats['domain_sizes'])
        stats['max_domain_size'] = max(stats['domain_sizes'])
    else:
        stats['avg_domain_size'] = 0
        stats['min_domain_size'] = 0
        stats['max_domain_size'] = 0
    
    if stats['coverage']:
        stats['avg_coverage'] = sum(stats['coverage']) / len(stats['coverage'])
        stats['min_coverage'] = min(stats['coverage'])
        stats['max_coverage'] = max(stats['coverage'])
    else:
        stats['avg_coverage'] = 0
        stats['min_coverage'] = 0
        stats['max_coverage'] = 0
    
    # Sort T-groups and H-groups by frequency
    stats['top_t_groups'] = sorted(stats['t_groups'].items(), key=lambda x: x[1], reverse=True)[:10]
    stats['top_h_groups'] = sorted(stats['h_groups'].items(), key=lambda x: x[1], reverse=True)[:10]
    
    # Print summary
    logger.info(f"Domain distribution summary:")
    logger.info(f"Total chains analyzed: {stats['total_chains']}")
    logger.info(f"Chains with domains: {stats['chains_with_domains']} ({stats['chains_with_domains']/stats['total_chains']*100:.1f}%)")
    logger.info(f"Unclassified chains: {stats['unclassified']} ({stats['unclassified']/stats['total_chains']*100:.1f}%)")
    logger.info(f"Single-domain chains: {stats['single_domain']} ({stats['single_domain']/stats['chains_with_domains']*100:.1f}% of classified)")
    logger.info(f"Multi-domain chains: {stats['multi_domain']} ({stats['multi_domain']/stats['chains_with_domains']*100:.1f}% of classified)")
    logger.info(f"Total domains: {stats['total_domains']}")
    logger.info(f"Average domains per chain: {stats['avg_domains_per_chain']:.2f}")
    logger.info(f"Discontinuous domains: {stats['discontinuous']} ({stats['discontinuous']/stats['total_domains']*100:.1f}%)")
    logger.info(f"Domain size: avg={stats['avg_domain_size']:.1f}, min={stats['min_domain_size']}, max={stats['max_domain_size']}")
    logger.info(f"Chain coverage: avg={stats['avg_coverage']:.1f}%, min={stats['min_coverage']:.1f}%, max={stats['max_coverage']:.1f}%")
    
    # Print domain count distribution
    logger.info("Domain count distribution:")
    for count, num_chains in sorted(stats['domain_counts'].items()):
        logger.info(f"  {count} domain(s): {num_chains} chains ({num_chains/stats['chains_with_domains']*100:.1f}%)")
    
    # Print top T-groups
    logger.info("Top 10 T-groups:")
    for t_group, count in stats['top_t_groups']:
        logger.info(f"  {t_group}: {count} domains ({count/stats['total_domains']*100:.1f}%)")
    
    # Print top H-groups
    logger.info("Top 10 H-groups:")
    for h_group, count in stats['top_h_groups']:
        logger.info(f"  {h_group}: {count} domains ({count/stats['total_domains']*100:.1f}%)")
    
    # Write statistics to JSON file if output file is specified
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(stats, f, indent=2)
        logger.info(f"Statistics written to {args.output}")
    
    return 0

def analyze_coverage_quality(args: argparse.Namespace) -> int:
    """
    Analyze domain coverage quality for a batch
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("ecod.analyze.coverage")
    
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
    
    logger.info(f"Analyzing domain coverage quality for batch {args.batch_id} ({batch_info[0]['batch_name']})")
    
    # Get statistics from process_status
    query = """
    SELECT 
        COUNT(*) as total_chains,
        SUM(CASE WHEN current_stage = 'domain_partition_complete' THEN 1 ELSE 0 END) as completed,
        SUM(CASE WHEN current_stage = 'domain_partition_failed' THEN 1 ELSE 0 END) as failed,
        SUM(CASE WHEN status = 'error' THEN 1 ELSE 0 END) as error_count
    FROM 
        ecod_schema.process_status
    WHERE 
        batch_id = %s
    """
    
    status_stats = context.db.execute_dict_query(query, (args.batch_id,))
    
    if status_stats:
        logger.info(f"Process status statistics:")
        logger.info(f"Total chains: {status_stats[0]['total_chains']}")
        logger.info(f"Completed: {status_stats[0]['completed']} ({status_stats[0]['completed']/status_stats[0]['total_chains']*100:.1f}%)")
        logger.info(f"Failed: {status_stats[0]['failed']} ({status_stats[0]['failed']/status_stats[0]['total_chains']*100:.1f}%)")
        logger.info(f"Error: {status_stats[0]['error_count']} ({status_stats[0]['error_count']/status_stats[0]['total_chains']*100:.1f}%)")
    
    # Get error messages
    if args.errors:
        error_query = """
        SELECT 
            p.pdb_id, p.chain_id, ps.error_message
        FROM 
            ecod_schema.process_status ps
        JOIN 
            ecod_schema.protein p ON ps.protein_id = p.id
        WHERE 
            ps.batch_id = %s
            AND ps.status = 'error'
        LIMIT 50
        """
        
        error_rows = context.db.execute_dict_query(error_query, (args.batch_id,))
        
        if error_rows:
            logger.info(f"Sample error messages ({len(error_rows)} of {status_stats[0]['error_count']}):")
            for row in error_rows:
                logger.info(f"  {row['pdb_id']}_{row['chain_id']}: {row['error_message']}")
    
    # Get file statistics
    file_query = """
    SELECT 
        file_type, 
        COUNT(*) as count,
        SUM(CASE WHEN file_exists = TRUE THEN 1 ELSE 0 END) as exist_count
    FROM 
        ecod_schema.process_file pf
    JOIN 
        ecod_schema.process_status ps ON pf.process_id = ps.id
    WHERE 
        ps.batch_id = %s
    GROUP BY 
        file_type
    """
    
    file_stats = context.db.execute_dict_query(file_query, (args.batch_id,))
    
    if file_stats:
        logger.info(f"File statistics:")
        for row in file_stats:
            logger.info(f"  {row['file_type']}: {row['exist_count']}/{row['count']} ({row['exist_count']/row['count']*100:.1f}%)")
    
    # Return success
    return 0

def analyze_mode(args: argparse.Namespace) -> int:
    """
    Run analyze mode with appropriate action
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("ecod.analyze")
    
    if args.action == "distribution":
        return analyze_domain_distribution(args)
    elif args.action == "coverage":
        return analyze_coverage_quality(args)
    else:
        logger.error(f"Unknown analysis action: {args.action}")
        return 1

#
# REPAIR MODE FUNCTIONS
#

def find_missing_partitions(args: argparse.Namespace) -> int:
    """
    Find missing domain partition files
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("ecod.repair.missing")
    
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
    
    logger.info(f"Finding missing domain partitions for batch {args.batch_id} ({batch_info[0]['batch_name']})")
    
    # Find proteins with domain summary but no domain partition
    query = """
    SELECT 
        p.id as protein_id, p.pdb_id, p.chain_id, ps.id as process_id,
        (SELECT pf.file_path FROM ecod_schema.process_file pf 
         WHERE pf.process_id = ps.id AND pf.file_type = 'domain_summary' AND pf.file_exists = TRUE
         LIMIT 1) as summary_path
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
        AND NOT EXISTS (
            SELECT 1 FROM ecod_schema.process_file pf 
            WHERE pf.process_id = ps.id AND pf.file_type = 'domain_partition' AND pf.file_exists = TRUE
        )
    """
    
    if args.reps_only:
        query += " AND ps.is_representative = TRUE"
    
    rows = context.db.execute_dict_query(query, (args.batch_id,))
    
    if not rows:
        logger.info("No missing domain partitions found")
        return 0
        
    logger.info(f"Found {len(rows)} proteins with missing domain partition files")
    
    # Summarize findings
    if args.summary_only:
        logger.info("Summary of missing partition files:")
        for i, row in enumerate(rows[:10]):  # Show first 10
            summary_path = os.path.join(batch_path, row['summary_path']) if row['summary_path'] else "Not found"
            logger.info(f"  {i+1}. {row['pdb_id']}_{row['chain_id']}: Summary at {summary_path}")
        
        if len(rows) > 10:
            logger.info(f"  ... and {len(rows) - 10} more")
        
        return 0
    
    # Initialize domain partition
    partition = DomainPartition(context)
    
    # Set force overwrite if needed
    if args.force:
        context.set_force_overwrite(True)
    
    # Process missing partitions
    success_count = 0
    error_count = 0
    
    # Build list of process IDs
    process_ids = [row['process_id'] for row in rows]
    
    # Apply limit if specified
    if args.limit and args.limit < len(process_ids):
        process_ids = process_ids[:args.limit]
        logger.info(f"Limited to processing {args.limit} proteins")
    
    # Process files
    result = partition.process_specific_ids(
        args.batch_id, 
        process_ids, 
        batch_path, 
        reference, 
        args.blast_only
    )
    
    return 0 if result else 1

def repair_db_file_paths(args: argparse.Namespace) -> int:
    """
    Repair database file paths for domain partitions
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("ecod.repair.db_paths")
    
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
    
    logger.info(f"Repairing database file paths for batch {args.batch_id} ({batch_info[0]['batch_name']})")
    
    # Find all domain partition files in the database
    query = """
    SELECT 
        pf.id, pf.process_id, pf.file_path, pf.file_exists,
        p.pdb_id, p.chain_id
    FROM 
        ecod_schema.process_file pf
    JOIN 
        ecod_schema.process_status ps ON pf.process_id = ps.id
    JOIN 
        ecod_schema.protein p ON ps.protein_id = p.id
    WHERE 
        ps.batch_id = %s
        AND pf.file_type = 'domain_partition'
    """
    
    rows = context.db.execute_dict_query(query, (args.batch_id,))
    
    if not rows:
        logger.info("No domain partition files found in database")
        return 0
        
    logger.info(f"Found {len(rows)} domain partition file records in database")
    
    # Initialize counters
    updated_count = 0
    error_count = 0
    
    # Process each file record
    for row in rows:
        file_id = row['id']
        pdb_id = row['pdb_id']
        chain_id = row['chain_id']
        current_path = row['file_path']
        file_exists = row['file_exists']
        
        # Calculate standard path
        standard_dir = os.path.join(batch_path, "domains")
        standard_filename = f"{pdb_id}_{chain_id}.{reference}.domains.xml"
        standard_path = os.path.join(standard_dir, standard_filename)
        
        # Convert to relative path for database
        relative_path = os.path.relpath(standard_path, batch_path)
        
        # Check if path needs updating
        if current_path != relative_path:
            logger.info(f"Path for {pdb_id}_{chain_id} needs updating:")
            logger.info(f"  Current: {current_path}")
            logger.info(f"  Standard: {relative_path}")
            
            # If file exists at current path but not at standard path, move it
            current_full_path = os.path.join(batch_path, current_path)
            
            if os.path.exists(current_full_path) and not os.path.exists(standard_path) and not args.dry_run:
                try:
                    # Create directory if needed
                    os.makedirs(os.path.dirname(standard_path), exist_ok=True)
                    # Copy file
                    import shutil
                    shutil.copy2(current_full_path, standard_path)
                    logger.info(f"  Copied file to standard location")
                except Exception as e:
                    logger.error(f"  Error copying file: {str(e)}")
            
            # Update database record
            if not args.dry_run:
                try:
                    # Check if file exists at standard path
                    file_exists_at_standard = os.path.exists(standard_path)
                    
                    context.db.update(
                        "ecod_schema.process_file",
                        {
                            "file_path": relative_path,
                            "file_exists": file_exists_at_standard,
                            "file_size": os.path.getsize(standard_path) if file_exists_at_standard else 0
                        },
                        "id = %s",
                        (file_id,)
                    )
                    
                    updated_count += 1
                    logger.info(f"  Updated database record")
                except Exception as e:
                    error_count += 1
                    logger.error(f"  Error updating database record: {str(e)}")
            else:
                logger.info(f"  Would update database record (dry run)")
                updated_count += 1
        
    logger.info(f"Repair summary: {updated_count} records updated, {error_count} errors")
    
    return 0

def repair_mode(args: argparse.Namespace) -> int:
    """
    Run repair mode with appropriate action
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("ecod.repair")
    
    if args.action == "missing":
        return find_missing_partitions(args)
    elif args.action == "db_paths":
        return repair_db_file_paths(args)
    else:
        logger.error(f"Unknown repair action: {args.action}")
        return 1

#
# VALIDATE MODE FUNCTIONS
#

def validate_summaries(args: argparse.Namespace) -> int:
    """
    Validate domain summary files for a batch
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("ecod.validate")
    
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
    
    # Find domains ready for partition
    query = """
    SELECT 
        p.id as protein_id, p.pdb_id, p.chain_id, ps.id as process_id,
        (SELECT pf.file_path FROM ecod_schema.process_file pf 
         WHERE pf.process_id = ps.id AND pf.file_type = 'domain_summary' AND pf.file_exists = TRUE
         LIMIT 1) as summary_path
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
        logger.error("No domain summaries found")
        return 1
        
    logger.info(f"Found {len(rows)} proteins with domain summaries")
    
    # Apply limit if specified
    if args.limit and args.limit < len(rows):
        rows = rows[:args.limit]
        logger.info(f"Limited to validating {args.limit} proteins")
    
    # Validate each summary
    valid_count = 0
    error_count = 0
    
    # Initialize statistics
    stats = {
        'total': len(rows),
        'valid': 0,
        'invalid': 0,
        'error_types': {},
        'blast_only': 0,
        'with_blast': 0,
        'with_hhsearch': 0,
        'with_domains': 0,
        'empty': 0
    }
    
    for row in rows:
        pdb_id = row['pdb_id']
        chain_id = row['chain_id']
        summary_path = row['summary_path']
        
        # Build full path
        if not os.path.isabs(summary_path):
            full_path = os.path.join(batch_path, summary_path)
        else:
            full_path = summary_path
        
        logger.debug(f"Validating {pdb_id}_{chain_id} summary at {full_path}")
        
        try:
            # Check if file exists
            if not os.path.exists(full_path):
                logger.error(f"Summary file not found: {full_path}")
                error_count += 1
                stats['invalid'] += 1
                
                if 'file_not_found' not in stats['error_types']:
                    stats['error_types']['file_not_found'] = 0
                stats['error_types']['file_not_found'] += 1
                
                continue
            
            # Parse XML
            tree = ET.parse(full_path)
            root = tree.getroot()
            
            # Basic validation of root element
            if root.tag != "domain_summ_doc":
                logger.error(f"Invalid root element: {root.tag}")
                error_count += 1
                stats['invalid'] += 1
                
                if 'invalid_root' not in stats['error_types']:
                    stats['error_types']['invalid_root'] = 0
                stats['error_types']['invalid_root'] += 1
                
                continue
            
            # Check for the main component tags expected in a summary
            if root.find("chain_blast_evidence") is None:
                if 'missing_chain_blast' not in stats['error_types']:
                    stats['error_types']['missing_chain_blast'] = 0
                stats['error_types']['missing_chain_blast'] += 1
            else:
                # Check if it has hits
                chain_blast_hits = root.findall(".//chain_blast_evidence/*/hit")
                if chain_blast_hits:
                    stats['with_blast'] += 1
            
            if root.find("domain_blast_evidence") is None:
                if 'missing_domain_blast' not in stats['error_types']:
                    stats['error_types']['missing_domain_blast'] = 0
                stats['error_types']['missing_domain_blast'] += 1
            
            if root.find("hhsearch_evidence") is None:
                # Might be blast-only summary
                if ".blast_only" in full_path:
                    stats['blast_only'] += 1
                else:
                    if 'missing_hhsearch' not in stats['error_types']:
                        stats['error_types']['missing_hhsearch'] = 0
                    stats['error_types']['missing_hhsearch'] += 1
            else:
                # Check if it has hits
                hhsearch_hits = root.findall(".//hhsearch_evidence/*/hh_hit")
                if hhsearch_hits:
                    stats['with_hhsearch'] += 1
            
            # Check for domain suggestions
            domain_suggestions = root.find("domain_suggestions")
            if domain_suggestions is not None:
                domains = domain_suggestions.findall("./domain")
                if domains:
                    stats['with_domains'] += 1
                else:
                    stats['empty'] += 1
            else:
                if 'missing_domain_suggestions' not in stats['error_types']:
                    stats['error_types']['missing_domain_suggestions'] = 0
                stats['error_types']['missing_domain_suggestions'] += 1
            
            # If we got here, summary is at least structurally valid
            valid_count += 1
            stats['valid'] += 1
            
        except ET.ParseError as e:
            logger.error(f"XML parsing error for {pdb_id}_{chain_id}: {str(e)}")
            error_count += 1
            stats['invalid'] += 1
            
            if 'parse_error' not in stats['error_types']:
                stats['error_types']['parse_error'] = 0
            stats['error_types']['parse_error'] += 1
            
        except Exception as e:
            logger.error(f"Error validating {pdb_id}_{chain_id}: {str(e)}")
            error_count += 1
            stats['invalid'] += 1
            
            if 'other_error' not in stats['error_types']:
                stats['error_types']['other_error'] = 0
            stats['error_types']['other_error'] += 1
    
    # Print summary statistics
    logger.info(f"Validation summary:")
    logger.info(f"Total summaries: {stats['total']}")
    logger.info(f"Valid: {stats['valid']} ({stats['valid']/stats['total']*100:.1f}%)")
    logger.info(f"Invalid: {stats['invalid']} ({stats['invalid']/stats['total']*100:.1f}%)")
    logger.info(f"BLAST-only: {stats['blast_only']} ({stats['blast_only']/stats['total']*100:.1f}%)")
    logger.info(f"With BLAST hits: {stats['with_blast']} ({stats['with_blast']/stats['total']*100:.1f}%)")
    logger.info(f"With HHSearch hits: {stats['with_hhsearch']} ({stats['with_hhsearch']/stats['total']*100:.1f}%)")
    logger.info(f"With domain suggestions: {stats['with_domains']} ({stats['with_domains']/stats['total']*100:.1f}%)")
    logger.info(f"Empty domain suggestions: {stats['empty']} ({stats['empty']/stats['total']*100:.1f}%)")
    
    if stats['error_types']:
        logger.info(f"Error types:")
        for error_type, count in stats['error_types'].items():
            logger.info(f"  {error_type}: {count} ({count/stats['total']*100:.1f}%)")
    
    # Write statistics to JSON file if output file is specified
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(stats, f, indent=2)
        logger.info(f"Statistics written to {args.output}")
    
    # Return success if more than half are valid
    return 0 if valid_count >= len(rows) / 2 else 1

def validate_mode(args: argparse.Namespace) -> int:
    """
    Run validate mode
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("ecod.validate")
    
    return validate_summaries(args)

#
# BATCHES MODE FUNCTIONS
#

def process_batch(batch_id: int, config_path: str, blast_only: bool = False, force: bool = False) -> Tuple[int, Dict[str, Any]]:
    """
    Process a single batch
    
    Args:
        batch_id: Batch ID to process
        config_path: Path to configuration file
        blast_only: Whether to use only BLAST results (no HHSearch)
        force: Force reprocessing of already processed results
        
    Returns:
        Tuple of (batch_id, stats dictionary)
    """
    logger = logging.getLogger(f"ecod.batches.{batch_id}")
    
    # Initialize application context
    context = ApplicationContext(config_path)
    
    # Set force overwrite if needed
    if force:
        context.set_force_overwrite(True)
    
    # Get batch info
    batch_info = context.db.execute_dict_query(
        "SELECT id, batch_name, base_path, ref_version FROM ecod_schema.batch WHERE id = %s", 
        (batch_id,)
    )
    
    if not batch_info:
        logger.error(f"Batch {batch_id} not found")
        return batch_id, {"error": "Batch not found"}
        
    batch_path = batch_info[0]['base_path']
    reference = batch_info[0]['ref_version']
    
    logger.info(f"Processing batch {batch_id} ({batch_info[0]['batch_name']})")
    
    # Initialize domain partition
    partition = DomainPartition(context)
    
    # Process batch
    try:
        start_time = datetime.now()
        
        # Process batch
        result = partition.process_batch(
            batch_id, 
            batch_path, 
            reference, 
            blast_only
        )
        
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        return batch_id, {
            "success": result,
            "duration": duration,
            "blast_only": blast_only
        }
    
    except Exception as e:
        logger.error(f"Error processing batch {batch_id}: {str(e)}")
        return batch_id, {
            "success": False,
            "error": str(e),
            "blast_only": blast_only
        }

def process_multiple_batches(args: argparse.Namespace) -> int:
    """
    Process multiple batches in parallel
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("ecod.batches")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get batch IDs to process
    if args.batch_ids:
        # Use specified batch IDs
        batch_ids = [int(b_id) for b_id in args.batch_ids]
        logger.info(f"Processing specified batches: {batch_ids}")
    else:
        # Get all batch IDs from database
        query = """
        SELECT 
            id, batch_name
        FROM 
            ecod_schema.batch
        ORDER BY 
            id
        """
        
        rows = context.db.execute_dict_query(query)
        
        if not rows:
            logger.error("No batches found in database")
            return 1
            
        batch_ids = [row['id'] for row in rows]
        logger.info(f"Found {len(batch_ids)} batches in database")
    
    # Filter out excluded batch IDs
    if args.exclude_batch_ids:
        exclude_ids = [int(b_id) for b_id in args.exclude_batch_ids]
        batch_ids = [b_id for b_id in batch_ids if b_id not in exclude_ids]
        logger.info(f"After exclusions: {len(batch_ids)} batches to process")
    
    # Process batches in parallel
    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        future_to_batch = {
            executor.submit(
                process_batch, 
                batch_id, 
                args.config, 
                args.blast_only, 
                args.force
            ): batch_id for batch_id in batch_ids
        }
        
        # Collect results
        results = {}
        for future in as_completed(future_to_batch):
            batch_id = future_to_batch[future]
            try:
                batch_id, stats = future.result()
                results[batch_id] = stats
                
                if stats.get("success", False):
                    logger.info(f"Batch {batch_id} completed successfully in {stats.get('duration', 0):.1f} seconds")
                else:
                    logger.error(f"Batch {batch_id} failed: {stats.get('error', 'Unknown error')}")
            except Exception as e:
                logger.error(f"Error processing batch {batch_id}: {str(e)}")
                results[batch_id] = {"success": False, "error": str(e)}
    
    # Count successful and failed batches
    success_count = sum(1 for stats in results.values() if stats.get("success", False))
    failed_count = len(results) - success_count
    
    logger.info(f"Processing summary: {success_count} batches succeeded, {failed_count} failed")
    
    # Write detailed results to output file if specified
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)
        logger.info(f"Detailed results written to {args.output}")
    
    return 0 if success_count > 0 else 1

def batches_mode(args: argparse.Namespace) -> int:
    """
    Run batches mode
    
    Args:
        args: Command line arguments
        
    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("ecod.batches")
    
    return process_multiple_batches(args)

#
# MAIN FUNCTION AND ARGUMENT PARSING
#

def main():
    """Main entry point"""
    # Create top-level parser
    parser = argparse.ArgumentParser(description='Domain Partition Tools - Unified domain partition toolkit')
    subparsers = parser.add_subparsers(dest="mode", help="Operating mode")

    # Partition mode parser
    partition_parser = subparsers.add_parser("partition", help="Run domain partition on summaries")
    partition_subparsers = partition_parser.add_subparsers(dest="backend", help="Backend to use")

    # Database backend subparser
    db_parser = partition_subparsers.add_parser("db", help="Use database backend")
    db_parser.add_argument('--config', type=str, default='config/config.yml',
                       help='Path to configuration file')
    db_parser.add_argument('--batch-id', type=int, required=True,
                       help='Batch ID to process')
    db_parser.add_argument('--process-ids', type=int, nargs='+',
                       help='Specific process IDs to process')
    db_parser.add_argument('--blast-only', action='store_true',
                       help='Use only BLAST results (no HHSearch)')
    db_parser.add_argument('--limit', type=int,
                       help='Limit number of proteins to process')
    db_parser.add_argument('--log-file', type=str,
                       help='Log file path')
    db_parser.add_argument('-v', '--verbose', action='store_true',
                       help='Enable verbose output')
    db_parser.add_argument('--force', action='store_true',
                       help='Force regeneration of domain files')
    db_parser.add_argument('--reps-only', action='store_true',
                       help='Process only representative proteins')

    # Filesystem backend subparser
    fs_parser = partition_subparsers.add_parser("fs", help="Use filesystem backend")
    fs_parser.add_argument('--config', type=str,
                       help='Path to configuration file')
    fs_parser.add_argument('--batch-path', type=str, required=True,
                       help='Path to batch directory')
    fs_parser.add_argument('--ref-version', type=str, default="develop291",
                       help='Reference version')
    fs_parser.add_argument('--blast-only', action='store_true',
                       help='Use only BLAST results (no HHSearch)')
    fs_parser.add_argument('--limit', type=int,
                       help='Limit number of proteins to process')
    fs_parser.add_argument('--log-file', type=str,
                       help='Log file path')
    fs_parser.add_argument('-v', '--verbose', action='store_true',
                       help='Enable verbose output')
    fs_parser.add_argument('--force', action='store_true',
                       help='Force regeneration of domain files')

    # Analyze mode parser
    analyze_parser = subparsers.add_parser("analyze", help="Analyze domain coverage and quality")
    analyze_subparsers = analyze_parser.add_subparsers(dest="action", help="Analysis action")

    # Domain distribution subparser
    distribution_parser = analyze_subparsers.add_parser("distribution", help="Analyze domain distribution")
    distribution_parser.add_argument('--config', type=str, default='config/config.yml',
                                 help='Path to configuration file')
    distribution_parser.add_argument('--batch-id', type=int, required=True,
                                 help='Batch ID to analyze')
    distribution_parser.add_argument('--output', type=str,
                                 help='Output JSON file for statistics')
    distribution_parser.add_argument('--log-file', type=str,
                                 help='Log file path')
    distribution_parser.add_argument('-v', '--verbose', action='store_true',
                                 help='Enable verbose output')

    # Coverage quality subparser
    coverage_parser = analyze_subparsers.add_parser("coverage", help="Analyze domain coverage quality")
    coverage_parser.add_argument('--config', type=str, default='config/config.yml',
                             help='Path to configuration file')
    coverage_parser.add_argument('--batch-id', type=int, required=True,
                             help='Batch ID to analyze')
    coverage_parser.add_argument('--errors', action='store_true',
                             help='Show error messages')
    coverage_parser.add_argument('--log-file', type=str,
                             help='Log file path')
    coverage_parser.add_argument('-v', '--verbose', action='store_true',
                             help='Enable verbose output')

    # Repair mode parser
    repair_parser = subparsers.add_parser("repair", help="Fix missing or problematic domain files")
    repair_subparsers = repair_parser.add_subparsers(dest="action", help="Repair action")

    # Find missing partitions subparser
    missing_parser = repair_subparsers.add_parser("missing", help="Find missing domain partition files")
    missing_parser.add_argument('--config', type=str, default='config/config.yml',
                            help='Path to configuration file')
    missing_parser.add_argument('--batch-id', type=int, required=True,
                            help='Batch ID to process')
    missing_parser.add_argument('--summary-only', action='store_true',
                            help='Only show summary, don\'t repair')
    missing_parser.add_argument('--blast-only', action='store_true',
                            help='Use only BLAST results (no HHSearch)')
    missing_parser.add_argument('--limit', type=int,
                            help='Limit number of proteins to process')
    missing_parser.add_argument('--log-file', type=str,
                            help='Log file path')
    missing_parser.add_argument('-v', '--verbose', action='store_true',
                            help='Enable verbose output')
    missing_parser.add_argument('--force', action='store_true',
                            help='Force regeneration of domain files')
    missing_parser.add_argument('--reps-only', action='store_true',
                            help='Process only representative proteins')

    # Repair database file paths subparser
    db_paths_parser = repair_subparsers.add_parser("db_paths", help="Repair database file paths")
    db_paths_parser.add_argument('--config', type=str, default='config/config.yml',
                              help='Path to configuration file')
    db_paths_parser.add_argument('--batch-id', type=int, required=True,
                              help='Batch ID to process')
    db_paths_parser.add_argument('--dry-run', action='store_true',
                              help='Don\'t actually update anything')
    db_paths_parser.add_argument('--log-file', type=str,
                              help='Log file path')
    db_paths_parser.add_argument('-v', '--verbose', action='store_true',
                              help='Enable verbose output')

    # Validate mode parser
    validate_parser = subparsers.add_parser("validate", help="Validate domain summaries")
    validate_parser.add_argument('--config', type=str, default='config/config.yml',
                             help='Path to configuration file')
    validate_parser.add_argument('--batch-id', type=int, required=True,
                             help='Batch ID to validate')
    validate_parser.add_argument('--limit', type=int,
                             help='Limit number of proteins to validate')
    validate_parser.add_argument('--output', type=str,
                             help='Output JSON file for statistics')
    validate_parser.add_argument('--log-file', type=str,
                             help='Log file path')
    validate_parser.add_argument('-v', '--verbose', action='store_true',
                             help='Enable verbose output')
    validate_parser.add_argument('--reps-only', action='store_true',
                             help='Validate only representative proteins')

    # Batches mode parser
    batches_parser = subparsers.add_parser("batches", help="Process multiple batches in parallel")
    batches_parser.add_argument('--config', type=str, default='config/config.yml',
                            help='Path to configuration file')
    batches_parser.add_argument('--batch-ids', type=int, nargs='+',
                            help='Specific batch IDs to process')
    batches_parser.add_argument('--exclude-batch-ids', type=int, nargs='+',
                            help='Batch IDs to exclude')
    batches_parser.add_argument('--blast-only', action='store_true',
                            help='Use only BLAST results (no HHSearch)')
    batches_parser.add_argument('--force', action='store_true',
                            help='Force regeneration of domain files')
    batches_parser.add_argument('--max-workers', type=int, default=4,
                            help='Maximum number of worker threads')
    batches_parser.add_argument('--output', type=str,
                            help='Output JSON file for detailed results')
    batches_parser.add_argument('--log-file', type=str,
                            help='Log file path')
    batches_parser.add_argument('-v', '--verbose', action='store_true',
                            help='Enable verbose output')

    # Parse arguments
    args = parser.parse_args()

    # Set up logging
    log_file = args.log_file if hasattr(args, 'log_file') and args.log_file else None
    verbose = args.verbose if hasattr(args, 'verbose') else False
    setup_logging(verbose, log_file)

    logger = logging.getLogger("domain_partition_tools")

    # Run appropriate mode
    if args.mode == "partition":
        return partition_mode(args)
    elif args.mode == "analyze":
        return analyze_mode(args)
    elif args.mode == "repair":
        return repair_mode(args)
    elif args.mode == "validate":
        return validate_mode(args)
    elif args.mode == "batches":
        return batches_mode(args)
    else:
        logger.error(f"Unknown mode: {args.mode}")
        parser.print_help()
        return 1

if __name__ == "__main__":
    sys.exit(main())
