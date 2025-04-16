#!/usr/bin/env python3
"""
batch_domain_analysis.py - Analyze domain generation results for a specific batch

This script analyzes the domain files from a batch to determine:
a) Quality of domain generation
b) Proteins with no domains and their blast hit status
c) Breakdown between chain and domain blast usage
d) Candidates for full HHsearch pipelines
"""

import os
import sys
import argparse
import json
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
from typing import Dict, List, Tuple, Any, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("domain_analysis")

def load_config(config_path: str) -> Dict:
    """Load ECOD configuration file
    
    This function tries to load a configuration file from the specified path.
    It also looks for a local config file that can override settings.
    """
    import yaml
    import os
    from pathlib import Path
    
    config = {}
    
    # Load main config
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
            logger.info(f"Loaded configuration from {config_path}")
    except Exception as e:
        logger.error(f"Error loading config from {config_path}: {str(e)}")
        sys.exit(1)
    
    # Check for local config file (same filename with .local suffix)
    config_dir = os.path.dirname(config_path)
    config_name = os.path.basename(config_path)
    local_config_path = os.path.join(config_dir, f"{os.path.splitext(config_name)[0]}.local.yml")
    
    if os.path.exists(local_config_path):
        try:
            with open(local_config_path, 'r') as f:
                local_config = yaml.safe_load(f)
                if local_config:
                    # Merge configs, with local taking precedence
                    logger.info(f"Merging settings from local config: {local_config_path}")
                    deep_update(config, local_config)
        except Exception as e:
            logger.warning(f"Error loading local config from {local_config_path}: {str(e)}")
    
    # Look for personal config file in user home directory
    user_config_path = Path.home() / ".ecod" / "config.yml"
    if user_config_path.exists():
        try:
            with open(user_config_path, 'r') as f:
                user_config = yaml.safe_load(f)
                if user_config:
                    # Merge configs, with user config taking precedence
                    logger.info(f"Merging settings from user config: {user_config_path}")
                    deep_update(config, user_config)
        except Exception as e:
            logger.warning(f"Error loading user config from {user_config_path}: {str(e)}")
            
    return config

def deep_update(original: Dict, update: Dict) -> Dict:
    """Deep update a nested dictionary
    
    Updates original dictionary with values from update dictionary.
    For nested dictionaries, updates keys recursively rather than replacing the dict.
    """
    for key, value in update.items():
        if key in original and isinstance(original[key], dict) and isinstance(value, dict):
            deep_update(original[key], value)
        else:
            original[key] = value
    return original

def connect_to_db(config: Dict, prompt_password: bool = False) -> Any:
    """Connect to the database using configuration
    
    Args:
        config: Configuration dictionary
        prompt_password: If True, prompt for password instead of using config
    """
    try:
        import psycopg2
        import psycopg2.extras
        import getpass
        from pathlib import Path
        import os
        
        # Get database connection parameters from config
        db_config = config.get('database', {})
        host = db_config.get('host', 'localhost')
        port = db_config.get('port', 5432)
        database = db_config.get('database', 'ecod')
        user = db_config.get('user', 'postgres')
        
        # Get password (3 options)
        password = None
        
        # Option 1: Check environment variable
        if 'ECOD_DB_PASSWORD' in os.environ:
            password = os.environ.get('ECOD_DB_PASSWORD')
            logger.debug("Using database password from environment variable")
        
        # Option 2: Check .pgpass file
        if password is None and not prompt_password:
            pgpass_path = Path.home() / '.pgpass'
            if pgpass_path.exists():
                try:
                    with open(pgpass_path, 'r') as f:
                        for line in f:
                            parts = line.strip().split(':')
                            if len(parts) == 5:
                                pg_host, pg_port, pg_db, pg_user, pg_pass = parts
                                if (pg_host == '*' or pg_host == host) and \
                                   (pg_port == '*' or pg_port == str(port)) and \
                                   (pg_db == '*' or pg_db == database) and \
                                   (pg_user == '*' or pg_user == user):
                                    password = pg_pass
                                    logger.debug("Using database password from .pgpass file")
                                    break
                except Exception as e:
                    logger.debug(f"Error reading .pgpass file: {str(e)}")
        
        # Option 3: Get from config if not found yet and not prompting
        if password is None and not prompt_password:
            password = db_config.get('password', '')
            logger.debug("Using database password from config file")
        
        # Option 4: Prompt user if requested or all else fails
        if password is None or prompt_password:
            password = getpass.getpass(f"Enter password for PostgreSQL user {user}: ")
            logger.debug("Using database password from user prompt")
        
        # Connect to database
        conn = psycopg2.connect(
            host=host,
            port=port,
            database=database,
            user=user,
            password=password
        )
        
        return conn
    except ImportError:
        logger.error("psycopg2 not installed. Please install it with: pip install psycopg2-binary")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Database connection error: {str(e)}")
        sys.exit(1)

def get_batch_info(conn, batch_id: int) -> Dict:
    """Get information about the specified batch"""
    query = """
    SELECT 
        b.id, 
        b.base_path, 
        b.ref_version,
        COUNT(ps.id) AS total_proteins
    FROM 
        ecod_schema.batch b
    LEFT JOIN 
        ecod_schema.process_status ps ON b.id = ps.batch_id
    WHERE 
        b.id = %s
    GROUP BY 
        b.id, b.base_path, b.ref_version
    """
    
    cursor = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    cursor.execute(query, (batch_id,))
    result = cursor.fetchone()
    
    if not result:
        logger.error(f"Batch {batch_id} not found")
        sys.exit(1)
    
    return dict(result)

def get_process_info(conn, batch_id: int) -> List[Dict]:
    """Get process information for all proteins in the batch"""
    query = """
    SELECT 
        ps.id AS process_id,
        p.id AS protein_id, 
        p.pdb_id, 
        p.chain_id,
        ps.status
    FROM 
        ecod_schema.process_status ps
    JOIN 
        ecod_schema.protein p ON ps.protein_id = p.id
    WHERE 
        ps.batch_id = %s
    ORDER BY 
        p.pdb_id, p.chain_id
    """
    
    cursor = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    cursor.execute(query, (batch_id,))
    return [dict(row) for row in cursor.fetchall()]

def get_file_paths(conn, batch_id: int) -> Dict[int, Dict[str, str]]:
    """Get file paths for all processes in the batch"""
    query = """
    SELECT 
        ps.id AS process_id,
        pf.file_type,
        pf.file_path,
        pf.file_exists
    FROM 
        ecod_schema.process_status ps
    JOIN 
        ecod_schema.process_file pf ON ps.id = pf.process_id
    WHERE 
        ps.batch_id = %s
    """
    
    cursor = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    cursor.execute(query, (batch_id,))
    
    file_paths = defaultdict(dict)
    for row in cursor.fetchall():
        process_id = row['process_id']
        file_type = row['file_type']
        file_paths[process_id][file_type] = {
            'path': row['file_path'],
            'exists': row['file_exists']
        }
    
    return file_paths

def load_domain_file(file_path: str) -> Dict:
    """Load domain JSON file"""
    try:
        with open(file_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        logger.warning(f"Error loading domain file {file_path}: {str(e)}")
        return {}

def load_blast_file(file_path: str) -> Dict:
    """Load BLAST results JSON file"""
    try:
        with open(file_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        logger.warning(f"Error loading BLAST file {file_path}: {str(e)}")
        return {}

def analyze_domains(domain_data: Dict) -> Dict:
    """Analyze domain data and return metrics"""
    if not domain_data:
        return {'domain_count': 0, 'source': None}
    
    domains = domain_data.get('domains', [])
    
    # Check source of domains (chain vs domain blast)
    source_counts = Counter()
    for domain in domains:
        source = domain.get('source', 'unknown')
        source_counts[source] += 1
    
    # Determine primary source
    primary_source = None
    if source_counts:
        primary_source = source_counts.most_common(1)[0][0]
    
    return {
        'domain_count': len(domains),
        'source': primary_source,
        'source_counts': dict(source_counts)
    }

def analyze_blast_hits(blast_data: Dict) -> Dict:
    """Analyze BLAST results and return metrics"""
    if not blast_data:
        return {'hit_count': 0, 'has_hits': False}
    
    hits = blast_data.get('hits', [])
    
    return {
        'hit_count': len(hits),
        'has_hits': len(hits) > 0
    }

def main():
    parser = argparse.ArgumentParser(description='Analyze domain generation results for a batch')
    parser.add_argument('--config', type=str, default='config/config.yml', 
                        help='Path to ECOD configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                        help='Batch ID to analyze')
    parser.add_argument('--output-dir', type=str, default='analysis_results',
                        help='Directory for output files')
    parser.add_argument('--prompt-password', action='store_true',
                        help='Prompt for database password instead of using config')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load configuration and connect to database
    config = load_config(args.config)
    conn = connect_to_db(config, prompt_password=args.prompt_password)
    
    try:
        # Get batch information
        batch_info = get_batch_info(conn, args.batch_id)
        base_path = batch_info['base_path']
        logger.info(f"Analyzing batch {args.batch_id} with {batch_info['total_proteins']} proteins")
        
        # Get process information and file paths
        processes = get_process_info(conn, args.batch_id)
        file_paths = get_file_paths(conn, args.batch_id)
        
        # Initialize results storage
        results = []
        
        # Process each protein
        for proc in processes:
            process_id = proc['process_id']
            pdb_id = proc['pdb_id']
            chain_id = proc['chain_id']
            protein_key = f"{pdb_id}_{chain_id}"
            
            # Skip if process failed
            if proc['status'] != 'success':
                logger.warning(f"Skipping {protein_key} (process {process_id}): Status is {proc['status']}")
                continue
            
            # Check for domain file
            domain_file_info = file_paths.get(process_id, {}).get('domain_file', {})
            domain_path = None
            if domain_file_info and domain_file_info.get('exists', False):
                domain_path = os.path.join(base_path, domain_file_info.get('path', ''))
            
            # Check for blast files
            chain_blast_info = file_paths.get(process_id, {}).get('chain_blast_result', {})
            domain_blast_info = file_paths.get(process_id, {}).get('domain_blast_result', {})
            
            chain_blast_path = None
            if chain_blast_info and chain_blast_info.get('exists', False):
                chain_blast_path = os.path.join(base_path, chain_blast_info.get('path', ''))
            
            domain_blast_path = None
            if domain_blast_info and domain_blast_info.get('exists', False):
                domain_blast_path = os.path.join(base_path, domain_blast_info.get('path', ''))
            
            # Process domain file
            domain_data = {}
            domain_analysis = {'domain_count': 0, 'source': None}
            if domain_path and os.path.exists(domain_path):
                domain_data = load_domain_file(domain_path)
                domain_analysis = analyze_domains(domain_data)
            
            # Process BLAST files
            chain_blast_data = {}
            chain_blast_analysis = {'hit_count': 0, 'has_hits': False}
            if chain_blast_path and os.path.exists(chain_blast_path):
                chain_blast_data = load_blast_file(chain_blast_path)
                chain_blast_analysis = analyze_blast_hits(chain_blast_data)
            
            domain_blast_data = {}
            domain_blast_analysis = {'hit_count': 0, 'has_hits': False}
            if domain_blast_path and os.path.exists(domain_blast_path):
                domain_blast_data = load_blast_file(domain_blast_path)
                domain_blast_analysis = analyze_blast_hits(domain_blast_data)
            
            # Determine HHsearch candidacy
            # Criteria: Has BLAST hits but no domains or few domains
            has_blast_hits = chain_blast_analysis['has_hits'] or domain_blast_analysis['has_hits']
            few_domains = 0 <= domain_analysis['domain_count'] <= 1
            hhsearch_candidate = has_blast_hits and few_domains
            
            # Store results
            results.append({
                'process_id': process_id,
                'protein_id': proc['protein_id'],
                'pdb_id': pdb_id,
                'chain_id': chain_id,
                'protein_key': protein_key,
                
                # Domain analysis
                'domain_count': domain_analysis['domain_count'],
                'domain_source': domain_analysis['source'],
                'domain_source_counts': domain_analysis.get('source_counts', {}),
                
                # BLAST analysis
                'chain_blast_hits': chain_blast_analysis['hit_count'],
                'domain_blast_hits': domain_blast_analysis['hit_count'],
                'has_chain_blast_hits': chain_blast_analysis['has_hits'],
                'has_domain_blast_hits': domain_blast_analysis['has_hits'],
                
                # HHsearch candidacy
                'hhsearch_candidate': hhsearch_candidate
            })
        
        # Convert results to DataFrame for analysis
        df = pd.DataFrame(results)
        
        # Save full results
        results_path = os.path.join(args.output_dir, f'batch_{args.batch_id}_domain_analysis.csv')
        df.to_csv(results_path, index=False)
        
        # Summary statistics
        total_proteins = len(df)
        proteins_with_domains = df[df['domain_count'] > 0].shape[0]
        proteins_no_domains = df[df['domain_count'] == 0].shape[0]
        
        # Domain source breakdown
        domain_source_counts = Counter()
        for _, row in df.iterrows():
            if row['domain_count'] > 0 and row['domain_source']:
                domain_source_counts[row['domain_source']] += 1
        
        # Blast hit statistics
        proteins_with_chain_hits = df[df['has_chain_blast_hits']].shape[0]
        proteins_with_domain_hits = df[df['has_domain_blast_hits']].shape[0]
        proteins_with_any_hits = df[df['has_chain_blast_hits'] | df['has_domain_blast_hits']].shape[0]
        proteins_with_no_hits = df[~df['has_chain_blast_hits'] & ~df['has_domain_blast_hits']].shape[0]
        
        # No domains but has hits
        no_domains_with_hits = df[(df['domain_count'] == 0) & 
                                  (df['has_chain_blast_hits'] | df['has_domain_blast_hits'])].shape[0]
        
        # HHsearch candidates
        hhsearch_candidates = df[df['hhsearch_candidate']].shape[0]
        
        # Create summary report
        summary = {
            'batch_id': args.batch_id,
            'total_proteins': total_proteins,
            'proteins_with_domains': proteins_with_domains,
            'proteins_no_domains': proteins_no_domains,
            'domain_source_breakdown': dict(domain_source_counts),
            'proteins_with_chain_hits': proteins_with_chain_hits,
            'proteins_with_domain_hits': proteins_with_domain_hits,
            'proteins_with_any_hits': proteins_with_any_hits,
            'proteins_with_no_hits': proteins_with_no_hits,
            'no_domains_with_hits': no_domains_with_hits,
            'hhsearch_candidates': hhsearch_candidates
        }
        
        # Save summary
        summary_path = os.path.join(args.output_dir, f'batch_{args.batch_id}_summary.json')
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        
        # Print summary
        print("\n===== DOMAIN ANALYSIS SUMMARY =====")
        print(f"Batch ID: {args.batch_id}")
        print(f"Total proteins analyzed: {total_proteins}")
        print(f"Proteins with domains: {proteins_with_domains} ({100*proteins_with_domains/total_proteins:.1f}%)")
        print(f"Proteins without domains: {proteins_no_domains} ({100*proteins_no_domains/total_proteins:.1f}%)")
        
        print("\nDomain source breakdown:")
        for source, count in domain_source_counts.most_common():
            print(f"  - {source}: {count} ({100*count/proteins_with_domains:.1f}%)")
        
        print("\nBLAST hit statistics:")
        print(f"Proteins with chain BLAST hits: {proteins_with_chain_hits} ({100*proteins_with_chain_hits/total_proteins:.1f}%)")
        print(f"Proteins with domain BLAST hits: {proteins_with_domain_hits} ({100*proteins_with_domain_hits/total_proteins:.1f}%)")
        print(f"Proteins with any BLAST hits: {proteins_with_any_hits} ({100*proteins_with_any_hits/total_proteins:.1f}%)")
        print(f"Proteins with no BLAST hits: {proteins_with_no_hits} ({100*proteins_with_no_hits/total_proteins:.1f}%)")
        
        print("\nQuality metrics:")
        print(f"Proteins with no domains despite BLAST hits: {no_domains_with_hits} ({100*no_domains_with_hits/total_proteins:.1f}%)")
        print(f"Potential HHsearch candidates: {hhsearch_candidates} ({100*hhsearch_candidates/total_proteins:.1f}%)")
        
        # Generate HHsearch candidate list
        hhsearch_df = df[df['hhsearch_candidate']][['pdb_id', 'chain_id', 'protein_key', 'domain_count', 'chain_blast_hits', 'domain_blast_hits']]
        hhsearch_path = os.path.join(args.output_dir, f'batch_{args.batch_id}_hhsearch_candidates.csv')
        hhsearch_df.to_csv(hhsearch_path, index=False)
        
        # Generate visualizations
        create_visualizations(df, args.output_dir, args.batch_id)
        
        print(f"\nResults saved to {args.output_dir}/")
        
    finally:
        # Close database connection
        conn.close()

def create_visualizations(df: pd.DataFrame, output_dir: str, batch_id: int):
    """Create visualizations of the analysis results"""
    # Set style
    plt.style.use('ggplot')
    
    # 1. Domain count distribution
    plt.figure(figsize=(10, 6))
    bins = np.arange(0, df['domain_count'].max() + 2) - 0.5
    plt.hist(df['domain_count'], bins=bins, alpha=0.7)
    plt.xlabel('Number of Domains')
    plt.ylabel('Number of Proteins')
    plt.title('Distribution of Domain Counts')
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(output_dir, f'batch_{batch_id}_domain_dist.png'), dpi=300, bbox_inches='tight')
    
    # 2. Domain source pie chart
    domain_sources = df[df['domain_count'] > 0]['domain_source'].value_counts()
    if not domain_sources.empty:
        plt.figure(figsize=(8, 8))
        plt.pie(domain_sources, labels=domain_sources.index, autopct='%1.1f%%', startangle=90)
        plt.axis('equal')
        plt.title('Domain Sources')
        plt.savefig(os.path.join(output_dir, f'batch_{batch_id}_domain_sources.png'), dpi=300, bbox_inches='tight')
    
    # 3. BLAST hit analysis
    plt.figure(figsize=(10, 6))
    categories = ['Chain Hits', 'Domain Hits', 'Any Hits', 'No Hits']
    values = [
        df['has_chain_blast_hits'].sum(),
        df['has_domain_blast_hits'].sum(),
        (df['has_chain_blast_hits'] | df['has_domain_blast_hits']).sum(),
        (~df['has_chain_blast_hits'] & ~df['has_domain_blast_hits']).sum()
    ]
    
    plt.bar(categories, values, alpha=0.7)
    plt.xlabel('BLAST Hit Type')
    plt.ylabel('Number of Proteins')
    plt.title('BLAST Hit Distribution')
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(output_dir, f'batch_{batch_id}_blast_hits.png'), dpi=300, bbox_inches='tight')
    
    # 4. HHsearch candidate analysis
    plt.figure(figsize=(10, 6))
    categories = ['Has Domains', 'No Domains', 'No Domains but Hits', 'HHsearch Candidates']
    values = [
        (df['domain_count'] > 0).sum(),
        (df['domain_count'] == 0).sum(),
        ((df['domain_count'] == 0) & (df['has_chain_blast_hits'] | df['has_domain_blast_hits'])).sum(),
        df['hhsearch_candidate'].sum()
    ]
    
    plt.bar(categories, values, alpha=0.7)
    plt.xlabel('Category')
    plt.ylabel('Number of Proteins')
    plt.title('Domain Generation Quality Metrics')
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(output_dir, f'batch_{batch_id}_quality_metrics.png'), dpi=300, bbox_inches='tight')

if __name__ == "__main__":
    main()