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
import json
import argparse
import logging
import getpass
from pathlib import Path
from collections import Counter, defaultdict
from typing import Dict, List, Tuple, Any, Optional

# Try importing optional dependencies
try:
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    HAS_PLOTTING = True
except ImportError:
    HAS_PLOTTING = False
    
try:
    import psycopg2
    import psycopg2.extras
    HAS_PSYCOPG2 = True
except ImportError:
    HAS_PSYCOPG2 = False
    
try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("domain_analysis")


def check_dependencies():
    """Check if required dependencies are installed"""
    missing_deps = []
    
    if not HAS_PSYCOPG2:
        missing_deps.append("psycopg2-binary")
    
    if not HAS_YAML:
        missing_deps.append("pyyaml")
    
    if not HAS_PLOTTING:
        missing_deps.append("pandas matplotlib numpy")
    
    if missing_deps:
        logger.error(f"Missing required dependencies: {', '.join(missing_deps)}")
        logger.error("Please install them with: pip install " + " ".join(missing_deps))
        sys.exit(1)


def load_config(config_path: str) -> Dict:
    """Load ECOD configuration file
    
    This function tries to load a configuration file from the specified path.
    It also looks for a local config file that can override settings.
    
    Args:
        config_path: Path to the main configuration file
        
    Returns:
        Dictionary containing configuration settings
    """
    if not HAS_YAML:
        logger.error("PyYAML not installed. Please install it with: pip install pyyaml")
        sys.exit(1)
    
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
    
    Args:
        original: Original dictionary to update
        update: Dictionary with values to update original with
        
    Returns:
        Updated original dictionary
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
        
    Returns:
        Database connection object
    """
    if not HAS_PSYCOPG2:
        logger.error("psycopg2 not installed. Please install it with: pip install psycopg2-binary")
        sys.exit(1)
        
    try:
        # Get database connection parameters from config
        db_config = config.get('database', {})
        host = db_config.get('host', 'localhost')
        port = db_config.get('port', 5432)
        database = db_config.get('database', 'ecod')
        user = db_config.get('user', 'postgres')
        
        # Get password (multiple options)
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
    except Exception as e:
        logger.error(f"Database connection error: {str(e)}")
        sys.exit(1)


def get_batch_info(conn, batch_id: int) -> Dict:
    """Get information about the specified batch
    
    Args:
        conn: Database connection
        batch_id: Batch ID to retrieve information for
        
    Returns:
        Dictionary containing batch information
    """
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
    """Get process information for all proteins in the batch with additional error details
    
    Args:
        conn: Database connection
        batch_id: Batch ID to analyze
        
    Returns:
        List of dictionaries containing process information
    """
    query = """
    SELECT 
        ps.id AS process_id,
        p.id AS protein_id, 
        p.pdb_id, 
        p.chain_id,
        ps.status,
        ps.error_message AS status_message,
        ps.current_stage,
        ps.updated_at AS status_time,
        0 AS retry_count,  -- Using 0 as default since retry_count isn't in the schema
        p.length AS chain_length,
        p.source_id,
        ps.is_representative,
        ps.relative_path,
        ps.created_at,
        ps.updated_at
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


def get_file_paths(conn, batch_id: int) -> Dict[int, Dict[str, Dict]]:
    """Get file paths for all processes in the batch
    
    Args:
        conn: Database connection
        batch_id: Batch ID to analyze
        
    Returns:
        Dictionary mapping process IDs to file information
    """
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


def detect_file_type(file_path: str) -> str:
    """Detect file type based on content rather than extension
    
    Args:
        file_path: Path to the file
        
    Returns:
        'json', 'xml', 'empty', 'missing', or 'unknown'
    """
    try:
        # Check if file exists
        if not os.path.exists(file_path):
            return 'missing'
            
        # Check if file is empty
        if os.path.getsize(file_path) == 0:
            return 'empty'
        
        # Read first few bytes to determine type
        with open(file_path, 'r') as f:
            start = f.read(1024).strip()
        
        # Check for JSON
        if start.startswith('{') or start.startswith('['):
            return 'json'
            
        # Check for XML
        if start.startswith('<?xml') or start.startswith('<'):
            return 'xml'
            
        # Unknown type
        return 'unknown'
        
    except Exception as e:
        logger.warning(f"Error detecting file type for {file_path}: {str(e)}")
        return 'error'


def load_domain_file(file_path: str) -> Dict:
    """Load domain JSON file with better error handling
    
    Args:
        file_path: Path to the domain file
        
    Returns:
        Dictionary containing domain data or empty dict on error
    """
    try:
        # Check file type
        file_type = detect_file_type(file_path)
        
        if file_type == 'missing':
            logger.warning(f"Domain file does not exist: {file_path}")
            return {}
            
        if file_type == 'empty':
            logger.warning(f"Domain file is empty: {file_path}")
            return {}
            
        if file_type != 'json':
            logger.warning(f"Domain file is not JSON format (detected: {file_type}): {file_path}")
            return {}
        
        # Load JSON file
        with open(file_path, 'r') as f:
            return json.load(f)
            
    except json.JSONDecodeError as e:
        logger.warning(f"Error parsing JSON in domain file {file_path}: {str(e)}")
        return {}
        
    except Exception as e:
        logger.warning(f"Error loading domain file {file_path}: {str(e)}")
        return {}


def load_blast_file(file_path: str) -> Dict:
    """Load BLAST results file - handles both JSON and XML formats
    
    Args:
        file_path: Path to the BLAST results file
        
    Returns:
        Dictionary containing BLAST results or empty dict on error
    """
    try:
        # Check if file exists and has content
        if not os.path.exists(file_path):
            logger.warning(f"BLAST file does not exist: {file_path}")
            return {}
            
        # Check file size
        if os.path.getsize(file_path) == 0:
            logger.warning(f"BLAST file is empty: {file_path}")
            return {}
        
        # Check file type
        file_type = detect_file_type(file_path)
        
        # Handle XML files
        if file_type == 'xml':
            import xml.etree.ElementTree as ET
            
            try:
                # Try to parse as XML
                tree = ET.parse(file_path)
                root = tree.getroot()
                
                # Extract hits from XML
                hits = []
                
                # Handle BLAST XML format
                if root.tag == 'BlastOutput':
                    iterations = root.findall('.//Iteration')
                    for iteration in iterations:
                        for hit in iteration.findall('.//Hit'):
                            hit_data = {
                                'hit_id': hit.find('Hit_id').text if hit.find('Hit_id') is not None else '',
                                'hit_def': hit.find('Hit_def').text if hit.find('Hit_def') is not None else '',
                                'hit_len': int(hit.find('Hit_len').text) if hit.find('Hit_len') is not None else 0,
                                'hsps': []
                            }
                            
                            for hsp in hit.findall('.//Hsp'):
                                hsp_data = {
                                    'bit_score': float(hsp.find('Hsp_bit-score').text) if hsp.find('Hsp_bit-score') is not None else 0,
                                    'evalue': float(hsp.find('Hsp_evalue').text) if hsp.find('Hsp_evalue') is not None else 0,
                                    'query_from': int(hsp.find('Hsp_query-from').text) if hsp.find('Hsp_query-from') is not None else 0,
                                    'query_to': int(hsp.find('Hsp_query-to').text) if hsp.find('Hsp_query-to') is not None else 0,
                                    'hit_from': int(hsp.find('Hsp_hit-from').text) if hsp.find('Hsp_hit-from') is not None else 0,
                                    'hit_to': int(hsp.find('Hsp_hit-to').text) if hsp.find('Hsp_hit-to') is not None else 0,
                                    'identity': int(hsp.find('Hsp_identity').text) if hsp.find('Hsp_identity') is not None else 0,
                                    'align_len': int(hsp.find('Hsp_align-len').text) if hsp.find('Hsp_align-len') is not None else 0
                                }
                                hit_data['hsps'].append(hsp_data)
                            
                            if hit_data['hsps']:  # Only add hits with HSPs
                                hits.append(hit_data)
                
                return {'hits': hits}
                
            except ET.ParseError as e:
                logger.warning(f"Error parsing XML in BLAST file {file_path}: {str(e)}")
                return {}
        
        # Handle JSON files (default)
        elif file_type == 'json':
            with open(file_path, 'r') as f:
                data = json.load(f)
                return data
        else:
            logger.warning(f"Unrecognized file type for BLAST file: {file_path} (detected: {file_type})")
            return {}
            
    except json.JSONDecodeError as e:
        # Check if file might be XML with wrong extension
        try:
            with open(file_path, 'r') as f:
                content = f.read().strip()
            
            if content.startswith('<') and ('xml' in content.lower() or '<BlastOutput>' in content):
                logger.warning(f"File {file_path} appears to be XML but has wrong extension. Trying XML parser...")
                
                # Create temporary symlink with .xml extension
                temp_path = file_path + ".xml"
                try:
                    os.symlink(file_path, temp_path)
                    result = load_blast_file(temp_path)
                    os.unlink(temp_path)
                    return result
                except Exception:
                    if os.path.exists(temp_path):
                        os.unlink(temp_path)
        except Exception:
            pass
            
        logger.warning(f"Error parsing JSON in BLAST file {file_path}: {str(e)}")
        return {}
        
    except Exception as e:
        logger.warning(f"Error loading BLAST file {file_path}: {str(e)}")
        return {}


def analyze_domains(domain_data: Dict) -> Dict:
    """Analyze domain data and return metrics
    
    Args:
        domain_data: Dictionary containing domain data
        
    Returns:
        Dictionary with domain analysis metrics
    """
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
    """Analyze BLAST results and return metrics
    
    Args:
        blast_data: Dictionary containing BLAST results
        
    Returns:
        Dictionary with BLAST analysis metrics
    """
    if not blast_data:
        return {'hit_count': 0, 'has_hits': False, 'error': True}
    
    hits = blast_data.get('hits', [])
    
    return {
        'hit_count': len(hits),
        'has_hits': len(hits) > 0,
        'error': False
    }


def classify_error(process_info: Dict) -> str:
    """Classify the type of error based on process information
    
    Args:
        process_info: Dictionary containing process information
        
    Returns:
        Error classification as a string
    """
    status = process_info.get('status', '')
    status_message = process_info.get('status_message', '')
    chain_length = process_info.get('chain_length', 0)
    current_stage = process_info.get('current_stage', '')
    retry_count = process_info.get('retry_count', 0)
    
    # Check for missing/incomplete chains
    if chain_length is None or chain_length == 0:
        return "missing_chain"
    
    if chain_length < 10:
        return "short_chain"
    
    # Check for specific stage failures
    if current_stage == 'blast_search':
        return "blast_error"
    
    if current_stage == 'domain_generation':
        return "domain_error"
    
    if current_stage == 'classification':
        return "classification_error"
    
    # Check for timeout errors
    if status_message and ('timeout' in status_message.lower() or 'timed out' in status_message.lower()):
        return "timeout"
    
    # Check for BLAST errors
    if status_message and ('blast' in status_message.lower() or 'search failed' in status_message.lower()):
        return "blast_error"
    
    # Check for file errors
    if status_message and ('file' in status_message.lower() or 'not found' in status_message.lower() or 'missing' in status_message.lower()):
        return "file_error"
    
    # Check for parsing errors
    if status_message and ('parse' in status_message.lower() or 'parsing' in status_message.lower() or 'invalid' in status_message.lower()):
        return "parse_error"
    
    # Check for domain errors
    if status_message and 'domain' in status_message.lower():
        return "domain_error"
    
    # Check for retry count
    if retry_count >= 3:
        return "max_retries"
    
    # Default case
    return "unknown_error"


def analyze_error_patterns(processes: List[Dict]) -> Dict:
    """Analyze patterns in error processes
    
    Args:
        processes: List of process information dictionaries
        
    Returns:
        Dictionary containing error pattern analysis
    """
    # Extract error processes
    error_processes = [p for p in processes if p.get('status') == 'error']
    
    # Count by PDB ID
    pdb_error_counts = Counter()
    for proc in error_processes:
        pdb_error_counts[proc.get('pdb_id')] += 1
    
    # Find PDBs with high error rates
    high_error_pdbs = {}
    for pdb_id, count in pdb_error_counts.items():
        # Count total chains for this PDB
        total_chains = sum(1 for p in processes if p.get('pdb_id') == pdb_id)
        error_rate = count / total_chains if total_chains > 0 else 0
        
        if error_rate > 0.5 and count >= 3:  # More than 50% failure and at least 3 chains
            high_error_pdbs[pdb_id] = {
                'error_count': count,
                'total_chains': total_chains,
                'error_rate': error_rate
            }
    
    # Classify errors
    error_classifications = Counter()
    error_details = {}
    
    # Track non-representative chains
    non_representative_errors = []
    
    # Track errors by current stage
    stage_errors = Counter()
    
    for proc in error_processes:
        error_class = classify_error(proc)
        error_classifications[error_class] += 1
        
        # Count errors by processing stage
        current_stage = proc.get('current_stage', 'unknown')
        stage_errors[current_stage] += 1
        
        # Check if non-representative chain
        if proc.get('is_representative') is False:
            non_representative_errors.append({
                'pdb_id': proc.get('pdb_id'),
                'chain_id': proc.get('chain_id'),
                'process_id': proc.get('process_id')
            })
        
        # Store details for each error type
        if error_class not in error_details:
            error_details[error_class] = []
        
        error_details[error_class].append({
            'pdb_id': proc.get('pdb_id'),
            'chain_id': proc.get('chain_id'),
            'process_id': proc.get('process_id'),
            'status_message': proc.get('status_message'),
            'chain_length': proc.get('chain_length'),
            'current_stage': proc.get('current_stage', ''),
            'is_representative': proc.get('is_representative', False),
            'source_id': proc.get('source_id', '')
        })
    
    # Check if errors are concentrated in specific chains types (by source_id pattern)
    source_patterns = Counter()
    for proc in error_processes:
        source_id = proc.get('source_id', '')
        # Extract pattern from source_id (e.g., PDB, AF, etc.)
        if source_id:
            pattern = source_id.split('_')[0] if '_' in source_id else source_id
            source_patterns[pattern] += 1
    
    return {
        'total_errors': len(error_processes),
        'error_classifications': dict(error_classifications),
        'high_error_pdbs': high_error_pdbs,
        'error_details': error_details,
        'non_representative_errors': non_representative_errors,
        'stage_errors': dict(stage_errors),
        'source_pattern_errors': dict(source_patterns)
    }


def check_table_exists(conn, table_name: str, schema: str = 'ecod_schema') -> bool:
    """Check if a table exists in the database
    
    Args:
        conn: Database connection
        table_name: Name of the table to check
        schema: Schema name (default: ecod_schema)
        
    Returns:
        True if table exists, False otherwise
    """
    try:
        query = """
        SELECT EXISTS (
            SELECT FROM information_schema.tables 
            WHERE table_schema = %s
            AND table_name = %s
        )
        """
        
        cursor = conn.cursor()
        cursor.execute(query, (schema, table_name))
        result = cursor.fetchone()[0]
        return result
    except Exception as e:
        logger.error(f"Error checking if table exists: {str(e)}")
        return False

def get_process_log_summary(conn, process_id: int) -> str:
    """Get summary of log entries for a specific process
    
    This function first checks if the process_log table exists.
    If not, it returns a message about using error_message from process_status instead.
    
    Args:
        conn: Database connection
        process_id: Process ID to get logs for
        
    Returns:
        Summary of log entries as a string
    """
    try:
        # Check if process_log table exists
        if not check_table_exists(conn, 'process_log'):
            # If table doesn't exist, get error_message from process_status
            query = """
            SELECT 
                error_message,
                updated_at,
                current_stage
            FROM 
                ecod_schema.process_status
            WHERE 
                id = %s
            """
            
            cursor = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
            cursor.execute(query, (process_id,))
            result = cursor.fetchone()
            
            if not result or not result['error_message']:
                return "No error message found in process_status"
            
            time_str = result['updated_at'].strftime('%Y-%m-%d %H:%M:%S')
            message = result['error_message'].replace('\n', ' ').strip()
            stage = result['current_stage']
            
            return f"{time_str} [ERROR] Stage: {stage}, Message: {message}"
        
        # If process_log table exists, use it
        query = """
        SELECT 
            log_level, 
            message,
            created_at
        FROM 
            ecod_schema.process_log
        WHERE 
            process_id = %s
        ORDER BY 
            created_at DESC
        LIMIT 10
        """
        
        cursor = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
        cursor.execute(query, (process_id,))
        logs = cursor.fetchall()
        
        if not logs:
            return "No log entries found"
        
        summary = []
        for log in logs:
            time_str = log['created_at'].strftime('%Y-%m-%d %H:%M:%S')
            level = log['log_level']
            message = log['message'].replace('\n', ' ').strip()
            if len(message) > 100:
                message = message[:97] + '...'
            
            summary.append(f"{time_str} [{level}] {message}")
        
        return '\n'.join(summary)
    except Exception as e:
        return f"Error retrieving logs: {str(e)}"


def create_visualizations(df, output_dir: str, batch_id: int):
    """Create visualizations of the analysis results
    
    Args:
        df: DataFrame containing analysis results
        output_dir: Directory to save visualizations
        batch_id: Batch ID for file naming
    """
    if not HAS_PLOTTING:
        logger.warning("Matplotlib and/or pandas not available. Skipping visualizations.")
        return
    
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
    plt.close()
    
    # 2. Domain source pie chart
    domain_sources = df[df['domain_count'] > 0]['domain_source'].value_counts()
    if not domain_sources.empty:
        plt.figure(figsize=(8, 8))
        plt.pie(domain_sources, labels=domain_sources.index, autopct='%1.1f%%', startangle=90)
        plt.axis('equal')
        plt.title('Domain Sources')
        plt.savefig(os.path.join(output_dir, f'batch_{batch_id}_domain_sources.png'), dpi=300, bbox_inches='tight')
        plt.close()
    
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
    plt.close()
    
    # 4. HHsearch candidate analysis
    plt.figure(figsize=(10, 6))
    categories = ['Has Domains', 'No Domains', 'No Domains but Hits', 'HHsearch Candidates']
    
    # Filter for valid BLAST files
    df_valid = df[~df['blast_file_errors']] if 'blast_file_errors' in df.columns else df
    
    values = [
        (df_valid['domain_count'] > 0).sum(),
        (df_valid['domain_count'] == 0).sum(),
        ((df_valid['domain_count'] == 0) & 
         (df_valid['has_chain_blast_hits'] | df_valid['has_domain_blast_hits'])).sum(),
        df_valid['hhsearch_candidate'].sum() if 'hhsearch_candidate' in df_valid.columns else 0
    ]
    
    plt.bar(categories, values, alpha=0.7)
    plt.xlabel('Category')
    plt.ylabel('Number of Proteins')
    plt.title('Domain Generation Quality Metrics')
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(output_dir, f'batch_{batch_id}_quality_metrics.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 5. Error analysis visualization (if error data is available)
    if 'error_class' in df.columns:
        error_counts = df['error_class'].value_counts()
        if not error_counts.empty and len(error_counts) > 0:
            plt.figure(figsize=(12, 6))
            error_counts.plot(kind='bar', alpha=0.7)
            plt.xlabel('Error Type')
            plt.ylabel('Number of Proteins')
            plt.title('Error Classification Distribution')
            plt.grid(True, alpha=0.3)
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'batch_{batch_id}_error_classes.png'), dpi=300, bbox_inches='tight')
            plt.close()


def main():
    """Main entry point for domain analysis script"""
    parser = argparse.ArgumentParser(description='Analyze domain generation results for a batch')
    parser.add_argument('--config', type=str, default='config/config.yml', 
                      help='Path to ECOD configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to analyze')
    parser.add_argument('--output-dir', type=str, default='analysis_results',
                      help='Directory for output files')
    parser.add_argument('--prompt-password', action='store_true',
                      help='Prompt for database password instead of using config')
    parser.add_argument('--error-analysis', action='store_true',
                      help='Perform detailed error analysis')
    parser.add_argument('--sample-logs', type=int, default=5,
                      help='Number of sample logs to retrieve per error type')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Check dependencies
    check_dependencies()
    
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
        
        # Analyze error patterns if requested
        if args.error_analysis:
            logger.info(f"Performing detailed error analysis for batch {args.batch_id}")
            error_analysis = analyze_error_patterns(processes)
            
            # Save error analysis
            error_analysis_path = os.path.join(args.output_dir, f'batch_{args.batch_id}_error_analysis.json')
            with open(error_analysis_path, 'w') as f:
                json.dump(error_analysis, f, indent=2)
            
            # Create DataFrame for error details
            error_rows = []
            
            for error_class, details in error_analysis['error_details'].items():
                for detail in details:
                    error_rows.append({
                        'error_class': error_class,
                        'pdb_id': detail['pdb_id'],
                        'chain_id': detail['chain_id'],
                        'process_id': detail['process_id'],
                        'chain_length': detail['chain_length'],
                        'status_message': detail['status_message']
                    })
            
            error_df = pd.DataFrame(error_rows)
            error_df_path = os.path.join(args.output_dir, f'batch_{args.batch_id}_error_details.csv')
            error_df.to_csv(error_df_path, index=False)
            
            # Sample logs for each error type
            if args.sample_logs > 0:
                log_samples = {}
                
                for error_class, details in error_analysis['error_details'].items():
                    # Select a sample of processes for this error type
                    sample_size = min(args.sample_logs, len(details))
                    sample_processes = details[:sample_size]
                    
                    log_samples[error_class] = {}
                    
                    for process in sample_processes:
                        process_id = process['process_id']
                        protein_key = f"{process['pdb_id']}_{process['chain_id']}"
                        
                        log_summary = get_process_log_summary(conn, process_id)
                        log_samples[error_class][protein_key] = {
                            'process_id': process_id,
                            'log_summary': log_summary
                        }
                
                # Save log samples
                log_samples_path = os.path.join(args.output_dir, f'batch_{args.batch_id}_error_logs.json')
                with open(log_samples_path, 'w') as f:
                    json.dump(log_samples, f, indent=2)
            
            # Print error analysis summary
            print("\n===== ERROR ANALYSIS SUMMARY =====")
            print(f"Total errors: {error_analysis['total_errors']}")
            
            print("\nError classifications:")
            for error_class, count in sorted(error_analysis['error_classifications'].items(), key=lambda x: x[1], reverse=True):
                print(f"  - {error_class}: {count} ({100 * count / error_analysis['total_errors']:.1f}%)")
            
            # Print errors by processing stage
            if 'stage_errors' in error_analysis and error_analysis['stage_errors']:
                print("\nErrors by processing stage:")
                for stage, count in sorted(error_analysis['stage_errors'].items(), key=lambda x: x[1], reverse=True):
                    print(f"  - {stage}: {count} ({100 * count / error_analysis['total_errors']:.1f}%)")
            
            # Print non-representative errors
            if 'non_representative_errors' in error_analysis and error_analysis['non_representative_errors']:
                non_rep_count = len(error_analysis['non_representative_errors'])
                print(f"\nNon-representative chains with errors: {non_rep_count} ({100 * non_rep_count / error_analysis['total_errors']:.1f}%)")
                if non_rep_count > 0 and non_rep_count <= 5:  # Show details if there are just a few
                    for chain in error_analysis['non_representative_errors']:
                        print(f"  - {chain['pdb_id']}_{chain['chain_id']} (process {chain['process_id']})")
            
            # Print errors by source pattern
            if 'source_pattern_errors' in error_analysis and error_analysis['source_pattern_errors']:
                print("\nErrors by source pattern:")
                for pattern, count in sorted(error_analysis['source_pattern_errors'].items(), key=lambda x: x[1], reverse=True):
                    print(f"  - {pattern}: {count} ({100 * count / error_analysis['total_errors']:.1f}%)")
            
            if error_analysis['high_error_pdbs']:
                print("\nPDBs with high error rates:")
                for pdb_id, info in error_analysis['high_error_pdbs'].items():
                    print(f"  - {pdb_id}: {info['error_count']}/{info['total_chains']} chains failed ({100 * info['error_rate']:.1f}%)")
            
            print(f"\nDetailed error information saved to:")
            print(f"  - {error_analysis_path}")
            print(f"  - {error_df_path}")
            if args.sample_logs > 0:
                print(f"  - {log_samples_path}")
        
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
                error_class = classify_error(proc) if proc['status'] == 'error' else proc['status']
                logger.warning(f"Skipping {protein_key} (process {process_id}): Status is {proc['status']} ({error_class})")
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
            domain_file_status = 'missing'
            
            if domain_path:
                domain_file_status = detect_file_type(domain_path)
                if domain_file_status == 'json':
                    domain_data = load_domain_file(domain_path)
                    domain_analysis = analyze_domains(domain_data)
                elif domain_file_status != 'missing':
                    logger.info(f"Domain file exists but is {domain_file_status}: {domain_path}")
            
            # Process BLAST files
            chain_blast_data = {}
            chain_blast_analysis = {'hit_count': 0, 'has_hits': False, 'error': True}
            chain_blast_status = 'missing'
            
            if chain_blast_path:
                chain_blast_status = detect_file_type(chain_blast_path)
                if chain_blast_status in ['json', 'xml']:
                    chain_blast_data = load_blast_file(chain_blast_path)
                    chain_blast_analysis = analyze_blast_hits(chain_blast_data)
                elif chain_blast_status != 'missing':
                    logger.info(f"Chain BLAST file exists but is {chain_blast_status}: {chain_blast_path}")
            
            domain_blast_data = {}
            domain_blast_analysis = {'hit_count': 0, 'has_hits': False, 'error': True}
            domain_blast_status = 'missing'
            
            if domain_blast_path:
                domain_blast_status = detect_file_type(domain_blast_path)
                if domain_blast_status in ['json', 'xml']:
                    domain_blast_data = load_blast_file(domain_blast_path)
                    domain_blast_analysis = analyze_blast_hits(domain_blast_data)
                elif domain_blast_status != 'missing':
                    logger.info(f"Domain BLAST file exists but is {domain_blast_status}: {domain_blast_path}")
            
            # Record file statuses
            file_statuses = {
                'domain_file': domain_file_status,
                'chain_blast_file': chain_blast_status,
                'domain_blast_file': domain_blast_status
            }
            
            # Check for empty BLAST files
            blast_file_errors = (
                chain_blast_analysis.get('error', True) or 
                domain_blast_analysis.get('error', True)
            )
            
            # Determine HHsearch candidacy with improved criteria
            # Consider blast file errors in the assessment
            has_blast_hits = (
                chain_blast_analysis.get('has_hits', False) or 
                domain_blast_analysis.get('has_hits', False)
            )
            few_domains = 0 <= domain_analysis.get('domain_count', 0) <= 1
            
            # Only consider as HHsearch candidate if BLAST files were properly parsed
            hhsearch_candidate = has_blast_hits and few_domains and not blast_file_errors
            
            # Store results with additional file status information
            results.append({
                'process_id': process_id,
                'protein_id': proc['protein_id'],
                'pdb_id': pdb_id,
                'chain_id': chain_id,
                'protein_key': protein_key,
                
                # Domain analysis
                'domain_count': domain_analysis.get('domain_count', 0),
                'domain_source': domain_analysis.get('source', None),
                'domain_source_counts': domain_analysis.get('source_counts', {}),
                
                # BLAST analysis
                'chain_blast_hits': chain_blast_analysis.get('hit_count', 0),
                'domain_blast_hits': domain_blast_analysis.get('hit_count', 0),
                'has_chain_blast_hits': chain_blast_analysis.get('has_hits', False),
                'has_domain_blast_hits': domain_blast_analysis.get('has_hits', False),
                
                # File status information
                'file_statuses': file_statuses,
                'blast_file_errors': blast_file_errors,
                
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
        
        # File status analysis
        file_status_counts = {
            'domain_file': Counter(),
            'chain_blast_file': Counter(),
            'domain_blast_file': Counter()
        }
        
        for _, row in df.iterrows():
            file_statuses = row.get('file_statuses', {})
            for file_type, status in file_statuses.items():
                file_status_counts[file_type][status] += 1
        
        proteins_with_blast_errors = df[df['blast_file_errors']].shape[0]
        
        # Domain source breakdown
        domain_source_counts = Counter()
        for _, row in df.iterrows():
            if row['domain_count'] > 0 and row['domain_source']:
                domain_source_counts[row['domain_source']] += 1
        
        # BLAST hit statistics
        proteins_with_chain_hits = df[df['has_chain_blast_hits']].shape[0]
        proteins_with_domain_hits = df[df['has_domain_blast_hits']].shape[0]
        proteins_with_any_hits = df[df['has_chain_blast_hits'] | df['has_domain_blast_hits']].shape[0]
        proteins_with_no_hits = df[~df['has_chain_blast_hits'] & ~df['has_domain_blast_hits']].shape[0]
        
        # Filter for proteins without BLAST file errors
        df_valid_blast = df[~df['blast_file_errors']]
        valid_total = len(df_valid_blast)
        
        # Valid proteins statistics (only those with properly parsed BLAST files)
        valid_with_domains = df_valid_blast[df_valid_blast['domain_count'] > 0].shape[0]
        valid_no_domains = df_valid_blast[df_valid_blast['domain_count'] == 0].shape[0]
        valid_with_hits = df_valid_blast[
            df_valid_blast['has_chain_blast_hits'] | df_valid_blast['has_domain_blast_hits']
        ].shape[0]
        valid_no_hits = df_valid_blast[
            ~df_valid_blast['has_chain_blast_hits'] & ~df_valid_blast['has_domain_blast_hits']
        ].shape[0]
        
        # No domains but has hits (valid files only)
        valid_no_domains_with_hits = df_valid_blast[
            (df_valid_blast['domain_count'] == 0) & 
            (df_valid_blast['has_chain_blast_hits'] | df_valid_blast['has_domain_blast_hits'])
        ].shape[0]
        
        # HHsearch candidates (already filtered for valid BLAST files)
        hhsearch_candidates = df[df['hhsearch_candidate']].shape[0]
        
        # Create summary report
        summary = {
            'batch_id': args.batch_id,
            'total_proteins': total_proteins,
            'proteins_with_domains': proteins_with_domains,
            'proteins_no_domains': proteins_no_domains,
            'domain_source_breakdown': dict(domain_source_counts),
            
            # File status information
            'file_status_counts': {k: dict(v) for k, v in file_status_counts.items()},
            'proteins_with_blast_errors': proteins_with_blast_errors,
            
            # BLAST hit statistics (all proteins)
            'proteins_with_chain_hits': proteins_with_chain_hits,
            'proteins_with_domain_hits': proteins_with_domain_hits,
            'proteins_with_any_hits': proteins_with_any_hits,
            'proteins_with_no_hits': proteins_with_no_hits,
            
            # Valid files statistics (excluding error files)
            'valid_total': valid_total,
            'valid_with_domains': valid_with_domains,
            'valid_no_domains': valid_no_domains,
            'valid_with_hits': valid_with_hits,
            'valid_no_hits': valid_no_hits,
            'valid_no_domains_with_hits': valid_no_domains_with_hits,
            
            # HHsearch candidates
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
        
        # File status summary
        print("\nFile Status Summary:")
        for file_type, counts in file_status_counts.items():
            print(f"  {file_type}:")
            for status, count in counts.most_common():
                print(f"    - {status}: {count} ({100*count/total_proteins:.1f}%)")
        
        print(f"\nProteins with BLAST file errors: {proteins_with_blast_errors} ({100*proteins_with_blast_errors/total_proteins:.1f}%)")
        print(f"Proteins with valid BLAST files: {valid_total} ({100*valid_total/total_proteins:.1f}%)")
        
        print("\n--- Overall Statistics (including error files) ---")
        print(f"Proteins with domains: {proteins_with_domains} ({100*proteins_with_domains/total_proteins:.1f}%)")
        print(f"Proteins without domains: {proteins_no_domains} ({100*proteins_no_domains/total_proteins:.1f}%)")
        
        print("\nDomain source breakdown:")
        for source, count in domain_source_counts.most_common():
            print(f"  - {source}: {count} ({100*count/proteins_with_domains:.1f}%)")
        
        print("\nBLAST hit statistics (all files):")
        print(f"Proteins with chain BLAST hits: {proteins_with_chain_hits} ({100*proteins_with_chain_hits/total_proteins:.1f}%)")
        print(f"Proteins with domain BLAST hits: {proteins_with_domain_hits} ({100*proteins_with_domain_hits/total_proteins:.1f}%)")
        print(f"Proteins with any BLAST hits: {proteins_with_any_hits} ({100*proteins_with_any_hits/total_proteins:.1f}%)")
        print(f"Proteins with no BLAST hits: {proteins_with_no_hits} ({100*proteins_with_no_hits/total_proteins:.1f}%)")
        
        if valid_total > 0:
            # Only show valid statistics if we have any valid files
            print("\n--- Valid Files Only (excluding error files) ---")
            print(f"Valid proteins with domains: {valid_with_domains} ({100*valid_with_domains/valid_total:.1f}%)")
            print(f"Valid proteins without domains: {valid_no_domains} ({100*valid_no_domains/valid_total:.1f}%)")
            
            print("\nValid BLAST hit statistics:")
            print(f"Valid proteins with BLAST hits: {valid_with_hits} ({100*valid_with_hits/valid_total:.1f}%)")
            print(f"Valid proteins with no BLAST hits: {valid_no_hits} ({100*valid_no_hits/valid_total:.1f}%)")
            
            print("\nQuality metrics (valid files only):")
            print(f"Valid proteins with no domains despite BLAST hits: {valid_no_domains_with_hits} ({100*valid_no_domains_with_hits/valid_total:.1f}%)")
        
        print("\nHHsearch candidacy:")
        print(f"Potential HHsearch candidates: {hhsearch_candidates} ({100*hhsearch_candidates/valid_total:.1f}% of valid proteins)")
        
        # Create a file with error cases for further investigation
        error_df = df[df['blast_file_errors']]
        if not error_df.empty:
            error_path = os.path.join(args.output_dir, f'batch_{args.batch_id}_file_errors.csv')
            error_df.to_csv(error_path, index=False)
            print(f"\nList of {len(error_df)} proteins with file errors saved to: {error_path}")
        
        # Generate visualizations
        create_visualizations(df, args.output_dir, args.batch_id)
        
        print(f"\nResults saved to {args.output_dir}/")
        
    finally:
        # Close database connection
        conn.close()


if __name__ == "__main__":
    main()