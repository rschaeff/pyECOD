#!/usr/bin/env python3
"""
Investigate Specific Protein Failures

This script takes specific protein IDs and examines their processing history,
files, and potential causes of failure in detail.

Usage:
    python scripts/investigate_specific_failures.py --config config.yml --proteins "8r0f_B,8r0g_A,8r15_A" [options]
"""

import os
import sys
import logging
import argparse
import json
import psycopg2
from psycopg2.extras import RealDictCursor
import yaml
import xml.etree.ElementTree as ET
from datetime import datetime
from typing import List, Dict, Any, Optional


def setup_logging(verbose=False):
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    format_str = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=level, format=format_str)
    return logging.getLogger(__name__)


def parse_config(config_path):
    """Parse configuration file."""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def get_db_connection(config):
    """Create database connection from config."""
    db_config = config.get('database', {})
    
    try:
        conn = psycopg2.connect(
            host=db_config.get('host', 'dione'),
            port=db_config.get('port', 45000),
            dbname=db_config.get('name', 'ecod_protein'),
            user=db_config.get('user', 'ecod'),
            password=db_config.get('password', '')
        )
        return conn
    except psycopg2.Error as e:
        logging.error(f"Database connection error: {e}")
        raise


def get_protein_details(conn, source_ids):
    """Get detailed information about specific proteins."""
    
    placeholders = ','.join(['%s'] * len(source_ids))
    query = f"""
    SELECT 
        ep.id,
        ep.source_id,
        ep.pdb_id,
        ep.chain_id,
        ep.name,
        ep.type,
        ep.tax_id,
        ep.length,
        ep.created_at,
        
        -- Process status
        ps.batch_id,
        ps.current_stage,
        ps.status,
        ps.error_message,
        ps.is_representative,
        ps.updated_at as process_updated,
        
        -- Batch info
        b.batch_name,
        b.ref_version,
        b.type as batch_type
        
    FROM ecod_schema.protein ep
    LEFT JOIN ecod_schema.process_status ps ON ep.id = ps.protein_id
    LEFT JOIN ecod_schema.batch b ON ps.batch_id = b.id
    WHERE ep.source_id IN ({placeholders})
    ORDER BY ep.source_id, ps.batch_id DESC
    """
    
    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(query, source_ids)
        return cur.fetchall()


def get_protein_files(conn, protein_ids):
    """Get file information for specific proteins."""
    
    query = """
    SELECT 
        ps.protein_id,
        ep.source_id,
        pf.file_type,
        pf.file_path,
        pf.file_exists,
        pf.file_size,
        pf.created_at,
        pf.updated_at
    FROM ecod_schema.process_file pf
    JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
    JOIN ecod_schema.protein ep ON ps.protein_id = ep.id
    WHERE ps.protein_id = ANY(%s)
    ORDER BY ep.source_id, pf.file_type, pf.created_at
    """
    
    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(query, (protein_ids,))
        return cur.fetchall()


def get_partition_attempts(conn, source_ids):
    """Check if there are any partition attempts for these proteins."""
    
    query = """
    SELECT 
        pp.pdb_id,
        pp.chain_id,
        pp.batch_id,
        pp.is_classified,
        pp.is_unclassified,
        pp.sequence_length,
        pp.coverage,
        pp.created_at,
        COUNT(pd.id) as domain_count
    FROM pdb_analysis.partition_proteins pp
    LEFT JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
    WHERE (pp.pdb_id || '_' || pp.chain_id) = ANY(%s)
    GROUP BY pp.id, pp.pdb_id, pp.chain_id, pp.batch_id, pp.is_classified, 
             pp.is_unclassified, pp.sequence_length, pp.coverage, pp.created_at
    ORDER BY pp.pdb_id, pp.chain_id, pp.batch_id
    """
    
    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(query, (source_ids,))
        return cur.fetchall()


def examine_file_content(file_path, file_type):
    """Examine the content of a file to understand why processing might have failed."""
    
    if not os.path.exists(file_path):
        return {"error": "File does not exist", "path": file_path}
    
    try:
        file_size = os.path.getsize(file_path)
        
        if file_size == 0:
            return {"error": "Empty file", "size": 0, "path": file_path}
        
        # For XML files, try to parse
        if file_type in ['domain_summary'] or file_path.endswith('.xml'):
            try:
                tree = ET.parse(file_path)
                root = tree.getroot()
                
                analysis = {
                    "valid_xml": True,
                    "root_tag": root.tag,
                    "size": file_size,
                    "path": file_path
                }
                
                # For domain summary files, check structure
                if file_type == 'domain_summary':
                    analysis.update({
                        "chain_blast_hits": len(root.findall(".//chain_blast_run//hit")),
                        "domain_blast_hits": len(root.findall(".//blast_run//hit")),
                        "hhsearch_hits": len(root.findall(".//hh_run//hit")),
                        "domain_suggestions": len(root.findall(".//domain_suggestions//domain"))
                    })
                
                return analysis
                
            except ET.ParseError as e:
                return {
                    "error": f"XML parse error: {str(e)}",
                    "size": file_size,
                    "path": file_path
                }
        
        # For other file types, just return basic info
        return {
            "valid_file": True,
            "size": file_size,
            "path": file_path
        }
        
    except Exception as e:
        return {
            "error": f"File access error: {str(e)}",
            "path": file_path
        }


def analyze_sequence_characteristics(conn, source_ids):
    """Get sequence characteristics if available."""
    
    # Try to get sequence from PDB analysis tables
    query = """
    SELECT 
        p.source_id,
        ps.sequence,
        LENGTH(ps.sequence) as sequence_length,
        -- Basic composition analysis
        (LENGTH(ps.sequence) - LENGTH(REPLACE(ps.sequence, 'G', ''))) as glycine_count,
        (LENGTH(ps.sequence) - LENGTH(REPLACE(ps.sequence, 'P', ''))) as proline_count,
        (LENGTH(ps.sequence) - LENGTH(REPLACE(ps.sequence, 'X', ''))) as unknown_count
    FROM pdb_analysis.protein p
    JOIN pdb_analysis.protein_sequence ps ON p.id = ps.protein_id
    WHERE p.source_id = ANY(%s)
    """
    
    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(query, (source_ids,))
        return cur.fetchall()


def investigate_protein(conn, source_id, base_data_dir="/data/ecod/pdb_updates"):
    """Thoroughly investigate a single protein failure."""
    
    print(f"\n{'='*80}")
    print(f"INVESTIGATING: {source_id}")
    print(f"{'='*80}")
    
    # Get protein details
    protein_details = get_protein_details(conn, [source_id])
    if not protein_details:
        print(f"❌ No protein data found for {source_id}")
        return
    
    protein = protein_details[0]  # Latest entry
    print(f"📋 Basic Info:")
    print(f"   PDB ID: {protein['pdb_id']}")
    print(f"   Chain: {protein['chain_id']}")  
    print(f"   Length: {protein['length']}")
    print(f"   Batch: {protein['batch_id']} ({protein['batch_name']})")
    print(f"   Status: {protein['status']}")
    print(f"   Stage: {protein['current_stage']}")
    print(f"   Representative: {protein['is_representative']}")
    
    if protein['error_message']:
        print(f"   ❌ Error: {protein['error_message']}")
    
    # Get files
    protein_files = get_protein_files(conn, [protein['id']])
    print(f"\n📁 File Status:")
    
    file_analysis = {}
    for file_info in protein_files:
        file_type = file_info['file_type']
        file_path = file_info['file_path']
        file_exists = file_info['file_exists']
        
        status = "✅" if file_exists else "❌"
        print(f"   {status} {file_type}: {file_path}")
        
        if file_exists and file_path:
            content_analysis = examine_file_content(file_path, file_type)
            file_analysis[file_type] = content_analysis
            
            if 'error' in content_analysis:
                print(f"      ⚠️  {content_analysis['error']}")
            elif file_type == 'domain_summary':
                hits = content_analysis.get('chain_blast_hits', 0) + content_analysis.get('domain_blast_hits', 0) + content_analysis.get('hhsearch_hits', 0)
                print(f"      📊 Evidence: {hits} total hits")
                if content_analysis.get('domain_suggestions', 0) > 0:
                    print(f"      🎯 Domain suggestions: {content_analysis['domain_suggestions']}")
    
    # Check for partition attempts
    partition_attempts = get_partition_attempts(conn, [source_id])
    if partition_attempts:
        print(f"\n🔄 Partition Attempts:")
        for attempt in partition_attempts:
            status = "✅ Classified" if attempt['is_classified'] else "❌ Failed"
            print(f"   Batch {attempt['batch_id']}: {status} (length: {attempt['sequence_length']}, domains: {attempt['domain_count']})")
    else:
        print(f"\n❌ No partition attempts found")
    
    # Get sequence characteristics
    sequences = analyze_sequence_characteristics(conn, [source_id])
    if sequences:
        seq = sequences[0]
        print(f"\n🧬 Sequence Analysis:")
        print(f"   Length: {seq['sequence_length']}")
        print(f"   Glycine %: {(seq['glycine_count']/seq['sequence_length']*100):.1f}%")
        print(f"   Proline %: {(seq['proline_count']/seq['sequence_length']*100):.1f}%")
        if seq['unknown_count'] > 0:
            print(f"   Unknown residues: {seq['unknown_count']}")
        
        # Check if it's likely a peptide
        if seq['sequence_length'] < 30:
            print(f"   ⚠️  PEPTIDE LENGTH - likely should be filtered")
    
    # Try to identify the specific failure mode
    print(f"\n🔍 Failure Analysis:")
    
    # Peptide detection issue
    if protein['length'] and protein['length'] < 30:
        print(f"   💡 HYPOTHESIS: Peptide not properly filtered (length: {protein['length']})")
    
    # File processing issues
    if 'domain_summary' in file_analysis:
        summary_analysis = file_analysis['domain_summary']
        if 'error' in summary_analysis:
            print(f"   💡 HYPOTHESIS: Domain summary file corrupted or malformed")
        elif summary_analysis.get('chain_blast_hits', 0) == 0 and summary_analysis.get('domain_blast_hits', 0) == 0:
            print(f"   💡 HYPOTHESIS: No BLAST evidence found - may be truly novel")
        elif summary_analysis.get('hhsearch_hits', 0) == 0:
            print(f"   💡 HYPOTHESIS: BLAST evidence exists but no HHSearch hits")
    
    # Processing pipeline issues
    has_summary = any(f['file_type'] == 'domain_summary' and f['file_exists'] for f in protein_files)
    has_partition = any(f['file_type'].endswith('partition') and f['file_exists'] for f in protein_files)
    
    if has_summary and not has_partition:
        print(f"   💡 HYPOTHESIS: Domain summary exists but partition failed - algorithm issue")
    elif not has_summary:
        print(f"   💡 HYPOTHESIS: Failed before domain summary stage - upstream pipeline issue")
    
    return {
        'source_id': source_id,
        'protein_details': dict(protein),
        'file_analysis': file_analysis,
        'partition_attempts': [dict(p) for p in partition_attempts],
        'sequence_analysis': dict(sequences[0]) if sequences else None
    }


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Investigate specific protein failures in detail')
    parser.add_argument('--config', required=True, help='Path to configuration file')
    parser.add_argument('--proteins', required=True, help='Comma-separated protein source IDs (e.g., "8r0f_B,8r0g_A")')
    parser.add_argument('--output', help='Output JSON file path')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')

    args = parser.parse_args()

    # Set up logging
    logger = setup_logging(args.verbose)

    # Parse config
    try:
        config = parse_config(args.config)
    except Exception as e:
        logger.error(f"Error parsing config file: {str(e)}")
        sys.exit(1)

    # Parse protein IDs
    try:
        source_ids = [pid.strip() for pid in args.proteins.split(',')]
        logger.info(f"Investigating proteins: {source_ids}")
    except ValueError:
        logger.error("Invalid protein IDs format. Use comma-separated source IDs like '8r0f_B,8r0g_A'")
        sys.exit(1)

    # Get database connection
    try:
        conn = get_db_connection(config)
    except Exception as e:
        logger.error(f"Error connecting to database: {str(e)}")
        sys.exit(1)

    try:
        # Investigate each protein
        investigation_results = []
        
        for source_id in source_ids:
            result = investigate_protein(conn, source_id)
            if result:
                investigation_results.append(result)
        
        # Prepare output
        output_data = {
            'timestamp': datetime.now().isoformat(),
            'investigated_proteins': source_ids,
            'total_investigated': len(investigation_results),
            'results': investigation_results
        }
        
        # Output results
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(output_data, f, indent=2, default=str)
            logger.info(f"Investigation results written to {args.output}")
        
        print(f"\n{'='*80}")
        print(f"INVESTIGATION COMPLETE")
        print(f"Investigated {len(investigation_results)} proteins")
        print(f"{'='*80}")

    except Exception as e:
        logger.error(f"Error during investigation: {str(e)}")
        if args.verbose:
            import traceback
            logger.error(traceback.format_exc())
        sys.exit(1)

    finally:
        conn.close()


if __name__ == "__main__":
    main()
