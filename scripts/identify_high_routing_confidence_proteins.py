#!/usr/bin/env python3
"""
This script identifies proteins with high confidence BLAST matches that can be
routed through the BLAST-only pipeline with the new 0.7 threshold.
"""

import argparse
import logging
import sys
import os
from collections import defaultdict
import math

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

# Import from ecod package
from ecod.core.context import ApplicationContext
from ecod.utils.file import safe_open

def evaluate_blast_confidence(domain_blast_results, chain_blast_results, sequence_length):
    """
    Evaluate confidence in BLAST results with the improved algorithm
    
    Args:
        domain_blast_results: List of domain BLAST hit dictionaries
        chain_blast_results: List of chain BLAST hit dictionaries
        sequence_length: Length of the query protein sequence
    
    Returns:
        float: Confidence score between 0 and 1
    """
    if not domain_blast_results and not chain_blast_results:
        return 0.0
    
    # --- Domain BLAST metrics ---
    # Consider top 5 domain hits
    domain_top_hits = domain_blast_results[:5] if len(domain_blast_results) >= 5 else domain_blast_results
    
    # Calculate domain metrics
    domain_min_evalue = 10.0
    domain_max_coverage = 0.0
    domain_max_identity = 0.0
    domain_hit_consistency = 0.0
    
    if domain_top_hits:
        # Extract minimum E-value
        all_evalues = []
        for hit in domain_top_hits:
            evalue_str = hit.get('evalue', '10.0')
            try:
                evalue = float(evalue_str)
                all_evalues.append(evalue)
            except (ValueError, TypeError):
                pass
        
        if all_evalues:
            domain_min_evalue = min(all_evalues)
        
        # Calculate coverage
        all_coverages = []
        for hit in domain_top_hits:
            query_regions = hit.get('query_regions', '').split(',')
            if not query_regions or query_regions[0] == '':
                continue
                
            positions = set()
            for region in query_regions:
                if '-' in region:
                    try:
                        start, end = map(int, region.split('-'))
                        positions.update(range(start, end + 1))
                    except (ValueError, TypeError):
                        continue
            
            if sequence_length > 0:
                coverage = len(positions) / sequence_length
                all_coverages.append(coverage)
        
        if all_coverages:
            domain_max_coverage = max(all_coverages)
        
        # Calculate identity
        all_identities = []
        for hit in domain_top_hits:
            identity_str = hit.get('identity', '0.0')
            try:
                identity = float(identity_str)
                all_identities.append(identity)
            except (ValueError, TypeError):
                pass
        
        if all_identities:
            domain_max_identity = max(all_identities)
        
        # Evaluate hit consistency based on T-group and H-group
        domain_groups = defaultdict(int)
        for hit in domain_top_hits:
            t_group = hit.get('t_group', '')
            h_group = hit.get('h_group', '')
            key = f"{t_group}_{h_group}"
            domain_groups[key] += 1
        
        if domain_groups:
            most_common_count = max(domain_groups.values())
            domain_hit_consistency = most_common_count / len(domain_top_hits)
    
    # --- Chain BLAST metrics ---
    # Consider top 5 chain hits
    chain_top_hits = chain_blast_results[:5] if len(chain_blast_results) >= 5 else chain_blast_results
    
    # Calculate chain metrics
    chain_min_evalue = 10.0
    chain_max_coverage = 0.0
    chain_max_identity = 0.0
    chain_hit_consistency = 0.0
    
    if chain_top_hits:
        # Extract minimum E-value
        all_evalues = []
        for hit in chain_top_hits:
            evalue_str = hit.get('evalue', '10.0')
            try:
                evalue = float(evalue_str)
                all_evalues.append(evalue)
            except (ValueError, TypeError):
                pass
        
        if all_evalues:
            chain_min_evalue = min(all_evalues)
        
        # Calculate coverage similar to domain hits
        all_coverages = []
        for hit in chain_top_hits:
            query_regions = hit.get('query_regions', '').split(',')
            if not query_regions or query_regions[0] == '':
                continue
                
            positions = set()
            for region in query_regions:
                if '-' in region:
                    try:
                        start, end = map(int, region.split('-'))
                        positions.update(range(start, end + 1))
                    except (ValueError, TypeError):
                        continue
            
            if sequence_length > 0:
                coverage = len(positions) / sequence_length
                all_coverages.append(coverage)
        
        if all_coverages:
            chain_max_coverage = max(all_coverages)
        
        # Calculate identity
        all_identities = []
        for hit in chain_top_hits:
            identity_str = hit.get('identity', '0.0')
            try:
                identity = float(identity_str)
                all_identities.append(identity)
            except (ValueError, TypeError):
                pass
        
        if all_identities:
            chain_max_identity = max(all_identities)
        
        # Simpler hit consistency for chains (based on PDB IDs)
        chain_ids = defaultdict(int)
        for hit in chain_top_hits:
            pdb_id = hit.get('pdb_id', '')
            chain_ids[pdb_id] += 1
        
        if chain_ids:
            most_common_count = max(chain_ids.values())
            chain_hit_consistency = most_common_count / len(chain_top_hits)
    
    # --- Combine metrics ---
    # Take best values from domain and chain BLAST
    min_evalue = min(domain_min_evalue, chain_min_evalue)
    max_coverage = max(domain_max_coverage, chain_max_coverage)
    max_identity = max(domain_max_identity, chain_max_identity)
    hit_consistency = max(domain_hit_consistency, chain_hit_consistency)
    
    # Calculate architecture coverage (how completely the hits cover the protein)
    all_positions = set()
    
    # Collect positions from domain hits
    for hit in domain_blast_results:
        query_regions = hit.get('query_regions', '').split(',')
        for region in query_regions:
            if '-' in region:
                try:
                    start, end = map(int, region.split('-'))
                    all_positions.update(range(start, end + 1))
                except (ValueError, TypeError):
                    continue
    
    # Calculate architecture coverage
    architecture_coverage = len(all_positions) / sequence_length if sequence_length > 0 else 0.0
    
    # Convert E-value to confidence score
    evalue_confidence = 0.0
    if min_evalue < 1e-300:
        evalue_confidence = 0.99
    else:
        # This gives higher confidence for lower E-values
        evalue_confidence = min(1.0, max(0.0, 1.0 - (1.0 / (1.0 + math.exp(-math.log10(min_evalue) - 3)))))
    
    # Calculate weighted confidence score
    confidence = (
        0.3 * evalue_confidence +        # E-value importance
        0.3 * max_coverage +             # Coverage is very important
        0.2 * max_identity +             # Identity is important
        0.1 * hit_consistency +          # Consistency is helpful
        0.1 * architecture_coverage      # Architecture completeness is helpful
    )
    
    return min(1.0, max(0.0, confidence))  # Ensure value between 0-1

def parse_blast_xml(blast_file):
    """
    Parse a BLAST XML file to extract hit information
    
    Args:
        blast_file: Path to BLAST XML file
        
    Returns:
        list: List of hit dictionaries
    """
    import xml.etree.ElementTree as ET
    
    hits = []
    try:
        tree = ET.parse(blast_file)
        root = tree.getroot()
        
        # Extract query length
        query_len = 0
        for iteration in root.findall(".//Iteration"):
            query_len_elem = iteration.find("Iteration_query-len")
            if query_len_elem is not None:
                query_len = int(query_len_elem.text)
                break
        
        # Process hits
        for iteration in root.findall(".//Iteration"):
            for hit in iteration.findall(".//Hit"):
                hit_num = hit.findtext("Hit_num", "")
                hit_id = hit.findtext("Hit_id", "")
                hit_def = hit.findtext("Hit_def", "")
                hit_len = int(hit.findtext("Hit_len", "0"))
                
                # Parse hit definition
                pdb_id = "unknown"
                chain_id = "unknown"
                domain_id = None
                
                # For domain BLAST hits
                import re
                domain_match = re.search(r"((d|g|e)(\d\w{3})\w+\d+)", hit_def)
                if domain_match:
                    domain_id = domain_match.group(1)
                    pdb_id = domain_match.group(3)
                else:
                    # For chain BLAST hits
                    pdb_match = re.search(r"(\d\w{3})\s+([A-Za-z0-9])", hit_def)
                    if pdb_match:
                        pdb_id = pdb_match.group(1)
                        chain_id = pdb_match.group(2)
                
                # Process HSPs
                evalues = []
                query_regions = []
                hit_regions = []
                identities = []
                
                for hsp in hit.findall(".//Hsp"):
                    hsp_evalue = float(hsp.findtext("Hsp_evalue", "999"))
                    
                    # Get alignment coordinates
                    hsp_query_from = int(hsp.findtext("Hsp_query-from", "0"))
                    hsp_query_to = int(hsp.findtext("Hsp_query-to", "0"))
                    hsp_hit_from = int(hsp.findtext("Hsp_hit-from", "0"))
                    hsp_hit_to = int(hsp.findtext("Hsp_hit-to", "0"))
                    
                    # Calculate identity percentage
                    hsp_identity = int(hsp.findtext("Hsp_identity", "0"))
                    hsp_align_len = int(hsp.findtext("Hsp_align-len", "0"))
                    identity_pct = hsp_identity / hsp_align_len if hsp_align_len > 0 else 0
                    
                    evalues.append(hsp_evalue)
                    query_regions.append(f"{hsp_query_from}-{hsp_query_to}")
                    hit_regions.append(f"{hsp_hit_from}-{hsp_hit_to}")
                    identities.append(identity_pct)
                
                # Skip hits with no HSPs
                if not evalues:
                    continue
                
                # Create hit dictionary
                hit_dict = {
                    'hit_num': hit_num,
                    'hit_id': hit_id,
                    'pdb_id': pdb_id,
                    'chain_id': chain_id,
                    'evalue': min(evalues),
                    'query_regions': ','.join(query_regions),
                    'hit_regions': ','.join(hit_regions),
                    'identity': max(identities) if identities else 0.0
                }
                
                if domain_id:
                    hit_dict['domain_id'] = domain_id
                
                hits.append(hit_dict)
                
        return hits, query_len
    except Exception as e:
        logging.error(f"Error parsing BLAST file {blast_file}: {e}")
        return [], 0

def main():
    parser = argparse.ArgumentParser(description='Identify high-confidence proteins for BLAST-only pipeline')
    parser.add_argument('--config', required=True, help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True, help='Batch ID to analyze')
    parser.add_argument('--threshold', type=float, default=0.7, help='Confidence threshold (default: 0.7)')
    parser.add_argument('--output', help='Output file for high-confidence proteins')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose output')
    args = parser.parse_args()
    
    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger('identify_high_confidence')
    
    # Initialize application context
    try:
        app_context = ApplicationContext(args.config)
        config = app_context.config_manager.config
        db = app_context.db
    except Exception as e:
        logger.error(f"Failed to initialize application: {e}")
        sys.exit(1)
    
    # Get proteins with BLAST results from the batch
    logger.info(f"Getting proteins with BLAST results for batch {args.batch_id}")
    
    query = """
        SELECT 
            p.id, p.source_id, ps.sequence, ps.sequence_length,
            b.base_path
        FROM 
            ecod_schema.protein p
        JOIN
            ecod_schema.protein_sequence ps ON p.id = ps.protein_id
        JOIN
            ecod_schema.batch_item bi ON p.id = bi.protein_id
        JOIN
            ecod_schema.batch b ON bi.batch_id = b.id
        WHERE 
            bi.batch_id = %s
    """
    
    try:
        proteins = db.execute_dict_query(query, (args.batch_id,))
        logger.info(f"Found {len(proteins)} proteins in batch {args.batch_id}")
    except Exception as e:
        logger.error(f"Database query failed: {e}")
        sys.exit(1)
    
    # Evaluate confidence for each protein
    high_confidence_proteins = []
    confidence_scores = []
    
    for protein in proteins:
        protein_id = protein['id']
        source_id = protein['source_id']
        sequence_length = protein['sequence_length']
        base_path = protein['base_path']
        
        # Construct paths to BLAST result files
        chain_blast_file = os.path.join(base_path, "chain_blast_results", f"{source_id}.chainwise_blast.xml")
        domain_blast_file = os.path.join(base_path, "domain_blast_results", f"{source_id}.domain_blast.xml")
        
        # Parse BLAST results
        domain_hits, domain_query_len = parse_blast_xml(domain_blast_file)
        chain_hits, chain_query_len = parse_blast_xml(chain_blast_file)
        
        if not domain_hits and not chain_hits:
            logger.debug(f"No BLAST results found for protein {source_id}")
            continue
        
        # Use the query length from BLAST or from database
        query_length = domain_query_len or chain_query_len or sequence_length or 0
        
        # Evaluate confidence
        confidence = evaluate_blast_confidence(domain_hits, chain_hits, query_length)
        confidence_scores.append((source_id, confidence))
        
        if confidence >= args.threshold:
            logger.debug(f"High confidence protein: {source_id} with score {confidence:.2f}")
            high_confidence_proteins.append({
                'id': protein_id,
                'source_id': source_id,
                'confidence': confidence
            })
    
    # Sort by confidence score
    high_confidence_proteins.sort(key=lambda x: x['confidence'], reverse=True)
    
    # Print summary
    logger.info(f"Found {len(high_confidence_proteins)} proteins with confidence >= {args.threshold}")
    logger.info(f"This represents {len(high_confidence_proteins) / len(proteins) * 100:.1f}% of the batch")
    
    # Print score distribution
    confidence_scores.sort(key=lambda x: x[1])
    percentiles = [0, 25, 50, 75, 90, 95, 99, 100]
    for p in percentiles:
        if confidence_scores:
            index = min(len(confidence_scores) - 1, int(p * len(confidence_scores) / 100))
            logger.info(f"{p}th percentile: {confidence_scores[index][1]:.2f}")
    
    # Write high confidence proteins to output file if specified
    if args.output:
        with safe_open(args.output, 'w') as f:
            f.write("protein_id,source_id,confidence\n")
            for protein in high_confidence_proteins:
                f.write(f"{protein['id']},{protein['source_id']},{protein['confidence']:.4f}\n")
        logger.info(f"Wrote {len(high_confidence_proteins)} proteins to {args.output}")
    
    # Print top 10 high confidence proteins
    if high_confidence_proteins:
        logger.info("Top 10 highest confidence proteins:")
        for i, protein in enumerate(high_confidence_proteins[:10]):
            logger.info(f"{i+1}. {protein['source_id']} - Score: {protein['confidence']:.4f}")

if __name__ == "__main__":
    main()