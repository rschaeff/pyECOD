#!/usr/bin/env python
"""
This script identifies proteins with high confidence BLAST matches that can be
routed through the BLAST-only pipeline with the new 0.7 threshold.
"""

import argparse
import logging
import sys
import os
from collections import defaultdict

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

# Import pyECOD modules
from ecod.core.application_context import ApplicationContext
from ecod.core.db.database import Database
from ecod.core.blast.blast_parser import BlastResultParser
from ecod.utils.logging_utils import setup_logging

def evaluate_blast_confidence(blast_hits, query_length):
    """
    Evaluate confidence in BLAST results with the improved algorithm
    
    Args:
        blast_hits: List of BLAST hit dictionaries
        query_length: Length of the query protein sequence
    
    Returns:
        float: Confidence score between 0 and 1
    """
    if not blast_hits or len(blast_hits) == 0:
        return 0.0
    
    # Factors for confidence calculation
    top_hits = blast_hits[:5] if len(blast_hits) >= 5 else blast_hits
    
    # Weight factors
    weights = {
        'evalue': 0.3,         # Lower E-value is better
        'coverage': 0.4,        # Higher coverage is better
        'identity': 0.2,        # Higher identity is better
        'consistency': 0.1      # More consistent hits are better
    }
    
    # Calculate E-value score (logarithmic scale, bounded between 0-1)
    min_evalue = min(hit.get('evalue', 1.0) for hit in top_hits)
    evalue_score = max(0, min(1, 1 - (min(100, -1 * (min_evalue + 1e-300).log10()) / 100)))
    
    # Calculate coverage score
    total_coverage = 0
    for hit in top_hits:
        query_start = hit.get('query_start', 0)
        query_end = hit.get('query_end', 0)
        coverage = (query_end - query_start + 1) / query_length if query_length > 0 else 0
        total_coverage += coverage
    coverage_score = min(1, total_coverage / len(top_hits))
    
    # Calculate identity score
    identity_score = sum(hit.get('identity', 0) / 100 for hit in top_hits) / len(top_hits)
    
    # Calculate consistency score (agreement among top hits)
    hit_domains = [hit.get('subject_id', '').split('|')[0] for hit in top_hits]
    domain_counts = defaultdict(int)
    for domain in hit_domains:
        domain_counts[domain] += 1
    
    most_common_domain_count = max(domain_counts.values()) if domain_counts else 0
    consistency_score = most_common_domain_count / len(top_hits) if top_hits else 0
    
    # Weighted sum of all factors
    confidence = (
        weights['evalue'] * evalue_score +
        weights['coverage'] * coverage_score +
        weights['identity'] * identity_score +
        weights['consistency'] * consistency_score
    )
    
    return confidence

def get_proteins_with_blast_results(db, batch_id):
    """
    Get all proteins in a batch that have BLAST results
    
    Args:
        db: Database connection
        batch_id: Batch ID to check
    
    Returns:
        list: List of protein dictionaries with IDs and sequences
    """
    query = """
        SELECT p.id, p.source_id, ps.sequence, ps.sequence_length 
        FROM ecod_schema.protein p
        JOIN ecod_schema.protein_sequence ps ON p.id = ps.protein_id
        JOIN ecod_schema.batch_item bi ON p.id = bi.protein_id
        WHERE bi.batch_id = %s
    """
    return db.execute_query(query, (batch_id,))

def get_blast_results(config, protein_id, source_id):
    """
    Get BLAST results for a specific protein
    
    Args:
        config: Application configuration
        protein_id: Protein ID
        source_id: Source ID (PDB ID + chain)
    
    Returns:
        list: List of BLAST hit dictionaries
    """
    blast_dir = os.path.join(config.get('data_dir'), 'blast_results')
    blast_file = os.path.join(blast_dir, f"{source_id}.blast")
    
    if not os.path.exists(blast_file):
        return []
    
    parser = BlastResultParser()
    return parser.parse_file(blast_file)

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
    setup_logging(log_level)
    logger = logging.getLogger('identify_high_confidence')
    
    # Initialize application context
    try:
        app_context = ApplicationContext(args.config)
        config = app_context.config
        db = Database(config)
    except Exception as e:
        logger.error(f"Failed to initialize application: {e}")
        sys.exit(1)
    
    # Get proteins with BLAST results
    logger.info(f"Getting proteins with BLAST results for batch {args.batch_id}")
    proteins = get_proteins_with_blast_results(db, args.batch_id)
    logger.info(f"Found {len(proteins)} proteins in batch {args.batch_id}")
    
    # Evaluate confidence for each protein
    high_confidence_proteins = []
    confidence_scores = []
    
    for protein in proteins:
        protein_id = protein['id']
        source_id = protein['source_id']
        sequence_length = protein['sequence_length']
        
        blast_hits = get_blast_results(config, protein_id, source_id)
        if not blast_hits:
            logger.debug(f"No BLAST results found for protein {source_id}")
            continue
        
        confidence = evaluate_blast_confidence(blast_hits, sequence_length)
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
        index = min(len(confidence_scores) - 1, int(p * len(confidence_scores) / 100))
        logger.info(f"{p}th percentile: {confidence_scores[index][1]:.2f}")
    
    # Write high confidence proteins to output file if specified
    if args.output:
        with open(args.output, 'w') as f:
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