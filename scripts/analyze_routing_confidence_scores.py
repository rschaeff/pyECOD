#!/usr/bin/env python3
"""
Analyze confidence score distribution for ECOD pipeline routing.
"""

import argparse
import logging
import sys, os
from typing import Dict, Any, List, Tuple
import matplotlib.pyplot as plt
import numpy as np

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
from ecod.core.context import ApplicationContext
from ecod.pipelines.domain_analysis.routing import ProcessingRouter
from ecod.db.manager import DBManager


def setup_logging(verbose: bool = False, log_file: str = None):
    """Set up logging configuration"""
    log_level = logging.DEBUG if verbose else logging.INFO
    
    handlers = [logging.StreamHandler()]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )


def analyze_confidence_distribution(batch_id: int, config_path: str, verbose: bool = False, log_file: str = None):
    """Analyze confidence score distribution for a batch"""
    # Setup logging
    setup_logging(verbose, log_file)
    logger = logging.getLogger("confidence_analysis")
    
    try:
        # Initialize application context
        logger.info(f"Initializing application context with config: {config_path}")
        context = ApplicationContext(config_path)
        
        # Initialize router
        logger.info("Initializing ProcessingRouter")
        router = ProcessingRouter(context)
        
        # Get database connection
        db = context.db
        
        # First check if confidence metrics exist in the database
        logger.info(f"Checking for existing confidence metrics for batch {batch_id}")
        existing_metrics = _get_existing_metrics(db, batch_id)
        
        if existing_metrics:
            logger.info(f"Found {len(existing_metrics)} existing confidence metrics")
            scores, components = existing_metrics
        else:
            logger.info("No existing metrics found, calculating confidence scores")
            # Calculate scores for all proteins in the batch
            scores, components = _calculate_confidence_scores(db, router, batch_id)
        
        if not scores:
            logger.error("No confidence scores could be calculated or retrieved")
            return 1
        
        # Analyze the scores
        _analyze_scores(scores, components, router.confidence_threshold, logger)
        
        return 0
        
    except Exception as e:
        logger.error(f"Error in confidence analysis: {str(e)}", exc_info=True)
        return 1


def _get_existing_metrics(db: DBManager, batch_id: int) -> Tuple[List[float], Dict[str, List[float]]]:
    """Get existing confidence metrics from database"""
    # Query to get proteins in the batch
    query = """
    SELECT protein_id 
    FROM ecod_schema.process_status
    WHERE batch_id = %s
    """
    
    try:
        protein_rows = db.execute_query(query, (batch_id,))
        protein_ids = [row[0] for row in protein_rows]
        
        if not protein_ids:
            return [], {}
        
        # Query for confidence metrics
        metrics_query = """
        SELECT 
            protein_id, 
            min_evalue, 
            max_coverage, 
            max_identity, 
            hit_consistency, 
            architecture_coverage, 
            overall_confidence
        FROM ecod_schema.blast_confidence_metrics
        WHERE protein_id IN %s
        """
        
        metrics_rows = db.execute_dict_query(metrics_query, (tuple(protein_ids),))
        
        if not metrics_rows:
            return [], {}
        
        # Organize the data
        scores = []
        components = {
            'evalue': [],
            'coverage': [],
            'identity': [],
            'consistency': [],
            'architecture': []
        }
        
        for row in metrics_rows:
            scores.append(float(row['overall_confidence']))
            
            # Store component values if they exist
            if 'min_evalue' in row:
                components['evalue'].append(float(row['min_evalue']))
            if 'max_coverage' in row:
                components['coverage'].append(float(row['max_coverage']))
            if 'max_identity' in row:
                components['identity'].append(float(row['max_identity']))
            if 'hit_consistency' in row:
                components['consistency'].append(float(row['hit_consistency']))
            if 'architecture_coverage' in row:
                components['architecture'].append(float(row['architecture_coverage']))
        
        return scores, components
        
    except Exception as e:
        logging.error(f"Error retrieving metrics: {e}")
        return [], {}


def _calculate_confidence_scores(db: DBManager, router: ProcessingRouter, batch_id: int) -> Tuple[List[float], Dict[str, List[float]]]:
    """Calculate confidence scores for all proteins in a batch"""
    # Query to get proteins in the batch
    query = """
    SELECT protein_id 
    FROM ecod_schema.process_status
    WHERE batch_id = %s
    """
    
    try:
        protein_rows = db.execute_query(query, (batch_id,))
        protein_ids = [row[0] for row in protein_rows]
        
        if not protein_ids:
            return [], {}
        
        # Calculate confidence for each protein
        scores = []
        components = {
            'evalue': [],
            'coverage': [],
            'identity': [],
            'consistency': [],
            'architecture': []
        }
        
        logger = logging.getLogger("confidence_analysis")
        
        for i, protein_id in enumerate(protein_ids):
            if i % 100 == 0:
                logger.info(f"Calculating confidence for protein {i+1}/{len(protein_ids)}")
            
            # Get BLAST results
            domain_results = router._get_domain_blast_results(protein_id)
            chain_results = router._get_chain_blast_results(protein_id)
            
            # Calculate metrics
            domain_metrics = router._calculate_domain_blast_metrics(domain_results)
            chain_metrics = router._calculate_chain_blast_metrics(chain_results)
            
            # Combine metrics
            min_evalue = min(
                domain_metrics.get('min_evalue', 10.0),
                chain_metrics.get('min_evalue', 10.0)
            )
            max_coverage = max(
                domain_metrics.get('max_coverage', 0.0),
                chain_metrics.get('max_coverage', 0.0)
            )
            max_identity = max(
                domain_metrics.get('max_identity', 0.0),
                chain_metrics.get('max_identity', 0.0)
            )
            hit_consistency = (
                domain_metrics.get('hit_consistency', 0.0) * 0.7 +
                chain_metrics.get('hit_consistency', 0.0) * 0.3
            )
            
            # Calculate domain architecture coverage
            architecture_coverage = router._calculate_architecture_coverage(
                protein_id, domain_results
            )
            
            # Calculate confidence
            evalue_score = router._evalue_to_confidence(min_evalue)
            confidence = (
                0.3 * evalue_score + 
                0.25 * max_coverage + 
                0.2 * max_identity + 
                0.15 * hit_consistency +
                0.1 * architecture_coverage
            )
            
            scores.append(confidence)
            
            # Save component values
            components['evalue'].append(min_evalue)
            components['coverage'].append(max_coverage)
            components['identity'].append(max_identity)
            components['consistency'].append(hit_consistency)
            components['architecture'].append(architecture_coverage)
        
        return scores, components
        
    except Exception as e:
        logging.error(f"Error calculating confidence scores: {e}")
        return [], {}


def _analyze_scores(scores: List[float], components: Dict[str, List[float]], threshold: float, logger: logging.Logger):
    """Analyze confidence scores and display distribution"""
    # Basic statistics
    scores_array = np.array(scores)
    mean_score = np.mean(scores_array)
    median_score = np.median(scores_array)
    min_score = np.min(scores_array)
    max_score = np.max(scores_array)
    std_dev = np.std(scores_array)
    
    # Count proteins above threshold
    above_threshold = np.sum(scores_array >= threshold)
    below_threshold = len(scores_array) - above_threshold
    
    # Log analysis results
    logger.info(f"Confidence score analysis:")
    logger.info(f"  Total proteins: {len(scores)}")
    logger.info(f"  Mean score: {mean_score:.4f}")
    logger.info(f"  Median score: {median_score:.4f}")
    logger.info(f"  Range: {min_score:.4f} - {max_score:.4f}")
    logger.info(f"  Standard deviation: {std_dev:.4f}")
    logger.info(f"  Above threshold ({threshold}): {above_threshold} ({above_threshold/len(scores)*100:.1f}%)")
    logger.info(f"  Below threshold ({threshold}): {below_threshold} ({below_threshold/len(scores)*100:.1f}%)")
    
    # Print values at different percentiles
    for percentile in [10, 25, 50, 75, 90, 95, 99]:
        value = np.percentile(scores_array, percentile)
        logger.info(f"  {percentile}th percentile: {value:.4f}")
        
    # Display recommended thresholds based on distribution
    logger.info("Recommended thresholds:")
    for target_percent in [10, 20, 30, 40, 50]:
        target_threshold = np.percentile(scores_array, 100 - target_percent)
        logger.info(f"  To route {target_percent}% to BLAST-only: {target_threshold:.4f}")
    
    # Output to CSV file
    try:
        with open('confidence_analysis.csv', 'w') as f:
            f.write("Score,E-value,Coverage,Identity,Consistency,Architecture\n")
            for i in range(len(scores)):
                row = [
                    str(scores[i]),
                    str(components['evalue'][i]) if i < len(components['evalue']) else '',
                    str(components['coverage'][i]) if i < len(components['coverage']) else '',
                    str(components['identity'][i]) if i < len(components['identity']) else '',
                    str(components['consistency'][i]) if i < len(components['consistency']) else '',
                    str(components['architecture'][i]) if i < len(components['architecture']) else ''
                ]
                f.write(','.join(row) + '\n')
        logger.info("Wrote detailed data to confidence_analysis.csv")
    except Exception as e:
        logger.error(f"Error writing CSV file: {e}")
    
    # Generate histogram
    try:
        plt.figure(figsize=(10, 6))
        plt.hist(scores_array, bins=20, alpha=0.7, color='blue')
        plt.axvline(x=threshold, color='red', linestyle='--', label=f'Threshold ({threshold})')
        plt.axvline(x=mean_score, color='green', linestyle='-', label=f'Mean ({mean_score:.4f})')
        plt.axvline(x=median_score, color='orange', linestyle='-', label=f'Median ({median_score:.4f})')
        plt.title('Confidence Score Distribution')
        plt.xlabel('Confidence Score')
        plt.ylabel('Number of Proteins')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig('confidence_histogram.png')
        logger.info("Generated histogram in confidence_histogram.png")
    except Exception as e:
        logger.error(f"Error generating histogram: {e}")


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description="Analyze confidence score distribution")
    parser.add_argument("--batch-id", type=int, required=True, help="Batch ID to analyze")
    parser.add_argument("--config", type=str, default="config/config.yml", help="Path to configuration file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    parser.add_argument("--log-file", type=str, help="Log to file in addition to stdout")
    
    args = parser.parse_args()
    
    return analyze_confidence_distribution(args.batch_id, args.config, args.verbose, args.log_file)


if __name__ == "__main__":
    sys.exit(main())