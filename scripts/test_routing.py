#!/usr/bin/env python3
"""
Test routing functionality for ECOD pipeline.
This script only tests the ProcessingRouter without reprocessing BLAST results.
"""

import argparse
import logging
import sys
from typing import Dict, Any, List

from ecod.core.context import ApplicationContext
from ecod.pipelines.domain_analysis.routing import ProcessingRouter
from ecod.exceptions import PipelineError


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


def test_routing(batch_id: int, config_path: str, verbose: bool = False, log_file: str = None):
    """Test routing for a specific batch"""
    # Setup logging
    setup_logging(verbose, log_file)
    logger = logging.getLogger("route_test")
    
    try:
        # Initialize application context
        logger.info(f"Initializing application context with config: {config_path}")
        context = ApplicationContext(config_path)
        
        # Initialize router
        logger.info("Initializing ProcessingRouter")
        router = ProcessingRouter(context)
        
        # Run routing
        logger.info(f"Running routing for batch {batch_id}")
        paths = router.assign_processing_paths(batch_id)
        
        # Display results
        path_counts = {k: len(v) for k, v in paths.items()}
        logger.info(f"Routing complete. Path assignments: {path_counts}")
        
        for path_type, protein_ids in paths.items():
            logger.info(f"{path_type} proteins: {len(protein_ids)}")
            if verbose and protein_ids:
                logger.debug(f"Sample protein IDs: {protein_ids[:5]}...")
        
        # Count proteins above threshold
        if 'blast_only' in paths:
            blast_only_count = len(paths['blast_only'])
            logger.info(f"BLAST-only proteins (confidence >= {router.confidence_threshold}): {blast_only_count}")
        
        if 'full_pipeline' in paths:
            full_pipeline_count = len(paths['full_pipeline'])
            logger.info(f"Full pipeline proteins (confidence < {router.confidence_threshold}): {full_pipeline_count}")
        
        return 0
        
    except Exception as e:
        logger.error(f"Error in routing test: {str(e)}", exc_info=True)
        return 1


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description="Test ECOD pipeline routing")
    parser.add_argument("--batch-id", type=int, required=True, help="Batch ID to process")
    parser.add_argument("--config", type=str, default="config/config.yml", help="Path to configuration file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    parser.add_argument("--log-file", type=str, help="Log to file in addition to stdout")
    
    args = parser.parse_args()
    
    return test_routing(args.batch_id, args.config, args.verbose, args.log_file)


if __name__ == "__main__":
    sys.exit(main())