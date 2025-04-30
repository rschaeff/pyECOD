#!/usr/bin/env python3
"""
register_hhsearch.py - Register HHR files and convert to XML format

This script provides a command-line interface to the HHResultRegistrar module,
which handles the registration and conversion of HHSearch result files for
the ECOD domain analysis pipeline.

Usage:
    python register_hhsearch.py --batch-id <batch_id> [--force] [--chains <pdb_id>_<chain_id> ...]
"""

import os
import sys
import logging
import argparse
from typing import List, Dict, Any, Optional, Tuple

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

# Import from pipeline module for consistency
from ecod.core.context import ApplicationContext
from ecod.pipelines.domain_analysis.hhresult_registrar import HHResultRegistrar


def setup_logging(verbose=False, log_file=None):
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


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Register and convert HHSearch results')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to process')
    parser.add_argument('--chains', nargs='+', default=None,
                      help='Specific chains to process (format: pdbid_chainid)')
    parser.add_argument('--force', action='store_true',
                      help='Force reprocessing of already processed results')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')

    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)

    logger = logging.getLogger("ecod.register_hhsearch")
    logger.info(f"Starting HHSearch registration and conversion for batch {args.batch_id}")

    # Use absolute path for config if provided
    if args.config and not os.path.isabs(args.config):
        args.config = os.path.abspath(args.config)

    # Initialize application context
    context = ApplicationContext(args.config)

    # Create registrar from the pipeline module
    registrar = HHResultRegistrar(context)

    try:
        # Process batch
        if args.chains:
            result = registrar.register_specific_chains(
                args.batch_id,
                args.chains,
                args.force
            )
        else:
            result = registrar.register_batch_results(
                args.batch_id,
                args.force
            )

        logger.info(f"Successfully registered {result} HHR files and converted them to XML")

        if result == 0:
            return 1

        return 0
    except Exception as e:
        logger.error(f"Error processing batch: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
