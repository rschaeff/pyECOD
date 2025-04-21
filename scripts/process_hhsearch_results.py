#!/usr/bin/env python3
"""
process_hhsearch_results.py - Process HHSearch results and integrate with BLAST evidence

This script processes HHSearch results for protein chains, converts them to XML,
and integrates them with BLAST evidence to create domain summaries.
"""

import os
import sys
import argparse
import logging
import glob
from typing import Optional

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext
from ecod.exceptions import PipelineError
from ecod.error_handlers import handle_exceptions
from ecod.pipelines.hhsearch.processor import HHSearchProcessor

def setup_logging(verbose: bool = False, log_file: Optional[str] = None):
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

def register_hhr_files(context, batch_id):
    """Register existing HHR files in the database
    
    Args:
        context: Application context
        batch_id: Batch ID
        
    Returns:
        Number of registered files
    """
    logger = logging.getLogger("ecod.hhsearch_processor")
    
    # Get batch info
    batch_query = """
    SELECT id, batch_name, base_path, ref_version 
    FROM ecod_schema.batch 
    WHERE id = %s
    """
    
    batch_info = context.db.execute_dict_query(batch_query, (batch_id,))
    if not batch_info:
        logger.error(f"Batch {batch_id} not found")
        return 0
    
    base_path = batch_info[0]['base_path']
    ref_version = batch_info[0]['ref_version']
    
    # Find all HHR files
    hhsearch_dir = os.path.join(base_path, "hhsearch")
    hhr_pattern = os.path.join(hhsearch_dir, f"*.{ref_version}.hhr")
    
    hhr_files = glob.glob(hhr_pattern)
    logger.info(f"Found {len(hhr_files)} HHR files on disk")
    
    if not hhr_files:
        logger.warning(f"No HHR files found in {hhsearch_dir}")
        return 0
    
    # Get all chains in this batch
    chains_query = """
    SELECT 
        p.id as protein_id,
        p.pdb_id,
        p.chain_id,
        ps.id as process_id
    FROM 
        ecod_schema.process_status ps
    JOIN
        ecod_schema.protein p ON ps.protein_id = p.id
    WHERE 
        ps.batch_id = %s
    """
    
    chains = context.db.execute_dict_query(chains_query, (batch_id,))
    logger.info(f"Found {len(chains)} chains in batch {batch_id}")
    
    # Create a mapping of pdb_chain to process_id
    chain_map = {}
    for chain in chains:
        pdb_chain = f"{chain['pdb_id']}_{chain['chain_id']}"
        chain_map[pdb_chain] = chain
    
    # Register each HHR file
    registered_count = 0
    for hhr_file in hhr_files:
        # Extract PDB and chain ID from filename
        filename = os.path.basename(hhr_file)
        parts = filename.split('.')
        pdb_chain = parts[0]  # Format: pdbid_chain
        
        if pdb_chain not in chain_map:
            logger.warning(f"No matching chain found for {pdb_chain}")
            continue
        
        chain_info = chain_map[pdb_chain]
        process_id = chain_info['process_id']
        
        # Check if already registered
        check_query = """
        SELECT id FROM ecod_schema.process_file
        WHERE process_id = %s AND file_type = 'hhr'
        """
        
        existing = context.db.execute_query(check_query, (process_id,))
        
        if existing:
            # Update existing record
            context.db.update(
                "ecod_schema.process_file",
                {
                    "file_path": hhr_file,
                    "file_exists": True,
                    "file_size": os.path.getsize(hhr_file)
                },
                "id = %s",
                (existing[0][0],)
            )
        else:
            # Insert new record
            context.db.insert(
                "ecod_schema.process_file",
                {
                    "process_id": process_id,
                    "file_type": "hhr",
                    "file_path": hhr_file,
                    "file_exists": True,
                    "file_size": os.path.getsize(hhr_file)
                }
            )
        
        # Update process status
        context.db.update(
            "ecod_schema.process_status",
            {
                "current_stage": "hhsearch_complete",
                "status": "success"
            },
            "id = %s",
            (process_id,)
        )
        
        registered_count += 1
        
        if registered_count % 50 == 0:
            logger.info(f"Registered {registered_count} HHR files so far")
    
    logger.info(f"Successfully registered {registered_count} HHR files in database")
    return registered_count

@handle_exceptions(exit_on_error=True)
def main():
    """Main entry point for HHSearch results processing script"""
    parser = argparse.ArgumentParser(description='Process ECOD HHSearch Results')
    parser.add_argument('--config', type=str, default='config/config.yml', 
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to process')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    parser.add_argument('--register-only', action='store_true',
                      help='Only register HHR files, do not process them')
    parser.add_argument('--force', action='store_true',
                      help='Force reprocessing of already processed results')

    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    
    # Initialize application context
    context = ApplicationContext(args.config)
    logger = logging.getLogger("ecod.hhsearch_processor")
    
    # Check if batch exists
    query = "SELECT id, batch_name FROM ecod_schema.batch WHERE id = %s"
    result = context.db.execute_dict_query(query, (args.batch_id,))
    
    if not result:
        logger.error(f"Batch {args.batch_id} not found")
        return 1
    
    batch_name = result[0]['batch_name']
    logger.info(f"Processing batch {args.batch_id} ({batch_name})")
    
    # Register HHR files
    registered = register_hhr_files(context, args.batch_id)
    logger.info(f"Registered {registered} HHR files in database")
    
    if args.register_only:
        logger.info("Skipping processing as --register-only was specified")
        return 0
    
    # Process HHSearch results
    try:
        logger.info(f"Starting HHSearch results processing for batch {args.batch_id}")
        
        processor = HHSearchProcessor(context)
        processed_count = processor.process_batch(args.batch_id, args.force)
        
        if processed_count > 0:
            logger.info(f"Successfully processed HHSearch results for {processed_count} chains in batch {args.batch_id}")
            return 0
        else:
            logger.warning(f"No chains were processed in batch {args.batch_id}")
            return 0
    
    except PipelineError as e:
        logger.error(f"Pipeline error: {str(e)}")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error: {str(e)}", exc_info=True)
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())#!/usr/bin/env python3
"""
process_hhsearch_results.py - Process HHSearch results and integrate with BLAST evidence

This script processes HHSearch results for protein chains, converts them to XML,
and integrates them with BLAST evidence to create domain summaries.
"""

import os
import sys
import argparse
import logging
import glob
from typing import Optional

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext
from ecod.exceptions import PipelineError
from ecod.error_handlers import handle_exceptions
from ecod.pipelines.hhsearch.processor import HHSearchProcessor

def setup_logging(verbose: bool = False, log_file: Optional[str] = None):
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

def register_hhr_files(context, batch_id):
    """Register existing HHR files in the database
    
    Args:
        context: Application context
        batch_id: Batch ID
        
    Returns:
        Number of registered files
    """
    logger = logging.getLogger("ecod.hhsearch_processor")
    
    # Get batch info
    batch_query = """
    SELECT id, batch_name, base_path, ref_version 
    FROM ecod_schema.batch 
    WHERE id = %s
    """
    
    batch_info = context.db.execute_dict_query(batch_query, (batch_id,))
    if not batch_info:
        logger.error(f"Batch {batch_id} not found")
        return 0
    
    base_path = batch_info[0]['base_path']
    ref_version = batch_info[0]['ref_version']
    
    # Find all HHR files
    hhsearch_dir = os.path.join(base_path, "hhsearch")
    hhr_pattern = os.path.join(hhsearch_dir, f"*.{ref_version}.hhr")
    
    hhr_files = glob.glob(hhr_pattern)
    logger.info(f"Found {len(hhr_files)} HHR files on disk")
    
    if not hhr_files:
        logger.warning(f"No HHR files found in {hhsearch_dir}")
        return 0
    
    # Get all chains in this batch
    chains_query = """
    SELECT 
        p.id as protein_id,
        p.pdb_id,
        p.chain_id,
        ps.id as process_id
    FROM 
        ecod_schema.process_status ps
    JOIN
        ecod_schema.protein p ON ps.protein_id = p.id
    WHERE 
        ps.batch_id = %s
    """
    
    chains = context.db.execute_dict_query(chains_query, (batch_id,))
    logger.info(f"Found {len(chains)} chains in batch {batch_id}")
    
    # Create a mapping of pdb_chain to process_id
    chain_map = {}
    for chain in chains:
        pdb_chain = f"{chain['pdb_id']}_{chain['chain_id']}"
        chain_map[pdb_chain] = chain
    
    # Register each HHR file
    registered_count = 0
    for hhr_file in hhr_files:
        # Extract PDB and chain ID from filename
        filename = os.path.basename(hhr_file)
        parts = filename.split('.')
        pdb_chain = parts[0]  # Format: pdbid_chain
        
        if pdb_chain not in chain_map:
            logger.warning(f"No matching chain found for {pdb_chain}")
            continue
        
        chain_info = chain_map[pdb_chain]
        process_id = chain_info['process_id']
        
        # Check if already registered
        check_query = """
        SELECT id FROM ecod_schema.process_file
        WHERE process_id = %s AND file_type = 'hhr'
        """
        
        existing = context.db.execute_query(check_query, (process_id,))
        
        if existing:
            # Update existing record
            context.db.update(
                "ecod_schema.process_file",
                {
                    "file_path": hhr_file,
                    "file_exists": True,
                    "file_size": os.path.getsize(hhr_file)
                },
                "id = %s",
                (existing[0][0],)
            )
        else:
            # Insert new record
            context.db.insert(
                "ecod_schema.process_file",
                {
                    "process_id": process_id,
                    "file_type": "hhr",
                    "file_path": hhr_file,
                    "file_exists": True,
                    "file_size": os.path.getsize(hhr_file)
                }
            )
        
        # Update process status
        context.db.update(
            "ecod_schema.process_status",
            {
                "current_stage": "hhsearch_complete",
                "status": "success"
            },
            "id = %s",
            (process_id,)
        )
        
        registered_count += 1
        
        if registered_count % 50 == 0:
            logger.info(f"Registered {registered_count} HHR files so far")
    
    logger.info(f"Successfully registered {registered_count} HHR files in database")
    return registered_count

@handle_exceptions(exit_on_error=True)
def main():
    """Main entry point for HHSearch results processing script"""
    parser = argparse.ArgumentParser(description='Process ECOD HHSearch Results')
    parser.add_argument('--config', type=str, default='config/config.yml', 
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to process')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    parser.add_argument('--register-only', action='store_true',
                      help='Only register HHR files, do not process them')
    parser.add_argument('--force', action='store_true',
                      help='Force reprocessing of already processed results')

    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    
    # Initialize application context
    context = ApplicationContext(args.config)
    logger = logging.getLogger("ecod.hhsearch_processor")
    
    # Check if batch exists
    query = "SELECT id, batch_name FROM ecod_schema.batch WHERE id = %s"
    result = context.db.execute_dict_query(query, (args.batch_id,))
    
    if not result:
        logger.error(f"Batch {args.batch_id} not found")
        return 1
    
    batch_name = result[0]['batch_name']
    logger.info(f"Processing batch {args.batch_id} ({batch_name})")
    
    # Register HHR files
    registered = register_hhr_files(context, args.batch_id)
    logger.info(f"Registered {registered} HHR files in database")
    
    if args.register_only:
        logger.info("Skipping processing as --register-only was specified")
        return 0
    
    # Process HHSearch results
    try:
        logger.info(f"Starting HHSearch results processing for batch {args.batch_id}")
        
        processor = HHSearchProcessor(context)
        processed_count = processor.process_batch(args.batch_id, args.force)
        
        if processed_count > 0:
            logger.info(f"Successfully processed HHSearch results for {processed_count} chains in batch {args.batch_id}")
            return 0
        else:
            logger.warning(f"No chains were processed in batch {args.batch_id}")
            return 0
    
    except PipelineError as e:
        logger.error(f"Pipeline error: {str(e)}")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error: {str(e)}", exc_info=True)
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())