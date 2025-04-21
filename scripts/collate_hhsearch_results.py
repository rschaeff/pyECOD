# Script to collate HHSearch results for batch 31 (representative processes only)
import os
import sys
import logging
import argparse

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))w
from ecod.core.context import ApplicationContext
from ecod.pipelines.domain_analysis.summary import DomainSummary

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

def collate_batch_results(batch_id, force=False):
    """Collate HHSearch results with BLAST results for representative processes in a batch"""
    logger = logging.getLogger("collator")
    
    # Initialize context and components
    context = ApplicationContext()
    domain_summary = DomainSummary(context)
    
    # Get batch info
    batch_query = """
    SELECT id, batch_name, base_path, ref_version 
    FROM ecod_schema.batch 
    WHERE id = %s
    """
    
    batch_info = context.db_manager.execute_dict_query(batch_query, (batch_id,))
    if not batch_info:
        logger.error(f"Batch {batch_id} not found")
        return False
    
    base_path = batch_info[0]['base_path']
    ref_version = batch_info[0]['ref_version']
    
    # Get representative proteins with HHSearch results
    protein_query = """
    SELECT 
        p.id, p.pdb_id, p.chain_id, ps.id as process_id
    FROM 
        ecod_schema.protein p
    JOIN
        ecod_schema.process_status ps ON p.id = ps.protein_id
    JOIN
        ecod_schema.process_file pf ON ps.id = pf.process_id
    WHERE 
        ps.batch_id = %s
        AND ps.is_representative = TRUE
        AND pf.file_type = 'hhr'
        AND pf.file_exists = TRUE
    """
    
    proteins = context.db_manager.execute_dict_query(protein_query, (batch_id,))
    logger.info(f"Found {len(proteins)} representative proteins with HHSearch results")
    
    success_count = 0
    for protein in proteins:
        pdb_id = protein['pdb_id']
        chain_id = protein['chain_id']
        process_id = protein['process_id']
        
        # Create domain summary with HHSearch evidence
        try:
            # Call DomainSummary's create_summary with blast_only=False
            summary_file = domain_summary.create_summary(
                pdb_id, 
                chain_id, 
                ref_version, 
                base_path, 
                blast_only=False  # Use HHSearch results
            )
            
            if summary_file:
                # Register summary file in database
                register_summary(context, process_id, summary_file, base_path)
                success_count += 1
                logger.info(f"Successfully created full domain summary for {pdb_id}_{chain_id}")
            else:
                logger.warning(f"Failed to create domain summary for {pdb_id}_{chain_id}")
        
        except Exception as e:
            logger.error(f"Error processing {pdb_id}_{chain_id}: {str(e)}")
    
    logger.info(f"Successfully collated results for {success_count}/{len(proteins)} proteins")
    return success_count > 0

def register_summary(context, process_id, summary_file, base_path):
    """Register domain summary in database"""
    logger = logging.getLogger("registrar")
    
    try:
        relative_path = os.path.relpath(summary_file, base_path)
        file_size = os.path.getsize(summary_file)
        
        # Check if summary already registered
        check_query = """
        SELECT id FROM ecod_schema.process_file
        WHERE process_id = %s AND file_type = 'domain_summary'
        """
        
        existing = context.db_manager.execute_query(check_query, (process_id,))
        
        if existing:
            # Update existing record
            context.db_manager.update(
                "ecod_schema.process_file",
                {
                    "file_path": relative_path,
                    "file_exists": True,
                    "file_size": file_size
                },
                "id = %s",
                (existing[0][0],)
            )
        else:
            # Insert new record
            context.db_manager.insert(
                "ecod_schema.process_file",
                {
                    "process_id": process_id,
                    "file_type": "domain_summary",
                    "file_path": relative_path,
                    "file_exists": True,
                    "file_size": file_size
                }
            )
        
        # Update process status
        context.db_manager.update(
            "ecod_schema.process_status",
            {
                "current_stage": "domain_summary_complete",
                "status": "success"
            },
            "id = %s",
            (process_id,)
        )
        return True
    except Exception as e:
        logger.error(f"Error registering summary: {str(e)}")
        return False

def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Collate ECOD HHSearch Results with BLAST')
    parser.add_argument('--batch-id', type=int, default=31,
                      help='Batch ID to process (default: 31)')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    parser.add_argument('--force', action='store_true',
                      help='Force reprocessing of already processed results')

    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    
    logger = logging.getLogger("main")
    logger.info(f"Starting collation of BLAST and HHSearch results for batch {args.batch_id}")
    
    success = collate_batch_results(args.batch_id, args.force)
    
    if success:
        logger.info(f"Successfully collated results for batch {args.batch_id}")
        return 0
    else:
        logger.error(f"Failed to collate results for batch {args.batch_id}")
        return 1

if __name__ == "__main__":
    sys.exit(main())