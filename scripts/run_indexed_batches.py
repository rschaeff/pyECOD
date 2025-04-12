#!/usr/bin/env python3
"""
run_indexed_batches.py - Run domain summary generation on indexed batches
"""

import os
import sys
import logging
import argparse
import subprocess
from datetime import datetime
from typing import List, Dict, Any

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

def setup_logging(verbose: bool = False, log_file: str = None):
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

def get_indexed_batches(config_path: str) -> List[Dict[str, Any]]:
    """Get list of indexed batches directly from the database"""
    logger = logging.getLogger("ecod.run_summaries")
    
    # Import db connection here to avoid circular imports
    sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
    from ecod.core.context import ApplicationContext
    
    context = ApplicationContext(config_path)
    
    query = """
    SELECT 
        id, batch_name, base_path, ref_version, total_items
    FROM 
        ecod_schema.batch
    WHERE 
        total_items = 5000 AND status = 'indexed'
    ORDER BY 
        id
    """
    
    result = context.db.execute_query(query)
    
    batches = []
    for row in result:
        batches.append({
            'id': row[0],
            'name': row[1],
            'path': row[2],
            'reference': row[3],
            'total_items': row[4]
        })
    
    logger.info(f"Found {len(batches)} indexed batches")
    
    return batches

def check_batch_summaries(batch_id: int, config_path: str) -> Dict[str, int]:
    """Check how many domain summaries exist for this batch"""
    logger = logging.getLogger("ecod.run_summaries")
    
    # Import db connection here to avoid circular imports
    sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
    from ecod.core.context import ApplicationContext
    
    context = ApplicationContext(config_path)
    
    query = """
    SELECT 
        COUNT(DISTINCT pf.id) as summary_count
    FROM 
        ecod_schema.process_status ps
    JOIN 
        ecod_schema.batch b ON ps.batch_id = b.id
    LEFT JOIN 
        ecod_schema.process_file pf ON ps.id = pf.process_id AND pf.file_type = 'domain_summary' AND pf.file_exists = TRUE
    WHERE 
        b.id = %s
    """
    
    result = context.db.execute_query(query, (batch_id,))
    
    summary_count = result[0][0] if result else 0
    
    return {
        'batch_id': batch_id,
        'summary_count': summary_count
    }

def run_domain_summary(batch_id: int, config_path: str, output_dir: str = None, 
                       blast_only: bool = True, threads: int = 8, log_dir: str = "logs",
                       force: bool = False, limit: int = None):
    """Run domain summary generation for a batch"""
    logger = logging.getLogger("ecod.run_summaries")
    
    # Ensure log directory exists
    os.makedirs(log_dir, exist_ok=True)
    
    # Build command
    cmd = [
        "python", "scripts/generate_batch_domain_summaries.py",
        "--config", config_path,
        "--batch-id", str(batch_id),
        "--threads", str(threads)
    ]
    
    if blast_only:
        cmd.append("--blast-only")
    
    if output_dir:
        cmd.extend(["--output-dir", output_dir])
    
    if force:
        cmd.append("--force")
    
    if limit:
        cmd.extend(["--limit", str(limit)])
    
    # Add log file
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f"domain_summary_batch_{batch_id}_{timestamp}.log")
    cmd.extend(["--log-file", log_file])
    
    # Execute command
    logger.info(f"Running domain summary for batch {batch_id}")
    logger.debug(f"Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd, 
            check=True, 
            capture_output=True, 
            text=True
        )
        logger.info(f"Domain summary completed for batch {batch_id}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running domain summary for batch {batch_id}: {e}")
        logger.error(f"stderr: {e.stderr}")
        return False

def main():
    """Main function to run domain summary generation on indexed batches"""
    parser = argparse.ArgumentParser(description='Run domain summary generation on indexed batches')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--output-dir', type=str,
                      help='Override output directory (default: from batch path)')
    parser.add_argument('--log-dir', type=str, default='logs',
                      help='Directory for log files')
    parser.add_argument('--batch-id', type=int,
                      help='Process specific batch ID only')
    parser.add_argument('--start-id', type=int,
                      help='Start processing from this batch ID')
    parser.add_argument('--end-id', type=int,
                      help='Stop processing at this batch ID')
    parser.add_argument('--blast-only', action='store_true', default=True,
                      help='Generate blast-only summaries (no HHSearch)')
    parser.add_argument('--threads', type=int, default=8,
                      help='Number of threads for parallel processing')
    parser.add_argument('--dry-run', action='store_true',
                      help='Show batches but do not run processing')
    parser.add_argument('--force', action='store_true',
                      help='Force regeneration of existing summaries')
    parser.add_argument('--limit', type=int,
                      help='Limit number of proteins per batch to process')
    parser.add_argument('--skip-complete', action='store_true',
                      help='Skip batches that already have all summaries')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    
    log_file = os.path.join(args.log_dir, f'run_indexed_batches_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
    setup_logging(args.verbose, log_file)
    logger = logging.getLogger("ecod.run_summaries")
    
    # Get indexed batches
    if args.batch_id:
        batches = [{'id': args.batch_id}]
        logger.info(f"Processing single batch ID: {args.batch_id}")
    else:
        batches = get_indexed_batches(args.config)
        
        # Filter by start/end ID if specified
        if args.start_id:
            batches = [b for b in batches if b['id'] >= args.start_id]
        if args.end_id:
            batches = [b for b in batches if b['id'] <= args.end_id]
        
    if not batches:
        logger.warning(f"No indexed batches found to process")
        return 0
    
    # Check current summary status for each batch
    for i, batch in enumerate(batches):
        batch_id = batch['id']
        summary_status = check_batch_summaries(batch_id, args.config)
        batch['summary_count'] = summary_status['summary_count']
        
        # Skip batches that already have all summaries
        if args.skip_complete and batch.get('total_items', 5000) == batch['summary_count']:
            logger.info(f"Batch {batch_id} already has all summaries ({batch['summary_count']}/{batch.get('total_items', 5000)})")
            batch['skip'] = True
        else:
            batch['skip'] = False
    
    # Display batch info
    logger.info(f"Found {len(batches)} batches to process:")
    for i, batch in enumerate(batches):
        batch_id = batch['id']
        if 'name' in batch:
            skip_status = " [SKIP - COMPLETE]" if batch.get('skip', False) else ""
            logger.info(f"  {i+1}. Batch ID: {batch_id} - {batch['name']}{skip_status}")
            logger.info(f"     Summary files: {batch['summary_count']}/{batch.get('total_items', 5000)}")
        else:
            logger.info(f"  {i+1}. Batch ID: {batch_id}")
    
    if args.dry_run:
        logger.info("Dry run completed. No processing performed.")
        return 0
    
    # Process each batch
    success_count = 0
    for batch in batches:
        batch_id = batch['id']
        
        # Skip if marked for skipping
        if batch.get('skip', False):
            logger.info(f"Skipping batch {batch_id} (already complete)")
            continue
        
        success = run_domain_summary(
            batch_id=batch_id, 
            config_path=args.config, 
            output_dir=args.output_dir,
            blast_only=args.blast_only,
            threads=args.threads,
            log_dir=args.log_dir,
            force=args.force,
            limit=args.limit
        )
        
        if success:
            success_count += 1
    
    logger.info(f"Processing complete: {success_count} of {len([b for b in batches if not b.get('skip', False)])} batches successful")
    
    return 0 if success_count == len([b for b in batches if not b.get('skip', False)]) else 1

if __name__ == "__main__":
    sys.exit(main())