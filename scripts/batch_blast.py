#!/usr/bin/env python3
import os
import sys
import argparse
import logging
from datetime import datetime
from typing import List, Dict, Any, Optional

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
# Import required modules from pyECOD
from ecod.core.context import ApplicationContext
from ecod.core.logging_config import LoggingManager
from ecod.pipelines.blast_pipeline import BlastPipeline
from ecod.exceptions import ConfigurationError, PipelineError

def setup_logging(verbose: bool, log_file: Optional[str]) -> logging.Logger:
    """Set up logging configuration"""
    return LoggingManager.configure(
        verbose=verbose,
        log_file=log_file,
        component="batch_blast",
    )

def find_proteins_with_db_fs_mismatch(context: ApplicationContext, batch_id: int) -> List[Dict]:
    """
    Find proteins with mismatches between database records and filesystem files.

    Args:
        context: Application context
        batch_id: Batch ID to check

    Returns:
        List of protein dictionaries with mismatches
    """
    logger = logging.getLogger("batch_blast.sync_check")

    # Get batch information
    batch_info = get_batch_statistics(context, batch_id)
    if not batch_info:
        logger.error(f"Batch {batch_id} not found")
        return []

    # Get base path for batch
    query = "SELECT base_path FROM ecod_schema.batch WHERE id = %s"
    result = context.db.execute_query(query, (batch_id,))
    if not result:
        logger.error(f"Batch {batch_id} not found or has no base_path")
        return []

    base_path = result[0][0]

    # Query for process files
    query = """
    SELECT
        p.id as protein_id,
        p.pdb_id,
        p.chain_id,
        ps.id as process_id,
        ps.relative_path,
        pf.file_type,
        pf.file_path,
        pf.file_exists,
        pf.id as file_id
    FROM
        ecod_schema.process_status ps
    JOIN
        ecod_schema.protein p ON ps.protein_id = p.id
    JOIN
        ecod_schema.process_file pf ON ps.id = pf.process_id
    WHERE
        ps.batch_id = %s
        AND pf.file_type IN ('chain_blast_result', 'domain_blast_result')
    ORDER BY
        p.pdb_id, p.chain_id, pf.file_type
    """

    file_records = context.db.execute_dict_query(query, (batch_id,))
    logger.info(f"Found {len(file_records)} BLAST file records")

    # Group by protein
    proteins_by_id = {}
    for record in file_records:
        protein_id = record['protein_id']
        if protein_id not in proteins_by_id:
            proteins_by_id[protein_id] = {
                'protein_id': protein_id,
                'pdb_id': record['pdb_id'],
                'chain_id': record['chain_id'],
                'process_id': record['process_id'],
                'relative_path': record['relative_path'],
                'files': []
            }

        proteins_by_id[protein_id]['files'].append({
            'file_type': record['file_type'],
            'file_path': record['file_path'],
            'file_exists_db': record['file_exists'],
            'file_id': record['file_id']
        })

    # Check filesystem against database records
    mismatched_proteins = []

    for protein_id, protein in proteins_by_id.items():
        has_mismatch = False

        for file_info in protein['files']:
            full_path = os.path.join(base_path, file_info['file_path']) if file_info['file_path'] else None
            file_exists_fs = False

            if full_path and os.path.exists(full_path):
                file_exists_fs = True
                file_info['file_exists_fs'] = True
                file_info['actual_file_size'] = os.path.getsize(full_path)
            else:
                file_info['file_exists_fs'] = False
                file_info['actual_file_size'] = None

            # Check for mismatch
            if file_info['file_exists_db'] != file_exists_fs:
                has_mismatch = True
                file_info['has_mismatch'] = True
            else:
                file_info['has_mismatch'] = False

        if has_mismatch:
            mismatched_proteins.append(protein)

    logger.info(f"Found {len(mismatched_proteins)} proteins with DB-FS mismatches")
    return mismatched_proteins

def find_potential_blast_paths(base_path: str, pdb_id: str, chain_id: str, file_type: str) -> List[str]:
    """
    Find potential paths where BLAST files might exist.

    Args:
        base_path: Base path for the batch
        pdb_id: PDB ID
        chain_id: Chain ID
        file_type: File type (chain_blast_result or domain_blast_result)

    Returns:
        List of potential file paths
    """
    patterns = []

    if file_type == 'chain_blast_result':
        # Common patterns for chain BLAST files
        patterns = [
            os.path.join(base_path, "chain_blast_results", f"{pdb_id}_{chain_id}.chainwise_blast.xml"),
            os.path.join(base_path, "chain_blast_results", f"{pdb_id.lower()}_{chain_id}.chainwise_blast.xml"),
            os.path.join(base_path, "blast", "chain", f"{pdb_id}_{chain_id}.xml"),
            os.path.join(base_path, "blast", "chain", f"{pdb_id}_{chain_id}.blast.xml"),
            # Legacy patterns
            os.path.join(base_path, "blast", f"{pdb_id}_{chain_id}.chainwise_blast.xml")
        ]
    elif file_type == 'domain_blast_result':
        # Common patterns for domain BLAST files
        patterns = [
            os.path.join(base_path, "domain_blast_results", f"{pdb_id}_{chain_id}.domain_blast.xml"),
            os.path.join(base_path, "blast", "domain", f"{pdb_id}_{chain_id}.xml"),
            os.path.join(base_path, "blast", "domain", f"{pdb_id}_{chain_id}.blast.xml"),
            # Legacy patterns
            os.path.join(base_path, "blast", f"{pdb_id}_{chain_id}.domain_blast.xml")
        ]

    # Check directories for matching files
    for dir_pattern in [
        os.path.join(base_path, "chain_blast_results"),
        os.path.join(base_path, "domain_blast_results"),
        os.path.join(base_path, "blast", "chain"),
        os.path.join(base_path, "blast", "domain"),
        os.path.join(base_path, "blast")
    ]:
        if os.path.exists(dir_pattern):
            for filename in os.listdir(dir_pattern):
                if f"{pdb_id}_{chain_id}" in filename and filename.endswith(".xml"):
                    patterns.append(os.path.join(dir_pattern, filename))

    return patterns

def fix_db_fs_sync(context: ApplicationContext, batch_id: int, dry_run: bool = False) -> Dict:
    """
    Fix synchronization issues between database and filesystem.

    Args:
        context: Application context
        batch_id: Batch ID to fix
        dry_run: Whether to perform a dry run (no actual changes)

    Returns:
        Dictionary with results
    """
    logger = logging.getLogger("batch_blast.fix_sync")

    # Get base path for batch
    query = "SELECT base_path FROM ecod_schema.batch WHERE id = %s"
    result = context.db.execute_query(query, (batch_id,))
    if not result:
        logger.error(f"Batch {batch_id} not found or has no base_path")
        return {"error": "Batch not found"}

    base_path = result[0][0]

    # Find mismatched proteins
    mismatched_proteins = find_proteins_with_db_fs_mismatch(context, batch_id)

    if not mismatched_proteins:
        logger.info(f"No synchronization issues found for batch {batch_id}")
        return {
            "batch_id": batch_id,
            "mismatched_proteins": 0,
            "fixed_count": 0
        }

    logger.info(f"Found {len(mismatched_proteins)} proteins with synchronization issues")

    # Fix mismatches
    fixed_count = 0

    for protein in mismatched_proteins:
        pdb_id = protein["pdb_id"]
        chain_id = protein["chain_id"]

        logger.info(f"Fixing synchronization for {pdb_id}_{chain_id}")

        for file_info in protein["files"]:
            if file_info.get("has_mismatch", False):
                file_type = file_info["file_type"]
                file_id = file_info["file_id"]

                # Case 1: File exists in DB but not on FS
                if file_info["file_exists_db"] and not file_info.get("file_exists_fs", False):
                    logger.info(f"  {file_type} exists in DB but not on FS")

                    # Look for the file in alternative locations
                    potential_paths = find_potential_blast_paths(base_path, pdb_id, chain_id, file_type)
                    found_path = None

                    for path in potential_paths:
                        if os.path.exists(path):
                            found_path = path
                            break

                    if found_path:
                        # File found in alternative location, update DB record
                        logger.info(f"  Found file at: {found_path}")
                        rel_path = os.path.relpath(found_path, base_path)

                        if not dry_run:
                            context.db.update(
                                "ecod_schema.process_file",
                                {
                                    "file_path": rel_path,
                                    "file_exists": True,
                                    "file_size": os.path.getsize(found_path)
                                },
                                "id = %s",
                                (file_id,)
                            )
                            fixed_count += 1
                        else:
                            logger.info(f"  Would update file path to: {rel_path}")
                    else:
                        # File truly doesn't exist, update DB record
                        logger.info(f"  File not found in any location, marking as non-existent")

                        if not dry_run:
                            context.db.update(
                                "ecod_schema.process_file",
                                {
                                    "file_exists": False,
                                    "file_size": 0
                                },
                                "id = %s",
                                (file_id,)
                            )
                            fixed_count += 1
                        else:
                            logger.info(f"  Would mark file as non-existent")

                # Case 2: File exists on FS but not in DB
                elif not file_info["file_exists_db"] and file_info.get("file_exists_fs", False):
                    logger.info(f"  {file_type} exists on FS but not in DB")

                    if not dry_run:
                        context.db.update(
                            "ecod_schema.process_file",
                            {
                                "file_exists": True,
                                "file_size": file_info.get("actual_file_size", 0)
                            },
                            "id = %s",
                            (file_id,)
                        )
                        fixed_count += 1
                    else:
                        logger.info(f"  Would mark file as existing")

    logger.info(f"Fixed {fixed_count} synchronization issues" if not dry_run else f"Would fix {fixed_count} synchronization issues")

    return {
        "batch_id": batch_id,
        "mismatched_proteins": len(mismatched_proteins),
        "fixed_count": fixed_count
    }

def generate_sync_report(context: ApplicationContext, batch_id: int, mismatched_proteins: List[Dict], output_path: str = None) -> str:
    """
    Generate a report on synchronization issues.

    Args:
        context: Application context
        batch_id: Batch ID
        mismatched_proteins: List of proteins with mismatches
        output_path: Path to save the report (optional)

    Returns:
        Path to the report file
    """
    logger = logging.getLogger("batch_blast.report")

    # Get batch information
    query = "SELECT batch_name FROM ecod_schema.batch WHERE id = %s"
    result = context.db.execute_query(query, (batch_id,))
    batch_name = result[0][0] if result else f"Batch {batch_id}"

    # Generate report content
    report_content = f"BLAST File Synchronization Report for {batch_name}\n"
    report_content += f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"

    report_content += f"Found {len(mismatched_proteins)} proteins with database-filesystem mismatches\n\n"

    if mismatched_proteins:
        # Group by mismatch type for better reporting
        db_only = []
        fs_only = []
        both_issues = []

        for protein in mismatched_proteins:
            db_only_files = []
            fs_only_files = []

            for file_info in protein["files"]:
                if file_info.get("has_mismatch", False):
                    if file_info["file_exists_db"] and not file_info.get("file_exists_fs", False):
                        db_only_files.append(file_info["file_type"])
                    elif not file_info["file_exists_db"] and file_info.get("file_exists_fs", False):
                        fs_only_files.append(file_info["file_type"])

            if db_only_files and fs_only_files:
                both_issues.append({
                    "protein": protein,
                    "db_only": db_only_files,
                    "fs_only": fs_only_files
                })
            elif db_only_files:
                db_only.append({
                    "protein": protein,
                    "file_types": db_only_files
                })
            elif fs_only_files:
                fs_only.append({
                    "protein": protein,
                    "file_types": fs_only_files
                })

        # Report on each category
        report_content += f"1. Files exist in database but not on filesystem: {len(db_only)} proteins\n"
        if db_only:
            for i, item in enumerate(db_only[:10]):  # Limit to first 10
                protein = item["protein"]
                report_content += f"   {i+1}. {protein['pdb_id']}_{protein['chain_id']}: {', '.join(item['file_types'])}\n"
            if len(db_only) > 10:
                report_content += f"   ... and {len(db_only) - 10} more proteins\n"
        report_content += "\n"

        report_content += f"2. Files exist on filesystem but not in database: {len(fs_only)} proteins\n"
        if fs_only:
            for i, item in enumerate(fs_only[:10]):  # Limit to first 10
                protein = item["protein"]
                report_content += f"   {i+1}. {protein['pdb_id']}_{protein['chain_id']}: {', '.join(item['file_types'])}\n"
            if len(fs_only) > 10:
                report_content += f"   ... and {len(fs_only) - 10} more proteins\n"
        report_content += "\n"

        report_content += f"3. Proteins with both types of issues: {len(both_issues)}\n"
        if both_issues:
            for i, item in enumerate(both_issues[:10]):  # Limit to first 10
                protein = item["protein"]
                report_content += f"   {i+1}. {protein['pdb_id']}_{protein['chain_id']}:\n"
                report_content += f"      - In DB only: {', '.join(item['db_only'])}\n"
                report_content += f"      - On FS only: {', '.join(item['fs_only'])}\n"
            if len(both_issues) > 10:
                report_content += f"   ... and {len(both_issues) - 10} more proteins\n"
        report_content += "\n"

        # Add recommendations
        report_content += "Recommendations:\n"
        report_content += "1. For files that exist in DB but not on FS:\n"
        report_content += "   - Update database records to mark files as non-existent\n"
        report_content += "   - Regenerate missing BLAST files if needed\n\n"

        report_content += "2. For files that exist on FS but not in DB:\n"
        report_content += "   - Update database records to register these files\n\n"

        report_content += "3. To fix these issues automatically, run:\n"
        report_content += f"   python batch_blast.py --config <config_file> --batch-id {batch_id} --fix-sync\n"

    # Write to file if output path provided
    if output_path:
        try:
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            with open(output_path, 'w') as f:
                f.write(report_content)
            logger.info(f"Report saved to {output_path}")
        except Exception as e:
            logger.error(f"Error writing report to {output_path}: {e}")
            output_path = None

    # If no output path provided or writing failed, use default path
    if not output_path:
        output_path = f"batch_{batch_id}_sync_report_{datetime.now().strftime('%Y%m%d_%H%M')}.txt"
        try:
            with open(output_path, 'w') as f:
                f.write(report_content)
            logger.info(f"Report saved to {output_path}")
        except Exception as e:
            logger.error(f"Error writing report to {output_path}: {e}")
            return None

    return output_path

def get_pending_batches(context: ApplicationContext) -> List[Dict[str, Any]]:
    """Get batches that need blast processing"""
    query = """
    SELECT 
        b.id, b.batch_name, b.base_path, b.type, b.ref_version
    FROM 
        ecod_schema.batch b
    WHERE 
        b.status = 'created'
        OR EXISTS (
            SELECT 1 FROM ecod_schema.process_status ps
            WHERE ps.batch_id = b.id
            AND ps.status IN ('pending', 'processing')
        )
    ORDER BY 
        b.id
    """
    
    return context.db.execute_dict_query(query)

def get_batch_statistics(context: ApplicationContext, batch_id: int) -> Dict[str, int]:
    """Get statistics for a batch"""
    query = """
    SELECT
        COUNT(*) as total,
        SUM(CASE WHEN ps.status = 'pending' THEN 1 ELSE 0 END) as pending,
        SUM(CASE WHEN ps.status = 'processing' THEN 1 ELSE 0 END) as processing,
        SUM(CASE WHEN ps.status = 'success' THEN 1 ELSE 0 END) as completed,
        SUM(CASE WHEN ps.status = 'error' THEN 1 ELSE 0 END) as failed
    FROM
        ecod_schema.process_status ps
    WHERE
        ps.batch_id = %s
    """
    
    result = context.db.execute_dict_query(query, (batch_id,))
    return result[0] if result else {}

def run_chain_blast(pipeline: BlastPipeline, batch_id: int, threads: int, 
                  memory: str, batch_size: int, force: bool
) -> List[str]:
    """Run chain blast for a batch"""
    try:
        return pipeline.run_chain_blast(
            batch_id=batch_id,
            batch_size=batch_size
        )
    except (ConfigurationError, PipelineError) as e:
        logger.error(f"Error running chain blast: {str(e)}")
        return []

def run_domain_blast(pipeline: BlastPipeline, batch_id: int, threads: int, 
                   memory: str, batch_size: int, force: bool
) -> List[str]:
    """Run domain blast for a batch"""
    try:
        return pipeline.run_domain_blast(
            batch_id=batch_id,
            batch_size=batch_size
        )
    except (ConfigurationError, PipelineError) as e:
        logger.error(f"Error running domain blast: {str(e)}")
        return []

def process_batch(context: ApplicationContext, batch_id: int, threads: int, 
                memory: str, batch_size: int, force: bool
) -> None:
    """Process a single batch"""
    logger.info(f"Processing batch {batch_id}")
    
    # Get batch statistics before processing
    stats_before = get_batch_statistics(context, batch_id)
    logger.info(f"Batch {batch_id} statistics before processing: {stats_before}")
    
    # Create pipeline
    pipeline = BlastPipeline(context)
    
    # Run chain blast
    logger.info(f"Running chain blast for batch {batch_id}")
    chain_jobs = run_chain_blast(pipeline, batch_id, threads, memory, batch_size, force)
    logger.info(f"Submitted {len(chain_jobs)} chain blast jobs")
    
    # Run domain blast
    logger.info(f"Running domain blast for batch {batch_id}")
    domain_jobs = run_domain_blast(pipeline, batch_id, threads, memory, batch_size, force)
    logger.info(f"Submitted {len(domain_jobs)} domain blast jobs")
    
    # Update batch status if no jobs were submitted
    if not chain_jobs and not domain_jobs:
        # Get current batch status
        query = "SELECT status FROM ecod_schema.batch WHERE id = %s"
        result = context.db.execute_query(query, (batch_id,))
        current_status = result[0][0] if result else None
        
        # Update status to 'processed' if currently 'created'
        if current_status == 'created':
            context.db.update(
                "ecod_schema.batch",
                {"status": "processed"},
                "id = %s",
                (batch_id,)
            )
            logger.info(f"Updated batch {batch_id} status to 'processed'")
    
    # Wait for jobs to complete if requested
    if chain_jobs or domain_jobs:
        logger.info(f"Checking job status for batch {batch_id}")
        pipeline.check_job_status(batch_id)

def main():
    parser = argparse.ArgumentParser(description="Run BLAST for ECOD batches")
    parser.add_argument("--config", default="config/config.yml", help="Path to config file")
    parser.add_argument("--batch-id", type=int, help="Process specific batch ID")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--memory", default="4G", help="Memory allocation")
    parser.add_argument("--batch-size", type=int, default=100, help="Batch size for job creation")
    parser.add_argument("--force", action="store_true", help="Force reprocessing")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--log-file", help="Log file path")
    parser.add_argument("--check-sync", action="store_true", help="Check synchronization between database and filesystem")
    parser.add_argument("--fix-sync", action="store_true", help="Fix synchronization issues between database and filesystem")
    parser.add_argument("--output-report", type=str, help="Path to save the synchronization report")
    args = parser.parse_args()
    
    # Setup logging
    global logger
    logger = setup_logging(args.verbose, args.log_file)
    
    # Create application context
    context = ApplicationContext(args.config)
    
    # Handle sync checking mode
    if args.check_sync or args.fix_sync:
        if not args.batch_id:
            logger.error("Batch ID is required for sync checking/fixing")
            return 1

        logger.info(f"{'Checking' if args.check_sync else 'Fixing'} synchronization for batch {args.batch_id}")

        # Find synchronization issues
        mismatched_proteins = find_proteins_with_db_fs_mismatch(context, args.batch_id)

        if not mismatched_proteins:
            logger.info(f"No synchronization issues found for batch {args.batch_id}")
            return 0

        logger.info(f"Found {len(mismatched_proteins)} proteins with synchronization issues")

        # Generate report
        report_path = generate_sync_report(context, args.batch_id, mismatched_proteins, args.output_report)
        if report_path:
            logger.info(f"Synchronization report saved to {report_path}")

        # Fix issues if requested
        if args.fix_sync:
            result = fix_db_fs_sync(context, args.batch_id, args.dry_run)
            if "error" in result:
                logger.error(f"Error fixing synchronization issues: {result['error']}")
                return 1

            action = "Would fix" if args.dry_run else "Fixed"
            logger.info(f"{action} {result['fixed_count']} synchronization issues for batch {args.batch_id}")

        return 0

    # Process specific batch or all pending batches
    if args.batch_id:
        process_batch(
            context=context,
            batch_id=args.batch_id,
            threads=args.threads,
            memory=args.memory,
            batch_size=args.batch_size,
            force=args.force
        )
    else:
        # Get all pending batches
        batches = get_pending_batches(context)
        logger.info(f"Found {len(batches)} pending batches")
        
        for batch in batches:
            process_batch(
                context=context,
                batch_id=batch['id'],
                threads=args.threads,
                memory=args.memory,
                batch_size=args.batch_size,
                force=args.force
            )

if __name__ == "__main__":
    main()
