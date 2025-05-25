#!/usr/bin/env python3
"""
Domain Partition Tool - Service-Based Implementation

This script has been updated to use the new service-oriented architecture:
- DomainPartitionService: High-level service interface
- PartitionOptions: Configuration management
- StatusTracker: Database status tracking
- BatchPartitionResults: Comprehensive result handling

All direct DomainPartition class usage has been replaced with service calls.
"""

import os
import sys
import logging
import argparse
import glob
import time
import re
import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple, Union, Set

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

# Core imports
from ecod.config import ConfigManager
from ecod.core.context import ApplicationContext
from ecod.db import DBManager

# NEW SERVICE-BASED IMPORTS
from ecod.pipelines.domain_analysis.partition import (
    DomainPartitionService,
    create_service,
    PartitionOptions,
    ValidationLevel,
    ProcessingMode,
    BatchPartitionResults
)
from ecod.models.pipeline.partition import DomainPartitionResult

# Utility imports
from ecod.utils.path_utils import get_standardized_paths, get_all_evidence_paths, resolve_file_path


def setup_logging(verbose: bool = False, log_file: Optional[str] = None) -> logging.Logger:
    """Configure logging with appropriate handlers and format"""
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

    return logging.getLogger("domain_partition")


#
# PROCESS MODE FUNCTIONS - UPDATED FOR SERVICE ARCHITECTURE
#

def process_batch(args: argparse.Namespace) -> int:
    """Process domain partition for a batch using the new service"""

    logger = logging.getLogger("domain_partition.process.batch")
    logger.info(f"Processing batch {args.batch_id} with service-based architecture")

    # Initialize context with config
    context = ApplicationContext(args.config)

    # Get batch information
    batch_info = get_batch_info(context, args.batch_id)
    if not batch_info:
        logger.error(f"Batch {args.batch_id} not found")
        return 1

    # Verify batch readiness if not forcing
    if not args.force and not verify_batch_readiness(context, args.batch_id, args.blast_only, args.reps_only):
        logger.error(f"Batch {args.batch_id} is not ready for domain partition")
        logger.error("Use --force to override or run domain summary first")
        return 1

    # Create service with custom configuration
    service_config = {
        'max_workers': getattr(args, 'workers', 4),
        'use_multiprocessing': getattr(args, 'use_multiprocessing', False),
        'save_intermediate': True,
        'track_status': True
    }

    service = DomainPartitionService(context, service_config)

    # Process the batch
    batch_path = batch_info["base_path"]

    try:
        logger.info(f"Running service-based partition for batch {args.batch_id}")

        # Create partition options
        options = {
            'blast_only': args.blast_only,
            'representatives_only': args.reps_only,
            'force_overwrite': args.force,
            'validation_level': ValidationLevel.NORMAL
        }

        # Process batch
        results = service.partition_batch(
            batch_id=args.batch_id,
            batch_path=batch_path,
            limit=args.limit,
            **options
        )

        # Log results
        logger.info(f"Domain partition complete: {results.success_count} succeeded, {results.failure_count} failed")

        # Show summary statistics
        if results.success_count > 0:
            logger.info("Classification summary:")
            logger.info(f"  Classified: {results.proteins_with_domains}")
            logger.info(f"  Peptides: {results.peptides_found}")
            logger.info(f"  Unclassified: {results.unclassified_proteins}")
            logger.info(f"  Total domains: {results.total_domains_found}")

            if results.proteins_with_domains > 0:
                avg_domains = results.total_domains_found / results.proteins_with_domains
                logger.info(f"  Average domains per classified protein: {avg_domains:.1f}")

        # Log failures if any
        if results.failure_count > 0:
            for i, (pdb_id, chain_id, error) in enumerate(results.failures[:3]):
                logger.error(f"Failure {i+1}: {pdb_id}_{chain_id}: {error}")

            if results.failure_count > 3:
                logger.error(f"... and {results.failure_count - 3} more failures")

        # Get service statistics
        stats = service.get_service_statistics()
        logger.info(f"Processing statistics: {stats['service']['proteins_processed']} proteins in "
                   f"{stats['service']['runtime_seconds']:.1f}s "
                   f"({stats['service']['proteins_per_minute']:.1f} proteins/min)")

        return 0 if results.success_count > 0 else 1

    except Exception as e:
        logger.error(f"Error processing batch: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return 1


def process_specific_proteins(args: argparse.Namespace) -> int:
    """Process domain partition for specific proteins using the service"""

    logger = logging.getLogger("domain_partition.process.specific")

    # Initialize context with config
    context = ApplicationContext(args.config)

    # Get batch ID if not specified
    batch_id = args.batch_id
    if not batch_id and args.process_ids:
        batch_id = get_batch_id_for_process(context, args.process_ids[0])

        if not batch_id:
            logger.error(f"Could not determine batch ID for process {args.process_ids[0]}")
            return 1

        logger.info(f"Using batch ID {batch_id} for process {args.process_ids[0]}")

    # Get batch information
    batch_info = get_batch_info(context, batch_id)
    if not batch_info:
        logger.error(f"Batch {batch_id} not found")
        return 1

    # Create service
    service = DomainPartitionService(context)

    # Get protein details for each process ID
    proteins = []
    for process_id in args.process_ids:
        protein_info = get_protein_for_process(context, process_id)
        if protein_info:
            proteins.append(protein_info)
        else:
            logger.warning(f"Process ID {process_id} not found")

    if not proteins:
        logger.error("No valid proteins found for specified process IDs")
        return 1

    # Process each protein
    success_count = 0
    failure_count = 0

    try:
        logger.info(f"Processing {len(proteins)} specific proteins from batch {batch_id}")

        for protein in proteins:
            pdb_id = protein['pdb_id']
            chain_id = protein['chain_id']
            process_id = protein['process_id']

            # Find domain summary
            summary_path = find_domain_summary(
                service, batch_info["base_path"], pdb_id, chain_id,
                batch_info["ref_version"], args.blast_only
            )

            if not summary_path:
                logger.error(f"Domain summary not found for {pdb_id}_{chain_id}")
                failure_count += 1
                continue

            # Process protein
            result = service.partition_protein(
                pdb_id=pdb_id,
                chain_id=chain_id,
                summary_path=summary_path,
                output_dir=batch_info["base_path"],
                process_id=process_id,
                blast_only=args.blast_only
            )

            if result.success:
                success_count += 1
                logger.info(f"Successfully processed {pdb_id}_{chain_id}")
            else:
                failure_count += 1
                logger.error(f"Failed to process {pdb_id}_{chain_id}: {result.error}")

        logger.info(f"Domain partition complete: {success_count} succeeded, {failure_count} failed")
        return 0 if success_count > 0 else 1

    except Exception as e:
        logger.error(f"Error processing specific proteins: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return 1


def process_single_protein(args: argparse.Namespace) -> int:
    """Process domain partition for a single protein using the service"""

    logger = logging.getLogger("domain_partition.process.single")

    if not args.pdb_id or not args.chain_id:
        logger.error("PDB ID and chain ID are required for single protein processing")
        return 1

    # Initialize context with config
    context = ApplicationContext(args.config)

    # Get batch information if batch ID provided
    batch_path = args.batch_path
    reference = args.reference
    process_id = None

    if args.batch_id:
        batch_info = get_batch_info(context, args.batch_id)
        if not batch_info:
            logger.error(f"Batch {args.batch_id} not found")
            return 1

        batch_path = batch_info["base_path"]
        reference = batch_info["ref_version"]

        # Look up process ID for status tracking
        process_id = get_process_id_for_protein(context, args.pdb_id, args.chain_id, args.batch_id)

    if not batch_path or not reference:
        logger.error("Batch path and reference are required - specify directly or provide batch ID")
        return 1

    # Create service
    service = DomainPartitionService(context)

    # Check if domain summary exists
    if not verify_domain_summary(batch_path, args.pdb_id, args.chain_id, reference, args.blast_only):
        logger.error(f"Domain summary not found for {args.pdb_id}_{args.chain_id}")
        logger.error("Run domain summary first or provide correct paths")
        return 1

    # Process single protein
    try:
        logger.info(f"Processing {args.pdb_id}_{args.chain_id} with service-based architecture")

        # Find domain summary path
        summary_path = find_domain_summary(
            service, batch_path, args.pdb_id, args.chain_id, reference, args.blast_only
        )

        if not summary_path:
            logger.error(f"Domain summary file not found for {args.pdb_id}_{args.chain_id}")
            return 1

        # Process using service
        result = service.partition_protein(
            pdb_id=args.pdb_id,
            chain_id=args.chain_id,
            summary_path=summary_path,
            output_dir=batch_path,
            process_id=process_id,
            blast_only=args.blast_only
        )

        # Report outcome
        if result.success:
            logger.info(f"Successfully processed {args.pdb_id}_{args.chain_id}")

            if result.domain_file:
                logger.info(f"Created domain file: {result.domain_file}")

            if result.is_peptide:
                logger.info("Classified as peptide")
            elif result.is_classified:
                logger.info(f"Classified with {len(result.domains)} domains")

                # Show domain details
                for i, domain in enumerate(result.domains):
                    logger.info(f"  Domain {i+1}: {domain.range} ({domain.get_classification_level()})")
            else:
                logger.info("Unclassified")

            # Show coverage statistics
            if result.sequence_length > 0:
                logger.info(f"Coverage: {result.coverage:.1%} ({result.residues_assigned}/{result.sequence_length} residues)")

        else:
            logger.error(f"Failed to process {args.pdb_id}_{args.chain_id}: {result.error}")

        return 0 if result.success else 1

    except Exception as e:
        logger.error(f"Error processing {args.pdb_id}_{args.chain_id}: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return 1


def process_all_batches(args: argparse.Namespace) -> int:
    """Process domain partition for all batches using the service"""

    logger = logging.getLogger("domain_partition.process.all")

    # Initialize context with config
    context = ApplicationContext(args.config)

    # Determine which batches to process
    batch_ids = args.batch_ids
    if not batch_ids:
        # Get all batches if not specified
        all_batches = get_all_batch_ids(context)
        logger.info(f"Found {len(all_batches)} total batches")

        # Filter by reference version if specified
        if args.reference:
            reference_batches = []
            for batch_id in all_batches:
                batch_info = get_batch_info(context, batch_id)
                if batch_info and batch_info.get("ref_version") == args.reference:
                    reference_batches.append(batch_id)

            batch_ids = reference_batches
            logger.info(f"Filtered to {len(batch_ids)} batches with reference {args.reference}")
        else:
            batch_ids = all_batches

    # Exclude specified batch IDs
    if args.exclude_batch_ids:
        batch_ids = [b_id for b_id in batch_ids if b_id not in args.exclude_batch_ids]
        logger.info(f"Excluding specified batches, {len(batch_ids)} batches remaining")

    # Sort batch IDs for consistent processing
    batch_ids.sort()

    if not batch_ids:
        logger.warning("No batches to process")
        return 0

    logger.info(f"Processing {len(batch_ids)} batches with service-based architecture: {batch_ids}")

    # Determine whether to use SLURM or process directly
    if args.use_slurm:
        return submit_batches_to_slurm(context, batch_ids, args)
    else:
        return process_batches_directly(context, batch_ids, args)


def submit_batches_to_slurm(context: ApplicationContext, batch_ids: List[int], args: argparse.Namespace) -> int:
    """Submit batch processing jobs to SLURM"""

    logger = logging.getLogger("domain_partition.process.slurm")
    logger.info(f"Submitting {len(batch_ids)} batches to SLURM")

    # Get job manager
    try:
        from ecod.jobs import SlurmJobManager
        job_manager = SlurmJobManager(context.config_manager.config)
    except ImportError:
        logger.error("Failed to import SlurmJobManager. Make sure SLURM integration is available.")
        return 1

    # Create temporary directory for job scripts
    temp_dir = os.path.join(context.config_manager.get_path('output_dir', '/tmp'),
                           f"domain_partition_jobs_{int(time.time())}")
    os.makedirs(temp_dir, exist_ok=True)
    logger.info(f"Created job directory: {temp_dir}")

    # Get paths
    config_path = os.path.abspath(args.config)
    script_path = os.path.abspath(sys.argv[0])

    # Process batches with appropriate batch size
    batch_groups = []
    for i in range(0, len(batch_ids), args.batch_size):
        batch_groups.append(batch_ids[i:i+args.batch_size])

    logger.info(f"Split into {len(batch_groups)} groups with max {args.batch_size} batches per group")

    # Submit jobs
    job_ids = []

    for group_idx, group in enumerate(batch_groups):
        group_dir = os.path.join(temp_dir, f"group_{group_idx}")
        os.makedirs(group_dir, exist_ok=True)

        for batch_id in group:
            # Create job name
            job_name = f"domain_partition_batch_{batch_id}"

            # Build command
            command = f"python {script_path} --config {config_path}"

            # Add global options
            if args.verbose:
                command += " --verbose"

            # Add subcommand and options
            command += f" process batch --batch-id {batch_id}"

            if args.blast_only:
                command += " --blast-only"
            if args.limit_per_batch:
                command += f" --limit {args.limit_per_batch}"
            if args.reps_only:
                command += " --reps-only"
            if args.force:
                command += " --force"

            # Create and submit job
            script_path_gen = job_manager.create_job_script(
                commands=[command],
                job_name=job_name,
                output_dir=group_dir,
                threads=args.slurm_threads,
                memory=args.slurm_memory,
                time=args.slurm_time
            )

            job_id = job_manager.submit_job(script_path_gen)

            if job_id:
                job_ids.append(job_id)
                logger.info(f"Submitted job for batch {batch_id}, SLURM job ID: {job_id}")
            else:
                logger.error(f"Failed to submit job for batch {batch_id}")

        # Wait between batch groups if specified
        if args.wait_between_groups and group_idx < len(batch_groups) - 1:
            logger.info(f"Waiting {args.wait_between_groups} seconds before next batch group")
            time.sleep(args.wait_between_groups)

    logger.info(f"Submitted {len(job_ids)} jobs to SLURM")

    # Monitor jobs if requested
    if args.wait_for_completion:
        return monitor_slurm_jobs(job_manager, job_ids, args, logger)

    return 0


def process_batches_directly(context: ApplicationContext, batch_ids: List[int], args: argparse.Namespace) -> int:
    """Process batches directly using the service"""

    logger = logging.getLogger("domain_partition.process.direct")

    # Create service with configuration for batch processing
    service_config = {
        'max_workers': getattr(args, 'workers', 1),
        'use_multiprocessing': False,  # Sequential for direct processing
        'batch_size': args.batch_size
    }

    service = DomainPartitionService(context, service_config)

    # Process batches with appropriate batch size
    batch_groups = []
    for i in range(0, len(batch_ids), args.batch_size):
        batch_groups.append(batch_ids[i:i+args.batch_size])

    logger.info(f"Split into {len(batch_groups)} groups with max {args.batch_size} batches per group")

    # Track overall statistics
    overall_results = BatchPartitionResults()
    success_count = 0
    failed_batches = []

    for group_idx, group in enumerate(batch_groups):
        logger.info(f"Processing batch group {group_idx+1}/{len(batch_groups)}: {group}")

        group_success = 0

        for batch_id in group:
            batch_info = get_batch_info(context, batch_id)
            if not batch_info:
                logger.error(f"Batch {batch_id} not found")
                failed_batches.append(batch_id)
                continue

            logger.info(f"Processing batch {batch_id} ({batch_info.get('batch_name', '')})")

            # Check batch readiness if not forcing
            if not args.force and not verify_batch_readiness(context, batch_id, args.blast_only, args.reps_only):
                logger.warning(f"Batch {batch_id} is not ready for domain partition, skipping")
                failed_batches.append(batch_id)
                continue

            # Process the batch
            try:
                batch_path = batch_info["base_path"]

                # Create options
                options = {
                    'blast_only': args.blast_only,
                    'representatives_only': args.reps_only,
                    'force_overwrite': args.force
                }

                results = service.partition_batch(
                    batch_id=batch_id,
                    batch_path=batch_path,
                    limit=args.limit_per_batch,
                    **options
                )

                # Merge results
                overall_results.total += results.total
                overall_results.success_count += results.success_count
                overall_results.failure_count += results.failure_count
                overall_results.proteins_with_domains += results.proteins_with_domains
                overall_results.total_domains_found += results.total_domains_found
                overall_results.peptides_found += results.peptides_found
                overall_results.unclassified_proteins += results.unclassified_proteins

                logger.info(f"Batch {batch_id} complete: {results.success_count} succeeded, {results.failure_count} failed")

                if results.success_count > 0:
                    group_success += 1
                    success_count += 1
                else:
                    failed_batches.append(batch_id)

            except Exception as e:
                logger.error(f"Error processing batch {batch_id}: {str(e)}")
                import traceback
                logger.error(traceback.format_exc())
                failed_batches.append(batch_id)

        logger.info(f"Batch group {group_idx+1} complete: {group_success}/{len(group)} succeeded")

        # Wait between batch groups if specified
        if args.wait_between_groups and group_idx < len(batch_groups) - 1:
            logger.info(f"Waiting {args.wait_between_groups} seconds before next batch group")
            time.sleep(args.wait_between_groups)

    # Log final summary
    overall_results.finalize()
    summary = overall_results.get_summary()

    logger.info(f"All batches processed: {success_count}/{len(batch_ids)} succeeded")
    logger.info("Overall statistics:")
    logger.info(f"  Total proteins processed: {summary['total_proteins']}")
    logger.info(f"  Successful: {summary['successful']}")
    logger.info(f"  Failed: {summary['failed']}")
    logger.info(f"  Classified: {summary['proteins_with_domains']}")
    logger.info(f"  Peptides: {summary['peptides']}")
    logger.info(f"  Unclassified: {summary['unclassified']}")
    logger.info(f"  Total domains: {summary['total_domains']}")
    logger.info(f"  Processing time: {summary['processing_time']:.1f}s")

    if summary['proteins_with_domains'] > 0:
        avg_domains = summary['total_domains'] / summary['proteins_with_domains']
        logger.info(f"  Average domains per classified protein: {avg_domains:.1f}")

    if failed_batches:
        logger.warning(f"Failed batches: {failed_batches}")

    # Clear caches after processing
    service.clear_all_caches()

    return 0 if len(failed_batches) == 0 else 1


def monitor_slurm_jobs(job_manager, job_ids: List[str], args: argparse.Namespace, logger) -> int:
    """Monitor SLURM jobs until completion"""

    logger.info(f"Waiting for jobs to complete, checking every {args.check_interval} seconds")

    start_time = time.time()
    completed_jobs = set()
    failed_jobs = set()

    while job_ids and (not args.timeout or time.time() - start_time < args.timeout):
        time.sleep(args.check_interval)

        # Check each job
        for job_id in list(job_ids):
            status = job_manager.check_job_status(job_id)

            if status in ["COMPLETED", "COMPLETING"]:
                logger.info(f"Job {job_id} completed successfully")
                job_ids.remove(job_id)
                completed_jobs.add(job_id)
            elif status in ["FAILED", "TIMEOUT", "CANCELLED", "NODE_FAIL"]:
                logger.error(f"Job {job_id} failed with status {status}")
                job_ids.remove(job_id)
                failed_jobs.add(job_id)

        # Log status
        if job_ids:
            logger.info(f"Waiting for {len(job_ids)} remaining jobs")

    # Log final status
    if not job_ids:
        logger.info("All jobs completed")
    else:
        logger.warning(f"{len(job_ids)} jobs still running at end of monitoring period")

    logger.info(f"Job completion summary: {len(completed_jobs)} completed, {len(failed_jobs)} failed")

    return 0 if not failed_jobs and not job_ids else 1


#
# ANALYZE MODE FUNCTIONS - UPDATED FOR SERVICE ARCHITECTURE
#

def analyze_batch_status(args: argparse.Namespace) -> int:
    """Analyze batch status for domain partition using service"""

    logger = logging.getLogger("domain_partition.analyze.batch")

    # Initialize context with config
    context = ApplicationContext(args.config)

    # Create service for status tracking
    service = DomainPartitionService(context)

    # Determine which batches to analyze
    batch_ids = args.batch_ids
    if not batch_ids:
        batch_ids = get_all_batch_ids(context)
        logger.info(f"Analyzing status for all {len(batch_ids)} batches")
    else:
        logger.info(f"Analyzing status for {len(batch_ids)} specified batches")

    # Create results table
    results = []

    for batch_id in batch_ids:
        batch_info = get_batch_info(context, batch_id)
        if not batch_info:
            logger.warning(f"Batch {batch_id} not found")
            continue

        # Get batch progress using service's tracker
        progress = service.tracker.get_batch_progress(batch_id)

        # Calculate readiness and completion
        total = progress.get("total", 0)
        ready = progress.get("files_created", 0)  # Files ready for partition
        complete = progress.get("complete", 0)

        readiness = (ready / total * 100) if total > 0 else 0
        completion = progress.get("complete_pct", 0)

        results.append({
            "batch_id": batch_id,
            "batch_name": batch_info.get("batch_name", ""),
            "total": total,
            "ready": ready,
            "complete": complete,
            "readiness": readiness,
            "completion": completion,
            "status": batch_info.get("status", ""),
            "errors": progress.get("errors", 0)
        })

    # Sort by batch ID
    results.sort(key=lambda x: x["batch_id"])

    # Display results
    logger.info(f"{'ID':5s} {'Name':20s} {'Total':8s} {'Ready':8s} {'Complete':8s} {'Ready%':8s} {'Complete%':8s} {'Errors':6s} {'Status':12s}")
    logger.info("-" * 100)

    for r in results:
        logger.info(f"{r['batch_id']:5d} {r['batch_name'][:20]:20s} {r['total']:8d} " +
                    f"{r['ready']:8d} {r['complete']:8d} {r['readiness']:7.1f}% " +
                    f"{r['completion']:8.1f}% {r['errors']:6d} {r['status'][:12]:12s}")

    logger.info("-" * 100)

    # Calculate totals
    total_proteins = sum(r["total"] for r in results)
    total_ready = sum(r["ready"] for r in results)
    total_complete = sum(r["complete"] for r in results)
    total_errors = sum(r["errors"] for r in results)

    logger.info(f"Total: {len(results)} batches, {total_proteins} proteins, " +
                f"{total_ready} ready ({total_ready/total_proteins*100:.1f}%), " +
                f"{total_complete} complete ({total_complete/total_proteins*100:.1f}%), " +
                f"{total_errors} errors")

    return 0


def analyze_protein_status(args: argparse.Namespace) -> int:
    """Analyze protein status using service"""

    logger = logging.getLogger("domain_partition.analyze.protein")

    if not args.pdb_id or not args.chain_id:
        logger.error("PDB ID and chain ID are required for protein status analysis")
        return 1

    # Initialize context with config
    context = ApplicationContext(args.config)

    # Find all occurrences of this protein
    proteins = find_protein_in_database(context, args.pdb_id, args.chain_id, args.batch_id)

    if not proteins:
        logger.error(f"Protein {args.pdb_id}_{args.chain_id} not found")
        return 1

    # Display protein information
    logger.info(f"Found {len(proteins)} instances of {args.pdb_id}_{args.chain_id}")

    for i, p in enumerate(proteins):
        batch_info = get_batch_info(context, p["batch_id"])
        batch_name = batch_info.get("batch_name", "unknown") if batch_info else "unknown"

        logger.info(f"Instance {i+1}:")
        logger.info(f"  Batch: {p['batch_id']} ({batch_name})")
        logger.info(f"  Process ID: {p['process_id']}")
        logger.info(f"  Status: {p['current_stage']} / {p['status']}")

        if p.get("is_representative"):
            logger.info(f"  Is Representative: Yes")

        if p["error_message"]:
            logger.info(f"  Error: {p['error_message']}")

        # Show file information
        files = p.get("files", {})
        logger.info("  Files:")
        for file_type, file_info in files.items():
            exists = "EXISTS" if file_info.get("file_exists", False) else "MISSING"
            logger.info(f"    {file_type:20s}: {exists}")

            if args.verbose and file_info.get("file_exists", False):
                file_path = file_info.get("file_path", "")
                if file_path and batch_info:
                    full_path = os.path.join(batch_info.get("base_path", ""), file_path)
                    logger.info(f"      Path: {full_path}")

        # Show domain information if available
        if batch_info and p.get("status") == "success" and "domain_partition" in files:
            try:
                domain_info = check_domain_result_with_models(
                    context, p["process_id"], args.pdb_id, args.chain_id, batch_info
                )

                if domain_info:
                    logger.info("  Domain Information:")
                    logger.info(f"    Is Classified: {domain_info.get('is_classified', False)}")
                    logger.info(f"    Is Peptide: {domain_info.get('is_peptide', False)}")
                    logger.info(f"    Domain Count: {len(domain_info.get('domains', []))}")
                    logger.info(f"    Coverage: {domain_info.get('coverage', 0):.1%}")
                    logger.info(f"    Processing Time: {domain_info.get('processing_time', 0):.2f}s")

                    # Show domain details in verbose mode
                    if args.verbose and domain_info.get("domains"):
                        logger.info("    Domains:")
                        for j, domain in enumerate(domain_info["domains"]):
                            classification = domain.get_classification_level() if hasattr(domain, 'get_classification_level') else "Unknown"
                            range_str = domain.range if hasattr(domain, 'range') else domain.get('range', 'unknown')
                            confidence = getattr(domain, 'confidence', domain.get('confidence', 0)) if domain else 0

                            logger.info(f"      Domain {j+1}: {range_str} ({classification}, conf={confidence:.3f})")

            except Exception as e:
                logger.error(f"Error reading domain information: {str(e)}")

    return 0


def check_domain_result_with_models(context: ApplicationContext, process_id: int, pdb_id: str, chain_id: str,
                                   batch_info: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """Check domain result file using model-based parsing"""

    # Get file path from database
    query = """
    SELECT file_path FROM ecod_schema.process_file
    WHERE process_id = %s AND file_type = 'domain_partition'
    """

    rows = context.db.execute_query(query, (process_id,))
    if not rows:
        return None

    file_path = rows[0][0]
    full_path = os.path.join(batch_info.get("base_path", ""), file_path)

    if not os.path.exists(full_path):
        return None

    # Parse domain file using model
    try:
        result = DomainPartitionResult.from_xml_file(full_path)

        return {
            "is_classified": result.is_classified,
            "is_unclassified": result.is_unclassified,
            "is_peptide": result.is_peptide,
            "coverage": result.coverage,
            "domains": result.domains,
            "sequence_length": result.sequence_length,
            "processing_time": result.processing_time
        }

    except Exception as e:
        logging.getLogger("domain_partition.analyze").error(f"Error parsing domain file: {str(e)}")
        return None


def analyze_domain_counts(args: argparse.Namespace) -> int:
    """Analyze domain count statistics using service"""

    logger = logging.getLogger("domain_partition.analyze.counts")

    # Initialize context with config
    context = ApplicationContext(args.config)

    # Determine which batches to analyze
    batch_ids = args.batch_ids
    if not batch_ids:
        batch_ids = get_all_batch_ids(context)
        logger.info(f"Analyzing domains for all {len(batch_ids)} batches")
    else:
        logger.info(f"Analyzing domains for {len(batch_ids)} specified batches")

    # Query database for domain files
    try:
        # Process each batch
        results = []

        for batch_id in batch_ids:
            batch_info = get_batch_info(context, batch_id)
            if not batch_info:
                continue

            # Get sample files for analysis
            sample_query = """
            SELECT p.pdb_id, p.chain_id, pf.file_path
            FROM ecod_schema.process_status ps
            JOIN ecod_schema.protein p ON ps.protein_id = p.id
            JOIN ecod_schema.process_file pf ON ps.id = pf.process_id
            WHERE pf.file_type = 'domain_partition'
              AND ps.batch_id = %s
              AND ps.status = 'success'
              AND pf.file_exists = TRUE
            ORDER BY RANDOM()
            LIMIT %s
            """

            samples = context.db.execute_dict_query(sample_query, (batch_id, args.sample_size))

            # Analyze using models
            domain_stats = {
                "total_domains": 0,
                "proteins_with_domains": 0,
                "multi_domain_proteins": 0,
                "single_domain_proteins": 0,
                "unclassified_proteins": 0,
                "peptide_proteins": 0,
                "domain_counts": {},
                "confidence_stats": [],
                "evidence_stats": []
            }

            for sample in samples:
                full_path = os.path.join(batch_info["base_path"], sample["file_path"])

                if os.path.exists(full_path):
                    try:
                        # Parse using DomainPartitionResult
                        result = DomainPartitionResult.from_xml_file(full_path)

                        if result.is_peptide:
                            domain_stats["peptide_proteins"] += 1
                        elif result.is_unclassified:
                            domain_stats["unclassified_proteins"] += 1
                        elif result.is_classified:
                            domain_count = len(result.domains)

                            if domain_count > 0:
                                domain_stats["proteins_with_domains"] += 1
                                domain_stats["total_domains"] += domain_count

                                # Track distribution
                                if domain_count not in domain_stats["domain_counts"]:
                                    domain_stats["domain_counts"][domain_count] = 0
                                domain_stats["domain_counts"][domain_count] += 1

                                # Track multi-domain proteins
                                if domain_count > 1:
                                    domain_stats["multi_domain_proteins"] += 1
                                else:
                                    domain_stats["single_domain_proteins"] += 1

                                # Collect confidence and evidence statistics
                                for domain in result.domains:
                                    if hasattr(domain, 'confidence'):
                                        domain_stats["confidence_stats"].append(domain.confidence)
                                    if hasattr(domain, 'evidence'):
                                        domain_stats["evidence_stats"].append(len(domain.evidence))

                    except Exception as e:
                        logger.debug(f"Error parsing domain file {full_path}: {str(e)}")

            # Calculate averages
            avg_domains = (domain_stats["total_domains"] / domain_stats["proteins_with_domains"]
                          if domain_stats["proteins_with_domains"] > 0 else 0)

            avg_confidence = (sum(domain_stats["confidence_stats"]) / len(domain_stats["confidence_stats"])
                            if domain_stats["confidence_stats"] else 0)

            avg_evidence = (sum(domain_stats["evidence_stats"]) / len(domain_stats["evidence_stats"])
                          if domain_stats["evidence_stats"] else 0)

            # Add to results
            results.append({
                "batch_id": batch_id,
                "batch_name": batch_info.get("batch_name", ""),
                "sample_size": len(samples),
                "proteins_with_domains": domain_stats["proteins_with_domains"],
                "multi_domain_proteins": domain_stats["multi_domain_proteins"],
                "single_domain_proteins": domain_stats["single_domain_proteins"],
                "unclassified_proteins": domain_stats["unclassified_proteins"],
                "peptide_proteins": domain_stats["peptide_proteins"],
                "avg_domains": avg_domains,
                "avg_confidence": avg_confidence,
                "avg_evidence": avg_evidence,
                "domain_counts": domain_stats["domain_counts"]
            })

        # Display results
        logger.info(f"{'ID':5s} {'Name':20s} {'Sample':6s} {'w/Domains':9s} {'Multi':5s} {'Peptides':8s} {'Avg Dom':7s} {'Avg Conf':8s} {'Avg Evid':8s}")
        logger.info("-" * 95)

        for r in results:
            logger.info(f"{r['batch_id']:5d} {r['batch_name'][:20]:20s} {r['sample_size']:6d} " +
                        f"{r['proteins_with_domains']:9d} " +
                        f"{r['multi_domain_proteins']:5d} {r['peptide_proteins']:8d} " +
                        f"{r['avg_domains']:7.1f} {r['avg_confidence']:8.3f} {r['avg_evidence']:8.1f}")

        logger.info("-" * 95)

        # Combined statistics
        if args.verbose:
            combined_counts = {}
            for r in results:
                for count, num in r["domain_counts"].items():
                    if count not in combined_counts:
                        combined_counts[count] = 0
                    combined_counts[count] += num

            logger.info("Domain count distribution:")
            for count in sorted(combined_counts.keys()):
                logger.info(f"  {count} domains: {combined_counts[count]} proteins")

        return 0

    except Exception as e:
        logger.error(f"Error analyzing domain counts: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return 1


#
# REPAIR MODE FUNCTIONS - UPDATED FOR SERVICE ARCHITECTURE
#

def repair_failed_processes(args: argparse.Namespace) -> int:
    """Reset and retry failed processes using service"""

    logger = logging.getLogger("domain_partition.repair.failed")

    # Initialize context with config
    context = ApplicationContext(args.config)

    # Validate batch ID
    batch_info = get_batch_info(context, args.batch_id)
    if not batch_info:
        logger.error(f"Batch {args.batch_id} not found")
        return 1

    # Reset failed processes
    from ecod.pipelines.domain_analysis.pipeline import DomainAnalysisPipeline
    pipeline = DomainAnalysisPipeline(context)

    summary_reset = pipeline.reset_failed_processes(args.batch_id, 'domain_summary_failed')
    partition_reset = pipeline.reset_failed_processes(args.batch_id, 'domain_partition_failed')

    total_reset = summary_reset + partition_reset

    logger.info(f"Reset {total_reset} failed processes:")
    logger.info(f"  Summary failures reset: {summary_reset}")
    logger.info(f"  Partition failures reset: {partition_reset}")

    if total_reset == 0:
        logger.info("No failed processes to reset")
        return 0

    # Re-run if requested
    if args.rerun:
        logger.info("Re-running partition for reset processes using service")

        # Create service
        service = DomainPartitionService(context)

        # Reprocess the batch
        try:
            results = service.reprocess_failed(
                batch_id=args.batch_id,
                batch_path=batch_info["base_path"],
                blast_only=args.blast_only
            )

            logger.info(f"Re-run complete: {results.success_count} succeeded, {results.failure_count} failed")
            return 0 if results.success_count > 0 else 1

        except Exception as e:
            logger.error(f"Error re-running processes: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
            return 1

    return 0


def repair_missing_files(args: argparse.Namespace) -> int:
    """Regenerate missing domain files using service"""

    logger = logging.getLogger("domain_partition.repair.missing")

    # Initialize context with config
    context = ApplicationContext(args.config)

    # Validate batch ID
    batch_info = get_batch_info(context, args.batch_id)
    if not batch_info:
        logger.error(f"Batch {args.batch_id} not found")
        return 1

    # Find processes with missing domain files
    query = """
    SELECT ps.id as process_id, p.pdb_id, p.chain_id
    FROM ecod_schema.process_status ps
    JOIN ecod_schema.protein p ON ps.protein_id = p.id
    LEFT JOIN ecod_schema.process_file pf_summary ON
        ps.id = pf_summary.process_id AND
        pf_summary.file_type = %s AND
        pf_summary.file_exists = TRUE
    LEFT JOIN ecod_schema.process_file pf_partition ON
        ps.id = pf_partition.process_id AND
        pf_partition.file_type = 'domain_partition' AND
        pf_partition.file_exists = TRUE
    WHERE ps.batch_id = %s
      AND pf_summary.id IS NOT NULL
      AND pf_partition.id IS NULL
    """

    summary_type = "blast_only_summary" if args.blast_only else "domain_summary"
    rows = context.db.execute_dict_query(query, (summary_type, args.batch_id))

    if not rows:
        logger.info(f"No missing domain files found for batch {args.batch_id}")
        return 0

    logger.info(f"Found {len(rows)} proteins with missing domain files")

    # Limit number of files to repair
    if args.limit and args.limit < len(rows):
        rows = rows[:args.limit]
        logger.info(f"Limiting repair to {args.limit} proteins")

    # Create service
    service = DomainPartitionService(context)

    # Process each missing file
    success_count = 0
    for row in rows:
        try:
            # Find summary path
            summary_path = find_domain_summary(
                service, batch_info["base_path"], row["pdb_id"], row["chain_id"],
                batch_info["ref_version"], args.blast_only
            )

            if not summary_path:
                logger.warning(f"No summary found for {row['pdb_id']}_{row['chain_id']}")
                continue

            # Process protein
            result = service.partition_protein(
                pdb_id=row["pdb_id"],
                chain_id=row["chain_id"],
                summary_path=summary_path,
                output_dir=batch_info["base_path"],
                process_id=row["process_id"],
                blast_only=args.blast_only
            )

            if result.success:
                success_count += 1
                logger.info(f"Regenerated domain file for {row['pdb_id']}_{row['chain_id']}")
            else:
                logger.error(f"Failed to regenerate for {row['pdb_id']}_{row['chain_id']}: {result.error}")

        except Exception as e:
            logger.error(f"Error processing {row['pdb_id']}_{row['chain_id']}: {e}")

    logger.info(f"Repair complete: {success_count}/{len(rows)} files regenerated")
    return 0 if success_count > 0 else 1


#
# UTILITY FUNCTIONS
#

def get_batch_info(context: ApplicationContext, batch_id: int) -> Optional[Dict[str, Any]]:
    """Get batch information from database"""
    query = """
    SELECT id, batch_name, base_path, ref_version, status
    FROM ecod_schema.batch
    WHERE id = %s
    """

    rows = context.db.execute_dict_query(query, (batch_id,))
    return rows[0] if rows else None


def get_all_batch_ids(context: ApplicationContext) -> List[int]:
    """Get all batch IDs from database"""
    query = "SELECT id FROM ecod_schema.batch ORDER BY id"
    rows = context.db.execute_query(query)
    return [row[0] for row in rows]


def get_batch_id_for_process(context: ApplicationContext, process_id: int) -> Optional[int]:
    """Get batch ID for a process"""
    query = "SELECT batch_id FROM ecod_schema.process_status WHERE id = %s"
    rows = context.db.execute_query(query, (process_id,))
    return rows[0][0] if rows else None


def get_process_id_for_protein(context: ApplicationContext, pdb_id: str, chain_id: str, batch_id: int) -> Optional[int]:
    """Get process ID for a specific protein"""
    query = """
    SELECT ps.id
    FROM ecod_schema.process_status ps
    JOIN ecod_schema.protein p ON ps.protein_id = p.id
    WHERE p.pdb_id = %s AND p.chain_id = %s AND ps.batch_id = %s
    """

    try:
        rows = context.db.execute_query(query, (pdb_id, chain_id, batch_id))
        return rows[0][0] if rows else None
    except Exception as e:
        logging.getLogger("domain_partition").error(f"Error getting process ID: {e}")
        return None


def get_protein_for_process(context: ApplicationContext, process_id: int) -> Optional[Dict[str, Any]]:
    """Get protein information for a process ID"""
    query = """
    SELECT p.pdb_id, p.chain_id, ps.id as process_id
    FROM ecod_schema.process_status ps
    JOIN ecod_schema.protein p ON ps.protein_id = p.id
    WHERE ps.id = %s
    """

    rows = context.db.execute_dict_query(query, (process_id,))
    return rows[0] if rows else None


def find_domain_summary(service: DomainPartitionService, batch_path: str,
                       pdb_id: str, chain_id: str, reference: str,
                       blast_only: bool = False) -> Optional[str]:
    """Find domain summary file using service utilities"""
    # Use the service's internal method through its analyzer
    file_type = 'blast_only_summary' if blast_only else 'domain_summary'

    try:
        evidence_paths = get_all_evidence_paths(batch_path, pdb_id, chain_id, reference)

        if file_type in evidence_paths and evidence_paths[file_type]['exists_at']:
            return evidence_paths[file_type]['exists_at']

    except Exception as e:
        logging.getLogger("domain_partition").warning(f"Error finding domain summary: {e}")

    return None


def verify_batch_readiness(context: ApplicationContext, batch_id: int, blast_only: bool = False,
                          reps_only: bool = False) -> bool:
    """Verify that a batch has the necessary domain summaries"""
    summary_type = "blast_only_summary" if blast_only else "domain_summary"
    query = """
    SELECT
        COUNT(*) as total,
        SUM(CASE WHEN EXISTS (
            SELECT 1 FROM ecod_schema.process_file pf
            WHERE pf.process_id = ps.id
            AND pf.file_type = %s
            AND pf.file_exists = TRUE
        ) THEN 1 ELSE 0 END) as ready_count
    FROM ecod_schema.process_status ps
    WHERE ps.batch_id = %s
    """

    params = [summary_type, batch_id]

    if reps_only:
        query += " AND ps.is_representative = TRUE"

    results = context.db.execute_dict_query(query, tuple(params))

    if not results:
        return False

    total = results[0]["total"]
    ready = results[0]["ready_count"]

    # Check if at least 90% of proteins are ready
    return total > 0 and ready / total >= 0.9


def verify_domain_summary(batch_path: str, pdb_id: str, chain_id: str, reference: str,
                         blast_only: bool = False) -> bool:
    """Verify that domain summary file exists for a protein"""
    all_paths = get_all_evidence_paths(batch_path, pdb_id, chain_id, reference)

    summary_type = "blast_only_summary" if blast_only else "domain_summary"

    return (summary_type in all_paths and all_paths[summary_type]["exists_at"])


def find_protein_in_database(context: ApplicationContext, pdb_id: str, chain_id: str,
                            batch_id: Optional[int] = None) -> List[Dict[str, Any]]:
    """Find a protein in the database"""
    query = """
    SELECT ps.id as process_id, ps.batch_id, ps.current_stage, ps.status,
           ps.error_message, ps.is_representative
    FROM ecod_schema.process_status ps
    JOIN ecod_schema.protein p ON ps.protein_id = p.id
    WHERE p.pdb_id = %s AND p.chain_id = %s
    """

    params = [pdb_id, chain_id]

    if batch_id is not None:
        query += " AND ps.batch_id = %s"
        params.append(batch_id)

    query += " ORDER BY ps.batch_id DESC"

    rows = context.db.execute_dict_query(query, tuple(params))

    # Enhance with file information
    for row in rows:
        process_id = row["process_id"]

        file_query = """
        SELECT file_type, file_path, file_exists, file_size
        FROM ecod_schema.process_file
        WHERE process_id = %s
        """

        files = context.db.execute_dict_query(file_query, (process_id,))
        row["files"] = {file["file_type"]: file for file in files}

    return rows


#
# MONITOR MODE FUNCTIONS - UPDATED FOR SERVICE ARCHITECTURE
#

def monitor_batch_status(args: argparse.Namespace) -> int:
    """Monitor batch processing status using service"""
    logger = logging.getLogger("domain_partition.monitor")

    # Initialize context with config
    context = ApplicationContext(args.config)

    # Create service
    service = DomainPartitionService(context)

    # Validate batch ID
    batch_info = get_batch_info(context, args.batch_id)
    if not batch_info:
        logger.error(f"Batch {args.batch_id} not found")
        return 1

    logger.info(f"Monitoring batch {args.batch_id} ({batch_info.get('batch_name', '')}) " +
               f"every {args.interval} seconds")

    # Monitor until completion or timeout
    start_time = time.time()
    last_progress = None

    while True:
        # Get current progress using service's tracker
        current_progress = service.tracker.get_batch_progress(args.batch_id)

        # Display progress if changed
        if current_progress != last_progress:
            total = current_progress.get("total", 0)
            complete = current_progress.get("complete", 0)
            errors = current_progress.get("errors", 0)
            processing = current_progress.get("processing", 0)

            logger.info(f"Progress: {complete}/{total} complete ({current_progress.get('complete_pct', 0):.1f}%), " +
                       f"{processing} processing, {errors} errors")

            last_progress = current_progress

            # Check if complete
            if complete + errors >= total and total > 0:
                logger.info(f"Batch {args.batch_id} processing complete")
                return 0

        # Check timeout
        if args.timeout and time.time() - start_time > args.timeout:
            logger.warning(f"Monitoring timeout reached after {args.timeout} seconds")
            return 0

        # Wait for next check
        time.sleep(args.interval)


#
# CLEANUP MODE FUNCTIONS
#

def cleanup_non_representative_status(args: argparse.Namespace) -> int:
    """Clean up status for non-representative proteins using service"""

    logger = logging.getLogger("domain_partition.cleanup")

    # Initialize context
    context = ApplicationContext(args.config)

    # Create service
    service = DomainPartitionService(context)

    # Determine which batches to clean up
    batch_ids = args.batch_ids if hasattr(args, 'batch_ids') and args.batch_ids else []

    if not batch_ids:
        # Get all batches between 13-61 (the ones mentioned with filtering issues)
        batch_ids = list(range(13, 62))
        logger.info(f"Cleaning up status for historical batches 13-61")

    total_updated = 0

    for batch_id in batch_ids:
        try:
            # Use service's tracker to update non-representative status
            updated = service.tracker.update_non_representative_status(batch_id)
            total_updated += updated

            if updated > 0:
                logger.info(f"Batch {batch_id}: Updated {updated} non-representative proteins")

        except Exception as e:
            logger.error(f"Error cleaning up batch {batch_id}: {e}")

    logger.info(f"Cleanup complete: Updated {total_updated} non-representative proteins across {len(batch_ids)} batches")
    return 0


#
# MAIN FUNCTION
#

def main():
    """Main entry point with service-based processing"""
    # Create top-level parser
    parser = argparse.ArgumentParser(
        description="Domain Partition Tool - Service-Based Implementation"
    )
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')

    # Create subparsers for different modes
    subparsers = parser.add_subparsers(dest="mode", help="Operating mode")

    # PROCESS mode - for running domain partition
    process_parser = subparsers.add_parser("process", help="Process domains using service")
    process_subparsers = process_parser.add_subparsers(dest="action", help="Process action")

    # Process batch subparser
    batch_parser = process_subparsers.add_parser("batch", help="Process a batch")
    batch_parser.add_argument('--batch-id', type=int, required=True,
                            help='Batch ID to process')
    batch_parser.add_argument('--blast-only', action='store_true',
                            help='Use only BLAST results (no HHSearch)')
    batch_parser.add_argument('--limit', type=int,
                            help='Maximum number of proteins to process')
    batch_parser.add_argument('--reps-only', action='store_true',
                            help='Process only representative proteins')
    batch_parser.add_argument('--force', action='store_true',
                            help='Force processing even if batch is not ready')
    batch_parser.add_argument('--workers', type=int, default=4,
                            help='Number of parallel workers')

    # Process all batches subparser
    all_parser = process_subparsers.add_parser("all", help="Process multiple batches")
    all_parser.add_argument('--batch-ids', type=int, nargs='+',
                          help='Batch IDs to process (default: all batches)')
    all_parser.add_argument('--exclude-batch-ids', type=int, nargs='+', default=[],
                          help='Batch IDs to exclude')
    all_parser.add_argument('--reference', type=str,
                          help='Filter batches by reference version')
    all_parser.add_argument('--blast-only', action='store_true',
                          help='Use only BLAST results (no HHSearch)')
    all_parser.add_argument('--limit-per-batch', type=int,
                          help='Maximum number of proteins to process per batch')
    all_parser.add_argument('--reps-only', action='store_true',
                          help='Process only representative proteins')
    all_parser.add_argument('--batch-size', type=int, default=5,
                          help='Number of batches to process simultaneously')
    all_parser.add_argument('--wait-between-groups', type=int, default=30,
                          help='Seconds to wait between batch groups')
    all_parser.add_argument('--force', action='store_true',
                          help='Force processing even if batches are not ready')
    all_parser.add_argument('--workers', type=int, default=1,
                          help='Number of parallel workers per batch')

    # SLURM-specific options
    all_parser.add_argument('--use-slurm', action='store_true',
                          help='Submit jobs to SLURM instead of running directly')
    all_parser.add_argument('--slurm-threads', type=int, default=8,
                          help='Number of threads to request per SLURM job')
    all_parser.add_argument('--slurm-memory', type=str, default='16G',
                          help='Memory to request per SLURM job')
    all_parser.add_argument('--slurm-time', type=str, default='12:00:00',
                          help='Time limit for SLURM jobs')
    all_parser.add_argument('--wait-for-completion', action='store_true',
                          help='Wait for SLURM jobs to complete')
    all_parser.add_argument('--check-interval', type=int, default=60,
                          help='Seconds between job status checks')
    all_parser.add_argument('--timeout', type=int,
                          help='Maximum time to wait for completion (seconds)')

    # Process specific proteins subparser
    specific_parser = process_subparsers.add_parser("specific", help="Process specific proteins")
    specific_parser.add_argument('--process-ids', type=int, nargs='+', required=True,
                               help='Process IDs to process')
    specific_parser.add_argument('--batch-id', type=int,
                               help='Batch ID (optional)')
    specific_parser.add_argument('--blast-only', action='store_true',
                               help='Use only BLAST results')

    # Process single protein subparser
    single_parser = process_subparsers.add_parser("single", help="Process a single protein")
    single_parser.add_argument('--pdb-id', type=str, required=True,
                             help='PDB ID')
    single_parser.add_argument('--chain-id', type=str, required=True,
                             help='Chain ID')
    single_parser.add_argument('--batch-id', type=int,
                             help='Batch ID (optional)')
    single_parser.add_argument('--batch-path', type=str,
                             help='Batch path (required if batch ID not provided)')
    single_parser.add_argument('--reference', type=str,
                             help='Reference version (required if batch ID not provided)')
    single_parser.add_argument('--blast-only', action='store_true',
                             help='Use only BLAST results')

    # ANALYZE mode
    analyze_parser = subparsers.add_parser("analyze", help="Examine domain partition results")
    analyze_subparsers = analyze_parser.add_subparsers(dest="action", help="Analysis action")

    # Analyze batch status subparser
    status_parser = analyze_subparsers.add_parser("status", help="Check batch status")
    status_parser.add_argument('--batch-ids', type=int, nargs='+',
                             help='Batch IDs to check (default: all batches)')
    status_parser.add_argument('--blast-only', action='store_true',
                             help='Check blast-only status')

    # Analyze protein status subparser
    protein_parser = analyze_subparsers.add_parser("protein", help="Check protein status")
    protein_parser.add_argument('--pdb-id', type=str, required=True,
                              help='PDB ID')
    protein_parser.add_argument('--chain-id', type=str, required=True,
                              help='Chain ID')
    protein_parser.add_argument('--batch-id', type=int,
                              help='Batch ID (optional)')

    # Analyze domain counts subparser
    counts_parser = analyze_subparsers.add_parser("counts", help="Analyze domain statistics")
    counts_parser.add_argument('--batch-ids', type=int, nargs='+',
                             help='Batch IDs to analyze (default: all batches)')
    counts_parser.add_argument('--sample-size', type=int, default=50,
                             help='Number of proteins to sample per batch')

    # REPAIR mode
    repair_parser = subparsers.add_parser("repair", help="Fix or regenerate problematic files")
    repair_subparsers = repair_parser.add_subparsers(dest="action", help="Repair action")

    # Repair failed processes subparser
    failed_parser = repair_subparsers.add_parser("failed", help="Reset and retry failed processes")
    failed_parser.add_argument('--batch-id', type=int, required=True,
                             help='Batch ID to repair')
    failed_parser.add_argument('--rerun', action='store_true',
                             help='Rerun partition after resetting')
    failed_parser.add_argument('--blast-only', action='store_true',
                             help='Use only BLAST results for rerun')

    # Repair missing files subparser
    missing_parser = repair_subparsers.add_parser("missing", help="Regenerate missing files")
    missing_parser.add_argument('--batch-id', type=int, required=True,
                              help='Batch ID to repair')
    missing_parser.add_argument('--blast-only', action='store_true',
                              help='Use only BLAST results')
    missing_parser.add_argument('--limit', type=int,
                              help='Maximum number of files to repair')

    # MONITOR mode
    monitor_parser = subparsers.add_parser("monitor", help="Monitor domain partition status")
    monitor_subparsers = monitor_parser.add_subparsers(dest="action", help="Monitor action")

    # Monitor batch status subparser
    batch_monitor_parser = monitor_subparsers.add_parser("batch", help="Monitor batch processing")
    batch_monitor_parser.add_argument('--batch-id', type=int, required=True,
                                    help='Batch ID to monitor')
    batch_monitor_parser.add_argument('--interval', type=int, default=60,
                                    help='Check interval in seconds')
    batch_monitor_parser.add_argument('--timeout', type=int,
                                    help='Timeout in seconds')

    # CLEANUP mode
    cleanup_parser = subparsers.add_parser("cleanup", help="Clean up historical status tracking issues")
    cleanup_subparsers = cleanup_parser.add_subparsers(dest="action", help="Cleanup action")

    # Cleanup non-representative status
    nonrep_parser = cleanup_subparsers.add_parser("non-representative",
                                                 help="Fix non-representative protein status")
    nonrep_parser.add_argument('--batch-ids', type=int, nargs='+',
                              help='Specific batch IDs to clean up (default: 13-61)')

    # Parse arguments and run appropriate function
    args = parser.parse_args()

    # Setup logging
    logger = setup_logging(args.verbose, args.log_file)
    logger.info("Starting domain partition tool with service-based implementation")

    # Use absolute path for config if provided
    if args.config and not os.path.isabs(args.config):
        args.config = os.path.abspath(args.config)

    # Handle different modes
    if args.mode == "process":
        if args.action == "batch":
            return process_batch(args)
        elif args.action == "specific":
            return process_specific_proteins(args)
        elif args.action == "single":
            return process_single_protein(args)
        elif args.action == "all":
            return process_all_batches(args)
        else:
            logger.error(f"Unknown process action: {args.action}")
            return 1

    elif args.mode == "analyze":
        if args.action == "status":
            return analyze_batch_status(args)
        elif args.action == "protein":
            return analyze_protein_status(args)
        elif args.action == "counts":
            return analyze_domain_counts(args)
        else:
            logger.error(f"Unknown analyze action: {args.action}")
            return 1

    elif args.mode == "repair":
        if args.action == "failed":
            return repair_failed_processes(args)
        elif args.action == "missing":
            return repair_missing_files(args)
        else:
            logger.error(f"Unknown repair action: {args.action}")
            return 1

    elif args.mode == "monitor":
        if args.action == "batch":
            return monitor_batch_status(args)
        else:
            logger.error(f"Unknown monitor action: {args.action}")
            return 1

    elif args.mode == "cleanup":
        if args.action == "non-representative":
            return cleanup_non_representative_status(args)
        else:
            logger.error(f"Unknown cleanup action: {args.action}")
            return 1

    else:
        logger.error(f"Unknown mode: {args.mode}")
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
