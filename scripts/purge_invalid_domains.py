#!/usr/bin/env python3
"""
purge_invalid_domains.py - Identify and remove invalid domain files

This script checks domain files in a batch for validity, identifying files with
non-standard roots or structural issues, and optionally purging them from the filesystem
and database.

Usage:
  python purge_invalid_domains.py --batch-id 68 --config config/config.yml --dry-run
  python purge_invalid_domains.py --batch-id 68 --config config/config.yml --purge
"""

import os
import sys
import argparse
import logging
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple, Set

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.config import ConfigManager
from ecod.core.context import ApplicationContext
from ecod.db import DBManager
from ecod.utils.path_utils import get_all_evidence_paths, find_files_with_legacy_paths

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
    
    return logging.getLogger("purge_domains")

def get_batch_info(db, batch_id):
    """Get batch information"""
    query = """
    SELECT id, batch_name, base_path, ref_version
    FROM ecod_schema.batch
    WHERE id = %s
    """
    
    batch_info = db.execute_dict_query(query, (batch_id,))
    return batch_info[0] if batch_info else None

def find_domains_files(batch_path: str, pattern="*.domains.xml") -> List[str]:
    """Find all domain files in the batch directory"""
    domains_dir = os.path.join(batch_path, "domains")
    if not os.path.exists(domains_dir):
        return []
    
    return [str(p) for p in Path(domains_dir).glob(pattern)]

def get_protein_info(db, pdb_id: str, chain_id: str, batch_id: int) -> Optional[Dict[str, Any]]:
    """Get protein and process information from the database"""
    query = """
    SELECT 
        p.id as protein_id, ps.id as process_id, ps.current_stage
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    WHERE 
        p.pdb_id = %s AND p.chain_id = %s AND ps.batch_id = %s
    """
    
    results = db.execute_dict_query(query, (pdb_id, chain_id, batch_id))
    return results[0] if results else None

def get_domain_file_record(db, process_id: int) -> Optional[Dict[str, Any]]:
    """Get domain partition file information from the database"""
    query = """
    SELECT 
        id, file_type, file_path, file_exists
    FROM 
        ecod_schema.process_file
    WHERE 
        process_id = %s AND 
        (file_type = 'domain_partition' OR file_type = 'blast_only_partition')
    """
    
    results = db.execute_dict_query(query, (process_id,))
    return results[0] if results else None

def update_file_record(db, file_id: int, file_exists: bool = False, dry_run: bool = True):
    """Update file record to mark as non-existent"""
    if dry_run:
        return True
    
    try:
        db.update(
            "ecod_schema.process_file",
            {"file_exists": file_exists},
            "id = %s",
            (file_id,)
        )
        return True
    except Exception as e:
        logging.error(f"Failed to update file record: {str(e)}")
        return False

def reset_process_status(db, process_id: int, stage: str = "domain_summary_complete", dry_run: bool = True):
    """Reset process status to a previous stage"""
    if dry_run:
        return True
    
    try:
        db.update(
            "ecod_schema.process_status",
            {"current_stage": stage},
            "id = %s",
            (process_id,)
        )
        return True
    except Exception as e:
        logging.error(f"Failed to update process status: {str(e)}")
        return False

def check_domain_file_validity(file_path: str) -> Tuple[bool, str, Dict[str, Any]]:
    """
    Check if a domain file is valid
    
    Returns:
        Tuple of (is_valid, error_message, details)
    """
    details = {
        "root_element": None,
        "domain_count": 0,
        "has_domain_list": False,
    }
    
    try:
        # Try to parse the XML
        tree = ET.parse(file_path)
        root = tree.getroot()
        
        # Store root element name
        details["root_element"] = root.tag
        
        # Check for expected root element
        if root.tag != "domain_doc":
            return False, f"Root element is not domain_doc (found {root.tag})", details
        
        # Check for domain_list element
        domain_list = root.find("domain_list")
        if domain_list is None:
            return False, "No domain_list element", details
        
        details["has_domain_list"] = True
        
        # Check for at least one domain
        domains = domain_list.findall("domain")
        domain_count = len(domains)
        details["domain_count"] = domain_count
        
        if domain_count == 0:
            return False, "No domains found in domain_list", details
        
        # Basic check for expected attributes in each domain
        for domain in domains:
            if 'range' not in domain.attrib:
                return False, f"Domain missing 'range' attribute", details
        
        # Check for statistics section
        stats = root.find("statistics")
        if stats is None:
            return False, "No statistics element", details
        
        # File passed all checks
        return True, "Valid domain file", details
        
    except ET.ParseError as e:
        return False, f"XML parsing error: {str(e)}", details
    except Exception as e:
        return False, f"Validation error: {str(e)}", details

def purge_invalid_domains(context, batch_id: int, criteria: List[str] = None, 
                        dry_run: bool = True) -> Dict[str, Any]:
    """
    Purge invalid domain files from a batch
    
    Args:
        context: Application context
        batch_id: Batch ID
        criteria: List of validation criteria to check 
                 (default: ['root_element', 'no_domains', 'parse_error'])
        dry_run: Whether to perform a dry run
        
    Returns:
        Dictionary with purge statistics
    """
    logger = logging.getLogger("purge_domains")
    
    # Default criteria if not specified
    if criteria is None:
        criteria = ['root_element', 'no_domains', 'parse_error']
    
    # Get batch info
    batch_info = get_batch_info(context.db, batch_id)
    if not batch_info:
        logger.error(f"Batch {batch_id} not found")
        return {"error": "Batch not found"}
    
    batch_path = batch_info['base_path']
    ref_version = batch_info['ref_version']
    
    logger.info(f"Looking for invalid domain files in batch {batch_id} ({batch_info['batch_name']})")
    logger.info(f"Base path: {batch_path}")
    logger.info(f"Using criteria: {', '.join(criteria)}")
    
    # Find all domain files
    domain_files = find_domains_files(batch_path)
    logger.info(f"Found {len(domain_files)} domain files")
    
    # Results storage
    results = {
        "total_files": len(domain_files),
        "invalid_files": 0,
        "purged_files": 0,
        "failed_purges": 0,
        "details": [],
        "reason_counts": {},
        "reset_process_status": 0
    }
    
    # Process each file
    for file_path in domain_files:
        # Extract PDB and chain ID from filename
        filename = os.path.basename(file_path)
        parts = filename.split(".")
        
        if len(parts) < 3 or "_" not in parts[0]:
            logger.warning(f"Could not parse filename: {filename}")
            continue
        
        pdb_chain = parts[0].split("_")
        if len(pdb_chain) != 2:
            logger.warning(f"Could not parse PDB/chain: {parts[0]}")
            continue
        
        pdb_id = pdb_chain[0]
        chain_id = pdb_chain[1]
        
        # Check file validity
        valid, message, details = check_domain_file_validity(file_path)
        
        # Determine if file should be purged based on criteria
        should_purge = False
        purge_reason = None
        
        if not valid:
            if 'root_element' in criteria and details["root_element"] != "domain_doc":
                should_purge = True
                purge_reason = f"Invalid root element: {details['root_element']}"
            elif 'no_domains' in criteria and details["has_domain_list"] and details["domain_count"] == 0:
                should_purge = True
                purge_reason = "No domains in domain_list"
            elif 'parse_error' in criteria and "parsing error" in message.lower():
                should_purge = True
                purge_reason = f"Parse error: {message}"
        
        # If file is invalid but no criteria match, log but don't purge
        if not valid and not should_purge:
            logger.info(f"Invalid file but criteria not met for purge: {file_path} - {message}")
            
        # If file should be purged, proceed with purging
        if should_purge:
            results["invalid_files"] += 1
            
            # Count reasons
            reason_type = purge_reason.split(":")[0] if ":" in purge_reason else purge_reason
            if reason_type not in results["reason_counts"]:
                results["reason_counts"][reason_type] = 0
            results["reason_counts"][reason_type] += 1
            
            # Get database information
            protein_info = get_protein_info(context.db, pdb_id, chain_id, batch_id)
            
            file_record = None
            if protein_info:
                file_record = get_domain_file_record(context.db, protein_info["process_id"])
            
            purge_result = {
                "pdb_id": pdb_id,
                "chain_id": chain_id,
                "file_path": file_path,
                "reason": purge_reason,
                "details": details,
                "db_updated": False,
                "status_reset": False,
                "file_deleted": False
            }
            
            # Perform the actual purge if not in dry run mode
            if not dry_run:
                # Delete file
                try:
                    os.remove(file_path)
                    purge_result["file_deleted"] = True
                    logger.info(f"Deleted file: {file_path}")
                except Exception as e:
                    logger.error(f"Failed to delete file {file_path}: {str(e)}")
                    purge_result["error"] = str(e)
                    results["failed_purges"] += 1
                    
                # Update database if record found
                if file_record:
                    # Mark file as non-existent
                    updated = update_file_record(context.db, file_record["id"], False, dry_run)
                    purge_result["db_updated"] = updated
                    
                    # Reset process status if it's domain_partition_complete
                    if protein_info["current_stage"] == "domain_partition_complete":
                        reset = reset_process_status(context.db, protein_info["process_id"], "domain_summary_complete", dry_run)
                        purge_result["status_reset"] = reset
                        
                        if reset:
                            results["reset_process_status"] += 1
            else:
                logger.info(f"Would purge: {file_path} - {purge_reason}")
            
            # Only count as purged if we would actually delete the file
            results["purged_files"] += 1
            results["details"].append(purge_result)
    
    # Log summary
    if dry_run:
        logger.info(f"\nDry run summary - would purge {results['purged_files']} invalid domain files")
    else:
        logger.info(f"\nPurged {results['purged_files']} invalid domain files")
    
    logger.info(f"Total domain files: {results['total_files']}")
    logger.info(f"Invalid files: {results['invalid_files']}")
    
    if results["reason_counts"]:
        logger.info("\nInvalid file reasons:")
        for reason, count in results["reason_counts"].items():
            logger.info(f"  {reason}: {count}")
    
    if not dry_run:
        logger.info(f"\nDatabase updates:")
        logger.info(f"  Process status resets: {results['reset_process_status']}")
    
    return results

def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Purge invalid domain files')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to process')
    
    # Purge operation options
    parser.add_argument('--dry-run', action='store_true',
                      help='Perform a dry run (no actual changes)')
    parser.add_argument('--purge', action='store_true',
                      help='Actually purge the invalid files')
    parser.add_argument('--criteria', type=str, nargs='+', 
                      choices=['root_element', 'no_domains', 'parse_error', 'all'],
                      default=['root_element', 'no_domains', 'parse_error'],
                      help='Criteria for identifying invalid files')
    
    # Output options
    parser.add_argument('--output', type=str,
                      help='Output JSON file for purge results')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')

    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(args.verbose, args.log_file)
    
    # Verify conflicting arguments
    if args.dry_run and args.purge:
        logger.error("Cannot specify both --dry-run and --purge. Use --purge to actually delete files.")
        return 1
    
    # Default to dry run if neither specified
    if not args.dry_run and not args.purge:
        args.dry_run = True
        logger.info("Neither --dry-run nor --purge specified, defaulting to --dry-run")
    
    # Use all criteria if 'all' is specified
    if 'all' in args.criteria:
        args.criteria = ['root_element', 'no_domains', 'parse_error']
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Run the purge operation
    dry_run = not args.purge
    results = purge_invalid_domains(context, args.batch_id, args.criteria, dry_run)
    
    # Write results to file if specified
    if args.output and 'error' not in results:
        import json
        try:
            with open(args.output, 'w') as f:
                json.dump(results, f, indent=2)
            logger.info(f"Wrote purge results to {args.output}")
        except Exception as e:
            logger.error(f"Failed to write output file: {str(e)}")
    
    # Return success if no error
    return 0 if 'error' not in results else 1

if __name__ == "__main__":
    sys.exit(main())
