#!/usr/bin/env python3
"""
validate_batch_directories.py - Validate and report on directory structure issues
in the ECOD batch directories, with options for remediation.
"""

import os
import sys
import argparse
import logging
import re
import json
import glob
import shutil
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional, Any
from collections import defaultdict, Counter
from datetime import datetime

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
# Uncomment to use the pyECOD utilities if needed
# from ecod.utils.path_utils import get_standardized_paths, get_file_db_path

# Standard directory structure definition
STANDARD_DIRS = {
    "fastas": {"required": True, "subdirs": []},
    "blast": {"required": True, "subdirs": ["chain", "domain"]},
    "hhsearch": {"required": True, "subdirs": ["profiles"]},
    "domains": {"required": True, "subdirs": []},
    "scripts": {"required": False, "subdirs": []},
}

# Standard file extension patterns
STANDARD_EXTS = {
    "fastas": [".fa"],
    "blast/chain": [".xml"],
    "blast/domain": [".xml"],
    "hhsearch": [".hhr", ".xml"],
    "hhsearch/profiles": [".a3m", ".hhm"],
    "domains": [".xml"],
}

# Non-standard directory names that should be normalized
NONSTANDARD_DIRS = {
    "query_fastas": "fastas",
    "fasta": "fastas",
    "blast_results": "blast",
    "hhsearch_results": "hhsearch",
    "domain_summaries": "domains",
}

def setup_logging(verbose: bool = False, log_file: Optional[str] = None) -> logging.Logger:
    """Set up logging with optional file output"""
    log_level = logging.DEBUG if verbose else logging.INFO
    
    logger = logging.getLogger("batch_validator")
    logger.setLevel(log_level)
    logger.handlers = []  # Clear any existing handlers
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)
    
    # File handler if specified
    if log_file:
        try:
            os.makedirs(os.path.dirname(log_file), exist_ok=True)
            file_handler = logging.FileHandler(log_file)
            file_handler.setLevel(log_level)
            file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            file_handler.setFormatter(file_formatter)
            logger.addHandler(file_handler)
        except Exception as e:
            logger.error(f"Failed to set up log file: {e}")
    
    return logger

def get_batch_directories(base_dir: str, pattern: Optional[str] = None) -> List[str]:
    """Find all batch directories in the base directory"""
    logger = logging.getLogger("batch_validator")
    
    if not os.path.exists(base_dir):
        logger.error(f"Base directory does not exist: {base_dir}")
        return []
    
    # Use the pattern if provided, otherwise find all possible batch directories
    if pattern:
        batch_pattern = os.path.join(base_dir, pattern)
        batch_dirs = glob.glob(batch_pattern)
    else:
        # Look for directories that match common batch naming patterns
        batch_dirs = []
        for item in os.listdir(base_dir):
            item_path = os.path.join(base_dir, item)
            if os.path.isdir(item_path) and (
                item.startswith("ecod_batch_") or
                item.startswith("batch_") or
                item.startswith("alt_rep_batch_")
            ):
                batch_dirs.append(item_path)
    
    logger.info(f"Found {len(batch_dirs)} batch directories")
    return sorted(batch_dirs)

def validate_directory_structure(batch_dir: str) -> Dict[str, Any]:
    """Validate the directory structure of a batch directory"""
    logger = logging.getLogger("batch_validator")
    
    logger.debug(f"Validating directory structure for: {batch_dir}")
    
    results = {
        "batch_dir": batch_dir,
        "batch_name": os.path.basename(batch_dir),
        "exists": os.path.exists(batch_dir),
        "is_dir": os.path.isdir(batch_dir),
        "standard_dirs": {},
        "nonstandard_dirs": [],
        "missing_dirs": [],
        "dir_count": 0,
        "file_count": 0,
        "total_size": 0,
        "file_extensions": defaultdict(list),
        "batch0_dirs": [],
        "errors": [],
    }
    
    if not results["exists"]:
        results["errors"].append(f"Batch directory does not exist: {batch_dir}")
        return results
    
    if not results["is_dir"]:
        results["errors"].append(f"Batch path is not a directory: {batch_dir}")
        return results
    
    # Get all directories in the batch directory
    try:
        all_items = os.listdir(batch_dir)
        
        dirs = []
        for item in all_items:
            item_path = os.path.join(batch_dir, item)
            if os.path.isdir(item_path):
                dirs.append(item)
                
                # Check for batch_0 subdirectories
                if item == "fastas":
                    batch0_path = os.path.join(batch_dir, item, "batch_0")
                    if os.path.isdir(batch0_path):
                        results["batch0_dirs"].append(f"fastas/batch_0")
        
        results["dir_count"] = len(dirs)
    except Exception as e:
        results["errors"].append(f"Failed to list directory contents: {str(e)}")
        return results
    
    # Check for standard directories
    for std_dir, specs in STANDARD_DIRS.items():
        if std_dir in dirs:
            # Directory exists
            results["standard_dirs"][std_dir] = {"exists": True, "subdirs": {}}
            
            # Check subdirectories
            for subdir in specs["subdirs"]:
                subdir_path = os.path.join(batch_dir, std_dir, subdir)
                results["standard_dirs"][std_dir]["subdirs"][subdir] = os.path.exists(subdir_path)
        else:
            # Directory doesn't exist
            results["standard_dirs"][std_dir] = {"exists": False, "subdirs": {}}
            if specs["required"]:
                results["missing_dirs"].append(std_dir)
    
    # Check for non-standard directories
    for dir_name in dirs:
        if dir_name not in STANDARD_DIRS and dir_name != "batch_0":
            # Check if it's a known non-standard directory
            if dir_name in NONSTANDARD_DIRS:
                results["nonstandard_dirs"].append({
                    "name": dir_name,
                    "standard_name": NONSTANDARD_DIRS[dir_name],
                    "type": "known"
                })
            else:
                results["nonstandard_dirs"].append({
                    "name": dir_name,
                    "standard_name": None,
                    "type": "unknown"
                })
    
    # Find all files and analyze extensions
    total_size = 0
    file_count = 0
    
    for root, _, files in os.walk(batch_dir):
        for file in files:
            file_path = os.path.join(root, file)
            try:
                file_size = os.path.getsize(file_path)
                total_size += file_size
                file_count += 1
                
                # Extract the relative path and file extension
                rel_path = os.path.relpath(os.path.dirname(file_path), batch_dir)
                _, ext = os.path.splitext(file)
                
                # Store in the correct category
                if rel_path == ".":
                    # Root directory
                    results["file_extensions"]["root"].append(ext)
                else:
                    results["file_extensions"][rel_path].append(ext)
            except Exception as e:
                results["errors"].append(f"Failed to process file {file_path}: {str(e)}")
    
    results["file_count"] = file_count
    results["total_size"] = total_size
    
    # Convert the file extension collections to Counters
    for key in results["file_extensions"]:
        results["file_extensions"][key] = dict(Counter(results["file_extensions"][key]))
    
    return results

def find_nonstandard_files(batch_dir: str) -> Dict[str, List[str]]:
    """Find files with non-standard extensions or in non-standard locations"""
    logger = logging.getLogger("batch_validator")
    
    logger.debug(f"Searching for non-standard files in: {batch_dir}")
    
    nonstandard_files = defaultdict(list)
    
    for root, _, files in os.walk(batch_dir):
        rel_dir = os.path.relpath(root, batch_dir)
        
        # Skip root directory since it has no standard files
        if rel_dir == ".":
            continue
        
        # Find the proper category for this directory
        category = None
        for std_dir in STANDARD_EXTS:
            if rel_dir == std_dir or rel_dir.startswith(f"{std_dir}/"):
                category = std_dir
                break
        
        # Skip if no category found
        if category is None:
            for file in files:
                nonstandard_files["unknown_location"].append(os.path.join(rel_dir, file))
            continue
        
        # Check file extensions
        standard_exts = STANDARD_EXTS.get(category, [])
        for file in files:
            _, ext = os.path.splitext(file)
            if ext not in standard_exts:
                nonstandard_files["wrong_extension"].append(os.path.join(rel_dir, file))
    
    return dict(nonstandard_files)

def check_batch_0_dirs(batch_dir: str) -> Dict[str, Any]:
    """Check for batch_0 subdirectories and analyze their contents"""
    batch0_path = os.path.join(batch_dir, "fastas", "batch_0")
    
    results = {
        "exists": os.path.exists(batch0_path),
        "file_count": 0,
        "extensions": Counter(),
        "fastas_files": [],
    }
    
    if not results["exists"]:
        return results
    
    # Count files and gather extensions
    for file in os.listdir(batch0_path):
        file_path = os.path.join(batch0_path, file)
        if os.path.isfile(file_path):
            results["file_count"] += 1
            _, ext = os.path.splitext(file)
            results["extensions"][ext] += 1
    
    # Check if there are any files in the main fastas directory
    fastas_dir = os.path.join(batch_dir, "fastas")
    if os.path.exists(fastas_dir):
        for file in os.listdir(fastas_dir):
            file_path = os.path.join(fastas_dir, file)
            if os.path.isfile(file_path) and not file.startswith("."):
                results["fastas_files"].append(file)
    
    return results

def generate_remediation_script(batch_dirs: List[str], output_path: str) -> str:
    """Generate a shell script to remediate directory structure issues"""
    logger = logging.getLogger("batch_validator")
    
    logger.info(f"Generating remediation script: {output_path}")
    
    # Create script content
    script_lines = [
        "#!/bin/bash",
        "",
        "# Auto-generated remediation script for batch directory structure",
        f"# Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "# Exit on error",
        "set -e",
        "",
        "# Dry run mode (set to 0 to actually perform changes)",
        "DRY_RUN=1",
        "",
        "# Function to display commands in dry run mode",
        "function run_cmd() {",
        "    if [ $DRY_RUN -eq 1 ]; then",
        "        echo \"[DRY RUN] Would execute: $@\"",
        "    else",
        "        echo \"Executing: $@\"",
        "        \"$@\"",
        "    fi",
        "}",
        "",
        "echo \"Starting batch directory structure remediation\"",
        "echo \"Dry run mode: $([ $DRY_RUN -eq 1 ] && echo 'ON' || echo 'OFF')\"",
        "echo",
        "",
    ]
    
    for batch_dir in batch_dirs:
        batch_name = os.path.basename(batch_dir)
        script_lines.extend([
            f"echo \"Processing batch: {batch_name}\"",
            f"cd \"{batch_dir}\"",
            "",
        ])
        
        # Check for and fix batch_0 directory
        batch0_check = f"if [ -d \"fastas/batch_0\" ]; then"
        script_lines.extend([
            "# Fix batch_0 directory (move files to fastas)",
            batch0_check,
            "    echo \"  Moving files from fastas/batch_0 to fastas\"",
            "    for file in fastas/batch_0/*; do",
            "        if [ -f \"$file\" ]; then",
            "            filename=$(basename \"$file\")",
            "            # Convert .fasta to .fa while moving",
            "            if [[ \"$filename\" == *.fasta ]]; then",
            "                newname=\"${filename%.fasta}.fa\"",
            "                run_cmd mv \"$file\" \"fastas/$newname\"",
            "            else",
            "                run_cmd mv \"$file\" \"fastas/\"",
            "            fi",
            "        fi",
            "    done",
            "    # Remove empty batch_0 directory",
            "    run_cmd rmdir \"fastas/batch_0\" 2>/dev/null || true",
            "fi",
            "",
        ])
        
        # Create standard directories if missing
        script_lines.extend([
            "# Create standard directories if missing",
            "for dir in fastas blast/chain blast/domain hhsearch/profiles domains; do",
            "    if [ ! -d \"$dir\" ]; then",
            "        run_cmd mkdir -p \"$dir\"",
            "        echo \"  Created directory: $dir\"",
            "    fi",
            "done",
            "",
        ])
        
        # Normalize non-standard directories
        script_lines.extend([
            "# Normalize non-standard directories",
            "if [ -d \"query_fastas\" ]; then",
            "    echo \"  Moving files from query_fastas to fastas\"",
            "    for file in query_fastas/*; do",
            "        if [ -f \"$file\" ]; then",
            "            filename=$(basename \"$file\")",
            "            # Convert .fasta to .fa while moving",
            "            if [[ \"$filename\" == *.fasta ]]; then",
            "                newname=\"${filename%.fasta}.fa\"",
            "                run_cmd mv \"$file\" \"fastas/$newname\"",
            "            else",
            "                run_cmd mv \"$file\" \"fastas/\"",
            "            fi",
            "        fi",
            "    done",
            "    run_cmd rmdir \"query_fastas\" 2>/dev/null || true",
            "fi",
            "",
            "if [ -d \"fasta\" ]; then",
            "    echo \"  Moving files from fasta to fastas\"",
            "    for file in fasta/*; do",
            "        if [ -f \"$file\" ]; then",
            "            filename=$(basename \"$file\")",
            "            # Convert .fasta to .fa while moving",
            "            if [[ \"$filename\" == *.fasta ]]; then",
            "                newname=\"${filename%.fasta}.fa\"",
            "                run_cmd mv \"$file\" \"fastas/$newname\"",
            "            else",
            "                run_cmd mv \"$file\" \"fastas/\"",
            "            fi",
            "        fi",
            "    done",
            "    run_cmd rmdir \"fasta\" 2>/dev/null || true",
            "fi",
            "",
        ])
        
        # Fix file extensions
        script_lines.extend([
            "# Fix file extensions",
            "echo \"  Fixing file extensions\"",
            "# Convert .fasta to .fa in fastas directory",
            "for file in fastas/*.fasta; do",
            "    if [ -f \"$file\" ]; then",
            "        newname=\"${file%.fasta}.fa\"",
            "        run_cmd mv \"$file\" \"$newname\"",
            "    fi",
            "done",
            "",
        ])
        
        script_lines.append("cd - > /dev/null\n")
    
    # Add final message
    script_lines.extend([
        "echo",
        "echo \"Remediation complete!\"",
        "echo \"Set DRY_RUN=0 at the top of this script to perform actual changes.\"",
        "",
    ])
    
    # Write script to file
    try:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, "w") as f:
            f.write("\n".join(script_lines))
        os.chmod(output_path, 0o755)  # Make executable
        logger.info(f"Remediation script created: {output_path}")
        return output_path
    except Exception as e:
        logger.error(f"Failed to write remediation script: {str(e)}")
        return None

def create_report(results: List[Dict[str, Any]], output_file: str) -> None:
    """Create a JSON report file with validation results"""
    logger = logging.getLogger("batch_validator")
    
    try:
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        with open(output_file, "w") as f:
            json.dump(results, f, indent=2)
        logger.info(f"Report written to: {output_file}")
    except Exception as e:
        logger.error(f"Failed to write report: {str(e)}")

def display_summary(results: List[Dict[str, Any]]) -> None:
    """Display a summary of validation results"""
    logger = logging.getLogger("batch_validator")
    
    total_batches = len(results)
    batches_with_errors = sum(1 for r in results if r["errors"])
    batches_with_missing_dirs = sum(1 for r in results if r["missing_dirs"])
    batches_with_nonstandard_dirs = sum(1 for r in results if r["nonstandard_dirs"])
    batches_with_batch0 = sum(1 for r in results if r["batch0_dirs"])
    
    total_files = sum(r["file_count"] for r in results)
    total_size = sum(r["total_size"] for r in results)
    
    logger.info("=" * 60)
    logger.info("BATCH STRUCTURE VALIDATION SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Total batches examined: {total_batches}")
    logger.info(f"Batches with errors: {batches_with_errors}")
    logger.info(f"Batches with missing directories: {batches_with_missing_dirs}")
    logger.info(f"Batches with non-standard directories: {batches_with_nonstandard_dirs}")
    logger.info(f"Batches with batch_0 subdirectories: {batches_with_batch0}")
    logger.info(f"Total files: {total_files}")
    logger.info(f"Total size: {format_size(total_size)}")
    logger.info("=" * 60)
    
    # List batches with issues
    if batches_with_batch0 > 0:
        logger.info("\nBatches with batch_0 subdirectories:")
        for r in results:
            if r["batch0_dirs"]:
                logger.info(f"  - {r['batch_name']}")
    
    if batches_with_missing_dirs > 0:
        logger.info("\nBatches with missing directories:")
        for r in results:
            if r["missing_dirs"]:
                logger.info(f"  - {r['batch_name']}: {', '.join(r['missing_dirs'])}")
    
    if batches_with_nonstandard_dirs > 0:
        logger.info("\nBatches with non-standard directories:")
        for r in results:
            if r["nonstandard_dirs"]:
                nonstandard = [d["name"] for d in r["nonstandard_dirs"]]
                logger.info(f"  - {r['batch_name']}: {', '.join(nonstandard)}")
    
    logger.info("\nTop 5 batches by file count:")
    top_by_files = sorted(results, key=lambda r: r["file_count"], reverse=True)[:5]
    for r in top_by_files:
        logger.info(f"  - {r['batch_name']}: {r['file_count']} files")
    
    logger.info("\nTop 5 batches by size:")
    top_by_size = sorted(results, key=lambda r: r["total_size"], reverse=True)[:5]
    for r in top_by_size:
        logger.info(f"  - {r['batch_name']}: {format_size(r['total_size'])}")

def format_size(size_bytes: int) -> str:
    """Format a size in bytes to a human-readable string"""
    if size_bytes < 1024:
        return f"{size_bytes} B"
    elif size_bytes < 1024 * 1024:
        return f"{size_bytes / 1024:.2f} KB"
    elif size_bytes < 1024 * 1024 * 1024:
        return f"{size_bytes / (1024 * 1024):.2f} MB"
    else:
        return f"{size_bytes / (1024 * 1024 * 1024):.2f} GB"

def main():
    parser = argparse.ArgumentParser(description="Validate ECOD batch directory structure")
    parser.add_argument("--base-dir", default="/data/ecod/pdb_updates/batches",
                      help="Base directory containing batch directories")
    parser.add_argument("--pattern", help="Pattern to match batch directories (e.g., 'ecod_batch_*')")
    parser.add_argument("--batch-id", type=str, help="Single batch ID or comma-separated list of batch IDs to validate")
    parser.add_argument("--report", default="validation_report.json", 
                      help="Path to output JSON report file")
    parser.add_argument("--script", default="remediate_batches.sh",
                      help="Path to output remediation script")
    parser.add_argument("--generate-script", action="store_true",
                      help="Generate a remediation script")
    parser.add_argument("--log-file", help="Path to log file")
    parser.add_argument("--verbose", "-v", action="store_true", 
                      help="Enable verbose output")
    
    args = parser.parse_args()
    
    # Set up logging
    logger = setup_logging(args.verbose, args.log_file)
    
    # Get batch directories
    if args.batch_id:
        # Process specific batch IDs
        batch_ids = args.batch_id.split(",")
        batch_dirs = []
        for bid in batch_ids:
            # Try to find matching batch directories
            pattern = f"*{bid}*" if bid.isdigit() else f"*{bid}*"
            matches = glob.glob(os.path.join(args.base_dir, pattern))
            if matches:
                batch_dirs.extend(matches)
            else:
                logger.warning(f"No batch directory found matching ID: {bid}")
    else:
        # Find all batch directories
        batch_dirs = get_batch_directories(args.base_dir, args.pattern)
    
    if not batch_dirs:
        logger.error("No batch directories found to validate")
        return 1
    
    logger.info(f"Validating {len(batch_dirs)} batch directories")
    
    # Validate each batch directory
    results = []
    for i, batch_dir in enumerate(batch_dirs):
        logger.info(f"[{i+1}/{len(batch_dirs)}] Validating: {os.path.basename(batch_dir)}")
        validation_result = validate_directory_structure(batch_dir)
        
        # Check for batch_0 directories
        batch0_result = check_batch_0_dirs(batch_dir)
        validation_result["batch0_details"] = batch0_result
        
        # Find non-standard files
        nonstandard_files = find_nonstandard_files(batch_dir)
        validation_result["nonstandard_files"] = nonstandard_files
        
        results.append(validation_result)
    
    # Generate report
    create_report(results, args.report)
    
    # Display summary
    display_summary(results)
    
    # Generate remediation script if requested
    if args.generate_script:
        script_path = generate_remediation_script(batch_dirs, args.script)
        if script_path:
            logger.info(f"\nRemediation script generated: {script_path}")
            logger.info("Review the script and run it with:")
            logger.info(f"  bash {script_path}")
            logger.info("Note: The script is in dry run mode by default. Set DRY_RUN=0 to perform actual changes.")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
