#!/usr/bin/env python3
"""
hhsearch_tools.py - Unified HHSearch processing and analysis toolkit

This script provides a comprehensive set of tools for working with HHSearch results:
- process: Convert HHR files to XML and create domain summaries
- collate: Combine HHSearch results with BLAST evidence
- analyze: Examine and diagnose HHSearch results
- repair: Fix missing or problematic files
- batches: Process multiple batches in parallel

Each mode has specific options and can work with different backends (database or filesystem).
"""

import os
import sys
import logging
import argparse
import glob
import random
import re
import xml.etree.ElementTree as ET
from xml.dom import minidom
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple, Union
from concurrent.futures import ThreadPoolExecutor, as_completed

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

# Import core modules
from ecod.core.context import ApplicationContext
from ecod.config import ConfigManager

# Import processing modules
from ecod.pipelines.hhsearch.processor import HHSearchProcessor, HHRToXMLConverter
from ecod.pipelines.domain_analysis.hhresult_registrar import HHResultRegistrar
from ecod.utils.hhsearch_utils import HHRParser

# Import path utilities
from ecod.utils.path_utils import (
    get_standardized_paths,
    get_all_evidence_paths,
    get_file_db_path,
    resolve_file_path,
    find_files_with_legacy_paths,
    migrate_file_to_standard_path,
    scan_batch_directory
)

def setup_logging(verbose: bool = False, log_file: Optional[str] = None) -> None:
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

#
# ANALYSIS MODE FUNCTIONS
#

def count_hits_in_hhr(hhr_file: str) -> int:
    """Count the number of hits in an HHR file"""
    logger = logging.getLogger("hhsearch.analyze.hhr_counter")

    try:
        with open(hhr_file, 'r') as f:
            content = f.read()

        lines = content.split('\n')
        hit_count = 0

        # Find the table header
        table_start = None
        for i, line in enumerate(lines):
            if line.startswith(' No Hit'):
                table_start = i + 1
                break

        if not table_start:
            return 0

        # Count hit entries
        for i in range(table_start, len(lines)):
            line = lines[i].strip()
            if line and line[0].isdigit() and not line.startswith('Q ') and not line.startswith('T '):
                hit_count += 1

        return hit_count
    except Exception as e:
        logger.warning(f"Error reading HHR file {hhr_file}: {str(e)}")
        return 0

def check_xml_content(xml_file: str) -> Tuple[int, List[str]]:
    """Check the content of an XML file and return hit count and info"""
    logger = logging.getLogger("hhsearch.analyze.xml_checker")

    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()

        # Look for hh_hit elements in hh_hit_list
        hit_list = root.find(".//hh_hit_list")
        if hit_list is None:
            return 0, []

        hits = hit_list.findall("hh_hit")

        hit_info = []
        for hit in hits:
            hit_num = hit.get("hit_num", "unknown")
            hit_id = hit.get("hit_id", "unknown")
            probability = hit.get("probability", "unknown")
            hit_info.append(f"{hit_num}: {hit_id} (prob: {probability})")

        return len(hits), hit_info
    except Exception as e:
        logger.warning(f"Error parsing XML file {xml_file}: {str(e)}")
        return 0, []

def check_summary_content(summary_file: str) -> Tuple[int, Dict[str, Any]]:
    """Check the content of a domain summary file"""
    logger = logging.getLogger("hhsearch.analyze.summary_checker")

    try:
        tree = ET.parse(summary_file)
        root = tree.getroot()

        result = {"has_chain_blast": False, "has_domain_blast": False, "has_hhsearch": False}

        # Look for chain blast hits
        chain_blast = root.find(".//chain_blast_evidence")
        if chain_blast is not None:
            hits = chain_blast.findall(".//hit")
            result["chain_blast_hits"] = len(hits)
            result["has_chain_blast"] = len(hits) > 0

        # Look for domain blast hits
        domain_blast = root.find(".//domain_blast_evidence")
        if domain_blast is not None:
            hits = domain_blast.findall(".//hit")
            result["domain_blast_hits"] = len(hits)
            result["has_domain_blast"] = len(hits) > 0

        # Look for hh_hit elements in hhsearch_evidence section
        hh_evidence = root.find(".//hhsearch_evidence")
        if hh_evidence is None:
            result["hhsearch_hits"] = 0
            return 0, result

        hit_list = hh_evidence.find(".//hh_hit_list")
        if hit_list is None:
            result["hhsearch_hits"] = 0
            return 0, result

        hits = hit_list.findall("hh_hit")
        result["hhsearch_hits"] = len(hits)
        result["has_hhsearch"] = len(hits) > 0

        return len(hits), result
    except Exception as e:
        logger.warning(f"Error parsing summary file {summary_file}: {str(e)}")
        return 0, {}

def analyze_files(batch_path: str, ref_version: str = "develop291", sample_size: int = 5,
                verbose: bool = False
) -> Dict[str, Any]:
    """
    Examine HHR, XML, and summary files to diagnose issues

    Args:
        batch_path: Path to batch directory
        ref_version: Reference version
        sample_size: Number of files to sample for detailed analysis
        verbose: Enable verbose logging

    Returns:
        Dictionary with analysis results
    """
    logger = logging.getLogger("hhsearch.analyze")

    # Scan batch directory to find all relevant files
    logger.info(f"Scanning batch directory: {batch_path}")
    files = scan_batch_directory(batch_path, ref_version)

    hhr_files = files['hhr']
    xml_files = files['hh_xml']
    summary_files = files['domain_summary'] + files['domain_partition']

    logger.info(f"Found {len(hhr_files)} HHR files")
    logger.info(f"Found {len(xml_files)} HHSearch XML files")
    logger.info(f"Found {len(summary_files)} domain summary files")

    # Check if counts match
    if len(hhr_files) != len(xml_files):
        logger.warning(f"Number of HHR files ({len(hhr_files)}) doesn't match number of XML files ({len(xml_files)})")
    else:
        logger.info(f"Number of HHR files matches number of XML files: {len(hhr_files)}")

    # Select random samples for detailed analysis
    if len(hhr_files) > sample_size:
        sample_hhr_files = random.sample(hhr_files, sample_size)
    else:
        sample_hhr_files = hhr_files

    analysis_results = {
        "total_files": {
            "hhr": len(hhr_files),
            "xml": len(xml_files),
            "summary": len(summary_files)
        },
        "hhr_xml_match": len(hhr_files) == len(xml_files),
        "detailed_analysis": []
    }

    logger.info(f"Performing detailed analysis on {len(sample_hhr_files)} sample files")

    # Analyze each sample
    for hhr_file in sample_hhr_files:
        filename = os.path.basename(hhr_file)
        pdb_chain = filename.split('.')[0]  # Format: pdbid_chain

        # Find corresponding XML and summary files
        xml_file = None
        summary_file = None

        for xml in xml_files:
            if os.path.basename(xml).startswith(pdb_chain):
                xml_file = xml
                break

        for summary in summary_files:
            if os.path.basename(summary).startswith(pdb_chain):
                summary_file = summary
                break

        # Analyze HHR file
        hhr_hits = count_hits_in_hhr(hhr_file)

        # Analyze XML file
        xml_hits = 0
        hit_info = []
        if xml_file:
            xml_hits, hit_info = check_xml_content(xml_file)

        # Analyze summary file
        summary_hits = 0
        summary_info = {}
        if summary_file:
            summary_hits, summary_info = check_summary_content(summary_file)

        # Log results
        logger.info(f"{pdb_chain}: {hhr_hits} HHR hits, {xml_hits} XML hits, {summary_hits} summary hits")

        if hhr_hits != xml_hits:
            logger.warning(f"{pdb_chain}: HHR and XML hit counts don't match ({hhr_hits} != {xml_hits})")

        if xml_hits != summary_hits and summary_file:
            logger.warning(f"{pdb_chain}: XML and summary hit counts don't match ({xml_hits} != {summary_hits})")

        # Print top 5 hits from XML if available
        if hit_info and verbose:
            top_hits = hit_info[:5]
            logger.info(f"{pdb_chain} top 5 hits: {', '.join(top_hits)}")

        # Add to analysis results
        analysis_results["detailed_analysis"].append({
            "pdb_chain": pdb_chain,
            "hhr_file": hhr_file,
            "xml_file": xml_file,
            "summary_file": summary_file,
            "hhr_hits": hhr_hits,
            "xml_hits": xml_hits,
            "summary_hits": summary_hits,
            "top_hits": hit_info[:5] if hit_info else [],
            "has_chain_blast": summary_info.get("has_chain_blast", False),
            "has_domain_blast": summary_info.get("has_domain_blast", False),
            "has_hhsearch": summary_info.get("has_hhsearch", False)
        })

    # Look for patterns in the results
    hhr_hits_avg = sum(r["hhr_hits"] for r in analysis_results["detailed_analysis"]) / len(analysis_results["detailed_analysis"]) if analysis_results["detailed_analysis"] else 0
    xml_hits_avg = sum(r["xml_hits"] for r in analysis_results["detailed_analysis"]) / len(analysis_results["detailed_analysis"]) if analysis_results["detailed_analysis"] else 0
    summary_hits_avg = sum(r["summary_hits"] for r in analysis_results["detailed_analysis"]) / len(analysis_results["detailed_analysis"]) if analysis_results["detailed_analysis"] else 0

    analysis_results["averages"] = {
        "hhr_hits": hhr_hits_avg,
        "xml_hits": xml_hits_avg,
        "summary_hits": summary_hits_avg
    }

    logger.info(f"Analysis complete. Average hits - HHR: {hhr_hits_avg:.1f}, XML: {xml_hits_avg:.1f}, Summary: {summary_hits_avg:.1f}")

    return analysis_results

def find_missing_files(batch_path: str, ref_version: str = "develop291") -> Dict[str, List[str]]:
    """
    Find missing HHR, XML, or summary files

    Args:
        batch_path: Path to batch directory
        ref_version: Reference version

    Returns:
        Dictionary with lists of missing files
    """
    logger = logging.getLogger("hhsearch.analyze.missing_files")

    # Scan batch directory to find all relevant files
    files = scan_batch_directory(batch_path, ref_version)

    # Extract PDB_CHAIN from filenames for each type
    hhr_chains = set()
    xml_chains = set()
    summary_chains = set()

    for hhr_file in files['hhr']:
        filename = os.path.basename(hhr_file)
        parts = filename.split('.')
        pdb_chain = parts[0]  # Format: pdbid_chain
        hhr_chains.add(pdb_chain)

    for xml_file in files['hh_xml']:
        filename = os.path.basename(xml_file)
        parts = filename.split('.')
        pdb_chain = parts[0]  # Format: pdbid_chain
        xml_chains.add(pdb_chain)

    for summary_file in files['domain_summary'] + files['domain_partition']:
        filename = os.path.basename(summary_file)
        parts = filename.split('.')
        pdb_chain = parts[0]  # Format: pdbid_chain
        summary_chains.add(pdb_chain)

    # Find chains with missing files
    missing_xml = hhr_chains - xml_chains
    missing_summary = hhr_chains - summary_chains

    # Find chains with XML but no HHR (unusual case)
    xml_no_hhr = xml_chains - hhr_chains

    # Log findings
    logger.info(f"Found {len(missing_xml)} chains with missing XML files")
    logger.info(f"Found {len(missing_summary)} chains with missing summary files")

    if xml_no_hhr:
        logger.warning(f"Found {len(xml_no_hhr)} chains with XML files but no HHR files")

    return {
        "missing_xml": sorted(list(missing_xml)),
        "missing_summary": sorted(list(missing_summary)),
        "xml_no_hhr": sorted(list(xml_no_hhr))
    }

def analyze_mode(args: argparse.Namespace) -> int:
    """
    Run analysis mode

    Args:
        args: Command line arguments

    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("hhsearch.analyze")
    logger.info(f"Running analysis mode for batch at {args.batch_path}")

    # Run appropriate analysis function based on subcommand
    if args.action == "content":
        analyze_files(args.batch_path, args.ref_version, args.sample_size, args.verbose)
    elif args.action == "missing":
        find_missing_files(args.batch_path, args.ref_version)
    else:
        logger.error(f"Unknown analysis action: {args.action}")
        return 1

    return 0

#
# PROCESS MODE FUNCTIONS
#

def process_via_database(args: argparse.Namespace) -> int:
    """
    Process HHSearch results using the database-backed pipeline

    Args:
        args: Command line arguments

    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("hhsearch.process.db")

    # Initialize application context with config path
    context = ApplicationContext(args.config)

    # Get batch info from database
    batch_query = """
    SELECT id, batch_name, base_path, ref_version
    FROM ecod_schema.batch
    WHERE id = %s
    """

    batch_info = context.db.execute_dict_query(batch_query, (args.batch_id,))
    if not batch_info:
        logger.error(f"Batch {args.batch_id} not found")
        return 1

    batch_path = batch_info[0]['base_path']
    ref_version = batch_info[0]['ref_version']
    batch_name = batch_info[0]['batch_name']

    logger.info(f"Processing batch {args.batch_id} ({batch_name}) with reference {ref_version}")

    # Create processor and run it
    processor = HHSearchProcessor(context)
    processed_count = processor.process_batch(args.batch_id, args.force)

    if processed_count > 0:
        logger.info(f"Successfully processed {processed_count} chains")
        return 0
    else:
        logger.warning("No chains were processed")
        return 1

def register_via_database(args: argparse.Namespace) -> int:
    """
    Register HHSearch results using HHResultRegistrar

    Args:
        args: Command line arguments

    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("hhsearch.register")

    # Initialize application context with config path
    context = ApplicationContext(args.config)

    # Create registrar
    registrar = HHResultRegistrar(context)

    try:
        # Process batch
        if hasattr(args, 'chains') and args.chains:
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

def parse_hhr_file(hhr_file: str, logger: logging.Logger) -> Optional[Dict[str, Any]]:
    """
    Parse an HHR file using the HHRParser class

    Args:
        hhr_file: Path to HHR file
        logger: Logger instance

    Returns:
        Dictionary with parsed HHR data or None if parsing failed
    """
    parser = HHRParser(logger)
    return parser.parse(hhr_file)

def process_via_filesystem(args: argparse.Namespace) -> int:
    """
    Process HHSearch results directly from the filesystem

    Args:
        args: Command line arguments

    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("hhsearch.process.fs")

    batch_path = args.batch_path
    ref_version = args.ref_version

    logger.info(f"Processing HHSearch results from filesystem at {batch_path}")

    # Find all HHR files
    hhr_pattern = os.path.join(batch_path, "hhsearch", f"*.{ref_version}.hhr")
    hhr_files = glob.glob(hhr_pattern)

    logger.info(f"Found {len(hhr_files)} HHR files on disk")

    if not hhr_files:
        logger.warning(f"No HHR files found matching pattern: {hhr_pattern}")
        return 1

    if args.limit:
        hhr_files = hhr_files[:args.limit]
        logger.info(f"Limited processing to {args.limit} files")

    # Create necessary directories
    hhsearch_dir = os.path.join(batch_path, "hhsearch")
    domains_dir = os.path.join(batch_path, "domains")

    os.makedirs(hhsearch_dir, exist_ok=True)
    os.makedirs(domains_dir, exist_ok=True)

    # Initialize parser and converter
    parser = HHRParser(logger)
    converter = HHRToXMLConverter(logger)
    #collator = DomainEvidenceCollator(logger)  #This is now deprecated!

    # Process each HHR file
    processed_count = 0
    for hhr_file in hhr_files:
        # Extract PDB and chain ID from filename
        filename = os.path.basename(hhr_file)
        parts = filename.split('.')
        pdb_chain = parts[0]  # Format: pdbid_chain
        pdb_parts = pdb_chain.split('_')

        if len(pdb_parts) != 2:
            logger.warning(f"Invalid filename format: {filename}")
            continue

        pdb_id = pdb_parts[0]
        chain_id = pdb_parts[1]

        # Get standardized file paths
        paths = get_standardized_paths(batch_path, pdb_id, chain_id, ref_version)

        # Skip if domain summary already exists and not forcing
        if os.path.exists(paths['domain_summary']) and not args.force:
            logger.debug(f"Domain summary already exists for {pdb_chain}, skipping")
            continue

        try:
            # Parse HHR file using the parser class method
            logger.debug(f"Parsing HHR file: {hhr_file}")
            hhr_data = parser.parse(hhr_file)

            if not hhr_data:
                logger.warning(f"Failed to parse HHR file: {hhr_file}")
                continue

            # Convert to XML using the converter class method
            logger.debug(f"Converting HHR data to XML for {pdb_chain}")
            xml_string = converter.convert(hhr_data, pdb_id, chain_id, ref_version)

            if not xml_string:
                logger.warning(f"Failed to convert HHR data to XML for {pdb_chain}")
                continue

            # Save HHSearch XML using the converter's save method
            logger.debug(f"Saving HHSearch XML: {paths['hh_xml']}")
            if not converter.save(xml_string, paths['hh_xml']):
                logger.warning(f"Failed to save HHSearch XML: {paths['hh_xml']}")
                continue

            # Create simple domain summary with HHSearch data only
            # (no BLAST evidence as we're working from filesystem)
            root = ET.Element("domain_summ_doc")

            # Add metadata
            metadata = ET.SubElement(root, "metadata")
            ET.SubElement(metadata, "pdb_id").text = pdb_id
            ET.SubElement(metadata, "chain_id").text = chain_id
            ET.SubElement(metadata, "reference").text = ref_version
            ET.SubElement(metadata, "creation_date").text = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

            # Add empty evidence sections
            ET.SubElement(root, "chain_blast_evidence")
            ET.SubElement(root, "domain_blast_evidence")

            # Add HHSearch evidence
            hhsearch_elem = ET.SubElement(root, "hhsearch_evidence")
            hhsearch_data = ET.parse(paths['hh_xml'])

            hh_doc = hhsearch_data.find("hh_summ_doc")
            if hh_doc is not None:
                for child in hh_doc:
                    hhsearch_elem.append(ET.fromstring(ET.tostring(child)))

            # Add empty domain suggestions section
            ET.SubElement(root, "domain_suggestions")

            # Convert to string with pretty formatting
            rough_string = ET.tostring(root, 'utf-8')
            reparsed = minidom.parseString(rough_string)
            pretty_xml = reparsed.toprettyxml(indent="  ")

            # Save domain summary
            logger.debug(f"Saving domain summary: {paths['domain_summary']}")
            try:
                with open(paths['domain_summary'], 'w', encoding='utf-8') as f:
                    f.write(pretty_xml)
            except Exception as e:
                logger.warning(f"Failed to save domain summary: {paths['domain_summary']}, error: {str(e)}")
                continue

            processed_count += 1

            if processed_count % 10 == 0:
                logger.info(f"Processed {processed_count} chains so far")

        except Exception as e:
            logger.error(f"Error processing {pdb_chain}: {str(e)}")
            continue

    logger.info(f"Successfully processed {processed_count} chains")
    return 0 if processed_count > 0 else 1

def process_mode(args: argparse.Namespace) -> int:
    """
    Run process mode with appropriate backend

    Args:
        args: Command line arguments

    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("hhsearch.process")

    # Choose between database or filesystem backend
    if args.backend == "db":
        return process_via_database(args)
    elif args.backend == "fs":
        return process_via_filesystem(args)
    elif args.backend == "register":
        return register_via_database(args)
    else:
        logger.error(f"Unknown backend: {args.backend}")
        return 1

#
# COLLATE MODE FUNCTIONS
#

def collate_with_blast(args: argparse.Namespace) -> int:
    """
    Collate HHSearch results with BLAST evidence

    Args:
        args: Command line arguments

    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("hhsearch.collate_w_blast")

    # Initialize application context
    context = ApplicationContext(args.config)

    # Get batch info
    batch_query = """
    SELECT id, batch_name, base_path, ref_version
    FROM ecod_schema.batch
    WHERE id = %s
    """

    batch_info = context.db.execute_dict_query(batch_query, (args.batch_id,))
    if not batch_info:
        logger.error(f"Batch {args.batch_id} not found")
        return 1

    base_path = batch_info[0]['base_path']
    ref_version = batch_info[0]['ref_version']
    batch_name = batch_info[0]['batch_name']

    logger.info(f"Collating HHSearch results with BLAST for batch {args.batch_id} ({batch_name})")

    # Get representative proteins with HHSearch results
    protein_query = """
    SELECT
        p.id as protein_id, p.pdb_id, p.chain_id, ps.id as process_id
    FROM
        ecod_schema.process_status ps
    JOIN
        ecod_schema.protein p ON ps.protein_id = p.id
    LEFT JOIN
        ecod_schema.process_file pf ON ps.id = pf.process_id AND pf.file_type = 'hhr'
    WHERE
        ps.batch_id = %s
        AND ps.is_representative = TRUE
        AND pf.file_exists = TRUE
    ORDER BY
        p.pdb_id, p.chain_id
    """

    proteins = context.db.execute_dict_query(protein_query, (args.batch_id,))
    logger.info(f"Found {len(proteins)} representative proteins with HHSearch results")

    if args.limit and args.limit < len(proteins):
        proteins = proteins[:args.limit]
        logger.info(f"Limited to {args.limit} proteins")

    # Initialize HHSearch processor
    processor = HHSearchProcessor(context)

    # Process each protein
    success_count = 0
    for protein in proteins:
        pdb_id = protein['pdb_id']
        chain_id = protein['chain_id']
        process_id = protein['process_id']

        logger.info(f"Processing {pdb_id}_{chain_id}")

        # Get standardized paths
        paths = get_standardized_paths(base_path, pdb_id, chain_id, ref_version, create_dirs=True)

        # Check if full pipeline domain partition exists
        if os.path.exists(paths['domain_partition']) and not args.force:
            # Verify this is a full pipeline summary (contains HHSearch evidence)
            try:
                tree = ET.parse(paths['domain_partition'])
                root = tree.getroot()

                # Look for HHSearch evidence section
                hhsearch_elem = root.find(".//hhsearch_evidence")
                if hhsearch_elem is not None:
                    # Check if it has any hit elements
                    hits = hhsearch_elem.findall(".//hh_hit")
                    if len(hits) > 0:
                        logger.info(f"Full pipeline domain summary already exists for {pdb_id}_{chain_id}, skipping")
                        success_count += 1
                        continue
            except Exception as e:
                logger.warning(f"Error checking domain partition: {str(e)}")

            logger.info(f"Found incomplete summary for {pdb_id}_{chain_id}, replacing with full pipeline version")

        # Process the chain using HHSearchProcessor
        result = processor._process_chain(
            pdb_id,
            chain_id,
            process_id,
            batch_info[0],
            ref_version,
            args.force
        )

        if result:
            success_count += 1
            logger.info(f"Successfully processed {pdb_id}_{chain_id}")
        else:
            logger.warning(f"Failed to process {pdb_id}_{chain_id}")

    logger.info(f"Successfully collated results for {success_count}/{len(proteins)} proteins")
    return 0 if success_count > 0 else 1

def collate_mode(args: argparse.Namespace) -> int:
    """
    Run collate mode

    Args:
        args: Command line arguments

    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("hhsearch.collate")
    logger.info("Running collate mode")

    return collate_with_blast(args)

#
# REPAIR MODE FUNCTIONS
#

def repair_missing_files(args: argparse.Namespace) -> int:
    """
    Repair missing HHR, XML, or summary files

    Args:
        args: Command line arguments

    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("hhsearch.repair")

    batch_path = args.batch_path
    ref_version = args.ref_version

    logger.info(f"Repairing missing files for batch at {batch_path}")

    # Find missing files
    missing = find_missing_files(batch_path, ref_version)

    if not missing['missing_xml'] and not missing['missing_summary']:
        logger.info("No missing files found to repair")
        return 0

    logger.info(f"Found {len(missing['missing_xml'])} chains with missing XML files")
    logger.info(f"Found {len(missing['missing_summary'])} chains with missing summary files")

    # Initialize parser and converter
    parser = HHRParser(logger)
    converter = HHRToXMLConverter(logger)

    # Repair missing XML files
    xml_repaired = 0
    if missing['missing_xml']:
        logger.info(f"Repairing missing XML files")

        for pdb_chain in missing['missing_xml']:
            if '_' not in pdb_chain:
                logger.warning(f"Invalid PDB chain format: {pdb_chain}")
                continue

            pdb_id, chain_id = pdb_chain.split('_')

            # Get standardized paths
            paths = get_standardized_paths(batch_path, pdb_id, chain_id, ref_version)

            # Find HHR file
            hhr_file = None
            for f in glob.glob(os.path.join(batch_path, "hhsearch", f"{pdb_chain}*.hhr")):
                hhr_file = f
                break

            if not hhr_file:
                logger.warning(f"Could not find HHR file for {pdb_chain}")
                continue

            try:
                # Parse HHR file
                hhr_data = parser.parse(hhr_file)
                if not hhr_data:
                    logger.warning(f"Failed to parse HHR file for {pdb_chain}")
                    continue

                # Convert to XML
                xml_string = converter.convert(hhr_data, pdb_id, chain_id, ref_version)
                if not xml_string:
                    logger.warning(f"Failed to convert HHR data to XML for {pdb_chain}")
                    continue

                # Save XML file
                if converter.save(xml_string, paths['hh_xml']):
                    xml_repaired += 1
                    logger.info(f"Repaired XML file for {pdb_chain}")
                else:
                    logger.warning(f"Failed to save XML file for {pdb_chain}")
            except Exception as e:
                logger.error(f"Error repairing XML for {pdb_chain}: {str(e)}")

    # Repair missing summary files
    summary_repaired = 0
    if missing['missing_summary']:
        logger.info(f"Repairing missing summary files")

        for pdb_chain in missing['missing_summary']:
            if '_' not in pdb_chain:
                logger.warning(f"Invalid PDB chain format: {pdb_chain}")
                continue

            pdb_id, chain_id = pdb_chain.split('_')

            # Get standardized paths
            paths = get_standardized_paths(batch_path, pdb_id, chain_id, ref_version)

            # Check if XML file exists
            if not os.path.exists(paths['hh_xml']):
                logger.warning(f"XML file missing for {pdb_chain}, skipping summary creation")
                continue

            try:
                # Create simple domain summary with HHSearch data only
                root = ET.Element("domain_summ_doc")

                # Add metadata
                metadata = ET.SubElement(root, "metadata")
                ET.SubElement(metadata, "pdb_id").text = pdb_id
                ET.SubElement(metadata, "chain_id").text = chain_id
                ET.SubElement(metadata, "reference").text = ref_version
                ET.SubElement(metadata, "creation_date").text = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

                # Add empty evidence sections
                ET.SubElement(root, "chain_blast_evidence")
                ET.SubElement(root, "domain_blast_evidence")

                # Add HHSearch evidence
                hhsearch_elem = ET.SubElement(root, "hhsearch_evidence")

                try:
                    hhsearch_data = ET.parse(paths['hh_xml'])
                    hh_doc = hhsearch_data.find("hh_summ_doc")
                    if hh_doc is not None:
                        for child in hh_doc:
                            hhsearch_elem.append(ET.fromstring(ET.tostring(child)))
                except Exception as e:
                    logger.warning(f"Error parsing XML for {pdb_chain}: {str(e)}")

                # Add empty domain suggestions section
                ET.SubElement(root, "domain_suggestions")

                # Convert to string with pretty formatting
                rough_string = ET.tostring(root, 'utf-8')
                reparsed = minidom.parseString(rough_string)
                pretty_xml = reparsed.toprettyxml(indent="  ")

                # Save domain summary
                with open(paths['domain_summary'], 'w', encoding='utf-8') as f:
                    f.write(pretty_xml)

                summary_repaired += 1
                logger.info(f"Repaired summary file for {pdb_chain}")
            except Exception as e:
                logger.error(f"Error repairing summary for {pdb_chain}: {str(e)}")

    logger.info(f"Repair summary: {xml_repaired} XML files and {summary_repaired} summary files repaired")
    return 0

def fix_database_paths(args: argparse.Namespace) -> int:
    """
    Fix file paths in the database to match standardized paths

    Args:
        args: Command line arguments

    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("hhsearch.repair.db_paths")

    # Initialize application context
    context = ApplicationContext(args.config)

    logger.info(f"Fixing database paths for batch {args.batch_id}")

    # Get batch info
    batch_query = """
    SELECT id, batch_name, base_path, ref_version
    FROM ecod_schema.batch
    WHERE id = %s
    """

    batch_info = context.db.execute_dict_query(batch_query, (args.batch_id,))
    if not batch_info:
        logger.error(f"Batch {args.batch_id} not found")
        return 1

    base_path = batch_info[0]['base_path']
    ref_version = batch_info[0]['ref_version']

    # Get all process files for this batch
    query = """
    SELECT
        pf.id, pf.process_id, pf.file_type, pf.file_path, pf.file_exists,
        p.pdb_id, p.chain_id
    FROM
        ecod_schema.process_file pf
    JOIN
        ecod_schema.process_status ps ON pf.process_id = ps.id
    JOIN
        ecod_schema.protein p ON ps.protein_id = p.id
    WHERE
        ps.batch_id = %s
    """

    process_files = context.db.execute_dict_query(query, (args.batch_id,))
    logger.info(f"Found {len(process_files)} files in database for batch {args.batch_id}")

    # Statistics
    stats = {
        'total': len(process_files),
        'updated': 0,
        'already_standard': 0,
        'errors': 0,
        'file_missing': 0
    }

    # Process each file
    for file in process_files:
        try:
            pdb_id = file['pdb_id']
            chain_id = file['chain_id']
            file_type = file['file_type']
            current_path = file['file_path']

            # Resolve current absolute path
            current_abs_path = resolve_file_path(base_path, current_path)

            # Get standardized paths
            paths = get_standardized_paths(base_path, pdb_id, chain_id, ref_version, create_dirs=False)

            # Skip if file type not in standardized paths
            if file_type not in paths:
                logger.warning(f"Unknown file type '{file_type}' for {pdb_id}_{chain_id}")
                stats['errors'] += 1
                continue

            # Get standard path
            standard_path = paths[file_type]
            standard_rel_path = get_file_db_path(base_path, standard_path)

            # Skip if already using standard path
            if current_path == standard_rel_path:
                logger.debug(f"Already using standard path: {current_path}")
                stats['already_standard'] += 1
                continue

            # Check if file exists at current path
            if not os.path.exists(current_abs_path):
                logger.warning(f"File does not exist at current path: {current_abs_path}")

                # Check if file exists at standard path
                if os.path.exists(standard_path):
                    logger.info(f"File exists at standard path: {standard_path}")
                else:
                    logger.warning(f"File missing at both current and standard paths")
                    stats['file_missing'] += 1
                    continue
            else:
                # Migrate file to standard path if needed
                if not os.path.exists(standard_path):
                    if not args.dry_run:
                        migrate_file_to_standard_path(current_abs_path, standard_path)

            # Update database path
            if not args.dry_run:
                context.db.update(
                    "ecod_schema.process_file",
                    {
                        "file_path": standard_rel_path
                    },
                    "id = %s",
                    (file['id'],)
                )
                logger.info(f"Updated path in database: {current_path} -> {standard_rel_path}")
            else:
                logger.info(f"Would update path: {current_path} -> {standard_rel_path}")

            stats['updated'] += 1

        except Exception as e:
            logger.error(f"Error processing file {file['id']}: {str(e)}")
            stats['errors'] += 1

    # Log statistics
    logger.info("Update Statistics:")
    logger.info(f"Total files processed: {stats['total']}")
    logger.info(f"Files already using standard paths: {stats['already_standard']}")
    logger.info(f"Files updated: {stats['updated']}")
    logger.info(f"Files missing: {stats['file_missing']}")
    logger.info(f"Errors: {stats['errors']}")

    return 0

def create_empty_summaries(args: argparse.Namespace) -> int:
    """
    Create empty summary files for peptides or no-hits proteins

    Args:
        args: Command line arguments

    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("hhsearch.repair.empty_summaries")

    # Initialize application context
    if hasattr(args, 'config') and args.config:
        context = ApplicationContext(args.config)

        # Get batch info from database
        batch_query = """
        SELECT id, batch_name, base_path, ref_version
        FROM ecod_schema.batch
        WHERE id = %s
        """

        batch_info = context.db.execute_dict_query(batch_query, (args.batch_id,))
        if not batch_info:
            logger.error(f"Batch {args.batch_id} not found")
            return 1

        batch_path = batch_info[0]['base_path']
        ref_version = batch_info[0]['ref_version']

        # Import the regenerate_missing_summaries functionality
        try:
            from regenerate_missing_summaries import process_missing_files

            # Run the process_missing_files function
            total, created, updated = process_missing_files(
                context,
                args.batch_id,
                args.dry_run,
                None  # Use default reference version
            )

            return 0 if created > 0 or args.dry_run else 1

        except ImportError:
            logger.error("Could not import regenerate_missing_summaries module")
            return 1
    else:
        logger.error("Config file required for creating empty summaries")
        return 1

def repair_mode(args: argparse.Namespace) -> int:
    """
    Run repair mode with appropriate action

    Args:
        args: Command line arguments

    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("hhsearch.repair")

    if args.action == "missing":
        return repair_missing_files(args)
    elif args.action == "db_paths":
        return fix_database_paths(args)
    elif args.action == "empty_summaries":
        return create_empty_summaries(args)
    elif args.action == "fix_summaries":
        return fix_domain_summaries(args)
    else:
        logger.error(f"Unknown repair action: {args.action}")
        return 1

def fix_domain_summaries(args: argparse.Namespace) -> int:
    """
    Fix domain summaries to include HHSearch hits

    Args:
        args: Command line arguments

    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("hhsearch.repair.fix_summaries")

    batch_path = args.batch_path
    ref_version = args.ref_version

    logger.info(f"Fixing domain summaries to include HHSearch hits in {batch_path}")

    # Find all HHSearch XML files
    hhsearch_dir = os.path.join(batch_path, "hhsearch")
    xml_pattern = os.path.join(hhsearch_dir, f"*.{ref_version}.hhsearch.xml")

    xml_files = glob.glob(xml_pattern)
    logger.info(f"Found {len(xml_files)} HHSearch XML files")

    if not xml_files:
        logger.warning(f"No HHSearch XML files found in {hhsearch_dir}")
        return 1

    # Limit processing if requested
    if args.limit and args.limit < len(xml_files):
        xml_files = xml_files[:args.limit]
        logger.info(f"Limited processing to {args.limit} files")

    # Process each XML file
    fixed_count = 0
    for xml_file in xml_files:
        # Extract pdb_chain from filename
        filename = os.path.basename(xml_file)
        parts = filename.split('.')
        pdb_chain = parts[0]  # Format: pdbid_chain

        # Fix domain summary
        success = fix_single_domain_summary(batch_path, pdb_chain, ref_version)

        if success:
            fixed_count += 1

        if fixed_count % 50 == 0 and fixed_count > 0:
            logger.info(f"Fixed {fixed_count} domain summaries so far")

    logger.info(f"Successfully fixed {fixed_count} domain summaries")
    return 0 if fixed_count > 0 else 1

def fix_single_domain_summary(batch_path: str, pdb_chain: str, ref_version: str) -> bool:
    """
    Fix a single domain summary to include HHSearch hits

    Args:
        batch_path: Path to batch directory
        pdb_chain: PDB chain ID (format: pdbid_chain)
        ref_version: Reference version

    Returns:
        True if successful, False otherwise
    """
    logger = logging.getLogger("hhsearch.repair.fix_single")

    hhsearch_dir = os.path.join(batch_path, "hhsearch")
    domains_dir = os.path.join(batch_path, "domains")

    # Define file paths
    xml_path = os.path.join(hhsearch_dir, f"{pdb_chain}.{ref_version}.hhsearch.xml")
    summary_path = os.path.join(domains_dir, f"{pdb_chain}.{ref_version}.domains.xml")

    if not os.path.exists(xml_path):
        logger.warning(f"HHSearch XML file not found: {xml_path}")
        return False

    if not os.path.exists(summary_path):
        # Try alternative summary path formats
        alt_paths = [
            os.path.join(domains_dir, f"{pdb_chain}.domains.xml"),
            os.path.join(domains_dir, f"{pdb_chain}.domain_summary.xml"),
            os.path.join(domains_dir, f"{pdb_chain}.{ref_version}.domain_summary.xml")
        ]

        for path in alt_paths:
            if os.path.exists(path):
                summary_path = path
                logger.debug(f"Found alternative summary path: {summary_path}")
                break

        if not os.path.exists(summary_path):
            logger.warning(f"Domain summary file not found: {summary_path}")
            return False

    try:
        # Parse XML files
        hh_tree = ET.parse(xml_path)
        hh_root = hh_tree.getroot()

        summary_tree = ET.parse(summary_path)
        summary_root = summary_tree.getroot()

        # Extract pdb_id and chain_id from summary
        pdb_id = summary_root.find(".//pdb_id").text if summary_root.find(".//pdb_id") is not None else "unknown"
        chain_id = summary_root.find(".//chain_id").text if summary_root.find(".//chain_id") is not None else "unknown"

        # Find hhsearch_evidence section in summary
        hhsearch_elem = summary_root.find(".//hhsearch_evidence")
        if hhsearch_elem is None:
            # Create hhsearch_evidence element if it doesn't exist
            logger.debug(f"Creating hhsearch_evidence section for {pdb_chain}")
            hhsearch_elem = ET.SubElement(summary_root, "hhsearch_evidence")
        else:
            # Clear existing content
            for child in list(hhsearch_elem):
                hhsearch_elem.remove(child)

        # Copy hit_list from HHSearch XML
        hit_list = hh_root.find(".//hh_hit_list")
        if hit_list is None:
            logger.warning(f"No hit_list found in HHSearch XML for {pdb_chain}")
            return False

        # Create a new hit_list in the hhsearch_evidence section
        new_hit_list = ET.SubElement(hhsearch_elem, "hh_hit_list")

        # Copy all hits
        hit_count = 0
        for hit in hit_list.findall("hh_hit"):
            # Create a copy of the hit element
            hit_string = ET.tostring(hit)
            new_hit = ET.fromstring(hit_string)
            new_hit_list.append(new_hit)
            hit_count += 1

        logger.info(f"Copied {hit_count} hits to domain summary for {pdb_chain}")

        # Save updated summary
        rough_string = ET.tostring(summary_root, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        pretty_xml = reparsed.toprettyxml(indent="  ")

        with open(summary_path, 'w', encoding='utf-8') as f:
            f.write(pretty_xml)

        logger.debug(f"Updated domain summary saved: {summary_path}")
        return True

    except Exception as e:
        logger.error(f"Error fixing domain summary for {pdb_chain}: {str(e)}")
        return False

def collate_single_batch(batch_id: int, config_path: str, force: bool = False) -> Tuple[int, int]:
    """
    Collate a single batch with BLAST evidence

    Args:
        batch_id: Batch ID to process
        config_path: Path to configuration file
        force: Force reprocessing of already processed results

    Returns:
        Tuple of (batch_id, number of proteins processed)
    """
    logger = logging.getLogger(f"collate_batch_{batch_id}")
    logger.info(f"Collating HHSearch results for batch {batch_id}")

    # Initialize application context
    context = ApplicationContext(config_path)

    # Create processor
    processor = HHSearchProcessor(context)

    try:
        # Get batch info
        batch_query = """
        SELECT id, batch_name, base_path, ref_version
        FROM ecod_schema.batch
        WHERE id = %s
        """
        batch_info = context.db.execute_dict_query(batch_query, (batch_id,))
        if not batch_info:
            logger.error(f"Batch {batch_id} not found")
            return batch_id, 0

        # Get proteins with HHSearch results
        protein_query = """
        SELECT
            p.id as protein_id, p.pdb_id, p.chain_id, ps.id as process_id
        FROM
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        LEFT JOIN
            ecod_schema.process_file pf ON ps.id = pf.process_id AND pf.file_type = 'hhr'
        WHERE
            ps.batch_id = %s
            AND ps.is_representative = TRUE
            AND pf.file_exists = TRUE
        """
        proteins = context.db.execute_dict_query(protein_query, (batch_id,))
        logger.info(f"Found {len(proteins)} proteins with HHSearch results")

        # Process each protein
        success_count = 0
        for protein in proteins:
            # Process the chain
            result = processor._process_chain(
                protein['pdb_id'],
                protein['chain_id'],
                protein['process_id'],
                batch_info[0],
                batch_info[0]['ref_version'],
                force
            )

            if result:
                success_count += 1

        logger.info(f"Successfully collated results for {success_count}/{len(proteins)} proteins in batch {batch_id}")
        return batch_id, success_count
    except Exception as e:
        logger.error(f"Error collating batch {batch_id}: {e}")
        logger.exception("Stack trace:")
        return batch_id, -1

def parallel_collate_batches(args: argparse.Namespace) -> int:
    """
    Collate HHSearch results with BLAST evidence across batches in parallel using SLURM

    Args:
        args: Command line arguments

    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("hhsearch.parallel_collate")
    logger.info(f"Parallel collating batches using SLURM")

    # Initialize application context
    context = ApplicationContext(args.config)

    # Get batch IDs to process
    if args.batch_ids:
        # Use specified batch IDs
        batch_ids = [b_id for b_id in args.batch_ids if b_id not in args.exclude_batch_ids]
    else:
        # Get all batch IDs from database
        all_batches = get_all_batch_ids(context)
        batch_ids = [b_id for b_id, _ in all_batches if b_id not in args.exclude_batch_ids]

    logger.info(f"Will collate {len(batch_ids)} batches: {batch_ids}")

    # Create a job manager for SLURM
    job_manager = context.job_manager

    # Create a temporary directory for job scripts
    temp_dir = os.path.join(context.config.get('output_dir', '/tmp'), f"collate_jobs_{int(time.time())}")
    os.makedirs(temp_dir, exist_ok=True)

    # Prepare job submissions
    job_ids = []

    for batch_id in batch_ids:
        # Create a command to run collation for this batch
        # Note: Use the collate mode of the same script for individual batch processing
        cmd = (
            f"python {os.path.abspath(sys.argv[0])} collate "
            f"--config {args.config} "
            f"--batch-id {batch_id} "
        )

        if args.force:
            cmd += "--force "

        if args.verbose:
            cmd += "--verbose "

        # Create job name and output path
        job_name = f"collate_batch_{batch_id}"
        output_dir = os.path.join(temp_dir, f"batch_{batch_id}")
        os.makedirs(output_dir, exist_ok=True)

        # Add log file if specified
        if args.log_file:
            batch_log = os.path.join(os.path.dirname(args.log_file), f"batch_{batch_id}_{os.path.basename(args.log_file)}")
            cmd += f"--log-file {batch_log} "

        # Create a job script
        script_path = job_manager.create_job_script(
            commands=[cmd],
            job_name=job_name,
            output_dir=output_dir,
            threads=args.threads or 4,
            memory=args.memory or "8G",
            time=args.time or "12:00:00"
        )

        # Submit the job
        job_id = job_manager.submit_job(script_path)

        if job_id:
            job_ids.append(job_id)
            logger.info(f"Submitted collation job for batch {batch_id}, SLURM job ID: {job_id}")
        else:
            logger.error(f"Failed to submit collation job for batch {batch_id}")

    logger.info(f"Submitted {len(job_ids)} collation jobs to SLURM")

    # Monitor jobs if requested
    if args.wait:
        logger.info(f"Waiting for jobs to complete, checking every {args.check_interval} seconds")

        # Record start time
        start_time = time.time()

        while job_ids:
            time.sleep(args.check_interval)

            # Check each job
            completed_jobs = []

            for job_id in job_ids:
                status = job_manager.check_job_status(job_id)

                if status in ["COMPLETED", "COMPLETING"]:
                    logger.info(f"Job {job_id} completed successfully")
                    completed_jobs.append(job_id)
                elif status in ["FAILED", "TIMEOUT", "CANCELLED", "NODE_FAIL"]:
                    logger.error(f"Job {job_id} failed with status {status}")
                    completed_jobs.append(job_id)

            # Remove completed jobs from tracking
            for job_id in completed_jobs:
                job_ids.remove(job_id)

            # Log status
            if job_ids:
                logger.info(f"Waiting for {len(job_ids)} remaining jobs")

            # Check timeout
            if args.timeout and time.time() - start_time > args.timeout:
                logger.warning(f"Timeout reached after {args.timeout} seconds")
                break

        # Final check
        if not job_ids:
            logger.info("All collation jobs completed")
        else:
            logger.warning(f"{len(job_ids)} jobs still running at end of monitoring period")
            return 1

    return 0 if job_ids else 1
#
# BATCH MODE FUNCTIONS
#

def get_all_batch_ids(context: ApplicationContext) -> List[Tuple[int, str]]:
    """Get all batch IDs from the database"""
    query = "SELECT id, batch_name FROM ecod_schema.batch ORDER BY id DESC"
    batches = context.db.execute_query(query)
    return [(row[0], row[1]) for row in batches]

def process_batch(batch_id: int, config_path: str, force: bool = False) -> Tuple[int, int]:
    """
    Process a single batch

    Args:
        batch_id: Batch ID to process
        config_path: Path to configuration file
        force: Force reprocessing of already processed results

    Returns:
        Tuple of (batch_id, number of files processed)
    """
    logger = logging.getLogger(f"process_batch_{batch_id}")
    logger.info(f"Processing HHSearch results for batch {batch_id}")

    # Initialize application context
    context = ApplicationContext(config_path)

    # Create registrar
    registrar = HHResultRegistrar(context)

    try:
        # Process batch
        result = registrar.register_batch_results(batch_id, force)

        logger.info(f"Successfully registered {result} HHR files for batch {batch_id}")
        return batch_id, result
    except Exception as e:
        logger.error(f"Error processing batch {batch_id}: {e}")
        logger.exception("Stack trace:")
        return batch_id, -1

def process_multiple_batches(args: argparse.Namespace) -> int:
    """
    Process multiple batches in parallel using the database backend

    Args:
        args: Command line arguments

    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("hhsearch.batches")
    logger.info(f"Processing multiple batches with {args.max_workers} workers")

    # Initialize application context
    context = ApplicationContext(args.config)

    # Get batch IDs to process
    if args.batch_ids:
        # Use specified batch IDs
        batch_ids = [b_id for b_id in args.batch_ids if b_id not in args.exclude_batch_ids]
    else:
        # Get all batch IDs from database
        all_batches = get_all_batch_ids(context)
        batch_ids = [b_id for b_id, _ in all_batches if b_id not in args.exclude_batch_ids]

    logger.info(f"Will process {len(batch_ids)} batches: {batch_ids}")

    # Process all batches using ThreadPoolExecutor
    results = {}
    failed_batches = []

    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        future_to_batch = {
            executor.submit(process_batch, batch_id, args.config, args.force): batch_id
            for batch_id in batch_ids
        }

        for future in as_completed(future_to_batch):
            batch_id = future_to_batch[future]
            try:
                batch_id, result = future.result()
                results[batch_id] = result

                if result < 0:
                    failed_batches.append(batch_id)
                    logger.error(f"Batch {batch_id} failed processing")
                else:
                    logger.info(f"Completed batch {batch_id}: registered {result} files")
            except Exception as e:
                failed_batches.append(batch_id)
                logger.error(f"Batch {batch_id} failed with exception: {e}")

    # Print summary
    total_processed = sum(count for count in results.values() if count > 0)
    logger.info(f"Completed processing all batches")
    logger.info(f"Total files processed: {total_processed}")

    if failed_batches:
        logger.error(f"Failed batches: {failed_batches}")
        return 1

    return 0

def batches_mode(args: argparse.Namespace) -> int:
    """
    Run batches mode

    Args:
        args: Command line arguments

    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("hhsearch.batches")
    logger.info("Running batches mode")

    return process_multiple_batches(args)

def collate_all_batches(args: argparse.Namespace) -> int:
    """
    Collate HHSearch results with BLAST evidence for all batches

    Args:
        args: Command line arguments

    Returns:
        Exit code (0 for success)
    """
    logger = logging.getLogger("hhsearch.collate.all_batches")
    logger.info("Starting collation across all batches")

    # Initialize application context
    context = ApplicationContext(args.config)

    # Get all batch IDs
    query = "SELECT id, batch_name FROM ecod_schema.batch ORDER BY id"
    batches = context.db.execute_query(query)

    if not batches:
        logger.error("No batches found in database")
        return 1

    logger.info(f"Found {len(batches)} batches to process")

    # Filter batches if exclude list is provided
    batch_ids = [batch[0] for batch in batches]
    if args.exclude_batch_ids:
        batch_ids = [b_id for b_id in batch_ids if b_id not in args.exclude_batch_ids]
        logger.info(f"Processing {len(batch_ids)} batches after exclusions")

    # Process each batch
    success_count = 0
    failed_batches = []

    for batch_id in batch_ids:
        batch_name = [b[1] for b in batches if b[0] == batch_id][0]
        logger.info(f"Processing batch {batch_id} ({batch_name})")

        # Create a separate log file for each batch if main log file is provided
        batch_log_file = None
        if args.log_file:
            log_dir = os.path.dirname(args.log_file)
            log_base = os.path.basename(args.log_file)
            batch_log_file = os.path.join(log_dir, f"batch_{batch_id}_{log_base}")

        # Create args for single batch collation
        batch_args = argparse.Namespace(
            config=args.config,
            batch_id=batch_id,
            force=args.force,
            limit=args.limit_per_batch,
            log_file=batch_log_file,
            verbose=args.verbose
        )

        # Process batch
        try:
            result = collate_with_blast(batch_args)
            if result == 0:
                success_count += 1
                logger.info(f"Successfully processed batch {batch_id} ({batch_name})")
            else:
                failed_batches.append(batch_id)
                logger.warning(f"Failed to process batch {batch_id} ({batch_name})")
        except Exception as e:
            failed_batches.append(batch_id)
            logger.error(f"Error processing batch {batch_id} ({batch_name}): {str(e)}")

    # Log final summary
    logger.info(f"Collation complete. Processed {len(batch_ids)} batches.")
    logger.info(f"Successful: {success_count}, Failed: {len(failed_batches)}")

    if failed_batches:
        logger.warning(f"Failed batches: {failed_batches}")
        return 1

    return 0

#
# MAIN FUNCTION AND ARGUMENT PARSING
#

def main():
    """Main entry point"""
    # Create top-level parser
    parser = argparse.ArgumentParser(description='HHSearch Tools - Unified HHSearch processing and analysis toolkit')
    subparsers = parser.add_subparsers(dest="mode", help="Operating mode")

    # Process mode parser
    process_parser = subparsers.add_parser("process", help="Convert HHR files to XML and create domain summaries")
    process_subparsers = process_parser.add_subparsers(dest="backend", help="Backend to use")

    # Database backend subparser
    db_parser = process_subparsers.add_parser("db", help="Use database backend")
    db_parser.add_argument('--config', type=str, default='config/config.yml',
                        help='Path to configuration file')
    db_parser.add_argument('--batch-id', type=int, required=True,
                        help='Batch ID to process')
    db_parser.add_argument('--limit', type=int,
                        help='Limit the number of files to process')
    db_parser.add_argument('--log-file', type=str,
                        help='Log file path')
    db_parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose output')
    db_parser.add_argument('--force', action='store_true',
                        help='Force reprocessing of already processed results')

    # Registration backend subparser (using HHResultRegistrar)
    register_parser = process_subparsers.add_parser("register", help="Register HHR files using HHResultRegistrar")
    register_parser.add_argument('--config', type=str, default='config/config.yml',
                              help='Path to configuration file')
    register_parser.add_argument('--batch-id', type=int, required=True,
                              help='Batch ID to process')
    register_parser.add_argument('--chains', nargs='+',
                              help='Specific chains to process (format: pdbid_chainid)')
    register_parser.add_argument('--log-file', type=str,
                              help='Log file path')
    register_parser.add_argument('-v', '--verbose', action='store_true',
                              help='Enable verbose output')
    register_parser.add_argument('--force', action='store_true',
                              help='Force reprocessing of already processed results')

    # Filesystem backend subparser
    fs_parser = process_subparsers.add_parser("fs", help="Use filesystem backend")
    fs_parser.add_argument('--batch-path', type=str, required=True,
                        help='Path to batch directory')
    fs_parser.add_argument('--ref-version', type=str, default="develop291",
                        help='Reference version')
    fs_parser.add_argument('--limit', type=int,
                        help='Limit the number of files to process')
    fs_parser.add_argument('--log-file', type=str,
                        help='Log file path')
    fs_parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose output')
    fs_parser.add_argument('--force', action='store_true',
                        help='Force reprocessing of already processed results')

    # Collate mode parser
    collate_parser = subparsers.add_parser("collate", help="Collate HHSearch results with BLAST evidence")
    collate_parser.add_argument('--config', type=str, default='config/config.yml',
                          help='Path to configuration file')
    collate_parser.add_argument('--batch-id', type=int, required=True,
                          help='Batch ID to process')
    collate_parser.add_argument('--limit', type=int,
                          help='Limit the number of proteins to process')
    collate_parser.add_argument('--log-file', type=str,
                          help='Log file path')
    collate_parser.add_argument('-v', '--verbose', action='store_true',
                          help='Enable verbose output')
    collate_parser.add_argument('--force', action='store_true',
                          help='Force reprocessing of already processed results')

    # Collate-all mode parser
    collate_all_parser = subparsers.add_parser("collate-all",
        help="Collate HHSearch results with BLAST evidence for all batches")
    collate_all_parser.add_argument('--config', type=str, default='config/config.yml',
                               help='Path to configuration file')
    collate_all_parser.add_argument('--exclude-batch-ids', nargs='+', type=int, default=[],
                               help='Batch IDs to exclude')
    collate_all_parser.add_argument('--limit-per-batch', type=int,
                               help='Limit the number of proteins to process per batch')
    collate_all_parser.add_argument('--log-file', type=str,
                               help='Log file path')
    collate_all_parser.add_argument('-v', '--verbose', action='store_true',
                               help='Enable verbose output')
    collate_all_parser.add_argument('--force', action='store_true',
                               help='Force reprocessing of already processed results')

    # Add parallel-collate mode parser
    parallel_collate_parser = subparsers.add_parser("parallel-collate",
        help="Collate HHSearch results with BLAST evidence in parallel using SLURM")
    parallel_collate_parser.add_argument('--config', type=str, default='config/config.yml',
                                    help='Path to configuration file')
    parallel_collate_parser.add_argument('--batch-ids', nargs='+', type=int, default=None,
                                    help='Specific batch IDs to process')
    parallel_collate_parser.add_argument('--exclude-batch-ids', nargs='+', type=int, default=[],
                                    help='Batch IDs to exclude')
    parallel_collate_parser.add_argument('--force', action='store_true',
                                    help='Force reprocessing of already processed results')
    parallel_collate_parser.add_argument('--threads', type=int,
                                    help='Number of threads per job')
    parallel_collate_parser.add_argument('--memory', type=str,
                                    help='Memory allocation per job (e.g. "8G")')
    parallel_collate_parser.add_argument('--time', type=str,
                                    help='Time limit per job (e.g. "12:00:00")')
    parallel_collate_parser.add_argument('--wait', action='store_true',
                                    help='Wait for jobs to complete')
    parallel_collate_parser.add_argument('--check-interval', type=int, default=60,
                                    help='Check interval when waiting for jobs (seconds)')
    parallel_collate_parser.add_argument('--timeout', type=int,
                                    help='Timeout when waiting for jobs (seconds)')
    parallel_collate_parser.add_argument('--log-file', type=str,
                                    help='Log file path')
    parallel_collate_parser.add_argument('-v', '--verbose', action='store_true',
                                    help='Enable verbose output')

    # Analyze mode parser
    analyze_parser = subparsers.add_parser("analyze", help="Examine and diagnose HHSearch results")
    analyze_subparsers = analyze_parser.add_subparsers(dest="action", help="Analysis action")

    # Content analysis subparser
    content_parser = analyze_subparsers.add_parser("content", help="Check content of HHR and XML files")
    content_parser.add_argument('--batch-path', type=str, required=True,
                           help='Path to batch directory')
    content_parser.add_argument('--ref-version', type=str, default="develop291",
                           help='Reference version')
    content_parser.add_argument('--sample-size', type=int, default=5,
                           help='Number of files to examine')
    content_parser.add_argument('-v', '--verbose', action='store_true',
                           help='Enable verbose output')
    content_parser.add_argument('--log-file', type=str,
                           help='Log file path')

    # Missing files analysis subparser
    missing_parser = analyze_subparsers.add_parser("missing", help="Find missing files")
    missing_parser.add_argument('--batch-path', type=str, required=True,
                            help='Path to batch directory')
    missing_parser.add_argument('--ref-version', type=str, default="develop291",
                            help='Reference version')
    missing_parser.add_argument('-v', '--verbose', action='store_true',
                            help='Enable verbose output')
    missing_parser.add_argument('--log-file', type=str,
                            help='Log file path')

    # Repair mode parser
    repair_parser = subparsers.add_parser("repair", help="Fix missing or problematic files")
    repair_subparsers = repair_parser.add_subparsers(dest="action", help="Repair action")

    # Repair missing files subparser
    missing_repair_parser = repair_subparsers.add_parser("missing", help="Repair missing files")
    missing_repair_parser.add_argument('--batch-path', type=str, required=True,
                                  help='Path to batch directory')
    missing_repair_parser.add_argument('--ref-version', type=str, default="develop291",
                                  help='Reference version')
    missing_repair_parser.add_argument('-v', '--verbose', action='store_true',
                                  help='Enable verbose output')
    missing_repair_parser.add_argument('--log-file', type=str,
                                  help='Log file path')

    # Add fix_summaries subparser
    fix_summaries_parser = repair_subparsers.add_parser("fix_summaries", help="Fix domain summaries to include HHSearch hits")
    fix_summaries_parser.add_argument('--batch-path', type=str, required=True,
                                 help='Path to batch directory')
    fix_summaries_parser.add_argument('--ref-version', type=str, default="develop291",
                                 help='Reference version')
    fix_summaries_parser.add_argument('--limit', type=int,
                                 help='Limit the number of files to process')
    fix_summaries_parser.add_argument('-v', '--verbose', action='store_true',
                                 help='Enable verbose output')
    fix_summaries_parser.add_argument('--log-file', type=str,
                                 help='Log file path')

    # Fix database paths subparser
    db_paths_parser = repair_subparsers.add_parser("db_paths", help="Fix file paths in database")
    db_paths_parser.add_argument('--config', type=str, default='config/config.yml',
                             help='Path to configuration file')
    db_paths_parser.add_argument('--batch-id', type=int, required=True,
                             help='Batch ID to process')
    db_paths_parser.add_argument('--dry-run', action='store_true',
                             help='Check but don\'t make changes')
    db_paths_parser.add_argument('-v', '--verbose', action='store_true',
                             help='Enable verbose output')
    db_paths_parser.add_argument('--log-file', type=str,
                             help='Log file path')

    # Create empty summaries subparser
    empty_parser = repair_subparsers.add_parser("empty_summaries", help="Create empty summary files for no-hits")
    empty_parser.add_argument('--config', type=str, default='config/config.yml',
                         help='Path to configuration file')
    empty_parser.add_argument('--batch-id', type=int, required=True,
                         help='Batch ID to process')
    empty_parser.add_argument('--dry-run', action='store_true',
                         help='Check but don\'t make changes')
    empty_parser.add_argument('-v', '--verbose', action='store_true',
                         help='Enable verbose output')
    empty_parser.add_argument('--log-file', type=str,
                         help='Log file path')

    # Batches mode parser
    batches_parser = subparsers.add_parser("batches", help="Process multiple batches in parallel")
    batches_parser.add_argument('--config', type=str, default='config/config.yml',
                            help='Path to configuration file')
    batches_parser.add_argument('--batch-ids', nargs='+', type=int, default=None,
                            help='Specific batch IDs to process')
    batches_parser.add_argument('--exclude-batch-ids', nargs='+', type=int, default=[],
                            help='Batch IDs to exclude')
    batches_parser.add_argument('--force', action='store_true',
                            help='Force reprocessing of already processed results')
    batches_parser.add_argument('--max-workers', type=int, default=4,
                            help='Maximum number of worker threads')
    batches_parser.add_argument('--log-file', type=str,
                            help='Log file path')
    batches_parser.add_argument('-v', '--verbose', action='store_true',
                            help='Enable verbose output')

    # Parse arguments
    args = parser.parse_args()

    # Set up logging
    log_file = args.log_file if hasattr(args, 'log_file') and args.log_file else None
    verbose = args.verbose if hasattr(args, 'verbose') else False
    setup_logging(verbose, log_file)

    logger = logging.getLogger("hhsearch_tools")

    # Run appropriate mode
    if args.mode == "process":
        return process_mode(args)
    elif args.mode == "collate":
        return collate_mode(args)
    elif args.mode == "collate-all":
        return collate_all_batches(args)
    elif args.mode == "analyze":
        return analyze_mode(args)
    elif args.mode == "repair":
        return repair_mode(args)
    elif args.mode == "batches":
        return batches_mode(args)
    elif args.mode == "parallel-collate":
        return parallel_collate_batches(args)
    else:
        logger.error(f"Unknown mode: {args.mode}")
        parser.print_help()
        return 1

if __name__ == "__main__":
    sys.exit(main())
