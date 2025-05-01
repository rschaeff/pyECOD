#!/usr/bin/env python3
"""
File system utilities for the ECOD pipeline.
Provides safe file operations with error handling.
"""
import os
import shutil
import logging
import tempfile
import hashlib
import re
from contextlib import contextmanager
from pathlib import Path
from typing import Optional, Any, Generator, BinaryIO, TextIO, Union, Dict, List, Tuple

from ecod.exceptions import FileOperationError

logger = logging.getLogger("ecod.utils.file")

def find_fasta_file(pdb_id: str, chain_id: str, base_path: str, db_manager=None) -> str:
    """Find FASTA file for a protein chain

    Args:
        pdb_id: PDB ID
        chain_id: Chain ID
        base_path: Base directory for job files
        db_manager: Optional database manager for db lookup

    Returns:
        Path to FASTA file or empty string if not found
    """
    pdb_chain = f"{pdb_id}_{chain_id}"
    logger.debug(f"Looking for FASTA file for {pdb_chain}")

    # Try standard locations for FASTA file - ordered by likelihood
    potential_paths = [
        os.path.join(base_path, "fastas", f"{pdb_chain}.fa"),
        os.path.join(base_path, "fastas", f"{pdb_chain}.fasta"),
        os.path.join(base_path, "fastas", "batch_0", f"{pdb_chain}.fa"),
        os.path.join(base_path, "fastas", "batch_0", f"{pdb_chain}.fasta"),
        os.path.join(base_path, "fastas", "batch_1", f"{pdb_chain}.fa"),
        os.path.join(base_path, "fastas", "batch_1", f"{pdb_chain}.fasta"),
        os.path.join(base_path, "fastas", "batch_2", f"{pdb_chain}.fa"),
        os.path.join(base_path, "fastas", "batch_2", f"{pdb_chain}.fasta"),
        os.path.join(base_path, pdb_chain, f"{pdb_chain}.fa"),
        os.path.join(base_path, pdb_chain, f"{pdb_chain}.fasta"),
        os.path.join(base_path, "ecod_dump", pdb_chain, f"{pdb_chain}.fa"),
        os.path.join(base_path, "ecod_dump", pdb_chain, f"{pdb_chain}.fasta")
    ]

    # First check standard locations
    for path in potential_paths:
        if os.path.exists(path) and os.path.getsize(path) > 0:
            logger.info(f"Found FASTA file at: {path}")
            return path

    # If not found and db_manager is provided, query database as last resort
    if db_manager:
        try:
            query = """
            SELECT pf.file_path
            FROM ecod_schema.process_file pf
            JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
            JOIN ecod_schema.protein p ON ps.protein_id = p.id
            WHERE p.pdb_id = %s AND p.chain_id = %s
            AND pf.file_type = 'fasta'
            AND pf.file_exists = TRUE
            LIMIT 1
            """
            rows = db_manager.execute_query(query, (pdb_id, chain_id))
            if rows:
                db_fasta_path = rows[0][0]
                # Handle relative paths
                if not os.path.isabs(db_fasta_path):
                    full_fasta_path = os.path.join(base_path, db_fasta_path)
                else:
                    full_fasta_path = db_fasta_path

                if os.path.exists(full_fasta_path) and os.path.getsize(full_fasta_path) > 0:
                    logger.info(f"Found FASTA file via database: {full_fasta_path}")
                    return full_fasta_path
        except Exception as e:
            logger.warning(f"Error querying database for FASTA file: {e}")

    logger.warning(f"No FASTA file found for {pdb_id}_{chain_id}")
    return ""

def find_domain_summary(pdb_id: str, chain_id: str, base_path: str,
                       reference: str, blast_only: bool = False,
                       db_manager=None) -> str:
    """Find domain summary file

    Args:
        pdb_id: PDB ID
        chain_id: Chain ID
        base_path: Base directory for job files
        reference: Reference version
        blast_only: Whether to look for blast-only summary
        db_manager: Optional database manager for db lookup

    Returns:
        Path to domain summary file or empty string if not found
    """
    pdb_chain = f"{pdb_id}_{chain_id}"
    logger.debug(f"Looking for domain summary for {pdb_chain}")

    # Determine file suffix
    suffix = ".blast_only" if blast_only else ""

    # Try standard locations for domain summary file
    standard_path = os.path.join(base_path, "domains", f"{pdb_chain}.{reference}.domain_summary{suffix}.xml")
    if os.path.exists(standard_path) and os.path.getsize(standard_path) > 0:
        logger.info(f"Found domain summary at standard path: {standard_path}")
        return standard_path

    # Try legacy locations
    legacy_paths = [
        os.path.join(base_path, "domains", f"{pdb_chain}.{reference}.domain_summ{suffix}.xml"),
        os.path.join(base_path, pdb_chain, f"{pdb_chain}.{reference}.domain_summary{suffix}.xml"),
        os.path.join(base_path, "ecod_dump", pdb_chain, f"{pdb_chain}.{reference}.domain_summary{suffix}.xml")
    ]

    for path in legacy_paths:
        if os.path.exists(path) and os.path.getsize(path) > 0:
            logger.info(f"Found domain summary at legacy path: {path}")
            return path

    # If not found and db_manager is provided, query database as last resort
    if db_manager:
        try:
            file_type = "blast_only_summary" if blast_only else "domain_summary"
            query = """
            SELECT pf.file_path
            FROM ecod_schema.process_file pf
            JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
            JOIN ecod_schema.protein p ON ps.protein_id = p.id
            WHERE p.pdb_id = %s AND p.chain_id = %s
            AND pf.file_type = %s
            AND pf.file_exists = TRUE
            LIMIT 1
            """
            rows = db_manager.execute_query(query, (pdb_id, chain_id, file_type))
            if rows:
                db_path = rows[0][0]
                # Handle relative paths
                if not os.path.isabs(db_path):
                    full_path = os.path.join(base_path, db_path)
                else:
                    full_path = db_path

                if os.path.exists(full_path) and os.path.getsize(full_path) > 0:
                    logger.info(f"Found domain summary via database: {full_path}")
                    return full_path
        except Exception as e:
            logger.warning(f"Error querying database for domain summary: {e}")

    logger.warning(f"No domain summary found for {pdb_id}_{chain_id}")
    return ""

def find_blast_results(pdb_id: str, chain_id: str, base_path: str,
                      reference: str, blast_type: str,
                      db_manager=None) -> str:
    """Find BLAST result file

    Args:
        pdb_id: PDB ID
        chain_id: Chain ID
        base_path: Base directory
        reference: Reference version
        blast_type: Type of BLAST results ('chain' or 'domain')
        db_manager: Optional database manager for db lookup

    Returns:
        Path to BLAST result file or empty string if not found
    """
    pdb_chain = f"{pdb_id}_{chain_id}"
    logger.debug(f"Looking for {blast_type} BLAST results for {pdb_chain}")

    # Verify blast_type is valid
    if blast_type not in ['chain', 'domain']:
        logger.error(f"Invalid blast_type: {blast_type} (must be 'chain' or 'domain')")
        return ""

    # Define standard and legacy patterns based on blast_type
    if blast_type == 'chain':
        standard_path = os.path.join(
            base_path, "blast", "chain", f"{pdb_chain}.{reference}.xml"
        )
        legacy_paths = [
            os.path.join(base_path, "blast", f"{pdb_chain}.{reference}.chainwise_blast.xml"),
            os.path.join(base_path, "blast", "chain", f"{pdb_chain}.chainwise.{reference}.xml"),
            os.path.join(base_path, "blast", f"{pdb_chain}.{reference}.chain_blast.xml"),
            os.path.join(base_path, "blast", "chain", "batch_0", f"{pdb_chain}.chainwise_blast.xml"),
            os.path.join(base_path, "blast", "chain", "batch_1", f"{pdb_chain}.chainwise_blast.xml"),
            os.path.join(base_path, "blast", "chain", "batch_2", f"{pdb_chain}.chainwise_blast.xml")
        ]
        db_file_type = "chain_blast_result"
    else:  # domain
        standard_path = os.path.join(
            base_path, "blast", "domain", f"{pdb_chain}.{reference}.xml"
        )
        legacy_paths = [
            os.path.join(base_path, "blast", f"{pdb_chain}.{reference}.blast.xml"),
            os.path.join(base_path, "blast", "domain", f"{pdb_chain}.{reference}.blast.xml"),
            os.path.join(base_path, "blast", f"{pdb_chain}.{reference}.domain_blast.xml"),
            os.path.join(base_path, "blast", "domain", "batch_0", f"{pdb_chain}.domain_blast.xml"),
            os.path.join(base_path, "blast", "domain", "batch_1", f"{pdb_chain}.domain_blast.xml"),
            os.path.join(base_path, "blast", "domain", "batch_2", f"{pdb_chain}.domain_blast.xml")
        ]
        db_file_type = "domain_blast_result"

    # Check standard path first
    if os.path.exists(standard_path) and os.path.getsize(standard_path) > 0:
        logger.info(f"Found {blast_type} BLAST results at standard path: {standard_path}")
        return standard_path

    # Check legacy paths
    for path in legacy_paths:
        if os.path.exists(path) and os.path.getsize(path) > 0:
            logger.info(f"Found {blast_type} BLAST results at legacy path: {path}")
            return path

    # If not found and db_manager is provided, query database as last resort
    if db_manager:
        try:
            query = """
            SELECT pf.file_path
            FROM ecod_schema.process_file pf
            JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
            JOIN ecod_schema.protein p ON ps.protein_id = p.id
            WHERE p.pdb_id = %s AND p.chain_id = %s
            AND pf.file_type = %s
            AND pf.file_exists = TRUE
            LIMIT 1
            """
            rows = db_manager.execute_query(query, (pdb_id, chain_id, db_file_type))
            if rows:
                db_path = rows[0][0]
                # Handle relative paths
                if not os.path.isabs(db_path):
                    full_path = os.path.join(base_path, db_path)
                else:
                    full_path = db_path

                if os.path.exists(full_path) and os.path.getsize(full_path) > 0:
                    logger.info(f"Found {blast_type} BLAST results via database: {full_path}")
                    return full_path
        except Exception as e:
            logger.warning(f"Error querying database for {blast_type} BLAST results: {e}")

    logger.warning(f"No {blast_type} BLAST results found for {pdb_id}_{chain_id}")
    return ""

def find_hhsearch_file(pdb_id: str, chain_id: str, base_path: str,
                      reference: str, db_manager=None) -> str:
    """Find HHSearch result file

    Args:
        pdb_id: PDB ID
        chain_id: Chain ID
        base_path: Base directory
        reference: Reference version
        db_manager: Optional database manager for db lookup

    Returns:
        Path to HHSearch result file or empty string if not found
    """
    pdb_chain = f"{pdb_id}_{chain_id}"
    logger.debug(f"Looking for HHSearch results for {pdb_chain}")

    # Define standard path (preferred format)
    standard_path = os.path.join(
        base_path, "hhsearch", f"{pdb_chain}.{reference}.hhsearch.xml"
    )

    # Check standard path first
    if os.path.exists(standard_path) and os.path.getsize(standard_path) > 0:
        logger.info(f"Found HHSearch results at standard path: {standard_path}")
        return standard_path

    # Define legacy paths
    legacy_paths = [
        os.path.join(base_path, "hhsearch", f"{pdb_chain}.{reference}.hhr"),
        os.path.join(base_path, "hhsearch", f"{pdb_chain}.hhsearch.{reference}.xml"),
        os.path.join(base_path, "hhsearch", f"{pdb_chain}.{reference}.hh_summ.xml"),
        os.path.join(base_path, pdb_chain, f"{pdb_chain}.{reference}.hhsearch.xml"),
        os.path.join(base_path, pdb_chain, f"{pdb_chain}.{reference}.hhr"),
        os.path.join(base_path, "ecod_dump", pdb_chain, f"{pdb_chain}.{reference}.hhsearch.xml"),
        os.path.join(base_path, "ecod_dump", pdb_chain, f"{pdb_chain}.{reference}.hhr")
    ]

    # Check legacy paths
    for path in legacy_paths:
        if os.path.exists(path) and os.path.getsize(path) > 0:
            logger.info(f"Found HHSearch results at legacy path: {path}")
            return path

    # If not found and db_manager is provided, query database as last resort
    if db_manager:
        try:
            # Check for various possible file types in DB
            for file_type in ["hhsearch_result", "hhr", "hhsearch_xml"]:
                query = """
                SELECT pf.file_path
                FROM ecod_schema.process_file pf
                JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
                JOIN ecod_schema.protein p ON ps.protein_id = p.id
                WHERE p.pdb_id = %s AND p.chain_id = %s
                AND pf.file_type = %s
                AND pf.file_exists = TRUE
                LIMIT 1
                """
                rows = db_manager.execute_query(query, (pdb_id, chain_id, file_type))
                if rows:
                    db_path = rows[0][0]
                    # Handle relative paths
                    if not os.path.isabs(db_path):
                        full_path = os.path.join(base_path, db_path)
                    else:
                        full_path = db_path

                    if os.path.exists(full_path) and os.path.getsize(full_path) > 0:
                        logger.info(f"Found HHSearch results via database: {full_path} ({file_type})")
                        return full_path
        except Exception as e:
            logger.warning(f"Error querying database for HHSearch results: {e}")

    logger.warning(f"No HHSearch results found for {pdb_id}_{chain_id}")
    return ""

def find_self_comparison(pdb_id: str, chain_id: str, base_path: str,
                        db_manager=None) -> str:
    """Find self-comparison file for a protein chain

    Args:
        pdb_id: PDB ID
        chain_id: Chain ID
        base_path: Base directory for job files
        db_manager: Optional database manager for db lookup

    Returns:
        Path to self-comparison file or empty string if not found
    """
    pdb_chain = f"{pdb_id}_{chain_id}"
    logger.debug(f"Looking for self-comparison for {pdb_chain}")

    # Check standard location
    standard_path = os.path.join(base_path, "self_comparisons", f"{pdb_chain}.self_comp.xml")
    if os.path.exists(standard_path) and os.path.getsize(standard_path) > 0:
        logger.info(f"Found self-comparison at standard location: {standard_path}")
        return standard_path

    # Check legacy locations
    legacy_paths = [
        os.path.join(base_path, pdb_chain, f"{pdb_chain}.self_comp.xml"),
        os.path.join(base_path, "ecod_dump", pdb_chain, f"{pdb_chain}.self_comp.xml")
    ]

    for path in legacy_paths:
        if os.path.exists(path) and os.path.getsize(path) > 0:
            logger.info(f"Found self-comparison at legacy location: {path}")
            return path

    # If not found and db_manager is provided, query database as last resort
    if db_manager:
        try:
            query = """
            SELECT pf.file_path
            FROM ecod_schema.process_file pf
            JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
            JOIN ecod_schema.protein p ON ps.protein_id = p.id
            WHERE p.pdb_id = %s AND p.chain_id = %s
            AND pf.file_type = 'self_comparison'
            AND pf.file_exists = TRUE
            LIMIT 1
            """
            rows = db_manager.execute_query(query, (pdb_id, chain_id))
            if rows:
                db_path = rows[0][0]
                # Handle relative paths
                if not os.path.isabs(db_path):
                    full_path = os.path.join(base_path, db_path)
                else:
                    full_path = db_path

                if os.path.exists(full_path) and os.path.getsize(full_path) > 0:
                    logger.info(f"Found self-comparison via database: {full_path}")
                    return full_path
        except Exception as e:
            logger.warning(f"Error querying database for self-comparison: {e}")

    logger.warning(f"No self-comparison found for {pdb_id}_{chain_id}")
    return ""

def read_sequence_from_fasta(fasta_path: str) -> Tuple[str, str]:
    """Read sequence from a FASTA file

    Args:
        fasta_path: Path to FASTA file

    Returns:
        Tuple of (header, sequence)

    Raises:
        FileOperationError: If file cannot be read
    """
    if not os.path.exists(fasta_path):
        error_msg = f"FASTA file does not exist: {fasta_path}"
        logger.error(error_msg)
        raise FileOperationError(error_msg)

    try:
        with safe_open(fasta_path, 'r') as f:
            lines = f.readlines()

        if not lines:
            error_msg = f"Empty FASTA file: {fasta_path}"
            logger.error(error_msg)
            raise FileOperationError(error_msg)

        # Get header (first line)
        header = lines[0].strip()
        if not header.startswith('>'):
            error_msg = f"Invalid FASTA format (header doesn't start with '>'): {fasta_path}"
            logger.error(error_msg)
            raise FileOperationError(error_msg)

        header = header[1:]  # Remove '>' character

        # Get sequence (remaining lines)
        sequence = ""
        for line in lines[1:]:
            sequence += line.strip()

        return header, sequence
    except (FileOperationError, UnicodeDecodeError) as e:
        # Re-raise FileOperationError, wrap other exceptions
        if isinstance(e, FileOperationError):
            raise
        error_msg = f"Error reading FASTA file {fasta_path}: {str(e)}"
        logger.error(error_msg)
        raise FileOperationError(error_msg) from e

def extract_pdb_chain_from_path(file_path: str) -> Tuple[Optional[str], Optional[str]]:
    """Extract PDB ID and chain ID from file path

    Args:
        file_path: Path to file

    Returns:
        Tuple of (pdb_id, chain_id) or (None, None) if not found
    """
    basename = os.path.basename(file_path)

    # Try different patterns
    # Pattern 1: pdbid_chainid.* (e.g., 1abc_A.fa, 1abc_A.develop291.hhr)
    match = re.match(r'([a-zA-Z0-9]{4})_([a-zA-Z0-9])\.', basename)
    if match:
        return match.group(1), match.group(2)

    # Pattern 2: pdbid_chainid_* (e.g., 1abc_A_domains.xml)
    match = re.match(r'([a-zA-Z0-9]{4})_([a-zA-Z0-9])_', basename)
    if match:
        return match.group(1), match.group(2)

    # Pattern 3: Directory name might be pdbid_chainid
    dir_name = os.path.basename(os.path.dirname(file_path))
    match = re.match(r'([a-zA-Z0-9]{4})_([a-zA-Z0-9])$', dir_name)
    if match:
        return match.group(1), match.group(2)

    return None, None

def check_input_files(paths: Dict[str, Dict[str, str]], required: List[str] = None,
                    optional: List[str] = None) -> Dict[str, bool]:
    """Check existence of input files with consolidated approach

    Args:
        paths: Dictionary of file paths from get_all_evidence_paths
        required: List of required file types
        optional: List of optional file types

    Returns:
        Dictionary of {file_type: exists} for all checked files
    """
    result = {}

    # Check required files
    for file_type in (required or []):
        file_path = paths.get(file_type, {}).get('exists_at')
        exists = file_path is not None and os.path.exists(file_path)
        result[file_type] = exists

        if not exists:
            logger.warning(f"Required file {file_type} not found")

    # Check optional files
    for file_type in (optional or []):
        file_path = paths.get(file_type, {}).get('exists_at')
        exists = file_path is not None and os.path.exists(file_path)
        result[file_type] = exists

    return result

def ensure_dir(directory: str) -> bool:
    """Ensure a directory exists, creating it if necessary

    Args:
        directory: Directory path

    Returns:
        True if successful

    Raises:
        FileOperationError: If directory cannot be created
    """
    try:
        if not os.path.exists(directory):
            logger.debug(f"Creating directory: {directory}")
            os.makedirs(directory, exist_ok=True)
        return True
    except OSError as e:
        error_msg = f"Error creating directory {directory}: {str(e)}"
        logger.error(error_msg)
        raise FileOperationError(error_msg, {"directory": directory}) from e


@contextmanager
def safe_open(file_path: str, mode: str = 'r', encoding: Optional[str] = None) -> Generator[Any, None, None]:
    """Safely open a file with error handling

    Args:
        file_path: Path to the file
        mode: File open mode
        encoding: File encoding

    Yields:
        Open file object

    Raises:
        FileOperationError: If file cannot be opened
    """
    try:
        # Ensure directory exists for write operations
        if ('w' in mode or 'a' in mode or '+' in mode) and not os.path.exists(os.path.dirname(file_path)):
            os.makedirs(os.path.dirname(file_path), exist_ok=True)
            logger.debug(f"Created directory: {os.path.dirname(file_path)}")

        # Open the file
        logger.debug(f"Opening file: {file_path} (mode: {mode})")
        if encoding and 'b' not in mode:
            file = open(file_path, mode, encoding=encoding)
        else:
            file = open(file_path, mode)

        try:
            yield file
        finally:
            file.close()
            logger.debug(f"Closed file: {file_path}")

    except (OSError, IOError) as e:
        error_msg = f"Error accessing file {file_path}: {str(e)}"
        logger.error(error_msg)
        raise FileOperationError(error_msg, {"file_path": file_path, "mode": mode}) from e


@contextmanager
def atomic_write(file_path: str, mode: str = 'w',
                encoding: Optional[str] = None) -> Generator[Union[TextIO, BinaryIO], None, None]:
    """Write to a file atomically using a temporary file

    Writes to a temporary file first, then renames it to the target file
    to ensure the operation is atomic and prevent partial writes.

    Args:
        file_path: Path to the file
        mode: File open mode (must be a write mode)
        encoding: File encoding

    Yields:
        Open temporary file object

    Raises:
        FileOperationError: If file operation fails
        ValueError: If mode is not a write mode
    """
    if 'w' not in mode and 'a' not in mode and '+' not in mode:
        error_msg = f"Invalid mode for atomic_write: {mode} (must be write mode)"
        logger.error(error_msg)
        raise ValueError(error_msg)

    # Create a temporary file in the same directory
    base_dir = os.path.dirname(file_path) or '.'

    try:
        ensure_dir(base_dir)

        temp_suffix = f".{os.path.basename(file_path)}.tmp"
        with tempfile.NamedTemporaryFile(mode='wb', suffix=temp_suffix,
                                       dir=base_dir, delete=False) as temp_file:
            temp_path = temp_file.name
            logger.debug(f"Created temporary file: {temp_path} for atomic write to {file_path}")

        # Reopen file with correct mode and encoding
        if 'b' in mode:
            final_file = open(temp_path, mode)
        else:
            final_file = open(temp_path, mode, encoding=encoding or 'utf-8')

        try:
            yield final_file

            # Close the file before renaming
            final_file.close()

            # Rename the temporary file to the target file
            shutil.move(temp_path, file_path)
            logger.debug(f"Atomically wrote to file: {file_path}")

        except Exception as e:
            # Make sure the file is closed
            if not final_file.closed:
                final_file.close()

            # Clean up the temporary file on error
            if os.path.exists(temp_path):
                os.unlink(temp_path)
                logger.debug(f"Deleted temporary file: {temp_path} after error")
            raise e

    except (OSError, IOError) as e:
        error_msg = f"Error during atomic write to {file_path}: {str(e)}"
        logger.error(error_msg)
        raise FileOperationError(error_msg, {"file_path": file_path}) from e

def verify_file_path(file_path: str, normalize: bool = True) -> str:
    """Verify and normalize a file path

    Args:
        file_path: Path to verify
        normalize: Whether to normalize the path (resolve '..')

    Returns:
        Normalized path if it exists, original path otherwise
    """
    if not file_path:
        return file_path

    if normalize:
        normalized_path = os.path.normpath(file_path)
        if os.path.exists(normalized_path):
            logger.debug(f"Normalized path {file_path} to {normalized_path}")
            return normalized_path

    if os.path.exists(file_path):
        return file_path

    logger.warning(f"File path does not exist: {file_path}")
    return file_path

def check_file_exists(file_path: str, min_size: int = 0,
                    expect_error: bool = False) -> bool:
    """Check if a file exists and has a minimum size

    Args:
        file_path: Path to the file
        min_size: Minimum file size in bytes
        expect_error: If True, don't log warning for missing file

    Returns:
        True if file exists and meets size requirement
    """
    try:
        if not os.path.exists(file_path):
            if not expect_error:
                logger.debug(f"File does not exist: {file_path}")
            return False

        if min_size > 0:
            size = os.path.getsize(file_path)
            if size < min_size:
                logger.debug(f"File too small: {file_path} ({size} bytes < {min_size} bytes)")
                return False

        return True
    except OSError as e:
        logger.warning(f"Error checking file {file_path}: {str(e)}")
        return False


def calculate_md5(file_path: str) -> str:
    """Calculate MD5 hash of a file

    Args:
        file_path: Path to the file

    Returns:
        MD5 hash as hexadecimal string

    Raises:
        FileOperationError: If file cannot be accessed
    """
    try:
        hash_md5 = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()
    except (OSError, IOError) as e:
        error_msg = f"Error calculating MD5 for file {file_path}: {str(e)}"
        logger.error(error_msg)
        raise FileOperationError(error_msg, {"file_path": file_path}) from e


def calculate_md5_string(data: str) -> str:
    """Calculate MD5 hash of a string

    Args:
        data: Input string

    Returns:
        MD5 hash as hexadecimal string
    """
    return hashlib.md5(data.encode('utf-8')).hexdigest()


def get_file_info(file_path: str) -> Dict[str, Any]:
    """Get file information

    Args:
        file_path: Path to the file

    Returns:
        Dictionary with file information
    """
    result = {
        'exists': False,
        'size': None,
        'is_dir': False,
        'modified': None,
        'created': None,
        'md5': None
    }

    try:
        if not os.path.exists(file_path):
            return result

        stat = os.stat(file_path)
        result['exists'] = True
        result['size'] = stat.st_size
        result['is_dir'] = os.path.isdir(file_path)
        result['modified'] = stat.st_mtime
        result['created'] = stat.st_ctime

        # Calculate MD5 for files smaller than 100MB
        if not result['is_dir'] and stat.st_size < 100 * 1024 * 1024:
            result['md5'] = calculate_md5(file_path)

        return result
    except OSError as e:
        logger.warning(f"Error getting file info for {file_path}: {str(e)}")
        return result


def read_text_file(file_path: str, encoding: str = 'utf-8') -> str:
    """Read text file with error handling

    Args:
        file_path: Path to the file
        encoding: File encoding

    Returns:
        File contents as string

    Raises:
        FileOperationError: If file cannot be read
    """
    try:
        with safe_open(file_path, 'r', encoding=encoding) as f:
            return f.read()
    except Exception as e:
        error_msg = f"Error reading file {file_path}: {str(e)}"
        logger.error(error_msg)
        raise FileOperationError(error_msg, {"file_path": file_path}) from e


def write_text_file(file_path: str, content: str, encoding: str = 'utf-8') -> None:
    """Write text file with error handling

    Args:
        file_path: Path to the file
        content: Content to write
        encoding: File encoding

    Raises:
        FileOperationError: If file cannot be written
    """
    try:
        with atomic_write(file_path, 'w', encoding=encoding) as f:
            f.write(content)
    except Exception as e:
        error_msg = f"Error writing file {file_path}: {str(e)}"
        logger.error(error_msg)
        raise FileOperationError(error_msg, {"file_path": file_path}) from e


def read_fasta(file_path: str) -> Dict[str, str]:
    """Read FASTA file

    Args:
        file_path: Path to the FASTA file

    Returns:
        Dictionary with sequence headers as keys and sequences as values

    Raises:
        FileOperationError: If file cannot be read
    """
    sequences = {}
    current_header = None
    current_sequence = []

    try:
        with safe_open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith('>'):
                    # Save previous sequence
                    if current_header is not None:
                        sequences[current_header] = ''.join(current_sequence)

                    # Start new sequence
                    current_header = line[1:]
                    current_sequence = []
                else:
                    current_sequence.append(line)

        # Save last sequence
        if current_header is not None:
            sequences[current_header] = ''.join(current_sequence)

        return sequences
    except Exception as e:
        error_msg = f"Error reading FASTA file {file_path}: {str(e)}"
        logger.error(error_msg)
        raise FileOperationError(error_msg, {"file_path": file_path}) from e


def write_fasta(file_path: str, sequences: Dict[str, str], line_width: int = 60) -> None:
    """Write FASTA file

    Args:
        file_path: Path to the FASTA file
        sequences: Dictionary with sequence headers as keys and sequences as values
        line_width: Width of sequence lines

    Raises:
        FileOperationError: If file cannot be written
    """
    try:
        with atomic_write(file_path, 'w') as f:
            for header, sequence in sequences.items():
                f.write(f">{header}\n")

                # Write sequence in chunks of line_width characters
                for i in range(0, len(sequence), line_width):
                    f.write(f"{sequence[i:i+line_width]}\n")
    except Exception as e:
        error_msg = f"Error writing FASTA file {file_path}: {str(e)}"
        logger.error(error_msg)
        raise FileOperationError(error_msg, {"file_path": file_path}) from e


def clean_directory(directory: str, pattern: Optional[str] = None,
                  recursive: bool = False) -> int:
    """Clean a directory by removing files (optionally matching a pattern)

    Args:
        directory: Directory path
        pattern: File pattern to match (glob style)
        recursive: Whether to clean subdirectories recursively

    Returns:
        Number of files removed

    Raises:
        FileOperationError: If directory cannot be cleaned
    """
    try:
        if not os.path.exists(directory):
            return 0

        count = 0
        path = Path(directory)

        if pattern:
            # Remove files matching the pattern
            for file_path in path.glob(pattern):
                if file_path.is_file():
                    file_path.unlink()
                    count += 1
        else:
            # Remove all files (and subdirectories if recursive)
            for item in path.iterdir():
                if item.is_file():
                    item.unlink()
                    count += 1
                elif item.is_dir() and recursive:
                    shutil.rmtree(item)
                    count += 1

        return count
    except Exception as e:
        error_msg = f"Error cleaning directory {directory}: {str(e)}"
        logger.error(error_msg)
        raise FileOperationError(error_msg, {"directory": directory}) from e


def list_files(directory: str, pattern: Optional[str] = None,
              recursive: bool = False) -> List[str]:
    """List files in a directory

    Args:
        directory: Directory path
        pattern: File pattern to match (glob style)
        recursive: Whether to search subdirectories recursively

    Returns:
        List of file paths
    """
    try:
        if not os.path.exists(directory):
            return []

        path = Path(directory)
        if pattern:
            if recursive:
                return [str(p) for p in path.glob(f"**/{pattern}") if p.is_file()]
            else:
                return [str(p) for p in path.glob(pattern) if p.is_file()]
        else:
            if recursive:
                return [str(p) for p in path.glob("**/*") if p.is_file()]
            else:
                return [str(p) for p in path.iterdir() if p.is_file()]
    except Exception as e:
        logger.warning(f"Error listing files in {directory}: {str(e)}")
        return []


def copy_file(source: str, destination: str, overwrite: bool = True) -> bool:
    """Copy a file with error handling

    Args:
        source: Source file path
        destination: Destination file path
        overwrite: Whether to overwrite existing files

    Returns:
        True if copy was successful

    Raises:
        FileOperationError: If file copy fails
    """
    try:
        if not os.path.exists(source):
            raise FileOperationError(f"Source file does not exist: {source}")

        if os.path.exists(destination) and not overwrite:
            logger.info(f"Destination file exists, not overwriting: {destination}")
            return False

        # Ensure destination directory exists
        ensure_dir(os.path.dirname(destination))

        # Copy file
        shutil.copy2(source, destination)
        logger.debug(f"Copied file: {source} -> {destination}")
        return True
    except FileOperationError:
        raise
    except Exception as e:
        error_msg = f"Error copying file {source} to {destination}: {str(e)}"
        logger.error(error_msg)
        raise FileOperationError(error_msg, {"source": source, "destination": destination}) from e
