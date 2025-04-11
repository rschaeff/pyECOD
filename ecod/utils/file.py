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
from contextlib import contextmanager
from pathlib import Path
from typing import Optional, Any, Generator, BinaryIO, TextIO, Union, Dict, List

from ecod.exceptions import FileOperationError

logger = logging.getLogger("ecod.utils.file")


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