# ecod/core/file_utils.py
import os
import logging
import shutil
import tempfile
from contextlib import contextmanager
from typing import Optional, Any, Generator
from .exceptions import FileOperationError

logger = logging.getLogger("ecod.file_utils")

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
        if 'w' in mode and not os.path.exists(os.path.dirname(file_path)):
            os.makedirs(os.path.dirname(file_path), exist_ok=True)
            logger.debug(f"Created directory: {os.path.dirname(file_path)}")
            
        # Open the file
        logger.debug(f"Opening file: {file_path} (mode: {mode})")
        if encoding:
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
        raise FileOperationError(error_msg) from e

@contextmanager
def atomic_write(file_path: str, mode: str = 'w', encoding: Optional[str] = None) -> Generator[Any, None, None]:
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
    """
    if 'w' not in mode and 'a' not in mode and '+' not in mode:
        error_msg = f"Invalid mode for atomic_write: {mode} (must be write mode)"
        logger.error(error_msg)
        raise ValueError(error_msg)
        
    # Create a temporary file in the same directory
    base_dir = os.path.dirname(file_path) or '.'
    try:
        os.makedirs(base_dir, exist_ok=True)
        temp_file = tempfile.NamedTemporaryFile(mode=mode.replace('b', '') + 'b',
                                              dir=base_dir,
                                              delete=False)
        temp_path = temp_file.name
        logger.debug(f"Created temporary file: {temp_path} for atomic write to {file_path}")
        
        try:
            if 'b' not in mode and encoding:
                # For text mode with encoding, we need to close and reopen
                temp_file.close()
                temp_file = open(temp_path, mode, encoding=encoding)
                
            yield temp_file
            
            # Close the file before renaming
            temp_file.close()
            
            # Rename the temporary file to the target file
            shutil.move(temp_path, file_path)
            logger.debug(f"Atomically wrote to file: {file_path}")
            
        except Exception as e:
            # Clean up the temporary file on error
            temp_file.close()
            if os.path.exists(temp_path):
                os.unlink(temp_path)
                logger.debug(f"Deleted temporary file: {temp_path} after error")
            raise e
            
    except (OSError, IOError) as e:
        error_msg = f"Error during atomic write to {file_path}: {str(e)}"
        logger.error(error_msg)
        raise FileOperationError(error_msg) from e

def ensure_dir(directory: str) -> bool:
    """Ensure a directory exists
    
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
        return True
    except OSError as e:
        error_msg = f"Error creating directory {directory}: {str(e)}"
        logger.error(error_msg)
        raise FileOperationError(error_msg) from e

def check_file_exists(file_path: str, min_size: int = 0) -> bool:
    """Check if a file exists and has a minimum size
    
    Args:
        file_path: Path to the file
        min_size: Minimum file size in bytes
        
    Returns:
        True if file exists and meets size requirement
    """
    try:
        if not os.path.exists(file_path):
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