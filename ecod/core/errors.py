#!/usr/bin/env python3
"""
Exception framework for the ECOD pipeline
Defines a hierarchy of exceptions for different error types
"""
import sys
import traceback
import logging
from typing import Optional, Dict, Any, Callable

logger = logging.getLogger("ecod.errors")

class ECODError(Exception):
    """Base exception for all ECOD pipeline errors"""
    def __init__(self, message: str, details: Optional[Dict[str, Any]] = None):
        self.message = message
        self.details = details or {}
        super().__init__(message)
        
    def log(self, level: int = logging.ERROR) -> None:
        """Log the error with additional details
        
        Args:
            level: Logging level (default: ERROR)
        """
        logger.log(level, f"{self.__class__.__name__}: {self.message}")
        if self.details:
            logger.log(level, f"Details: {self.details}")

class ConfigurationError(ECODError):
    """Error related to configuration issues"""
    pass

class DatabaseError(ECODError):
    """Base class for database-related errors"""
    pass

class ConnectionError(DatabaseError):
    """Error connecting to the database"""
    pass

class QueryError(DatabaseError):
    """Error executing database query"""
    pass

class ValidationError(ECODError):
    """Error validating data or configurations"""
    pass

class ProcessingError(ECODError):
    """Error during data processing"""
    pass

class PipelineError(ECODError):
    """Base class for pipeline process errors"""
    pass

class JobSubmissionError(PipelineError):
    """Error submitting a job to the cluster"""
    pass

class JobExecutionError(PipelineError):
    """Error during job execution"""
    pass

class FileOperationError(ECODError):
    """Error during file operations"""
    pass

class DataError(ECODError):
    """Error related to data integrity or format"""
    pass

def format_error(error: Exception, verbose: bool = False) -> str:
    """Format an error message for display
    
    Args:
        error: Exception object
        verbose: Whether to include detailed information
        
    Returns:
        Formatted error message
    """
    if isinstance(error, ECODError):
        if verbose and error.details:
            return f"{error.__class__.__name__}: {error.message}\nDetails: {error.details}"
        return f"{error.__class__.__name__}: {error.message}"
    elif verbose:
        return f"Unexpected Error ({error.__class__.__name__}): {str(error)}\n{traceback.format_exc()}"
    else:
        return f"Unexpected Error: {str(error)}"

def handle_exceptions(func: Callable) -> Callable:
    """Decorator to handle exceptions in CLI commands
    
    Args:
        func: Function to wrap
        
    Returns:
        Wrapped function with exception handling
    """
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except KeyboardInterrupt:
            logger.info("Operation cancelled by user")
            print("\nOperation cancelled by user", file=sys.stderr)
            return 130  # Standard exit code for SIGINT
        except ECODError as e:
            # Known application error
            e.log()
            print(format_error(e), file=sys.stderr)
            return 1
        except Exception as e:
            # Unexpected error
            logger.error(f"Unexpected error: {str(e)}", exc_info=True)
            print(format_error(e, verbose=False), file=sys.stderr)
            print("See log for details. Run with --verbose for more information.", file=sys.stderr)
            return 2
    return wrapper

def print_error_and_exit(error: Exception, exit_code: int = 1, verbose: bool = False) -> None:
    """Print error message and exit
    
    Args:
        error: Exception object
        exit_code: Exit code
        verbose: Whether to include detailed information
    """
    print(format_error(error, verbose), file=sys.stderr)
    sys.exit(exit_code)