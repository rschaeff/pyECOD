#!/usr/bin/env python3
"""
Error handling utilities for the ECOD pipeline.
Provides formatting and decorators for consistent error handling.
"""
import sys
import traceback
import logging
from functools import wraps
from typing import Callable, TypeVar, Any, Dict, Optional, Union

from .exceptions import ECODError

# Type variable for decorators
T = TypeVar('T')


def format_error(error: Exception, verbose: bool = False) -> str:
    """Format an error message for display
    
    Args:
        error: Exception object
        verbose: Whether to include detailed information
        
    Returns:
        Formatted error message
    """
    if isinstance(error, ECODError):
        msg = f"{error.__class__.__name__}: {error.message}"
        if verbose and error.details:
            msg += f"\nDetails: {error.details}"
        return msg
    else:
        if verbose:
            return f"Unexpected Error ({error.__class__.__name__}): {str(error)}\n{traceback.format_exc()}"
        else:
            return f"Unexpected Error: {str(error)}"


def handle_exceptions(exit_on_error: bool = False) -> Callable[[Callable[..., T]], Callable[..., Union[T, int]]]:
    """Decorator to handle exceptions in functions
    
    Args:
        exit_on_error: Whether to exit the program on error
        
    Returns:
        Decorated function
    """
    def decorator(func: Callable[..., T]) -> Callable[..., Union[T, int]]:
        @wraps(func)
        def wrapper(*args, **kwargs) -> Union[T, int]:
            logger = logging.getLogger(func.__module__)
            try:
                return func(*args, **kwargs)
            except KeyboardInterrupt:
                logger.info("Operation cancelled by user")
                print("\nOperation cancelled by user", file=sys.stderr)
                if exit_on_error:
                    sys.exit(130)  # Standard exit code for SIGINT
                return 130
            except ECODError as e:
                # Known application error
                logger.error(str(e))
                print(format_error(e), file=sys.stderr)
                if exit_on_error:
                    sys.exit(1)
                return 1
            except Exception as e:
                # Unexpected error
                logger.error(f"Unexpected error: {str(e)}", exc_info=True)
                print(format_error(e, verbose=False), file=sys.stderr)
                print("See log for details. Run with --verbose for more information.", file=sys.stderr)
                if exit_on_error:
                    sys.exit(2)
                return 2
        return wrapper
    return decorator


def cli_error_handler(func: Callable[..., T]) -> Callable[..., T]:
    """Specialized decorator for CLI commands
    
    Args:
        func: CLI command function
        
    Returns:
        Decorated function
    """
    return handle_exceptions(exit_on_error=True)(func)


def log_exception(logger: logging.Logger, 
                 error: Exception, 
                 level: int = logging.ERROR, 
                 context: Optional[Dict[str, Any]] = None) -> None:
    """Log an exception with context
    
    Args:
        logger: Logger instance
        error: Exception object
        level: Logging level
        context: Additional context for the log
    """
    if isinstance(error, ECODError):
        # Include any details from the error
        ctx = {**(error.details or {}), **(context or {})}
        logger.log(level, f"{error.__class__.__name__}: {error.message}", 
                  extra={"context": ctx} if ctx else None, 
                  exc_info=True)
    else:
        # Generic exception
        logger.log(level, f"Unexpected error: {str(error)}", 
                  extra={"context": context} if context else None, 
                  exc_info=True)