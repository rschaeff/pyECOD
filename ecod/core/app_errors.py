# ecod/core/app_errors.py
import sys
import traceback
import logging
from typing import Optional, Dict, Any, Callable
from ecod.exceptions import ECODError, ConfigurationError, DatabaseError, PipelineError

def format_error(error: Exception, verbose: bool = False) -> str:
    """Format an error message for display
    
    Args:
        error: Exception object
        verbose: Whether to include detailed information
        
    Returns:
        Formatted error message
    """
    if isinstance(error, ConfigurationError):
        return f"Configuration Error: {str(error)}"
    elif isinstance(error, DatabaseError):
        return f"Database Error: {str(error)}"
    elif isinstance(error, PipelineError):
        return f"Pipeline Error: {str(error)}"
    elif isinstance(error, ECODError):
        return f"ECOD Error: {str(error)}"
    else:
        if verbose:
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
        logger = logging.getLogger("ecod.cli")
        try:
            return func(*args, **kwargs)
        except KeyboardInterrupt:
            logger.info("Operation cancelled by user")
            print("\nOperation cancelled by user", file=sys.stderr)
            return 130  # Standard exit code for SIGINT
        except ECODError as e:
            # Known application error
            logger.error(str(e))
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