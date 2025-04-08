# ecod/core/cli_utils.py
import sys
import logging
import traceback
import json
from typing import Optional, Dict, Any, Callable, TypeVar, Union
from ecod.exceptions import ECODError

T = TypeVar('T')

def handle_command_errors(func: Callable[..., T]) -> Callable[..., Union[T, int]]:
    """Decorator to handle errors in CLI commands
    
    Args:
        func: Command function to wrap
        
    Returns:
        Wrapped function that handles errors
    """
    def wrapper(*args, **kwargs) -> Union[T, int]:
        logger = logging.getLogger("ecod.cli")
        try:
            return func(*args, **kwargs)
        except ECODError as e:
            # Known application error
            logger.error(str(e))
            print(f"Error: {str(e)}", file=sys.stderr)
            return 1
        except Exception as e:
            # Unexpected error
            logger.error(f"Unexpected error: {str(e)}", exc_info=True)
            print(f"Unexpected error: {str(e)}", file=sys.stderr)
            print(f"See log for details", file=sys.stderr)
            return 2
    return wrapper

def format_error_response(error: Exception, include_traceback: bool = False) -> Dict[str, Any]:
    """Format an error into a standard response dictionary
    
    Args:
        error: Exception object
        include_traceback: Whether to include traceback
        
    Returns:
        Error response dictionary
    """
    response = {
        'success': False,
        'error': {
            'type': error.__class__.__name__,
            'message': str(error),
        }
    }
    
    if include_traceback:
        response['error']['traceback'] = traceback.format_exc()
        
    return response

def exit_with_error(message: str, error_code: int = 1) -> None:
    """Print error message and exit with error code
    
    Args:
        message: Error message
        error_code: Exit code
    """
    print(f"Error: {message}", file=sys.stderr)
    sys.exit(error_code)

def print_response(response: Dict[str, Any], json_output: bool = False) -> None:
    """Print response in human-readable or JSON format
    
    Args:
        response: Response dictionary
        json_output: Whether to output as JSON
    """
    if json_output:
        print(json.dumps(response, indent=2))
    else:
        if response.get('success', False):
            if 'message' in response:
                print(response['message'])
            else:
                print("Operation completed successfully")
                
            # Print results if present
            if 'data' in response:
                data = response['data']
                if isinstance(data, dict):
                    for key, value in data.items():
                        print(f"{key}: {value}")
                elif isinstance(data, list):
                    for item in data:
                        print(f"- {item}")
                else:
                    print(data)
        else:
            error = response.get('error', {})
            print(f"Error ({error.get('type', 'Unknown')}): {error.get('message', 'Unknown error')}", file=sys.stderr)