"""
Command-line interface for the ECOD pipeline.

This package provides a unified command-line interface for all ECOD functionality,
organizing commands into logical groups for different aspects of the pipeline.
"""

import logging
from importlib import import_module
from typing import Dict, List, Any, Optional

__all__ = ['get_command_groups', 'get_commands']

# Define command groups and their descriptions
COMMAND_GROUPS = {
    'blast': 'Commands for running and managing BLAST searches',
    'hhsearch': 'Commands for running and managing HHSearch analyses',
    'domain': 'Commands for domain boundary detection and classification',
    'db': 'Database management commands',
    'batch': 'Batch processing commands'
}

# Logger for CLI operations
logger = logging.getLogger("ecod.cli")

def get_command_groups() -> Dict[str, str]:
    """Return all available command groups and their descriptions"""
    return COMMAND_GROUPS

def get_commands(group: str) -> List[str]:
    """Return all commands available in a specific group"""
    try:
        # Import the module for the requested command group
        module = import_module(f"ecod.cli.{group}")
        
        # Get the commands from the module (should be defined in a COMMANDS dict)
        if hasattr(module, 'COMMANDS'):
            return module.COMMANDS
        else:
            logger.warning(f"Command group '{group}' does not define a COMMANDS dictionary")
            return []
    except ImportError as e:
        logger.error(f"Failed to import command group '{group}': {e}")
        return []

def setup_logging(verbosity: int = 0, log_file: Optional[str] = None) -> None:
    """
    Configure logging for CLI commands
    
    Args:
        verbosity: Verbosity level (0=warning, 1=info, 2+=debug)
        log_file: Optional path to log file
    """
    log_level = {
        0: logging.WARNING,
        1: logging.INFO,
        2: logging.DEBUG
    }.get(verbosity, logging.DEBUG)
    
    handlers = [logging.StreamHandler()]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )