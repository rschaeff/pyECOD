# ecod/cli/base_command.py
import logging
import argparse
from typing import Optional
from functools import wraps

from ecod.config import ConfigManager
from ecod.db import DBManager
from ecod.core.context import ApplicationContext
from ecod.error_handlers import cli_error_handler
from ecod.exceptions import ConfigurationError

# Add this decorator function before the BaseCommand class
def handle_command_errors(func):
    """Decorator to handle command errors gracefully"""
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        try:
            return func(self, *args, **kwargs)
        except Exception as e:
            self.logger.error(f"Command error: {str(e)}", exc_info=True)
            return 1
    return wrapper

class BaseCommand:
    """Base class for all CLI commands"""
    
def __init__(self, config_path: Optional[str] = None):
    """Initialize base command
    
    Args:
        config_path: Optional path to configuration file
    """
    try:
        self.context = ApplicationContext(config_path)
        self.db = self.context.db
        # Change this line to access the config through config_manager
        self.config = self.context.config_manager.config
        self.logger = logging.getLogger(f"ecod.cli.{self.__class__.__name__.lower()}")
    except Exception as e:
        # Create a logger if context initialization fails
        self.logger = logging.getLogger(f"ecod.cli.{self.__class__.__name__.lower()}")
        self.logger.error(f"Failed to initialize command: {str(e)}")
        raise ConfigurationError(f"Failed to initialize command: {str(e)}") from ee
    
    def setup_parser(self, parser: argparse.ArgumentParser) -> None:
        """Set up the argument parser for this command
        
        Args:
            parser: Argument parser to configure
            
        Note:
            Must be implemented by subclasses
        """
        raise NotImplementedError("Subclasses must implement setup_parser")
    
    @handle_command_errors
    def run_command(self, args: argparse.Namespace) -> int:
        """Run the command with the provided arguments
        
        Args:
            args: Command line arguments
            
        Returns:
            Exit code (0 for success, non-zero for errors)
            
        Note:
            Must be implemented by subclasses
        """
        raise NotImplementedError("Subclasses must implement run_command")