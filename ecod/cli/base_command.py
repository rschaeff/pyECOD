# ecod/cli/base_command.py
import logging
import argparse
from typing import Optional

from ecod.config import ConfigManager
from ecod.db import DBManager
from ecod.core.context import ApplicationContext
from ecod.core.cli_utils import handle_command_errors

class BaseCommand:
    """Base class for all CLI commands"""
    
    def __init__(self, config_path: Optional[str] = None):
        """Initialize base command
        
        Args:
            config_path: Optional path to configuration file
        """
        self.context = ApplicationContext(config_path)
        self.db = self.context.db
        self.config = self.context.config
        self.logger = logging.getLogger(f"ecod.cli.{self.__class__.__name__.lower()}")
    
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