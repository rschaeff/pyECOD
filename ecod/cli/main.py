#!/usr/bin/env python3
"""
Main entry point for the ECOD command-line interface
"""

import argparse
import sys
import logging
from typing import List, Optional

from ecod.cli import get_command_groups, setup_logging
from ecod.core.config import ConfigManager

def create_parser() -> argparse.ArgumentParser:
    """Create the main argument parser"""
    parser = argparse.ArgumentParser(
        description="ECOD (Evolutionary Classification Of protein Domains) Pipeline",
        prog="ecod"
    )
    
    # Global options
    parser.add_argument('--config', type=str, help='Path to configuration file')
    parser.add_argument('-v', '--verbose', action='count', default=0, 
                       help='Increase verbosity (can be used multiple times)')
    parser.add_argument('--log-file', type=str, 
                       help='Log to file in addition to stdout')
    parser.add_argument('--version', action='store_true',
                       help='Show version information and exit')
    
    # Create subparsers for command groups
    subparsers = parser.add_subparsers(dest='command_group', help='Command group')
    
    # Add each command group
    for group, description in get_command_groups().items():
        group_parser = subparsers.add_parser(group, help=description)
        
        # Import the module for this command group
        try:
            module = __import__(f"ecod.cli.{group}", fromlist=['setup_parser'])
            if hasattr(module, 'setup_parser'):
                module.setup_parser(group_parser)
            else:
                logging.warning(f"Command group '{group}' does not define setup_parser function")
        except ImportError as e:
            logging.error(f"Failed to import command group '{group}': {e}")
    
    return parser

def main(args: Optional[List[str]] = None) -> int:
    """
    Main entry point for the ECOD CLI
    
    Args:
        args: Command line arguments (uses sys.argv if None)
        
    Returns:
        Exit code (0 for success, non-zero for errors)
    """
    parser = create_parser()
    parsed_args = parser.parse_args(args)
    
    # Handle version request
    if parsed_args.version:
        from ecod import __version__
        print(f"ECOD Pipeline version {__version__}")
        return 0
    
    # Setup logging
    setup_logging(parsed_args.verbose, parsed_args.log_file)
    
    # If no command group specified, show help and exit
    if not parsed_args.command_group:
        parser.print_help()
        return 1
    
    # Import the command group module
    try:
        module = __import__(f"ecod.cli.{parsed_args.command_group}", fromlist=['run_command'])
        if hasattr(module, 'run_command'):
            return module.run_command(parsed_args)
        else:
            logging.error(f"Command group '{parsed_args.command_group}' does not define run_command function")
            return 1
    except ImportError as e:
        logging.error(f"Failed to import command group '{parsed_args.command_group}': {e}")
        return 1
    except Exception as e:
        logging.error(f"Error executing command: {e}", exc_info=True)
        return 1

if __name__ == "__main__":
    sys.exit(main())