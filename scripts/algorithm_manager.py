#!/usr/bin/env python3
"""
Standalone Algorithm Version Management Script

This script provides algorithm version management functionality without
relying on the CLI framework. Can be used directly or imported.

Usage:
    python scripts/algorithm_manager.py list
    python scripts/algorithm_manager.py register config/algorithms/v1_baseline.yml
    python scripts/algorithm_manager.py show v1.0_baseline
    python scripts/algorithm_manager.py promote v1.0_baseline testing
"""

import os
import sys
import json
import argparse
import logging
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext
from ecod.evaluation.algorithm_versions.manager import (
    AlgorithmVersionManager, 
    AlgorithmVersion, 
    AlgorithmStatus
)


def setup_logging(verbose=False):
    """Setup logging"""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )


def list_algorithms(manager: AlgorithmVersionManager, args):
    """List algorithm versions"""
    status_filter = AlgorithmStatus(args.status) if args.status else None
    algorithms = manager.list_versions(status=status_filter)
    
    if args.format == 'json':
        data = []
        for algo in algorithms:
            algo_dict = {
                'version_id': algo.version_id,
                'name': algo.name,
                'description': algo.description,
                'status': algo.status.value,
                'parent_version': algo.parent_version,
                'created_at': algo.created_at.isoformat() if algo.created_at else None,
                'created_by': algo.created_by
            }
            data.append(algo_dict)
        print(json.dumps(data, indent=2))
    else:
        if not algorithms:
            print("No algorithm versions found.")
            return
        
        print(f"{'Version ID':25} {'Name':30} {'Status':12} {'Parent':20} {'Created':20}")
        print("-" * 107)
        
        for algo in algorithms:
            created_str = algo.created_at.strftime('%Y-%m-%d %H:%M') if algo.created_at else 'Unknown'
            parent_str = algo.parent_version or 'None'
            
            print(f"{algo.version_id[:24]:25} {algo.name[:29]:30} {algo.status.value:12} "
                  f"{parent_str[:19]:20} {created_str:20}")


def register_algorithm(manager: AlgorithmVersionManager, args):
    """Register new algorithm version"""
    config_path = Path(args.config_file)
    
    if not config_path.exists():
        print(f"ERROR: Configuration file not found: {config_path}")
        return 1
    
    try:
        # Load algorithm from config file
        algorithm = AlgorithmVersion.from_config_file(str(config_path))
        
        # Check if version already exists
        existing = manager.get_version(algorithm.version_id)
        if existing and not args.force:
            print(f"ERROR: Algorithm version {algorithm.version_id} already exists. Use --force to overwrite.")
            return 1
        
        algorithm_id = manager.register_version(algorithm)
        print(f"Successfully registered algorithm version {algorithm.version_id} with ID {algorithm_id}")
        return 0
        
    except Exception as e:
        print(f"ERROR: Failed to register algorithm: {e}")
        return 1


def show_algorithm(manager: AlgorithmVersionManager, args):
    """Show algorithm version details"""
    algorithm = manager.get_version(args.version_id)
    
    if not algorithm:
        print(f"ERROR: Algorithm version {args.version_id} not found")
        return 1
    
    if args.format == 'json':
        config_dict = algorithm.to_config_dict()
        config_dict.update({
            'metadata': {
                'database_id': algorithm.database_id,
                'status': algorithm.status.value,
                'parent_version': algorithm.parent_version,
                'created_at': algorithm.created_at.isoformat() if algorithm.created_at else None,
                'created_by': algorithm.created_by,
                'notes': algorithm.notes
            }
        })
        print(json.dumps(config_dict, indent=2))
    else:
        try:
            import yaml
            config_dict = algorithm.to_config_dict()
            config_dict.update({
                'metadata': {
                    'database_id': algorithm.database_id,
                    'status': algorithm.status.value,
                    'parent_version': algorithm.parent_version,
                    'created_at': algorithm.created_at.isoformat() if algorithm.created_at else None,
                    'created_by': algorithm.created_by,
                    'notes': algorithm.notes
                }
            })
            print(yaml.dump(config_dict, default_flow_style=False))
        except ImportError:
            print("ERROR: PyYAML not available, try --format json")
            return 1
    
    return 0


def promote_algorithm(manager: AlgorithmVersionManager, args):
    """Promote algorithm version"""
    try:
        new_status = AlgorithmStatus(args.status)
        
        if not args.force:
            # Show current status and ask for confirmation
            algorithm = manager.get_version(args.version_id)
            if not algorithm:
                print(f"ERROR: Algorithm version {args.version_id} not found")
                return 1
            
            print(f"Current status: {algorithm.status.value}")
            print(f"New status: {new_status.value}")
            
            confirm = input("Proceed with promotion? (y/N): ")
            if confirm.lower() != 'y':
                print("Promotion cancelled.")
                return 0
        
        success = manager.promote_version(args.version_id, new_status)
        
        if success:
            print(f"Successfully promoted {args.version_id} to {new_status.value}")
            return 0
        else:
            print("ERROR: Promotion failed")
            return 1
            
    except Exception as e:
        print(f"ERROR: Failed to promote algorithm: {e}")
        return 1


def export_algorithm(manager: AlgorithmVersionManager, args):
    """Export algorithm configuration"""
    try:
        manager.export_version(args.version_id, args.output_file)
        print(f"Exported algorithm {args.version_id} to {args.output_file}")
        return 0
    except Exception as e:
        print(f"ERROR: Failed to export algorithm: {e}")
        return 1


def import_algorithm(manager: AlgorithmVersionManager, args):
    """Import algorithm from configuration file"""
    try:
        algorithm_id = manager.import_version(args.config_file)
        print(f"Successfully imported algorithm with ID {algorithm_id}")
        return 0
    except Exception as e:
        print(f"ERROR: Failed to import algorithm: {e}")
        return 1


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Algorithm Version Management for pyECOD",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    %(prog)s list
    %(prog)s list --status production
    %(prog)s register config/algorithms/v1_baseline.yml
    %(prog)s show v1.0_baseline
    %(prog)s promote v1.0_baseline testing
    %(prog)s export v1.0_baseline /tmp/exported.yml
        """
    )
    
    parser.add_argument('--config', type=str, default='config/config.yml',
                       help='Configuration file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Verbose output')
    
    subparsers = parser.add_subparsers(dest='action', help='Available actions')
    
    # List command
    list_parser = subparsers.add_parser('list', help='List algorithm versions')
    list_parser.add_argument('--status', type=str, 
                           choices=['development', 'testing', 'production', 'deprecated'],
                           help='Filter by status')
    list_parser.add_argument('--format', type=str, choices=['table', 'json'], 
                           default='table', help='Output format')
    
    # Register command
    register_parser = subparsers.add_parser('register', help='Register new algorithm version')
    register_parser.add_argument('config_file', type=str, help='Algorithm configuration file')
    register_parser.add_argument('--force', action='store_true',
                               help='Force registration even if version exists')
    
    # Show command
    show_parser = subparsers.add_parser('show', help='Show algorithm version details')
    show_parser.add_argument('version_id', type=str, help='Algorithm version ID')
    show_parser.add_argument('--format', type=str, choices=['yaml', 'json'], 
                           default='yaml', help='Output format')
    
    # Promote command
    promote_parser = subparsers.add_parser('promote', help='Promote algorithm version')
    promote_parser.add_argument('version_id', type=str, help='Algorithm version ID')
    promote_parser.add_argument('status', type=str, 
                              choices=['testing', 'production', 'deprecated'],
                              help='New status')
    promote_parser.add_argument('--force', action='store_true',
                              help='Force promotion without confirmation')
    
    # Export command
    export_parser = subparsers.add_parser('export', help='Export algorithm configuration')
    export_parser.add_argument('version_id', type=str, help='Algorithm version ID')
    export_parser.add_argument('output_file', type=str, help='Output file path')
    
    # Import command
    import_parser = subparsers.add_parser('import', help='Import algorithm from file')
    import_parser.add_argument('config_file', type=str, help='Configuration file path')
    
    args = parser.parse_args()
    
    if not args.action:
        parser.print_help()
        return 1
    
    setup_logging(args.verbose)
    
    try:
        # Initialize context and manager
        context = ApplicationContext(args.config)
        manager = AlgorithmVersionManager(context)
        
        # Execute the requested action
        if args.action == 'list':
            return list_algorithms(manager, args)
        elif args.action == 'register':
            return register_algorithm(manager, args)
        elif args.action == 'show':
            return show_algorithm(manager, args)
        elif args.action == 'promote':
            return promote_algorithm(manager, args)
        elif args.action == 'export':
            return export_algorithm(manager, args)
        elif args.action == 'import':
            return import_algorithm(manager, args)
        else:
            print(f"ERROR: Unknown action: {args.action}")
            return 1
            
    except Exception as e:
        print(f"ERROR: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
