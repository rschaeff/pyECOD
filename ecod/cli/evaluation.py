#!/usr/bin/env python3
"""
Functional CLI module for algorithm evaluation and testing

This implements the actual functionality for the evaluation CLI commands.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Optional

from ecod.cli.base_command import BaseCommand, handle_command_errors


class EvaluationCommand(BaseCommand):
    """Evaluation command group for algorithm testing and validation"""

    def setup_parser(self, parser: argparse.ArgumentParser) -> None:
        """Set up the evaluation command parser"""
        subparsers = parser.add_subparsers(dest='eval_command', help='Evaluation commands')

        # Algorithm version management
        self._setup_algorithm_parser(subparsers)

        # Test set management (placeholder for now)
        self._setup_testset_parser(subparsers)

        # Evaluation runs (placeholder for now)
        self._setup_run_parser(subparsers)

    def _setup_algorithm_parser(self, subparsers):
        """Set up algorithm version management commands"""
        algo_parser = subparsers.add_parser('algorithm', help='Algorithm version management')
        algo_subs = algo_parser.add_subparsers(dest='algo_action', help='Algorithm actions')

        # List algorithms
        list_parser = algo_subs.add_parser('list', help='List algorithm versions')
        list_parser.add_argument('--status', type=str,
                               choices=['development', 'testing', 'production', 'deprecated'],
                               help='Filter by status')
        list_parser.add_argument('--format', type=str, choices=['table', 'json'],
                               default='table', help='Output format')

        # Register new algorithm
        register_parser = algo_subs.add_parser('register', help='Register new algorithm version')
        register_parser.add_argument('config_file', type=str,
                                   help='Algorithm configuration file')
        register_parser.add_argument('--force', action='store_true',
                                   help='Force registration even if version exists')

        # Show algorithm details
        show_parser = algo_subs.add_parser('show', help='Show algorithm version details')
        show_parser.add_argument('version_id', type=str, help='Algorithm version ID')
        show_parser.add_argument('--format', type=str, choices=['yaml', 'json'],
                               default='yaml', help='Output format')

        # Promote algorithm
        promote_parser = algo_subs.add_parser('promote', help='Promote algorithm version')
        promote_parser.add_argument('version_id', type=str, help='Algorithm version ID')
        promote_parser.add_argument('status', type=str,
                                  choices=['testing', 'production', 'deprecated'],
                                  help='New status')
        promote_parser.add_argument('--force', action='store_true',
                                  help='Force promotion')

        # Export algorithm
        export_parser = algo_subs.add_parser('export', help='Export algorithm configuration')
        export_parser.add_argument('version_id', type=str, help='Algorithm version ID')
        export_parser.add_argument('output_file', type=str, help='Output file path')

        # Import algorithm
        import_parser = algo_subs.add_parser('import', help='Import algorithm from file')
        import_parser.add_argument('config_file', type=str, help='Configuration file path')

    def _setup_testset_parser(self, subparsers):
        """Set up test set management commands (placeholder)"""
        testset_parser = subparsers.add_parser('testset', help='Test set management')
        testset_subs = testset_parser.add_subparsers(dest='testset_action', help='Test set actions')

        # Placeholder commands
        list_parser = testset_subs.add_parser('list', help='List test sets (not implemented)')

    def _setup_run_parser(self, subparsers):
        """Set up evaluation run commands (placeholder)"""
        run_parser = subparsers.add_parser('run', help='Run evaluations')
        run_subs = run_parser.add_subparsers(dest='run_action', help='Run types')

        # Placeholder commands
        quick_parser = run_subs.add_parser('quick', help='Quick test (not implemented)')

    @handle_command_errors
    def run_command(self, args: argparse.Namespace) -> int:
        """Run the evaluation command"""
        if not args.eval_command:
            self.logger.error("No evaluation command specified")
            return 1

        if args.eval_command == 'algorithm':
            return self._handle_algorithm_commands(args)
        elif args.eval_command == 'testset':
            return self._handle_testset_commands(args)
        elif args.eval_command == 'run':
            return self._handle_run_commands(args)
        else:
            self.logger.error(f"Unknown evaluation command: {args.eval_command}")
            return 1

    def _handle_algorithm_commands(self, args) -> int:
        """Handle algorithm version management commands"""
        if not args.algo_action:
            self.logger.error("No algorithm action specified")
            return 1

        # Import here to avoid circular imports
        from ecod.evaluation.algorithm_versions import AlgorithmVersionManager, AlgorithmStatus

        manager = AlgorithmVersionManager(self.context)

        if args.algo_action == 'list':
            return self._list_algorithms(manager, args)
        elif args.algo_action == 'register':
            return self._register_algorithm(manager, args)
        elif args.algo_action == 'show':
            return self._show_algorithm(manager, args)
        elif args.algo_action == 'promote':
            return self._promote_algorithm(manager, args)
        elif args.algo_action == 'export':
            return self._export_algorithm(manager, args)
        elif args.algo_action == 'import':
            return self._import_algorithm(manager, args)
        else:
            self.logger.error(f"Unknown algorithm action: {args.algo_action}")
            return 1

    def _list_algorithms(self, manager, args) -> int:
        """List algorithm versions"""
        from ecod.evaluation.algorithm_versions import AlgorithmStatus

        status_filter = AlgorithmStatus(args.status) if args.status else None
        algorithms = manager.list_versions(status=status_filter)

        if args.format == 'json':
            # Convert to JSON-serializable format
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
            # Table format
            if not algorithms:
                print("No algorithm versions found.")
                return 0

            print(f"{'Version ID':25} {'Name':30} {'Status':12} {'Parent':20} {'Created':20}")
            print("-" * 107)

            for algo in algorithms:
                created_str = algo.created_at.strftime('%Y-%m-%d %H:%M') if algo.created_at else 'Unknown'
                parent_str = algo.parent_version or 'None'

                print(f"{algo.version_id[:24]:25} {algo.name[:29]:30} {algo.status.value:12} "
                      f"{parent_str[:19]:20} {created_str:20}")

        return 0

    def _register_algorithm(self, manager, args) -> int:
        """Register new algorithm version"""
        config_path = Path(args.config_file)

        if not config_path.exists():
            self.logger.error(f"Configuration file not found: {config_path}")
            return 1

        try:
            from ecod.evaluation.algorithm_versions import AlgorithmVersion

            # Load algorithm from config file
            algorithm = AlgorithmVersion.from_config_file(str(config_path))

            # Check if version already exists
            existing = manager.get_version(algorithm.version_id)
            if existing and not args.force:
                self.logger.error(f"Algorithm version {algorithm.version_id} already exists. Use --force to overwrite.")
                return 1

            if existing and args.force:
                self.logger.warning(f"Version {algorithm.version_id} already exists, but --force specified")
                # For now, we don't support overwriting - would need delete functionality
                self.logger.error("Overwriting existing versions not yet supported")
                return 1

            algorithm_id = manager.register_version(algorithm)
            print(f"Successfully registered algorithm version {algorithm.version_id} with ID {algorithm_id}")

            return 0

        except Exception as e:
            self.logger.error(f"Failed to register algorithm: {e}")
            return 1

    def _show_algorithm(self, manager, args) -> int:
        """Show algorithm version details"""
        algorithm = manager.get_version(args.version_id)

        if not algorithm:
            self.logger.error(f"Algorithm version {args.version_id} not found")
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
            # YAML format
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

        return 0

    def _promote_algorithm(self, manager, args) -> int:
        """Promote algorithm version"""
        try:
            from ecod.evaluation.algorithm_versions import AlgorithmStatus

            new_status = AlgorithmStatus(args.status)

            if not args.force:
                # Show current status and ask for confirmation
                algorithm = manager.get_version(args.version_id)
                if not algorithm:
                    self.logger.error(f"Algorithm version {args.version_id} not found")
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
                self.logger.error("Promotion failed")
                return 1

        except Exception as e:
            self.logger.error(f"Failed to promote algorithm: {e}")
            return 1

    def _export_algorithm(self, manager, args) -> int:
        """Export algorithm configuration"""
        try:
            manager.export_version(args.version_id, args.output_file)
            print(f"Exported algorithm {args.version_id} to {args.output_file}")
            return 0
        except Exception as e:
            self.logger.error(f"Failed to export algorithm: {e}")
            return 1

    def _import_algorithm(self, manager, args) -> int:
        """Import algorithm from configuration file"""
        try:
            algorithm_id = manager.import_version(args.config_file)
            print(f"Successfully imported algorithm with ID {algorithm_id}")
            return 0
        except Exception as e:
            self.logger.error(f"Failed to import algorithm: {e}")
            return 1

    def _handle_testset_commands(self, args) -> int:
        """Handle test set management commands"""
        self.logger.info("Test set commands not yet implemented")
        return 1

    def _handle_run_commands(self, args) -> int:
        """Handle evaluation run commands"""
        self.logger.info("Run commands not yet implemented")
        return 1


def setup_parser(parser: argparse.ArgumentParser) -> None:
    """
    Set up the evaluation command parser

    This function is called by the main CLI system to configure
    the evaluation command group.
    """
    command = EvaluationCommand()
    command.setup_parser(parser)


def run_command(args: argparse.Namespace) -> int:
    """
    Run evaluation command

    This function is called by the main CLI system to execute
    evaluation commands.
    """
    try:
        command = EvaluationCommand(args.config)
        return command.run_command(args)
    except Exception as e:
        import logging
        logger = logging.getLogger(__name__)
        logger.error(f"Evaluation command failed: {e}", exc_info=True)
        return 1
