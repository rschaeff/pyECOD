# ecod/cli/evaluation.py
from ecod.cli.base_command import BaseCommand
from ecod.evaluation.test_sets import TestSetManager
from ecod.evaluation.algorithm_versions import AlgorithmVersionManager

class EvaluationCommand(BaseCommand):
    """Evaluation command group for algorithm testing and validation"""
    
    def setup_parser(self, parser):
        subparsers = parser.add_subparsers(dest='eval_command')
        
        # Algorithm version management
        version_parser = subparsers.add_parser('algorithm', help='Algorithm version management')
        version_subs = version_parser.add_subparsers(dest='version_action')
        
        # Test set management  
        testset_parser = subparsers.add_parser('testset', help='Test set management')
        testset_subs = testset_parser.add_subparsers(dest='testset_action')
        
        # Evaluation runs
        eval_parser = subparsers.add_parser('run', help='Run evaluations')
        
    def run_command(self, args):
        if args.eval_command == 'algorithm':
            return self._handle_algorithm_commands(args)
        elif args.eval_command == 'testset':
            return self._handle_testset_commands(args)
        elif args.eval_command == 'run':
            return self._handle_evaluation_runs(args)
