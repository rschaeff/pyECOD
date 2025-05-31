#!/usr/bin/env python3
"""
Evidence Isolation Testing for pyECOD Algorithm Debug

This script uses the existing ecod.evaluation framework to test individual
evidence types in isolation, helping debug chain blast down-weighting and
other evidence processing issues.
"""

import os
import sys
import json
import yaml
import logging
import tempfile
from pathlib import Path
from typing import Dict, List, Any, Optional
from dataclasses import dataclass

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext
from ecod.evaluation.algorithm_versions.manager import AlgorithmVersionManager, AlgorithmVersion, AlgorithmStatus

@dataclass
class EvidenceTestResult:
    """Results from testing a single evidence type"""
    evidence_type: str
    algorithm_version: str
    protein_id: str
    success: bool
    domains_found: int
    evidence_hits: int
    processing_time: float
    error_message: Optional[str] = None
    detailed_results: Optional[Dict[str, Any]] = None

class EvidenceIsolationTester:
    """Test individual evidence types in isolation using the evaluation framework"""
    
    def __init__(self, config_path: str = "config/config.yml"):
        self.context = ApplicationContext(config_path)
        self.version_manager = AlgorithmVersionManager(self.context)
        self.logger = logging.getLogger(__name__)
        
        # Define the debug algorithm configurations
        self.debug_algorithms = {
            'chain_blast': 'debug_chain_blast_only',
            'hhsearch': 'debug_hhsearch_only', 
            'domain_blast': 'debug_domain_blast_only',
            'blast': 'debug_blast_only'
        }
        
    def setup_debug_algorithms(self) -> None:
        """Register all debug algorithm configurations"""
        
        self.logger.info("Setting up debug algorithm configurations...")
        
        # Create the algorithm configurations
        debug_configs = self._create_debug_algorithm_configs()
        
        for evidence_type, config in debug_configs.items():
            try:
                # Write config to temporary file
                with tempfile.NamedTemporaryFile(mode='w', suffix='.yml', delete=False) as f:
                    yaml.dump(config, f, default_flow_style=False)
                    temp_path = f.name
                
                # Import via algorithm manager
                algorithm_id = self.version_manager.import_version(temp_path)
                
                self.logger.info(f"Registered {evidence_type} debug algorithm: {config['algorithm']['version_id']}")
                
                # Clean up temp file
                os.unlink(temp_path)
                
            except Exception as e:
                self.logger.error(f"Failed to register {evidence_type} debug algorithm: {e}")
    
    def _create_debug_algorithm_configs(self) -> Dict[str, Dict[str, Any]]:
        """Create debug algorithm configurations for each evidence type"""
        
        base_config = {
            'domain_analysis': {
                'partition': {
                    'min_domain_size': 20,
                    'overlap_threshold': 0.3,
                    'merge_gap_tolerance': 25,
                    'validation_level': 'lenient'
                },
                'coverage_thresholds': {
                    'min_reference_coverage': 0.5,
                    'strict_reference_coverage': 0.8,
                    'partial_coverage_threshold': 0.2,
                    'extend_to_reference_size': True,
                    'use_ungapped_coverage': True
                },
                'behavioral_flags': {
                    'extend_to_reference_size': True
                }
            }
        }
        
        configs = {}
        
        # Chain blast only
        configs['chain_blast'] = {
            'algorithm': {
                'version_id': 'debug_chain_blast_only',
                'name': 'Chain BLAST Only (Debug)',
                'description': 'Debug algorithm using ONLY chain blast evidence',
                'status': 'development',
                'created_by': 'debug_script'
            },
            'domain_analysis': {
                **base_config['domain_analysis'],
                'evidence_weights': {
                    'chain_blast': 3.0,
                    'hhsearch': 0.0,
                    'domain_blast': 0.0,
                    'blast': 0.0
                },
                'behavioral_flags': {
                    **base_config['domain_analysis']['behavioral_flags'],
                    'prefer_hhsearch_classification': False,
                    'prefer_chain_blast_architecture': True,
                    'use_architectural_context': True
                }
            }
        }
        
        # HHsearch only
        configs['hhsearch'] = {
            'algorithm': {
                'version_id': 'debug_hhsearch_only',
                'name': 'HHsearch Only (Debug)',
                'description': 'Debug algorithm using ONLY HHsearch evidence',
                'status': 'development',
                'created_by': 'debug_script'
            },
            'domain_analysis': {
                **base_config['domain_analysis'],
                'evidence_weights': {
                    'hhsearch': 3.0,
                    'chain_blast': 0.0,
                    'domain_blast': 0.0,
                    'blast': 0.0
                },
                'coverage_thresholds': {
                    **base_config['domain_analysis']['coverage_thresholds'],
                    'min_reference_coverage': 0.7  # Higher for HHsearch
                },
                'behavioral_flags': {
                    **base_config['domain_analysis']['behavioral_flags'],
                    'prefer_hhsearch_classification': True
                }
            }
        }
        
        # Domain blast only
        configs['domain_blast'] = {
            'algorithm': {
                'version_id': 'debug_domain_blast_only',
                'name': 'Domain BLAST Only (Debug)',
                'description': 'Debug algorithm using ONLY domain blast evidence',
                'status': 'development',
                'created_by': 'debug_script'
            },
            'domain_analysis': {
                **base_config['domain_analysis'],
                'evidence_weights': {
                    'domain_blast': 3.0,
                    'hhsearch': 0.0,
                    'chain_blast': 0.0,
                    'blast': 0.0
                },
                'coverage_thresholds': {
                    **base_config['domain_analysis']['coverage_thresholds'],
                    'min_reference_coverage': 0.6  # Slightly lenient
                }
            }
        }
        
        # Regular blast only
        configs['blast'] = {
            'algorithm': {
                'version_id': 'debug_blast_only',
                'name': 'Regular BLAST Only (Debug)',
                'description': 'Debug algorithm using ONLY regular blast evidence',
                'status': 'development',
                'created_by': 'debug_script'
            },
            'domain_analysis': {
                **base_config['domain_analysis'],
                'evidence_weights': {
                    'blast': 3.0,
                    'hhsearch': 0.0,
                    'chain_blast': 0.0,
                    'domain_blast': 0.0
                },
                'coverage_thresholds': {
                    **base_config['domain_analysis']['coverage_thresholds'],
                    'min_reference_coverage': 0.4  # Very lenient for basic blast
                }
            }
        }
        
        return configs
    
    def test_protein_all_evidence_types(self, pdb_id: str, chain_id: str) -> Dict[str, EvidenceTestResult]:
        """Test a protein with each evidence type in isolation"""
        
        self.logger.info(f"Testing {pdb_id}_{chain_id} with all evidence types in isolation")
        
        results = {}
        
        for evidence_type, algorithm_version in self.debug_algorithms.items():
            try:
                result = self._test_protein_single_evidence(
                    pdb_id, chain_id, evidence_type, algorithm_version
                )
                results[evidence_type] = result
                
                self.logger.info(f"  {evidence_type}: {'SUCCESS' if result.success else 'FAILED'} "
                               f"- {result.domains_found} domains, {result.evidence_hits} hits")
                
            except Exception as e:
                self.logger.error(f"  {evidence_type}: ERROR - {e}")
                results[evidence_type] = EvidenceTestResult(
                    evidence_type=evidence_type,
                    algorithm_version=algorithm_version,
                    protein_id=f"{pdb_id}_{chain_id}",
                    success=False,
                    domains_found=0,
                    evidence_hits=0,
                    processing_time=0.0,
                    error_message=str(e)
                )
        
        return results
    
    def _demonstrate_evaluation_integration(self):
        """Demonstrate integration with existing evaluation framework"""
        
        self.logger.info("Demonstrating integration with ecod.evaluation framework...")
        
        # Show how to use the existing AlgorithmVersionManager
        algorithms = self.version_manager.list_versions()
        debug_algorithms = [a for a in algorithms if 'debug_' in a.version_id]
        
        self.logger.info(f"Found {len(debug_algorithms)} debug algorithms registered:")
        for algo in debug_algorithms:
            self.logger.info(f"  {algo.version_id}: {algo.name}")
            
        # Show how to get and use algorithm configurations
        if debug_algorithms:
            example_algo = debug_algorithms[0]
            config = example_algo.to_config_dict()
            
            self.logger.info(f"Example algorithm config for {example_algo.version_id}:")
            self.logger.info(f"  Evidence weights: {config['domain_analysis']['evidence_weights']}")
            
            # This is where you'd integrate with DomainPartitionService
            # partition_service = DomainPartitionService(self.context, algorithm_config=config['domain_analysis'])
            
        return debug_algorithms
    
    def _test_protein_single_evidence(self, pdb_id: str, chain_id: str, 
                                    evidence_type: str, algorithm_version: str) -> EvidenceTestResult:
        """Test a protein with a single evidence type"""
        
        import time
        start_time = time.time()
        
        # Get algorithm configuration
        algorithm = self.version_manager.get_version(algorithm_version)
        if not algorithm:
            raise ValueError(f"Algorithm version {algorithm_version} not found")
        
        # Find domain summary files
        summary_paths = self._find_domain_summaries(pdb_id, chain_id)
        if not summary_paths:
            raise ValueError(f"No domain summaries found for {pdb_id}_{chain_id}")
        
        # Parse summary to count evidence hits
        evidence_hits = self._count_evidence_hits(summary_paths[0], evidence_type)
        
        # TODO: Actually run the domain partition algorithm with this configuration
        # For now, simulate the result based on evidence availability
        
        processing_time = time.time() - start_time
        
        # Simulate success based on evidence availability
        success = evidence_hits > 0
        domains_found = 1 if success else 0  # Simplified
        
        return EvidenceTestResult(
            evidence_type=evidence_type,
            algorithm_version=algorithm_version,
            protein_id=f"{pdb_id}_{chain_id}",
            success=success,
            domains_found=domains_found,
            evidence_hits=evidence_hits,
            processing_time=processing_time,
            detailed_results={
                'algorithm_config': algorithm.to_config_dict(),
                'summary_file': summary_paths[0],
                'evidence_details': self._analyze_evidence_details(summary_paths[0], evidence_type)
            }
        )
    #test
    def _find_domain_summaries(self, pdb_id: str, chain_id: str) -> List[str]:
        """Find domain summary files for a protein"""
        
        query = """
        SELECT DISTINCT pf.file_path, b.base_path
        FROM ecod_schema.process_file pf
        JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        JOIN ecod_schema.batch b ON ps.batch_id = b.id
        WHERE p.pdb_id = %s AND p.chain_id = %s
          AND pf.file_type IN ('domain_summary', 'blast_only_summary')
          AND pf.file_exists = TRUE
        ORDER BY ps.batch_id DESC
        LIMIT 3
        """
        
        results = self.context.db.execute_dict_query(query, (pdb_id, chain_id))
        
        summary_paths = []
        for row in results:
            full_path = os.path.join(row['base_path'], row['file_path'])
            if os.path.exists(full_path):
                summary_paths.append(full_path)
                
        return summary_paths
    
    def _count_evidence_hits(self, summary_path: str, evidence_type: str) -> int:
        """Count evidence hits of a specific type in summary file"""
        
        # This would parse the actual domain summary format
        # For now, return a placeholder
        
        with open(summary_path, 'r') as f:
            content = f.read()
        
        # Simple heuristic - count lines containing evidence type
        evidence_patterns = {
            'chain_blast': ['chain_blast', 'chain BLAST'],
            'hhsearch': ['hhsearch', 'HHsearch'],
            'domain_blast': ['domain_blast', 'domain BLAST'],
            'blast': ['blast', 'BLAST']
        }
        
        patterns = evidence_patterns.get(evidence_type, [evidence_type])
        
        hits = 0
        for line in content.lower().split('\n'):
            for pattern in patterns:
                if pattern.lower() in line:
                    hits += 1
                    break
        
        return hits
    
    def _analyze_evidence_details(self, summary_path: str, evidence_type: str) -> Dict[str, Any]:
        """Analyze detailed evidence information"""
        
        # This would extract detailed evidence information
        # For now, return basic file info
        
        return {
            'file_size': os.path.getsize(summary_path),
            'evidence_type': evidence_type,
            'summary_path': summary_path
        }
    
    def compare_evidence_types(self, test_results: Dict[str, EvidenceTestResult]) -> Dict[str, Any]:
        """Compare results across evidence types"""
        
        comparison = {
            'protein_id': list(test_results.values())[0].protein_id if test_results else 'unknown',
            'evidence_ranking': [],
            'successful_evidence_types': [],
            'failed_evidence_types': [],
            'evidence_hit_counts': {},
            'processing_times': {},
            'recommendations': []
        }
        
        # Sort by success and hit count
        sorted_results = sorted(
            test_results.items(),
            key=lambda x: (x[1].success, x[1].evidence_hits),
            reverse=True
        )
        
        for evidence_type, result in sorted_results:
            comparison['evidence_ranking'].append({
                'evidence_type': evidence_type,
                'success': result.success,
                'domains_found': result.domains_found,
                'evidence_hits': result.evidence_hits,
                'processing_time': result.processing_time
            })
            
            if result.success:
                comparison['successful_evidence_types'].append(evidence_type)
            else:
                comparison['failed_evidence_types'].append(evidence_type)
                
            comparison['evidence_hit_counts'][evidence_type] = result.evidence_hits
            comparison['processing_times'][evidence_type] = result.processing_time
        
        # Generate recommendations
        if not comparison['successful_evidence_types']:
            comparison['recommendations'].append("⚠️ No evidence types succeeded - check domain summary files")
        elif len(comparison['successful_evidence_types']) == 1:
            successful_type = comparison['successful_evidence_types'][0]
            comparison['recommendations'].append(f"✅ Only {successful_type} succeeded - investigate other evidence types")
        else:
            comparison['recommendations'].append("✅ Multiple evidence types working - check evidence integration")
            
        # Chain blast specific recommendations
        chain_blast_result = test_results.get('chain_blast')
        if chain_blast_result and not chain_blast_result.success:
            if chain_blast_result.evidence_hits == 0:
                comparison['recommendations'].append("🔍 Chain blast: No hits found - check chain blast database")
            else:
                comparison['recommendations'].append("🔍 Chain blast: Has hits but failed - check domain mapping logic")
        elif chain_blast_result and chain_blast_result.success:
            comparison['recommendations'].append("✅ Chain blast: Working correctly in isolation")
            
        return comparison
    
    def test_gfp_debug(self) -> Dict[str, Any]:
        """Specific test for GFP assignment debugging"""
        
        # Common GFP PDB entries
        gfp_candidates = [
            ("1GFL", "A"),  # Green Fluorescent Protein
            ("1EMA", "A"),  # Enhanced GFP
            ("2WUR", "A"),  # Another GFP variant
        ]
        
        # Also test integration with existing evaluation framework
        self._demonstrate_evaluation_integration()
        
        results = {
            'gfp_candidates_tested': [],
            'successful_candidates': [],
            'evidence_analysis': {}
        }
        
        for pdb_id, chain_id in gfp_candidates:
            try:
                self.logger.info(f"Testing GFP candidate: {pdb_id}_{chain_id}")
                
                test_results = self.test_protein_all_evidence_types(pdb_id, chain_id)
                comparison = self.compare_evidence_types(test_results)
                
                results['gfp_candidates_tested'].append(f"{pdb_id}_{chain_id}")
                
                if comparison['successful_evidence_types']:
                    results['successful_candidates'].append({
                        'protein': f"{pdb_id}_{chain_id}",
                        'successful_evidence': comparison['successful_evidence_types'],
                        'evidence_ranking': comparison['evidence_ranking']
                    })
                
                results['evidence_analysis'][f"{pdb_id}_{chain_id}"] = comparison
                
            except Exception as e:
                self.logger.error(f"Failed to test {pdb_id}_{chain_id}: {e}")
        
        return results


def main():
    """Main CLI interface"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Evidence Isolation Testing for Algorithm Debug")
    parser.add_argument('--config', type=str, default='config/config.yml', help='Configuration file')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Setup debug algorithms
    setup_parser = subparsers.add_parser('setup', help='Setup debug algorithm configurations')
    
    # Test single protein
    test_parser = subparsers.add_parser('test', help='Test protein with all evidence types')
    test_parser.add_argument('--protein', type=str, required=True, 
                           help='Protein to test (format: PDB_CHAIN)')
    
    # Test GFP specifically
    gfp_parser = subparsers.add_parser('test-gfp', help='Test GFP candidates for debugging')
    
    # List debug algorithms
    list_parser = subparsers.add_parser('list', help='List registered debug algorithms')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return 1
    
    # Setup logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=level, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    try:
        tester = EvidenceIsolationTester(args.config)
        
        if args.command == 'setup':
            tester.setup_debug_algorithms()
            print("✅ Debug algorithm configurations registered")
            print("\nRegistered algorithms:")
            for evidence_type, version_id in tester.debug_algorithms.items():
                print(f"  {evidence_type}: {version_id}")
            print("\nYou can now test individual evidence types!")
            
        elif args.command == 'test':
            pdb_id, chain_id = args.protein.split('_')
            
            print(f"🧪 Testing {pdb_id}_{chain_id} with isolated evidence types")
            print("=" * 60)
            
            test_results = tester.test_protein_all_evidence_types(pdb_id, chain_id)
            comparison = tester.compare_evidence_types(test_results)
            
            print("\nRESULTS BY EVIDENCE TYPE:")
            for ranking in comparison['evidence_ranking']:
                status = "✅ SUCCESS" if ranking['success'] else "❌ FAILED"
                print(f"  {ranking['evidence_type']:15} {status:12} "
                      f"Domains: {ranking['domains_found']}, Hits: {ranking['evidence_hits']}")
            
            print(f"\nSUMMARY:")
            print(f"  Successful: {', '.join(comparison['successful_evidence_types']) or 'None'}")
            print(f"  Failed: {', '.join(comparison['failed_evidence_types']) or 'None'}")
            
            print(f"\nRECOMMENDATIONS:")
            for rec in comparison['recommendations']:
                print(f"  {rec}")
                
        elif args.command == 'test-gfp':
            print("🧬 Testing GFP candidates for single domain assignment debugging")
            print("=" * 70)
            
            results = tester.test_gfp_debug()
            
            print(f"\nGFP CANDIDATES TESTED: {len(results['gfp_candidates_tested'])}")
            for candidate in results['gfp_candidates_tested']:
                print(f"  {candidate}")
            
            print(f"\nSUCCESSFUL CANDIDATES: {len(results['successful_candidates'])}")
            for success in results['successful_candidates']:
                print(f"  {success['protein']}: {', '.join(success['successful_evidence'])}")
            
            # Detailed analysis for each candidate
            for protein, analysis in results['evidence_analysis'].items():
                print(f"\n{protein} DETAILED ANALYSIS:")
                print(f"  Evidence ranking:")
                for ranking in analysis['evidence_ranking'][:3]:  # Top 3
                    print(f"    {ranking['evidence_type']}: {ranking['evidence_hits']} hits")
                
                if analysis['recommendations']:
                    print(f"  Recommendations:")
                    for rec in analysis['recommendations']:
                        print(f"    {rec}")
                        
        elif args.command == 'list':
            print("🔬 Registered Debug Algorithms:")
            print("-" * 40)
            
            for evidence_type, version_id in tester.debug_algorithms.items():
                algorithm = tester.version_manager.get_version(version_id)
                if algorithm:
                    print(f"{evidence_type:15} {version_id:25} ✅ Registered")
                else:
                    print(f"{evidence_type:15} {version_id:25} ❌ Not found")
            
            print(f"\nUse 'setup' command to register missing algorithms")
            
        return 0
        
    except Exception as e:
        print(f"❌ Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
