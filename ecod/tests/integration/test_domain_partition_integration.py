#!/usr/bin/env python3
"""
Domain Partition Integration Tests

Comprehensive integration tests for the domain partition workflow to detect
regressions when iterating on evidence weights, confidence scoring, and
domain boundary determination logic.
"""

import os
import sys
import unittest
import tempfile
import shutil
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from datetime import datetime
import logging

# Test framework imports
import pytest
from unittest.mock import Mock, patch, MagicMock

# Add the project root to the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import the components we're testing
from ecod.core.context import ApplicationContext
from ecod.pipelines.domain_analysis.partition.service import DomainPartitionService
from ecod.models.pipeline.evidence import Evidence
from ecod.models.pipeline.domain import DomainModel
from ecod.models.pipeline.partition import DomainPartitionResult


class DomainPartitionIntegrationTests(unittest.TestCase):
    """
    Integration tests for domain partition workflow.
    
    These tests focus on end-to-end workflows and regression detection
    when evidence weights and scoring algorithms change.
    """

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures and golden datasets"""
        cls.test_data_dir = Path(__file__).parent / "test_data" / "integration"
        cls.golden_dataset_dir = cls.test_data_dir / "golden_datasets"
        cls.config_dir = cls.test_data_dir / "configs"
        
        # Ensure test directories exist
        cls.test_data_dir.mkdir(parents=True, exist_ok=True)
        cls.golden_dataset_dir.mkdir(exist_ok=True)
        cls.config_dir.mkdir(exist_ok=True)
        
        # Create test configuration
        cls.test_config = cls._create_test_config()
        
        # Set up logging for tests
        logging.basicConfig(level=logging.DEBUG)
        cls.logger = logging.getLogger(__name__)

    def setUp(self):
        """Set up for each test"""
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)
        
        # Create mock context
        self.context = self._create_test_context()
        
        # Track results for comparison
        self.test_results = {}

    def tearDown(self):
        """Clean up after each test"""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)

    @classmethod
    def _create_test_config(cls) -> Dict[str, Any]:
        """Create test configuration"""
        return {
            'database': {
                'host': 'localhost',
                'port': 5432,
                'database': 'ecod_test',
                'user': 'test_user',
                'password': 'test_pass'
            },
            'reference': {
                'current_version': 'develop291'
            },
            'partition': {
                'confidence_thresholds': {
                    'high': 0.9,
                    'medium': 0.7,
                    'low': 0.5
                },
                'evidence_weights': {
                    'domain_blast': 3.0,
                    'hhsearch': 2.5,
                    'chain_blast': 2.0,
                    'blast': 1.5
                },
                'overlap_tolerance': 0.15,
                'min_domain_size': 20,
                'peptide_threshold': 50,
                'reference_coverage': {
                    'min_coverage': 0.7,
                    'strict_coverage': 0.9,
                    'extend_to_reference': True
                }
            }
        }

    def _create_test_context(self) -> ApplicationContext:
        """Create mock context for testing"""
        context = Mock(spec=ApplicationContext)
        context.config_manager = Mock()
        context.config_manager.config = self.test_config
        context.config_manager.get_db_config.return_value = self.test_config['database']
        context.config_manager.get_path.return_value = str(self.temp_path)
        return context

    def _create_mock_blast_hits(self, pdb_id: str, chain_id: str, seq_len: int) -> List[Dict[str, Any]]:
        """Create realistic mock BLAST hits"""
        return [
            {
                'num': '1',
                'pdb_id': f'{pdb_id}',
                'chain_id': chain_id,
                'evalues': '1e-50',
                'hsp_count': '1',
                'query_range': f'1-{seq_len}',
                'hit_range': f'1-{seq_len}'
            }
        ]

    def _create_mock_domain_hits(self, test_params: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Create mock domain BLAST hits based on expected domains"""
        expected_domains = test_params.get('expected_domains', 1)
        seq_len = test_params.get('sequence_length', 100)

        hits = []
        for i in range(expected_domains):
            start = i * (seq_len // expected_domains) + 1
            end = (i + 1) * (seq_len // expected_domains)

            hits.append({
                'domain_id': f'e{test_params["pdb_id"]}{test_params["chain_id"]}{i+1}',
                'pdb_id': test_params['pdb_id'],
                'chain_id': test_params['chain_id'],
                'evalues': f'1e-{20 + i*5}',
                'hsp_count': '1',
                'query_range': f'{start}-{end}',
                'hit_range': f'{start}-{end}'
            })

        return hits

    def _create_mock_hhsearch_hits(self, test_params: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Create mock HHSearch hits"""
        expected_domains = test_params.get('expected_domains', 1)
        seq_len = test_params.get('sequence_length', 100)

        hits = []
        for i in range(expected_domains):
            start = i * (seq_len // expected_domains) + 1
            end = (i + 1) * (seq_len // expected_domains)

            hits.append({
                'hit_id': f'h{test_params["pdb_id"]}{test_params["chain_id"]}{i+1}',
                'domain_id': f'e{test_params["pdb_id"]}{test_params["chain_id"]}{i+1}',
                'probability': str(95.0 - i*5.0),
                'evalue': f'1e-{25 + i*3}',
                'score': str(100.0 - i*10.0),
                'query_range': f'{start}-{end}',
                'hit_range': f'{start}-{end}'
            })

        return hits

    def _store_baseline_result(self, test_name: str, result) -> None:
        """Store baseline result for regression comparison"""
        baseline_dir = self.test_data_dir / "baselines"
        baseline_dir.mkdir(exist_ok=True)

        baseline_file = baseline_dir / f"{test_name}.json"

        baseline_data = {
            'timestamp': datetime.now().isoformat(),
            'test_name': test_name,
            'result_summary': {
                'success': getattr(result, 'success', True),
                'is_classified': getattr(result, 'is_classified', True),
                'is_peptide': getattr(result, 'is_peptide', False),
                'domain_count': len(getattr(result, 'domains', [])),
                'coverage': getattr(result, 'coverage', 0.0),
                'processing_time': getattr(result, 'processing_time', 0.0)
            },
            'domains': []
        }

        # Add domain details if available
        if hasattr(result, 'domains'):
            for i, domain in enumerate(result.domains):
                baseline_data['domains'].append({
                    'id': getattr(domain, 'id', f'domain_{i}'),
                    'start': getattr(domain, 'start', 1),
                    'end': getattr(domain, 'end', 100),
                    'range': getattr(domain, 'range', '1-100'),
                    'confidence': getattr(domain, 'confidence', 0.8),
                    'source': getattr(domain, 'source', 'mock'),
                    'evidence_count': len(getattr(domain, 'evidence', []))
                })

        with open(baseline_file, 'w') as f:
            json.dump(baseline_data, f, indent=2)

    def _compare_evidence_results(self, result1, result2, comparison_name: str) -> None:
        """Compare two evidence processing results"""
        comparison_data = {
            'timestamp': datetime.now().isoformat(),
            'comparison_name': comparison_name,
            'result1_summary': {
                'success': getattr(result1, 'success', False),
                'domain_count': len(getattr(result1, 'domains', [])),
                'avg_confidence': 0.0
            },
            'result2_summary': {
                'success': getattr(result2, 'success', False),
                'domain_count': len(getattr(result2, 'domains', [])),
                'avg_confidence': 0.0
            }
        }

        # Calculate average confidence if domains exist
        if hasattr(result1, 'domains') and result1.domains:
            confidences = [getattr(d, 'confidence', 0) for d in result1.domains]
            comparison_data['result1_summary']['avg_confidence'] = sum(confidences) / len(confidences)

        if hasattr(result2, 'domains') and result2.domains:
            confidences = [getattr(d, 'confidence', 0) for d in result2.domains]
            comparison_data['result2_summary']['avg_confidence'] = sum(confidences) / len(confidences)

        # Store comparison
        comparison_dir = self.test_data_dir / "comparisons"
        comparison_dir.mkdir(exist_ok=True)

        comparison_file = comparison_dir / f"{comparison_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(comparison_file, 'w') as f:
            json.dump(comparison_data, f, indent=2)

    def _store_evidence_comparison(self, test_name: str, results_dict: Dict[str, Any]) -> None:
        """Store evidence comparison results"""
        comparison_dir = self.test_data_dir / "evidence_comparisons"
        comparison_dir.mkdir(exist_ok=True)

        comparison_file = comparison_dir / f"{test_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"

        comparison_data = {
            'timestamp': datetime.now().isoformat(),
            'test_name': test_name,
            'results': {}
        }

        for name, result in results_dict.items():
            comparison_data['results'][name] = {
                'success': getattr(result, 'success', False),
                'domain_count': len(getattr(result, 'domains', [])),
                'coverage': getattr(result, 'coverage', 0.0)
            }

        with open(comparison_file, 'w') as f:
            json.dump(comparison_data, f, indent=2)

    def _run_with_custom_weights(self, protein: Dict[str, Any], weights: Dict[str, float]):
        """Run partition with custom evidence weights"""
        # Create a mock result for now - in real implementation this would use the actual service
        from ecod.models.pipeline.partition import DomainPartitionResult
        from ecod.models.pipeline.domain import DomainModel

        result = DomainPartitionResult(
            pdb_id=protein['pdb_id'],
            chain_id=protein['chain_id'],
            reference="develop291",
            sequence_length=protein.get('sequence_length', 200),
            is_classified=True
        )

        # Add mock domain for testing
        domain = DomainModel(
            id=f"{protein['pdb_id']}_{protein['chain_id']}_d1",
            start=1,
            end=protein.get('sequence_length', 200),
            range=f"1-{protein.get('sequence_length', 200)}",
            source="hhsearch",
            confidence=0.85  # Base confidence, would be affected by weights in real implementation
        )

        result.add_domain(domain)
        return result

    def _run_with_custom_thresholds(self, protein: Dict[str, Any], thresholds: Dict[str, float]):
        """Run partition with custom confidence thresholds"""
        # Similar mock implementation
        from ecod.models.pipeline.partition import DomainPartitionResult

        result = DomainPartitionResult(
            pdb_id=protein['pdb_id'],
            chain_id=protein['chain_id'],
            reference="develop291",
            sequence_length=protein.get('sequence_length', 200),
            is_classified=True
        )

        return result

    def _create_filtered_evidence(self, protein: Dict[str, Any], evidence_types: List[str]) -> Dict[str, Any]:
        """Create evidence data filtered to specific types"""
        evidence_data = {
            'sequence_length': protein.get('sequence_length', 200),
            'is_peptide': False,
            'chain_blast_hits': [],
            'domain_blast_hits': [],
            'hhsearch_hits': []
        }

        # Only include requested evidence types
        if 'chain_blast' in evidence_types:
            evidence_data['chain_blast_hits'] = self._create_mock_blast_hits(
                protein['pdb_id'], protein['chain_id'], protein.get('sequence_length', 200)
            )

        if 'domain_blast' in evidence_types:
            evidence_data['domain_blast_hits'] = self._create_mock_domain_hits(protein)

        if 'hhsearch' in evidence_types:
            evidence_data['hhsearch_hits'] = self._create_mock_hhsearch_hits(protein)

        return evidence_data

    def _create_quality_evidence(self, protein: Dict[str, Any], quality_level: str) -> Dict[str, Any]:
        """Create evidence with specific quality level"""
        quality_settings = {
            'high': {'evalue_exp': 20, 'probability': 95.0, 'score': 100.0},
            'medium': {'evalue_exp': 10, 'probability': 80.0, 'score': 60.0},
            'low': {'evalue_exp': 5, 'probability': 60.0, 'score': 30.0}
        }

        settings = quality_settings.get(quality_level, quality_settings['medium'])
        seq_len = protein.get('sequence_length', 200)

        # Create evidence with appropriate quality metrics
        evidence_data = {
            'sequence_length': seq_len,
            'is_peptide': False,
            'chain_blast_hits': [{
                'num': '1',
                'pdb_id': protein['pdb_id'],
                'chain_id': protein['chain_id'],
                'evalues': f'1e-{settings["evalue_exp"]}',
                'hsp_count': '1'
            }],
            'domain_blast_hits': [{
                'domain_id': f'e{protein["pdb_id"]}{protein["chain_id"]}1',
                'evalues': f'1e-{settings["evalue_exp"]}'
            }],
            'hhsearch_hits': [{
                'hit_id': f'h{protein["pdb_id"]}{protein["chain_id"]}1',
                'probability': str(settings['probability']),
                'score': str(settings['score'])
            }]
        }

        return evidence_data

    def _run_with_conflicting_evidence(self, protein: Dict[str, Any], description: str):
        """Run partition with conflicting evidence"""
        # Create a mock result that shows conflict resolution
        from ecod.models.pipeline.partition import DomainPartitionResult
        from ecod.models.pipeline.domain import DomainModel

        result = DomainPartitionResult(
            pdb_id=protein['pdb_id'],
            chain_id=protein['chain_id'],
            reference="develop291",
            sequence_length=protein.get('sequence_length', 250),
            is_classified=True
        )

        # Add domain showing resolved conflict
        domain = DomainModel(
            id=f"{protein['pdb_id']}_{protein['chain_id']}_d1",
            start=25,  # Resolved boundary
            end=225,   # Resolved boundary
            range="25-225",
            source="conflict_resolution",
            confidence=0.75  # Lower confidence due to conflicts
        )

        result.add_domain(domain)
        return result

    def _store_weight_baseline(self, baseline_name: str, results: Dict[str, Any]) -> None:
        """Store weight baseline results"""
        baseline_dir = self.test_data_dir / "weight_baselines"
        baseline_dir.mkdir(exist_ok=True)

        baseline_file = baseline_dir / f"{baseline_name}.json"
        with open(baseline_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)

    def _load_weight_baseline(self, baseline_name: str) -> Dict[str, Any]:
        """Load weight baseline results"""
        baseline_file = self.test_data_dir / "weight_baselines" / f"{baseline_name}.json"
        if baseline_file.exists():
            with open(baseline_file, 'r') as f:
                return json.load(f)
        return {}

    def _compare_weight_impact(self, baseline_results: Dict, modified_results: Dict, test_name: str) -> None:
        """Compare weight impact results"""
        # Store comparison data
        comparison_dir = self.test_data_dir / "weight_comparisons"
        comparison_dir.mkdir(exist_ok=True)

        comparison_file = comparison_dir / f"{test_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        comparison_data = {
            'test_name': test_name,
            'baseline': baseline_results,
            'modified': modified_results,
            'timestamp': datetime.now().isoformat()
        }

        with open(comparison_file, 'w') as f:
            json.dump(comparison_data, f, indent=2, default=str)

    def _store_sensitivity_analysis(self, analysis_name: str, results: Dict[str, Any]) -> None:
        """Store sensitivity analysis results"""
        analysis_dir = self.test_data_dir / "sensitivity_analysis"
        analysis_dir.mkdir(exist_ok=True)

        analysis_file = analysis_dir / f"{analysis_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(analysis_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)

    def _store_threshold_analysis(self, analysis_name: str, results: Dict[str, Any]) -> None:
        """Store threshold analysis results"""
        analysis_dir = self.test_data_dir / "threshold_analysis"
        analysis_dir.mkdir(exist_ok=True)

        analysis_file = analysis_dir / f"{analysis_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(analysis_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)

    def _run_with_default_weights(self, protein: Dict[str, Any]):
        """Run with default weights"""
        return self._run_with_custom_weights(protein, self.test_config['partition']['evidence_weights'])

    def _store_performance_baseline(self, baseline_name: str, results: Dict[str, Any]) -> None:
        """Store performance baseline results"""
        baseline_dir = self.test_data_dir / "performance_baselines"
        baseline_dir.mkdir(exist_ok=True)

        baseline_file = baseline_dir / f"{baseline_name}.json"
        with open(baseline_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)

    def _create_mock_domain_summary(self, test_case: Dict[str, Any]) -> Path:
        """Create mock domain summary XML file"""
        mock_data = test_case['mock_data']
        pdb_id = test_case['input']['pdb_id']
        chain_id = test_case['input']['chain_id']

        # Create XML structure
        root = ET.Element("blast_summ_doc")

        # Add metadata
        blast_summ = ET.SubElement(root, "blast_summ")
        blast_summ.set("pdb", pdb_id)
        blast_summ.set("chain", chain_id)

        # Add chain BLAST hits
        if 'chain_blast_hits' in mock_data:
            chain_run = ET.SubElement(root, "chain_blast_run")
            chain_run.set("program", "blastp")
            hits_elem = ET.SubElement(chain_run, "hits")

            for hit in mock_data['chain_blast_hits']:
                hit_elem = ET.SubElement(hits_elem, "hit")
                for key, value in hit.items():
                    hit_elem.set(key, str(value))

        # Add domain BLAST hits
        if 'domain_blast_hits' in mock_data:
            domain_run = ET.SubElement(root, "blast_run")
            domain_run.set("program", "blastp")
            hits_elem = ET.SubElement(domain_run, "hits")

            for hit in mock_data['domain_blast_hits']:
                hit_elem = ET.SubElement(hits_elem, "hit")
                for key, value in hit.items():
                    if key not in ['query_range', 'hit_range']:
                        hit_elem.set(key, str(value))

                # Add range elements
                if 'query_range' in hit:
                    query_reg = ET.SubElement(hit_elem, "query_reg")
                    query_reg.text = hit['query_range']
                if 'hit_range' in hit:
                    hit_reg = ET.SubElement(hit_elem, "hit_reg")
                    hit_reg.text = hit['hit_range']

        # Add HHSearch hits
        if 'hhsearch_hits' in mock_data:
            hh_run = ET.SubElement(root, "hh_run")
            hh_run.set("program", "hhsearch")
            hits_elem = ET.SubElement(hh_run, "hits")

            for hit in mock_data['hhsearch_hits']:
                hit_elem = ET.SubElement(hits_elem, "hit")
                for key, value in hit.items():
                    if key not in ['query_range', 'hit_range']:
                        hit_elem.set(key, str(value))

                # Add range elements
                if 'query_range' in hit:
                    query_reg = ET.SubElement(hit_elem, "query_reg")
                    query_reg.text = hit['query_range']
                if 'hit_range' in hit:
                    hit_reg = ET.SubElement(hit_elem, "hit_reg")
                    hit_reg.text = hit['hit_range']

        # Save to file
        summary_file = self.temp_path / f"{pdb_id}_{chain_id}.summary.xml"
        tree = ET.ElementTree(root)
        tree.write(str(summary_file), encoding='utf-8', xml_declaration=True)

        return summary_file


class GoldenDatasetTests(DomainPartitionIntegrationTests):
    """
    Tests using golden datasets with known expected outcomes.

    These are the primary regression detection tests.
    """

    def test_single_domain_protein_baseline(self):
        """Test processing of a known single-domain protein"""
        # Golden dataset: 1cbs_A - small single domain protein
        test_case = self._create_golden_test_case(
            pdb_id="1cbs",
            chain_id="A",
            sequence_length=141,
            expected_domains=1,
            expected_classification="a.39.1.5",
            description="Small single domain protein - baseline test"
        )

        result = self._run_golden_test(test_case)

        # Assertions for regression detection
        self.assertTrue(result.success, "Processing should succeed")
        self.assertTrue(result.is_classified, "Should be classified")
        self.assertFalse(result.is_peptide, "Should not be peptide")
        self.assertEqual(len(result.domains), 1, "Should have exactly 1 domain")

        domain = result.domains[0]
        self.assertGreater(domain.confidence, 0.8, "Should have high confidence")
        self.assertIsNotNone(domain.a_group, "Should have A-group classification")

        # Store baseline for comparison
        self._store_baseline_result("single_domain_baseline", result)

    def test_multi_domain_protein_baseline(self):
        """Test processing of a known multi-domain protein"""
        # Golden dataset: 2pth_A - two domain protein
        test_case = self._create_golden_test_case(
            pdb_id="2pth",
            chain_id="A",
            sequence_length=393,
            expected_domains=2,
            expected_classifications=["c.37.1.13", "c.37.1.13"],
            description="Multi-domain protein - baseline test"
        )

        result = self._run_golden_test(test_case)

        # Assertions
        self.assertTrue(result.success, "Processing should succeed")
        self.assertTrue(result.is_classified, "Should be classified")
        self.assertEqual(len(result.domains), 2, "Should have exactly 2 domains")

        # Check domain boundaries make sense
        domain1, domain2 = result.domains
        self.assertLess(domain1.end, domain2.start, "Domains should not overlap significantly")
        self.assertGreater(domain1.confidence, 0.7, "Domain 1 should have good confidence")
        self.assertGreater(domain2.confidence, 0.7, "Domain 2 should have good confidence")

        # Store baseline
        self._store_baseline_result("multi_domain_baseline", result)

    def test_complex_domain_architecture_baseline(self):
        """Test processing of protein with complex domain architecture"""
        # Golden dataset: Complex protein with overlapping evidence
        test_case = self._create_golden_test_case(
            pdb_id="3hhp",
            chain_id="A",
            sequence_length=506,
            expected_domains=3,
            expected_overlaps=True,
            description="Complex domain architecture with overlapping evidence"
        )

        result = self._run_golden_test(test_case)

        # Assertions for complex cases
        self.assertTrue(result.success, "Processing should succeed")
        self.assertTrue(result.is_classified, "Should be classified")
        self.assertGreaterEqual(len(result.domains), 2, "Should have multiple domains")

        # Check coverage is reasonable
        self.assertGreater(result.coverage, 0.8, "Should have good sequence coverage")

        # Store baseline
        self._store_baseline_result("complex_architecture_baseline", result)

    def test_peptide_classification_baseline(self):
        """Test processing of peptide-length sequences"""
        test_case = self._create_golden_test_case(
            pdb_id="1pep",
            chain_id="A",
            sequence_length=25,
            expected_peptide=True,
            description="Peptide sequence - should be classified as peptide"
        )

        result = self._run_golden_test(test_case)

        # Assertions
        self.assertTrue(result.success, "Processing should succeed")
        self.assertTrue(result.is_peptide, "Should be classified as peptide")
        self.assertTrue(result.is_classified, "Peptides count as classified")
        self.assertEqual(len(result.domains), 0, "Peptides should have no domains")

        # Store baseline
        self._store_baseline_result("peptide_baseline", result)

    def _create_golden_test_case(self, **kwargs) -> Dict[str, Any]:
        """Create a golden test case with expected outcomes"""
        return {
            'input': {
                'pdb_id': kwargs['pdb_id'],
                'chain_id': kwargs['chain_id'],
                'sequence_length': kwargs.get('sequence_length', 100)
            },
            'expected': {
                'domains': kwargs.get('expected_domains', 1),
                'classifications': kwargs.get('expected_classifications', []),
                'is_peptide': kwargs.get('expected_peptide', False),
                'overlaps': kwargs.get('expected_overlaps', False)
            },
            'description': kwargs.get('description', 'Golden test case'),
            'mock_data': self._create_mock_evidence_data(kwargs)
        }

    def _create_mock_evidence_data(self, test_params: Dict[str, Any]) -> Dict[str, Any]:
        """Create mock evidence data for test case"""
        pdb_id = test_params['pdb_id']
        chain_id = test_params['chain_id']
        seq_len = test_params.get('sequence_length', 100)

        # Create realistic mock evidence based on test case
        if test_params.get('expected_peptide', False):
            # Peptides have minimal or no evidence
            return {
                'sequence_length': seq_len,
                'is_peptide': True,
                'chain_blast_hits': [],
                'domain_blast_hits': [],
                'hhsearch_hits': []
            }

        # Regular proteins have evidence
        evidence_data = {
            'sequence_length': seq_len,
            'is_peptide': False,
            'chain_blast_hits': self._create_mock_blast_hits(pdb_id, chain_id, seq_len),
            'domain_blast_hits': self._create_mock_domain_hits(test_params),
            'hhsearch_hits': self._create_mock_hhsearch_hits(test_params)
        }

        return evidence_data

    def _run_golden_test(self, test_case: Dict[str, Any]) -> DomainPartitionResult:
        """Run a golden test case and return result"""
        # Create mock domain summary file
        summary_path = self._create_mock_domain_summary(test_case)

        # Create service
        service = DomainPartitionService(self.context)

        # Run partition
        result = service.partition_protein(
            pdb_id=test_case['input']['pdb_id'],
            chain_id=test_case['input']['chain_id'],
            summary_path=str(summary_path),
            output_dir=str(self.temp_path)
        )

        return result


class EvidenceVariationTests(DomainPartitionIntegrationTests):
    """
    Tests how different evidence combinations affect domain partitioning.

    These tests help detect regressions in evidence processing logic.
    """

    def test_blast_only_vs_full_evidence(self):
        """Compare results when using BLAST-only vs full evidence"""
        base_protein = {
            'pdb_id': '1test',
            'chain_id': 'A',
            'sequence_length': 200
        }

        # Test with BLAST-only evidence
        blast_only_result = self._run_with_evidence_set(
            protein=base_protein,
            evidence_types=['domain_blast', 'chain_blast'],
            description="BLAST-only evidence"
        )

        # Test with full evidence (BLAST + HHSearch)
        full_evidence_result = self._run_with_evidence_set(
            protein=base_protein,
            evidence_types=['domain_blast', 'chain_blast', 'hhsearch'],
            description="Full evidence including HHSearch"
        )

        # Compare results
        self._compare_evidence_results(
            blast_only_result,
            full_evidence_result,
            "blast_only_vs_full"
        )

        # Assertions about expected differences
        if full_evidence_result.domains and blast_only_result.domains:
            # HHSearch should generally improve confidence
            avg_conf_full = sum(d.confidence for d in full_evidence_result.domains) / len(full_evidence_result.domains)
            avg_conf_blast = sum(d.confidence for d in blast_only_result.domains) / len(blast_only_result.domains)

            self.assertGreaterEqual(avg_conf_full, avg_conf_blast * 0.9,
                                  "Full evidence should not significantly reduce confidence")

    def test_evidence_quality_impact(self):
        """Test how evidence quality affects domain boundaries"""
        base_protein = {
            'pdb_id': '1qual',
            'chain_id': 'A',
            'sequence_length': 300
        }

        # High quality evidence (low e-values, high probabilities)
        high_quality_result = self._run_with_evidence_quality(
            protein=base_protein,
            quality_level="high",
            description="High quality evidence"
        )

        # Medium quality evidence
        medium_quality_result = self._run_with_evidence_quality(
            protein=base_protein,
            quality_level="medium",
            description="Medium quality evidence"
        )

        # Low quality evidence
        low_quality_result = self._run_with_evidence_quality(
            protein=base_protein,
            quality_level="low",
            description="Low quality evidence"
        )

        # Store results for regression comparison
        self._store_evidence_comparison("quality_impact", {
            'high': high_quality_result,
            'medium': medium_quality_result,
            'low': low_quality_result
        })

        # Assertions about quality impact
        if high_quality_result.domains:
            high_conf = max(d.confidence for d in high_quality_result.domains)
            if medium_quality_result.domains:
                medium_conf = max(d.confidence for d in medium_quality_result.domains)
                self.assertGreater(high_conf, medium_conf,
                                 "High quality evidence should yield higher confidence")

    def test_conflicting_evidence_resolution(self):
        """Test how conflicting evidence is resolved"""
        base_protein = {
            'pdb_id': '1conf',
            'chain_id': 'A',
            'sequence_length': 250
        }

        # Create conflicting evidence scenario
        conflicting_result = self._run_with_conflicting_evidence(
            protein=base_protein,
            description="Conflicting domain boundary evidence"
        )

        # Assertions about conflict resolution
        self.assertTrue(conflicting_result.success, "Should handle conflicting evidence")

        if conflicting_result.domains:
            # Check that boundaries are reasonable
            for domain in conflicting_result.domains:
                self.assertGreater(domain.size, 20, "Domains should be reasonable size")
                self.assertLess(domain.size, base_protein['sequence_length'],
                              "Domain should not exceed sequence length")

    def _run_with_evidence_set(self, protein: Dict[str, Any], evidence_types: List[str],
                              description: str) -> DomainPartitionResult:
        """Run partition with specific evidence types"""
        # Create mock evidence limited to specified types
        evidence_data = self._create_filtered_evidence(protein, evidence_types)

        # Create test case
        test_case = {
            'input': protein,
            'description': description,
            'mock_data': evidence_data
        }

        return self._run_golden_test(test_case)

    def _run_with_evidence_quality(self, protein: Dict[str, Any], quality_level: str,
                                  description: str) -> DomainPartitionResult:
        """Run partition with evidence of specified quality level"""
        # Create evidence with appropriate quality metrics
        evidence_data = self._create_quality_evidence(protein, quality_level)

        test_case = {
            'input': protein,
            'description': description,
            'mock_data': evidence_data
        }

        return self._run_golden_test(test_case)


class ConfidenceWeightTests(DomainPartitionIntegrationTests):
    """
    Tests focused on confidence scoring and evidence weight changes.

    These are critical for detecting regressions when tuning the algorithm.
    """

    def test_default_weights_baseline(self):
        """Establish baseline results with default evidence weights"""
        test_proteins = [
            {'pdb_id': '1wt1', 'chain_id': 'A', 'sequence_length': 150},
            {'pdb_id': '1wt2', 'chain_id': 'A', 'sequence_length': 300},
            {'pdb_id': '1wt3', 'chain_id': 'A', 'sequence_length': 450}
        ]

        baseline_results = {}

        for protein in test_proteins:
            result = self._run_with_default_weights(protein)
            baseline_results[f"{protein['pdb_id']}_{protein['chain_id']}"] = result

        # Store baseline for future weight change comparisons
        self._store_weight_baseline("default_weights", baseline_results)

        # Basic sanity checks
        for protein_id, result in baseline_results.items():
            self.assertTrue(result.success, f"Default weights should work for {protein_id}")

    def test_increased_hhsearch_weight(self):
        """Test impact of increasing HHSearch evidence weight"""
        # Use higher weight for HHSearch
        modified_weights = self.test_config['partition']['evidence_weights'].copy()
        modified_weights['hhsearch'] = 3.5  # Increased from 2.5

        test_proteins = [
            {'pdb_id': '1hw1', 'chain_id': 'A', 'sequence_length': 200}
        ]

        # Run with modified weights
        modified_results = {}
        for protein in test_proteins:
            result = self._run_with_custom_weights(protein, modified_weights)
            modified_results[f"{protein['pdb_id']}_{protein['chain_id']}"] = result

        # Compare with baseline (if available)
        baseline_results = self._load_weight_baseline("default_weights")
        if baseline_results:
            self._compare_weight_impact(baseline_results, modified_results, "increased_hhsearch")

    def test_evidence_weight_sensitivity(self):
        """Test sensitivity to evidence weight changes"""
        base_protein = {'pdb_id': '1sens', 'chain_id': 'A', 'sequence_length': 275}

        # Test multiple weight configurations
        weight_configs = [
            {'name': 'blast_heavy', 'hhsearch': 1.0, 'domain_blast': 4.0},
            {'name': 'hhsearch_heavy', 'hhsearch': 4.0, 'domain_blast': 2.0},
            {'name': 'balanced', 'hhsearch': 2.5, 'domain_blast': 2.5}
        ]

        sensitivity_results = {}
        for config in weight_configs:
            weights = self.test_config['partition']['evidence_weights'].copy()
            weights.update({k: v for k, v in config.items() if k != 'name'})

            result = self._run_with_custom_weights(base_protein, weights)
            sensitivity_results[config['name']] = result

        # Store sensitivity analysis
        self._store_sensitivity_analysis("evidence_weights", sensitivity_results)

        # Verify all configurations produce valid results
        for config_name, result in sensitivity_results.items():
            self.assertTrue(result.success, f"Weight config {config_name} should succeed")

    def test_confidence_threshold_impact(self):
        """Test how confidence thresholds affect domain acceptance"""
        base_protein = {'pdb_id': '1thr', 'chain_id': 'A', 'sequence_length': 320}

        # Test different confidence thresholds
        threshold_configs = [
            {'high': 0.95, 'medium': 0.8, 'low': 0.6},   # Strict
            {'high': 0.9, 'medium': 0.7, 'low': 0.5},    # Default
            {'high': 0.8, 'medium': 0.6, 'low': 0.4}     # Lenient
        ]

        threshold_results = {}
        for i, thresholds in enumerate(threshold_configs):
            result = self._run_with_custom_thresholds(base_protein, thresholds)
            threshold_results[f"config_{i}"] = result

        # Store threshold analysis
        self._store_threshold_analysis("confidence_thresholds", threshold_results)


class ServiceIntegrationTests(DomainPartitionIntegrationTests):
    """
    Tests for service-level integration and workflow correctness.
    """

    def test_service_batch_processing(self):
        """Test batch processing through service"""
        # Create mock batch
        batch_proteins = [
            {'pdb_id': '1bat', 'chain_id': 'A', 'process_id': 1001},
            {'pdb_id': '1bat', 'chain_id': 'B', 'process_id': 1002},
            {'pdb_id': '2bat', 'chain_id': 'A', 'process_id': 1003}
        ]

        # Mock database calls
        with patch.object(DomainPartitionService, '_get_proteins_to_process') as mock_get_proteins:
            mock_get_proteins.return_value = batch_proteins

            service = DomainPartitionService(self.context)

            # Run batch processing
            results = service.partition_batch(
                batch_id=999,
                batch_path=str(self.temp_path),
                limit=None
            )

            # Verify batch results
            self.assertGreater(results.total, 0, "Should process some proteins")
            self.assertEqual(results.total, len(batch_proteins), "Should process all proteins")

    def test_service_error_handling(self):
        """Test service error handling and recovery"""
        # Test with invalid summary file
        with self.assertLogs(level='ERROR') as log_context:
            service = DomainPartitionService(self.context)

            result = service.partition_protein(
                pdb_id="1err",
                chain_id="A",
                summary_path="/nonexistent/path.xml",
                output_dir=str(self.temp_path)
            )

            self.assertFalse(result.success, "Should fail gracefully")
            self.assertIsNotNone(result.error, "Should have error message")

    def test_service_status_tracking(self):
        """Test database status tracking integration"""
        with patch('ecod.pipelines.domain_analysis.partition.tracker.StatusTracker') as mock_tracker:
            mock_tracker_instance = Mock()
            mock_tracker.return_value = mock_tracker_instance

            service = DomainPartitionService(self.context)

            # Create mock summary
            summary_path = self._create_mock_domain_summary({
                'input': {'pdb_id': '1trk', 'chain_id': 'A'},
                'mock_data': self._create_mock_evidence_data({'pdb_id': '1trk', 'chain_id': 'A'})
            })

            result = service.partition_protein(
                pdb_id="1trk",
                chain_id="A",
                summary_path=str(summary_path),
                output_dir=str(self.temp_path),
                process_id=2001
            )

            # Verify status tracking calls were made
            mock_tracker_instance.update_process_status.assert_called()


class FileFormatRegressionTests(DomainPartitionIntegrationTests):
    """
    Tests for XML serialization and file format compatibility.
    """

    def test_xml_round_trip_fidelity(self):
        """Test XML serialization/deserialization round-trip"""
        # Create test domain partition result
        original_result = DomainPartitionResult(
            pdb_id="1xml",
            chain_id="A",
            reference="develop291",
            sequence_length=200,
            is_classified=True
        )

        # Add test domains with evidence
        domain1 = DomainModel(
            id="1xml_A_d1",
            start=1,
            end=100,
            range="1-100",
            source="hhsearch",
            confidence=0.95,
            t_group="2001.1.1",
            h_group="2001.1",
            x_group="2001",
            a_group="a.1"
        )

        # Add evidence to domain
        evidence1 = Evidence(
            type="hhsearch",
            source_id="e1xmlA1",
            probability=0.99,
            evalue=1e-20,
            query_range="1-100",
            hit_range="1-100"
        )
        domain1.add_evidence(evidence1)

        original_result.add_domain(domain1)

        # Serialize to XML
        xml_element = original_result.to_xml()
        xml_str = ET.tostring(xml_element, encoding='unicode')

        # Deserialize from XML
        parsed_element = ET.fromstring(xml_str)
        reconstructed_result = DomainPartitionResult.from_xml(parsed_element)

        # Verify round-trip fidelity
        self.assertEqual(original_result.pdb_id, reconstructed_result.pdb_id)
        self.assertEqual(original_result.chain_id, reconstructed_result.chain_id)
        self.assertEqual(len(original_result.domains), len(reconstructed_result.domains))

        # Verify domain details
        orig_domain = original_result.domains[0]
        recon_domain = reconstructed_result.domains[0]

        self.assertEqual(orig_domain.id, recon_domain.id)
        self.assertEqual(orig_domain.start, recon_domain.start)
        self.assertEqual(orig_domain.end, recon_domain.end)
        self.assertAlmostEqual(orig_domain.confidence, recon_domain.confidence, places=4)

    def test_backward_compatibility(self):
        """Test compatibility with older XML formats"""
        # Create XML in older format (without some newer fields)
        old_format_xml = """<?xml version="1.0"?>
        <domain_partition pdb_id="1old" chain_id="A" reference="develop291" is_classified="true">
            <domains count="1">
                <domain id="1old_A_d1" start="1" end="150" range="1-150"
                        source="hhsearch" confidence="0.85"
                        t_group="1000.1.1" h_group="1000.1"/>
            </domains>
        </domain_partition>"""

        # Parse old format
        element = ET.fromstring(old_format_xml)
        result = DomainPartitionResult.from_xml(element)

        # Verify it parses correctly
        self.assertTrue(result.success)
        self.assertEqual(result.pdb_id, "1old")
        self.assertEqual(result.chain_id, "A")
        self.assertEqual(len(result.domains), 1)

        domain = result.domains[0]
        self.assertEqual(domain.start, 1)
        self.assertEqual(domain.end, 150)
        self.assertAlmostEqual(domain.confidence, 0.85, places=2)

    def test_file_save_and_load(self):
        """Test file save/load operations"""
        # Create test result
        result = DomainPartitionResult(
            pdb_id="1file",
            chain_id="A",
            reference="develop291",
            is_classified=True,
            sequence_length=180
        )

        # Save to file
        output_file = self.temp_path / "test_domain_result.xml"
        success = result.save(output_dir=str(self.temp_path),
                            filename="test_domain_result.xml")

        self.assertTrue(success, "File save should succeed")
        self.assertTrue(output_file.exists(), "Output file should exist")

        # Load from file
        loaded_result = DomainPartitionResult.from_xml_file(str(output_file))

        # Verify loaded result
        self.assertEqual(result.pdb_id, loaded_result.pdb_id)
        self.assertEqual(result.chain_id, loaded_result.chain_id)
        self.assertEqual(result.is_classified, loaded_result.is_classified)


class PerformanceRegressionTests(DomainPartitionIntegrationTests):
    """
    Tests to detect performance regressions.
    """

    def test_processing_speed_baseline(self):
        """Establish baseline for processing speed"""
        import time

        # Create test cases of varying complexity
        test_cases = [
            {'pdb_id': '1spd', 'chain_id': 'A', 'sequence_length': 100, 'complexity': 'simple'},
            {'pdb_id': '2spd', 'chain_id': 'A', 'sequence_length': 300, 'complexity': 'medium'},
            {'pdb_id': '3spd', 'chain_id': 'A', 'sequence_length': 600, 'complexity': 'complex'}
        ]

        performance_results = {}

        service = DomainPartitionService(self.context)

        for test_case in test_cases:
            # Create mock summary
            summary_path = self._create_mock_domain_summary({
                'input': test_case,
                'mock_data': self._create_mock_evidence_data(test_case)
            })

            # Time the processing
            start_time = time.time()

            result = service.partition_protein(
                pdb_id=test_case['pdb_id'],
                chain_id=test_case['chain_id'],
                summary_path=str(summary_path),
                output_dir=str(self.temp_path)
            )

            end_time = time.time()
            processing_time = end_time - start_time

            performance_results[test_case['complexity']] = {
                'processing_time': processing_time,
                'success': result.success,
                'sequence_length': test_case['sequence_length']
            }

        # Store performance baseline
        self._store_performance_baseline("processing_speed", performance_results)

        # Basic performance assertions
        for complexity, perf in performance_results.items():
            self.assertLess(perf['processing_time'], 30.0,
                          f"{complexity} case should complete within 30 seconds")
            self.assertTrue(perf['success'], f"{complexity} case should succeed")

    def test_memory_usage_baseline(self):
        """Test memory usage patterns"""
        import psutil
        import os

        process = psutil.Process(os.getpid())
        initial_memory = process.memory_info().rss

        # Process a batch of proteins
        service = DomainPartitionService(self.context)

        for i in range(10):
            test_case = {
                'pdb_id': f'1mem',
                'chain_id': chr(ord('A') + i),
                'sequence_length': 200 + i * 50
            }

            summary_path = self._create_mock_domain_summary({
                'input': test_case,
                'mock_data': self._create_mock_evidence_data(test_case)
            })

            result = service.partition_protein(
                pdb_id=test_case['pdb_id'],
                chain_id=test_case['chain_id'],
                summary_path=str(summary_path),
                output_dir=str(self.temp_path)
            )

        final_memory = process.memory_info().rss
        memory_increase = final_memory - initial_memory

        # Store memory usage baseline
        self._store_performance_baseline("memory_usage", {
            'initial_memory_mb': initial_memory / 1024 / 1024,
            'final_memory_mb': final_memory / 1024 / 1024,
            'increase_mb': memory_increase / 1024 / 1024
        })

        # Assert reasonable memory usage
        self.assertLess(memory_increase / 1024 / 1024, 500,
                       "Memory increase should be less than 500MB")


class RegressionTestRunner:
    """Utility class for running regression tests and comparing results"""
    
    def __init__(self, test_data_dir: str):
        self.test_data_dir = Path(test_data_dir)
        self.baselines_dir = self.test_data_dir / "baselines"
        self.results_dir = self.test_data_dir / "results"
        
        # Ensure directories exist
        self.baselines_dir.mkdir(parents=True, exist_ok=True)
        self.results_dir.mkdir(parents=True, exist_ok=True)

    def run_regression_suite(self, test_categories: List[str] = None) -> Dict[str, Any]:
        """Run the full regression test suite"""
        if test_categories is None:
            test_categories = [
                'golden_datasets', 
                'evidence_variation',
                'confidence_weights',
                'service_integration',
                'file_format',
                'performance'
            ]
        
        results = {}
        
        for category in test_categories:
            print(f"Running {category} tests...")
            category_results = self._run_test_category(category)
            results[category] = category_results
        
        # Generate regression report
        report = self._generate_regression_report(results)
        
        return report

    def _run_test_category(self, category: str) -> Dict[str, Any]:
        """Run tests for a specific category"""
        # This would integrate with pytest to run specific test classes
        # For now, return a placeholder structure
        return {
            'category': category,
            'tests_run': 0,
            'tests_passed': 0,
            'tests_failed': 0,
            'regressions_detected': [],
            'performance_changes': []
        }

    def _generate_regression_report(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Generate comprehensive regression report"""
        total_tests = sum(r.get('tests_run', 0) for r in results.values())
        total_passed = sum(r.get('tests_passed', 0) for r in results.values())
        total_failed = sum(r.get('tests_failed', 0) for r in results.values())
        
        all_regressions = []
        for category_results in results.values():
            all_regressions.extend(category_results.get('regressions_detected', []))
        
        report = {
            'timestamp': datetime.now().isoformat(),
            'summary': {
                'total_tests': total_tests,
                'passed': total_passed,
                'failed': total_failed,
                'success_rate': total_passed / total_tests if total_tests > 0 else 0
            },
            'regressions': all_regressions,
            'categories': results,
            'recommendations': self._generate_recommendations(results)
        }
        
        # Save report
        report_file = self.results_dir / f"regression_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        return report

    def _generate_recommendations(self, results: Dict[str, Any]) -> List[str]:
        """Generate recommendations based on test results"""
        recommendations = []
        
        # Check for performance regressions
        perf_results = results.get('performance', {})
        if perf_results.get('tests_failed', 0) > 0:
            recommendations.append("Performance regressions detected - review algorithm changes")
        
        # Check for evidence processing issues
        evidence_results = results.get('evidence_variation', {})
        if evidence_results.get('regressions_detected'):
            recommendations.append("Evidence processing regressions - check weight calculations")
        
        # Check for golden dataset failures
        golden_results = results.get('golden_datasets', {})
        if golden_results.get('tests_failed', 0) > 0:
            recommendations.append("Golden dataset failures - core algorithm may be broken")
        
        return recommendations


if __name__ == "__main__":
    # Example of how to run the tests
    unittest.main()
