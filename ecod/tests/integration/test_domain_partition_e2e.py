#!/usr/bin/env python3
"""
End-to-End Integration Tests for Domain Partition Pipeline

Updated to work with the synthesized conftest.py that combines robust database setup,
golden datasets, baseline testing, and performance monitoring.
"""

import os
import sys
import pytest
import subprocess
import json
import time
from pathlib import Path
from typing import Dict, Any, List, Optional

# Import utilities from conftest
from conftest import (
    run_script_subprocess,
    temporary_environment,
    PROJECT_ROOT
)

# Add project root to path for local imports
sys.path.insert(0, str(PROJECT_ROOT))

from ecod.models.pipeline.partition import DomainPartitionResult


@pytest.mark.integration
@pytest.mark.e2e
@pytest.mark.database
class TestDomainPartitionEndToEnd:
    """End-to-end tests using the actual domain_partition_run.py script with enhanced fixtures"""

    def test_script_help_functionality(self, script_paths):
        """Test that the script provides proper help information"""
        script_path = script_paths['domain_partition']

        # Test main help
        result = subprocess.run(
            [sys.executable, str(script_path), "--help"],
            capture_output=True,
            text=True,
            cwd=str(PROJECT_ROOT)
        )

        assert result.returncode == 0
        assert "Domain Partition Tool" in result.stdout
        assert "process" in result.stdout
        assert "analyze" in result.stdout

        # Test subcommand help
        result = subprocess.run(
            [sys.executable, str(script_path), "process", "--help"],
            capture_output=True,
            text=True,
            cwd=str(PROJECT_ROOT)
        )

        assert result.returncode == 0
        assert "batch" in result.stdout
        assert "single" in result.stdout

    @pytest.mark.skip(reason="ProteinRepository not used in production - architectural stub")
    def test_single_protein_processing_with_evidence(self, script_paths, test_config_file,
                                                   temp_test_dir, test_data_creator,
                                                   evidence_data_factory,
                                                   mock_domain_summary_factory):
        """Test single protein processing with realistic evidence data"""
        script_path = script_paths['domain_partition']

        # Create test batch with realistic protein
        protein_data = {
            "pdb_id": "1TST",
            "chain_id": "A",
            "length": 150,
            "is_rep": True,
            "expected_domains": 1
        }

        batch_data = test_data_creator.create_test_batch(
            temp_test_dir / "test_batch",
            [protein_data]
        )

        # Generate evidence and domain summary
        evidence = evidence_data_factory(protein_data, quality='high')
        summary_file = mock_domain_summary_factory(
            protein_data["pdb_id"],
            protein_data["chain_id"],
            evidence
        )

        # Register summary file in database
        test_data_creator.db.insert(
            "ecod_schema.process_file",
            {
                "process_id": batch_data["proteins"][0]["process_id"],
                "file_type": "domain_summary",
                "file_path": str(summary_file.relative_to(temp_test_dir / "test_batch")),
                "file_exists": True,
                "file_size": summary_file.stat().st_size
            }
        )

        # Run single protein processing
        result = run_script_subprocess(
            str(script_path), test_config_file,
            [
                "process", "single",
                "--pdb-id", protein_data["pdb_id"],
                "--chain-id", protein_data["chain_id"],
                "--batch-id", str(batch_data["batch_id"])
            ]
        )

        print(f"STDOUT:\n{result.stdout}")
        if result.stderr:
            print(f"STDERR:\n{result.stderr}")

        if result.returncode != 0:
            pytest.fail(f"Script failed with return code {result.returncode}\nSTDERR: {result.stderr}")

        # Verify processing success
        assert "Successfully processed" in result.stdout or "complete" in result.stdout.lower()

        # Check domain file creation
        expected_domain_file = (temp_test_dir / "test_batch" / "domains" /
                              f"{protein_data['pdb_id']}_{protein_data['chain_id']}.develop291.partition_v14.xml")
        assert expected_domain_file.exists(), f"Domain partition file should be created at {expected_domain_file}"

        # Validate domain file content
        result = DomainPartitionResult.from_xml_file(str(expected_domain_file))
        assert result.pdb_id == protein_data["pdb_id"]
        assert result.chain_id == protein_data["chain_id"]
        assert result.is_classified or result.is_unclassified

    @pytest.mark.skip(reason="ProteinRepository not used in production - architectural stub")
    def test_batch_processing_with_mixed_proteins(self, script_paths, test_config_file,
                                                 temp_test_dir, test_data_creator,
                                                 evidence_data_factory,
                                                 mock_domain_summary_factory):
        """Test batch processing with mixed protein types (regular, multi-domain, peptide)"""
        script_path = script_paths['domain_partition']

        # Create mixed protein dataset
        mixed_proteins = [
            {"pdb_id": "1TST", "chain_id": "A", "length": 150, "is_rep": True, "expected_domains": 1},
            {"pdb_id": "2TST", "chain_id": "A", "length": 400, "is_rep": False, "expected_domains": 2},
            {"pdb_id": "3TST", "chain_id": "B", "length": 35, "is_rep": True, "expected_domains": 0}  # Peptide
        ]

        batch_data = test_data_creator.create_test_batch(
            temp_test_dir / "mixed_batch",
            mixed_proteins
        )

        # Generate evidence and summaries for each protein
        for i, protein_data in enumerate(mixed_proteins):
            evidence = evidence_data_factory(protein_data, quality='medium')
            summary_file = mock_domain_summary_factory(
                protein_data["pdb_id"],
                protein_data["chain_id"],
                evidence
            )

            # Register in database
            test_data_creator.db.insert(
                "ecod_schema.process_file",
                {
                    "process_id": batch_data["proteins"][i]["process_id"],
                    "file_type": "domain_summary",
                    "file_path": str(summary_file.relative_to(temp_test_dir / "mixed_batch")),
                    "file_exists": True,
                    "file_size": summary_file.stat().st_size
                }
            )

        # Run batch processing
        result = run_script_subprocess(
            str(script_path), test_config_file,
            [
                "process", "batch",
                "--batch-id", str(batch_data["batch_id"]),
                "--force"
            ]
        )

        print(f"STDOUT:\n{result.stdout}")
        if result.stderr:
            print(f"STDERR:\n{result.stderr}")

        if result.returncode != 0:
            pytest.fail(f"Batch processing failed: {result.stderr}")

        # Verify processing summary
        assert "Domain partition complete" in result.stdout or "complete" in result.stdout.lower()

        # Check files for non-peptide proteins
        domains_dir = temp_test_dir / "mixed_batch" / "domains"

        for protein_data in mixed_proteins:
            domain_file = domains_dir / f"{protein_data['pdb_id']}_{protein_data['chain_id']}.develop291.partition_v14.xml"

            if protein_data["length"] >= 50:  # Non-peptide
                assert domain_file.exists(), f"Domain file should exist for {protein_data['pdb_id']}_{protein_data['chain_id']}"

                # Verify content
                result = DomainPartitionResult.from_xml_file(str(domain_file))
                assert result.pdb_id == protein_data["pdb_id"]
                assert result.chain_id == protein_data["chain_id"]

    @pytest.mark.skip(reason="ProteinRepository not used in production - architectural stub")
    def test_representatives_only_processing(self, script_paths, test_config_file,
                                           temp_test_dir, test_data_creator,
                                           evidence_data_factory,
                                           mock_domain_summary_factory):
        """Test processing only representative proteins"""
        script_path = script_paths['domain_partition']

        # Create proteins with mixed representative status
        proteins = [
            {"pdb_id": "1REP", "chain_id": "A", "length": 200, "is_rep": True, "expected_domains": 1},
            {"pdb_id": "2NON", "chain_id": "A", "length": 180, "is_rep": False, "expected_domains": 1},
            {"pdb_id": "3REP", "chain_id": "B", "length": 250, "is_rep": True, "expected_domains": 1}
        ]

        batch_data = test_data_creator.create_test_batch(
            temp_test_dir / "rep_batch",
            proteins
        )

        # Generate summaries for all proteins
        for i, protein_data in enumerate(proteins):
            evidence = evidence_data_factory(protein_data, quality='high')
            summary_file = mock_domain_summary_factory(
                protein_data["pdb_id"],
                protein_data["chain_id"],
                evidence
            )

            test_data_creator.db.insert(
                "ecod_schema.process_file",
                {
                    "process_id": batch_data["proteins"][i]["process_id"],
                    "file_type": "domain_summary",
                    "file_path": str(summary_file.relative_to(temp_test_dir / "rep_batch")),
                    "file_exists": True,
                    "file_size": summary_file.stat().st_size
                }
            )

        # Run representatives-only processing
        result = run_script_subprocess(
            str(script_path), test_config_file,
            [
                "process", "batch",
                "--batch-id", str(batch_data["batch_id"]),
                "--reps-only",
                "--force"
            ]
        )

        print(f"STDOUT:\n{result.stdout}")
        if result.returncode != 0:
            pytest.fail(f"Representatives-only processing failed: {result.stderr}")

        # Verify only representatives were processed
        domains_dir = temp_test_dir / "rep_batch" / "domains"

        for protein_data in proteins:
            domain_file = domains_dir / f"{protein_data['pdb_id']}_{protein_data['chain_id']}.develop291.partition_v14.xml"

            if protein_data["is_rep"]:
                assert domain_file.exists(), f"Domain file should exist for representative {protein_data['pdb_id']}_{protein_data['chain_id']}"

    @pytest.mark.skip(reason="ProteinRepository not used in production - architectural stub")
    @pytest.mark.golden
    def test_golden_dataset_processing(self, script_paths, test_config_file, golden_datasets,
                                     temp_test_dir, test_data_creator,
                                     evidence_data_factory, mock_domain_summary_factory):
        """Test processing against golden datasets if available"""
        if not golden_datasets:
            pytest.skip("No golden datasets available")

        script_path = script_paths['domain_partition']

        # Use single_domain dataset if available
        if 'single_domain' not in golden_datasets:
            pytest.skip("single_domain golden dataset not available")

        golden_proteins = golden_datasets['single_domain']

        # Create batch from golden dataset
        batch_data = test_data_creator.create_test_batch(
            temp_test_dir / "golden_batch",
            golden_proteins
        )

        # Generate summaries based on expected results
        for i, protein_data in enumerate(golden_proteins):
            evidence = evidence_data_factory(protein_data, quality='high')
            summary_file = mock_domain_summary_factory(
                protein_data["pdb_id"],
                protein_data["chain_id"],
                evidence
            )

            test_data_creator.db.insert(
                "ecod_schema.process_file",
                {
                    "process_id": batch_data["proteins"][i]["process_id"],
                    "file_type": "domain_summary",
                    "file_path": str(summary_file.relative_to(temp_test_dir / "golden_batch")),
                    "file_exists": True,
                    "file_size": summary_file.stat().st_size
                }
            )

        # Process golden dataset
        result = run_script_subprocess(
            str(script_path), test_config_file,
            [
                "process", "batch",
                "--batch-id", str(batch_data["batch_id"]),
                "--force"
            ]
        )

        if result.returncode != 0:
            pytest.fail(f"Golden dataset processing failed: {result.stderr}")

        # Verify all expected files were created
        domains_dir = temp_test_dir / "golden_batch" / "domains"

        created_files = []
        for protein_data in golden_proteins:
            if protein_data.get("length", 100) >= 50:  # Non-peptide
                domain_file = domains_dir / f"{protein_data['pdb_id']}_{protein_data['chain_id']}.develop291.partition_v14.xml"
                if domain_file.exists():
                    created_files.append(domain_file)

        assert len(created_files) > 0, "No domain files were created for golden dataset"

    @pytest.mark.skip(reason="ProteinRepository not used in production - architectural stub")
    @pytest.mark.performance
    def test_batch_processing_performance(self, script_paths, test_config_file,
                                        temp_test_dir, test_data_creator,
                                        evidence_data_factory, mock_domain_summary_factory,
                                        performance_monitor):
        """Test batch processing performance with monitoring"""
        script_path = script_paths['domain_partition']

        # Create larger dataset for performance testing
        batch_size = 15
        perf_proteins = []

        for i in range(batch_size):
            perf_proteins.append({
                "pdb_id": f"P{i:03d}",
                "chain_id": "A",
                "length": 200 + i * 20,
                "is_rep": i % 4 == 0,  # Every 4th is representative
                "expected_domains": 1 + (i % 3)  # 1-3 domains
            })

        batch_data = test_data_creator.create_test_batch(
            temp_test_dir / "perf_batch",
            perf_proteins
        )

        # Generate summaries
        for i, protein_data in enumerate(perf_proteins):
            evidence = evidence_data_factory(protein_data, quality='medium')
            summary_file = mock_domain_summary_factory(
                protein_data["pdb_id"],
                protein_data["chain_id"],
                evidence
            )

            test_data_creator.db.insert(
                "ecod_schema.process_file",
                {
                    "process_id": batch_data["proteins"][i]["process_id"],
                    "file_type": "domain_summary",
                    "file_path": str(summary_file.relative_to(temp_test_dir / "perf_batch")),
                    "file_exists": True,
                    "file_size": summary_file.stat().st_size
                }
            )

        # Monitor performance
        test_name = "batch_processing_performance"
        performance_monitor.start_monitoring(test_name)

        # Run batch processing
        result = run_script_subprocess(
            str(script_path), test_config_file,
            [
                "process", "batch",
                "--batch-id", str(batch_data["batch_id"]),
                "--force",
                "--workers", "2"
            ],
            timeout=180
        )

        performance_monitor.stop_monitoring(test_name)

        if result.returncode != 0:
            pytest.fail(f"Performance test failed: {result.stderr}")

        # Analyze performance
        perf_results = performance_monitor.get_results()
        test_results = perf_results[test_name]

        # Calculate proteins per second
        if test_results['duration'] > 0:
            test_results['proteins_per_second'] = batch_size / test_results['duration']

        print(f"Performance Results:")
        print(f"  Processed {batch_size} proteins in {test_results['duration']:.2f} seconds")
        print(f"  Rate: {test_results.get('proteins_per_second', 0):.2f} proteins/second")
        print(f"  Memory usage: {test_results['memory_delta_mb']:.1f} MB")

        # Performance assertions (adjust thresholds as needed)
        assert test_results['duration'] < 120, f"Processing took too long: {test_results['duration']}s"
        assert test_results.get('proteins_per_second', 0) > 0.1, "Processing rate too slow"

        # Save performance results for tracking
        perf_file = temp_test_dir / "performance_results.json"
        performance_monitor.save_results(perf_file)

    @pytest.mark.skip(reason="ProteinRepository not used in production - architectural stub")
    def test_analyze_batch_status(self, script_paths, test_config_file, temp_test_dir, test_data_creator):
        """Test batch status analysis functionality"""
        script_path = script_paths['domain_partition']

        # Create simple batch for status testing
        proteins = [{"pdb_id": "1STA", "chain_id": "A", "length": 150, "is_rep": True}]
        batch_data = test_data_creator.create_test_batch(temp_test_dir / "status_batch", proteins)

        # Run status analysis
        result = run_script_subprocess(
            str(script_path), test_config_file,
            [
                "analyze", "status",
                "--batch-ids", str(batch_data["batch_id"])
            ]
        )

        print(f"STDOUT:\n{result.stdout}")

        if result.returncode != 0:
            pytest.fail(f"Status analysis failed: {result.stderr}")

        # Verify status output
        assert str(batch_data["batch_id"]) in result.stdout
        assert "integration_test_batch" in result.stdout

    def test_error_handling_invalid_batch(self, script_paths, test_config_file):
        """Test error handling for invalid batch ID"""
        script_path = script_paths['domain_partition']

        # Run with non-existent batch ID
        result = run_script_subprocess(
            str(script_path), test_config_file,
            [
                "process", "batch",
                "--batch-id", "99999"
            ]
        )

        print(f"STDOUT:\n{result.stdout}")
        print(f"STDERR:\n{result.stderr}")

        # Should fail gracefully
        assert result.returncode != 0
        assert "not found" in result.stderr.lower() or "not found" in result.stdout.lower()

    def test_error_handling_missing_config(self, script_paths):
        """Test error handling for missing configuration"""
        script_path = script_paths['domain_partition']

        # Run with non-existent config file
        result = run_script_subprocess(
            str(script_path), "/nonexistent/config.yml",
            ["analyze", "status"]
        )

        # Should fail due to missing config
        assert result.returncode != 0

    @pytest.mark.baseline
    def test_baseline_comparison(self, script_paths, test_config_file, baseline_results_dir,
                                temp_test_dir, test_data_creator, evidence_data_factory,
                                mock_domain_summary_factory, request):
        """Test baseline comparison functionality if baseline exists"""
        if request.config.getoption("--establish-baseline"):
            pytest.skip("Establishing baseline, not comparing")

        script_path = script_paths['domain_partition']

        # Create test batch with known proteins
        baseline_proteins = [
            {"pdb_id": "1BAS", "chain_id": "A", "length": 200, "is_rep": True, "expected_domains": 1},
            {"pdb_id": "2BAS", "chain_id": "A", "length": 350, "is_rep": True, "expected_domains": 2}
        ]

        batch_data = test_data_creator.create_test_batch(
            temp_test_dir / "baseline_batch",
            baseline_proteins
        )

        # Generate consistent evidence for baseline comparison
        for i, protein_data in enumerate(baseline_proteins):
            evidence = evidence_data_factory(protein_data, quality='high')
            summary_file = mock_domain_summary_factory(
                protein_data["pdb_id"],
                protein_data["chain_id"],
                evidence
            )

            test_data_creator.db.insert(
                "ecod_schema.process_file",
                {
                    "process_id": batch_data["proteins"][i]["process_id"],
                    "file_type": "domain_summary",
                    "file_path": str(summary_file.relative_to(temp_test_dir / "baseline_batch")),
                    "file_exists": True,
                    "file_size": summary_file.stat().st_size
                }
            )

        # Process batch
        result = run_script_subprocess(
            str(script_path), test_config_file,
            [
                "process", "batch",
                "--batch-id", str(batch_data["batch_id"]),
                "--force"
            ]
        )

        if result.returncode != 0:
            pytest.fail(f"Baseline test processing failed: {result.stderr}")

        # If establishing baseline, save results
        if request.config.getoption("--establish-baseline"):
            results_file = baseline_results_dir / "baseline_results.json"
            baseline_data = {
                "proteins_processed": len(baseline_proteins),
                "timestamp": time.time(),
                "config": request.config.getoption("--config")
            }

            with open(results_file, 'w') as f:
                json.dump(baseline_data, f, indent=2)

        # Otherwise compare against baseline
        else:
            baseline_file = baseline_results_dir / "baseline_results.json"
            if baseline_file.exists():
                with open(baseline_file, 'r') as f:
                    baseline_data = json.load(f)

                # Simple comparison - could be enhanced
                assert baseline_data["proteins_processed"] == len(baseline_proteins)


@pytest.mark.integration
class TestDomainPartitionConfiguration:
    """Test configuration handling and validation with enhanced fixtures"""

    def test_configuration_loading(self, test_config, test_database_config):
        """Test that configuration loads correctly"""
        assert 'database' in test_config
        assert 'reference' in test_config
        assert 'partition' in test_config

        # Verify database config has required fields
        required_db_fields = ['host', 'port', 'database', 'user', 'password']
        for field in required_db_fields:
            assert field in test_database_config

    def test_service_configuration_loading(self, test_config):
        """Test that service configuration is present"""
        assert 'services' in test_config
        assert 'domain_partition' in test_config['services']

        service_config = test_config['services']['domain_partition']
        assert 'max_workers' in service_config
        assert 'track_status' in service_config


# Utility function for running test suites
def run_integration_test_suite():
    """Run the complete integration test suite"""
    return pytest.main([
        __file__,
        "-v",
        "--tb=short",
        "-m", "integration and not slow"
    ])


if __name__ == "__main__":
    # Run tests when script is executed directly
    pytest.main([__file__, "-v", "--tb=short"])
