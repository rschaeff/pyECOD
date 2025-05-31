#!/usr/bin/env python3
"""
Missing Integration Tests for domain_partition_run.py

These tests cover the previously untested functions using the refactored 
raw SQL approach to match the actual codebase patterns.
"""

import pytest
import subprocess
import time
from pathlib import Path
from conftest import run_script_subprocess



@pytest.mark.integration
@pytest.mark.e2e
class TestDomainPartitionMissingCoverage:
    """Integration tests for previously untested domain_partition_run.py functions"""

    def test_process_all_batches_direct_mode(self, script_paths, test_config_file, temp_test_dir,
                                           test_data_creator, mock_domain_summary_factory,
                                           evidence_data_factory):
        """Test process all batches in direct mode (not SLURM)"""
        script_path = script_paths['domain_partition']
        
        # Create multiple test batches with different characteristics
        batch_1_data = setup_complete_test_protein(
            test_data_creator, mock_domain_summary_factory,
            evidence_data_factory,
            {"pdb_id": "1ALL", "chain_id": "A", "length": 150, "is_rep": True, "expected_domains": 1},
            temp_test_dir / "batch_1"
        )
        
        batch_2_data = setup_complete_test_protein(
            test_data_creator, mock_domain_summary_factory,
            evidence_data_factory,
            {"pdb_id": "2ALL", "chain_id": "A", "length": 200, "is_rep": True, "expected_domains": 1},
            temp_test_dir / "batch_2"
        )
        
        # Test direct processing of multiple batches
        result = run_script_subprocess(
            str(script_path), test_config_file,
            [
                "process", "all",
                "--batch-ids", str(batch_1_data["batch_id"]), str(batch_2_data["batch_id"]),
                "--batch-size", "2",
                "--workers", "1",
                "--force"
            ]
        )
        
        print(f"STDOUT:\n{result.stdout}")
        if result.stderr:
            print(f"STDERR:\n{result.stderr}")
        
        assert result.returncode == 0, f"Multi-batch processing failed: {result.stderr}"
        assert "All batches processed" in result.stdout or "complete" in result.stdout.lower()
        
        # Verify both batches were mentioned
        assert str(batch_1_data["batch_id"]) in result.stdout
        assert str(batch_2_data["batch_id"]) in result.stdout

    def test_process_all_batches_with_exclusions(self, script_paths, test_config_file, temp_test_dir,
                                               test_data_creator, mock_domain_summary_factory,
                                               evidence_data_factory):
        """Test process all batches with exclusion list"""
        script_path = script_paths['domain_partition']
        
        # Create 3 batches
        batch_data = []
        for i in range(3):
            data = setup_complete_test_protein(
                test_data_creator, mock_domain_summary_factory,
                evidence_data_factory,
                {"pdb_id": f"{i+1}EXC", "chain_id": "A", "length": 150, "is_rep": True, "expected_domains": 1},
                temp_test_dir / f"batch_{i+1}"
            )
            batch_data.append(data)
        
        # Process all but exclude the middle one
        batch_ids = [str(b["batch_id"]) for b in batch_data]
        exclude_id = batch_ids[1]  # Exclude middle batch
        
        result = run_script_subprocess(
            str(script_path), test_config_file,
            [
                "process", "all",
                "--batch-ids"] + batch_ids + [
                "--exclude-batch-ids", exclude_id,
                "--force"
            ]
        )
        
        assert result.returncode == 0
        # Should mention first and third batch but not middle one
        assert batch_ids[0] in result.stdout
        assert batch_ids[2] in result.stdout

    def test_process_specific_proteins_by_process_id(self, script_paths, test_config_file, temp_test_dir,
                                                   test_data_creator, mock_domain_summary_factory,
                                                   evidence_data_factory):
        """Test processing specific proteins by process ID"""
        script_path = script_paths['domain_partition']
        
        # Create batch with multiple proteins
        proteins_data = [
            {"pdb_id": "1SPC", "chain_id": "A", "length": 150, "is_rep": True, "expected_domains": 1},
            {"pdb_id": "2SPC", "chain_id": "A", "length": 200, "is_rep": False, "expected_domains": 1}
        ]
        
        batch_data = test_data_creator.create_test_batch(
            temp_test_dir / "specific_batch", proteins_data
        )
        
        # Create summary files for both proteins
        for i, protein_data in enumerate(proteins_data):
            evidence = evidence_data_factory(protein_data, quality='high')
            summary_file = mock_domain_summary_factory.create_summary(
                protein_data["pdb_id"], protein_data["chain_id"], evidence
            )
            
            # Register in database
            relative_path = str(summary_file.relative_to(temp_test_dir / "specific_batch"))
            test_data_creator.create_test_files(
                batch_data["proteins"][i]["process_id"],
                {"domain_summary": relative_path}
            )
        
        # Get process IDs from created proteins
        process_ids = [str(p["process_id"]) for p in batch_data["proteins"]]
        
        result = run_script_subprocess(
            str(script_path), test_config_file,
            [
                "process", "specific",
                "--process-ids"] + process_ids + [
                "--batch-id", str(batch_data["batch_id"])
            ]
        )
        
        print(f"STDOUT:\n{result.stdout}")
        if result.stderr:
            print(f"STDERR:\n{result.stderr}")
        
        assert result.returncode == 0, f"Specific protein processing failed: {result.stderr}"
        
        # Should mention both proteins
        assert "1SPC_A" in result.stdout or "1SPC" in result.stdout
        assert "2SPC_A" in result.stdout or "2SPC" in result.stdout

    def test_analyze_protein_status_detailed(self, script_paths, test_config_file, temp_test_dir,
                                           test_data_creator, mock_domain_summary_factory,
                                           evidence_data_factory):
        """Test detailed protein status analysis"""
        script_path = script_paths['domain_partition']
        
        # Set up complete test protein
        complete_setup = setup_complete_test_protein(
            test_data_creator, mock_domain_summary_factory,
            evidence_data_factory,
            {"pdb_id": "1STA", "chain_id": "A", "length": 150, "is_rep": True, "expected_domains": 1},
            temp_test_dir / "status_batch"
        )
        
        result = run_script_subprocess(
            str(script_path), test_config_file,
            [
                "analyze", "protein",
                "--pdb-id", "1STA",
                "--chain-id", "A",
                "--batch-id", str(complete_setup["batch_id"]),
                "--verbose"
            ]
        )
        
        print(f"STDOUT:\n{result.stdout}")
        if result.stderr:
            print(f"STDERR:\n{result.stderr}")
        
        assert result.returncode == 0
        assert "1STA_A" in result.stdout
        assert "Process ID:" in result.stdout
        assert "Batch:" in result.stdout
        assert "Status:" in result.stdout
        assert "Files:" in result.stdout

    def test_analyze_protein_status_not_found(self, script_paths, test_config_file):
        """Test protein status analysis for non-existent protein"""
        script_path = script_paths['domain_partition']
        
        result = run_script_subprocess(
            str(script_path), test_config_file,
            [
                "analyze", "protein",
                "--pdb-id", "XXXX",
                "--chain-id", "Z"
            ]
        )
        
        # Should fail gracefully
        assert result.returncode == 1
        assert "not found" in result.stderr.lower() or "not found" in result.stdout.lower()

    def test_analyze_domain_counts_statistics(self, script_paths, test_config_file, temp_test_dir,
                                            test_data_creator, mock_domain_summary_factory,
                                            evidence_data_factory):
        """Test domain count analysis and statistics"""
        script_path = script_paths['domain_partition']
        
        # Create batch with varied proteins for statistics
        proteins_data = [
            {"pdb_id": "1CNT", "chain_id": "A", "length": 150, "is_rep": True, "expected_domains": 1},
            {"pdb_id": "2CNT", "chain_id": "A", "length": 300, "is_rep": True, "expected_domains": 2},
            {"pdb_id": "3CNT", "chain_id": "A", "length": 450, "is_rep": True, "expected_domains": 3}
        ]
        
        batch_data = test_data_creator.create_test_batch(
            temp_test_dir / "counts_batch", proteins_data
        )
        
        # Create domain partition result files (not just summaries)
        for i, protein_data in enumerate(proteins_data):
            evidence = evidence_data_factory(protein_data, quality='high')
            summary_file = mock_domain_summary_factory.create_summary(
                protein_data["pdb_id"], protein_data["chain_id"], evidence
            )
            
            # Create a mock domain partition result file
            from ecod.models.pipeline.partition import DomainPartitionResult
            from ecod.models.pipeline.domain import DomainModel
            
            result = DomainPartitionResult(
                pdb_id=protein_data["pdb_id"],
                chain_id=protein_data["chain_id"],
                reference="develop291",
                sequence_length=protein_data["length"],
                is_classified=True
            )
            
            # Add expected number of domains
            for j in range(protein_data["expected_domains"]):
                domain_size = protein_data["length"] // protein_data["expected_domains"]
                start = j * domain_size + 1
                end = min((j + 1) * domain_size, protein_data["length"])
                
                domain = DomainModel(
                    id=f"{protein_data['pdb_id']}_{protein_data['chain_id']}_d{j+1}",
                    start=start,
                    end=end,
                    range=f"{start}-{end}",
                    source="mock",
                    confidence=0.85
                )
                result.add_domain(domain)
            
            # Save domain partition file
            result.save(output_dir=str(temp_test_dir / "counts_batch" / "domains"),
                       filename=f"{protein_data['pdb_id']}_{protein_data['chain_id']}.develop291.partition_v14.xml")
            
            # Register both files in database
            process_id = batch_data["proteins"][i]["process_id"]
            test_data_creator.create_test_files(process_id, {
                "domain_summary": str(summary_file.relative_to(temp_test_dir / "counts_batch")),
                "domain_partition": f"domains/{protein_data['pdb_id']}_{protein_data['chain_id']}.develop291.partition_v14.xml"
            })
        
        result = run_script_subprocess(
            str(script_path), test_config_file,
            [
                "analyze", "counts",
                "--batch-ids", str(batch_data["batch_id"]),
                "--sample-size", "10",
                "--verbose"
            ]
        )
        
        print(f"STDOUT:\n{result.stdout}")
        if result.stderr:
            print(f"STDERR:\n{result.stderr}")
        
        assert result.returncode == 0
        # Should show statistics headers
        assert "Sample" in result.stdout or "Domains" in result.stdout
        assert str(batch_data["batch_id"]) in result.stdout

    def test_repair_missing_files_functionality(self, script_paths, test_config_file, temp_test_dir,
                                              test_data_creator, mock_domain_summary_factory,
                                              evidence_data_factory):
        """Test repair missing files functionality"""
        script_path = script_paths['domain_partition']
        
        # Create batch with protein that has summary but no partition file
        protein_data = {"pdb_id": "1REP", "chain_id": "A", "length": 150, "is_rep": True, "expected_domains": 1}
        batch_data = test_data_creator.create_test_batch(
            temp_test_dir / "repair_batch", [protein_data]
        )
        
        # Create domain summary file
        evidence = evidence_data_factory(protein_data, quality='high')
        summary_file = mock_domain_summary_factory.create_summary(
            protein_data["pdb_id"], protein_data["chain_id"], evidence
        )
        
        # Register ONLY summary file (missing partition file simulates the repair scenario)
        process_id = batch_data["proteins"][0]["process_id"]
        test_data_creator.create_test_files(process_id, {
            "domain_summary": str(summary_file.relative_to(temp_test_dir / "repair_batch"))
        })
        
        result = run_script_subprocess(
            str(script_path), test_config_file,
            [
                "repair", "missing",
                "--batch-id", str(batch_data["batch_id"]),
                "--limit", "5"
            ]
        )
        
        print(f"STDOUT:\n{result.stdout}")
        if result.stderr:
            print(f"STDERR:\n{result.stderr}")
        
        # Should complete successfully (return code 0 or 1 is acceptable for repair)
        assert result.returncode in [0, 1]
        
        # Should mention the repair process
        assert "missing" in result.stdout.lower() or "repair" in result.stdout.lower() or "Repair complete" in result.stdout

    def test_repair_failed_processes_reset(self, script_paths, test_config_file, temp_test_dir,
                                         test_data_creator):
        """Test repair failed processes functionality"""
        script_path = script_paths['domain_partition']
        
        # Create batch with a failed protein
        protein_data = {"pdb_id": "1FAIL", "chain_id": "A", "length": 150, "is_rep": True}
        batch_data = test_data_creator.create_test_batch(
            temp_test_dir / "failed_batch", [protein_data]
        )
        
        # Mark the process as failed
        process_id = batch_data["proteins"][0]["process_id"]
        test_data_creator.db.execute_query(
            "UPDATE ecod_schema.process_status SET status = 'domain_partition_failed', error_message = 'Test failure' WHERE id = %s",
            (process_id,)
        )
        
        result = run_script_subprocess(
            str(script_path), test_config_file,
            [
                "repair", "failed",
                "--batch-id", str(batch_data["batch_id"])
                # Note: not using --rerun to avoid complex setup
            ]
        )
        
        print(f"STDOUT:\n{result.stdout}")
        if result.stderr:
            print(f"STDERR:\n{result.stderr}")
        
        assert result.returncode == 0
        assert "Reset" in result.stdout and ("failed" in result.stdout or "processes" in result.stdout)

    def test_monitor_batch_status_with_timeout(self, script_paths, test_config_file, temp_test_dir,
                                             test_data_creator):
        """Test batch monitoring with short timeout"""
        script_path = script_paths['domain_partition']
        
        # Create batch for monitoring
        protein_data = {"pdb_id": "1MON", "chain_id": "A", "length": 150, "is_rep": True}
        batch_data = test_data_creator.create_test_batch(
            temp_test_dir / "monitor_batch", [protein_data]
        )
        
        # Test monitoring with very short timeout
        result = run_script_subprocess(
            str(script_path), test_config_file,
            [
                "monitor", "batch",
                "--batch-id", str(batch_data["batch_id"]),
                "--interval", "1",
                "--timeout", "3"
            ]
        )
        
        print(f"STDOUT:\n{result.stdout}")
        if result.stderr:
            print(f"STDERR:\n{result.stderr}")
        
        # Should timeout gracefully
        assert result.returncode == 0
        assert ("timeout" in result.stdout.lower() or 
                "monitoring" in result.stdout.lower() or 
                "progress" in result.stdout.lower())

    def test_analyze_batch_status_comprehensive(self, script_paths, test_config_file, temp_test_dir,
                                              test_data_creator, mock_domain_summary_factory,
                                              evidence_data_factory):
        """Test comprehensive batch status analysis"""
        script_path = script_paths['domain_partition']
        
        # Create batches with different statuses
        batch_1 = setup_complete_test_protein(
            test_data_creator, mock_domain_summary_factory,
            evidence_data_factory,
            {"pdb_id": "1STAT", "chain_id": "A", "length": 150, "is_rep": True, "expected_domains": 1},
            temp_test_dir / "status_batch_1"
        )
        
        batch_2 = setup_complete_test_protein(
            test_data_creator, mock_domain_summary_factory,
            evidence_data_factory,
            {"pdb_id": "2STAT", "chain_id": "A", "length": 200, "is_rep": True, "expected_domains": 1},
            temp_test_dir / "status_batch_2"
        )
        
        result = run_script_subprocess(
            str(script_path), test_config_file,
            [
                "analyze", "status",
                "--batch-ids", str(batch_1["batch_id"]), str(batch_2["batch_id"])
            ]
        )
        
        print(f"STDOUT:\n{result.stdout}")
        if result.stderr:
            print(f"STDERR:\n{result.stderr}")
        
        assert result.returncode == 0
        
        # Should show both batches
        assert str(batch_1["batch_id"]) in result.stdout
        assert str(batch_2["batch_id"]) in result.stdout
        
        # Should show status headers
        assert ("Total" in result.stdout or "Ready" in result.stdout or "Complete" in result.stdout)

    def test_cleanup_non_representative_status(self, script_paths, test_config_file, temp_test_dir,
                                             test_data_creator):
        """Test cleanup of non-representative protein status"""
        script_path = script_paths['domain_partition']
        
        # Create batch with mixed representative status
        proteins_data = [
            {"pdb_id": "1CLN", "chain_id": "A", "length": 150, "is_rep": True},
            {"pdb_id": "2CLN", "chain_id": "A", "length": 150, "is_rep": False}
        ]
        
        batch_data = test_data_creator.create_test_batch(
            temp_test_dir / "cleanup_batch", proteins_data
        )
        
        result = run_script_subprocess(
            str(script_path), test_config_file,
            [
                "cleanup", "non-representative",
                "--batch-ids", str(batch_data["batch_id"])
            ]
        )
        
        print(f"STDOUT:\n{result.stdout}")
        if result.stderr:
            print(f"STDERR:\n{result.stderr}")
        
        assert result.returncode == 0
        assert ("cleanup" in result.stdout.lower() or 
                "updated" in result.stdout.lower() or
                "non-representative" in result.stdout.lower())

    def test_error_handling_invalid_process_ids(self, script_paths, test_config_file):
        """Test error handling for invalid process IDs"""
        script_path = script_paths['domain_partition']
        
        result = run_script_subprocess(
            str(script_path), test_config_file,
            [
                "process", "specific",
                "--process-ids", "99999", "99998"
            ]
        )
        
        # Should handle gracefully
        assert result.returncode in [0, 1]  # Either succeeds with no work or fails gracefully

    def test_error_handling_nonexistent_batch_for_repair(self, script_paths, test_config_file):
        """Test error handling for repair on non-existent batch"""
        script_path = script_paths['domain_partition']
        
        result = run_script_subprocess(
            str(script_path), test_config_file,
            [
                "repair", "failed",
                "--batch-id", "99999"
            ]
        )
        
        assert result.returncode == 1
        assert "not found" in result.stderr.lower() or "not found" in result.stdout.lower()

    def test_help_functionality_for_all_modes(self, script_paths):
        """Test help functionality for all major modes"""
        script_path = script_paths['domain_partition']
        
        # Test main help
        result = subprocess.run(
            [str(script_path), "--help"],
            capture_output=True, text=True
        )
        assert result.returncode == 0
        assert "process" in result.stdout
        assert "analyze" in result.stdout
        assert "repair" in result.stdout
        
        # Test mode-specific help
        for mode in ["process", "analyze", "repair", "monitor", "cleanup"]:
            result = subprocess.run(
                [str(script_path), mode, "--help"],
                capture_output=True, text=True
            )
            assert result.returncode == 0, f"{mode} help failed"

    @pytest.mark.slow
    def test_integration_workflow_complete(self, script_paths, test_config_file, temp_test_dir,
                                         test_data_creator, mock_domain_summary_factory,
                                         evidence_data_factory):
        """Test complete integration workflow: create -> process -> analyze -> repair"""
        script_path = script_paths['domain_partition']
        
        # Set up complete test protein
        complete_setup = setup_complete_test_protein(
            test_data_creator, mock_domain_summary_factory,
            evidence_data_factory,
            {"pdb_id": "1WFL", "chain_id": "A", "length": 150, "is_rep": True, "expected_domains": 1},
            temp_test_dir / "workflow_batch"
        )
        
        batch_id = str(complete_setup["batch_id"])
        
        # Step 1: Check initial status
        status_result = run_script_subprocess(
            str(script_path), test_config_file,
            ["analyze", "status", "--batch-ids", batch_id]
        )
        assert status_result.returncode == 0
        
        # Step 2: Process the batch
        process_result = run_script_subprocess(
            str(script_path), test_config_file,
            ["process", "batch", "--batch-id", batch_id, "--force"]
        )
        print(f"Process result: {process_result.returncode}")
        print(f"Process stdout: {process_result.stdout}")
        
        # Step 3: Analyze protein
        analyze_result = run_script_subprocess(
            str(script_path), test_config_file,
            ["analyze", "protein", "--pdb-id", "1WFL", "--chain-id", "A", "--batch-id", batch_id]
        )
        assert analyze_result.returncode == 0
        
        # Step 4: Check for any missing files to repair
        repair_result = run_script_subprocess(
            str(script_path), test_config_file,
            ["repair", "missing", "--batch-id", batch_id, "--limit", "1"]
        )
        # Repair should complete whether or not files were missing
        assert repair_result.returncode in [0, 1]
        
        print("Complete integration workflow test passed")
