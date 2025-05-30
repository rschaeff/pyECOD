#!/usr/bin/env python3
"""
End-to-End Integration Tests for Domain Partition Pipeline

Tests the complete domain partition workflow using the actual script entry point
to ensure the entire system works as users would experience it.
"""

import os
import sys
import pytest
import tempfile
import subprocess
import yaml
import json
import shutil
import time
from pathlib import Path
from typing import Dict, Any, List, Optional
from unittest.mock import patch

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ecod.db import DBManager
from ecod.db.migration_manager import MigrationManager
from ecod.db.repositories.protein_repository import ProteinRepository
from ecod.models.protein import Protein, ProteinSequence

class TestDomainPartitionEndToEnd:
    """End-to-end tests using the actual domain_partition_run.py script"""
    
    @pytest.fixture(scope="class")
    def test_db_config(self):
        """Test database configuration"""
        return {
            'host': 'localhost',
            'port': 5432,
            'database': 'ecod_test',
            'user': 'test_user',
            'password': 'test_pass'
        }
    
    @pytest.fixture(scope="class")
    def test_db_manager(self, test_db_config):
        """Set up test database with schema"""
        db_manager = DBManager(test_db_config)
        
        # Apply migrations to ensure schema exists
        migration_manager = MigrationManager(test_db_config, "ecod/db/migrations")
        try:
            migration_manager.apply_migrations()
        except Exception as e:
            pytest.skip(f"Could not set up test database: {e}")
        
        yield db_manager
        
        # Cleanup - truncate test data but keep schema
        with db_manager.get_connection() as conn:
            with conn.cursor() as cursor:
                # Truncate tables in correct order (respecting foreign keys)
                tables = [
                    'ecod_schema.job_item',
                    'ecod_schema.job', 
                    'ecod_schema.process_file',
                    'ecod_schema.process_status',
                    'ecod_schema.protein_sequence',
                    'ecod_schema.protein',
                    'ecod_schema.batch',
                ]
                for table in tables:
                    cursor.execute(f"TRUNCATE TABLE {table} CASCADE")
    
    @pytest.fixture
    def test_config_file(self, test_db_config, tmp_path):
        """Create test configuration file"""
        config = {
            'database': test_db_config,
            'reference': {'current_version': 'develop291'},
            'partition': {
                'confidence_thresholds': {'high': 0.9, 'medium': 0.7, 'low': 0.5},
                'evidence_weights': {'domain_blast': 3.0, 'hhsearch': 2.5, 'chain_blast': 2.0},
                'overlap_tolerance': 0.15,
                'min_domain_size': 20,
                'peptide_threshold': 50
            },
            'job_manager': {'type': 'local'},
            'paths': {
                'output_dir': str(tmp_path / "output"),
                'temp_dir': str(tmp_path / "temp")
            }
        }
        
        config_file = tmp_path / "test_config.yml"
        with open(config_file, 'w') as f:
            yaml.dump(config, f)
        
        return str(config_file)
    
    @pytest.fixture
    def script_path(self):
        """Path to the domain partition script"""
        return os.path.join(os.path.dirname(__file__), "..", "scripts", "domain_partition_run.py")
    
    @pytest.fixture
    def test_batch_data(self, test_db_manager, tmp_path):
        """Create test batch with realistic data"""
        batch_dir = tmp_path / "test_batch"
        batch_dir.mkdir()
        
        # Create batch in database
        batch_id = test_db_manager.insert(
            "ecod_schema.batch",
            {
                "batch_name": "integration_test_batch",
                "base_path": str(batch_dir),
                "type": "domain_analysis",
                "ref_version": "develop291",
                "total_items": 3,
                "status": "processing"
            },
            "id"
        )
        
        # Create test proteins and their data
        protein_repo = ProteinRepository(test_db_manager)
        proteins_data = []
        
        test_proteins = [
            {"pdb_id": "1TST", "chain_id": "A", "length": 150, "is_rep": True},
            {"pdb_id": "2TST", "chain_id": "A", "length": 300, "is_rep": False}, 
            {"pdb_id": "3TST", "chain_id": "B", "length": 45, "is_rep": True}  # Peptide
        ]
        
        for prot_data in test_proteins:
            # Create protein
            protein = Protein(
                pdb_id=prot_data["pdb_id"],
                chain_id=prot_data["chain_id"],
                source_id=f"{prot_data['pdb_id']}_{prot_data['chain_id']}",
                length=prot_data["length"]
            )
            
            # Add sequence
            sequence = "M" + "AKVLTKSPG" * (prot_data["length"] // 9)
            sequence = sequence[:prot_data["length"]]
            protein.sequence = ProteinSequence(
                sequence=sequence,
                sequence_md5=f"test_md5_{prot_data['pdb_id']}"
            )
            
            protein_id = protein_repo.create(protein)
            
            # Create process status
            process_id = test_db_manager.insert(
                "ecod_schema.process_status",
                {
                    "protein_id": protein_id,
                    "batch_id": batch_id,
                    "current_stage": "domain_summary_complete",
                    "status": "ready",
                    "is_representative": prot_data["is_rep"]
                },
                "id"
            )
            
            # Create domain summary files
            self._create_test_domain_summary(
                batch_dir, prot_data["pdb_id"], prot_data["chain_id"], 
                prot_data["length"]
            )
            
            # Register summary files in database
            summary_path = f"domains/{prot_data['pdb_id']}_{prot_data['chain_id']}.develop291.domains_v14.xml"
            test_db_manager.insert(
                "ecod_schema.process_file",
                {
                    "process_id": process_id,
                    "file_type": "domain_summary",
                    "file_path": summary_path,
                    "file_exists": True,
                    "file_size": 1024
                }
            )
            
            proteins_data.append({
                "protein_id": protein_id,
                "process_id": process_id,
                **prot_data
            })
        
        return {
            "batch_id": batch_id,
            "batch_dir": str(batch_dir),
            "proteins": proteins_data
        }
    
    def _create_test_domain_summary(self, batch_dir: Path, pdb_id: str, chain_id: str, length: int):
        """Create realistic test domain summary file"""
        domains_dir = batch_dir / "domains"
        domains_dir.mkdir(exist_ok=True)
        
        # Create XML content based on protein length
        if length < 50:  # Peptide
            xml_content = f"""<?xml version="1.0" encoding="utf-8"?>
<blast_summ_doc>
    <blast_summ pdb="{pdb_id}" chain="{chain_id}"/>
    <!-- No hits for peptide -->
</blast_summ_doc>"""
        else:
            # Regular protein with hits
            xml_content = f"""<?xml version="1.0" encoding="utf-8"?>
<blast_summ_doc>
    <blast_summ pdb="{pdb_id}" chain="{chain_id}"/>
    <blast_run program="blastp">
        <hits>
            <hit domain_id="e{pdb_id}{chain_id}1" evalues="1e-25" hsp_count="1">
                <query_reg>1-{length}</query_reg>
                <hit_reg>1-{length}</hit_reg>
            </hit>
        </hits>
    </blast_run>
    <hh_run program="hhsearch">
        <hits>
            <hit domain_id="e{pdb_id}{chain_id}1" probability="95.0" evalue="1e-20" score="85.0">
                <query_reg>1-{length}</query_reg>
                <hit_reg>1-{length}</hit_reg>
            </hit>
        </hits>
    </hh_run>
</blast_summ_doc>"""
        
        summary_file = domains_dir / f"{pdb_id}_{chain_id}.develop291.domains_v14.xml"
        summary_file.write_text(xml_content)
    
    def _run_script(self, script_path: str, config_file: str, args: List[str], 
                   timeout: int = 60) -> subprocess.CompletedProcess:
        """Run the domain partition script with given arguments"""
        cmd = [
            sys.executable, script_path,
            "--config", config_file,
            "--verbose"
        ] + args
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout
        )
        
        return result
    
    def test_single_protein_processing(self, script_path, test_config_file, test_batch_data):
        """Test processing a single protein through the script"""
        batch_id = test_batch_data["batch_id"]
        protein = test_batch_data["proteins"][0]  # Use first protein
        
        # Run single protein processing
        result = self._run_script(
            script_path, test_config_file,
            [
                "process", "single",
                "--pdb-id", protein["pdb_id"],
                "--chain-id", protein["chain_id"], 
                "--batch-id", str(batch_id)
            ]
        )
        
        # Check execution success
        print(f"STDOUT:\n{result.stdout}")
        if result.stderr:
            print(f"STDERR:\n{result.stderr}")
        
        assert result.returncode == 0, f"Script failed with return code {result.returncode}"
        
        # Verify output mentions processing
        assert f"Processing {protein['pdb_id']}_{protein['chain_id']}" in result.stdout
        assert "Successfully processed" in result.stdout
        
        # Check that domain file was created
        expected_domain_file = Path(test_batch_data["batch_dir"]) / "domains" / f"{protein['pdb_id']}_{protein['chain_id']}.develop291.partition_v14.xml"
        assert expected_domain_file.exists(), "Domain partition file should be created"
    
    def test_batch_processing(self, script_path, test_config_file, test_batch_data):
        """Test processing an entire batch through the script"""
        batch_id = test_batch_data["batch_id"]
        
        # Run batch processing
        result = self._run_script(
            script_path, test_config_file,
            [
                "process", "batch",
                "--batch-id", str(batch_id),
                "--force"  # Force even if not 100% ready
            ]
        )
        
        print(f"STDOUT:\n{result.stdout}")
        if result.stderr:
            print(f"STDERR:\n{result.stderr}")
        
        assert result.returncode == 0, f"Batch processing failed with return code {result.returncode}"
        
        # Verify processing summary
        assert "Domain partition complete" in result.stdout
        assert "succeeded" in result.stdout
        
        # Verify classification summary appears
        assert "Classification summary:" in result.stdout
        
        # Check that domain files were created for non-peptide proteins
        domains_dir = Path(test_batch_data["batch_dir"]) / "domains"
        
        non_peptide_proteins = [p for p in test_batch_data["proteins"] if p["length"] >= 50]
        for protein in non_peptide_proteins:
            domain_file = domains_dir / f"{protein['pdb_id']}_{protein['chain_id']}.develop291.partition_v14.xml"
            assert domain_file.exists(), f"Domain file should exist for {protein['pdb_id']}_{protein['chain_id']}"
    
    def test_representatives_only_processing(self, script_path, test_config_file, test_batch_data):
        """Test processing only representative proteins"""
        batch_id = test_batch_data["batch_id"]
        
        # Run representatives-only processing
        result = self._run_script(
            script_path, test_config_file,
            [
                "process", "batch", 
                "--batch-id", str(batch_id),
                "--reps-only",
                "--force"
            ]
        )
        
        print(f"STDOUT:\n{result.stdout}")
        assert result.returncode == 0
        
        # Should process fewer proteins (only representatives)
        rep_count = sum(1 for p in test_batch_data["proteins"] if p["is_rep"])
        
        # Verify only representatives were processed
        assert "Domain partition complete" in result.stdout
        
        # Check that only representative proteins have domain files
        domains_dir = Path(test_batch_data["batch_dir"]) / "domains"
        
        for protein in test_batch_data["proteins"]:
            domain_file = domains_dir / f"{protein['pdb_id']}_{protein['chain_id']}.develop291.partition_v14.xml"
            
            if protein["is_rep"] and protein["length"] >= 50:
                assert domain_file.exists(), f"Domain file should exist for representative {protein['pdb_id']}_{protein['chain_id']}"
            elif not protein["is_rep"]:
                # Non-representative proteins should not be processed in reps-only mode
                pass  # They might or might not have files from previous tests
    
    def test_blast_only_processing(self, script_path, test_config_file, test_batch_data):
        """Test processing with BLAST-only mode"""
        batch_id = test_batch_data["batch_id"]
        
        # Create BLAST-only summary files
        for protein in test_batch_data["proteins"]:
            if protein["length"] >= 50:  # Skip peptides
                self._create_blast_only_summary(
                    Path(test_batch_data["batch_dir"]), 
                    protein["pdb_id"], 
                    protein["chain_id"],
                    protein["length"]
                )
        
        # Run BLAST-only processing
        result = self._run_script(
            script_path, test_config_file,
            [
                "process", "batch",
                "--batch-id", str(batch_id),
                "--blast-only",
                "--force"
            ]
        )
        
        print(f"STDOUT:\n{result.stdout}")
        assert result.returncode == 0
        
        # Verify BLAST-only processing
        assert "Domain partition complete" in result.stdout
    
    def _create_blast_only_summary(self, batch_dir: Path, pdb_id: str, chain_id: str, length: int):
        """Create BLAST-only summary file"""
        domains_dir = batch_dir / "domains"
        
        xml_content = f"""<?xml version="1.0" encoding="utf-8"?>
<blast_summ_doc>
    <blast_summ pdb="{pdb_id}" chain="{chain_id}"/>
    <blast_run program="blastp">
        <hits>
            <hit domain_id="e{pdb_id}{chain_id}1" evalues="1e-25" hsp_count="1">
                <query_reg>1-{length}</query_reg>
                <hit_reg>1-{length}</hit_reg>
            </hit>
        </hits>
    </blast_run>
</blast_summ_doc>"""
        
        blast_file = domains_dir / f"{pdb_id}_{chain_id}.develop291.blast_domains_v14.xml"
        blast_file.write_text(xml_content)
    
    def test_analyze_batch_status(self, script_path, test_config_file, test_batch_data):
        """Test batch status analysis"""
        batch_id = test_batch_data["batch_id"]
        
        # Run status analysis
        result = self._run_script(
            script_path, test_config_file,
            [
                "analyze", "status",
                "--batch-ids", str(batch_id)
            ]
        )
        
        print(f"STDOUT:\n{result.stdout}")
        assert result.returncode == 0
        
        # Verify status output format
        assert "ID" in result.stdout and "Name" in result.stdout
        assert str(batch_id) in result.stdout
        assert "integration_test_batch" in result.stdout
    
    def test_analyze_protein_status(self, script_path, test_config_file, test_batch_data):
        """Test individual protein status analysis"""
        protein = test_batch_data["proteins"][0]
        
        # Run protein analysis
        result = self._run_script(
            script_path, test_config_file,
            [
                "analyze", "protein",
                "--pdb-id", protein["pdb_id"],
                "--chain-id", protein["chain_id"]
            ]
        )
        
        print(f"STDOUT:\n{result.stdout}")
        assert result.returncode == 0
        
        # Verify protein information is displayed
        assert f"{protein['pdb_id']}_{protein['chain_id']}" in result.stdout
        assert "Process ID:" in result.stdout
        assert "Status:" in result.stdout
    
    def test_analyze_domain_counts(self, script_path, test_config_file, test_batch_data):
        """Test domain count analysis after processing"""
        batch_id = test_batch_data["batch_id"]
        
        # First process the batch to generate domain files
        self._run_script(
            script_path, test_config_file,
            [
                "process", "batch",
                "--batch-id", str(batch_id),
                "--force"
            ]
        )
        
        # Then analyze domain counts
        result = self._run_script(
            script_path, test_config_file,
            [
                "analyze", "counts",
                "--batch-ids", str(batch_id),
                "--sample-size", "10"
            ]
        )
        
        print(f"STDOUT:\n{result.stdout}")
        assert result.returncode == 0
        
        # Verify count analysis output
        assert "Sample" in result.stdout
        assert "w/Domains" in result.stdout
        assert str(batch_id) in result.stdout
    
    def test_specific_proteins_processing(self, script_path, test_config_file, test_batch_data):
        """Test processing specific proteins by process ID"""
        # Get process IDs for first two proteins
        process_ids = [str(p["process_id"]) for p in test_batch_data["proteins"][:2]]
        
        # Run specific protein processing
        result = self._run_script(
            script_path, test_config_file,
            [
                "process", "specific",
                "--process-ids"] + process_ids
        )
        
        print(f"STDOUT:\n{result.stdout}")
        assert result.returncode == 0
        
        # Verify processing of specific proteins
        assert f"Processing {len(process_ids)} specific proteins" in result.stdout
        assert "Domain partition complete" in result.stdout
    
    def test_repair_missing_files(self, script_path, test_config_file, test_batch_data, test_db_manager):
        """Test repairing missing domain files"""
        batch_id = test_batch_data["batch_id"]
        
        # First, process batch to create domain files
        self._run_script(
            script_path, test_config_file,
            [
                "process", "batch",
                "--batch-id", str(batch_id),
                "--force"
            ]
        )
        
        # Remove a domain file to simulate missing file
        protein = test_batch_data["proteins"][0]
        domain_file = Path(test_batch_data["batch_dir"]) / "domains" / f"{protein['pdb_id']}_{protein['chain_id']}.develop291.partition_v14.xml"
        
        if domain_file.exists():
            domain_file.unlink()
            
            # Update database to reflect missing file
            test_db_manager.update(
                "ecod_schema.process_file",
                {"file_exists": False},
                "process_id = %s AND file_type = 'domain_partition'",
                (protein["process_id"],)
            )
        
        # Run repair
        result = self._run_script(
            script_path, test_config_file,
            [
                "repair", "missing",
                "--batch-id", str(batch_id),
                "--limit", "5"
            ]
        )
        
        print(f"STDOUT:\n{result.stdout}")
        assert result.returncode == 0
        
        # Verify repair process
        assert "missing domain files" in result.stdout or "No missing domain files found" in result.stdout
    
    def test_error_handling_invalid_batch(self, script_path, test_config_file):
        """Test error handling for invalid batch ID"""
        # Run with non-existent batch ID
        result = self._run_script(
            script_path, test_config_file,
            [
                "process", "batch",
                "--batch-id", "99999"
            ]
        )
        
        print(f"STDOUT:\n{result.stdout}")
        print(f"STDERR:\n{result.stderr}")
        
        # Should fail gracefully
        assert result.returncode != 0
        assert "not found" in result.stderr or "not found" in result.stdout
    
    def test_error_handling_missing_config(self, script_path):
        """Test error handling for missing configuration"""
        # Run with non-existent config file
        result = self._run_script(
            script_path, "/nonexistent/config.yml",
            ["analyze", "status"]
        )
        
        # Should fail due to missing config
        assert result.returncode != 0
    
    @pytest.mark.slow
    def test_monitor_batch_processing(self, script_path, test_config_file, test_batch_data):
        """Test batch monitoring functionality"""
        batch_id = test_batch_data["batch_id"]
        
        # Start monitoring with short timeout
        result = self._run_script(
            script_path, test_config_file,
            [
                "monitor", "batch",
                "--batch-id", str(batch_id),
                "--interval", "5",
                "--timeout", "10"  # Short timeout for testing
            ]
        )
        
        print(f"STDOUT:\n{result.stdout}")
        
        # Should complete within timeout or show monitoring output
        assert result.returncode == 0
        assert f"Monitoring batch {batch_id}" in result.stdout
    
    def test_configuration_validation(self, script_path, tmp_path):
        """Test script behavior with invalid configuration"""
        # Create invalid config (missing required fields)
        invalid_config = {
            'database': {
                'host': 'localhost'
                # Missing other required database fields
            }
        }
        
        config_file = tmp_path / "invalid_config.yml"
        with open(config_file, 'w') as f:
            yaml.dump(invalid_config, f)
        
        # Should handle configuration errors gracefully
        result = self._run_script(
            script_path, str(config_file),
            ["analyze", "status"]
        )
        
        assert result.returncode != 0
        # Should have meaningful error message
        assert "error" in result.stderr.lower() or "error" in result.stdout.lower()


class TestDomainPartitionScriptIntegration:
    """Additional integration tests for script functionality"""
    
    def test_help_and_usage(self, tmp_path):
        """Test that script provides proper help and usage information"""
        script_path = os.path.join(os.path.dirname(__file__), "..", "scripts", "domain_partition_run.py")
        
        # Test main help
        result = subprocess.run(
            [sys.executable, script_path, "--help"],
            capture_output=True,
            text=True
        )
        
        assert result.returncode == 0
        assert "Domain Partition Tool" in result.stdout
        assert "process" in result.stdout
        assert "analyze" in result.stdout
        
        # Test subcommand help
        result = subprocess.run(
            [sys.executable, script_path, "process", "--help"],
            capture_output=True,
            text=True
        )
        
        assert result.returncode == 0
        assert "batch" in result.stdout
        assert "single" in result.stdout
    
    @pytest.mark.parametrize("mode,action", [
        ("process", "batch"),
        ("analyze", "status"), 
        ("repair", "missing"),
        ("monitor", "batch")
    ])
    def test_subcommand_help(self, mode, action):
        """Test help for various subcommands"""
        script_path = os.path.join(os.path.dirname(__file__), "..", "scripts", "domain_partition_run.py")
        
        result = subprocess.run(
            [sys.executable, script_path, mode, action, "--help"],
            capture_output=True,
            text=True
        )
        
        assert result.returncode == 0
        assert "--" in result.stdout  # Should show command-line options


@pytest.mark.performance
class TestDomainPartitionPerformance:
    """Performance tests for the domain partition script"""
    
    def test_batch_processing_performance(self, test_db_manager, tmp_path):
        """Test performance with larger batch sizes"""
        # Create larger test batch
        batch_size = 20
        
        # This would be similar to the other tests but with more data
        # Implementation depends on performance requirements
        pass
    
    def test_memory_usage_monitoring(self):
        """Test memory usage during processing"""
        # Could use psutil to monitor memory usage
        pass


if __name__ == "__main__":
    # Run with: python -m pytest test_domain_partition_e2e.py -v
    pytest.main([__file__, "-v", "--tb=short"])
