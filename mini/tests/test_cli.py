#!/usr/bin/env python3
"""
CLI tests for mini_pyecod

Basic smoke tests for the command-line interface.
"""

import pytest
from pathlib import Path
import subprocess
import sys

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))


class TestCLIBasics:
    """Test basic CLI functionality"""
    
    @pytest.mark.unit
    def test_cli_help(self):
        """Test that help works"""
        result = subprocess.run(
            [sys.executable, "pyecod_mini.py", "--help"],
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent
        )
        
        assert result.returncode == 0
        assert "pyECOD Mini" in result.stdout
        assert "--verbose" in result.stdout
        assert "--batch-id" in result.stdout
    
    @pytest.mark.unit
    def test_cli_validate(self):
        """Test configuration validation"""
        result = subprocess.run(
            [sys.executable, "pyecod_mini.py", "--validate"],
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent
        )
        
        # Should run without crashing
        assert result.returncode == 0
        assert "Configuration" in result.stdout
    
    @pytest.mark.unit
    def test_cli_list_batches(self):
        """Test listing batches"""
        result = subprocess.run(
            [sys.executable, "pyecod_mini.py", "--list-batches"],
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent
        )
        
        # Should run without crashing
        assert result.returncode == 0
        # Either shows batches or says none found
        assert "batch" in result.stdout.lower()
    
    @pytest.mark.unit
    def test_cli_missing_protein_id(self):
        """Test error when protein ID is missing"""
        result = subprocess.run(
            [sys.executable, "pyecod_mini.py"],
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent
        )
        
        # Should show error
        assert result.returncode != 0
        assert "protein_id is required" in result.stderr or "protein_id is required" in result.stdout
    
    @pytest.mark.integration
    @pytest.mark.slow
    def test_cli_with_protein(self, stable_batch_dir):
        """Test running with a real protein ID"""
        # Only run if we have batch data
        if not Path(stable_batch_dir).exists():
            pytest.skip("No batch directory available")
        
        result = subprocess.run(
            [sys.executable, "pyecod_mini.py", "8ovp_A", "--batch-id", "036"],
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent,
            timeout=60  # Timeout after 60 seconds
        )
        
        # Check that it ran (might fail if data missing, but shouldn't crash)
        output = result.stdout + result.stderr
        
        # Should either succeed or give meaningful error
        assert "RESULTS:" in output or "ERROR:" in output or "not found" in output


class TestBatchFinder:
    """Test the BatchFinder functionality"""
    
    @pytest.mark.unit
    def test_batch_finder_import(self):
        """Test that BatchFinder can be imported and instantiated"""
        from mini.pyecod_mini import BatchFinder
        
        finder = BatchFinder("/tmp/test")
        assert finder.base_dir == Path("/tmp/test")
        assert isinstance(finder.stable_batches, dict)
    
    @pytest.mark.unit
    def test_config_import(self):
        """Test that PyEcodMiniConfig can be imported"""
        from mini.pyecod_mini import PyEcodMiniConfig
        
        config = PyEcodMiniConfig()
        assert config.base_dir == Path("/data/ecod/pdb_updates/batches")
        assert config.test_data_dir.name == "test_data"


class TestMainFunctions:
    """Test main module functions"""
    
    @pytest.mark.unit
    def test_partition_protein_import(self):
        """Test that main functions can be imported"""
        from mini.pyecod_mini import partition_protein, analyze_protein_batches
        
        # Functions should exist
        assert callable(partition_protein)
        assert callable(analyze_protein_batches)


if __name__ == "__main__":
    # Allow running tests directly
    pytest.main([__file__, "-v"])
