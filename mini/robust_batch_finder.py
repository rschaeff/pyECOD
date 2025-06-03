#!/usr/bin/env python3
"""
Robust batch detection that handles proteins in multiple batches

Key improvements:
1. Warns when proteins exist in multiple batches
2. Has configurable batch preferences
3. Allows pinning specific proteins to specific batches
4. Provides detailed information about batch differences
"""

import os
from pathlib import Path
from typing import Optional, List, Tuple, Dict
import json

class RobustBatchFinder:
    """Enhanced batch finder that handles multi-batch proteins properly"""
    
    def __init__(self, base_dir: str):
        self.base_dir = Path(base_dir)
        self._batch_cache = {}
        
        # Configurable preferences
        self.batch_preferences = self._load_batch_preferences()
        self.protein_pins = self._load_protein_pins()
    
    def find_batch_for_protein(self, protein_id: str, verbose: bool = False) -> Optional[str]:
        """
        Find the best batch for a protein with multi-batch awareness.
        
        Strategy:
        1. Check if protein is pinned to a specific batch
        2. Find all batches containing the protein
        3. Apply batch preferences to select the best one
        4. Warn if multiple batches found
        """
        # Check for explicit pin first
        if protein_id in self.protein_pins:
            pinned_batch = self.protein_pins[protein_id]
            if self._protein_exists_in_batch(protein_id, pinned_batch):
                if verbose:
                    print(f"Using pinned batch for {protein_id}: {pinned_batch}")
                return pinned_batch
            else:
                if verbose:
                    print(f"WARNING: Pinned batch {pinned_batch} doesn't contain {protein_id}")
        
        # Find all batches containing this protein
        available_batches = self._get_available_batches()
        found_batches = []
        
        if verbose:
            print(f"Searching for {protein_id} across {len(available_batches)} batches...")
        
        for batch_name in available_batches:
            if self._protein_exists_in_batch(protein_id, batch_name):
                found_batches.append(batch_name)
                if verbose:
                    print(f"  ✓ Found in {batch_name}")
        
        if not found_batches:
            if verbose:
                print(f"  ✗ {protein_id} not found in any batch")
            return None
        
        if len(found_batches) == 1:
            return found_batches[0]
        
        # Multiple batches found - apply selection logic
        if verbose:
            print(f"  ⚠️  Found in {len(found_batches)} batches: {found_batches}")
        
        selected_batch = self._select_best_batch(found_batches, protein_id, verbose)
        
        if verbose:
            print(f"  → Selected: {selected_batch}")
            print(f"  💡 To pin this choice: echo '{protein_id}: {selected_batch}' >> mini/.protein_pins")
        
        return selected_batch
    
    def _select_best_batch(self, found_batches: List[str], protein_id: str, verbose: bool = False) -> str:
        """
        Select the best batch from multiple options using preferences.
        
        Selection criteria (in order):
        1. Explicit batch preferences
        2. Known stable batches for test cases
        3. Most recent batch
        """
        # Apply explicit preferences
        for preferred_batch in self.batch_preferences.get('preferred_order', []):
            if preferred_batch in found_batches:
                if verbose:
                    print(f"    Using preferred batch: {preferred_batch}")
                return preferred_batch
        
        # Check for known test case batches
        test_case_batches = self.batch_preferences.get('test_case_batches', {})
        if protein_id in test_case_batches:
            test_batch = test_case_batches[protein_id]
            if test_batch in found_batches:
                if verbose:
                    print(f"    Using test case batch: {test_batch}")
                return test_batch
        
        # Fall back to most recent (with warning)
        most_recent = found_batches[-1]  # Assumes sorted order
        if verbose:
            print(f"    Defaulting to most recent: {most_recent}")
            print(f"    ⚠️  This may give inconsistent results across runs")
        
        return most_recent
    
    def analyze_multi_batch_protein(self, protein_id: str) -> Dict[str, any]:
        """
        Analyze a protein that exists in multiple batches to help user decide.
        
        Returns information about differences between batches.
        """
        available_batches = self._get_available_batches()
        found_batches = []
        
        for batch_name in available_batches:
            if self._protein_exists_in_batch(protein_id, batch_name):
                found_batches.append(batch_name)
        
        if len(found_batches) <= 1:
            return {"multi_batch": False, "batches": found_batches}
        
        # Analyze differences between batches
        batch_analysis = {
            "multi_batch": True,
            "batches": found_batches,
            "analysis": {}
        }
        
        for batch_name in found_batches:
            batch_info = self._analyze_protein_in_batch(protein_id, batch_name)
            batch_analysis["analysis"][batch_name] = batch_info
        
        return batch_analysis
    
    def _analyze_protein_in_batch(self, protein_id: str, batch_name: str) -> Dict[str, any]:
        """Analyze a protein's data in a specific batch"""
        batch_dir = self.base_dir / batch_name
        
        info = {
            "has_domain_summary": False,
            "has_domains": False,
            "has_blast": False,
            "file_sizes": {},
            "file_dates": {}
        }
        
        # Check domain summary
        domain_summary_file = batch_dir / "domains" / f"{protein_id}.develop291.domain_summary.xml"
        if domain_summary_file.exists():
            info["has_domain_summary"] = True
            info["file_sizes"]["domain_summary"] = domain_summary_file.stat().st_size
            info["file_dates"]["domain_summary"] = domain_summary_file.stat().st_mtime
        
        # Check domains file
        domains_file = batch_dir / "domains" / f"{protein_id}.develop291.domains.xml"
        if domains_file.exists():
            info["has_domains"] = True
            info["file_sizes"]["domains"] = domains_file.stat().st_size
            info["file_dates"]["domains"] = domains_file.stat().st_mtime
        
        # Check BLAST
        blast_file = batch_dir / "blast" / "chain" / f"{protein_id}.develop291.xml"
        if blast_file.exists():
            info["has_blast"] = True
            info["file_sizes"]["blast"] = blast_file.stat().st_size
            info["file_dates"]["blast"] = blast_file.stat().st_mtime
        
        return info
    
    def _load_batch_preferences(self) -> Dict[str, any]:
        """Load batch preferences from config file"""
        config_file = Path(__file__).parent / ".batch_preferences"
        
        default_preferences = {
            "preferred_order": [
                # Add specific batch names in preference order
                # e.g., "ecod_batch_036_20250406_1424"  # Known stable batch
            ],
            "test_case_batches": {
                # Pin specific test cases to specific batches for reproducibility
                "8ovp_A": "ecod_batch_036_20250406_1424",  # Our validated test case
            },
            "avoid_batches": [
                # List any batches known to have issues
            ]
        }
        
        if config_file.exists():
            try:
                with open(config_file, 'r') as f:
                    user_preferences = json.load(f)
                # Merge with defaults
                for key, value in user_preferences.items():
                    if key in default_preferences:
                        if isinstance(value, dict):
                            default_preferences[key].update(value)
                        else:
                            default_preferences[key] = value
            except Exception as e:
                print(f"Warning: Error loading batch preferences: {e}")
        
        return default_preferences
    
    def _load_protein_pins(self) -> Dict[str, str]:
        """Load protein->batch pins from config file"""
        pins_file = Path(__file__).parent / ".protein_pins"
        pins = {}
        
        if pins_file.exists():
            try:
                with open(pins_file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if line and ':' in line and not line.startswith('#'):
                            protein, batch = line.split(':', 1)
                            pins[protein.strip()] = batch.strip()
            except Exception as e:
                print(f"Warning: Error loading protein pins: {e}")
        
        return pins
    
    def suggest_similar_proteins(self, protein_id: str, max_suggestions: int = 5) -> List[str]:
        """Suggest similar protein IDs that exist (unchanged from original)"""
        pdb_id = protein_id.split('_')[0] if '_' in protein_id else protein_id[:4]
        
        suggestions = []
        for batch_name in self._get_available_batches():
            proteins = self._get_proteins_in_batch(batch_name)
            same_pdb = [p for p in proteins if p.startswith(pdb_id)]
            suggestions.extend(same_pdb)
            
            if len(suggestions) >= max_suggestions:
                break
        
        unique_suggestions = list(dict.fromkeys(suggestions))[:max_suggestions]
        return [s for s in unique_suggestions if s != protein_id]
    
    def _get_available_batches(self) -> List[str]:
        """Get list of available batch directories"""
        if not self.base_dir.exists():
            return []
        
        batch_dirs = [d.name for d in self.base_dir.iterdir() 
                     if d.is_dir() and d.name.startswith("ecod_batch_")]
        
        return sorted(batch_dirs)
    
    def _protein_exists_in_batch(self, protein_id: str, batch_name: str) -> bool:
        """Check if protein exists in a specific batch"""
        batch_dir = self.base_dir / batch_name
        domain_file = batch_dir / "domains" / f"{protein_id}.develop291.domain_summary.xml"
        return domain_file.exists()
    
    def _get_proteins_in_batch(self, batch_name: str) -> List[str]:
        """Get list of proteins in a batch (cached)"""
        if batch_name in self._batch_cache:
            return self._batch_cache[batch_name]
        
        batch_dir = self.base_dir / batch_name
        domains_dir = batch_dir / "domains"
        
        if not domains_dir.exists():
            self._batch_cache[batch_name] = []
            return []
        
        proteins = []
        for domain_file in domains_dir.glob("*.develop291.domain_summary.xml"):
            protein_id = domain_file.stem.replace('.develop291.domain_summary', '')
            proteins.append(protein_id)
        
        self._batch_cache[batch_name] = sorted(proteins)
        return self._batch_cache[batch_name]

# Utility functions for configuration management

def pin_protein_to_batch(protein_id: str, batch_name: str):
    """Pin a protein to always use a specific batch"""
    pins_file = Path(__file__).parent / ".protein_pins"
    
    # Read existing pins
    existing_pins = {}
    if pins_file.exists():
        with open(pins_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and ':' in line and not line.startswith('#'):
                    protein, batch = line.split(':', 1)
                    existing_pins[protein.strip()] = batch.strip()
    
    # Add/update pin
    existing_pins[protein_id] = batch_name
    
    # Write back
    with open(pins_file, 'w') as f:
        f.write("# Protein to batch pins for reproducible results\n")
        f.write("# Format: protein_id: batch_name\n")
        f.write("\n")
        for protein, batch in sorted(existing_pins.items()):
            f.write(f"{protein}: {batch}\n")
    
    print(f"Pinned {protein_id} to batch {batch_name}")

def analyze_protein_batches(protein_id: str, base_dir: str = "/data/ecod/pdb_updates/batches"):
    """Analyze a protein across multiple batches"""
    finder = RobustBatchFinder(base_dir)
    analysis = finder.analyze_multi_batch_protein(protein_id)
    
    if not analysis["multi_batch"]:
        print(f"{protein_id} exists in only one batch: {analysis['batches']}")
        return
    
    print(f"{protein_id} exists in {len(analysis['batches'])} batches:")
    print("=" * 60)
    
    for batch_name in analysis["batches"]:
        info = analysis["analysis"][batch_name]
        print(f"\n{batch_name}:")
        print(f"  Domain summary: {'✓' if info['has_domain_summary'] else '✗'}")
        print(f"  Domains file: {'✓' if info['has_domains'] else '✗'}")
        print(f"  BLAST file: {'✓' if info['has_blast'] else '✗'}")
        
        if info.get("file_sizes"):
            print(f"  File sizes:")
            for file_type, size in info["file_sizes"].items():
                print(f"    {file_type}: {size:,} bytes")
    
    print(f"\nRecommendation:")
    print(f"  1. Test with different batches to compare results")
    print(f"  2. Pin to a specific batch for reproducibility:")
    print(f"     python -c \"from mini.robust_batch_detection import pin_protein_to_batch; pin_protein_to_batch('{protein_id}', 'CHOSEN_BATCH')\"")

if __name__ == "__main__":
    # Test the robust batch detection
    import sys
    
    if len(sys.argv) > 1:
        protein_id = sys.argv[1]
        print(f"Analyzing {protein_id} across batches...")
        analyze_protein_batches(protein_id)
    else:
        print("Usage: python robust_batch_detection.py PROTEIN_ID")
        print("Example: python robust_batch_detection.py 8ovp_A")
