#!/usr/bin/env python3
"""
Smart batch detection that finds the right batch for a protein

Instead of blindly picking the "latest" batch, this searches across
available batches to find where the protein actually exists.
"""

import os
from pathlib import Path
from typing import Optional, List, Tuple, Dict

class SmartBatchFinder:
    """Find the right batch for a protein without database access"""

    def __init__(self, base_dir: str = "/data/ecod/pdb_updates/batches"):
        self.base_dir = Path(base_dir)
        self._batch_cache = {}  # Cache batch contents

    def find_batch_for_protein(self, protein_id: str, verbose: bool = False) -> Optional[str]:
        """
        Find which batch contains the given protein.

        Args:
            protein_id: Protein ID (e.g., "8ovp_A")
            verbose: Print search details

        Returns:
            Batch name that contains the protein, or None if not found
        """
        available_batches = self._get_available_batches()

        if verbose:
            print(f"Searching for {protein_id} across {len(available_batches)} batches...")

        found_batches = []

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

        # Multiple batches have this protein - pick the most recent
        if verbose:
            print(f"  Found in multiple batches: {found_batches}")
            print(f"  Using most recent: {found_batches[-1]}")

        return found_batches[-1]  # Most recent

    def suggest_similar_proteins(self, protein_id: str, max_suggestions: int = 5) -> List[str]:
        """
        Suggest similar protein IDs that exist in available batches.

        Args:
            protein_id: Target protein ID
            max_suggestions: Maximum number of suggestions

        Returns:
            List of similar protein IDs that exist
        """
        # Extract PDB ID from protein_id
        pdb_id = protein_id.split('_')[0] if '_' in protein_id else protein_id[:4]

        suggestions = []
        available_batches = self._get_available_batches()

        for batch_name in available_batches:
            proteins = self._get_proteins_in_batch(batch_name)

            # Look for proteins with same PDB ID
            same_pdb = [p for p in proteins if p.startswith(pdb_id)]
            suggestions.extend(same_pdb)

            if len(suggestions) >= max_suggestions:
                break

        # Remove duplicates and limit
        unique_suggestions = list(dict.fromkeys(suggestions))[:max_suggestions]

        return [s for s in unique_suggestions if s != protein_id]

    def get_batch_info(self, batch_name: str) -> Dict[str, any]:
        """Get information about a batch"""
        batch_dir = self.base_dir / batch_name

        if not batch_dir.exists():
            return {'exists': False}

        domains_dir = batch_dir / "domains"
        blast_dir = batch_dir / "blast" / "chain"

        info = {
            'exists': True,
            'path': str(batch_dir),
            'has_domains': domains_dir.exists(),
            'has_blast': blast_dir.exists(),
        }

        if domains_dir.exists():
            domain_files = list(domains_dir.glob("*.domain_summary.xml"))
            info['protein_count'] = len(domain_files)
            info['sample_proteins'] = [
                f.stem.replace('.develop291.domain_summary', '')
                for f in domain_files[:5]
            ]

        return info

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

class EnhancedPyEcodMiniConfig:
    """Enhanced configuration with smart batch detection"""

    def __init__(self):
        self.base_dir = Path("/data/ecod/pdb_updates/batches")
        self.test_data_dir = Path(__file__).parent / "test_data"
        self.output_dir = Path("/tmp")

        # Reference files
        self.domain_lengths_file = self.test_data_dir / "domain_lengths.csv"
        self.protein_lengths_file = self.test_data_dir / "protein_lengths.csv"
        self.domain_definitions_file = self.test_data_dir / "domain_definitions.csv"

        # Smart batch finder
        self.batch_finder = SmartBatchFinder(str(self.base_dir))

        # Default batch strategy: None (force smart detection)
        self.default_batch = None

        # Visualization settings
        self.pdb_repo_path = "/usr2/pdb/data"
        self.visualization_output_dir = "/tmp/pymol_comparison"

    def get_batch_for_protein(self, protein_id: str, batch_id: Optional[str] = None,
                             verbose: bool = False) -> Optional[str]:
        """
        Get the right batch for a protein using smart detection.

        Args:
            protein_id: Protein ID
            batch_id: Optional explicit batch ID
            verbose: Print search details

        Returns:
            Batch name or None if not found
        """
        if batch_id is not None:
            # Explicit batch specified - validate it exists
            if self._validate_batch_id(batch_id):
                return self._resolve_batch_name(batch_id)
            else:
                raise ValueError(f"Specified batch not found: {batch_id}")

        # Smart detection: find batch that contains this protein
        found_batch = self.batch_finder.find_batch_for_protein(protein_id, verbose)

        if found_batch is None:
            # Provide helpful error message with suggestions
            suggestions = self.batch_finder.suggest_similar_proteins(protein_id)
            error_msg = f"Protein {protein_id} not found in any batch"

            if suggestions:
                error_msg += f"\\nSimilar proteins available: {suggestions[:3]}"

            available_batches = self.batch_finder._get_available_batches()
            if available_batches:
                error_msg += f"\\nAvailable batches: {available_batches[-3:]}"

            raise FileNotFoundError(error_msg)

        return found_batch

    def _validate_batch_id(self, batch_id: str) -> bool:
        """Check if batch_id corresponds to an existing batch"""
        try:
            batch_name = self._resolve_batch_name(batch_id)
            batch_dir = self.base_dir / batch_name
            return batch_dir.exists()
        except:
            return False

    def _resolve_batch_name(self, batch_id: str) -> str:
        """Convert batch_id to full batch name"""
        if batch_id.isdigit():
            # Number like "036" -> find "ecod_batch_036_*"
            pattern = f"ecod_batch_{batch_id.zfill(3)}_*"
            matches = list(self.base_dir.glob(pattern))
            if matches:
                return matches[0].name
            else:
                raise ValueError(f"No batch found matching pattern: {pattern}")
        else:
            # Assume it's already a full batch name
            return batch_id

    def get_batch_dir(self, batch_name: str) -> Path:
        """Get batch directory path from batch name"""
        return self.base_dir / batch_name

    def get_paths_for_protein(self, protein_id: str, batch_id: Optional[str] = None,
                             verbose: bool = False) -> Dict[str, Path]:
        """Get all file paths for a protein with smart batch detection"""

        # Find the right batch for this protein
        batch_name = self.get_batch_for_protein(protein_id, batch_id, verbose)
        batch_dir = self.get_batch_dir(batch_name)

        if verbose:
            print(f"Using batch: {batch_name}")

        return {
            'batch_dir': batch_dir,
            'batch_name': batch_name,
            'domain_summary': batch_dir / "domains" / f"{protein_id}.develop291.domain_summary.xml",
            'blast_xml': batch_dir / "blast" / "chain" / f"{protein_id}.develop291.xml",
            'blast_dir': batch_dir / "blast" / "chain",
            'domain_lengths': self.domain_lengths_file,
            'protein_lengths': self.protein_lengths_file,
            'domain_definitions': self.domain_definitions_file,
            'output': self.output_dir / f"{protein_id}_mini.domains.xml",
            'old_domains': batch_dir / "domains" / f"{protein_id}.develop291.domains.xml"
        }

    def list_available_batches(self) -> List[str]:
        """List all available batch directories with protein counts"""
        batches = self.batch_finder._get_available_batches()

        batch_info = []
        for batch_name in batches:
            info = self.batch_finder.get_batch_info(batch_name)
            protein_count = info.get('protein_count', 0)
            batch_info.append((batch_name, protein_count))

        return batch_info

    def validate_setup(self) -> Tuple[bool, List[str]]:
        """Validate configuration"""
        issues = []

        # Check base directory
        if not self.base_dir.exists():
            issues.append(f"Base directory not found: {self.base_dir}")
            return False, issues

        # Check if any batches exist
        available_batches = self.batch_finder._get_available_batches()
        if not available_batches:
            issues.append("No batch directories found")

        # Check test data files
        for name, path in [
            ("domain lengths", self.domain_lengths_file),
            ("protein lengths", self.protein_lengths_file),
        ]:
            if not path.exists():
                issues.append(f"{name.title()} file not found: {path}")

        # Domain definitions are optional
        if not self.domain_definitions_file.exists():
            issues.append(f"Domain definitions file not found (chain BLAST decomposition disabled): {self.domain_definitions_file}")

        return len(issues) == 0, issues

# Example usage and testing
def test_smart_batch_detection():
    """Test the smart batch detection"""
    finder = SmartBatchFinder()

    # Test cases
    test_proteins = ["8ovp_A", "1ubq_A", "nonexistent_protein"]

    print("Testing smart batch detection:")
    print("=" * 40)

    for protein in test_proteins:
        print(f"\\nSearching for {protein}:")
        batch = finder.find_batch_for_protein(protein, verbose=True)

        if batch:
            print(f"  → Found in: {batch}")

            # Show batch info
            info = finder.get_batch_info(batch)
            print(f"  → Batch has {info.get('protein_count', 0)} proteins")
        else:
            print(f"  → Not found")
            suggestions = finder.suggest_similar_proteins(protein)
            if suggestions:
                print(f"  → Similar proteins: {suggestions}")

if __name__ == "__main__":
    test_smart_batch_detection()
