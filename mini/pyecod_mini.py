#!/usr/bin/env python3
"""
pyecod_mini - Clean Minimal Domain Partitioning Tool

A standalone domain partitioning tool with smart batch detection,
visualization, and testing capabilities.

Usage:
    pyecod_mini 8ovp_A                    # Basic domain partitioning
    pyecod_mini 8ovp_A --visualize        # With PyMOL comparison
    pyecod_mini --test-suite               # Run formal test cases
    pyecod_mini --setup-references         # Generate reference files
"""

import sys
import os
from pathlib import Path
from typing import Optional, Dict, List, Tuple

# Add parent directory to path so we can import ecod
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.parser import parse_domain_summary, load_reference_lengths, load_protein_lengths
from mini.partitioner import partition_domains
from mini.writer import write_domain_partition
from mini.decomposer import load_domain_definitions
from mini.blast_parser import load_chain_blast_alignments

class BatchFinder:
    """Simple, robust batch finder for proteins"""

    def __init__(self, base_dir: str):
        self.base_dir = Path(base_dir)
        self._batch_cache = {}

        # Known stable batches for test cases
        self.stable_batches = {
            "8ovp_A": "ecod_batch_036_20250406_1424",  # Our validated test case
        }

    def find_batch_for_protein(self, protein_id: str, verbose: bool = False) -> Optional[str]:
        """Find which batch contains the protein"""

        # Check for known stable test cases first
        if protein_id in self.stable_batches:
            stable_batch = self.stable_batches[protein_id]
            if self._protein_exists_in_batch(protein_id, stable_batch):
                if verbose:
                    print(f"Using stable batch for {protein_id}: {stable_batch}")
                return stable_batch
            else:
                if verbose:
                    print(f"WARNING: Stable batch {stable_batch} doesn't contain {protein_id}")

        # Search all available batches
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

        # Multiple batches found - warn and use most recent
        if verbose:
            print(f"  ⚠️  Found in {len(found_batches)} batches: {found_batches}")
            print(f"  → Using most recent: {found_batches[-1]}")
            print(f"  💡 Use --batch-id to specify a particular batch")

        return found_batches[-1]  # Most recent

    def suggest_similar_proteins(self, protein_id: str, max_suggestions: int = 5) -> List[str]:
        """Suggest similar protein IDs that exist"""
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

    def analyze_protein_batches(self, protein_id: str) -> Dict[str, any]:
        """Analyze a protein across multiple batches"""
        available_batches = self._get_available_batches()
        found_batches = []

        for batch_name in available_batches:
            if self._protein_exists_in_batch(protein_id, batch_name):
                found_batches.append(batch_name)

        if len(found_batches) <= 1:
            return {"multi_batch": False, "batches": found_batches}

        return {"multi_batch": True, "batches": found_batches}

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

class PyEcodMiniConfig:
    """Configuration manager with integrated batch detection"""

    def __init__(self):
        self.base_dir = Path("/data/ecod/pdb_updates/batches")
        self.test_data_dir = Path(__file__).parent / "test_data"
        self.output_dir = Path("/tmp")

        # Default reference files
        self.domain_lengths_file = self.test_data_dir / "domain_lengths.csv"
        self.protein_lengths_file = self.test_data_dir / "protein_lengths.csv"
        self.domain_definitions_file = self.test_data_dir / "domain_definitions.csv"
        self.reference_blacklist_file = self.test_data_dir / "reference_blacklist.csv"

        # Batch finder
        self.batch_finder = BatchFinder(str(self.base_dir))

        # Visualization settings
        self.pdb_repo_path = "/usr2/pdb/data"
        self.visualization_output_dir = "/tmp/pymol_comparison"

    def get_batch_for_protein(self, protein_id: str, batch_id: Optional[str] = None,
                             verbose: bool = False) -> str:
        """Get the right batch for a protein"""
        if batch_id is not None:
            # Explicit batch specified - validate and resolve
            return self._resolve_batch_name(batch_id)

        # Smart detection
        found_batch = self.batch_finder.find_batch_for_protein(protein_id, verbose)

        if found_batch is None:
            # Provide helpful error message
            suggestions = self.batch_finder.suggest_similar_proteins(protein_id)
            error_msg = f"Protein {protein_id} not found in any batch"

            if suggestions:
                error_msg += f"\nSimilar proteins available: {suggestions[:3]}"

            available_batches = self.batch_finder._get_available_batches()
            if available_batches:
                error_msg += f"\nAvailable batches: {available_batches[-3:]}"
                error_msg += f"\nTo specify a batch: pyecod_mini {protein_id} --batch-id BATCH_NAME"

            raise FileNotFoundError(error_msg)

        return found_batch

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
            if (self.base_dir / batch_id).exists():
                return batch_id
            else:
                raise ValueError(f"Batch directory not found: {batch_id}")

    def get_batch_dir(self, batch_name: str) -> Path:
        """Get batch directory path from batch name"""
        return self.base_dir / batch_name

    def get_paths_for_protein(self, protein_id: str, batch_id: Optional[str] = None,
                             verbose: bool = False) -> Dict[str, Path]:
        """Get all file paths for a protein"""

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

    def list_available_batches(self) -> List[Tuple[str, int]]:
        """List all available batch directories with protein counts"""
        batches = self.batch_finder._get_available_batches()

        batch_info = []
        for batch_name in batches:
            proteins = self.batch_finder._get_proteins_in_batch(batch_name)
            batch_info.append((batch_name, len(proteins)))

        return batch_info

    def validate_setup(self) -> Tuple[bool, List[str]]:
        """Validate that the configuration is usable"""
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

        # Domain definitions are optional but recommended
        if not self.domain_definitions_file.exists():
            issues.append(f"Domain definitions file not found (chain BLAST decomposition disabled): {self.domain_definitions_file}")

        # Blacklist is optional
        if not self.reference_blacklist_file.exists():
            if verbose:
                print(f"No reference blacklist found: {self.reference_blacklist_file}")

            return len(issues) == 0, issues

def partition_protein(protein_id: str, config: PyEcodMiniConfig,
                     batch_id: Optional[str] = None,
                     verbose: bool = False,
                     visualize: bool = False) -> Optional[List]:
    """Partition domains for a single protein"""

    try:
        # Get all paths with smart batch detection
        paths = config.get_paths_for_protein(protein_id, batch_id, verbose)

        if verbose:
            print(f"Processing: {protein_id}")
            print(f"Batch: {paths['batch_name']}")
            print(f"Domain summary: {paths['domain_summary']}")
            print(f"BLAST XML: {paths['blast_xml']}")

        # Validate required files exist
        required_files = [
            ('domain_summary', 'Domain summary'),
            ('domain_lengths', 'Domain lengths'),
            ('protein_lengths', 'Protein lengths')
        ]

        for key, description in required_files:
            if not paths[key].exists():
                print(f"ERROR: {description} file not found: {paths[key]}")
                return None

        # Load reference data
        if verbose:
            print("\nLoading reference data...")

        domain_lengths = load_reference_lengths(str(paths['domain_lengths']))
        protein_lengths = load_protein_lengths(str(paths['protein_lengths']))

        # Load domain definitions (optional but recommended)
        domain_definitions = {}
        if paths['domain_definitions'].exists():
            blacklist_path = str(config.reference_blacklist_file) if config.reference_blacklist_file.exists() else None
            domain_definitions = load_domain_definitions(
                str(paths['domain_definitions']),
                verbose=verbose,
                blacklist_path=blacklist_path  # PASS BLACKLIST FILE PATH
            )
            if verbose:
                print(f"Loaded domain definitions for {len(domain_definitions)} protein chains")
        else:
            print("WARNING: Domain definitions not found - chain BLAST decomposition disabled")

        # Load BLAST alignments (optional but needed for accurate decomposition)
        blast_alignments = {}
        if paths['blast_xml'].exists():
            # Parse protein ID
            parts = protein_id.split('_')
            pdb_id, chain_id = parts[0], parts[1] if len(parts) > 1 else 'A'

            blast_alignments = load_chain_blast_alignments(
                str(paths['blast_dir']), pdb_id, chain_id, verbose=verbose
            )
            if verbose:
                print(f"Loaded {len(blast_alignments)} BLAST alignments")
        else:
            if verbose:
                print("WARNING: BLAST XML not found - alignment-based decomposition disabled")

        # Parse evidence
        if verbose:
            print(f"\nParsing evidence from {paths['domain_summary'].name}...")

        evidence = parse_domain_summary(
            str(paths['domain_summary']),
            reference_lengths=domain_lengths,
            protein_lengths=protein_lengths,
            blast_alignments=blast_alignments,
            require_reference_lengths=True,
            verbose=verbose
        )

        if not evidence:
            print("ERROR: No evidence with reference lengths found")
            return None

        # Show evidence summary
        evidence_by_type = {}
        for ev in evidence:
            evidence_by_type[ev.type] = evidence_by_type.get(ev.type, 0) + 1

        print(f"\nFound {len(evidence)} evidence items:")
        for etype, count in sorted(evidence_by_type.items()):
            print(f"  {etype}: {count}")

        # Estimate sequence length
        max_pos = max(ev.query_range.segments[-1].end for ev in evidence)
        sequence_length = int(max_pos * 1.1)

        if verbose:
            print(f"Estimated sequence length: {sequence_length}")

        # Show decomposition readiness
        chain_blast_count = evidence_by_type.get('chain_blast', 0)
        with_alignments = sum(1 for e in evidence if e.type == 'chain_blast' and e.alignment)

        print(f"\nDecomposition status:")
        print(f"  Chain BLAST evidence: {chain_blast_count}")
        print(f"  With alignment data: {with_alignments}")
        print(f"  Domain definitions: {'✓' if domain_definitions else '✗'}")

        # Partition domains
        print(f"\nPartitioning domains...")

        domains = partition_domains(
            evidence,
            sequence_length=sequence_length,
            domain_definitions=domain_definitions if domain_definitions else None,
            verbose=verbose
        )

        if not domains:
            print("No domains found")
            return None

        # Show results
        print(f"\n{'='*50}")
        print(f"RESULTS: {len(domains)} domains found")
        print(f"{'='*50}")

        total_coverage = sum(d.range.total_length for d in domains)
        coverage_pct = (total_coverage / sequence_length) * 100 if sequence_length > 0 else 0

        for i, domain in enumerate(domains, 1):
            print(f"\n{i}. Domain {domain.id}:")
            print(f"   Family: {domain.family}")
            print(f"   Range: {domain.range}")
            print(f"   Size: {domain.range.total_length} residues")
            print(f"   Source: {domain.source}")
            if domain.range.is_discontinuous:
                print(f"   Discontinuous: {len(domain.range.segments)} segments")

        print(f"\nTotal coverage: {total_coverage}/{sequence_length} residues ({coverage_pct:.1f}%)")

        # Write output
        parts = protein_id.split('_')
        pdb_id, chain_id = parts[0], parts[1] if len(parts) > 1 else 'A'

        write_domain_partition(domains, pdb_id, chain_id, str(paths['output']))
        print(f"\n✓ Output written to: {paths['output']}")

        # Generate visualization if requested
        if visualize:
            print(f"\nGenerating PyMOL comparison...")
            try:
                from mini.visualization import quick_comparison
                script_path = quick_comparison(protein_id, str(paths['batch_dir']))
                print(f"✓ PyMOL script: {script_path}")
                print(f"  Run: pymol {script_path}")
            except Exception as e:
                print(f"⚠️  Visualization failed: {e}")

        return domains

    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        return None
    except Exception as e:
        print(f"ERROR processing {protein_id}: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        return None

def analyze_protein_batches(protein_id: str, config: PyEcodMiniConfig) -> bool:
    """Analyze a protein across multiple batches"""

    analysis = config.batch_finder.analyze_protein_batches(protein_id)

    if not analysis["multi_batch"]:
        if analysis["batches"]:
            print(f"{protein_id} exists in only one batch: {analysis['batches'][0]}")
        else:
            print(f"{protein_id} not found in any batch")
        return True

    print(f"{protein_id} exists in {len(analysis['batches'])} batches:")
    print("=" * 60)

    for batch_name in analysis["batches"]:
        print(f"  {batch_name}")

    print(f"\n💡 Recommendations:")
    print(f"  • Test different batches: pyecod_mini {protein_id} --batch-id BATCH_NAME")
    print(f"  • Compare results to choose the best batch")
    print(f"  • For test cases, we use: ecod_batch_036_20250406_1424 (known stable)")

    return True

def setup_references(cache_file: str = None, output_dir: str = None) -> bool:
    """Setup reference files from ECOD range cache"""

    if cache_file is None:
        cache_file = "/data/ecod/database_versions/v291/ecod.develop291.range_cache.txt"

    if output_dir is None:
        output_dir = str(Path(__file__).parent / "test_data")

    print(f"Setting up reference files...")
    print(f"  Cache file: {cache_file}")
    print(f"  Output directory: {output_dir}")

    if not os.path.exists(cache_file):
        print(f"ERROR: Cache file not found: {cache_file}")
        return False

    try:
        from mini.range_cache_parser import (
            create_domain_lengths_from_cache,
            create_domain_definitions_from_cache,
            extract_protein_lengths_from_cache
        )

        os.makedirs(output_dir, exist_ok=True)

        print("Generating domain lengths...")
        create_domain_lengths_from_cache(cache_file, f"{output_dir}/domain_lengths.csv")

        print("Generating domain definitions...")
        create_domain_definitions_from_cache(cache_file, f"{output_dir}/domain_definitions.csv")

        print("Generating protein lengths...")
        extract_protein_lengths_from_cache(cache_file, f"{output_dir}/protein_lengths.csv")

        print("✅ Reference files generated successfully")
        return True

    except Exception as e:
        print(f"ERROR: Failed to generate reference files: {e}")
        return False

def run_test_suite(config: PyEcodMiniConfig, verbose: bool = False) -> bool:
    """Run the formal test suite"""

    print("Running formal test suite...")

    try:
        from mini.tests.test_cases import TestRunner

        # Use the first available batch for the test runner
        available_batches = config.batch_finder._get_available_batches()
        if not available_batches:
            print("ERROR: No batches available for testing")
            return False

        batch_dir = str(config.get_batch_dir(available_batches[0]))
        runner = TestRunner(batch_dir, str(config.test_data_dir))
        summary = runner.run_all_tests(verbose)

        if summary['failed'] > 0:
            print(f"\n❌ {summary['failed']} test(s) failed")
            return False
        else:
            print(f"\n✅ All {summary['passed']} test(s) passed")
            return True

    except Exception as e:
        print(f"ERROR: Failed to run test suite: {e}")
        return False

def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description='pyECOD Mini - Clean Domain Partitioning Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  pyecod_mini 8ovp_A                    # Basic domain partitioning
  pyecod_mini 8ovp_A --visualize        # With PyMOL comparison
  pyecod_mini 8ovp_A --batch-id 036     # Use specific batch
  pyecod_mini 8ovp_A --analyze-batches  # Check if protein exists in multiple batches
  pyecod_mini --test-suite               # Run formal test cases
  pyecod_mini --setup-references         # Generate reference files
  pyecod_mini --list-batches             # Show available batches
        """
    )

    # Main action
    parser.add_argument('protein_id', nargs='?',
                        help='Protein ID to process (e.g., 8ovp_A)')

    # Options
    parser.add_argument('--batch-id',
                        help='Batch ID (number like "036" or full name)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Show detailed output')
    parser.add_argument('--visualize', action='store_true',
                        help='Generate PyMOL comparison visualization')

    # Utility commands
    parser.add_argument('--list-batches', action='store_true',
                        help='List available batch directories')
    parser.add_argument('--validate', action='store_true',
                        help='Validate configuration and show status')
    parser.add_argument('--test-suite', action='store_true',
                        help='Run formal test suite')
    parser.add_argument('--setup-references', action='store_true',
                        help='Generate reference files from ECOD cache')
    parser.add_argument('--analyze-batches', action='store_true',
                        help='Analyze protein across multiple batches')
    parser.add_argument('--cache-file',
                        default='/data/ecod/database_versions/v291/ecod.develop291.range_cache.txt',
                        help='ECOD range cache file (for --setup-references)')

    args = parser.parse_args()

    # Initialize configuration
    config = PyEcodMiniConfig()

    # Handle utility commands
    if args.list_batches:
        batch_info = config.list_available_batches()
        if batch_info:
            print("Available batches:")
            for batch_name, protein_count in batch_info:
                print(f"  {batch_name} ({protein_count:,} proteins)")
        else:
            print("No batch directories found")
        return

    if args.validate:
        print("Validating pyECOD Mini configuration...")
        print(f"Base directory: {config.base_dir}")
        print(f"Test data directory: {config.test_data_dir}")

        # Show some available batches
        batch_info = config.list_available_batches()
        if batch_info:
            print(f"Available batches: {len(batch_info)} found")
            # Show most recent few
            for batch_name, protein_count in batch_info[-3:]:
                print(f"  {batch_name} ({protein_count:,} proteins)")

        is_valid, issues = config.validate_setup()
        if is_valid:
            print("✓ Configuration is valid")
        else:
            print("✗ Configuration issues found:")
            for issue in issues:
                print(f"  - {issue}")
        return

    if args.setup_references:
        success = setup_references(args.cache_file, str(config.test_data_dir))
        if not success:
            sys.exit(1)
        return

    if args.test_suite:
        success = run_test_suite(config, args.verbose)
        if not success:
            sys.exit(1)
        return

    if args.analyze_batches:
        if not args.protein_id:
            print("ERROR: --analyze-batches requires a protein_id")
            print("Usage: pyecod_mini PROTEIN_ID --analyze-batches")
            sys.exit(1)

        success = analyze_protein_batches(args.protein_id, config)
        if not success:
            sys.exit(1)
        return

    # Main processing
    if not args.protein_id:
        parser.print_help()
        print("\nERROR: protein_id is required for domain partitioning")
        print("Use --test-suite to run tests or --setup-references to prepare data")
        sys.exit(1)

    # Validate configuration
    is_valid, issues = config.validate_setup()
    if not is_valid:
        print("Configuration issues found:")
        for issue in issues:
            print(f"  - {issue}")
        print("\nUse --validate for more details or --setup-references to fix")
        sys.exit(1)

    # Process protein
    result = partition_protein(args.protein_id, config, args.batch_id,
                              args.verbose, args.visualize)

    if result is None:
        sys.exit(1)
    else:
        print(f"\n✅ Successfully processed {args.protein_id}")
        if args.visualize:
            print("✅ Visualization generated")

if __name__ == "__main__":
    main()
