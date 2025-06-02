#!/usr/bin/env python3
"""
pyecod_mini - Enhanced Minimal Domain Partitioning Tool

A standalone domain partitioning tool with intelligent defaults,
automatic discovery of data files, visualization, and testing capabilities.

Usage:
    pyecod_mini 8ovp_A                    # Basic domain partitioning
    pyecod_mini 8ovp_A --visualize        # With PyMOL comparison
    pyecod_mini --test-suite               # Run formal test cases
    pyecod_mini --setup-references         # Generate reference files
"""

import sys
import os
import json
from pathlib import Path
from typing import Optional, Dict, List, Tuple

# Add parent directory to path so we can import ecod
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.parser import parse_domain_summary, load_reference_lengths, load_protein_lengths
from mini.partitioner import partition_domains
from mini.writer import write_domain_partition
from mini.decomposer import load_domain_definitions
from mini.blast_parser import load_chain_blast_alignments

class PyEcodMiniConfig:
    """Enhanced configuration manager with visualization and testing support"""

    def __init__(self):
        self.base_dir = Path("/data/ecod/pdb_updates/batches")
        self.test_data_dir = Path(__file__).parent / "test_data"
        self.output_dir = Path("/tmp")

        # Default reference files
        self.domain_lengths_file = self.test_data_dir / "domain_lengths.csv"
        self.protein_lengths_file = self.test_data_dir / "protein_lengths.csv"
        self.domain_definitions_file = self.test_data_dir / "domain_definitions.csv"

        # Default batch (latest)
        self.default_batch = self._find_latest_batch()

        # Visualization settings
        self.pdb_repo_path = "/usr2/pdb/data"
        self.visualization_output_dir = "/tmp/pymol_comparison"

    def _find_latest_batch(self) -> Optional[str]:
        """Find the most recent batch directory"""
        if not self.base_dir.exists():
            return None

        batch_dirs = [d for d in self.base_dir.iterdir()
                     if d.is_dir() and d.name.startswith("ecod_batch_")]

        if not batch_dirs:
            return None

        # Sort by name (which includes timestamp)
        latest = sorted(batch_dirs)[-1]
        return latest.name

    def get_batch_dir(self, batch_id: Optional[str] = None) -> Path:
        """Get batch directory path"""
        if batch_id is None:
            batch_id = self.default_batch

        if batch_id is None:
            raise ValueError("No batch ID specified and no default batch found")

        # Handle different batch ID formats
        if batch_id.isdigit():
            # Just number: "036" -> find "ecod_batch_036_*"
            pattern = f"ecod_batch_{batch_id.zfill(3)}_*"
            matches = list(self.base_dir.glob(pattern))
            if matches:
                return matches[0]  # Take first match
            else:
                raise ValueError(f"No batch found matching pattern: {pattern}")
        else:
            # Full name or partial name
            batch_dir = self.base_dir / batch_id
            if batch_dir.exists():
                return batch_dir
            else:
                raise ValueError(f"Batch directory not found: {batch_dir}")

    def get_paths_for_protein(self, protein_id: str, batch_id: Optional[str] = None) -> Dict[str, Path]:
        """Get all file paths for a protein"""
        batch_dir = self.get_batch_dir(batch_id)

        return {
            'batch_dir': batch_dir,
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
        """List all available batch directories"""
        if not self.base_dir.exists():
            return []

        batch_dirs = [d.name for d in self.base_dir.iterdir()
                     if d.is_dir() and d.name.startswith("ecod_batch_")]

        return sorted(batch_dirs)

    def validate_setup(self) -> Tuple[bool, List[str]]:
        """Validate that the configuration is usable"""
        issues = []

        # Check base directory
        if not self.base_dir.exists():
            issues.append(f"Base directory not found: {self.base_dir}")

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

        # Check if any batches exist
        if not self.list_available_batches():
            issues.append("No batch directories found")

        return len(issues) == 0, issues

def partition_protein(protein_id: str, config: PyEcodMiniConfig,
                     batch_id: Optional[str] = None,
                     verbose: bool = False,
                     visualize: bool = False) -> Optional[List]:
    """
    Partition domains for a single protein with full decomposition support.

    Args:
        protein_id: Protein to process (e.g., "8ovp_A")
        config: Configuration object
        batch_id: Optional batch ID (uses default if None)
        verbose: Enable detailed output
        visualize: Generate PyMOL comparison

    Returns:
        List of domains or None if failed
    """

    try:
        # Get all paths
        paths = config.get_paths_for_protein(protein_id, batch_id)

        if verbose:
            print(f"Processing: {protein_id}")
            print(f"Batch: {paths['batch_dir'].name}")
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
            domain_definitions = load_domain_definitions(str(paths['domain_definitions']), verbose=verbose)
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

    except Exception as e:
        print(f"ERROR processing {protein_id}: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        return None

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

        runner = TestRunner(str(config.get_batch_dir()), str(config.test_data_dir))
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
    """Main entry point with enhanced functionality"""
    import argparse

    parser = argparse.ArgumentParser(
        description='pyECOD Mini - Enhanced Domain Partitioning Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  pyecod_mini 8ovp_A                    # Basic domain partitioning
  pyecod_mini 8ovp_A --visualize        # With PyMOL comparison
  pyecod_mini 8ovp_A --batch-id 036     # Use specific batch
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
    parser.add_argument('--cache-file',
                        default='/data/ecod/database_versions/v291/ecod.develop291.range_cache.txt',
                        help='ECOD range cache file (for --setup-references)')

    args = parser.parse_args()

    # Initialize configuration
    config = PyEcodMiniConfig()

    # Handle utility commands
    if args.list_batches:
        batches = config.list_available_batches()
        if batches:
            print("Available batches:")
            for batch in batches:
                marker = " (default)" if batch == config.default_batch else ""
                print(f"  {batch}{marker}")
        else:
            print("No batch directories found")
        return

    if args.validate:
        print("Validating pyECOD Mini configuration...")
        print(f"Base directory: {config.base_dir}")
        print(f"Test data directory: {config.test_data_dir}")
        print(f"Default batch: {config.default_batch}")

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
