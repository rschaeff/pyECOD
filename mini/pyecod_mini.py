#!/usr/bin/env python3
"""
pyecod_mini - Minimal Domain Partitioning Tool

A standalone domain partitioning tool with intelligent defaults and 
automatic discovery of data files within the pyECOD batch structure.

Usage:
    pyecod_mini 8ovp_A                    # Use all defaults
    pyecod_mini 8ovp_A --verbose          # Show detailed output
    pyecod_mini 8ovp_A --batch-id 036     # Use specific batch
    pyecod_mini --list-batches             # Show available batches
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
    """Configuration manager with intelligent defaults"""
    
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
            'output': self.output_dir / f"{protein_id}_mini.domains.xml"
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
            ("domain definitions", self.domain_definitions_file)
        ]:
            if not path.exists():
                issues.append(f"{name.title()} file not found: {path}")
        
        # Check if any batches exist
        if not self.list_available_batches():
            issues.append("No batch directories found")
        
        return len(issues) == 0, issues

def partition_protein(protein_id: str, config: PyEcodMiniConfig, 
                     batch_id: Optional[str] = None, verbose: bool = False) -> Optional[List]:
    """
    Partition domains for a single protein with full decomposition support.
    
    Args:
        protein_id: Protein to process (e.g., "8ovp_A")
        config: Configuration object
        batch_id: Optional batch ID (uses default if None)
        verbose: Enable detailed output
        
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
        
        return domains
        
    except Exception as e:
        print(f"ERROR processing {protein_id}: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        return None

def main():
    """Main entry point"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='pyECOD Mini - Minimal Domain Partitioning Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  pyecod_mini 8ovp_A                    # Process with all defaults
  pyecod_mini 8ovp_A --verbose          # Show detailed output  
  pyecod_mini 8ovp_A --batch-id 036     # Use specific batch
  pyecod_mini --list-batches             # Show available batches
  pyecod_mini --validate                 # Check configuration
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
    
    # Information commands
    parser.add_argument('--list-batches', action='store_true',
                        help='List available batch directories')
    parser.add_argument('--validate', action='store_true',
                        help='Validate configuration and show status')
    
    args = parser.parse_args()
    
    # Initialize configuration
    config = PyEcodMiniConfig()
    
    # Handle information commands
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
    
    # Main processing
    if not args.protein_id:
        parser.print_help()
        print("\nERROR: protein_id is required")
        print("Use --list-batches to see available data")
        sys.exit(1)
    
    # Validate configuration
    is_valid, issues = config.validate_setup()
    if not is_valid:
        print("Configuration issues found:")
        for issue in issues:
            print(f"  - {issue}")
        print("\nUse --validate for more details")
        sys.exit(1)
    
    # Process protein
    result = partition_protein(args.protein_id, config, args.batch_id, args.verbose)
    
    if result is None:
        sys.exit(1)
    else:
        print(f"\n✓ Successfully processed {args.protein_id}")

if __name__ == "__main__":
    main()
