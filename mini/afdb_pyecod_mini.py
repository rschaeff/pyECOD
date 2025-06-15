#!/usr/bin/env python3
"""
afdb_pyecod_mini - AFDB-specific Domain Partitioning Tool

Modified version of pyecod_mini for processing AlphaFold Database (AFDB) proteins.
Key differences:
- Handles UniProt IDs (P21580_F1) instead of PDB IDs (8ovp_A)
- Uses /data/ecod/pyecod_mini_bench/batches/ directory structure
- Outputs TSV format for bulk processing
- No database integration requirements

Usage:
    afdb_pyecod_mini P21580_F1                    # Process single AFDB protein
    afdb_pyecod_mini --batch afdb_human_pyecod_batch_001  # Process entire batch
    afdb_pyecod_mini --validate                   # Validate AFDB data files
"""

import sys
import os
import csv
from pathlib import Path
from typing import Optional, Dict, List, Tuple
from dataclasses import dataclass

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.core.parser import parse_domain_summary, load_reference_lengths, load_protein_lengths
from mini.core.partitioner import partition_domains
from mini.core.decomposer import load_domain_definitions
from mini.core.blast_parser import load_chain_blast_alignments

@dataclass
class AFDBPartitionResult:
    """Results from AFDB domain partitioning"""
    protein_id: str
    success: bool
    domain_count: int = 0
    domains: List[Dict] = None
    error: str = ""
    coverage_percent: float = 0.0
    evidence_counts: Dict[str, int] = None

class AFDBBatchFinder:
    """Batch finder specifically for AFDB data structure"""
    
    def __init__(self, base_dir: str = "/data/ecod/pyecod_mini_bench/batches"):
        self.base_dir = Path(base_dir)
    
    def find_batch_for_protein(self, protein_id: str, verbose: bool = False) -> Optional[str]:
        """Find which AFDB batch contains the protein"""
        
        available_batches = self._get_available_batches()
        found_batches = []
        
        if verbose:
            print(f"Searching for {protein_id} across {len(available_batches)} AFDB batches...")
        
        for batch_name in available_batches:
            if self._protein_exists_in_batch(protein_id, batch_name):
                found_batches.append(batch_name)
                if verbose:
                    print(f"  ✓ Found in {batch_name}")
        
        if not found_batches:
            if verbose:
                print(f"  ✗ {protein_id} not found in any AFDB batch")
            return None
        
        if len(found_batches) == 1:
            return found_batches[0]
        
        # Multiple batches found - use most recent
        if verbose:
            print(f"  ⚠️  Found in {len(found_batches)} batches, using most recent: {found_batches[-1]}")
        
        return found_batches[-1]
    
    def _get_available_batches(self) -> List[str]:
        """Get list of available AFDB batch directories"""
        if not self.base_dir.exists():
            return []
        
        batch_dirs = [d.name for d in self.base_dir.iterdir()
                     if d.is_dir() and d.name.startswith("afdb_human_pyecod_batch_")]
        
        return sorted(batch_dirs)
    
    def _protein_exists_in_batch(self, protein_id: str, batch_name: str) -> bool:
        """Check if protein exists in a specific AFDB batch"""
        batch_dir = self.base_dir / batch_name
        domain_file = batch_dir / "summaries" / f"{protein_id}.develop291.domain_summary.xml"
        return domain_file.exists()
    
    def get_batch_proteins(self, batch_name: str) -> List[str]:
        """Get all proteins in an AFDB batch"""
        batch_dir = self.base_dir / batch_name
        summaries_dir = batch_dir / "summaries"
        
        if not summaries_dir.exists():
            return []
        
        proteins = []
        for summary_file in summaries_dir.glob("*.develop291.domain_summary.xml"):
            protein_id = summary_file.stem.replace('.develop291.domain_summary', '')
            proteins.append(protein_id)
        
        return sorted(proteins)

class AFDBConfig:
    """Configuration for AFDB processing"""
    
    def __init__(self, base_dir: str = "/data/ecod/pyecod_mini_bench/batches"):
        self.base_dir = Path(base_dir)
        self.test_data_dir = Path(__file__).parent / "test_data"
        self.output_dir = Path("/tmp/afdb_domains")
        
        # Reference files (same as regular mini)
        self.domain_lengths_file = self.test_data_dir / "domain_lengths.csv"
        self.protein_lengths_file = self.test_data_dir / "protein_lengths.csv"
        self.domain_definitions_file = self.test_data_dir / "domain_definitions.csv"
        self.reference_blacklist_file = self.test_data_dir / "reference_blacklist.csv"
        
        # AFDB-specific batch finder
        self.batch_finder = AFDBBatchFinder(str(self.base_dir))
        
        # Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def get_paths_for_protein(self, protein_id: str, batch_name: str) -> Dict[str, Path]:
        """Get file paths for an AFDB protein"""
        batch_dir = self.base_dir / batch_name
        
        return {
            'batch_dir': batch_dir,
            'batch_name': batch_name,
            'domain_summary': batch_dir / "summaries" / f"{protein_id}.develop291.domain_summary.xml",
            'blast_dir': batch_dir / "blast",  # AFDB may not have chain subdirectory
            'domain_lengths': self.domain_lengths_file,
            'protein_lengths': self.protein_lengths_file,
            'domain_definitions': self.domain_definitions_file,
            'output_xml': self.output_dir / f"{protein_id}_mini.domains.xml",
            'output_tsv': self.output_dir / f"{protein_id}_mini.domains.tsv"
        }
    
    def validate_setup(self) -> Tuple[bool, List[str]]:
        """Validate AFDB configuration"""
        issues = []
        
        if not self.base_dir.exists():
            issues.append(f"AFDB base directory not found: {self.base_dir}")
            return False, issues
        
        available_batches = self.batch_finder._get_available_batches()
        if not available_batches:
            issues.append("No AFDB batch directories found")
        
        # Check reference files
        for name, path in [
            ("domain lengths", self.domain_lengths_file),
            ("protein lengths", self.protein_lengths_file),
        ]:
            if not path.exists():
                issues.append(f"{name.title()} file not found: {path}")
        
        return len(issues) == 0, issues

def parse_afdb_protein_id(protein_id: str) -> Tuple[str, str]:
    """Parse AFDB protein ID into UniProt accession and fragment"""
    # Format: P21580_F1 -> ('P21580', 'F1')
    if '_' in protein_id:
        uniprot_acc, fragment = protein_id.split('_', 1)
        return uniprot_acc, fragment
    else:
        # Fallback: treat as single part
        return protein_id, 'F1'

def write_tsv_output(domains: List, protein_id: str, output_path: str, 
                    sequence_length: int = None, evidence_counts: Dict[str, int] = None):
    """Write domain results to TSV format for bulk analysis"""
    
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        
        # Header
        writer.writerow([
            'protein_id', 'domain_num', 'domain_range', 'domain_size', 
            'family', 't_group', 'h_group', 'x_group', 'source', 
            'is_discontinuous', 'evidence_count'
        ])
        
        if not domains:
            # Write empty row for proteins with no domains
            writer.writerow([
                protein_id, 0, '', 0, 'none', '', '', '', 'none', False, 0
            ])
        else:
            for i, domain in enumerate(domains, 1):
                writer.writerow([
                    protein_id,
                    i,
                    str(domain.range),
                    domain.range.total_length,
                    domain.family,
                    domain.t_group or '',
                    domain.h_group or '',
                    domain.x_group or '',
                    domain.source,
                    domain.range.is_discontinuous,
                    domain.evidence_count
                ])

def partition_afdb_protein(protein_id: str, config: AFDBConfig, 
                          batch_name: str = None, verbose: bool = False) -> AFDBPartitionResult:
    """Partition domains for a single AFDB protein"""
    
    try:
        # Find batch if not specified
        if batch_name is None:
            batch_name = config.batch_finder.find_batch_for_protein(protein_id, verbose)
            if batch_name is None:
                return AFDBPartitionResult(
                    protein_id=protein_id,
                    success=False,
                    error=f"Protein {protein_id} not found in any AFDB batch"
                )
        
        # Get file paths
        paths = config.get_paths_for_protein(protein_id, batch_name)
        
        if verbose:
            print(f"Processing AFDB protein: {protein_id}")
            print(f"Batch: {batch_name}")
            print(f"Domain summary: {paths['domain_summary']}")
        
        # Validate required files
        if not paths['domain_summary'].exists():
            return AFDBPartitionResult(
                protein_id=protein_id,
                success=False,
                error=f"Domain summary not found: {paths['domain_summary']}"
            )
        
        # Load reference data
        domain_lengths = load_reference_lengths(str(paths['domain_lengths']))
        protein_lengths = load_protein_lengths(str(paths['protein_lengths']))
        
        # Load domain definitions (optional)
        domain_definitions = {}
        if paths['domain_definitions'].exists():
            blacklist_path = str(config.reference_blacklist_file) if config.reference_blacklist_file.exists() else None
            domain_definitions = load_domain_definitions(
                str(paths['domain_definitions']),
                verbose=verbose,
                blacklist_path=blacklist_path
            )
        
        # Parse evidence
        evidence = parse_domain_summary(
            str(paths['domain_summary']),
            reference_lengths=domain_lengths,
            protein_lengths=protein_lengths,
            require_reference_lengths=True,
            verbose=verbose
        )
        
        # Count evidence by type
        evidence_counts = {}
        for ev in evidence:
            evidence_counts[ev.type] = evidence_counts.get(ev.type, 0) + 1
        
        if not evidence:
            # No evidence found - this is valid (0 domains)
            write_tsv_output([], protein_id, str(paths['output_tsv']))
            return AFDBPartitionResult(
                protein_id=protein_id,
                success=True,
                domain_count=0,
                domains=[],
                evidence_counts=evidence_counts
            )
        
        # Estimate sequence length from evidence
        max_pos = max(ev.query_range.segments[-1].end for ev in evidence)
        sequence_length = int(max_pos * 1.1)
        
        # Partition domains
        domains = partition_domains(
            evidence,
            sequence_length=sequence_length,
            domain_definitions=domain_definitions if domain_definitions else None,
            verbose=verbose
        )
        
        # Calculate coverage
        total_coverage = sum(d.range.total_length for d in domains) if domains else 0
        coverage_percent = (total_coverage / sequence_length * 100) if sequence_length > 0 else 0
        
        # Convert domains to simple dict format for TSV output
        domain_dicts = []
        for domain in domains:
            domain_dicts.append({
                'range': str(domain.range),
                'size': domain.range.total_length,
                'family': domain.family,
                't_group': domain.t_group,
                'h_group': domain.h_group,
                'x_group': domain.x_group,
                'source': domain.source,
                'is_discontinuous': domain.range.is_discontinuous,
                'evidence_count': domain.evidence_count
            })
        
        # Write outputs
        write_tsv_output(domains, protein_id, str(paths['output_tsv']), 
                        sequence_length, evidence_counts)
        
        # Also write XML for compatibility
        if domains:
            from mini.core.writer import write_domain_partition
            uniprot_acc, fragment = parse_afdb_protein_id(protein_id)
            write_domain_partition(domains, uniprot_acc, 'A', str(paths['output_xml']))
        
        return AFDBPartitionResult(
            protein_id=protein_id,
            success=True,
            domain_count=len(domains),
            domains=domain_dicts,
            coverage_percent=coverage_percent,
            evidence_counts=evidence_counts
        )
        
    except Exception as e:
        return AFDBPartitionResult(
            protein_id=protein_id,
            success=False,
            error=str(e)
        )

def process_afdb_batch(batch_name: str, config: AFDBConfig, 
                      max_proteins: int = None, verbose: bool = False) -> Dict[str, any]:
    """Process an entire AFDB batch"""
    
    proteins = config.batch_finder.get_batch_proteins(batch_name)
    
    if max_proteins:
        proteins = proteins[:max_proteins]
    
    print(f"Processing {len(proteins)} proteins from batch {batch_name}")
    
    results = []
    success_count = 0
    total_domains = 0
    
    # Batch output file
    batch_output = config.output_dir / f"{batch_name}_domains.tsv"
    
    with open(batch_output, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([
            'protein_id', 'domain_num', 'domain_range', 'domain_size', 
            'family', 't_group', 'h_group', 'x_group', 'source', 
            'is_discontinuous', 'evidence_count'
        ])
        
        for i, protein_id in enumerate(proteins, 1):
            if verbose:
                print(f"  {i}/{len(proteins)}: {protein_id}")
            
            result = partition_afdb_protein(protein_id, config, batch_name, verbose=False)
            results.append(result)
            
            if result.success:
                success_count += 1
                total_domains += result.domain_count
                
                # Write to batch file
                if result.domains:
                    for j, domain in enumerate(result.domains, 1):
                        writer.writerow([
                            protein_id, j, domain['range'], domain['size'],
                            domain['family'], domain['t_group'], domain['h_group'],
                            domain['x_group'], domain['source'], domain['is_discontinuous'],
                            domain['evidence_count']
                        ])
                else:
                    # No domains found
                    writer.writerow([
                        protein_id, 0, '', 0, 'none', '', '', '', 'none', False, 0
                    ])
            else:
                if verbose:
                    print(f"    ❌ {result.error}")
    
    print(f"\nBatch {batch_name} results:")
    print(f"  Processed: {len(proteins)} proteins")
    print(f"  Success: {success_count}")
    print(f"  Failed: {len(proteins) - success_count}")
    print(f"  Total domains: {total_domains}")
    print(f"  Avg domains/protein: {total_domains/success_count:.2f}" if success_count > 0 else "")
    print(f"  Output: {batch_output}")
    
    return {
        'batch_name': batch_name,
        'total_proteins': len(proteins),
        'success_count': success_count,
        'failed_count': len(proteins) - success_count,
        'total_domains': total_domains,
        'avg_domains': total_domains/success_count if success_count > 0 else 0,
        'output_file': str(batch_output),
        'results': results
    }

def validate_afdb_setup(config: AFDBConfig) -> bool:
    """Validate AFDB setup and data files"""
    
    print("Validating AFDB setup...")
    
    # Basic configuration validation
    is_valid, issues = config.validate_setup()
    
    if not is_valid:
        print("❌ Configuration issues:")
        for issue in issues:
            print(f"  - {issue}")
        return False
    
    # Check available batches
    available_batches = config.batch_finder._get_available_batches()
    print(f"✅ Found {len(available_batches)} AFDB batches:")
    for batch_name in available_batches:
        protein_count = len(config.batch_finder.get_batch_proteins(batch_name))
        print(f"  - {batch_name}: {protein_count:,} proteins")
    
    # Validate a sample of domain summary files
    if available_batches:
        print(f"\nValidating sample domain summaries...")
        
        # Import the validator
        sys.path.insert(0, str(Path(__file__).parent))
        try:
            from afdb_validator import AFDBDomainSummaryValidator
            
            validator = AFDBDomainSummaryValidator()
            sample_batch = available_batches[0]
            batch_dir = config.base_dir / sample_batch
            
            result = validator.validate_batch(str(batch_dir), sample_size=10)
            
            if result.get('batch_valid'):
                print(f"✅ Sample validation passed ({result['valid_files']}/{result['total_files']} files)")
            else:
                print(f"❌ Sample validation failed")
                if result.get('unique_errors'):
                    print("Errors found:")
                    for error in result['unique_errors'][:3]:
                        print(f"  - {error}")
                return False
        except ImportError:
            print("⚠️  Could not import validator - skipping domain summary validation")
    
    print("✅ AFDB setup validation complete")
    return True

def main():
    """Main entry point for AFDB processing"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='AFDB PyECOD Mini - AlphaFold Domain Partitioning',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  afdb_pyecod_mini P21580_F1                           # Process single protein
  afdb_pyecod_mini --batch afdb_human_pyecod_batch_001  # Process entire batch
  afdb_pyecod_mini --batch afdb_human_pyecod_batch_001 --max-proteins 100  # Process first 100
  afdb_pyecod_mini --validate                          # Validate setup
  afdb_pyecod_mini --list-batches                      # Show available batches
        """
    )
    
    # Main actions
    parser.add_argument('protein_id', nargs='?', help='AFDB protein ID (e.g., P21580_F1)')
    parser.add_argument('--batch', help='Process entire AFDB batch')
    parser.add_argument('--max-proteins', type=int, help='Limit number of proteins to process')
    
    # Options
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    parser.add_argument('--base-dir', default='/data/ecod/pyecod_mini_bench/batches',
                       help='Base directory for AFDB batches')
    parser.add_argument('--output-dir', default='/tmp/afdb_domains',
                       help='Output directory for results')
    
    # Utility commands
    parser.add_argument('--validate', action='store_true', help='Validate AFDB setup')
    parser.add_argument('--list-batches', action='store_true', help='List available batches')
    
    args = parser.parse_args()
    
    # Initialize configuration
    config = AFDBConfig(args.base_dir)
    if args.output_dir:
        config.output_dir = Path(args.output_dir)
        config.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Handle utility commands
    if args.validate:
        success = validate_afdb_setup(config)
        sys.exit(0 if success else 1)
    
    if args.list_batches:
        available_batches = config.batch_finder._get_available_batches()
        if available_batches:
            print("Available AFDB batches:")
            for batch_name in available_batches:
                protein_count = len(config.batch_finder.get_batch_proteins(batch_name))
                print(f"  {batch_name}: {protein_count:,} proteins")
        else:
            print("No AFDB batches found")
        return
    
    # Validate setup before processing
    if not validate_afdb_setup(config):
        print("❌ Setup validation failed")
        sys.exit(1)
    
    # Main processing
    if args.batch:
        # Process entire batch
        result = process_afdb_batch(args.batch, config, args.max_proteins, args.verbose)
        
        if result['success_count'] > 0:
            print(f"✅ Successfully processed {result['success_count']} proteins")
        if result['failed_count'] > 0:
            print(f"❌ Failed to process {result['failed_count']} proteins")
        
        sys.exit(0 if result['failed_count'] == 0 else 1)
    
    elif args.protein_id:
        # Process single protein
        result = partition_afdb_protein(args.protein_id, config, verbose=args.verbose)
        
        if result.success:
            print(f"✅ Successfully processed {args.protein_id}")
            print(f"   Domains found: {result.domain_count}")
            print(f"   Coverage: {result.coverage_percent:.1f}%")
            print(f"   Output: {config.output_dir / f'{args.protein_id}_mini.domains.tsv'}")
        else:
            print(f"❌ Failed to process {args.protein_id}: {result.error}")
            sys.exit(1)
    
    else:
        parser.print_help()
        print("\nERROR: Must specify either protein_id or --batch")
        sys.exit(1)

if __name__ == "__main__":
    main()
