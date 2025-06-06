#!/usr/bin/env python3
"""
Domain Architecture Analysis Tool for Mini PyECOD Production Results

Analyzes the distribution of domain architectures (N-C ordering of T-groups) 
across quality tiers to identify the most successful structural patterns.
Generates PyMOL superposition scripts for structural validation.

Usage:
    python analyze_domain_architectures.py --analyze-tiers
    python analyze_domain_architectures.py --architecture-examples "1.10.8 → 1.25.40"
    python analyze_domain_architectures.py --superposition-candidates 10
    python analyze_domain_architectures.py --generate-superposition "1.10.8 → 1.25.40"
"""

import os
import sys
import argparse
import yaml
import xml.etree.ElementTree as ET
import gzip
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass, asdict
from collections import defaultdict, Counter
import json
import logging
import re

# Import coordinate translation from visualization module
sys.path.insert(0, str(Path(__file__).parent / ".." / "mini" / "core"))
from visualization import CoordinateTranslator

logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class DomainInfo:
    """Domain information parsed from XML"""
    id: str
    range_str: str
    family: str
    t_group: Optional[str]
    h_group: Optional[str]
    source: str
    segments: List[Tuple[int, int]]  # (start, end) pairs

@dataclass
class ProteinResult:
    """Protein result with domains and quality assessment"""
    protein_id: str
    pdb_id: str
    chain_id: str
    batch_name: str
    domains: List[DomainInfo]
    tier: str
    coverage_percentage: float
    total_domains: int
    classified_domains: int
    
    @property
    def is_two_domain(self) -> bool:
        """Check if this is a two-domain protein"""
        return self.total_domains == 2
    
    @property
    def architecture_string(self) -> str:
        """Get architecture as string"""
        if not self.domains:
            return "no_domains"
        
        t_groups = []
        for domain in self.domains:
            if domain.t_group and domain.t_group.strip():
                t_groups.append(domain.t_group.strip())
            else:
                t_groups.append("unclassified")
        
        return " → ".join(t_groups)

@dataclass
class DomainArchitecture:
    """Domain architecture with T-group sequence"""
    t_group_sequence: Tuple[str, ...]  # N to C terminal order
    architecture_string: str  # Human-readable representation
    domain_count: int
    examples: List[ProteinResult]  # Protein results with this architecture
    
    @property
    def complexity(self) -> str:
        """Classify architecture complexity"""
        if self.domain_count == 1:
            return "single_domain"
        elif self.domain_count == 2:
            return "two_domain"
        elif self.domain_count <= 4:
            return "multi_domain"
        else:
            return "complex"

class PyMOLSuperpositionGenerator:
    """Generate PyMOL scripts for structural superposition of similar architectures"""
    
    def __init__(self, pdb_repo_path: str = "/usr2/pdb/data"):
        self.pdb_repo_path = pdb_repo_path
        self.colors = [
            "red", "blue", "green", "cyan", "magenta", "orange",
            "salmon", "lightblue", "lightgreen", "lightcyan", "pink", "gold"
        ]
    
    def find_structure(self, pdb_id: str) -> Optional[str]:
        """Find structure file in PDB repository"""
        pdb_id = pdb_id.lower()
        
        possible_paths = [
            f"{self.pdb_repo_path}/structures/divided/mmCIF/{pdb_id[1:3]}/{pdb_id}.cif.gz",
            f"{self.pdb_repo_path}/structures/divided/mmCIF/{pdb_id[1:3]}/{pdb_id}.cif",
            f"{self.pdb_repo_path}/mmCIF/{pdb_id}.cif.gz",
            f"{self.pdb_repo_path}/mmCIF/{pdb_id}.cif",
        ]
        
        for path in possible_paths:
            if os.path.exists(path):
                return path
        return None
    
    def parse_range_segments(self, range_str: str) -> List[Tuple[int, int]]:
        """Parse range string into list of (start, end) tuples"""
        segments = []
        if not range_str:
            return segments
        
        try:
            for segment in range_str.split(','):
                segment = segment.strip()
                if '-' in segment:
                    start, end = segment.split('-')
                    segments.append((int(start), int(end)))
                else:
                    pos = int(segment)
                    segments.append((pos, pos))
        except ValueError as e:
            logger.warning(f"Could not parse range '{range_str}': {e}")
        
        return segments
    
    def generate_superposition_script(self, architecture: DomainArchitecture, 
                                    max_structures: int = 10,
                                    output_dir: str = "/tmp/pymol_superposition") -> str:
        """Generate PyMOL script for superposing proteins with identical architecture"""
        
        # Limit to max_structures for performance
        proteins = architecture.examples[:max_structures]
        
        if len(proteins) < 2:
            raise ValueError(f"Need at least 2 proteins for superposition, got {len(proteins)}")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Generate script filename
        arch_clean = re.sub(r'[^a-zA-Z0-9._-]', '_', architecture.architecture_string)
        script_path = os.path.join(output_dir, f"superposition_{arch_clean}.pml")
        
        lines = [
            f"# PyMOL Superposition Script",
            f"# Architecture: {architecture.architecture_string}",
            f"# Domain count: {architecture.domain_count}",
            f"# Proteins: {len(proteins)}",
            f"# Generated for structural validation",
            "",
            "# Clear workspace",
            "delete all",
            "bg_color white",
            "",
            "# Load all structures"
        ]
        
        # Load structures and validate availability
        valid_proteins = []
        for i, protein in enumerate(proteins):
            pdb_id = protein.pdb_id
            chain_id = protein.chain_id
            structure_path = self.find_structure(pdb_id)
            
            if structure_path:
                object_name = f"{pdb_id}_{chain_id}"
                lines.extend([
                    f"# Protein {i+1}: {protein.protein_id}",
                    f"load {structure_path}, {object_name}",
                    f"# Extract chain {chain_id}",
                    f"create {object_name}_chain, {object_name} and chain {chain_id}",
                    f"delete {object_name}",
                    ""
                ])
                valid_proteins.append(protein)
            else:
                lines.append(f"# WARNING: Structure not found for {pdb_id}")
        
        if len(valid_proteins) < 2:
            raise FileNotFoundError(f"Insufficient structures found for superposition")
        
        # Set up reference structure (first one)
        ref_protein = valid_proteins[0]
        ref_name = f"{ref_protein.pdb_id}_{ref_protein.chain_id}_chain"
        
        lines.extend([
            f"# Reference structure: {ref_protein.protein_id}",
            f"show cartoon, {ref_name}",
            "",
            "# Superpose all other structures to reference"
        ])

        # Superpose each structure
        for i, protein in enumerate(valid_proteins[1:], 1):
            target_name = f"{protein.pdb_id}_{protein.chain_id}_chain"

            lines.extend([
                f"# Superpose {protein.protein_id}",
                f"align {target_name}, {ref_name}",
                f"show cartoon, {target_name}",
                ""
            ])

        # Color domains by architecture
        lines.extend([
            "",
            "# Color domains by architecture position",
            "# This helps validate domain boundary consistency"
        ])

        # Color everything gray80 first to show unassigned regions
        lines.extend([
            "",
            "# Color all structures gray80 first (shows unassigned regions)",
            "color gray80, all",
            ""
        ])
        
        # Color each domain position across all structures
        domain_colors = ["red", "blue", "green", "yellow", "magenta"][:architecture.domain_count]
        
        for domain_idx in range(architecture.domain_count):
            domain_color = domain_colors[domain_idx]
            t_group = architecture.t_group_sequence[domain_idx]
            
            lines.append(f"# Domain {domain_idx + 1}: {t_group}")
            
            for protein in valid_proteins:
                if domain_idx < len(protein.domains):
                    domain = protein.domains[domain_idx]
                    object_name = f"{protein.pdb_id}_{protein.chain_id}_chain"
                    
                    # Parse domain range and create selection
                    segments = self.parse_range_segments(domain.range_str)
                    for start, end in segments:
                        if start == end:
                            selection = f"{object_name} and resi {start}"
                        else:
                            selection = f"{object_name} and resi {start}-{end}"
                        lines.append(f"color {domain_color}, {selection}")
            
            lines.append("")
        
        # Final setup
        lines.extend([
            "# Final visualization setup",
            "set cartoon_fancy_helices, 1",
            "set cartoon_fancy_sheets, 1",
            "orient",
            "zoom all",
            "",
            "# Create selection for each domain across all structures",
        ])
        
        for domain_idx in range(architecture.domain_count):
            t_group = architecture.t_group_sequence[domain_idx]
            selection_parts = []
            
            for protein in valid_proteins:
                if domain_idx < len(protein.domains):
                    domain = protein.domains[domain_idx]
                    object_name = f"{protein.pdb_id}_{protein.chain_id}_chain"
                    segments = self.parse_range_segments(domain.range_str)
                    
                    for start, end in segments:
                        if start == end:
                            selection_parts.append(f"({object_name} and resi {start})")
                        else:
                            selection_parts.append(f"({object_name} and resi {start}-{end})")
            
            if selection_parts:
                selection_name = f"domain_{domain_idx + 1}_{t_group.replace('.', '_')}"
                selection_expr = " or ".join(selection_parts)
                lines.extend([
                    f"select {selection_name}, {selection_expr}",
                    f"# Use 'show sticks, {selection_name}' to highlight this domain"
                ])
        
        lines.extend([
            "",
            f"# Save superposition session",
            f"save {output_dir}/superposition_{arch_clean}.pse",
            "",
            "# Analysis commands:",
            "# rmsd_cur - show RMSD between structures",
            "# select domain_X_Y, [domain selection] - select specific domains",
            "# show sticks, [selection] - highlight specific regions",
            f"# Architecture validated: {architecture.architecture_string}",
            f"# Total structures: {len(valid_proteins)}",
        ])
        
        # Write script
        with open(script_path, 'w') as f:
            f.write('\n'.join(lines))
        
        logger.info(f"Generated superposition script: {script_path}")
        logger.info(f"  Architecture: {architecture.architecture_string}")
        logger.info(f"  Structures: {len(valid_proteins)}")
        
        return script_path
    class DomainArchitectureAnalyzer:
    """Analyze domain architectures from mini PyECOD XML results"""
    
    def __init__(self, config_path: str = "config/config.local.yml"):
        self.config = self._load_config(config_path)
        self.pymol_generator = PyMOLSuperpositionGenerator()
        
    def _load_config(self, config_path: str) -> Dict:
        """Load configuration"""
        try:
            with open(config_path, 'r') as f:
                return yaml.safe_load(f)
        except FileNotFoundError:
            return {"paths": {"batch_base_dir": "/data/ecod/pdb_updates/batches"}}
    
    def parse_domain_xml(self, xml_path: Path) -> Optional[ProteinResult]:
        """Parse mini domain XML file and assess quality"""
        
        try:
            tree = ET.parse(xml_path)
            root = tree.getroot()
            
            # Extract protein info
            protein_id = xml_path.stem.replace('.mini.domains', '')
            pdb_id = protein_id.split('_')[0]
            chain_id = protein_id.split('_')[1] if '_' in protein_id else 'A'
            batch_name = xml_path.parent.parent.name
            
            # Parse domains
            domains = []
            domain_elements = root.findall(".//domain")
            
            total_length = 0
            classified_count = 0
            
            for domain_elem in domain_elements:
                domain_id = domain_elem.get('id', 'unknown')
                range_str = domain_elem.get('range', '')
                family = domain_elem.get('family', 'unknown')
                t_group = domain_elem.get('t_group')
                h_group = domain_elem.get('h_group')
                source = domain_elem.get('source', 'unknown')
                
                # Parse range segments
                segments = self.pymol_generator.parse_range_segments(range_str)
                segment_length = sum(end - start + 1 for start, end in segments)
                total_length += segment_length
                
                if t_group and t_group.strip():
                    classified_count += 1
                
                domain_info = DomainInfo(
                    id=domain_id,
                    range_str=range_str,
                    family=family,
                    t_group=t_group,
                    h_group=h_group,
                    source=source,
                    segments=segments
                )
                domains.append(domain_info)
            
            # Estimate sequence length (rough approximation)
            if domains:
                max_end = max(max(end for start, end in domain.segments) 
                            for domain in domains if domain.segments)
                sequence_length = int(max_end * 1.1)  # Add 10% buffer
            else:
                sequence_length = 1
            
            # Calculate coverage and quality tier
            coverage_percentage = (total_length / sequence_length) * 100 if sequence_length > 0 else 0
            total_domains = len(domains)
            
            # Simplified tier assignment (focusing on excellent tier)
            if (total_domains > 0 and 
                classified_count >= total_domains * 0.8 and 
                coverage_percentage >= 60):
                tier = "excellent"
            elif (total_domains > 0 and 
                  classified_count >= total_domains * 0.6 and 
                  coverage_percentage >= 40):
                tier = "good"
            elif total_domains > 0 and coverage_percentage >= 20:
                tier = "acceptable"
            elif total_domains > 0:
                tier = "poor"
            else:
                tier = "failed"
            
            return ProteinResult(
                protein_id=protein_id,
                pdb_id=pdb_id,
                chain_id=chain_id,
                batch_name=batch_name,
                domains=domains,
                tier=tier,
                coverage_percentage=coverage_percentage,
                total_domains=total_domains,
                classified_domains=classified_count
            )
            
        except Exception as e:
            logger.warning(f"Failed to parse {xml_path}: {e}")
            return None
    
    def scan_mini_results(self, batch_filter: Optional[str] = None) -> List[ProteinResult]:
        """Scan filesystem for mini domain XML results"""
        
        batch_base = Path(self.config["paths"]["batch_base_dir"])
        
        if not batch_base.exists():
            raise FileNotFoundError(f"Batch base directory not found: {batch_base}")
        
        print(f"🔍 Scanning mini results in {batch_base}")
        
        results = []
        
        # Find batch directories
        if batch_filter:
            batch_dirs = [batch_base / batch_filter]
        else:
            batch_dirs = [d for d in batch_base.iterdir() if d.is_dir()]
        
        for batch_dir in batch_dirs:
            mini_domains_dir = batch_dir / "mini_domains"
            
            if not mini_domains_dir.exists():
                continue
            
            xml_files = list(mini_domains_dir.glob("*.mini.domains.xml"))
            print(f"📁 Batch {batch_dir.name}: {len(xml_files)} results")
            
            for xml_file in xml_files:
                result = self.parse_domain_xml(xml_file)
                if result:
                    results.append(result)
        
        print(f"✓ Loaded {len(results)} protein results")
        return results
    
    def get_two_domain_architectures(self, results: List[ProteinResult], 
                                   tier_filter: str = "excellent") -> Dict[str, DomainArchitecture]:
        """Get two-domain architectures from specified tier"""
        
        # Filter to two-domain proteins in excellent tier
        two_domain_results = [
            r for r in results 
            if r.is_two_domain and r.tier == tier_filter
        ]
        
        print(f"🎯 Found {len(two_domain_results)} two-domain {tier_filter} proteins")
        
        # Group by architecture
        architecture_map = defaultdict(list)
        
        for result in two_domain_results:
            arch_str = result.architecture_string
            architecture_map[arch_str].append(result)
        
        # Create DomainArchitecture objects
        architectures = {}
        for arch_str, protein_list in architecture_map.items():
            if len(protein_list) >= 2:  # Need at least 2 for superposition
                first_protein = protein_list[0]
                t_groups = tuple(d.t_group or "unclassified" for d in first_protein.domains)
                
                architecture = DomainArchitecture(
                    t_group_sequence=t_groups,
                    architecture_string=arch_str,
                    domain_count=2,
                    examples=protein_list
                )
                architectures[arch_str] = architecture
        
        # Sort by frequency
        sorted_architectures = dict(sorted(architectures.items(), 
                                         key=lambda x: len(x[1].examples), 
                                         reverse=True))
        
        print(f"✓ Found {len(sorted_architectures)} two-domain architectures with ≥2 examples")
        
        return sorted_architectures
    
    def analyze_tier_architectures(self, results: List[ProteinResult], tier: str = "excellent") -> Dict[str, DomainArchitecture]:
        """Analyze all architectures in a specific tier"""
        
        # Filter to specified tier
        tier_results = [r for r in results if r.tier == tier]
        
        print(f"🎯 Found {len(tier_results)} {tier} proteins")
        
        # Group by architecture
        architecture_map = defaultdict(list)
        
        for result in tier_results:
            arch_str = result.architecture_string
            architecture_map[arch_str].append(result)
        
        # Create DomainArchitecture objects
        architectures = {}
        for arch_str, protein_list in architecture_map.items():
            first_protein = protein_list[0]
            if first_protein.domains:
                t_groups = tuple(d.t_group or "unclassified" for d in first_protein.domains)
                
                architecture = DomainArchitecture(
                    t_group_sequence=t_groups,
                    architecture_string=arch_str,
                    domain_count=len(first_protein.domains),
                    examples=protein_list
                )
                architectures[arch_str] = architecture
        
        # Sort by frequency
        sorted_architectures = dict(sorted(architectures.items(), 
                                         key=lambda x: len(x[1].examples), 
                                         reverse=True))
        
        print(f"✓ Found {len(sorted_architectures)} unique architectures in {tier} tier")
        
        return sorted_architectures


    def print_two_domain_summary(self, architectures: Dict[str, DomainArchitecture]):
        """Print two-domain architecture summary"""
        
        print(f"\n🏗️  TWO-DOMAIN ARCHITECTURE ANALYSIS (EXCELLENT TIER)")
        print("=" * 80)
        
        if not architectures:
            print("No two-domain architectures found in excellent tier")
            return
        
        total_proteins = sum(len(arch.examples) for arch in architectures.values())
        print(f"Total two-domain excellent proteins: {total_proteins:,}")
        print(f"Unique two-domain architectures: {len(architectures):,}")
        print()
        
        print(f"📈 Two-Domain Architectures (sorted by frequency):")
        print(f"{'Rank':<6} {'Architecture':<50} {'Count':<8} {'Examples'}")
        print("-" * 90)
        
        for i, (arch_str, arch) in enumerate(architectures.items(), 1):
            count = len(arch.examples)
            example_ids = [p.protein_id for p in arch.examples[:3]]
            examples_str = ", ".join(example_ids)
            if len(arch.examples) > 3:
                examples_str += f" ... (+{len(arch.examples) - 3} more)"
            
            print(f"{i:<6} {arch_str:<50} {count:<8} {examples_str}")
            
            if i >= 20:  # Show top 20
                remaining = len(architectures) - 20
                if remaining > 0:
                    print(f"... and {remaining} more architectures")
                break


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Analyze Two-Domain Architectures in Mini PyECOD Production Results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze all architectures in excellent tier
  python analyze_domain_architectures.py --analyze-excellent
  
  # Focus on two-domain architectures in excellent tier
  python analyze_domain_architectures.py --two-domain-analysis
  
  # Generate superposition script for specific architecture
  python analyze_domain_architectures.py --generate-superposition "1.10.8 → 1.25.40"
  
  # Find examples of specific architecture
  python analyze_domain_architectures.py --architecture-examples "1.10.8 → 1.25.40"
  
  # Get top candidates for superposition
  python analyze_domain_architectures.py --superposition-candidates 10
        """
    )
    
    parser.add_argument('--two-domain-analysis', action='store_true',
                       help='Analyze two-domain architectures in excellent tier')
    parser.add_argument('--analyze-excellent', action='store_true',
                       help='Analyze all architectures in excellent tier')
    parser.add_argument('--architecture-examples', type=str,
                       help='Find examples of specific architecture (e.g., "1.10.8 → 1.25.40")')
    parser.add_argument('--generate-superposition', type=str,
                       help='Generate PyMOL superposition script for architecture')
    parser.add_argument('--superposition-candidates', type=int, metavar='N',
                       help='Show top N two-domain architectures for superposition')
    parser.add_argument('--batch-filter', type=str,
                       help='Filter to specific batch (e.g., "ecod_batch_036_20250406_1424")')
    parser.add_argument('--max-structures', type=int, default=10,
                       help='Maximum structures per superposition script')
    parser.add_argument('--min-examples', type=int, default=2,
                       help='Minimum examples required for superposition')
    parser.add_argument('--output-dir', type=str, default='/tmp/pymol_superposition',
                       help='Output directory for PyMOL scripts')
    parser.add_argument('--config', type=str, default='config/config.local.yml',
                       help='Config file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Enable verbose output')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Initialize analyzer
    analyzer = DomainArchitectureAnalyzer(args.config)
    
    # Load all results
    print("🔍 Scanning filesystem for mini PyECOD results...")
    results = analyzer.scan_mini_results(args.batch_filter)
    
    if not results:
        print("❌ No results found to analyze")
        return
    
    print(f"✓ Loaded {len(results):,} protein results")
    
    # Get two-domain architectures from excellent tier
    two_domain_architectures = analyzer.get_two_domain_architectures(results, "excellent")
    
    if args.two_domain_analysis:
        analyzer.print_two_domain_summary(two_domain_architectures)
    
    if args.analyze_excellent:
        print(f"\n🏆 EXCELLENT TIER ARCHITECTURE ANALYSIS")
        print("=" * 80)
        
        excellent_architectures = analyzer.analyze_tier_architectures(results, "excellent")
        
        if excellent_architectures:
            total_excellent = sum(len(arch.examples) for arch in excellent_architectures.values())
            print(f"Total excellent proteins: {total_excellent:,}")
            print(f"Unique architectures: {len(excellent_architectures):,}")
            print()
            
            # Group by complexity
            by_complexity = defaultdict(list)
            for arch_str, arch in excellent_architectures.items():
                by_complexity[arch.complexity].append((arch_str, arch))
            
            for complexity in ["single_domain", "two_domain", "multi_domain", "complex"]:
                if complexity not in by_complexity:
                    continue
                
                complexity_archs = by_complexity[complexity]
                complexity_total = sum(len(arch[1].examples) for arch in complexity_archs)
                
                print(f"📊 {complexity.replace('_', ' ').title()} ({len(complexity_archs)} architectures, {complexity_total} proteins):")
                
                # Show top 5 in this complexity
                sorted_complexity = sorted(complexity_archs, key=lambda x: len(x[1].examples), reverse=True)
                for i, (arch_str, arch) in enumerate(sorted_complexity[:5], 1):
                    count = len(arch.examples)
                    pct = count / total_excellent * 100
                    print(f"  {i}. {arch_str:<45} {count:>6} ({pct:>5.1f}%)")
                
                if len(complexity_archs) > 5:
                    print(f"     ... and {len(complexity_archs) - 5} more")
                print()
        else:
            print("No architectures found in excellent tier")
    
    if args.architecture_examples:
        print(f"\n🔍 EXAMPLES OF ARCHITECTURE: {args.architecture_examples}")
        print("=" * 80)
        
        # First check in two-domain architectures
        found_examples = []
        source = "two-domain excellent"
        
        if args.architecture_examples in two_domain_architectures:
            arch = two_domain_architectures[args.architecture_examples]
            found_examples = arch.examples
        else:
            # Check in all excellent tier architectures
            excellent_architectures = analyzer.analyze_tier_architectures(results, "excellent")
            if args.architecture_examples in excellent_architectures:
                arch = excellent_architectures[args.architecture_examples]
                found_examples = arch.examples
                source = "excellent tier"
        
        if found_examples:
            examples = found_examples[:20]  # Show up to 20 examples
            
            print(f"Found {len(found_examples)} total examples in {source} (showing first {len(examples)}):")
            print(f"{'#':<4} {'Protein ID':<15} {'PDB':<6} {'Chain':<6} {'Coverage':<10} {'Domains'}")
            print("-" * 60)
            
            for i, example in enumerate(examples, 1):
                print(f"{i:<4} {example.protein_id:<15} {example.pdb_id:<6} {example.chain_id:<6} "
                      f"{example.coverage_percentage:>6.1f}%   {example.total_domains}")
        else:
            print(f"Architecture '{args.architecture_examples}' not found in excellent tier")
            print("Available two-domain architectures:")
            for arch_str in list(two_domain_architectures.keys())[:10]:
                print(f"  {arch_str}")
            
            if len(two_domain_architectures) > 10:
                print(f"  ... and {len(two_domain_architectures) - 10} more")
    
    if args.generate_superposition:
        print(f"\n🧬 GENERATING SUPERPOSITION SCRIPT: {args.generate_superposition}")
        print("=" * 80)
        
        # First check in two-domain architectures
        target_arch = None
        if args.generate_superposition in two_domain_architectures:
            target_arch = two_domain_architectures[args.generate_superposition]
        else:
            # Check in all excellent tier architectures
            excellent_architectures = analyzer.analyze_tier_architectures(results, "excellent")
            if args.generate_superposition in excellent_architectures:
                target_arch = excellent_architectures[args.generate_superposition]
        
        if target_arch:
            try:
                script_path = analyzer.pymol_generator.generate_superposition_script(
                    target_arch, args.max_structures, args.output_dir)
                
                print(f"✅ PyMOL superposition script generated!")
                print(f"   Script: {script_path}")
                print(f"   Architecture: {target_arch.architecture_string}")
                print(f"   Structures: {min(len(target_arch.examples), args.max_structures)}")
                print(f"\nTo run:")
                print(f"   pymol {script_path}")
                
            except Exception as e:
                print(f"❌ Failed to generate script: {e}")
        else:
            print(f"Architecture '{args.generate_superposition}' not found in excellent tier")
    
    if args.superposition_candidates:
        print(f"\n🎯 TOP {args.superposition_candidates} SUPERPOSITION CANDIDATES")
        print("=" * 80)
        
        # Get all excellent tier architectures for candidates
        excellent_architectures = analyzer.analyze_tier_architectures(results, "excellent")
        
        # Filter architectures with sufficient examples
        candidates = {
            arch_str: arch for arch_str, arch in excellent_architectures.items()
            if len(arch.examples) >= args.min_examples
        }
        
        if not candidates:
            print(f"No architectures found with ≥{args.min_examples} examples")
            return
        
        # Sort by example count
        sorted_candidates = sorted(candidates.items(), 
                                 key=lambda x: len(x[1].examples), 
                                 reverse=True)
        
        print(f"{'Rank':<6} {'Architecture':<50} {'Examples':<10} {'Domains':<8} {'T-groups'}")
        print("-" * 100)
        
        for i, (arch_str, arch) in enumerate(sorted_candidates[:args.superposition_candidates], 1):
            count = len(arch.examples)
            domain_count = arch.domain_count
            t_groups = " + ".join(arch.t_group_sequence)
            print(f"{i:<6} {arch_str:<50} {count:<10} {domain_count:<8} {t_groups}")
        
        print(f"\nTo generate superposition script:")
        print(f'python {sys.argv[0]} --generate-superposition "ARCHITECTURE_STRING"')
    
    # Print overall summary
    excellent_count = len([r for r in results if r.tier == "excellent"])
    two_domain_excellent = len([r for r in results if r.tier == "excellent" and r.is_two_domain])
    
    print(f"\n📊 OVERALL SUMMARY")
    print("=" * 50)
    print(f"Total proteins analyzed:        {len(results):,}")
    print(f"Excellent tier:                 {excellent_count:,}")
    print(f"Two-domain excellent:           {two_domain_excellent:,}")
    print(f"Unique two-domain architectures: {len(two_domain_architectures):,}")
    
    if two_domain_architectures:
        most_common = max(two_domain_architectures.values(), key=lambda x: len(x.examples))
        print(f"Most common architecture:       {most_common.architecture_string} ({len(most_common.examples)} examples)")


if __name__ == "__main__":
    main()
