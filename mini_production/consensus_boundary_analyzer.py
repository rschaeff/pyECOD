#!/usr/bin/env python3
"""
Consensus Boundary Analysis for Domain Architectures

Analyzes boundary consistency across proteins with identical architectures.
Identifies cases where boundaries should be spatially/sequentially aligned
but are inconsistent across similar proteins.

Usage:
    python consensus_boundary_analyzer.py --architecture "11.1.1 → 11.1.1"
    python consensus_boundary_analyzer.py --analyze-all-architectures
    python consensus_boundary_analyzer.py --boundary-variance-analysis
"""

import os
import sys
import argparse
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass
from collections import defaultdict
import statistics
import json

@dataclass
class ProteinBoundary:
    """Boundary information for a single protein"""
    protein_id: str
    pdb_id: str
    chain_id: str
    domain_count: int
    boundary_positions: List[int]  # Residue numbers where domains end
    domain_ranges: List[str]
    domain_lengths: List[int]
    sequence_length: int
    
    @property
    def normalized_boundaries(self) -> List[float]:
        """Boundary positions as fraction of sequence length"""
        if self.sequence_length == 0:
            return []
        return [pos / self.sequence_length for pos in self.boundary_positions]

@dataclass
class ArchitectureBoundaries:
    """Boundary analysis for all proteins with same architecture"""
    architecture: str
    proteins: List[ProteinBoundary]
    
    @property
    def boundary_variance(self) -> List[float]:
        """Calculate variance in boundary positions across proteins"""
        if len(self.proteins) < 2:
            return []
        
        variances = []
        max_boundaries = max(len(p.normalized_boundaries) for p in self.proteins)
        
        for boundary_idx in range(max_boundaries):
            positions = []
            for protein in self.proteins:
                if boundary_idx < len(protein.normalized_boundaries):
                    positions.append(protein.normalized_boundaries[boundary_idx])
            
            if len(positions) >= 2:
                variances.append(statistics.variance(positions))
            else:
                variances.append(0.0)
        
        return variances
    
    @property
    def consensus_boundaries(self) -> List[float]:
        """Calculate consensus boundary positions"""
        if len(self.proteins) < 2:
            return []
        
        consensus = []
        max_boundaries = max(len(p.normalized_boundaries) for p in self.proteins)
        
        for boundary_idx in range(max_boundaries):
            positions = []
            for protein in self.proteins:
                if boundary_idx < len(protein.normalized_boundaries):
                    positions.append(protein.normalized_boundaries[boundary_idx])
            
            if positions:
                consensus.append(statistics.mean(positions))
        
        return consensus
    
    @property
    def inconsistency_score(self) -> float:
        """Overall inconsistency score (higher = more inconsistent)"""
        variances = self.boundary_variance
        if not variances:
            return 0.0
        return sum(variances) / len(variances)

class ConsensusBoundaryAnalyzer:
    """Analyze boundary consistency across proteins with identical architectures"""
    
    def __init__(self, batch_base: str = "/data/ecod/pdb_updates/batches"):
        self.batch_base = Path(batch_base)
    
    def find_mini_domain_files(self) -> List[Path]:
        """Find all mini domain XML files"""
        xml_files = []
        
        for batch_dir in self.batch_base.iterdir():
            if not batch_dir.is_dir():
                continue
                
            mini_domains_dir = batch_dir / "mini_domains"
            if mini_domains_dir.exists():
                xml_files.extend(mini_domains_dir.glob("*.mini.domains.xml"))
        
        return xml_files
    
    def parse_range_segments(self, range_str: str) -> List[Tuple[int, int]]:
        """Parse range string into segments"""
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
        except ValueError:
            pass
        
        return segments
    
    def extract_protein_boundaries(self, xml_file: Path) -> Optional[ProteinBoundary]:
        """Extract boundary information from a mini domain XML file"""
        
        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()
            
            # Extract protein info
            protein_id = xml_file.stem.replace('.mini.domains', '')
            pdb_id = protein_id.split('_')[0]
            chain_id = protein_id.split('_')[1] if '_' in protein_id else 'A'
            
            # Parse domains
            domain_elements = root.findall(".//domain")
            
            if len(domain_elements) < 2:
                return None  # Only analyze multi-domain proteins
            
            domain_ranges = []
            domain_lengths = []
            all_positions = set()
            
            for domain_elem in domain_elements:
                range_str = domain_elem.get('range', '')
                segments = self.parse_range_segments(range_str)
                domain_length = sum(end - start + 1 for start, end in segments)
                
                domain_ranges.append(range_str)
                domain_lengths.append(domain_length)
                
                # Collect all positions
                for start, end in segments:
                    all_positions.update(range(start, end + 1))
            
            if not all_positions:
                return None
            
            sequence_length = max(all_positions)
            
            # Calculate boundary positions (end of each domain except the last)
            boundary_positions = []
            for i in range(len(domain_ranges) - 1):
                segments = self.parse_range_segments(domain_ranges[i])
                if segments:
                    domain_end = max(end for start, end in segments)
                    boundary_positions.append(domain_end)
            
            return ProteinBoundary(
                protein_id=protein_id,
                pdb_id=pdb_id,
                chain_id=chain_id,
                domain_count=len(domain_elements),
                boundary_positions=boundary_positions,
                domain_ranges=domain_ranges,
                domain_lengths=domain_lengths,
                sequence_length=sequence_length
            )
            
        except Exception as e:
            print(f"Warning: Failed to parse {xml_file}: {e}")
            return None
    
    def get_architecture_string(self, xml_file: Path) -> str:
        """Get architecture string for a protein"""
        
        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()
            
            domain_elements = root.findall(".//domain")
            t_groups = []
            
            for domain_elem in domain_elements:
                t_group = domain_elem.get('t_group', '').strip()
                if t_group:
                    t_groups.append(t_group)
                else:
                    t_groups.append('unclassified')
            
            return ' → '.join(t_groups)
            
        except Exception:
            return 'unknown'
    
    def analyze_architecture_boundaries(self, architecture: str, max_proteins: int = 50) -> Optional[ArchitectureBoundaries]:
        """Analyze boundary consistency for a specific architecture"""
        
        print(f"🔍 Analyzing boundaries for architecture: {architecture}")
        
        xml_files = self.find_mini_domain_files()
        matching_proteins = []
        
        for xml_file in xml_files:
            if len(matching_proteins) >= max_proteins:
                break
                
            arch = self.get_architecture_string(xml_file)
            if arch == architecture:
                boundary_info = self.extract_protein_boundaries(xml_file)
                if boundary_info:
                    matching_proteins.append(boundary_info)
        
        if len(matching_proteins) < 2:
            print(f"❌ Need at least 2 proteins for consensus analysis, found {len(matching_proteins)}")
            return None
        
        print(f"✓ Found {len(matching_proteins)} proteins with architecture {architecture}")
        
        return ArchitectureBoundaries(
            architecture=architecture,
            proteins=matching_proteins
        )
    
    def print_boundary_analysis(self, arch_boundaries: ArchitectureBoundaries):
        """Print detailed boundary analysis"""
        
        print(f"\n📊 BOUNDARY CONSENSUS ANALYSIS")
        print("=" * 100)
        print(f"Architecture: {arch_boundaries.architecture}")
        print(f"Proteins analyzed: {len(arch_boundaries.proteins)}")
        print(f"Overall inconsistency score: {arch_boundaries.inconsistency_score:.4f}")
        
        # Print consensus boundaries
        consensus = arch_boundaries.consensus_boundaries
        variances = arch_boundaries.boundary_variance
        
        print(f"\n🎯 CONSENSUS BOUNDARY POSITIONS:")
        for i, (cons_pos, variance) in enumerate(zip(consensus, variances)):
            print(f"  Boundary {i+1}: {cons_pos:.3f} ± {variance:.4f} (position in sequence)")
        
        # Print individual protein boundaries
        print(f"\n📋 INDIVIDUAL PROTEIN BOUNDARIES:")
        print(f"{'Protein':<15} {'Seq Len':<8} {'Boundaries (residue)':<25} {'Boundaries (normalized)':<30} {'Deviation'}")
        print("-" * 100)
        
        for protein in arch_boundaries.proteins:
            boundaries_str = ', '.join(map(str, protein.boundary_positions))
            normalized_str = ', '.join(f"{pos:.3f}" for pos in protein.normalized_boundaries)
            
            # Calculate deviation from consensus
            if consensus and protein.normalized_boundaries:
                deviations = [abs(actual - expected) for actual, expected in zip(protein.normalized_boundaries, consensus)]
                avg_deviation = sum(deviations) / len(deviations) if deviations else 0.0
                deviation_str = f"{avg_deviation:.4f}"
            else:
                deviation_str = "N/A"
            
            print(f"{protein.protein_id:<15} {protein.sequence_length:<8} {boundaries_str:<25} {normalized_str:<30} {deviation_str}")
        
        # Identify outliers
        print(f"\n⚠️  BOUNDARY OUTLIERS (high deviation from consensus):")
        outlier_threshold = 0.05  # 5% of sequence length
        
        outliers = []
        for protein in arch_boundaries.proteins:
            if consensus and protein.normalized_boundaries:
                deviations = [abs(actual - expected) for actual, expected in zip(protein.normalized_boundaries, consensus)]
                max_deviation = max(deviations) if deviations else 0.0
                if max_deviation > outlier_threshold:
                    outliers.append((protein, max_deviation))
        
        if outliers:
            outliers.sort(key=lambda x: x[1], reverse=True)
            for protein, deviation in outliers[:5]:  # Show top 5 outliers
                print(f"  {protein.protein_id}: max deviation {deviation:.4f}")
                print(f"    Expected: {[f'{p:.3f}' for p in consensus]}")
                print(f"    Actual:   {[f'{p:.3f}' for p in protein.normalized_boundaries]}")
        else:
            print("  No significant outliers detected")
    
    def analyze_all_two_domain_architectures(self, min_proteins: int = 5) -> Dict[str, ArchitectureBoundaries]:
        """Analyze boundary consistency for all two-domain architectures"""
        
        print(f"🔍 Finding all two-domain architectures with ≥{min_proteins} examples...")
        
        # First pass: count architectures
        architecture_counts = defaultdict(int)
        xml_files = self.find_mini_domain_files()
        
        for xml_file in xml_files:
            arch = self.get_architecture_string(xml_file)
            if ' → ' in arch and arch.count(' → ') == 1:  # Two domains
                architecture_counts[arch] += 1
        
        # Filter to architectures with enough examples
        candidate_architectures = {arch: count for arch, count in architecture_counts.items() 
                                 if count >= min_proteins}
        
        print(f"✓ Found {len(candidate_architectures)} two-domain architectures with ≥{min_proteins} examples")
        
        # Analyze each architecture
        results = {}
        for arch in sorted(candidate_architectures.keys(), key=lambda x: candidate_architectures[x], reverse=True)[:10]:
            result = self.analyze_architecture_boundaries(arch, max_proteins=20)
            if result:
                results[arch] = result
        
        return results
    
    def export_consensus_boundaries(self, results: Dict[str, ArchitectureBoundaries], output_file: str):
        """Export consensus boundary data for template-based boundary assignment"""
        
        export_data = {}
        
        for arch, arch_boundaries in results.items():
            consensus = arch_boundaries.consensus_boundaries
            variances = arch_boundaries.boundary_variance
            
            export_data[arch] = {
                'consensus_boundaries': consensus,
                'boundary_variances': variances,
                'inconsistency_score': arch_boundaries.inconsistency_score,
                'protein_count': len(arch_boundaries.proteins),
                'example_proteins': [p.protein_id for p in arch_boundaries.proteins[:5]]
            }
        
        with open(output_file, 'w') as f:
            json.dump(export_data, f, indent=2)
        
        print(f"✓ Exported consensus boundaries to {output_file}")


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Analyze Boundary Consensus Across Identical Domain Architectures'
    )
    
    parser.add_argument('--architecture', type=str,
                       help='Specific architecture to analyze (e.g., "11.1.1 → 11.1.1")')
    parser.add_argument('--analyze-all-architectures', action='store_true',
                       help='Analyze all two-domain architectures')
    parser.add_argument('--boundary-variance-analysis', action='store_true',
                       help='Focus on architectures with high boundary variance')
    parser.add_argument('--min-proteins', type=int, default=5,
                       help='Minimum proteins required for analysis')
    parser.add_argument('--export-consensus', type=str,
                       help='Export consensus boundaries to JSON file')
    
    args = parser.parse_args()
    
    analyzer = ConsensusBoundaryAnalyzer()
    
    if args.architecture:
        # Analyze specific architecture
        result = analyzer.analyze_architecture_boundaries(args.architecture)
        if result:
            analyzer.print_boundary_analysis(result)
        
    elif args.analyze_all_architectures:
        # Analyze all two-domain architectures
        results = analyzer.analyze_all_two_domain_architectures(args.min_proteins)
        
        if results:
            print(f"\n🏆 BOUNDARY CONSISTENCY SUMMARY")
            print("=" * 80)
            print(f"{'Architecture':<50} {'Proteins':<10} {'Inconsistency':<12} {'Status'}")
            print("-" * 80)
            
            for arch, arch_boundaries in sorted(results.items(), 
                                              key=lambda x: x[1].inconsistency_score, 
                                              reverse=True):
                protein_count = len(arch_boundaries.proteins)
                inconsistency = arch_boundaries.inconsistency_score
                
                if inconsistency > 0.01:
                    status = "⚠️  High variance"
                elif inconsistency > 0.005:
                    status = "⚡ Medium variance"
                else:
                    status = "✓ Consistent"
                
                print(f"{arch:<50} {protein_count:<10} {inconsistency:<12.4f} {status}")
            
            # Show detailed analysis for most problematic
            if results:
                worst_arch = max(results.items(), key=lambda x: x[1].inconsistency_score)
                print(f"\n🔍 DETAILED ANALYSIS - MOST INCONSISTENT ARCHITECTURE:")
                analyzer.print_boundary_analysis(worst_arch[1])
            
            if args.export_consensus:
                analyzer.export_consensus_boundaries(results, args.export_consensus)
        
    else:
        # Default: analyze the double Ig domain architecture
        result = analyzer.analyze_architecture_boundaries("11.1.1 → 11.1.1")
        if result:
            analyzer.print_boundary_analysis(result)
    
    print(f"\n💡 INSIGHTS:")
    print(f"High boundary variance indicates inconsistent domain splitting")
    print(f"Proteins with high deviation may need boundary refinement")
    print(f"Consensus boundaries can be used as templates for consistent assignment")


if __name__ == "__main__":
    main()
