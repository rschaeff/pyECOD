#!/usr/bin/env python3
"""
Analyze Anomalous Double Ig Domain Predictions - OVERSEGMENTATION & UNRESOLVED REGIONS

Investigates why certain proteins classified as "11.1.1 → 11.1.1" 
are structurally anomalous, focusing on:
1. Oversegmentation of Ig domains (single domain split into pieces)
2. Domains assigned to unresolved/missing structural regions

Usage:
    python anomaly_analyzer.py --analyze-proteins "8g8c_H,8g8d_H,8gav_H"
    python anomaly_analyzer.py --oversegmentation-analysis
    python anomaly_analyzer.py --structure-resolution-check
"""

import os
import sys
import argparse
import xml.etree.ElementTree as ET
import gzip
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass
import re

# Add the domain architecture analyzer
sys.path.insert(0, str(Path(__file__).parent))

@dataclass
class StructuralResolution:
    """Track which residues are resolved in the crystal structure"""
    resolved_residues: Set[int]
    missing_residues: Set[int]
    total_sequence_length: int
    resolution_percentage: float
    
@dataclass
class DomainAnalysis:
    """Enhanced domain analysis for oversegmentation and resolution issues"""
    protein_id: str
    pdb_id: str
    chain_id: str
    total_domains: int
    domain_ranges: List[str]
    domain_lengths: List[int]
    domain_positions: List[Set[int]]  # Actual positions for each domain
    sequence_length: int
    total_coverage: int
    coverage_percentage: float
    gaps: List[Tuple[int, int]]  # (start, end) of gaps
    largest_gap: int
    domain_sources: List[str]
    evidence_types: List[str]
    
    # NEW: Oversegmentation analysis
    undersized_domains: List[int]  # Domain indices that are too small
    proximal_domains: List[Tuple[int, int, int]]  # (domain1_idx, domain2_idx, gap_size)
    
    # NEW: Structural resolution analysis
    structural_resolution: Optional[StructuralResolution] = None
    domains_in_unresolved: List[Tuple[int, Set[int]]] = None  # (domain_idx, unresolved_positions)
    
    @property
    def is_antibody_like(self) -> bool:
        """Check if this looks like an antibody chain"""
        return (self.chain_id in ['H', 'L'] or 
                self.protein_id.endswith('_H') or 
                self.protein_id.endswith('_L'))
    
    @property
    def typical_ig_domain_size(self) -> Tuple[int, int]:
        """Typical Ig domain size range"""
        return (90, 130)  # Typical Ig domain is ~110 residues
    
    @property
    def oversegmentation_score(self) -> float:
        """Score for likely oversegmentation (higher = more likely oversegmented)"""
        score = 0.0
        min_size, max_size = self.typical_ig_domain_size
        
        # Count undersized domains
        undersized_count = len(self.undersized_domains)
        if undersized_count > 0:
            score += undersized_count * 2.0
        
        # Check for proximal domains (likely fragments of one domain)
        for domain1_idx, domain2_idx, gap_size in self.proximal_domains:
            if gap_size <= 20:  # Very close domains
                combined_size = self.domain_lengths[domain1_idx] + self.domain_lengths[domain2_idx] + gap_size
                if min_size <= combined_size <= max_size * 1.2:
                    score += 3.0  # Strong evidence of oversegmentation
        
        return score
    
    @property
    def unresolved_region_score(self) -> float:
        """Score for domains in unresolved regions"""
        if not self.domains_in_unresolved or not self.structural_resolution:
            return 0.0
        
        score = 0.0
        for domain_idx, unresolved_positions in self.domains_in_unresolved:
            domain_size = self.domain_lengths[domain_idx]
            unresolved_count = len(unresolved_positions)
            unresolved_fraction = unresolved_count / domain_size
            
            if unresolved_fraction > 0.3:  # >30% of domain is unresolved
                score += unresolved_fraction * 5.0
        
        return score

class StructuralAnalyzer:
    """Analyze structural resolution from mmCIF files"""
    
    def __init__(self, pdb_repo_path: str = "/usr2/pdb/data"):
        self.pdb_repo_path = pdb_repo_path
    
    def find_structure_file(self, pdb_id: str) -> Optional[str]:
        """Find mmCIF structure file"""
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
    
    def parse_resolved_residues(self, structure_path: str, chain_id: str) -> Optional[StructuralResolution]:
        """Parse which residues are resolved in the structure"""
        
        if structure_path.endswith('.gz'):
            opener = gzip.open
            mode = 'rt'
        else:
            opener = open
            mode = 'r'
        
        resolved_residues = set()
        sequence_length = 0
        
        try:
            with opener(structure_path, mode) as f:
                in_atom_site = False
                header_processed = False
                
                for line in f:
                    line = line.strip()
                    
                    # Look for atom_site section
                    if line.startswith('loop_'):
                        # Check if next few lines contain _atom_site
                        continue
                    elif line.startswith('_atom_site.'):
                        in_atom_site = True
                        continue
                    elif in_atom_site and line.startswith('_'):
                        continue
                    elif in_atom_site and (line.startswith('loop_') or line.startswith('#')):
                        break
                    elif in_atom_site and line:
                        # Parse atom site line
                        parts = line.split()
                        if len(parts) >= 8:  # Minimum required fields
                            try:
                                atom_chain_id = parts[6]  # auth_asym_id (chain)
                                auth_seq_id = parts[8]    # auth_seq_id (residue number)
                                
                                if atom_chain_id == chain_id:
                                    try:
                                        residue_num = int(auth_seq_id)
                                        resolved_residues.add(residue_num)
                                        sequence_length = max(sequence_length, residue_num)
                                    except ValueError:
                                        # Skip residues with insertion codes or non-numeric
                                        continue
                            except (IndexError, ValueError):
                                continue
                        
            # Calculate missing residues (gaps in sequence)
            if resolved_residues and sequence_length > 0:
                all_expected = set(range(1, sequence_length + 1))
                missing_residues = all_expected - resolved_residues
                resolution_percentage = len(resolved_residues) / sequence_length * 100
                
                return StructuralResolution(
                    resolved_residues=resolved_residues,
                    missing_residues=missing_residues,
                    total_sequence_length=sequence_length,
                    resolution_percentage=resolution_percentage
                )
            
        except Exception as e:
            print(f"Warning: Could not parse structure {structure_path}: {e}")
        
        return None

class AnomalyAnalyzer:
    """Analyze anomalous domain predictions with focus on oversegmentation and unresolved regions"""
    
    def __init__(self, config_path: str = "config/config.local.yml"):
        self.batch_base = Path("/data/ecod/pdb_updates/batches")
        self.structural_analyzer = StructuralAnalyzer()
        
    def find_protein_xml(self, protein_id: str) -> Optional[Path]:
        """Find the mini domain XML file for a protein"""
        
        # Search across all batch directories
        for batch_dir in self.batch_base.iterdir():
            if not batch_dir.is_dir():
                continue
                
            mini_domains_dir = batch_dir / "mini_domains"
            if not mini_domains_dir.exists():
                continue
            
            xml_file = mini_domains_dir / f"{protein_id}.mini.domains.xml"
            if xml_file.exists():
                return xml_file
        
        return None
    
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
        except ValueError as e:
            print(f"Warning: Could not parse range '{range_str}': {e}")
        
        return segments
    
    def analyze_oversegmentation(self, domain_ranges: List[str], domain_lengths: List[int]) -> Tuple[List[int], List[Tuple[int, int, int]]]:
        """Analyze potential oversegmentation"""
        
        min_size, max_size = 90, 130  # Typical Ig domain size
        
        # Find undersized domains
        undersized_domains = []
        for i, length in enumerate(domain_lengths):
            if length < min_size:
                undersized_domains.append(i)
        
        # Find proximal domains (close together, might be fragments)
        proximal_domains = []
        domain_segments = [self.parse_range_segments(range_str) for range_str in domain_ranges]
        
        for i in range(len(domain_segments) - 1):
            # Get end of domain i and start of domain i+1
            if domain_segments[i] and domain_segments[i+1]:
                domain_i_end = max(end for start, end in domain_segments[i])
                domain_j_start = min(start for start, end in domain_segments[i+1])
                
                gap_size = domain_j_start - domain_i_end - 1
                
                if gap_size <= 30:  # Domains very close together
                    proximal_domains.append((i, i+1, gap_size))
        
        return undersized_domains, proximal_domains
    
    def analyze_unresolved_regions(self, domain_ranges: List[str], structural_resolution: StructuralResolution) -> List[Tuple[int, Set[int]]]:
        """Find domains assigned to unresolved structural regions"""
        
        domains_in_unresolved = []
        
        for i, range_str in enumerate(domain_ranges):
            segments = self.parse_range_segments(range_str)
            domain_positions = set()
            
            for start, end in segments:
                domain_positions.update(range(start, end + 1))
            
            # Find overlap with unresolved regions
            unresolved_overlap = domain_positions.intersection(structural_resolution.missing_residues)
            
            if unresolved_overlap:
                domains_in_unresolved.append((i, unresolved_overlap))
        
        return domains_in_unresolved
        
    def find_protein_xml(self, protein_id: str) -> Optional[Path]:
        """Find the mini domain XML file for a protein"""
        
        # Search across all batch directories
        for batch_dir in self.batch_base.iterdir():
            if not batch_dir.is_dir():
                continue
                
            mini_domains_dir = batch_dir / "mini_domains"
            if not mini_domains_dir.exists():
                continue
            
            xml_file = mini_domains_dir / f"{protein_id}.mini.domains.xml"
            if xml_file.exists():
                return xml_file
        
        return None
    
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
        except ValueError as e:
            print(f"Warning: Could not parse range '{range_str}': {e}")
        
        return segments
    
    def analyze_protein(self, protein_id: str, include_structure: bool = True) -> Optional[DomainAnalysis]:
        """Analyze a single protein for oversegmentation and unresolved region issues"""
        
        xml_file = self.find_protein_xml(protein_id)
        if not xml_file:
            print(f"❌ XML file not found for {protein_id}")
            return None
        
        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()
            
            # Extract basic info
            pdb_id = protein_id.split('_')[0]
            chain_id = protein_id.split('_')[1] if '_' in protein_id else 'A'
            
            # Parse domains
            domain_elements = root.findall(".//domain")
            domain_ranges = []
            domain_lengths = []
            domain_positions = []
            domain_sources = []
            evidence_types = []
            all_positions = set()
            
            for domain_elem in domain_elements:
                range_str = domain_elem.get('range', '')
                source = domain_elem.get('source', 'unknown')
                
                segments = self.parse_range_segments(range_str)
                domain_length = sum(end - start + 1 for start, end in segments)
                
                # Get all positions for this domain
                positions = set()
                for start, end in segments:
                    positions.update(range(start, end + 1))
                
                domain_ranges.append(range_str)
                domain_lengths.append(domain_length)
                domain_positions.append(positions)
                domain_sources.append(source)
                evidence_types.append(source)
                
                # Collect all assigned positions
                all_positions.update(positions)
            
            # Calculate sequence metrics
            if all_positions:
                sequence_length = max(all_positions)
                total_coverage = len(all_positions)
                coverage_percentage = (total_coverage / sequence_length) * 100
                
                # Find gaps
                gaps = []
                if domain_ranges:
                    all_segments = []
                    for range_str in domain_ranges:
                        all_segments.extend(self.parse_range_segments(range_str))
                    
                    # Sort segments by start position
                    all_segments.sort()
                    
                    # Find gaps between segments
                    for i in range(len(all_segments) - 1):
                        current_end = all_segments[i][1]
                        next_start = all_segments[i + 1][0]
                        
                        if next_start > current_end + 1:
                            gap_start = current_end + 1
                            gap_end = next_start - 1
                            gaps.append((gap_start, gap_end))
                    
                    # Check for N-terminal gap
                    if all_segments and all_segments[0][0] > 1:
                        gaps.insert(0, (1, all_segments[0][0] - 1))
                    
                    # Check for C-terminal gap
                    if all_segments and all_segments[-1][1] < sequence_length:
                        gaps.append((all_segments[-1][1] + 1, sequence_length))
                
                largest_gap = max((end - start + 1 for start, end in gaps), default=0)
            else:
                sequence_length = 0
                total_coverage = 0
                coverage_percentage = 0
                gaps = []
                largest_gap = 0
            
            # NEW: Analyze oversegmentation
            undersized_domains, proximal_domains = self.analyze_oversegmentation(domain_ranges, domain_lengths)
            
            # NEW: Analyze structural resolution
            structural_resolution = None
            domains_in_unresolved = []
            
            if include_structure:
                structure_path = self.structural_analyzer.find_structure_file(pdb_id)
                if structure_path:
                    structural_resolution = self.structural_analyzer.parse_resolved_residues(structure_path, chain_id)
                    if structural_resolution:
                        domains_in_unresolved = self.analyze_unresolved_regions(domain_ranges, structural_resolution)
            
            return DomainAnalysis(
                protein_id=protein_id,
                pdb_id=pdb_id,
                chain_id=chain_id,
                total_domains=len(domain_elements),
                domain_ranges=domain_ranges,
                domain_lengths=domain_lengths,
                domain_positions=domain_positions,
                sequence_length=sequence_length,
                total_coverage=total_coverage,
                coverage_percentage=coverage_percentage,
                gaps=gaps,
                largest_gap=largest_gap,
                domain_sources=domain_sources,
                evidence_types=list(set(evidence_types)),
                undersized_domains=undersized_domains,
                proximal_domains=proximal_domains,
                structural_resolution=structural_resolution,
                domains_in_unresolved=domains_in_unresolved
            )
            
        except Exception as e:
            print(f"❌ Failed to analyze {protein_id}: {e}")
            return None
    
    def analyze_protein_list(self, protein_ids: List[str], include_structure: bool = True) -> List[DomainAnalysis]:
        """Analyze a list of proteins"""
        
        analyses = []
        for protein_id in protein_ids:
            analysis = self.analyze_protein(protein_id, include_structure)
            if analysis:
                analyses.append(analysis)
        
        return analyses
    
    def print_detailed_analysis(self, analyses: List[DomainAnalysis]):
        """Print detailed analysis focusing on oversegmentation and unresolved regions"""
        
        print(f"\n🔍 OVERSEGMENTATION & UNRESOLVED REGION ANALYSIS")
        print("=" * 120)
        print(f"{'Protein':<12} {'Type':<6} {'Doms':<5} {'Overseg':<8} {'Unresolved':<10} {'SmallDoms':<10} {'Resolution':<11}")
        print("-" * 120)
        
        for analysis in sorted(analyses, key=lambda x: x.oversegmentation_score + x.unresolved_region_score, reverse=True):
            protein_type = "Heavy" if analysis.protein_id.endswith('_H') else "Light" if analysis.protein_id.endswith('_L') else "Other"
            
            overseg_score = analysis.oversegmentation_score
            unresolved_score = analysis.unresolved_region_score
            small_domain_count = len(analysis.undersized_domains)
            
            resolution_info = "N/A"
            if analysis.structural_resolution:
                resolution_info = f"{analysis.structural_resolution.resolution_percentage:.1f}%"
            
            print(f"{analysis.protein_id:<12} {protein_type:<6} {analysis.total_domains:<5} "
                  f"{overseg_score:>6.1f}   {unresolved_score:>8.1f}   {small_domain_count:<10} {resolution_info:<11}")
        
        print(f"\n📊 DETAILED BREAKDOWN:")
        print("-" * 100)
        
        for analysis in analyses:
            print(f"\n🔬 {analysis.protein_id} ({analysis.pdb_id} chain {analysis.chain_id}):")
            print(f"   Sequence length: {analysis.sequence_length}")
            print(f"   Domains found: {analysis.total_domains}")
            print(f"   Coverage: {analysis.total_coverage}/{analysis.sequence_length} residues ({analysis.coverage_percentage:.1f}%)")
            
            # Structural resolution analysis
            if analysis.structural_resolution:
                sr = analysis.structural_resolution
                print(f"   📡 Structural Resolution: {sr.resolution_percentage:.1f}% ({len(sr.resolved_residues)}/{sr.total_sequence_length} residues)")
                if sr.missing_residues:
                    missing_ranges = self._format_residue_ranges(sr.missing_residues)
                    print(f"       Missing/unresolved: {missing_ranges}")
            else:
                print(f"   📡 Structural Resolution: Could not analyze (structure file not found)")
            
            # Oversegmentation analysis
            print(f"   🧩 Oversegmentation Analysis:")
            if analysis.undersized_domains:
                print(f"       Small domains (likely fragments): {len(analysis.undersized_domains)}")
                for domain_idx in analysis.undersized_domains:
                    length = analysis.domain_lengths[domain_idx]
                    range_str = analysis.domain_ranges[domain_idx]
                    print(f"         Domain {domain_idx+1}: {range_str} ({length} residues - expect ~110)")
            
            if analysis.proximal_domains:
                print(f"       Proximal domain pairs (likely oversegmented):")
                for d1_idx, d2_idx, gap_size in analysis.proximal_domains:
                    combined_size = analysis.domain_lengths[d1_idx] + analysis.domain_lengths[d2_idx] + gap_size
                    print(f"         Domains {d1_idx+1} & {d2_idx+1}: gap={gap_size}, combined_size={combined_size}")
                    print(f"           {analysis.domain_ranges[d1_idx]} + {analysis.domain_ranges[d2_idx]}")
            
            if not analysis.undersized_domains and not analysis.proximal_domains:
                print(f"       No obvious oversegmentation detected")
            
            # Unresolved region analysis
            print(f"   👻 Unresolved Region Analysis:")
            if analysis.domains_in_unresolved:
                print(f"       Domains with unresolved residues: {len(analysis.domains_in_unresolved)}")
                for domain_idx, unresolved_positions in analysis.domains_in_unresolved:
                    unresolved_count = len(unresolved_positions)
                    domain_size = analysis.domain_lengths[domain_idx]
                    unresolved_fraction = unresolved_count / domain_size
                    range_str = analysis.domain_ranges[domain_idx]
                    
                    print(f"         Domain {domain_idx+1}: {range_str}")
                    print(f"           {unresolved_count}/{domain_size} residues unresolved ({unresolved_fraction:.1%})")
                    
                    if unresolved_fraction > 0.3:
                        print(f"           ⚠️  >30% unresolved - likely spurious domain assignment!")
            else:
                if analysis.structural_resolution:
                    print(f"       All domains appear to be in resolved regions")
                else:
                    print(f"       Could not check (no structural data)")
            
            print(f"   📈 Anomaly Scores:")
            print(f"       Oversegmentation: {analysis.oversegmentation_score:.2f}")
            print(f"       Unresolved regions: {analysis.unresolved_region_score:.2f}")
    
    def _format_residue_ranges(self, residues: Set[int]) -> str:
        """Format a set of residues into range strings"""
        if not residues:
            return "None"
        
        sorted_residues = sorted(residues)
        ranges = []
        start = sorted_residues[0]
        end = start
        
        for residue in sorted_residues[1:]:
            if residue == end + 1:
                end = residue
            else:
                if start == end:
                    ranges.append(str(start))
                else:
                    ranges.append(f"{start}-{end}")
                start = end = residue
        
        # Add the last range
        if start == end:
            ranges.append(str(start))
        else:
            ranges.append(f"{start}-{end}")
        
        return ", ".join(ranges)
    
    def suggest_specific_hypotheses(self, analyses: List[DomainAnalysis]):
        """Suggest specific hypotheses based on analysis results"""
        
        print(f"\n🎯 SPECIFIC ANOMALY HYPOTHESES")
        print("=" * 80)
        
        # Classify problems
        oversegmented_proteins = [a for a in analyses if a.oversegmentation_score > 1.0]
        unresolved_domain_proteins = [a for a in analyses if a.unresolved_region_score > 1.0]
        small_domain_proteins = [a for a in analyses if len(a.undersized_domains) > 0]
        proximal_domain_proteins = [a for a in analyses if len(a.proximal_domains) > 0]
        
        print(f"🧩 HYPOTHESIS 1: OVERSEGMENTATION ({len(oversegmented_proteins)} proteins)")
        print(f"   Single Ig domains are being split into multiple pieces")
        print(f"   Evidence:")
        print(f"     - {len(small_domain_proteins)} proteins with undersized domains (<90 residues)")
        print(f"     - {len(proximal_domain_proteins)} proteins with domains very close together (<30 residues)")
        
        if proximal_domain_proteins:
            print(f"   Likely oversegmented candidates:")
            for analysis in proximal_domain_proteins[:3]:  # Show top 3
                for d1_idx, d2_idx, gap_size in analysis.proximal_domains:
                    combined_size = analysis.domain_lengths[d1_idx] + analysis.domain_lengths[d2_idx] + gap_size
                    if 90 <= combined_size <= 150:  # Would make a normal Ig domain
                        print(f"     {analysis.protein_id}: domains {d1_idx+1}&{d2_idx+1} → {combined_size} residues (normal Ig size)")
        
        print(f"\n👻 HYPOTHESIS 2: DOMAINS IN UNRESOLVED REGIONS ({len(unresolved_domain_proteins)} proteins)")
        print(f"   Domains assigned to crystal structure regions with missing density")
        print(f"   Evidence:")
        
        if unresolved_domain_proteins:
            for analysis in unresolved_domain_proteins:
                if analysis.structural_resolution:
                    print(f"     {analysis.protein_id}: {analysis.structural_resolution.resolution_percentage:.1f}% resolved")
                    for domain_idx, unresolved_positions in analysis.domains_in_unresolved:
                        unresolved_count = len(unresolved_positions)
                        domain_size = analysis.domain_lengths[domain_idx]
                        fraction = unresolved_count / domain_size
                        if fraction > 0.3:
                            print(f"       Domain {domain_idx+1}: {fraction:.1%} unresolved (likely spurious)")
        
        print(f"\n🎯 RECOMMENDED FIXES:")
        print(f"   1. For oversegmentation: Increase minimum domain size threshold")
        print(f"   2. For oversegmentation: Merge domains separated by <20 residues")
        print(f"   3. For unresolved regions: Filter domains with >30% unresolved residues")
        print(f"   4. For unresolved regions: Check ECOD references for structural completeness")


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Analyze Anomalous Double Ig Domain Predictions - Focus on Oversegmentation & Unresolved Regions'
    )
    
    parser.add_argument('--analyze-proteins', type=str,
                       help='Comma-separated list of protein IDs to analyze')
    parser.add_argument('--oversegmentation-analysis', action='store_true',
                       help='Focus on oversegmentation patterns')
    parser.add_argument('--structure-resolution-check', action='store_true',
                       help='Check domains against structural resolution')
    parser.add_argument('--no-structure', action='store_true',
                       help='Skip structural analysis (faster)')
    
    args = parser.parse_args()
    
    # Default problematic proteins if none specified
    if not args.analyze_proteins:
        problematic_proteins = [
            "8g8c_H", "8g8d_H", "8gav_H", "8gau_H", "8g8a_A", 
            "8gav_L", "8g5b_H", "8g2m_L", "8g5b_I", "8g5b_M"
        ]
    else:
        problematic_proteins = [p.strip() for p in args.analyze_proteins.split(',')]
    
    print(f"🔍 Analyzing {len(problematic_proteins)} potentially anomalous proteins...")
    print(f"Focus: Oversegmentation & Unresolved Regions")
    
    analyzer = AnomalyAnalyzer()
    include_structure = not args.no_structure
    
    if not include_structure:
        print("⚠️  Skipping structural analysis (--no-structure)")
    
    analyses = analyzer.analyze_protein_list(problematic_proteins, include_structure)
    
    if not analyses:
        print("❌ No proteins could be analyzed")
        return
    
    print(f"✓ Successfully analyzed {len(analyses)} proteins")
    
    # Detailed analysis with focus on new hypotheses
    analyzer.print_detailed_analysis(analyses)
    
    # Specific hypotheses for these types of anomalies
    analyzer.suggest_specific_hypotheses(analyses)
    
    # Summary statistics
    print(f"\n📊 SUMMARY STATISTICS")
    print("=" * 60)
    
    oversegmented_count = len([a for a in analyses if a.oversegmentation_score > 1.0])
    unresolved_count = len([a for a in analyses if a.unresolved_region_score > 1.0])
    small_domain_count = len([a for a in analyses if len(a.undersized_domains) > 0])
    proximal_count = len([a for a in analyses if len(a.proximal_domains) > 0])
    
    print(f"Proteins with oversegmentation evidence: {oversegmented_count}/{len(analyses)} ({oversegmented_count/len(analyses)*100:.1f}%)")
    print(f"Proteins with domains in unresolved regions: {unresolved_count}/{len(analyses)} ({unresolved_count/len(analyses)*100:.1f}%)")
    print(f"Proteins with undersized domains: {small_domain_count}/{len(analyses)} ({small_domain_count/len(analyses)*100:.1f}%)")
    print(f"Proteins with proximal domains: {proximal_count}/{len(analyses)} ({proximal_count/len(analyses)*100:.1f}%)")
    
    if include_structure:
        structural_data_count = len([a for a in analyses if a.structural_resolution is not None])
        print(f"Proteins with structural data: {structural_data_count}/{len(analyses)} ({structural_data_count/len(analyses)*100:.1f}%)")
        
        if structural_data_count > 0:
            avg_resolution = sum(a.structural_resolution.resolution_percentage 
                               for a in analyses if a.structural_resolution) / structural_data_count
            print(f"Average structural resolution: {avg_resolution:.1f}%")


if __name__ == "__main__":
    main()
