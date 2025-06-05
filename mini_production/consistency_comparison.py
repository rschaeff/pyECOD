#!/usr/bin/env python3
"""
Consistency-Based Domain Comparison

Evaluates domain partitions based on consistency with established ECOD patterns
rather than absolute truth. The database becomes the reference standard.

Key principle: Good annotations are consistent with existing high-quality annotations
for the same families.
"""

import sys
import yaml
import psycopg2
import psycopg2.extras
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict, Counter
import numpy as np
from dataclasses import dataclass
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class DomainAnnotation:
    """Represents a domain annotation for consistency analysis"""
    protein_id: str
    domain_range: str
    t_group: str
    h_group: str
    x_group: str
    length: int
    start_pos: int
    end_pos: int
    confidence: float = 0.0

@dataclass
class FamilyPattern:
    """Represents consistent patterns within a domain family"""
    family_id: str  # t_group or h_group
    typical_length_range: Tuple[int, int]
    common_boundaries: List[int]  # Common start/end positions (normalized)
    split_frequency: float  # How often this family is split into multiple domains
    representative_examples: List[str]  # High-quality example protein IDs

class ConsistencyComparator:
    """Compare domain partitions based on database consistency"""
    
    def __init__(self, config_path: str = "config/config.local.yml"):
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)
        
        self.db_conn = psycopg2.connect(**self.config['database'])
        
        # Cache for family patterns to avoid repeated queries
        self.family_patterns = {}
        
    def get_family_patterns(self, family_ids: List[str], family_type: str = "t_group") -> Dict[str, FamilyPattern]:
        """Extract consistency patterns for domain families from ECOD database"""
        
        patterns = {}
        
        for family_id in family_ids:
            if family_id in self.family_patterns:
                patterns[family_id] = self.family_patterns[family_id]
                continue
                
            try:
                with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
                    # Get all domains in this family from high-quality sources
                    cursor.execute("""
                        SELECT 
                            p.source_id,
                            d.range,
                            d.t_group,
                            d.h_group, 
                            d.x_group,
                            d.length,
                            ps.resolution
                        FROM pdb_analysis.domain d
                        JOIN pdb_analysis.protein p ON d.protein_id = p.id
                        LEFT JOIN pdb_analysis.protein_structure ps ON p.id = ps.protein_id
                        WHERE d.{} = %s
                        AND d.length BETWEEN 30 AND 1000  -- Reasonable domain sizes
                        AND (ps.resolution IS NULL OR ps.resolution < 3.0)  -- High-quality structures
                        ORDER BY ps.resolution ASC NULLS LAST
                        LIMIT 200  -- Don't overload with too many examples
                    """.format(family_type), (family_id,))
                    
                    rows = cursor.fetchall()
                    
                    if not rows:
                        continue
                    
                    # Analyze patterns
                    lengths = [row['length'] for row in rows if row['length']]
                    ranges = [row['range'] for row in rows if row['range']]
                    
                    # Length distribution
                    if lengths:
                        length_mean = np.mean(lengths)
                        length_std = np.std(lengths)
                        typical_range = (
                            max(30, int(length_mean - 2 * length_std)),
                            int(length_mean + 2 * length_std)
                        )
                    else:
                        typical_range = (50, 300)  # Default range
                    
                    # Boundary analysis (simplified - would need more sophisticated analysis)
                    boundaries = []
                    for range_str in ranges:
                        try:
                            if '-' in range_str:
                                parts = range_str.split(',')[0]  # Take first segment for now
                                if '-' in parts:
                                    start, end = map(int, parts.split('-'))
                                    boundaries.extend([start % 10, end % 10])  # Modulo for pattern analysis
                        except (ValueError, IndexError):
                            continue
                    
                    common_boundaries = []
                    if boundaries:
                        boundary_counts = Counter(boundaries)
                        # Find boundaries that appear in >20% of cases
                        threshold = len(ranges) * 0.2
                        common_boundaries = [b for b, count in boundary_counts.items() if count > threshold]
                    
                    # How often is this family split across multiple domains in one protein?
                    cursor.execute("""
                        SELECT p.source_id, COUNT(d.id) as domain_count
                        FROM pdb_analysis.domain d
                        JOIN pdb_analysis.protein p ON d.protein_id = p.id
                        WHERE d.{} = %s
                        GROUP BY p.source_id
                        HAVING COUNT(d.id) > 1
                    """.format(family_type), (family_id,))
                    
                    split_cases = cursor.fetchall()
                    
                    # FIXED: Calculate split frequency properly
                    total_proteins_with_family = len(set(row['source_id'] for row in rows))
                    proteins_with_splits = len(split_cases)
                    
                    split_frequency = proteins_with_splits / max(1, total_proteins_with_family)
                    split_frequency = max(0.0, min(1.0, split_frequency))  # Ensure 0-1 range
                    
                    # Representative examples (highest quality)
                    representatives = [row['source_id'] for row in rows[:10]]
                    
                    pattern = FamilyPattern(
                        family_id=family_id,
                        typical_length_range=typical_range,
                        common_boundaries=common_boundaries,
                        split_frequency=split_frequency,
                        representative_examples=representatives
                    )
                    
                    patterns[family_id] = pattern
                    self.family_patterns[family_id] = pattern
                    
            except Exception as e:
                logger.warning(f"Could not analyze pattern for family {family_id}: {e}")
                continue
        
        return patterns
    
    def parse_partition_domains(self, xml_file: Path) -> List[DomainAnnotation]:
        """Parse domain partition into DomainAnnotation objects"""
        
        domains = []
        
        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()
            
            protein_id = xml_file.stem.replace('.mini.domains', '').replace('.develop291.domains', '')
            
            for domain_elem in root.findall(".//domain"):
                range_str = domain_elem.get('range', '')
                t_group = domain_elem.get('t_group', '') or domain_elem.get('family', '')
                h_group = domain_elem.get('h_group', '')
                x_group = domain_elem.get('x_group', '')
                
                # FIXED: Better length parsing
                length = int(domain_elem.get('length', '0') or domain_elem.get('size', '0'))
                confidence = float(domain_elem.get('confidence', '0.0'))
                
                # Parse start/end positions and calculate length if missing
                start_pos, end_pos = 0, 0
                calculated_length = 0
                
                if range_str and '-' in range_str:
                    try:
                        # Parse all segments to get total length
                        segments = []
                        total_length = 0
                        
                        for segment in range_str.split(','):
                            segment = segment.strip()
                            if '-' in segment:
                                s, e = map(int, segment.split('-'))
                                segments.extend([s, e])
                                total_length += (e - s + 1)
                        
                        if segments:
                            start_pos = min(segments)
                            end_pos = max(segments)
                            calculated_length = total_length
                    except (ValueError, IndexError):
                        pass
                
                # Use calculated length if XML length is missing
                if length == 0 and calculated_length > 0:
                    length = calculated_length
                
                if t_group:  # Only include classified domains
                    domains.append(DomainAnnotation(
                        protein_id=protein_id,
                        domain_range=range_str,
                        t_group=t_group,
                        h_group=h_group,
                        x_group=x_group,
                        length=length,
                        start_pos=start_pos,
                        end_pos=end_pos,
                        confidence=confidence
                    ))
        
        except Exception as e:
            logger.error(f"Error parsing {xml_file}: {e}")
        
        return domains
    
    def calculate_family_consistency_score(self, domains: List[DomainAnnotation]) -> float:
        """Calculate how consistent domains are with their family patterns"""
        
        if not domains:
            return 0.0
        
        # Get family patterns for all families present
        family_ids = list(set(d.t_group for d in domains if d.t_group))
        if not family_ids:
            return 0.0
        
        family_patterns = self.get_family_patterns(family_ids, "t_group")
        
        family_scores = []
        
        for domain in domains:
            if domain.t_group not in family_patterns:
                family_scores.append(0.5)  # Unknown family - neutral score
                continue
            
            pattern = family_patterns[domain.t_group]
            
            # FIXED: Ensure individual scores are properly bounded
            score_components = []
            
            # 1. Length consistency
            if domain.length == 0:
                length_score = 0.0  # Can't evaluate if no length
            elif pattern.typical_length_range[0] <= domain.length <= pattern.typical_length_range[1]:
                length_score = 1.0
            else:
                # Penalize based on how far outside typical range
                if domain.length < pattern.typical_length_range[0]:
                    ratio = domain.length / max(1, pattern.typical_length_range[0])
                else:
                    ratio = pattern.typical_length_range[1] / max(1, domain.length)
                length_score = max(0.0, min(1.0, ratio))
            
            score_components.append(length_score * 0.4)
            
            # 2. Split consistency
            family_count = sum(1 for d in domains if d.t_group == domain.t_group)
            
            if family_count == 1:
                # Single domain of this family - good if family usually appears once
                split_score = max(0.0, 1.0 - pattern.split_frequency)
            else:
                # Multiple domains of this family - good if family usually gets split
                split_score = min(1.0, pattern.split_frequency)
            
            score_components.append(split_score * 0.3)
            
            # 3. Confidence component
            confidence_score = max(0.0, min(1.0, domain.confidence))
            score_components.append(confidence_score * 0.3)
            
            # FIXED: Sum components properly
            total_score = sum(score_components)
            family_scores.append(max(0.0, min(1.0, total_score)))
        
        return max(0.0, min(1.0, np.mean(family_scores)))
    
    def calculate_boundary_consistency_score(self, domains: List[DomainAnnotation]) -> float:
        """Calculate consistency of domain boundaries"""
        
        if len(domains) < 2:
            return 1.0  # Single domain - no boundary issues
        
        # Sort domains by position
        sorted_domains = sorted(domains, key=lambda d: d.start_pos)
        
        boundary_scores = []
        
        for i in range(len(sorted_domains) - 1):
            current = sorted_domains[i]
            next_domain = sorted_domains[i + 1]
            
            # Check for gaps or overlaps
            gap = next_domain.start_pos - current.end_pos - 1
            
            if gap == 0:
                # Perfect boundary - domains are adjacent
                boundary_scores.append(1.0)
            elif 1 <= gap <= 10:
                # Small gap - reasonable linker region
                boundary_scores.append(0.8)
            elif -10 <= gap < 0:
                # Small overlap - might be annotation uncertainty
                boundary_scores.append(0.6)
            elif gap > 10:
                # Large gap - unusual but not necessarily wrong
                boundary_scores.append(0.4)
            else:
                # Large overlap - likely problematic
                boundary_scores.append(0.2)
        
        return np.mean(boundary_scores) if boundary_scores else 1.0
    
    def calculate_overall_consistency_score(self, domains: List[DomainAnnotation]) -> Dict[str, float]:
        """Calculate overall consistency score with breakdown"""
        
        family_score = self.calculate_family_consistency_score(domains)
        boundary_score = self.calculate_boundary_consistency_score(domains)
        
        # Coverage score - do domains cover reasonable portion of protein?
        if domains:
            total_coverage = sum(d.length for d in domains)
            # Estimate protein length from domain span
            positions = [d.start_pos for d in domains if d.start_pos > 0] + [d.end_pos for d in domains if d.end_pos > 0]
            if positions:
                min_pos = min(positions)
                max_pos = max(positions)
                estimated_length = max(max_pos - min_pos + 1, total_coverage)
            else:
                estimated_length = total_coverage
            
            if estimated_length > 0:
                coverage_ratio = total_coverage / estimated_length
                coverage_score = min(1.0, coverage_ratio / 0.8)  # 80% coverage = perfect
            else:
                coverage_score = 0.0
        else:
            coverage_score = 0.0
        
        # FIXED: Ensure all scores are bounded 0-1
        family_score = max(0.0, min(1.0, family_score))
        boundary_score = max(0.0, min(1.0, boundary_score))
        coverage_score = max(0.0, min(1.0, coverage_score))
        
        # Weighted combination
        overall_score = (
            family_score * 0.5 +      # Family consistency most important
            boundary_score * 0.3 +    # Boundary quality important
            coverage_score * 0.2      # Coverage less important
        )
        
        # FIXED: Ensure overall score is bounded
        overall_score = max(0.0, min(1.0, overall_score))
        
        return {
            'overall': overall_score,
            'family_consistency': family_score,
            'boundary_quality': boundary_score,
            'coverage': coverage_score,
            'domain_count': len(domains)
        }
    
    def compare_partitions_consistency(self, mini_file: Path, main_file: Path) -> Dict:
        """Compare two partitions based on consistency scores"""
        
        # Parse both partitions
        mini_domains = self.parse_partition_domains(mini_file)
        main_domains = self.parse_partition_domains(main_file) if main_file.exists() else []
        
        # Calculate consistency scores
        mini_scores = self.calculate_overall_consistency_score(mini_domains)
        main_scores = self.calculate_overall_consistency_score(main_domains) if main_domains else {
            'overall': 0.0, 'family_consistency': 0.0, 'boundary_quality': 0.0, 
            'coverage': 0.0, 'domain_count': 0
        }
        
        # Determine winner
        score_diff = mini_scores['overall'] - main_scores['overall']
        
        if abs(score_diff) < 0.05:
            winner = "tie"
            reason = "similar_consistency"
        elif score_diff > 0:
            winner = "mini"
            # Determine why mini wins
            if mini_scores['family_consistency'] > main_scores['family_consistency'] + 0.1:
                reason = "better_family_consistency"
            elif mini_scores['boundary_quality'] > main_scores['boundary_quality'] + 0.1:
                reason = "better_boundaries"
            else:
                reason = "overall_consistency"
        else:
            winner = "main"
            if main_scores['family_consistency'] > mini_scores['family_consistency'] + 0.1:
                reason = "better_family_consistency"
            elif main_scores['boundary_quality'] > mini_scores['boundary_quality'] + 0.1:
                reason = "better_boundaries"
            else:
                reason = "overall_consistency"
        
        return {
            'protein_id': mini_domains[0].protein_id if mini_domains else "unknown",
            'winner': winner,
            'reason': reason,
            'score_difference': score_diff,
            'mini_scores': mini_scores,
            'main_scores': main_scores,
            'mini_domains': len(mini_domains),
            'main_domains': len(main_domains),
            'mini_families': list(set(d.t_group for d in mini_domains if d.t_group)),
            'main_families': list(set(d.t_group for d in main_domains if d.t_group))
        }

def debug_consistency_comparison(protein_id: str = "8ovp_A"):
    """Debug consistency-based comparison for specific protein"""
    
    print(f"🔍 CONSISTENCY-BASED Analysis for {protein_id}")
    print("=" * 70)
    
    comparator = ConsistencyComparator()
    
    # Find files
    batch_base = Path("../data/ecod/pdb_updates/batches")
    mini_file = None
    main_file = None
    
    for batch_dir in batch_base.iterdir():
        if not batch_dir.is_dir():
            continue
        
        # Look for mini file
        mini_domains_dir = batch_dir / "mini_domains"
        if mini_domains_dir.exists():
            test_file = mini_domains_dir / f"{protein_id}.mini.domains.xml"
            if test_file.exists():
                mini_file = test_file
        
        # Look for main file
        domains_dir = batch_dir / "domains"
        if domains_dir.exists():
            test_file = domains_dir / f"{protein_id}.develop291.domains.xml"
            if test_file.exists():
                main_file = test_file
    
    if not mini_file:
        print("❌ No mini file found!")
        return
    
    # Parse domains
    mini_domains = comparator.parse_partition_domains(mini_file)
    main_domains = comparator.parse_partition_domains(main_file) if main_file and main_file.exists() else []
    
    print(f"📊 DOMAIN ANNOTATIONS:")
    print(f"  Mini: {len(mini_domains)} domains")
    for i, domain in enumerate(mini_domains):
        print(f"    D{i+1}: {domain.t_group} | {domain.domain_range} | length={domain.length}")
    
    print(f"  Main: {len(main_domains)} domains")
    for i, domain in enumerate(main_domains):
        print(f"    D{i+1}: {domain.t_group} | {domain.domain_range} | length={domain.length}")
    
    # Get family patterns
    all_families = list(set(d.t_group for d in mini_domains + main_domains if d.t_group))
    family_patterns = comparator.get_family_patterns(all_families)
    
    print(f"\n📚 FAMILY PATTERNS:")
    for family_id, pattern in family_patterns.items():
        print(f"  {family_id}:")
        print(f"    Typical length: {pattern.typical_length_range[0]}-{pattern.typical_length_range[1]}")
        print(f"    Split frequency: {pattern.split_frequency:.2f}")
        print(f"    Examples: {', '.join(pattern.representative_examples[:3])}")
    
    # Calculate consistency scores
    mini_scores = comparator.calculate_overall_consistency_score(mini_domains)
    main_scores = comparator.calculate_overall_consistency_score(main_domains)
    
    print(f"\n🎯 CONSISTENCY SCORES:")
    print(f"  Mini Algorithm:")
    print(f"    Overall: {mini_scores['overall']:.3f}")
    print(f"    Family consistency: {mini_scores['family_consistency']:.3f}")
    print(f"    Boundary quality: {mini_scores['boundary_quality']:.3f}")
    print(f"    Coverage: {mini_scores['coverage']:.3f}")
    
    print(f"  Main Algorithm:")
    print(f"    Overall: {main_scores['overall']:.3f}")
    print(f"    Family consistency: {main_scores['family_consistency']:.3f}")
    print(f"    Boundary quality: {main_scores['boundary_quality']:.3f}")
    print(f"    Coverage: {main_scores['coverage']:.3f}")
    
    # Winner
    score_diff = mini_scores['overall'] - main_scores['overall']
    if score_diff > 0.05:
        winner = "MINI"
    elif score_diff < -0.05:
        winner = "MAIN"
    else:
        winner = "TIE"
    
    print(f"\n🏆 CONSISTENCY WINNER: {winner}")
    print(f"  Score difference: {score_diff:+.3f}")
    
    # Analysis
    print(f"\n💡 CONSISTENCY ANALYSIS:")
    
    for domain in mini_domains:
        if domain.t_group in family_patterns:
            pattern = family_patterns[domain.t_group]
            in_range = pattern.typical_length_range[0] <= domain.length <= pattern.typical_length_range[1]
            print(f"  Mini {domain.t_group}: length {domain.length} "
                  f"{'✅' if in_range else '⚠️'} (typical: {pattern.typical_length_range[0]}-{pattern.typical_length_range[1]})")
    
    for domain in main_domains:
        if domain.t_group in family_patterns:
            pattern = family_patterns[domain.t_group]
            in_range = pattern.typical_length_range[0] <= domain.length <= pattern.typical_length_range[1]
            print(f"  Main {domain.t_group}: length {domain.length} "
                  f"{'✅' if in_range else '⚠️'} (typical: {pattern.typical_length_range[0]}-{pattern.typical_length_range[1]})")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Consistency-based domain comparison')
    parser.add_argument('--protein', default='8ovp_A', help='Protein to analyze')
    
    args = parser.parse_args()
    
    debug_consistency_comparison(args.protein)
