#!/usr/bin/env python3
"""
Overnight Consistency Analysis - Production Ready

Scans all available mini results and identifies:
1. Family preservation patterns  
2. Spurious tiny domains
3. Over-segmentation cases
4. Best/worst examples for slide deck

Run: python overnight_consistency_scan.py --scan-all --output-report
"""

import os
import sys
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Tuple
from collections import defaultdict, Counter
from dataclasses import dataclass, asdict
import argparse

@dataclass
class DomainIssue:
    protein_id: str
    issue_type: str  # 'tiny_domain', 'family_fragmentation', 'over_segmentation'
    details: str
    severity: str    # 'high', 'medium', 'low'

@dataclass
class FamilyAnalysis:
    family_id: str
    total_occurrences: int
    fragmentation_cases: int
    typical_length_range: Tuple[int, int]
    suspicious_proteins: List[str]

class ConsistencyScanner:
    """Fast consistency scanner for batch analysis"""
    
    def __init__(self, batch_base_dir: str = "../data/ecod/pdb_updates/batches"):
        self.batch_base = Path(batch_base_dir)
        self.issues = []
        self.family_stats = defaultdict(lambda: {
            'occurrences': [],
            'lengths': [],
            'proteins': [],
            'fragmentations': []
        })
        
    def scan_all_mini_results(self) -> Dict:
        """Scan all available mini results for consistency issues"""
        
        print("🔍 Scanning all mini results for consistency issues...")
        
        total_proteins = 0
        total_domains = 0
        
        # Find all mini result files
        for batch_dir in self.batch_base.iterdir():
            if not batch_dir.is_dir():
                continue
                
            mini_domains_dir = batch_dir / "mini_domains"
            if not mini_domains_dir.exists():
                continue
                
            xml_files = list(mini_domains_dir.glob("*.mini.domains.xml"))
            print(f"  Processing {batch_dir.name}: {len(xml_files)} files")
            
            for xml_file in xml_files:
                protein_data = self.analyze_protein_file(xml_file)
                if protein_data:
                    total_proteins += 1
                    total_domains += protein_data['domain_count']
        
        print(f"✅ Processed {total_proteins} proteins, {total_domains} domains")
        
        # Analyze family patterns
        family_analysis = self.analyze_family_patterns()
        
        # Generate summary
        summary = {
            'scan_summary': {
                'total_proteins': total_proteins,
                'total_domains': total_domains,
                'total_issues': len(self.issues),
                'families_analyzed': len(family_analysis)
            },
            'issue_breakdown': self.categorize_issues(),
            'family_analysis': family_analysis,
            'worst_examples': self.get_worst_examples(),
            'best_examples': self.get_best_examples()
        }
        
        return summary
    
    def analyze_protein_file(self, xml_file: Path) -> Dict:
        """Analyze single protein for consistency issues"""
        
        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()
            
            protein_id = xml_file.stem.replace('.mini.domains', '')
            domains = []
            
            for domain_elem in root.findall(".//domain"):
                range_str = domain_elem.get('range', '')
                family = domain_elem.get('t_group', '') or domain_elem.get('family', '')
                size = int(domain_elem.get('size', '0') or domain_elem.get('length', '0'))
                
                # Calculate size from range if missing
                if size == 0 and range_str:
                    try:
                        total_length = 0
                        for segment in range_str.split(','):
                            if '-' in segment:
                                start, end = map(int, segment.split('-'))
                                total_length += (end - start + 1)
                        size = total_length
                    except:
                        pass
                
                if family:  # Only analyze classified domains
                    domains.append({
                        'family': family,
                        'range': range_str,
                        'size': size
                    })
                    
                    # Track family statistics
                    self.family_stats[family]['occurrences'].append(protein_id)
                    self.family_stats[family]['lengths'].append(size)
                    self.family_stats[family]['proteins'].append(protein_id)
            
            # Check for issues
            self.check_protein_issues(protein_id, domains)
            
            return {
                'protein_id': protein_id,
                'domain_count': len(domains),
                'families': [d['family'] for d in domains]
            }
            
        except Exception as e:
            print(f"Error processing {xml_file}: {e}")
            return None
    
    def check_protein_issues(self, protein_id: str, domains: List[Dict]):
        """Check for consistency issues in a protein"""
        
        # 1. Tiny domains (likely spurious)
        for domain in domains:
            if domain['size'] < 30:
                self.issues.append(DomainIssue(
                    protein_id=protein_id,
                    issue_type='tiny_domain',
                    details=f"Domain {domain['family']} only {domain['size']} residues",
                    severity='high'
                ))
        
        # 2. Family fragmentation (same family multiple times)
        family_counts = Counter(d['family'] for d in domains)
        for family, count in family_counts.items():
            if count > 1:
                self.family_stats[family]['fragmentations'].append(protein_id)
                
                # High fragmentation is suspicious
                if count > 3:
                    self.issues.append(DomainIssue(
                        protein_id=protein_id,
                        issue_type='family_fragmentation',
                        details=f"Family {family} split into {count} pieces",
                        severity='high'
                    ))
                elif count == 2:
                    self.issues.append(DomainIssue(
                        protein_id=protein_id,
                        issue_type='family_fragmentation',
                        details=f"Family {family} split into {count} pieces",
                        severity='medium'
                    ))
        
        # 3. Over-segmentation (too many domains for protein size)
        total_coverage = sum(d['size'] for d in domains)
        if len(domains) > 5 and total_coverage < 300:
            self.issues.append(DomainIssue(
                protein_id=protein_id,
                issue_type='over_segmentation',
                details=f"{len(domains)} domains in {total_coverage} residues",
                severity='medium'
            ))
    
    def analyze_family_patterns(self) -> Dict[str, FamilyAnalysis]:
        """Analyze patterns within each family"""
        
        family_analysis = {}
        
        for family_id, stats in self.family_stats.items():
            lengths = stats['lengths']
            if not lengths:
                continue
                
            # Calculate typical length range (remove outliers)
            sorted_lengths = sorted(lengths)
            if len(sorted_lengths) >= 3:
                # Use 10th to 90th percentile as "typical" range
                p10_idx = max(0, len(sorted_lengths) // 10)
                p90_idx = min(len(sorted_lengths) - 1, len(sorted_lengths) * 9 // 10)
                typical_range = (sorted_lengths[p10_idx], sorted_lengths[p90_idx])
            else:
                typical_range = (min(lengths), max(lengths))
            
            # Identify suspicious proteins
            suspicious = []
            
            # Proteins with unusually small domains in this family
            for i, length in enumerate(lengths):
                protein_id = stats['proteins'][i]
                if length < typical_range[0] * 0.5:  # Less than half typical size
                    suspicious.append(f"{protein_id}(size:{length})")
            
            # Proteins with high fragmentation
            fragmentation_counts = Counter(stats['fragmentations'])
            for protein_id, frag_count in fragmentation_counts.items():
                if frag_count > 2:
                    suspicious.append(f"{protein_id}(frags:{frag_count})")
            
            family_analysis[family_id] = FamilyAnalysis(
                family_id=family_id,
                total_occurrences=len(stats['occurrences']),
                fragmentation_cases=len(stats['fragmentations']),
                typical_length_range=typical_range,
                suspicious_proteins=suspicious[:10]  # Top 10 most suspicious
            )
        
        return family_analysis
    
    def categorize_issues(self) -> Dict:
        """Categorize and count issues"""
        
        breakdown = {
            'tiny_domain': {'high': 0, 'medium': 0, 'low': 0},
            'family_fragmentation': {'high': 0, 'medium': 0, 'low': 0},
            'over_segmentation': {'high': 0, 'medium': 0, 'low': 0}
        }
        
        for issue in self.issues:
            breakdown[issue.issue_type][issue.severity] += 1
        
        return breakdown
    
    def get_worst_examples(self) -> List[Dict]:
        """Get worst examples for detailed analysis"""
        
        # Group issues by protein
        protein_issues = defaultdict(list)
        for issue in self.issues:
            protein_issues[issue.protein_id].append(issue)
        
        # Sort by number of high severity issues
        worst_proteins = []
        for protein_id, issues in protein_issues.items():
            high_severity_count = sum(1 for i in issues if i.severity == 'high')
            total_issues = len(issues)
            
            worst_proteins.append({
                'protein_id': protein_id,
                'high_severity_issues': high_severity_count,
                'total_issues': total_issues,
                'issue_details': [asdict(i) for i in issues]
            })
        
        # Return top 20 worst
        return sorted(worst_proteins, key=lambda x: (x['high_severity_issues'], x['total_issues']), reverse=True)[:20]
    
    def get_best_examples(self) -> List[str]:
        """Get proteins with no consistency issues (good examples)"""
        
        proteins_with_issues = set(issue.protein_id for issue in self.issues)
        
        # Find proteins with multiple domains but no issues
        all_proteins = set()
        for stats in self.family_stats.values():
            all_proteins.update(stats['proteins'])
        
        good_proteins = []
        for protein_id in all_proteins:
            if protein_id not in proteins_with_issues:
                # Count domains for this protein
                domain_count = sum(1 for stats in self.family_stats.values() 
                                 if protein_id in stats['proteins'])
                if domain_count >= 2:  # Multi-domain proteins with no issues
                    good_proteins.append(protein_id)
        
        return good_proteins[:20]  # Top 20 good examples

def generate_slide_ready_report(summary: Dict, output_file: str):
    """Generate slide-ready summary report"""
    
    lines = [
        "# Mini Algorithm Consistency Analysis Report",
        "",
        f"**Analysis Date:** {os.popen('date').read().strip()}",
        f"**Proteins Analyzed:** {summary['scan_summary']['total_proteins']:,}",
        f"**Domains Analyzed:** {summary['scan_summary']['total_domains']:,}",
        "",
        "## 🚨 Issue Summary",
        "",
        "### High-Priority Issues",
        f"- **Tiny Domains:** {summary['issue_breakdown']['tiny_domain']['high']} cases",
        f"- **Family Fragmentation:** {summary['issue_breakdown']['family_fragmentation']['high']} cases", 
        f"- **Over-segmentation:** {summary['issue_breakdown']['over_segmentation']['high']} cases",
        "",
        "### Most Problematic Families",
        ""
    ]
    
    # Top families with issues
    family_issues = []
    for family_id, analysis in summary['family_analysis'].items():
        if analysis.fragmentation_cases > 0:
            fragmentation_rate = analysis.fragmentation_cases / analysis.total_occurrences
            family_issues.append((family_id, fragmentation_rate, analysis.fragmentation_cases))
    
    family_issues.sort(key=lambda x: x[1], reverse=True)
    
    for family_id, rate, count in family_issues[:10]:
        lines.append(f"- **{family_id}**: {count} fragmentations in {summary['family_analysis'][family_id].total_occurrences} occurrences ({rate:.1%})")
    
    lines.extend([
        "",
        "## 🎯 Worst Examples (for detailed analysis)",
        ""
    ])
    
    for example in summary['worst_examples'][:10]:
        lines.append(f"- **{example['protein_id']}**: {example['high_severity_issues']} high-severity issues")
    
    lines.extend([
        "",
        "## ✅ Best Examples (for positive showcase)", 
        ""
    ])
    
    for protein_id in summary['best_examples'][:10]:
        lines.append(f"- **{protein_id}**: Multi-domain, no consistency issues")
    
    lines.extend([
        "",
        "## 📊 Key Statistics for Slides",
        "",
        f"- **Issue Rate:** {len(summary.get('worst_examples', [])) / summary['scan_summary']['total_proteins'] * 100:.1f}% of proteins have consistency issues",
        f"- **Family Preservation:** {summary['scan_summary']['families_analyzed']} families analyzed",
        f"- **Spurious Domains:** {summary['issue_breakdown']['tiny_domain']['high']} domains < 30 residues (likely spurious)",
        ""
    ])
    
    with open(output_file, 'w') as f:
        f.write('\n'.join(lines))

def main():
    parser = argparse.ArgumentParser(description='Overnight Consistency Analysis')
    parser.add_argument('--scan-all', action='store_true', help='Scan all available mini results')
    parser.add_argument('--output-report', action='store_true', help='Generate slide-ready report')
    parser.add_argument('--batch-dir', default='../data/ecod/pdb_updates/batches', help='Batch directory')
    parser.add_argument('--output-dir', default='/tmp/consistency_analysis', help='Output directory')
    
    args = parser.parse_args()
    
    if not args.scan_all:
        parser.print_help()
        print("\n💡 Run: python overnight_consistency_scan.py --scan-all --output-report")
        return
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Run analysis
    scanner = ConsistencyScanner(args.batch_dir)
    summary = scanner.scan_all_mini_results()
    
    # Save detailed results
    with open(output_dir / 'detailed_analysis.json', 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    
    if args.output_report:
        generate_slide_ready_report(summary, str(output_dir / 'slide_ready_report.md'))
        print(f"📋 Slide-ready report: {output_dir / 'slide_ready_report.md'}")
    
    print(f"✅ Analysis complete!")
    print(f"📁 Results: {output_dir}")
    print(f"🔍 Found {len(scanner.issues)} consistency issues across {summary['scan_summary']['total_proteins']} proteins")

if __name__ == "__main__":
    main()
