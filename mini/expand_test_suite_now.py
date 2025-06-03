#!/usr/bin/env python3
"""
Quick script to expand mini_pyecod test suite with proteins from actual batch

This analyzes proteins in your batch to find diverse test cases based on
their T-group hits and evidence characteristics.

Usage:
    python expand_test_suite_now.py
    python expand_test_suite_now.py --sample-size 50
    python expand_test_suite_now.py --generate-tests
"""

import os
import sys
import argparse
from pathlib import Path
from typing import List, Dict, Tuple, Set
import xml.etree.ElementTree as ET
from collections import defaultdict

# Add parent directory
sys.path.insert(0, str(Path(__file__).parent.parent))


class BatchProteinAnalyzer:
    """Analyze proteins in batch to find test candidates"""
    
    def __init__(self, batch_dir: str = None):
        self.batch_dir = Path(batch_dir or "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424")
        self.domains_dir = self.batch_dir / "domains"
        
        # Track what we find
        self.proteins_analyzed = []
        self.categorized = defaultdict(list)
        
        # Categories based on actual batch content
        self.categories = {
            "validated": "Already validated test case",
            "single_domain_small": "Small proteins with single domain signatures",
            "single_domain_medium": "Medium proteins with single domain signatures", 
            "multi_domain_clear": "Clear multi-domain proteins",
            "chain_blast_multi": "Multi-domain via chain BLAST (good for decomposition)",
            "large_complex": "Large proteins with complex architecture",
            "high_evidence": "Proteins with abundant evidence",
            "diverse_tgroups": "Proteins hitting multiple T-groups",
            "minimal_evidence": "Proteins with minimal but clear evidence"
        }
    
    def analyze_batch_sample(self, sample_size: int = 100) -> Dict[str, List]:
        """Analyze a sample of proteins from the batch"""
        
        print(f"Analyzing proteins in: {self.domains_dir}")
        
        if not self.domains_dir.exists():
            print(f"ERROR: Domains directory not found: {self.domains_dir}")
            return {}
        
        # Get domain summary files
        domain_summaries = list(self.domains_dir.glob("*.develop291.domain_summary.xml"))
        print(f"Found {len(domain_summaries)} domain summary files")
        
        # Always include 8ovp_A if it exists
        ovp_file = self.domains_dir / "8ovp_A.develop291.domain_summary.xml"
        if ovp_file.exists():
            self.categorized['validated'].append({
                'protein_id': '8ovp_A',
                'pdb_id': '8ovp',
                'chain_id': 'A',
                'length': 569,  # Known from your tests
                'category': 'validated',
                'description': 'GFP-PBP fusion - validated test case',
                'evidence': {
                    'chain_blast': 1,  # Known to have chain BLAST
                    'domain_blast': 1,
                    'hhsearch': 1,
                    'total': 3
                },
                'tgroups': ['5025.1.1', '7579.1.1'],  # GFP-like and PBP-like
                'estimated_domains': 3
            })

        # Sample proteins
        sample_files = domain_summaries[:sample_size]

        for summary_file in sample_files:
            protein_id = summary_file.stem.replace('.develop291.domain_summary', '')

            if protein_id == '8ovp_A':  # Already handled
                continue

            analysis = self.analyze_protein(summary_file, protein_id)
            if analysis:
                self.proteins_analyzed.append(analysis)
                self.categorize_protein(analysis)

        return dict(self.categorized)

    def analyze_protein(self, summary_file: Path, protein_id: str) -> Dict:
        """Analyze a single protein's domain summary"""

        try:
            tree = ET.parse(summary_file)
            root = tree.getroot()

            # Parse protein ID
            parts = protein_id.split('_')
            pdb_id = parts[0]
            chain_id = parts[1] if len(parts) > 1 else 'A'

            # Count evidence
            chain_blast = len(root.findall(".//chain_blast_run/hits/hit"))
            domain_blast = len(root.findall(".//blast_run/hits/hit"))
            hhsearch = len(root.findall(".//hh_run/hits/hit"))
            total_evidence = chain_blast + domain_blast + hhsearch

            if total_evidence == 0:
                return None  # Skip proteins with no evidence

            # Extract T-groups
            tgroups = set()
            domain_ids = set()

            # From domain BLAST
            for hit in root.findall(".//blast_run/hits/hit"):
                t_group = hit.get("t_group")
                if t_group:
                    tgroups.add(t_group)
                domain_id = hit.get("domain_id")
                if domain_id:
                    domain_ids.add(domain_id)

            # Check for high-confidence hits to estimate domains
            high_conf_domains = set()

            for hit in root.findall(".//blast_run/hits/hit"):
                evalue = float(hit.get("evalues", "999"))
                if evalue < 1e-10:
                    domain_id = hit.get("domain_id", "")
                    if domain_id:
                        # Group by source PDB to avoid counting same domain multiple times
                        high_conf_domains.add(domain_id[1:5] if len(domain_id) > 4 else domain_id)

            for hit in root.findall(".//hh_run/hits/hit"):
                prob = float(hit.get("probability", "0"))
                if prob > 90:
                    hit_id = hit.get("hit_id", "")
                    if hit_id and hit_id.startswith("e"):
                        high_conf_domains.add(hit_id[1:5])

            # Estimate protein length from max hit position
            estimated_length = 0
            for hit in root.findall(".//*/hit"):
                query_reg = hit.find("query_reg")
                if query_reg is not None and query_reg.text:
                    # Parse range to get max position
                    try:
                        parts = query_reg.text.strip().split('-')
                        if len(parts) == 2:
                            end_pos = int(parts[1])
                            estimated_length = max(estimated_length, end_pos)
                    except:
                        pass

            return {
                'protein_id': protein_id,
                'pdb_id': pdb_id,
                'chain_id': chain_id,
                'file': summary_file,
                'evidence': {
                    'chain_blast': chain_blast,
                    'domain_blast': domain_blast,
                    'hhsearch': hhsearch,
                    'total': total_evidence
                },
                'tgroups': tgroups,
                'domain_ids': domain_ids,
                'estimated_domains': len(high_conf_domains) if high_conf_domains else 1,
                'estimated_length': estimated_length
            }

        except Exception as e:
            print(f"Error analyzing {protein_id}: {e}")
            return None

    def categorize_protein(self, analysis: Dict):
        """Categorize protein based on its characteristics"""

        protein_info = {
            'protein_id': analysis['protein_id'],
            'pdb_id': analysis['pdb_id'],
            'chain_id': analysis['chain_id'],
            'length': analysis['estimated_length'],
            'evidence': analysis['evidence'],
            'tgroups': list(analysis['tgroups']),
            'estimated_domains': analysis['estimated_domains']
        }

        # Categorize by evidence pattern and size
        length = analysis['estimated_length']
        domains = analysis['estimated_domains']
        evidence = analysis['evidence']

        # High evidence proteins
        if evidence['total'] > 100:
            protein_info['category'] = 'high_evidence'
            protein_info['description'] = f"High evidence protein ({evidence['total']} hits)"
            self.categorized['high_evidence'].append(protein_info)

        # Multi-domain candidates
        if domains > 1:
            if evidence['chain_blast'] > 0:
                protein_info['category'] = 'chain_blast_multi'
                protein_info['description'] = f"Multi-domain with chain BLAST ({domains} domains)"
                self.categorized['chain_blast_multi'].append(protein_info)
            else:
                protein_info['category'] = 'multi_domain_clear'
                protein_info['description'] = f"Clear multi-domain ({domains} domains)"
                self.categorized['multi_domain_clear'].append(protein_info)

        # Single domain by size
        elif domains == 1:
            if length < 150:
                protein_info['category'] = 'single_domain_small'
                protein_info['description'] = f"Small single domain ({length} residues)"
                self.categorized['single_domain_small'].append(protein_info)
            elif length < 400:
                protein_info['category'] = 'single_domain_medium'
                protein_info['description'] = f"Medium single domain ({length} residues)"
                self.categorized['single_domain_medium'].append(protein_info)

        # Large complex proteins
        if length > 800:
            protein_info['category'] = 'large_complex'
            protein_info['description'] = f"Large protein ({length} residues, {domains} domains)"
            self.categorized['large_complex'].append(protein_info)

        # Diverse T-groups
        if len(analysis['tgroups']) > 2:
            protein_info['category'] = 'diverse_tgroups'
            protein_info['description'] = f"Hits {len(analysis['tgroups'])} different T-groups"
            self.categorized['diverse_tgroups'].append(protein_info)

        # Minimal evidence
        if 0 < evidence['total'] < 20:
            protein_info['category'] = 'minimal_evidence'
            protein_info['description'] = f"Minimal evidence ({evidence['total']} hits)"
            self.categorized['minimal_evidence'].append(protein_info)

    def select_test_proteins(self, max_count: int = 15) -> List[Dict]:
        """Select diverse test proteins from categorized results"""

        selected = []
        used_proteins = set()

        # Selection priorities
        priorities = [
            ('validated', 1),           # Keep 8ovp_A
            ('chain_blast_multi', 3),   # Good decomposition tests
            ('single_domain_small', 2), # Simple cases
            ('single_domain_medium', 2), # Simple cases
            ('multi_domain_clear', 2),  # Multi-domain without chain BLAST
            ('large_complex', 1),       # Large protein test
            ('high_evidence', 2),       # Well-supported proteins
            ('diverse_tgroups', 1),     # Complex classification
            ('minimal_evidence', 1)     # Edge case
        ]

        for category, target in priorities:
            if category in self.categorized:
                candidates = self.categorized[category]

                if not candidates:  # Skip empty categories
                    continue

                # Sort by evidence quality
                sorted_candidates = sorted(
                    candidates,
                    key=lambda x: x.get('evidence', {}).get('total', 0),
                    reverse=True
                )

                added = 0
                for candidate in sorted_candidates:
                    if candidate['protein_id'] not in used_proteins and added < target:
                        selected.append(candidate)
                        used_proteins.add(candidate['protein_id'])
                        added += 1

                if len(selected) >= max_count:
                    break

        return selected[:max_count]

    def generate_test_config(self, proteins: List[Dict], output_file: str):
        """Generate test configuration"""

        lines = [
            '"""',
            'Test proteins discovered from batch analysis',
            'Auto-generated by expand_test_suite_now.py',
            '"""',
            '',
            '# Proteins discovered from actual batch content',
            'BATCH_TEST_PROTEINS = {'
        ]

        for protein in proteins:
            evidence = protein.get('evidence', {})
            lines.extend([
                f'    "{protein["protein_id"]}": {{',
                f'        "category": "{protein.get("category", "unknown")}",',
                f'        "description": "{protein.get("description", "")}",',
                f'        "estimated_length": {protein.get("length", 0)},',
                f'        "estimated_domains": {protein.get("estimated_domains", 1)},',
                f'        "evidence_total": {evidence.get("total", 0)},',
                f'        "chain_blast": {evidence.get("chain_blast", 0)},',
                f'        "domain_blast": {evidence.get("domain_blast", 0)},',
                f'        "hhsearch": {evidence.get("hhsearch", 0)},',
                f'        "tgroups": {protein.get("tgroups", [])}',
                '    },'
            ])

        lines.extend([
            '}',
            '',
            '# Categories found in batch',
            'BATCH_CATEGORIES = {'
        ])

        # Group by category
        by_category = defaultdict(list)
        for protein in proteins:
            by_category[protein.get('category', 'unknown')].append(protein['protein_id'])

        for category, protein_ids in by_category.items():
            lines.append(f'    "{category}": {protein_ids},')

        lines.append('}')

        with open(output_file, 'w') as f:
            f.write('\n'.join(lines))

        print(f"\nGenerated test configuration: {output_file}")

    def generate_test_cases(self, proteins: List[Dict], output_file: str):
        """Generate test cases for discovered proteins"""

        lines = [
            '#!/usr/bin/env python3',
            '"""',
            'Test cases for proteins discovered from batch analysis',
            '"""',
            '',
            'import pytest',
            'from pathlib import Path',
            'import sys',
            '',
            'sys.path.insert(0, str(Path(__file__).parent.parent.parent))',
            '',
            'from mini.batch_test_proteins import BATCH_TEST_PROTEINS',
            '',
            '',
            'class TestBatchProteins:',
            '    """Test cases for proteins from batch"""',
            '',
            '    @pytest.fixture',
            '    def run_mini_algorithm(self, stable_batch_dir, real_reference_data, blast_alignments):',
            '        """Run mini algorithm on a protein"""',
            '        def _run(protein_id):'
        ]

        # Add the algorithm runner (similar to your existing fixture)
        lines.extend([
            '            from mini.parser import parse_domain_summary',
            '            from mini.partitioner import partition_domains',
            '            import os',
            '            ',
            '            xml_path = os.path.join(stable_batch_dir, "domains", ',
            '                                   f"{protein_id}.develop291.domain_summary.xml")',
            '            ',
            '            if not os.path.exists(xml_path):',
            '                return {"success": False, "error": "Domain summary not found"}',
            '            ',
            '            try:',
            '                evidence = parse_domain_summary(',
            '                    xml_path,',
            '                    reference_lengths=real_reference_data.get("domain_lengths", {}),',
            '                    protein_lengths=real_reference_data.get("protein_lengths", {}),',
            '                    blast_alignments=blast_alignments,',
            '                    require_reference_lengths=False',
            '                )',
            '                ',
            '                if not evidence:',
            '                    return {"success": False, "error": "No evidence found"}',
            '                ',
            '                max_pos = max(ev.query_range.segments[-1].end for ev in evidence)',
            '                sequence_length = int(max_pos * 1.1)',
            '                ',
            '                domains = partition_domains(',
            '                    evidence,',
            '                    sequence_length=sequence_length,',
            '                    domain_definitions=real_reference_data.get("domain_definitions", {})',
            '                )',
            '                ',
            '                return {',
            '                    "success": True,',
            '                    "domains": domains,',
            '                    "sequence_length": sequence_length,',
            '                    "evidence_count": len(evidence)',
            '                }',
            '            except Exception as e:',
            '                return {"success": False, "error": str(e)}',
            '        ',
            '        return _run',
            ''
        ])

        # Group proteins by category for parametrized tests
        by_category = defaultdict(list)
        for protein in proteins:
            if protein.get('category') != 'validated':  # Skip 8ovp_A
                by_category[protein.get('category', 'unknown')].append(protein['protein_id'])

        # Generate parametrized tests
        for category, protein_ids in by_category.items():
            if not protein_ids:
                continue

            lines.extend([
                f'    @pytest.mark.integration',
                f'    @pytest.mark.parametrize("protein_id", {protein_ids})',
                f'    def test_{category}_proteins(self, protein_id, run_mini_algorithm):',
                f'        """Test {category.replace("_", " ")} proteins"""',
                '        info = BATCH_TEST_PROTEINS[protein_id]',
                '        result = run_mini_algorithm(protein_id)',
                '        ',
                '        # Basic success check',
                '        assert result["success"], f"Failed: {result.get(\'error\')}"',
                '        assert len(result["domains"]) > 0, "Should find at least one domain"',
                '        ',
                '        # Category-specific checks',
            ])

            if category == 'single_domain_small' or category == 'single_domain_medium':
                lines.extend([
                    '        # Single domain proteins should have 1 domain',
                    '        assert len(result["domains"]) == 1, \\',
                    '            f"Expected 1 domain, found {len(result[\'domains\'])}"'
                ])
            elif 'multi_domain' in category:
                lines.extend([
                    '        # Multi-domain proteins should have multiple domains',
                    '        assert len(result["domains"]) >= 2, \\',
                    '            f"Expected multiple domains, found {len(result[\'domains\'])}"'
                ])

            lines.append('')

        with open(output_file, 'w') as f:
            f.write('\n'.join(lines))

        print(f"Generated test cases: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Quick expansion of test suite from batch content'
    )

    parser.add_argument('--batch-dir',
                       default='/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424',
                       help='Batch directory to analyze')
    parser.add_argument('--sample-size', type=int, default=100,
                       help='Number of proteins to analyze')
    parser.add_argument('--generate-tests', action='store_true',
                       help='Generate test configuration and cases')

    args = parser.parse_args()

    analyzer = BatchProteinAnalyzer(args.batch_dir)

    # Analyze batch
    print(f"Analyzing proteins from batch...")
    categorized = analyzer.analyze_batch_sample(args.sample_size)

    if not categorized:
        print("No suitable proteins found!")
        return

    # Print summary
    print(f"\nAnalyzed {len(analyzer.proteins_analyzed)} proteins")
    print("\nProteins by category:")
    for category, proteins in categorized.items():
        print(f"  {category}: {len(proteins)} proteins")

    # Select test proteins
    selected = analyzer.select_test_proteins()

    print(f"\nSelected {len(selected)} test proteins:")
    for protein in selected:
        evidence = protein.get('evidence', {})
        print(f"\n{protein['protein_id']}:")
        print(f"  Category: {protein.get('category', 'unknown')}")
        print(f"  Description: {protein.get('description', '')}")
        print(f"  Evidence: CB={evidence.get('chain_blast', 0)}, "
              f"DB={evidence.get('domain_blast', 0)}, "
              f"HH={evidence.get('hhsearch', 0)}")
    
    # Generate test files if requested
    if args.generate_tests:
        print("\nGenerating test files...")
        
        # Create output directory
        output_dir = Path("mini")
        output_dir.mkdir(exist_ok=True)
        
        # Generate config
        analyzer.generate_test_config(
            selected,
            output_dir / "batch_test_proteins.py"
        )
        
        # Generate tests
        test_dir = output_dir / "tests"
        test_dir.mkdir(exist_ok=True)
        
        analyzer.generate_test_cases(
            selected,
            test_dir / "test_batch_proteins.py"
        )
        
        print("\n✅ Test files generated!")
        print("\nTo run tests:")
        print("  cd mini")
        print("  python -m pytest tests/test_batch_proteins.py -v")


if __name__ == "__main__":
    main()
