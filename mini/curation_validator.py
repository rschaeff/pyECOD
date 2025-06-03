#!/usr/bin/env python3
"""
Domain Curation Validator for mini_pyecod

This script helps validate manually curated domain boundaries against algorithm results
and generates regression tests for the validated proteins.

Usage:
    python curation_validator.py --parse-curation        # Parse curation files
    python curation_validator.py --validate-boundaries   # Compare algorithm vs manual
    python curation_validator.py --generate-tests        # Create regression tests
    python curation_validator.py --check-progress        # Show curation progress
"""

import os
import sys
import re
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple, NamedTuple
from dataclasses import dataclass
import xml.etree.ElementTree as ET

# Add paths
sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.batch_test_proteins import BATCH_TEST_PROTEINS

class ManualDomain(NamedTuple):
    """Manually curated domain"""
    range_str: str
    family: str
    notes: str

@dataclass
class CurationData:
    """Parsed curation data for a protein"""
    protein_id: str
    curation_status: str  # 'pending', 'curated', 'validated'
    manual_domains: List[ManualDomain]
    algorithm_domains: List[Dict]
    validation_notes: str = ""
    boundaries_acceptable: bool = False
    ready_for_testing: bool = False

@dataclass
class BoundaryComparison:
    """Comparison between algorithm and manual boundaries"""
    protein_id: str
    domain_count_match: bool
    domain_count_algorithm: int
    domain_count_manual: int
    boundary_accuracy: float  # 0.0 to 1.0
    details: List[Dict]

class CurationValidator:
    """Validates manually curated domain boundaries"""
    
    def __init__(self, curation_dir: str = "/tmp/test_suite_expansion/curation"):
        self.curation_dir = Path(curation_dir)
        self.output_dir = Path(curation_dir).parent / "validation"
        self.output_dir.mkdir(exist_ok=True)
        
        # Parsed curation data
        self.curations: Dict[str, CurationData] = {}
        
    def parse_curation_files(self) -> Dict[str, CurationData]:
        """Parse all curation markdown files"""
        
        print("📖 Parsing curation files...")
        
        if not self.curation_dir.exists():
            print(f"ERROR: Curation directory not found: {self.curation_dir}")
            return {}
        
        curation_files = list(self.curation_dir.glob("*_curation.md"))
        
        if not curation_files:
            print(f"No curation files found in {self.curation_dir}")
            return {}
        
        parsed = {}
        
        for curation_file in curation_files:
            protein_id = curation_file.stem.replace('_curation', '')
            
            try:
                curation_data = self._parse_single_curation(curation_file, protein_id)
                if curation_data:
                    parsed[protein_id] = curation_data
                    status = "curated" if curation_data.manual_domains else "pending"
                    print(f"  ✓ {protein_id}: {status} ({len(curation_data.manual_domains)} manual domains)")
                else:
                    print(f"  ⚠️  {protein_id}: Could not parse")
            except Exception as e:
                print(f"  ✗ {protein_id}: Error - {e}")
        
        self.curations = parsed
        print(f"\nParsed {len(parsed)} curation files")
        return parsed
    
    def _parse_single_curation(self, file_path: Path, protein_id: str) -> Optional[CurationData]:
        """Parse a single curation markdown file"""
        
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Parse manual domains from markdown
        manual_domains = []
        
        # Look for manual domain sections
        manual_section = re.search(r'### Expected Domains.*?(?=###|$)', content, re.DOTALL)
        if manual_section:
            domain_matches = re.finditer(
                r'\d+\.\s*\*\*Domain\s+\d+\*\*.*?Range:\s*([^\n]+).*?Family/Function:\s*([^\n]+).*?Notes:\s*([^\n]+)',
                manual_section.group(0),
                re.DOTALL
            )
            
            for match in domain_matches:
                range_str = match.group(1).strip()
                family = match.group(2).strip()
                notes = match.group(3).strip()
                
                # Skip TODO entries
                if range_str != "TODO" and family != "TODO":
                    manual_domains.append(ManualDomain(range_str, family, notes))
        
        # Parse validation status checkboxes
        boundaries_acceptable = '- [x] Algorithm boundaries acceptable' in content
        ready_for_testing = '- [x] Test case ready for regression testing' in content
        
        # Parse algorithm domains (from the file content)
        algorithm_domains = []
        alg_section = re.search(r'## Domain Boundaries \(Algorithm\).*?(?=##|$)', content, re.DOTALL)
        if alg_section:
            domain_sections = re.finditer(
                r'### Domain \d+: (\w+).*?Range\*\*:\s*([^\n]+).*?Family\*\*:\s*([^\n]+).*?Source\*\*:\s*([^\n]+).*?Size\*\*:\s*(\d+)',
                alg_section.group(0),
                re.DOTALL
            )
            
            for match in domain_sections:
                algorithm_domains.append({
                    'id': match.group(1),
                    'range': match.group(2).strip(),
                    'family': match.group(3).strip(),
                    'source': match.group(4).strip(),
                    'size': int(match.group(5))
                })
        
        # Determine curation status
        if ready_for_testing and manual_domains:
            status = 'validated'
        elif manual_domains:
            status = 'curated'
        else:
            status = 'pending'
        
        return CurationData(
            protein_id=protein_id,
            curation_status=status,
            manual_domains=manual_domains,
            algorithm_domains=algorithm_domains,
            boundaries_acceptable=boundaries_acceptable,
            ready_for_testing=ready_for_testing
        )
    
    def validate_boundaries(self) -> Dict[str, BoundaryComparison]:
        """Compare algorithm boundaries with manual curation"""
        
        if not self.curations:
            self.parse_curation_files()
        
        print("🔍 Validating boundaries against manual curation...")
        
        comparisons = {}
        
        for protein_id, curation in self.curations.items():
            if curation.curation_status in ['curated', 'validated']:
                comparison = self._compare_boundaries(curation)
                comparisons[protein_id] = comparison
                
                accuracy = comparison.boundary_accuracy
                match_indicator = "✓" if comparison.domain_count_match else "✗"
                
                print(f"  {match_indicator} {protein_id}: {accuracy:.1%} boundary accuracy "
                      f"({comparison.domain_count_algorithm} alg vs {comparison.domain_count_manual} manual)")
        
        # Save comparison results
        self._save_boundary_comparisons(comparisons)
        
        return comparisons
    
    def _compare_boundaries(self, curation: CurationData) -> BoundaryComparison:
        """Compare algorithm vs manual boundaries for one protein"""
        
        alg_domains = curation.algorithm_domains
        manual_domains = curation.manual_domains
        
        domain_count_match = len(alg_domains) == len(manual_domains)
        
        details = []
        total_accuracy = 0.0
        
        # Compare each algorithm domain with closest manual domain
        for i, alg_domain in enumerate(alg_domains):
            alg_range = alg_domain['range']
            alg_positions = self._parse_range_positions(alg_range)
            
            best_match = None
            best_overlap = 0.0
            
            for j, manual_domain in enumerate(manual_domains):
                manual_positions = self._parse_range_positions(manual_domain.range_str)
                
                if alg_positions and manual_positions:
                    overlap = len(alg_positions & manual_positions)
                    union = len(alg_positions | manual_positions)
                    jaccard = overlap / union if union > 0 else 0.0
                    
                    if jaccard > best_overlap:
                        best_overlap = jaccard
                        best_match = j
            
            details.append({
                'algorithm_domain': i,
                'algorithm_range': alg_range,
                'manual_match': best_match,
                'manual_range': manual_domains[best_match].range_str if best_match is not None else None,
                'jaccard_similarity': best_overlap
            })
            
            total_accuracy += best_overlap
        
        # Average accuracy across domains
        boundary_accuracy = total_accuracy / len(alg_domains) if alg_domains else 0.0
        
        return BoundaryComparison(
            protein_id=curation.protein_id,
            domain_count_match=domain_count_match,
            domain_count_algorithm=len(alg_domains),
            domain_count_manual=len(manual_domains),
            boundary_accuracy=boundary_accuracy,
            details=details
        )
    
    def _parse_range_positions(self, range_str: str) -> set:
        """Parse range string to set of positions"""
        positions = set()
        
        try:
            for segment in range_str.split(','):
                segment = segment.strip()
                if '-' in segment:
                    start, end = map(int, segment.split('-'))
                    positions.update(range(start, end + 1))
                else:
                    positions.add(int(segment))
        except:
            pass  # Invalid range format
        
        return positions
    
    def _save_boundary_comparisons(self, comparisons: Dict[str, BoundaryComparison]):
        """Save boundary comparison results"""
        
        # Convert to JSON-serializable format
        serializable = {}
        for protein_id, comparison in comparisons.items():
            serializable[protein_id] = {
                'protein_id': comparison.protein_id,
                'domain_count_match': comparison.domain_count_match,
                'domain_count_algorithm': comparison.domain_count_algorithm,
                'domain_count_manual': comparison.domain_count_manual,
                'boundary_accuracy': comparison.boundary_accuracy,
                'details': comparison.details
            }
        
        output_file = self.output_dir / "boundary_comparisons.json"
        with open(output_file, 'w') as f:
            json.dump(serializable, f, indent=2)
        
        print(f"\n💾 Boundary comparisons saved to: {output_file}")
    
    def check_curation_progress(self):
        """Show curation progress summary"""
        
        if not self.curations:
            self.parse_curation_files()
        
        print("📊 CURATION PROGRESS SUMMARY")
        print("=" * 60)
        
        total = len(self.curations)
        pending = sum(1 for c in self.curations.values() if c.curation_status == 'pending')
        curated = sum(1 for c in self.curations.values() if c.curation_status == 'curated')
        validated = sum(1 for c in self.curations.values() if c.curation_status == 'validated')
        
        print(f"Total proteins: {total}")
        print(f"Pending curation: {pending}")
        print(f"Manually curated: {curated}")
        print(f"Validated & test-ready: {validated}")
        
        if total > 0:
            print(f"Progress: {(curated + validated)/total*100:.1f}% curated")
            print(f"Test-ready: {validated/total*100:.1f}%")
        
        # Show by category
        print(f"\nBy category:")
        
        by_category = {}
        for protein_id, curation in self.curations.items():
            info = BATCH_TEST_PROTEINS.get(protein_id, {})
            category = info.get('category', 'unknown')
            if category not in by_category:
                by_category[category] = {'pending': 0, 'curated': 0, 'validated': 0}
            by_category[category][curation.curation_status] += 1
        
        for category, counts in by_category.items():
            total_cat = sum(counts.values())
            curated_cat = counts['curated'] + counts['validated']
            print(f"  {category}: {curated_cat}/{total_cat} curated")
        
        # Next steps
        if pending > 0:
            print(f"\n🎯 Next Steps:")
            print(f"1. Complete curation of {pending} pending proteins")
            print(f"2. Fill in manual domain boundaries in curation files")
            print(f"3. Mark proteins as test-ready when boundaries are validated")
        
        if validated > 0:
            print(f"\n✅ Ready to generate regression tests for {validated} validated proteins")
    
    def generate_regression_tests(self):
        """Generate regression test files for validated proteins"""
        
        if not self.curations:
            self.parse_curation_files()
        
        validated_proteins = [
            (protein_id, curation) for protein_id, curation in self.curations.items()
            if curation.curation_status == 'validated' and curation.ready_for_testing
        ]
        
        if not validated_proteins:
            print("No validated proteins ready for testing!")
            print("Complete manual curation and mark proteins as test-ready first.")
            return
        
        print(f"🧪 Generating regression tests for {len(validated_proteins)} validated proteins...")
        
        # Generate test file
        test_file = self.output_dir / "test_curated_regression.py"
        self._create_regression_test_file(validated_proteins, test_file)
        
        # Generate expected results file
        expected_file = self.output_dir / "expected_domain_boundaries.json"
        self._create_expected_boundaries_file(validated_proteins, expected_file)
        
        print(f"✅ Generated regression test infrastructure:")
        print(f"  Test file: {test_file}")
        print(f"  Expected boundaries: {expected_file}")
        print(f"\nTo run regression tests:")
        print(f"  cd mini")
        print(f"  python -m pytest {test_file} -v")
    
    def _create_regression_test_file(self, validated_proteins: List[Tuple[str, CurationData]], output_file: Path):
        """Create pytest regression test file"""
        
        protein_ids = [protein_id for protein_id, _ in validated_proteins]
        
        lines = [
            '#!/usr/bin/env python3',
            '"""',
            'Regression tests for manually curated domain boundaries',
            '',
            'These tests validate that the mini_pyecod algorithm produces',
            'domain boundaries that match manually curated expectations.',
            '"""',
            '',
            'import pytest',
            'import json',
            'from pathlib import Path',
            'import sys',
            '',
            'sys.path.insert(0, str(Path(__file__).parent.parent))',
            '',
            'from mini.parser import parse_domain_summary',
            'from mini.partitioner import partition_domains',
            'from mini.decomposer import load_domain_definitions',
            '',
            '',
            'class TestCuratedRegression:',
            '    """Regression tests against manually curated boundaries"""',
            '',
            '    @pytest.fixture(scope="class")',
            '    def expected_boundaries(self):',
            '        """Load expected domain boundaries from curation"""',
            f'        expected_file = Path(__file__).parent / "expected_domain_boundaries.json"',
            '        with open(expected_file, "r") as f:',
            '            return json.load(f)',
            '',
            '    @pytest.fixture',
            '    def algorithm_runner(self, stable_batch_dir, real_reference_data, blast_alignments):',
            '        """Run the algorithm on a protein"""',
            '        def _run(protein_id):',
            '            from pathlib import Path',
            '            import os',
            '            ',
            '            xml_path = os.path.join(stable_batch_dir, "domains",',
            '                                   f"{protein_id}.develop291.domain_summary.xml")',
            '            ',
            '            evidence = parse_domain_summary(',
            '                xml_path,',
            '                reference_lengths=real_reference_data.get("domain_lengths", {}),',
            '                protein_lengths=real_reference_data.get("protein_lengths", {}),',
            '                blast_alignments=blast_alignments,',
            '                require_reference_lengths=True',
            '            )',
            '            ',
            '            max_pos = max(ev.query_range.segments[-1].end for ev in evidence)',
            '            sequence_length = int(max_pos * 1.1)',
            '            ',
            '            domains = partition_domains(',
            '                evidence,',
            '                sequence_length=sequence_length,',
            '                domain_definitions=real_reference_data.get("domain_definitions", {})',
            '            )',
            '            ',
            '            return domains',
            '        return _run',
            '',
            f'    @pytest.mark.integration',
            f'    @pytest.mark.parametrize("protein_id", {protein_ids})',
            '    def test_curated_domain_boundaries(self, protein_id, expected_boundaries, algorithm_runner):',
            '        """Test that algorithm matches manually curated boundaries"""',
            '        ',
            '        # Get expected boundaries',
            '        expected = expected_boundaries[protein_id]',
            '        expected_count = len(expected["domains"])',
            '        ',
            '        # Run algorithm',
            '        algorithm_domains = algorithm_runner(protein_id)',
            '        algorithm_count = len(algorithm_domains)',
            '        ',
            '        # Test domain count',
            '        assert algorithm_count == expected_count, \\',
            '            f"Domain count mismatch: algorithm={algorithm_count}, expected={expected_count}"',
            '        ',
            '        # Test boundary accuracy',
            '        total_accuracy = 0.0',
            '        ',
            '        for i, alg_domain in enumerate(algorithm_domains):',
            '            alg_positions = set(alg_domain.range.to_positions_simple())',
            '            ',
            '            best_jaccard = 0.0',
            '            for exp_domain in expected["domains"]:',
            '                exp_positions = self._parse_range_positions(exp_domain["range"])',
            '                ',
            '                if alg_positions and exp_positions:',
            '                    overlap = len(alg_positions & exp_positions)',
            '                    union = len(alg_positions | exp_positions)',
            '                    jaccard = overlap / union if union > 0 else 0.0',
            '                    best_jaccard = max(best_jaccard, jaccard)',
            '            ',
            '            total_accuracy += best_jaccard',
            '        ',
            '        # Average accuracy should be high',
            '        average_accuracy = total_accuracy / algorithm_count if algorithm_count > 0 else 0.0',
            '        min_accuracy = expected.get("min_boundary_accuracy", 0.8)',
            '        ',
            '        assert average_accuracy >= min_accuracy, \\',
            '            f"Boundary accuracy {average_accuracy:.2%} below threshold {min_accuracy:.2%}"',
            '',
            '    def _parse_range_positions(self, range_str: str) -> set:',
            '        """Parse range string to set of positions"""',
            '        positions = set()',
            '        try:',
            '            for segment in range_str.split(","):',
            '                segment = segment.strip()',
            '                if "-" in segment:',
            '                    start, end = map(int, segment.split("-"))',
            '                    positions.update(range(start, end + 1))',
            '                else:',
            '                    positions.add(int(segment))',
            '        except:',
            '            pass',
            '        return positions',
            '',
            '    @pytest.mark.integration',
            '    def test_overall_regression_performance(self, expected_boundaries, algorithm_runner):',
            '        """Test overall performance across all curated proteins"""',
            '        ',
            '        correct_counts = 0',
            '        total_accuracy = 0.0',
            '        total_proteins = len(expected_boundaries)',
            '        ',
            '        for protein_id, expected in expected_boundaries.items():',
            '            algorithm_domains = algorithm_runner(protein_id)',
            '            ',
            '            # Check domain count accuracy',
            '            if len(algorithm_domains) == len(expected["domains"]):',
            '                correct_counts += 1',
            '            ',
            '            # Calculate boundary accuracy for this protein',
            '            protein_accuracy = 0.0',
            '            if algorithm_domains:',
            '                for alg_domain in algorithm_domains:',
            '                    alg_positions = set(alg_domain.range.to_positions_simple())',
            '                    best_jaccard = 0.0',
            '                    ',
            '                    for exp_domain in expected["domains"]:',
            '                        exp_positions = self._parse_range_positions(exp_domain["range"])',
            '                        if alg_positions and exp_positions:',
            '                            overlap = len(alg_positions & exp_positions)',
            '                            union = len(alg_positions | exp_positions)',
            '                            jaccard = overlap / union if union > 0 else 0.0',
            '                            best_jaccard = max(best_jaccard, jaccard)',
            '                    ',
            '                    protein_accuracy += best_jaccard',
            '                ',
            '                protein_accuracy /= len(algorithm_domains)',
            '            ',
            '            total_accuracy += protein_accuracy',
            '        ',
            '        # Overall metrics',
            '        count_accuracy = correct_counts / total_proteins',
            '        boundary_accuracy = total_accuracy / total_proteins',
            '        ',
            '        print(f"\\nRegression Test Summary:")',
            '        print(f"  Proteins tested: {total_proteins}")',
            '        print(f"  Domain count accuracy: {count_accuracy:.1%}")',
            '        print(f"  Average boundary accuracy: {boundary_accuracy:.1%}")',
            '        ',
            '        # Performance thresholds',
            '        assert count_accuracy >= 0.8, f"Domain count accuracy {count_accuracy:.1%} too low"',
            '        assert boundary_accuracy >= 0.8, f"Boundary accuracy {boundary_accuracy:.1%} too low"'
        ]
        
        with open(output_file, 'w') as f:
            f.write('\n'.join(lines))
    
    def _create_expected_boundaries_file(self, validated_proteins: List[Tuple[str, CurationData]], output_file: Path):
        """Create expected boundaries JSON file"""
        
        expected = {}
        
        for protein_id, curation in validated_proteins:
            domains = []
            for manual_domain in curation.manual_domains:
                domains.append({
                    'range': manual_domain.range_str,
                    'family': manual_domain.family,
                    'notes': manual_domain.notes
                })
            
            expected[protein_id] = {
                'domains': domains,
                'min_boundary_accuracy': 0.8,  # Default threshold
                'notes': f'Manually curated boundaries for {protein_id}'
            }
        
        with open(output_file, 'w') as f:
            json.dump(expected, f, indent=2)


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Validate manually curated domain boundaries',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Workflow:
1. python curation_validator.py --check-progress     # See curation status
2. [Complete manual curation of domain boundaries]
3. python curation_validator.py --validate-boundaries # Check accuracy
4. python curation_validator.py --generate-tests      # Create regression tests

Examples:
    python curation_validator.py --parse-curation     # Parse curation files
    python curation_validator.py --check-progress     # Show progress summary
    python curation_validator.py --validate-boundaries # Compare boundaries
    python curation_validator.py --generate-tests     # Create test files
        """
    )
    
    parser.add_argument('--parse-curation', action='store_true',
                       help='Parse curation markdown files')
    parser.add_argument('--check-progress', action='store_true',
                       help='Show curation progress summary')
    parser.add_argument('--validate-boundaries', action='store_true',
                       help='Compare algorithm vs manual boundaries')
    parser.add_argument('--generate-tests', action='store_true',
                       help='Generate regression test files')
    parser.add_argument('--curation-dir', 
                       default='/tmp/test_suite_expansion/curation',
                       help='Directory containing curation files')
    
    args = parser.parse_args()
    
    validator = CurationValidator(args.curation_dir)
    
    if args.parse_curation:
        validator.parse_curation_files()
        return
        
    if args.check_progress:
        validator.check_curation_progress()
        return
        
    if args.validate_boundaries:
        comparisons = validator.validate_boundaries()
        
        if comparisons:
            # Summary statistics
            total = len(comparisons)
            exact_matches = sum(1 for c in comparisons.values() if c.domain_count_match)
            avg_accuracy = sum(c.boundary_accuracy for c in comparisons.values()) / total
            
            print(f"\n📊 VALIDATION SUMMARY")
            print(f"Proteins validated: {total}")
            print(f"Exact domain count matches: {exact_matches}/{total} ({exact_matches/total:.1%})")
            print(f"Average boundary accuracy: {avg_accuracy:.1%}")
            
            # Show proteins ready for testing
            high_accuracy = [p for p, c in comparisons.items() if c.boundary_accuracy >= 0.8]
            print(f"High accuracy proteins (≥80%): {len(high_accuracy)}")
            
            if high_accuracy:
                print(f"\n✅ Proteins ready for regression testing:")
                for protein_id in high_accuracy:
                    accuracy = comparisons[protein_id].boundary_accuracy
                    print(f"  {protein_id}: {accuracy:.1%} accuracy")
        
        return
        
    if args.generate_tests:
        validator.generate_regression_tests()
        return
    
    # No action specified
    parser.print_help()
    print(f"\n💡 Tip: Start with --check-progress to see current status")

if __name__ == "__main__":
    main()
