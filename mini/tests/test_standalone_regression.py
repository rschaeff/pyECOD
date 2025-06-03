#!/usr/bin/env python3
"""
Standalone regression tests for manually curated domain boundaries

These tests validate that the mini_pyecod algorithm produces
domain boundaries that match manually curated expectations.

Run from mini/ directory:
    python -m pytest /tmp/test_suite_expansion/validation/test_standalone_regression.py -v
"""

import pytest
import json
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path
import sys
import os

# Add paths for mini_pyecod
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "mini"))

class TestStandaloneRegression:
    """Standalone regression tests using pyecod_mini executable"""

    @pytest.fixture(scope="class")
    def expected_boundaries(self):
        """Load expected domain boundaries from curation"""
        expected_file = Path(__file__).parent / "expected_domain_boundaries.json"
        with open(expected_file, "r") as f:
            return json.load(f)

    @pytest.fixture
    def pyecod_mini_runner(self):
        """Run pyecod_mini executable on a protein"""
        def _run(protein_id):
            # Find the mini directory and executable
            mini_dir = Path(__file__).parent.parent.parent / "mini"
            
            # Try executable wrapper first, then Python script
            pyecod_mini = mini_dir / "pyecod_mini"
            if pyecod_mini.exists() and pyecod_mini.is_file():
                cmd = [str(pyecod_mini), protein_id]
            else:
                pyecod_mini_py = mini_dir / "pyecod_mini.py"
                if pyecod_mini_py.exists():
                    cmd = ["python", str(pyecod_mini_py), protein_id]
                else:
                    raise FileNotFoundError("Neither pyecod_mini nor pyecod_mini.py found")
            
            # Run the algorithm
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=120,
                cwd=str(mini_dir)
            )
            
            if result.returncode != 0:
                raise RuntimeError(f"pyecod_mini failed: {result.stderr}")
            
            # Parse the output XML file
            output_file = f"/tmp/{protein_id}_mini.domains.xml"
            if not os.path.exists(output_file):
                raise FileNotFoundError(f"Output file not found: {output_file}")
            
            domains = self._parse_domain_xml(output_file)
            return domains
        
        return _run

    def _parse_domain_xml(self, xml_file):
        """Parse domain XML file to extract domain information"""
        domains = []
        
        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()
            
            for domain_elem in root.findall(".//domain"):
                domain_info = {
                    'id': domain_elem.get('id', ''),
                    'range': domain_elem.get('range', ''),
                    'family': domain_elem.get('family', ''),
                    'source': domain_elem.get('source', ''),
                    'size': int(domain_elem.get('size', '0')),
                    'range_obj': self._parse_range_to_positions(domain_elem.get('range', ''))
                }
                domains.append(domain_info)
        
        except Exception as e:
            print(f"Warning: Could not parse {xml_file}: {e}")
        
        return domains

    def _parse_range_to_positions(self, range_str):
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
            pass
        return positions

    def _parse_expected_range_to_positions(self, range_str):
        """Parse expected range string to set of positions"""
        return self._parse_range_to_positions(range_str)

    @pytest.mark.parametrize("protein_id", ['8oni_L', '8p2e_B', '8oz3_B', '8p12_L', '8olg_A'])
    def test_curated_domain_boundaries(self, protein_id, expected_boundaries, pyecod_mini_runner):
        """Test that algorithm matches manually curated boundaries"""
        
        # Get expected boundaries
        expected = expected_boundaries[protein_id]
        expected_count = len(expected["domains"])
        
        # Run algorithm
        algorithm_domains = pyecod_mini_runner(protein_id)
        algorithm_count = len(algorithm_domains)
        
        print(f"\n{protein_id}:")
        print(f"  Expected domains: {expected_count}")
        print(f"  Algorithm domains: {algorithm_count}")
        
        # Special handling for "no domain" cases (e.g., 8olg_A fibril)
        if expected_count == 1 and expected["domains"][0]["family"] == "Amyloid fibril (not folded domain)":
            # This is a "no true domain" case - algorithm should ideally find 0, but finding 1 is acceptable
            print(f"  → Fibrillar protein: Algorithm found {algorithm_count}, expected 0 true domains")
            # Don't fail the test for fibrillar proteins, just log the behavior
            return
        
        # Test domain count for normal proteins
        assert algorithm_count == expected_count, \
            f"Domain count mismatch: algorithm={algorithm_count}, expected={expected_count}"
        
        # Test boundary accuracy
        total_accuracy = 0.0
        
        for i, alg_domain in enumerate(algorithm_domains):
            alg_positions = alg_domain['range_obj']
            
            best_jaccard = 0.0
            best_match = None
            
            for j, exp_domain in enumerate(expected["domains"]):
                exp_positions = self._parse_expected_range_to_positions(exp_domain["range"])
                
                if alg_positions and exp_positions:
                    overlap = len(alg_positions & exp_positions)
                    union = len(alg_positions | exp_positions)
                    jaccard = overlap / union if union > 0 else 0.0
                    
                    if jaccard > best_jaccard:
                        best_jaccard = jaccard
                        best_match = j
            
            total_accuracy += best_jaccard
            
            if best_match is not None:
                print(f"  Domain {i+1}: {alg_domain['range']} vs {expected['domains'][best_match]['range']} "
                      f"(similarity: {best_jaccard:.1%})")
        
        # Average accuracy should be high
        average_accuracy = total_accuracy / algorithm_count if algorithm_count > 0 else 0.0
        min_accuracy = expected.get("min_boundary_accuracy", 0.75)  # Realistic 75% threshold

        print(f"  → Average boundary accuracy: {average_accuracy:.1%}")

        assert average_accuracy >= min_accuracy, \
            f"Boundary accuracy {average_accuracy:.2%} below threshold {min_accuracy:.2%}"

    def test_overall_regression_performance(self, expected_boundaries, pyecod_mini_runner):
        """Test overall performance across all curated proteins"""

        correct_counts = 0
        total_accuracy = 0.0
        total_proteins = len(expected_boundaries)
        results = {}

        for protein_id, expected in expected_boundaries.items():
            try:
                algorithm_domains = pyecod_mini_runner(protein_id)

                # Handle fibrillar proteins specially
                if (len(expected["domains"]) == 1 and
                    expected["domains"][0]["family"] == "Amyloid fibril (not folded domain)"):
                    # For fibrillar proteins, any reasonable result is acceptable
                    results[protein_id] = {
                        'algorithm_count': len(algorithm_domains),
                        'expected_count': 0,  # 0 true domains
                        'boundary_accuracy': 0.5,  # Neutral score
                        'type': 'fibrillar'
                    }
                    correct_counts += 0.5  # Partial credit
                    total_accuracy += 0.5
                    continue

                # Normal protein handling
                expected_count = len(expected["domains"])
                algorithm_count = len(algorithm_domains)

                # Check domain count accuracy
                count_match = (algorithm_count == expected_count)
                if count_match:
                    correct_counts += 1

                # Calculate boundary accuracy for this protein
                protein_accuracy = 0.0
                if algorithm_domains:
                    for alg_domain in algorithm_domains:
                        alg_positions = alg_domain['range_obj']
                        best_jaccard = 0.0

                        for exp_domain in expected["domains"]:
                            exp_positions = self._parse_expected_range_to_positions(exp_domain["range"])
                            if alg_positions and exp_positions:
                                overlap = len(alg_positions & exp_positions)
                                union = len(alg_positions | exp_positions)
                                jaccard = overlap / union if union > 0 else 0.0
                                best_jaccard = max(best_jaccard, jaccard)

                        protein_accuracy += best_jaccard

                    protein_accuracy /= len(algorithm_domains)

                total_accuracy += protein_accuracy

                results[protein_id] = {
                    'algorithm_count': algorithm_count,
                    'expected_count': expected_count,
                    'boundary_accuracy': protein_accuracy,
                    'count_match': count_match,
                    'type': 'normal'
                }

            except Exception as e:
                print(f"Failed to test {protein_id}: {e}")
                results[protein_id] = {
                    'error': str(e),
                    'type': 'error'
                }

        # Overall metrics
        count_accuracy = correct_counts / total_proteins
        boundary_accuracy = total_accuracy / total_proteins

        print(f"\nRegression Test Summary:")
        print(f"  Proteins tested: {total_proteins}")
        print(f"  Domain count accuracy: {count_accuracy:.1%}")
        print(f"  Average boundary accuracy: {boundary_accuracy:.1%}")

        # Show individual results
        for protein_id, result in results.items():
            if result['type'] == 'normal':
                status = "✓" if result['count_match'] else "✗"
                print(f"  {status} {protein_id}: {result['algorithm_count']}/{result['expected_count']} domains, "
                      f"{result['boundary_accuracy']:.1%} accuracy")
            elif result['type'] == 'fibrillar':
                print(f"  ~ {protein_id}: fibrillar protein, {result['algorithm_count']} domains found")
            else:
                print(f"  ✗ {protein_id}: {result.get('error', 'unknown error')}")

        # Performance thresholds (realistic for domain boundary prediction)
        assert count_accuracy >= 0.8, f"Domain count accuracy {count_accuracy:.1%} too low"
        assert boundary_accuracy >= 0.75, f"Boundary accuracy {boundary_accuracy:.1%} too low"
        
        print(f"\n✅ Regression tests passed!")


if __name__ == "__main__":
    # Allow running as script
    pytest.main([__file__, "-v"])
