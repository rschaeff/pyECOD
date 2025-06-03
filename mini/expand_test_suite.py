#!/usr/bin/env python3
"""
Test Suite Expansion for mini_pyecod

This script systematically tests all candidate proteins from batch_test_proteins.py,
captures results, and prepares them for manual curation.

Usage:
    python expand_test_suite.py --run-all         # Test all candidates  
    python expand_test_suite.py --validate-setup  # Check that everything is ready
    python expand_test_suite.py --curate          # Generate curation templates
    python expand_test_suite.py --summary         # Show current status
"""

import os
import sys
import json
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, asdict
import xml.etree.ElementTree as ET
from datetime import datetime

# Add paths for imports
sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.batch_test_proteins import BATCH_TEST_PROTEINS, BATCH_CATEGORIES

@dataclass
class TestResult:
    """Result of running mini_pyecod on a protein"""
    protein_id: str
    success: bool
    domains_found: int
    sequence_length: int
    runtime_seconds: float
    output_file: Optional[str] = None
    error_message: Optional[str] = None
    evidence_count: Optional[int] = None
    domain_details: List[Dict] = None

    def __post_init__(self):
        if self.domain_details is None:
            self.domain_details = []

@dataclass  
class DomainBoundary:
    """Domain boundary information for curation"""
    id: str
    range_str: str
    family: str
    source: str
    size: int
    start: int
    end: int
    is_discontinuous: bool = False
    segments: List[Tuple[int, int]] = None

    def __post_init__(self):
        if self.segments is None:
            self.segments = []

class TestSuiteExpander:
    """Manages expansion and testing of the protein test suite"""
    
    def __init__(self, output_dir: str = "/tmp/test_suite_expansion"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Results storage
        self.results_file = self.output_dir / "test_results.json"
        self.curation_dir = self.output_dir / "curation"
        self.curation_dir.mkdir(exist_ok=True)
        
        # Mini pyecod executable
        self.mini_executable = Path(__file__).parent / "pyecod_mini"
        
        # Test results
        self.test_results: Dict[str, TestResult] = {}
        
    def validate_setup(self) -> Tuple[bool, List[str]]:
        """Validate that everything is ready for testing"""
        issues = []
        
        print("🔍 Validating test suite expansion setup...")
        
        # Check mini_pyecod executable
        if not self.mini_executable.exists():
            issues.append(f"pyecod_mini executable not found: {self.mini_executable}")
        else:
            print(f"✓ Found pyecod_mini: {self.mini_executable}")
            
        # Check that pyecod_mini works
        try:
            result = subprocess.run(
                [str(self.mini_executable), "--validate"],
                capture_output=True,
                text=True,
                timeout=30
            )
            if result.returncode != 0:
                issues.append(f"pyecod_mini validation failed: {result.stderr}")
            else:
                print("✓ pyecod_mini validation passed")
        except Exception as e:
            issues.append(f"Failed to run pyecod_mini: {e}")
            
        # Check test proteins
        print(f"✓ Found {len(BATCH_TEST_PROTEINS)} test protein candidates")
        for category, proteins in BATCH_CATEGORIES.items():
            print(f"  {category}: {len(proteins)} proteins")
            
        # Check batch directory availability
        try:
            result = subprocess.run(
                [str(self.mini_executable), "--list-batches"],
                capture_output=True,
                text=True,
                timeout=10
            )
            if result.returncode == 0:
                print("✓ Batch directories accessible")
            else:
                issues.append("Cannot access batch directories")
        except Exception as e:
            issues.append(f"Failed to check batch directories: {e}")
            
        return len(issues) == 0, issues
    
    def run_single_test(self, protein_id: str, verbose: bool = False) -> TestResult:
        """Run mini_pyecod on a single protein and capture results"""
        
        print(f"Testing {protein_id}...")
        
        start_time = datetime.now()
        
        try:
            # Run mini_pyecod
            cmd = [str(self.mini_executable), protein_id]
            if verbose:
                cmd.append("--verbose")
                
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=120  # 2 minute timeout
            )
            
            end_time = datetime.now()
            runtime = (end_time - start_time).total_seconds()
            
            if result.returncode == 0:
                # Success - parse output file
                output_file = f"/tmp/{protein_id}_mini.domains.xml"
                
                domains = []
                sequence_length = 0
                evidence_count = None
                
                if Path(output_file).exists():
                    domains, sequence_length = self._parse_domain_output(output_file)
                
                # Try to extract evidence count from stdout
                for line in result.stdout.split('\n'):
                    if "Found" in line and "evidence items" in line:
                        try:
                            evidence_count = int(line.split()[1])
                        except:
                            pass
                
                test_result = TestResult(
                    protein_id=protein_id,
                    success=True,
                    domains_found=len(domains),
                    sequence_length=sequence_length,
                    runtime_seconds=runtime,
                    output_file=output_file,
                    evidence_count=evidence_count,
                    domain_details=domains
                )
                
                print(f"  ✓ Success: {len(domains)} domains, {sequence_length} residues")
                
            else:
                # Failure
                test_result = TestResult(
                    protein_id=protein_id,
                    success=False,
                    domains_found=0,
                    sequence_length=0,
                    runtime_seconds=runtime,
                    error_message=result.stderr or result.stdout
                )
                
                print(f"  ✗ Failed: {result.stderr[:100]}...")
                
        except subprocess.TimeoutExpired:
            test_result = TestResult(
                protein_id=protein_id,
                success=False,
                domains_found=0,
                sequence_length=0,
                runtime_seconds=120.0,
                error_message="Timeout after 2 minutes"
            )
            print(f"  ✗ Timeout")
            
        except Exception as e:
            test_result = TestResult(
                protein_id=protein_id,
                success=False,
                domains_found=0,
                sequence_length=0,
                runtime_seconds=0.0,
                error_message=str(e)
            )
            print(f"  ✗ Error: {e}")
            
        return test_result
    
    def _parse_domain_output(self, xml_file: str) -> Tuple[List[Dict], int]:
        """Parse domain XML output to extract domain details"""
        domains = []
        sequence_length = 0
        
        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()
            
            for domain_elem in root.findall(".//domain"):
                domain_id = domain_elem.get("id", "")
                range_str = domain_elem.get("range", "")
                family = domain_elem.get("family", "")
                source = domain_elem.get("source", "")
                size = int(domain_elem.get("size", "0"))
                is_disc = domain_elem.get("is_discontinuous", "false").lower() == "true"
                
                # Parse range to get start/end
                start, end = 0, 0
                segments = []
                
                if range_str:
                    try:
                        for segment in range_str.split(','):
                            if '-' in segment:
                                s, e = map(int, segment.split('-'))
                                segments.append((s, e))
                                if not start or s < start:
                                    start = s
                                if e > end:
                                    end = e
                    except:
                        pass
                
                domains.append({
                    "id": domain_id,
                    "range": range_str,
                    "family": family,
                    "source": source,
                    "size": size,
                    "start": start,
                    "end": end,
                    "is_discontinuous": is_disc,
                    "segments": segments
                })
                
                # Track max position for sequence length estimate
                if end > sequence_length:
                    sequence_length = end
                    
        except Exception as e:
            print(f"Warning: Could not parse {xml_file}: {e}")
            
        return domains, sequence_length
    
    def run_all_tests(self, verbose: bool = False, 
                     categories: List[str] = None,
                     limit_per_category: int = None) -> Dict[str, TestResult]:
        """Run tests on all candidate proteins"""
        
        print("🚀 Running mini_pyecod on all test candidates...")
        print("=" * 60)
        
        # Determine which proteins to test
        proteins_to_test = []
        
        if categories:
            for category in categories:
                if category in BATCH_CATEGORIES:
                    cat_proteins = BATCH_CATEGORIES[category]
                    if limit_per_category:
                        cat_proteins = cat_proteins[:limit_per_category]
                    proteins_to_test.extend(cat_proteins)
        else:
            # Test all proteins, with optional limits per category
            for category, cat_proteins in BATCH_CATEGORIES.items():
                if limit_per_category:
                    proteins_to_test.extend(cat_proteins[:limit_per_category])
                else:
                    proteins_to_test.extend(cat_proteins)
        
        # Remove duplicates while preserving order
        proteins_to_test = list(dict.fromkeys(proteins_to_test))
        
        print(f"Testing {len(proteins_to_test)} proteins...")
        
        results = {}
        successful = 0
        failed = 0
        
        for i, protein_id in enumerate(proteins_to_test, 1):
            print(f"\n[{i}/{len(proteins_to_test)}] ", end="")
            
            result = self.run_single_test(protein_id, verbose)
            results[protein_id] = result
            
            if result.success:
                successful += 1
            else:
                failed += 1
        
        # Save results
        self.test_results.update(results)
        self._save_results()
        
        # Summary
        print(f"\n{'='*60}")
        print(f"TEST SUMMARY")
        print(f"{'='*60}")
        print(f"Total tested: {len(proteins_to_test)}")
        print(f"Successful: {successful}")
        print(f"Failed: {failed}")
        print(f"Success rate: {successful/len(proteins_to_test)*100:.1f}%")
        
        return results
    
    def _save_results(self):
        """Save test results to JSON file"""
        # Convert to JSON-serializable format
        serializable = {}
        for protein_id, result in self.test_results.items():
            serializable[protein_id] = asdict(result)
            
        with open(self.results_file, 'w') as f:
            json.dump(serializable, f, indent=2)
            
        print(f"\n💾 Results saved to: {self.results_file}")
    
    def load_results(self) -> bool:
        """Load existing test results"""
        if not self.results_file.exists():
            return False
            
        try:
            with open(self.results_file, 'r') as f:
                data = json.load(f)
                
            self.test_results = {}
            for protein_id, result_data in data.items():
                self.test_results[protein_id] = TestResult(**result_data)
                
            print(f"📖 Loaded {len(self.test_results)} previous results")
            return True
            
        except Exception as e:
            print(f"Warning: Could not load results: {e}")
            return False
    
    def show_summary(self):
        """Show summary of current test results"""
        if not self.test_results:
            if not self.load_results():
                print("No test results available. Run --run-all first.")
                return
        
        print("📊 TEST SUITE SUMMARY")
        print("=" * 60)
        
        # Overall stats
        total = len(self.test_results)
        successful = sum(1 for r in self.test_results.values() if r.success)
        failed = total - successful
        
        print(f"Total proteins tested: {total}")
        print(f"Successful: {successful} ({successful/total*100:.1f}%)")
        print(f"Failed: {failed} ({failed/total*100:.1f}%)")
        
        # By category
        print(f"\nResults by category:")
        for category, proteins in BATCH_CATEGORIES.items():
            tested = [p for p in proteins if p in self.test_results]
            if tested:
                success_count = sum(1 for p in tested if self.test_results[p].success)
                print(f"  {category}: {success_count}/{len(tested)} successful")
        
        # Successful proteins details
        print(f"\nSuccessful proteins:")
        for protein_id, result in self.test_results.items():
            if result.success:
                info = BATCH_TEST_PROTEINS.get(protein_id, {})
                category = info.get('category', 'unknown')
                print(f"  {protein_id}: {result.domains_found} domains, "
                      f"{result.sequence_length} residues ({category})")
        
        # Failed proteins
        failures = [p for p, r in self.test_results.items() if not r.success]
        if failures:
            print(f"\nFailed proteins:")
            for protein_id in failures:
                result = self.test_results[protein_id]
                error = result.error_message[:50] if result.error_message else "Unknown error"
                print(f"  {protein_id}: {error}...")
    
    def generate_curation_templates(self):
        """Generate templates for manual curation of successful proteins"""
        
        if not self.test_results:
            if not self.load_results():
                print("No test results available. Run --run-all first.")
                return
        
        print("📝 Generating curation templates...")
        
        successful_proteins = [p for p, r in self.test_results.items() if r.success]
        
        if not successful_proteins:
            print("No successful proteins to curate!")
            return
        
        # Create individual curation files
        for protein_id in successful_proteins:
            self._create_curation_template(protein_id)
        
        # Create master curation file
        self._create_master_curation_file(successful_proteins)
        
        print(f"✅ Created curation templates in: {self.curation_dir}")
        print(f"📋 Master file: {self.curation_dir / 'master_curation.md'}")
    
    def _create_curation_template(self, protein_id: str):
        """Create individual curation template for a protein"""
        
        result = self.test_results[protein_id]
        info = BATCH_TEST_PROTEINS.get(protein_id, {})
        
        template_file = self.curation_dir / f"{protein_id}_curation.md"
        
        lines = [
            f"# Domain Curation for {protein_id}",
            f"",
            f"## Protein Information",
            f"- **Category**: {info.get('category', 'unknown')}",
            f"- **Description**: {info.get('description', '')}",
            f"- **Estimated Length**: {info.get('estimated_length', 'unknown')} residues",
            f"- **Evidence Count**: {result.evidence_count or 'unknown'}",
            f"",
            f"## Algorithm Results",
            f"- **Domains Found**: {result.domains_found}",
            f"- **Sequence Length**: {result.sequence_length} residues",
            f"- **Runtime**: {result.runtime_seconds:.2f} seconds",
            f"- **Output File**: {result.output_file}",
            f"",
            f"## Domain Boundaries (Algorithm)",
            f""
        ]
        
        if result.domain_details:
            for i, domain in enumerate(result.domain_details, 1):
                lines.extend([
                    f"### Domain {i}: {domain['id']}",
                    f"- **Range**: {domain['range']}",
                    f"- **Family**: {domain['family']}",
                    f"- **Source**: {domain['source']}",
                    f"- **Size**: {domain['size']} residues",
                    f"- **Discontinuous**: {domain['is_discontinuous']}",
                    f""
                ])
        
        lines.extend([
            f"## Manual Curation",
            f"",
            f"### Expected Domains",
            f"<!-- Fill in your manual domain assignments -->",
            f"",
            f"1. **Domain 1**",
            f"   - Range: TODO",
            f"   - Family/Function: TODO",
            f"   - Notes: TODO",
            f"",
            f"2. **Domain 2** (if applicable)",
            f"   - Range: TODO", 
            f"   - Family/Function: TODO",
            f"   - Notes: TODO",
            f"",
            f"### Validation Status",
            f"- [ ] Domains manually curated",
            f"- [ ] Algorithm boundaries acceptable",
            f"- [ ] Algorithm boundaries need adjustment",
            f"- [ ] Test case ready for regression testing",
            f"",
            f"### Notes",
            f"<!-- Add any notes about this protein's domain architecture -->",
            f"",
            f"### PyMOL Visualization Commands",
            f"```pymol",
            f"# Load structure",
            f"fetch {protein_id.split('_')[0]}",
            f"",
            f"# Algorithm domains",
        ])
        
        # Add PyMOL coloring commands for algorithm domains
        colors = ["red", "blue", "green", "yellow", "cyan", "magenta"]
        if result.domain_details:
            for i, domain in enumerate(result.domain_details):
                color = colors[i % len(colors)]
                chain = protein_id.split('_')[1] if '_' in protein_id else 'A'
                lines.append(f"color {color}, chain {chain} and resi {domain['range']}")
        
        lines.extend([
            f"```",
            f"",
            f"---",
            f"*Generated by test suite expansion on {datetime.now().strftime('%Y-%m-%d %H:%M')}*"
        ])
        
        with open(template_file, 'w') as f:
            f.write('\n'.join(lines))
    
    def _create_master_curation_file(self, proteins: List[str]):
        """Create master curation tracking file"""
        
        master_file = self.curation_dir / "master_curation.md"
        
        lines = [
            "# Test Suite Curation Master List",
            "",
            f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
            f"Total proteins: {len(proteins)}",
            "",
            "## Curation Progress",
            ""
        ]
        
        # Group by category
        by_category = {}
        for protein_id in proteins:
            info = BATCH_TEST_PROTEINS.get(protein_id, {})
            category = info.get('category', 'unknown')
            if category not in by_category:
                by_category[category] = []
            by_category[category].append(protein_id)
        
        for category, cat_proteins in by_category.items():
            lines.extend([
                f"### {category.replace('_', ' ').title()} ({len(cat_proteins)} proteins)",
                ""
            ])
            
            for protein_id in cat_proteins:
                result = self.test_results[protein_id]
                lines.append(f"- [ ] **{protein_id}**: {result.domains_found} domains, "
                           f"{result.sequence_length} residues")
            
            lines.append("")
        
        lines.extend([
            "## Curation Workflow",
            "",
            "1. **Review Algorithm Results**: Check each protein's domain assignments",
            "2. **Manual Validation**: Use PyMOL, literature, or domain databases",
            "3. **Curate Boundaries**: Define expected domain boundaries",
            "4. **Mark Complete**: Check off completed proteins above",
            "5. **Generate Tests**: Create regression tests for curated proteins",
            "",
            "## Files",
            ""
        ])
        
        for protein_id in proteins:
            lines.append(f"- [{protein_id}_curation.md](./{protein_id}_curation.md)")
        
        with open(master_file, 'w') as f:
            f.write('\n'.join(lines))

def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Expand mini_pyecod test suite',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python expand_test_suite.py --validate-setup     # Check readiness
    python expand_test_suite.py --run-all            # Test all proteins
    python expand_test_suite.py --run-all --limit 2  # Test 2 per category
    python expand_test_suite.py --categories validated single_domain_small
    python expand_test_suite.py --summary            # Show current results
    python expand_test_suite.py --curate             # Generate curation templates
        """
    )
    
    parser.add_argument('--validate-setup', action='store_true',
                       help='Validate that everything is ready for testing')
    parser.add_argument('--run-all', action='store_true',
                       help='Run tests on all candidate proteins')
    parser.add_argument('--categories', nargs='+',
                       choices=list(BATCH_CATEGORIES.keys()),
                       help='Test only specific categories')
    parser.add_argument('--limit', type=int,
                       help='Limit number of proteins per category')
    parser.add_argument('--summary', action='store_true',
                       help='Show summary of current test results')
    parser.add_argument('--curate', action='store_true',
                       help='Generate curation templates for successful proteins')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Verbose output during testing')
    parser.add_argument('--output-dir', default='/tmp/test_suite_expansion',
                       help='Output directory for results and curation files')
    
    args = parser.parse_args()
    
    expander = TestSuiteExpander(args.output_dir)
    
    if args.validate_setup:
        is_valid, issues = expander.validate_setup()
        if is_valid:
            print("\n✅ Setup validation passed - ready to run tests!")
        else:
            print("\n❌ Setup validation failed:")
            for issue in issues:
                print(f"  - {issue}")
            sys.exit(1)
        return
    
    if args.summary:
        expander.show_summary()
        return
        
    if args.curate:
        expander.generate_curation_templates()
        return
    
    if args.run_all:
        results = expander.run_all_tests(
            verbose=args.verbose,
            categories=args.categories,
            limit_per_category=args.limit
        )
        
        # Suggest next steps
        successful = sum(1 for r in results.values() if r.success)
        if successful > 0:
            print(f"\n🎯 Next Steps:")
            print(f"1. Review results: python expand_test_suite.py --summary")
            print(f"2. Generate curation templates: python expand_test_suite.py --curate")
            print(f"3. Manually curate domain boundaries in {expander.curation_dir}")
            print(f"4. Create regression tests for curated proteins")
        
        return
    
    # No action specified
    parser.print_help()
    print(f"\n💡 Tip: Start with --validate-setup, then --run-all")

if __name__ == "__main__":
    main()
