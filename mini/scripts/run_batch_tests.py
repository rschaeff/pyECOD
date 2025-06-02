#!/usr/bin/env python3
"""Batch test runner for mini_pyecod with result analysis"""

import sys
import os
from pathlib import Path
import json
from datetime import datetime

# Add parent directory to path so we can import ecod
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.parser import parse_domain_summary
from mini.partitioner import partition_domains
from mini.writer import write_domain_partition

def analyze_test_results(protein_id, domains, evidence, sequence_length):
    """Analyze test results and generate metrics"""
    
    # Basic metrics
    metrics = {
        'protein_id': protein_id,
        'domain_count': len(domains),
        'sequence_length': sequence_length,
        'evidence_count': len(evidence),
        'coverage': 0.0,
        'domains': []
    }
    
    # Calculate coverage
    covered_positions = set()
    for domain in domains:
        positions = set(domain.range.get_positions())
        covered_positions.update(positions)
        
        # Domain details
        domain_info = {
            'id': domain.id,
            'family': domain.family,
            'range': str(domain.range),
            'size': domain.range.size,
            'source': domain.source,
            'discontinuous': domain.range.is_discontinuous
        }
        metrics['domains'].append(domain_info)
    
    metrics['coverage'] = len(covered_positions) / sequence_length if sequence_length > 0 else 0
    
    # Evidence breakdown
    evidence_types = {}
    for ev in evidence:
        evidence_types[ev.type] = evidence_types.get(ev.type, 0) + 1
    metrics['evidence_types'] = evidence_types
    
    return metrics

def run_single_test(protein_id, batch_dir, verbose=False):
    """Run test for a single protein and return results"""
    
    # Parse PDB and chain from protein_id
    parts = protein_id.split('_')
    if len(parts) >= 2:
        pdb_id = parts[0]
        chain_id = parts[1]
    else:
        return None, f"Invalid protein ID format: {protein_id}"
    
    # Construct file path
    xml_path = os.path.join(batch_dir, "domains", f"{protein_id}.develop291.domain_summary.xml")
    
    if not os.path.exists(xml_path):
        return None, f"Domain summary file not found: {xml_path}"
    
    try:
        # Parse evidence
        evidence = parse_domain_summary(xml_path, verbose=verbose)
        
        # Get sequence length from evidence
        max_pos = 0
        for ev in evidence:
            for start, end in ev.query_range.segments:
                max_pos = max(max_pos, end)
        sequence_length = int(max_pos * 1.1)
        
        # Filter to high-quality evidence
        good_evidence = [e for e in evidence if e.confidence > 0.5 or (e.evalue and e.evalue < 1e-3)]
        
        # Partition domains (suppress verbose output)
        domains = partition_domains(good_evidence, sequence_length=sequence_length, verbose=False)
        
        # Analyze results
        metrics = analyze_test_results(protein_id, domains, good_evidence, sequence_length)
        
        # Write output
        output_path = f"/tmp/batch_test/{protein_id}_mini.domains.xml"
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        write_domain_partition(domains, pdb_id, chain_id, output_path)
        metrics['output_file'] = output_path
        
        return metrics, None
        
    except Exception as e:
        import traceback
        error_msg = f"Error processing {protein_id}: {str(e)}\n{traceback.format_exc()}"
        return None, error_msg

def run_batch_tests(test_proteins, batch_dir, verbose=False):
    """Run tests on multiple proteins and generate report"""
    
    print(f"\n{'='*60}")
    print(f"Running batch tests on {len(test_proteins)} proteins")
    print(f"{'='*60}\n")
    
    results = []
    errors = []
    
    for i, protein_id in enumerate(test_proteins, 1):
        print(f"[{i}/{len(test_proteins)}] Testing {protein_id}...", end='', flush=True)
        
        metrics, error = run_single_test(protein_id, batch_dir, verbose=verbose)
        
        if metrics:
            results.append(metrics)
            print(f" ✓ {metrics['domain_count']} domains found")
        else:
            errors.append({'protein_id': protein_id, 'error': error})
            print(f" ✗ Failed")
    
    return results, errors

def generate_report(results, errors, output_dir):
    """Generate comprehensive test report"""
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_file = os.path.join(output_dir, f"test_report_{timestamp}.txt")
    json_file = os.path.join(output_dir, f"test_results_{timestamp}.json")
    
    # Save raw results as JSON
    with open(json_file, 'w') as f:
        json.dump({'results': results, 'errors': errors}, f, indent=2)
    
    # Generate human-readable report
    with open(report_file, 'w') as f:
        f.write("="*60 + "\n")
        f.write("pyECOD Mini Test Report\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("="*60 + "\n\n")
        
        # Summary statistics
        f.write("SUMMARY\n")
        f.write("-"*30 + "\n")
        f.write(f"Total tests run: {len(results) + len(errors)}\n")
        f.write(f"Successful: {len(results)}\n")
        f.write(f"Failed: {len(errors)}\n")
        
        if results:
            # Domain count distribution
            domain_counts = {}
            for r in results:
                count = r['domain_count']
                domain_counts[count] = domain_counts.get(count, 0) + 1
            
            f.write("\nDomain count distribution:\n")
            for count in sorted(domain_counts.keys()):
                f.write(f"  {count} domains: {domain_counts[count]} proteins\n")
            
            # Coverage statistics
            coverages = [r['coverage'] for r in results]
            avg_coverage = sum(coverages) / len(coverages)
            f.write(f"\nAverage sequence coverage: {avg_coverage:.1%}\n")
        
        # Detailed results
        f.write("\n" + "="*60 + "\n")
        f.write("DETAILED RESULTS\n")
        f.write("="*60 + "\n")
        
        # Group by domain count
        by_domain_count = {}
        for r in results:
            count = r['domain_count']
            if count not in by_domain_count:
                by_domain_count[count] = []
            by_domain_count[count].append(r)
        
        for count in sorted(by_domain_count.keys()):
            f.write(f"\n{count}-DOMAIN PROTEINS ({len(by_domain_count[count])} proteins)\n")
            f.write("-"*40 + "\n")
            
            for r in by_domain_count[count]:
                f.write(f"\n{r['protein_id']}:\n")
                f.write(f"  Sequence length: {r['sequence_length']}\n")
                f.write(f"  Coverage: {r['coverage']:.1%}\n")
                f.write(f"  Evidence: {r['evidence_count']} items")
                
                # Evidence breakdown
                if r['evidence_types']:
                    types_str = ", ".join([f"{k}:{v}" for k, v in r['evidence_types'].items()])
                    f.write(f" ({types_str})")
                f.write("\n")
                
                # Domain details
                for d in r['domains']:
                    f.write(f"  Domain {d['id']}: {d['family']} @ {d['range']} ")
                    f.write(f"({d['size']} residues, source: {d['source']})\n")
        
        # Errors
        if errors:
            f.write("\n" + "="*60 + "\n")
            f.write("ERRORS\n")
            f.write("="*60 + "\n")
            
            for e in errors:
                f.write(f"\n{e['protein_id']}:\n")
                f.write(f"  {e['error']}\n")
    
    print(f"\n{'='*60}")
    print(f"Report generated:")
    print(f"  Text report: {report_file}")
    print(f"  JSON results: {json_file}")
    print(f"{'='*60}")

def main():
    """Main entry point"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Run batch tests for mini_pyecod')
    parser.add_argument('--batch-dir', default="/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424",
                        help='Batch directory containing domain files')
    parser.add_argument('--output-dir', default="/tmp/batch_test",
                        help='Output directory for results and reports')
    parser.add_argument('--test-list', help='File containing list of protein IDs to test')
    parser.add_argument('--auto-find', type=int, metavar='N',
                        help='Automatically find and test N proteins')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose output')
    parser.add_argument('proteins', nargs='*', help='Protein IDs to test')
    
    args = parser.parse_args()
    
    # Determine which proteins to test
    test_proteins = []
    
    if args.test_list:
        # Read from file
        with open(args.test_list, 'r') as f:
            test_proteins = [line.strip() for line in f if line.strip()]
    elif args.auto_find:
        # Auto-find test cases
        print("Auto-finding test cases...")
        from mini.find_real_test_cases import analyze_existing_files, verify_test_suite
        test_suite = analyze_existing_files(args.batch_dir, sample_size=args.auto_find * 2)
        verified = verify_test_suite(args.batch_dir, test_suite)
        test_proteins = [t['protein'] for t in verified[:args.auto_find]]
    elif args.proteins:
        # Use command line arguments
        test_proteins = args.proteins
    else:
        # Default test set
        test_proteins = [
            '8ovp_A',    # Known 3-domain case
            '8opd_Aa',   # 2-domain case
            '9nlb_E',    # From batch
        ]
    
    if not test_proteins:
        print("No proteins to test!")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Run tests
    results, errors = run_batch_tests(test_proteins, args.batch_dir, verbose=args.verbose)
    
    # Generate report
    generate_report(results, errors, args.output_dir)

if __name__ == "__main__":
    main()
