#!/usr/bin/env python3
"""
GFP Beta Barrel Domain Assignment Analysis

Find and analyze GFP proteins (ECOD T_id 271.1.1) in curation sets
to debug why they're not being assigned as single domains.
"""

import os
import sys
import json
import logging
from typing import Dict, List, Any, Optional

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext
from ecod.evaluation.algorithm_versions.manager import AlgorithmVersionManager

class GFPAnalyzer:
    """Analyze GFP domain assignment issues"""
    
    def __init__(self, config_path: str = "config/config.yml"):
        self.context = ApplicationContext(config_path)
        self.db = self.context.db
        self.logger = logging.getLogger(__name__)
        
    def find_gfp_proteins_in_curation(self) -> List[Dict[str, Any]]:
        """Find GFP proteins (T_id 271.1.1) in curation sets"""
        
        query = """
        SELECT DISTINCT 
            p.pdb_id, 
            p.chain_id,
            p.id as protein_id,
            cd.has_domain,
            cd.domain_assigned_correctly,
            cd.boundaries_correct,
            cd.confidence_level,
            cd.curator_name,
            cs.session_id
        FROM pdb_analysis.protein p
        JOIN pdb_analysis.curation_decision cd ON p.id = cd.protein_id
        JOIN pdb_analysis.curation_session cs ON cd.session_id = cs.id
        WHERE EXISTS (
            -- Find proteins with ECOD assignments to GFP topology
            SELECT 1 FROM pdb_analysis.partition_domains pd
            JOIN pdb_analysis.partition_proteins pp ON pd.protein_id = pp.id
            WHERE pp.pdb_id = p.pdb_id AND pp.chain_id = p.chain_id
              AND pd.t_group = '271.1.1'
        )
        AND cs.status = 'committed'
        ORDER BY p.pdb_id, p.chain_id
        """
        
        results = self.db.execute_dict_query(query)
        
        self.logger.info(f"Found {len(results)} GFP proteins in curation sets")
        
        return results
    
    def analyze_gfp_domain_assignment(self, pdb_id: str, chain_id: str) -> Dict[str, Any]:
        """Analyze how a specific GFP protein is currently assigned domains"""
        
        # Get current partition results
        partition_query = """
        SELECT 
            pp.id as partition_protein_id,
            pp.is_classified,
            pp.coverage,
            pp.sequence_length,
            pp.residues_assigned,
            pp.domains_with_evidence,
            pp.fully_classified_domains,
            pp.timestamp
        FROM pdb_analysis.partition_proteins pp
        JOIN pdb_analysis.protein p ON pp.pdb_id = p.pdb_id AND pp.chain_id = p.chain_id
        WHERE p.pdb_id = %s AND p.chain_id = %s
        ORDER BY pp.timestamp DESC
        LIMIT 1
        """
        
        partition_results = self.db.execute_dict_query(partition_query, (pdb_id, chain_id))
        
        if not partition_results:
            return {
                'error': f'No partition results found for {pdb_id}_{chain_id}',
                'protein': f'{pdb_id}_{chain_id}'
            }
        
        partition_info = partition_results[0]
        partition_protein_id = partition_info['partition_protein_id']
        
        # Get domain assignments
        domains_query = """
        SELECT 
            domain_number,
            start_pos,
            end_pos,
            range,
            source,
            source_id,
            confidence,
            t_group,
            h_group,
            x_group,
            a_group,
            pdb_range,
            length
        FROM pdb_analysis.partition_domains
        WHERE protein_id = %s
        ORDER BY domain_number
        """
        
        domains = self.db.execute_dict_query(domains_query, (partition_protein_id,))
        
        # Analyze the assignment
        analysis = {
            'protein': f'{pdb_id}_{chain_id}',
            'partition_info': partition_info,
            'domain_count': len(domains),
            'domains': domains,
            'is_single_domain': len(domains) == 1,
            'gfp_domains': [d for d in domains if d.get('t_group') == '271.1.1'],
            'non_gfp_domains': [d for d in domains if d.get('t_group') != '271.1.1'],
            'issues': []
        }
        
        # Identify issues
        if len(domains) == 0:
            analysis['issues'].append('Not classified as having any domains')
        elif len(domains) > 1:
            analysis['issues'].append(f'Over-segmented into {len(domains)} domains (should be 1)')
            
            # Check if multiple domains have GFP assignment
            gfp_domain_count = len(analysis['gfp_domains'])
            if gfp_domain_count > 1:
                analysis['issues'].append(f'Multiple domains ({gfp_domain_count}) assigned to GFP topology')
            elif gfp_domain_count == 0:
                analysis['issues'].append('No domains assigned to GFP topology despite known GFP structure')
            
        elif len(domains) == 1:
            domain = domains[0]
            if domain.get('t_group') != '271.1.1':
                analysis['issues'].append(f'Single domain but wrong topology: {domain.get("t_group")} (should be 271.1.1)')
            else:
                analysis['issues'].append('Correctly assigned as single GFP domain')
        
        return analysis
    
    def find_domain_summaries_for_gfp(self, pdb_id: str, chain_id: str) -> List[str]:
        """Find domain summary files for GFP protein"""
        
        query = """
        SELECT DISTINCT pf.file_path, b.base_path, ps.batch_id
        FROM ecod_schema.process_file pf
        JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        JOIN ecod_schema.batch b ON ps.batch_id = b.id
        WHERE p.pdb_id = %s AND p.chain_id = %s
          AND pf.file_type IN ('domain_summary', 'blast_only_summary')
          AND pf.file_exists = TRUE
        ORDER BY ps.batch_id DESC
        LIMIT 5
        """
        
        results = self.db.execute_dict_query(query, (pdb_id, chain_id))
        
        summary_paths = []
        for row in results:
            full_path = os.path.join(row['base_path'], row['file_path'])
            if os.path.exists(full_path):
                summary_paths.append({
                    'path': full_path,
                    'batch_id': row['batch_id']
                })
        
        return summary_paths
    
    def analyze_gfp_evidence(self, pdb_id: str, chain_id: str) -> Dict[str, Any]:
        """Analyze evidence types available for GFP protein"""
        
        summary_files = self.find_domain_summaries_for_gfp(pdb_id, chain_id)
        
        if not summary_files:
            return {
                'error': f'No domain summary files found for {pdb_id}_{chain_id}',
                'protein': f'{pdb_id}_{chain_id}'
            }
        
        # Analyze the most recent summary file
        latest_summary = summary_files[0]['path']
        
        try:
            with open(latest_summary, 'r') as f:
                content = f.read()
                
            # Count evidence types (this would need to match your actual format)
            evidence_counts = {
                'hhsearch': content.lower().count('hhsearch'),
                'chain_blast': content.lower().count('chain_blast') + content.lower().count('chain blast'),
                'domain_blast': content.lower().count('domain_blast') + content.lower().count('domain blast'),
                'blast': content.lower().count('blast') - content.lower().count('chain_blast') - content.lower().count('domain_blast'),
                'total_lines': len(content.split('\n'))
            }
            
            # Look for GFP-specific evidence
            gfp_evidence = {
                'gfp_mentions': content.lower().count('fluorescent') + content.lower().count('gfp'),
                't_group_271': content.count('271.1.1'),
                'beta_barrel_mentions': content.lower().count('beta barrel') + content.lower().count('barrel')
            }
            
            return {
                'protein': f'{pdb_id}_{chain_id}',
                'summary_file': latest_summary,
                'batch_id': summary_files[0]['batch_id'],
                'evidence_counts': evidence_counts,
                'gfp_evidence': gfp_evidence,
                'file_size': len(content),
                'available_summaries': len(summary_files)
            }
            
        except Exception as e:
            return {
                'error': f'Failed to read summary file: {e}',
                'protein': f'{pdb_id}_{chain_id}',
                'summary_file': latest_summary
            }
    
    def comprehensive_gfp_analysis(self) -> Dict[str, Any]:
        """Comprehensive analysis of all GFP proteins in curation"""
        
        gfp_proteins = self.find_gfp_proteins_in_curation()
        
        analysis = {
            'total_gfp_proteins': len(gfp_proteins),
            'proteins_analyzed': [],
            'assignment_summary': {
                'correctly_single_domain': 0,
                'over_segmented': 0,
                'not_classified': 0,
                'wrong_topology': 0
            },
            'curation_summary': {
                'has_domain_true': 0,
                'domain_assigned_correctly_true': 0,
                'boundaries_correct_true': 0,
                'high_confidence': 0  # confidence >= 4
            },
            'evidence_summary': {
                'proteins_with_summaries': 0,
                'avg_evidence_counts': {},
                'gfp_specific_evidence': 0
            },
            'problematic_cases': []
        }
        
        for protein_info in gfp_proteins:
            pdb_id = protein_info['pdb_id']
            chain_id = protein_info['chain_id']
            
            self.logger.info(f"Analyzing {pdb_id}_{chain_id}")
            
            # Analyze domain assignment
            domain_analysis = self.analyze_gfp_domain_assignment(pdb_id, chain_id)
            
            # Analyze evidence
            evidence_analysis = self.analyze_gfp_evidence(pdb_id, chain_id)
            
            protein_analysis = {
                'protein': f'{pdb_id}_{chain_id}',
                'curation_info': protein_info,
                'domain_analysis': domain_analysis,
                'evidence_analysis': evidence_analysis
            }
            
            analysis['proteins_analyzed'].append(protein_analysis)
            
            # Update summaries
            if protein_info['has_domain']:
                analysis['curation_summary']['has_domain_true'] += 1
            if protein_info['domain_assigned_correctly']:
                analysis['curation_summary']['domain_assigned_correctly_true'] += 1
            if protein_info['boundaries_correct']:
                analysis['curation_summary']['boundaries_correct_true'] += 1
            if protein_info['confidence_level'] >= 4:
                analysis['curation_summary']['high_confidence'] += 1
                
            # Update assignment summary
            if 'error' not in domain_analysis:
                if domain_analysis['is_single_domain'] and len(domain_analysis['gfp_domains']) == 1:
                    analysis['assignment_summary']['correctly_single_domain'] += 1
                elif domain_analysis['domain_count'] > 1:
                    analysis['assignment_summary']['over_segmented'] += 1
                elif domain_analysis['domain_count'] == 0:
                    analysis['assignment_summary']['not_classified'] += 1
                elif domain_analysis['is_single_domain'] and len(domain_analysis['gfp_domains']) == 0:
                    analysis['assignment_summary']['wrong_topology'] += 1
                    
            # Update evidence summary
            if 'error' not in evidence_analysis:
                analysis['evidence_summary']['proteins_with_summaries'] += 1
                
            # Identify problematic cases
            if domain_analysis.get('issues') and any('over-segmented' in issue.lower() or 'wrong topology' in issue.lower() for issue in domain_analysis['issues']):
                analysis['problematic_cases'].append({
                    'protein': f'{pdb_id}_{chain_id}',
                    'issues': domain_analysis['issues'],
                    'domain_count': domain_analysis['domain_count'],
                    'confidence': protein_info['confidence_level']
                })
        
        return analysis


def main():
    """Main CLI interface"""
    import argparse
    
    parser = argparse.ArgumentParser(description="GFP Beta Barrel Domain Assignment Analysis")
    parser.add_argument('--config', type=str, default='config/config.yml', help='Configuration file')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Find GFP proteins
    find_parser = subparsers.add_parser('find', help='Find GFP proteins in curation sets')
    
    # Analyze specific protein
    analyze_parser = subparsers.add_parser('analyze', help='Analyze specific GFP protein')
    analyze_parser.add_argument('--protein', type=str, required=True, 
                               help='Protein to analyze (format: PDB_CHAIN)')
    
    # Comprehensive analysis
    comprehensive_parser = subparsers.add_parser('comprehensive', help='Comprehensive GFP analysis')
    
    # List problematic cases
    problems_parser = subparsers.add_parser('problems', help='List problematic GFP cases')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return 1
    
    # Setup logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=level, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    try:
        analyzer = GFPAnalyzer(args.config)
        
        if args.command == 'find':
            gfp_proteins = analyzer.find_gfp_proteins_in_curation()
            
            print(f"🧬 Found {len(gfp_proteins)} GFP proteins (T_id 271.1.1) in curation sets:")
            print("-" * 70)
            print(f"{'PDB_Chain':15} {'Has Domain':12} {'Correct':10} {'Boundaries':12} {'Confidence':10} {'Curator':15}")
            print("-" * 70)
            
            for protein in gfp_proteins:
                pdb_chain = f"{protein['pdb_id']}_{protein['chain_id']}"
                has_domain = "✅" if protein['has_domain'] else "❌"
                correct = "✅" if protein['domain_assigned_correctly'] else "❌" if protein['domain_assigned_correctly'] is not None else "❓"
                boundaries = "✅" if protein['boundaries_correct'] else "❌" if protein['boundaries_correct'] is not None else "❓"
                confidence = protein['confidence_level']
                curator = protein['curator_name'][:14]
                
                print(f"{pdb_chain:15} {has_domain:12} {correct:10} {boundaries:12} {confidence:10} {curator:15}")
                
        elif args.command == 'analyze':
            pdb_id, chain_id = args.protein.split('_')
            
            print(f"🔍 Analyzing GFP protein: {pdb_id}_{chain_id}")
            print("=" * 50)
            
            domain_analysis = analyzer.analyze_gfp_domain_assignment(pdb_id, chain_id)
            evidence_analysis = analyzer.analyze_gfp_evidence(pdb_id, chain_id)
            
            # Domain assignment analysis
            print("\nDOMAIN ASSIGNMENT:")
            if 'error' in domain_analysis:
                print(f"  ❌ {domain_analysis['error']}")
            else:
                print(f"  Domain count: {domain_analysis['domain_count']}")
                print(f"  Is single domain: {'✅' if domain_analysis['is_single_domain'] else '❌'}")
                print(f"  GFP domains: {len(domain_analysis['gfp_domains'])}")
                print(f"  Non-GFP domains: {len(domain_analysis['non_gfp_domains'])}")
                
                if domain_analysis['issues']:
                    print(f"  Issues:")
                    for issue in domain_analysis['issues']:
                        status = "✅" if "correctly assigned" in issue.lower() else "⚠️"
                        print(f"    {status} {issue}")
                
                if domain_analysis['domains']:
                    print(f"  Domain details:")
                    for i, domain in enumerate(domain_analysis['domains'], 1):
                        print(f"    Domain {i}: {domain['range']} (T:{domain.get('t_group', 'N/A')})")
            
            # Evidence analysis
            print(f"\nEVIDENCE ANALYSIS:")
            if 'error' in evidence_analysis:
                print(f"  ❌ {evidence_analysis['error']}")
            else:
                print(f"  Summary file: {evidence_analysis['summary_file']}")
                print(f"  Batch ID: {evidence_analysis['batch_id']}")
                print(f"  Evidence counts:")
                for evidence_type, count in evidence_analysis['evidence_counts'].items():
                    print(f"    {evidence_type}: {count}")
                print(f"  GFP-specific evidence:")
                for evidence_type, count in evidence_analysis['gfp_evidence'].items():
                    print(f"    {evidence_type}: {count}")
                    
        elif args.command == 'comprehensive':
            print("🧬 Comprehensive GFP Analysis")
            print("=" * 50)
            
            analysis = analyzer.comprehensive_gfp_analysis()
            
            print(f"\nOVERALL SUMMARY:")
            print(f"  Total GFP proteins in curation: {analysis['total_gfp_proteins']}")
            
            print(f"\nCURATION SUMMARY:")
            curation = analysis['curation_summary']
            print(f"  Marked as having domain: {curation['has_domain_true']}/{analysis['total_gfp_proteins']}")
            print(f"  Domain assigned correctly: {curation['domain_assigned_correctly_true']}/{analysis['total_gfp_proteins']}")
            print(f"  Boundaries correct: {curation['boundaries_correct_true']}/{analysis['total_gfp_proteins']}")
            print(f"  High confidence (≥4): {curation['high_confidence']}/{analysis['total_gfp_proteins']}")
            
            print(f"\nAUTOMATED ASSIGNMENT SUMMARY:")
            assignment = analysis['assignment_summary']
            print(f"  Correctly single domain: {assignment['correctly_single_domain']}")
            print(f"  Over-segmented: {assignment['over_segmented']}")
            print(f"  Not classified: {assignment['not_classified']}")
            print(f"  Wrong topology: {assignment['wrong_topology']}")
            
            print(f"\nEVIDENCE SUMMARY:")
            evidence = analysis['evidence_summary']
            print(f"  Proteins with summary files: {evidence['proteins_with_summaries']}/{analysis['total_gfp_proteins']}")
            
        elif args.command == 'problems':
            analysis = analyzer.comprehensive_gfp_analysis()
            
            print("🚨 Problematic GFP Cases")
            print("=" * 40)
            
            if not analysis['problematic_cases']:
                print("✅ No problematic cases found!")
            else:
                for case in analysis['problematic_cases']:
                    print(f"\n{case['protein']}:")
                    print(f"  Domain count: {case['domain_count']} (should be 1)")
                    print(f"  Confidence: {case['confidence']}")
                    print(f"  Issues:")
                    for issue in case['issues']:
                        print(f"    - {issue}")
                        
                print(f"\nSUMMARY: {len(analysis['problematic_cases'])} problematic cases found")
                
                # Show most common issues
                all_issues = []
                for case in analysis['problematic_cases']:
                    all_issues.extend(case['issues'])
                    
                from collections import Counter
                issue_counts = Counter(all_issues)
                
                print(f"\nMOST COMMON ISSUES:")
                for issue, count in issue_counts.most_common(5):
                    print(f"  {count}x: {issue}")
        
        return 0
        
    except Exception as e:
        print(f"❌ Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
