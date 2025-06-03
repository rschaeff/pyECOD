#!/usr/bin/env python3
"""
ECOD T-group validation test for 8ovp_A

This test validates that mini correctly assigns domains to ECOD T-groups
using the actual ECOD domains.txt classifications.
"""

import pytest
from pathlib import Path
from typing import Dict, List, Set

# Add parent directory for imports
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.parser import parse_domain_summary, load_reference_lengths, load_protein_lengths
from mini.decomposer import load_domain_definitions
from mini.blast_parser import load_chain_blast_alignments
from mini.partitioner import partition_domains
from mini.ecod_domains_parser import load_ecod_classifications

class TestEcodTGroupValidation:
    """ECOD T-group validation tests"""
    
    @pytest.mark.integration
    @pytest.mark.slow
    def test_8ovp_A_tgroup_assignment(self, stable_batch_dir, real_reference_data, 
                                     blast_alignments, temp_output_dir):
        """
        PRODUCTION TEST: Validate ECOD T-group assignments for 8ovp_A
        
        This tests what actually matters for ECOD production:
        - Correct T-group classification (not superficial family names)
        - Biological domain architecture accuracy
        - ECOD database consistency
        """
        protein_id = "8ovp_A"
        
        # Run the algorithm
        result = self._run_mini_algorithm(protein_id, stable_batch_dir, real_reference_data, blast_alignments)
        assert result['success'], f"Algorithm failed: {result.get('error', 'Unknown error')}"
        
        domains = result['domains']
        
        # Load ECOD classifications for validation
        ecod_file = Path("test_data/ecod_classifications.csv")
        if ecod_file.exists():
            ecod_classifications = load_ecod_classifications(str(ecod_file))
            tgroup_validation = self._validate_tgroup_assignments(domains, ecod_classifications)
        else:
            # Fallback validation without ECOD data
            tgroup_validation = self._validate_domain_characteristics(domains)
        
        # Print detailed results
        print(f"\\n=== ECOD T-GROUP VALIDATION RESULTS ===")
        print(f"Protein: {protein_id}")
        print(f"Domains found: {len(domains)}")
        
        for i, domain in enumerate(domains, 1):
            disc = " (discontinuous)" if domain.range.is_discontinuous else ""
            print(f"  {i}. {domain.family}: {domain.range} ({domain.range.total_length} res){disc}")
        
        print(f"\\nValidation Results:")
        for check, passed in tgroup_validation['checks'].items():
            status = "✅" if passed else "❌"
            print(f"  {status} {check}")
        
        if tgroup_validation.get('tgroup_details'):
            print(f"\\nT-group Details:")
            for detail in tgroup_validation['tgroup_details']:
                print(f"  {detail}")
        
        # Core assertions for production readiness
        assert tgroup_validation['checks']['domain_count'], "Incorrect domain count"
        assert tgroup_validation['checks']['gfp_domain'], "GFP domain not found"
        assert tgroup_validation['checks']['pbp_tgroup'], "PBP T-group not correctly assigned"
        assert tgroup_validation['checks']['discontinuous_architecture'], "Missing discontinuous architecture"
        assert tgroup_validation['checks']['coverage'], "Insufficient coverage"
        
        print(f"\\n✅ ECOD T-GROUP VALIDATION PASSED")
        print(f"   Mini correctly assigns domains to ECOD T-groups")
        print(f"   Biological architecture properly detected")
        print(f"   Production-quality results achieved")
    
    def _run_mini_algorithm(self, protein_id: str, batch_dir: str, reference_data: Dict, blast_alignments: Dict) -> Dict:
        """Run mini algorithm and return results"""
        import os
        
        parts = protein_id.split('_')
        pdb_id, chain_id = parts[0], parts[1] if len(parts) > 1 else 'A'
        
        xml_path = os.path.join(batch_dir, "domains", f"{protein_id}.develop291.domain_summary.xml")
        if not os.path.exists(xml_path):
            return {'success': False, 'error': f"Domain summary not found: {xml_path}"}
        
        try:
            # Parse evidence
            evidence = parse_domain_summary(
                xml_path,
                reference_lengths=reference_data.get('domain_lengths', {}),
                protein_lengths=reference_data.get('protein_lengths', {}),
                blast_alignments=blast_alignments,
                require_reference_lengths=True,
                verbose=False
            )
            
            if not evidence:
                return {'success': False, 'error': 'No evidence with reference lengths found'}
            
            # Estimate sequence length
            max_pos = max(ev.query_range.segments[-1].end for ev in evidence)
            sequence_length = int(max_pos * 1.1)
            
            # Partition domains
            domains = partition_domains(
                evidence,
                sequence_length=sequence_length,
                domain_definitions=reference_data.get('domain_definitions', {}),
                verbose=False
            )
            
            return {'success': True, 'domains': domains, 'sequence_length': sequence_length}
            
        except Exception as e:
            return {'success': False, 'error': str(e)}
    
    def _validate_tgroup_assignments(self, domains: List, ecod_classifications: Dict) -> Dict:
        """Validate T-group assignments using actual ECOD T-groups"""

        validation = {
            'checks': {},
            'tgroup_details': []
        }

        # Expected T-groups for 8ovp
        expected_pbp_tgroup = "7523.1.1"  # Periplasmic binding protein-like II

        for domain in domains:
            # Get the actual T-group from the evidence
            if hasattr(domain, 'evidence_items') and domain.evidence_items:
                t_group = domain.evidence_items[0].t_group if domain.evidence_items[0].t_group else "unknown"
            else:
                t_group = "unknown"

            validation['tgroup_details'].append(
                f"{domain.family}: T-group={t_group}, range={domain.range}"
            )

            # Check if it's the expected PBP T-group
            if t_group.startswith(expected_pbp_tgroup):
                validation['pbp_domains_found'] = validation.get('pbp_domains_found', 0) + 1

        # Validate results
        validation['checks']['correct_pbp_tgroup'] = validation.get('pbp_domains_found', 0) >= 2
        
        return validation

        def _validate_domain_characteristics(self, domains: List) -> Dict:
            """Fallback validation without ECOD data"""

            validation = {
                'checks': {},
                'tgroup_details': ['Fallback validation (no ECOD data available)']
            }

            # Basic structural validation
            validation['checks']['domain_count'] = len(domains) == 3
            validation['checks']['gfp_domain'] = any('6dgv' in d.family.lower() for d in domains)
            validation['checks']['pbp_tgroup'] = any('2vha' in d.family.lower() or '2ia4' in d.family.lower() for d in domains)
            validation['checks']['discontinuous_architecture'] = any(d.range.is_discontinuous for d in domains)
            validation['checks']['coverage'] = sum(d.range.total_length for d in domains) >= 450
            validation['checks']['decomposition_success'] = sum(1 for d in domains if d.source == 'chain_blast_decomposed') >= 1

            return validation

if __name__ == "__main__":
    # Run the test directly
    print("Running ECOD T-group validation test...")
    pytest.main([__file__, "-v"])
