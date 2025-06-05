#!/usr/bin/env python3
"""
Coordinate Translation Validation Script

Validates coordinate translation against known good examples from ECOD database.
Tests both the translation accuracy and identifies potential systematic issues.

Usage:
    python validate_coordinates.py --test-proteins 8ovp_A,2ia4_A,6dgv_A
    python validate_coordinates.py --random-sample 20
    python validate_coordinates.py --batch-dir /path/to/batch --all-proteins
"""

import sys
import os
import random
import psycopg2
import yaml
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
from collections import defaultdict

# Import our coordinate translation
sys.path.insert(0, str(Path(__file__).parent / "../mini/core"))
from visualization import CoordinateTranslator

@dataclass
class ValidationResult:
    """Result of coordinate validation for one protein"""
    protein_id: str
    pdb_id: str
    chain_id: str
    status: str  # 'success', 'failed', 'no_structure', 'no_domains'
    total_residues_mapped: int
    ecod_domains_tested: int
    translation_accuracy: float  # 0.0 - 1.0
    issues: List[str]
    sample_translations: List[Tuple[int, str]]  # (seqid, pdb_num) pairs

@dataclass
class ECODDomain:
    """ECOD domain from database"""
    ecod_uid: int
    protein_id: str
    range_str: str
    t_group: str
    family: str

class CoordinateValidator:
    """Validates coordinate translation against ECOD database"""
    
    def __init__(self, config_path: str = "../config.yml"):
        self.config = self._load_config(config_path)
        self.db_conn = None
        self._connect_database()
    
    def _load_config(self, config_path: str) -> dict:
        """Load database configuration"""
        try:
            with open(config_path, 'r') as f:
                return yaml.safe_load(f)
        except Exception as e:
            print(f"Warning: Could not load config {config_path}: {e}")
            # Default config for testing
            return {
                'database': {
                    'host': 'localhost',
                    'port': 5432,
                    'dbname': 'ecod_protein',
                    'user': 'ecod',
                    'password': 'your_password'
                }
            }
    
    def _connect_database(self):
        """Connect to PostgreSQL database"""
        try:
            db_config = self.config['database']
            self.db_conn = psycopg2.connect(
                host=db_config['host'],
                port=db_config['port'],
                dbname=db_config['dbname'],
                user=db_config['user'],
                password=db_config['password']
            )
            print("✅ Connected to ECOD database")
        except Exception as e:
            print(f"⚠️ Could not connect to database: {e}")
            print("   Validation will be limited to structure-based testing")
    
    def get_test_proteins(self, protein_list: List[str] = None, 
                         random_sample: int = None,
                         batch_dir: str = None) -> List[str]:
        """Get list of proteins to test"""
        
        if protein_list:
            return protein_list
        
        if random_sample and self.db_conn:
            # Get random sample from database
            cursor = self.db_conn.cursor()
            cursor.execute("""
                SELECT DISTINCT p.source_id 
                FROM pdb_analysis.protein p
                JOIN pdb_analysis.domain d ON p.id = d.protein_id
                WHERE p.source_id LIKE '%_%'
                ORDER BY RANDOM()
                LIMIT %s
            """, (random_sample,))
            
            return [row[0] for row in cursor.fetchall()]
        
        if batch_dir:
            # Get proteins from batch directory
            batch_path = Path(batch_dir)
            if batch_path.exists():
                mini_domains = batch_path / "mini_domains"
                if mini_domains.exists():
                    proteins = []
                    for xml_file in mini_domains.glob("*.mini.domains.xml"):
                        protein_id = xml_file.stem.replace('.mini.domains', '')
                        proteins.append(protein_id)
                    return proteins
        
        # Default test set
        return [
            "8ovp_A",  # Known good fusion protein
            "2ia4_A",  # Classic test case
            "6dgv_A",  # Another test case
            "1xyz_A",  # May not exist - test error handling
        ]
    
    def get_ecod_domains(self, protein_id: str) -> List[ECODDomain]:
        """Get ECOD domains for a protein from database"""
        
        if not self.db_conn:
            return []
        
        cursor = self.db_conn.cursor()
        cursor.execute("""
            SELECT d.ecod_uid, p.source_id, d.range, d.t_group, 
                   COALESCE(d.t_group, d.h_group, 'unknown') as family
            FROM pdb_analysis.protein p
            JOIN pdb_analysis.domain d ON p.id = d.protein_id
            WHERE p.source_id = %s
            ORDER BY d.id
        """, (protein_id,))
        
        domains = []
        for row in cursor.fetchall():
            domains.append(ECODDomain(
                ecod_uid=row[0],
                protein_id=row[1],
                range_str=row[2],
                t_group=row[3],
                family=row[4]
            ))
        
        return domains
    
    def validate_protein(self, protein_id: str, 
                        pdb_repo_path: str = "/usr2/pdb/data") -> ValidationResult:
        """Validate coordinate translation for one protein"""
        
        print(f"🧪 Validating {protein_id}...")
        
        pdb_id = protein_id.split('_')[0]
        chain_id = protein_id.split('_')[1] if '_' in protein_id else 'A'
        
        # Find structure file
        structure_path = self._find_structure(pdb_id, pdb_repo_path)
        if not structure_path:
            return ValidationResult(
                protein_id=protein_id,
                pdb_id=pdb_id,
                chain_id=chain_id,
                status='no_structure',
                total_residues_mapped=0,
                ecod_domains_tested=0,
                translation_accuracy=0.0,
                issues=[f"Structure file not found for {pdb_id}"],
                sample_translations=[]
            )
        
        # Test coordinate translation
        try:
            translator = CoordinateTranslator(structure_path, chain_id)
            mapping_info = translator.get_mapping_info()
            
            if mapping_info['total_residues_mapped'] == 0:
                return ValidationResult(
                    protein_id=protein_id,
                    pdb_id=pdb_id,
                    chain_id=chain_id,
                    status='failed',
                    total_residues_mapped=0,
                    ecod_domains_tested=0,
                    translation_accuracy=0.0,
                    issues=['Coordinate mapping failed - no residues mapped'],
                    sample_translations=[]
                )
            
            # Get ECOD domains for validation
            ecod_domains = self.get_ecod_domains(protein_id)
            
            # Test domain range translations
            translation_accuracy = 1.0
            issues = []
            domains_tested = 0
            
            for domain in ecod_domains:
                domains_tested += 1
                accuracy, domain_issues = self._test_domain_translation(
                    domain, translator
                )
                translation_accuracy = min(translation_accuracy, accuracy)
                issues.extend(domain_issues)
            
            # Get sample translations for inspection
            sample_translations = []
            sample_seqids = list(translator.seqid_to_pdb.keys())[:10]
            for seqid in sample_seqids:
                pdb_num = translator.seqid_to_pdb[seqid]
                sample_translations.append((seqid, pdb_num))
            
            # Check for structural issues
            struct_issues = self._check_structural_consistency(translator)
            issues.extend(struct_issues)
            
            status = 'success' if translation_accuracy > 0.8 and len(issues) < 3 else 'failed'
            
            return ValidationResult(
                protein_id=protein_id,
                pdb_id=pdb_id,
                chain_id=chain_id,
                status=status,
                total_residues_mapped=mapping_info['total_residues_mapped'],
                ecod_domains_tested=domains_tested,
                translation_accuracy=translation_accuracy,
                issues=issues,
                sample_translations=sample_translations
            )
            
        except Exception as e:
            return ValidationResult(
                protein_id=protein_id,
                pdb_id=pdb_id,
                chain_id=chain_id,
                status='failed',
                total_residues_mapped=0,
                ecod_domains_tested=0,
                translation_accuracy=0.0,
                issues=[f"Translation error: {str(e)}"],
                sample_translations=[]
            )
    
    def _test_domain_translation(self, domain: ECODDomain, 
                               translator: CoordinateTranslator) -> Tuple[float, List[str]]:
        """Test translation accuracy for a specific domain"""
        
        issues = []
        
        try:
            # Parse domain range
            segments = self._parse_range(domain.range_str)
            if not segments:
                return 0.0, [f"Could not parse range: {domain.range_str}"]
            
            # Test each segment translation
            successful_translations = 0
            total_segments = len(segments)
            
            for start, end in segments:
                pdb_start, pdb_end = translator.translate_seqid_range(start, end)
                
                # Check if translation succeeded
                if pdb_start and pdb_end:
                    successful_translations += 1
                    
                    # Validate translation makes sense
                    if not self._validate_pdb_range(pdb_start, pdb_end):
                        issues.append(f"Invalid PDB range: {pdb_start}-{pdb_end}")
                else:
                    issues.append(f"Translation failed for seqid {start}-{end}")
            
            accuracy = successful_translations / total_segments if total_segments > 0 else 0.0
            return accuracy, issues
            
        except Exception as e:
            return 0.0, [f"Domain translation error: {str(e)}"]
    
    def _parse_range(self, range_str: str) -> List[Tuple[int, int]]:
        """Parse range string into segments"""
        segments = []
        
        try:
            for segment in range_str.split(','):
                segment = segment.strip()
                if '-' in segment:
                    start, end = segment.split('-')
                    segments.append((int(start), int(end)))
                else:
                    pos = int(segment)
                    segments.append((pos, pos))
        except Exception:
            return []
        
        return segments
    
    def _validate_pdb_range(self, pdb_start: str, pdb_end: str) -> bool:
        """Validate that PDB range makes sense"""
        
        import re
        
        try:
            # Extract numeric parts
            start_num = int(re.sub(r'[A-Z]', '', pdb_start))
            end_num = int(re.sub(r'[A-Z]', '', pdb_end))
            
            # End should be >= start
            return end_num >= start_num
            
        except Exception:
            return False
    
    def _check_structural_consistency(self, translator: CoordinateTranslator) -> List[str]:
        """Check for structural consistency issues"""
        
        issues = []
        
        # Check for large gaps in seqid sequence
        seqids = sorted(translator.seqid_to_pdb.keys())
        large_gaps = []
        
        for i in range(len(seqids) - 1):
            gap = seqids[i+1] - seqids[i]
            if gap > 10:  # Gap of more than 10 residues
                large_gaps.append((seqids[i], seqids[i+1], gap))
        
        if large_gaps:
            issues.append(f"Large gaps in sequence: {len(large_gaps)} gaps > 10 residues")
        
        # Check for insertion codes (expected in some structures)
        has_insertions = any(not pdb_num.isdigit() 
                           for pdb_num in translator.seqid_to_pdb.values())
        if has_insertions:
            issues.append("Structure has insertion codes (expected for some PDB entries)")
        
        return issues
    
    def _find_structure(self, pdb_id: str, pdb_repo_path: str) -> Optional[str]:
        """Find structure file"""
        pdb_id = pdb_id.lower()
        
        possible_paths = [
            f"{pdb_repo_path}/structures/divided/mmCIF/{pdb_id[1:3]}/{pdb_id}.cif.gz",
            f"{pdb_repo_path}/structures/divided/mmCIF/{pdb_id[1:3]}/{pdb_id}.cif",
            f"{pdb_repo_path}/mmCIF/{pdb_id}.cif.gz",
            f"{pdb_repo_path}/mmCIF/{pdb_id}.cif",
        ]
        
        for path in possible_paths:
            if os.path.exists(path):
                return path
        
        return None
    
    def run_validation(self, proteins: List[str]) -> Dict:
        """Run validation on list of proteins"""
        
        print(f"🎯 Validating coordinate translation for {len(proteins)} proteins")
        
        results = []
        summary = {
            'total_proteins': len(proteins),
            'successful': 0,
            'failed': 0,
            'no_structure': 0,
            'no_domains': 0,
            'avg_accuracy': 0.0,
            'common_issues': defaultdict(int)
        }
        
        for protein_id in proteins:
            result = self.validate_protein(protein_id)
            results.append(result)
            
            # Update summary
            summary[result.status] += 1
            
            if result.status == 'success':
                summary['avg_accuracy'] += result.translation_accuracy
            
            # Track common issues
            for issue in result.issues:
                # Categorize issues
                if 'gap' in issue.lower():
                    summary['common_issues']['sequence_gaps'] += 1
                elif 'insertion' in issue.lower():
                    summary['common_issues']['insertion_codes'] += 1
                elif 'translation failed' in issue.lower():
                    summary['common_issues']['translation_failures'] += 1
                elif 'structure' in issue.lower():
                    summary['common_issues']['structure_issues'] += 1
                else:
                    summary['common_issues']['other'] += 1
        
        # Calculate average accuracy
        if summary['successful'] > 0:
            summary['avg_accuracy'] /= summary['successful']
        
        return {
            'summary': summary,
            'results': results
        }
    
    def generate_report(self, validation_data: Dict, output_file: str):
        """Generate validation report"""
        
        summary = validation_data['summary']
        results = validation_data['results']
        
        lines = [
            "# Coordinate Translation Validation Report",
            "",
            f"**Validation Date:** {os.popen('date').read().strip()}",
            f"**Proteins Tested:** {summary['total_proteins']}",
            "",
            "## 📊 Summary Statistics",
            "",
            f"- **Successful:** {summary['successful']} ({summary['successful']/summary['total_proteins']*100:.1f}%)",
            f"- **Failed:** {summary['failed']} ({summary['failed']/summary['total_proteins']*100:.1f}%)",
            f"- **No Structure:** {summary['no_structure']} ({summary['no_structure']/summary['total_proteins']*100:.1f}%)",
            f"- **Average Accuracy:** {summary['avg_accuracy']:.1%}",
            "",
            "## 🚨 Common Issues",
            ""
        ]
        
        for issue_type, count in summary['common_issues'].items():
            percentage = count / summary['total_proteins'] * 100
            lines.append(f"- **{issue_type.replace('_', ' ').title()}:** {count} cases ({percentage:.1f}%)")
        
        lines.extend([
            "",
            "## ✅ Successful Cases",
            ""
        ])
        
        successful_cases = [r for r in results if r.status == 'success']
        for result in successful_cases[:10]:  # Top 10
            lines.append(f"- **{result.protein_id}**: {result.total_residues_mapped} residues mapped, "
                        f"{result.ecod_domains_tested} domains tested, "
                        f"{result.translation_accuracy:.1%} accuracy")
        
        lines.extend([
            "",
            "## ❌ Failed Cases",
            ""
        ])
        
        failed_cases = [r for r in results if r.status == 'failed']
        for result in failed_cases[:10]:  # Top 10 failures
            lines.append(f"- **{result.protein_id}**: {', '.join(result.issues[:2])}")
        
        lines.extend([
            "",
            "## 🔍 Sample Coordinate Mappings",
            ""
        ])
        
        # Show sample translations from successful cases
        for result in successful_cases[:3]:
            if result.sample_translations:
                lines.append(f"### {result.protein_id}")
                lines.append("| SeqID | PDB Residue |")
                lines.append("|-------|-------------|")
                for seqid, pdb_num in result.sample_translations[:5]:
                    lines.append(f"| {seqid} | {pdb_num} |")
                lines.append("")
        
        lines.extend([
            "",
            "## 🎯 Recommendations",
            ""
        ])
        
        if summary['successful'] / summary['total_proteins'] > 0.8:
            lines.append("✅ **Coordinate translation is working well** (>80% success rate)")
        else:
            lines.append("⚠️ **Coordinate translation needs improvement** (<80% success rate)")
        
        if summary['common_issues']['structure_issues'] > summary['total_proteins'] * 0.2:
            lines.append("- Consider improving structure file detection/parsing")
        
        if summary['common_issues']['translation_failures'] > summary['total_proteins'] * 0.1:
            lines.append("- Review mmCIF parsing logic for edge cases")
        
        if summary['no_structure'] > summary['total_proteins'] * 0.3:
            lines.append("- Check PDB repository path and file organization")
        
        with open(output_file, 'w') as f:
            f.write('\n'.join(lines))

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Validate coordinate translation')
    parser.add_argument('--test-proteins', help='Comma-separated list of protein IDs')
    parser.add_argument('--random-sample', type=int, help='Test N random proteins from database')
    parser.add_argument('--batch-dir', help='Test all proteins from batch directory')
    parser.add_argument('--all-proteins', action='store_true', help='Test all proteins in batch-dir')
    parser.add_argument('--config', default='../config.yml', help='Database config file')
    parser.add_argument('--output', default='/tmp/coordinate_validation.md', help='Output report file')
    parser.add_argument('--pdb-repo', default='/usr2/pdb/data', help='PDB repository path')
    
    args = parser.parse_args()
    
    # Determine test set
    proteins = []
    if args.test_proteins:
        proteins = args.test_proteins.split(',')
    elif args.random_sample:
        proteins = None  # Will be determined by validator
    elif args.batch_dir:
        proteins = None  # Will be determined by validator
    else:
        proteins = ["8ovp_A", "2ia4_A", "6dgv_A"]  # Default test set
    
    # Run validation
    validator = CoordinateValidator(args.config)
    
    if proteins is None:
        proteins = validator.get_test_proteins(
            random_sample=args.random_sample,
            batch_dir=args.batch_dir
        )
    
    validation_data = validator.run_validation(proteins)
    
    # Generate report
    validator.generate_report(validation_data, args.output)
    
    # Print summary
    summary = validation_data['summary']
    print(f"\n📋 Validation Complete!")
    print(f"   Tested: {summary['total_proteins']} proteins")
    print(f"   Success: {summary['successful']} ({summary['successful']/summary['total_proteins']*100:.1f}%)")
    print(f"   Failed: {summary['failed']} ({summary['failed']/summary['total_proteins']*100:.1f}%)")
    print(f"   Avg Accuracy: {summary['avg_accuracy']:.1%}")
    print(f"   Report: {args.output}")
    
    if validator.db_conn:
        validator.db_conn.close()

if __name__ == "__main__":
    main()
