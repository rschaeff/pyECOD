#!/usr/bin/env python3
"""
Find diverse test proteins for mini_pyecod test suite expansion

This script queries ecod_schema to find proteins with domain summaries
and analyzes their T-group hits to identify diverse test cases.

Usage:
    python find_test_proteins.py --config config.yml --batch-id 36
    python find_test_proteins.py --config config.yml --analyze-tgroups
    python find_test_proteins.py --config config.yml --output test_proteins.csv
"""

import psycopg2
import yaml
import argparse
import csv
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass
from collections import defaultdict


# Common ECOD T-groups from classic proteins for reference
REFERENCE_TGROUPS = {
    "221.1.1": "Ubiquitin-like",
    "109.1.1": "Thioredoxin-like", 
    "2484.1.1": "TIM barrel",
    "1.1.1": "Globin-like",
    "2.1.1": "SH3-like barrel",
    "8040.1.1": "SH2 domain",
    "218.1.1": "Ankyrin repeat",
    "7579.1.1": "Periplasmic binding protein-like II",  # PBP
    "5025.1.1": "GFP-like"  # GFP
}


@dataclass
class TestProteinCandidate:
    """Candidate protein for test suite"""
    protein_id: str
    pdb_id: str
    chain_id: str
    length: int
    batch_id: int
    process_id: int
    # From domain summary analysis
    evidence_counts: Dict[str, int]
    tgroup_hits: Set[str]
    domain_count_estimate: int
    has_discontinuous: bool = False
    has_chain_blast_decomposition: bool = False
    description: str = ""
    category: str = ""
    priority: int = 0


class EcodSchemaTestFinder:
    """Find test proteins from ecod_schema that have domain summaries"""
    
    def __init__(self, db_config: Dict, batch_id: int = 36):
        self.db_config = db_config
        self.batch_id = batch_id
        self.conn = None
        
        # Categories based on what we find
        self.categories = {
            "single_tgroup": "Proteins with hits to single T-group",
            "multi_tgroup": "Proteins with hits to multiple T-groups", 
            "ubiquitin_like": "Proteins with ubiquitin-like hits (221.1.1)",
            "tim_barrel": "Proteins with TIM barrel hits (2484.1.1)",
            "pbp_like": "Proteins with PBP-like hits (7579.1.1)",
            "repeat_proteins": "Proteins with repeat domain hits",
            "large_proteins": "Large proteins (>800 residues)",
            "small_proteins": "Small proteins (<150 residues)",
            "decomposition_candidates": "Good chain BLAST decomposition candidates"
        }
    
    def connect(self):
        """Connect to database"""
        try:
            self.conn = psycopg2.connect(
                host=self.db_config['host'],
                database=self.db_config['database'],
                user=self.db_config['user'],
                password=self.db_config['password']
            )
            print(f"Connected to database: {self.db_config['database']}")
        except Exception as e:
            print(f"ERROR: Failed to connect to database: {e}")
            raise
    
    def find_proteins_with_domain_summaries(self) -> List[Dict]:
        """Find proteins in batch that have domain summary files"""
        
        query = """
        WITH batch_proteins AS (
            SELECT 
                p.id as protein_id,
                p.pdb_id,
                p.chain_id,
                p.source_id,
                p.length,
                ps.id as process_id,
                ps.batch_id,
                ps.current_stage,
                ps.status,
                ps.relative_path
            FROM ecod_schema.protein p
            JOIN ecod_schema.process_status ps ON p.id = ps.protein_id
            WHERE ps.batch_id = %s
        ),
        domain_files AS (
            SELECT 
                pf.process_id,
                pf.file_type,
                pf.file_path,
                pf.file_exists
            FROM ecod_schema.process_file pf
            WHERE pf.file_type = 'domain_summary'
                AND pf.file_exists = true
        )
        SELECT 
            bp.protein_id,
            bp.pdb_id,
            bp.chain_id,
            bp.source_id,
            bp.length,
            bp.process_id,
            bp.batch_id,
            bp.current_stage,
            bp.status,
            df.file_path
        FROM batch_proteins bp
        JOIN domain_files df ON bp.process_id = df.process_id
        WHERE bp.length IS NOT NULL
        ORDER BY bp.length DESC;
        """
        
        proteins = []
        with self.conn.cursor() as cur:
            cur.execute(query, (self.batch_id,))
            
            for row in cur.fetchall():
                proteins.append({
                    'protein_id': row[0],
                    'pdb_id': row[1],
                    'chain_id': row[2],
                    'source_id': row[3],
                    'length': row[4],
                    'process_id': row[5],
                    'batch_id': row[6],
                    'current_stage': row[7],
                    'status': row[8],
                    'domain_summary_path': row[9]
                })
        
        print(f"Found {len(proteins)} proteins with domain summaries in batch {self.batch_id}")
        return proteins
    
    def analyze_domain_summary(self, xml_path: str) -> Dict:
        """Analyze a domain summary to extract T-group hits and evidence"""
        
        result = {
            'evidence_counts': {
                'chain_blast': 0,
                'domain_blast': 0,
                'hhsearch': 0
            },
            'tgroup_hits': set(),
            'domain_ids': set(),
            'has_discontinuous': False,
            'estimated_domains': 0
        }
        
        try:
            # Handle both direct paths and paths that need batch prefix
            if not xml_path.startswith('/'):
                # Relative path - need to construct full path
                batch_base = f"/data/ecod/pdb_updates/batches/ecod_batch_{str(self.batch_id).zfill(3)}_*"
                # This is a simplification - in practice you'd resolve the full batch name
                xml_path = f"/data/ecod/pdb_updates/batches/{xml_path}"
            
            tree = ET.parse(xml_path)
            root = tree.getroot()
            
            # Count evidence
            result['evidence_counts']['chain_blast'] = len(root.findall(".//chain_blast_run/hits/hit"))
            result['evidence_counts']['domain_blast'] = len(root.findall(".//blast_run/hits/hit"))
            result['evidence_counts']['hhsearch'] = len(root.findall(".//hh_run/hits/hit"))
            
            # Extract T-groups from hits
            # From domain BLAST hits
            for hit in root.findall(".//blast_run/hits/hit"):
                t_group = hit.get("t_group")
                if t_group:
                    result['tgroup_hits'].add(t_group)
                
                domain_id = hit.get("domain_id")
                if domain_id:
                    result['domain_ids'].add(domain_id)
            
            # From HHsearch hits - might have t_group attribute or need to parse from hit_id
            for hit in root.findall(".//hh_run/hits/hit"):
                t_group = hit.get("t_group")
                if t_group:
                    result['tgroup_hits'].add(t_group)
                
                # Sometimes T-group info is embedded in the hit description
                hit_id = hit.get("hit_id", "")
                if hit_id and hit_id.startswith("e"):
                    result['domain_ids'].add(hit_id)
            
            # Check for discontinuous ranges
            for hit in root.findall(".//*/hit"):
                query_reg = hit.find("query_reg")
                if query_reg is not None and query_reg.text and ',' in query_reg.text:
                    result['has_discontinuous'] = True
                    break
            
            # Estimate domain count from unique significant hits
            high_conf_domains = set()
            
            # High confidence domain BLAST hits
            for hit in root.findall(".//blast_run/hits/hit"):
                evalue = float(hit.get("evalues", "999"))
                if evalue < 1e-10:
                    domain_id = hit.get("domain_id", "")
                    if domain_id:
                        high_conf_domains.add(domain_id.split('_')[0])  # Group by PDB
            
            # High confidence HHsearch hits
            for hit in root.findall(".//hh_run/hits/hit"):
                prob = float(hit.get("probability", "0"))
                if prob > 90:
                    hit_id = hit.get("hit_id", "")
                    if hit_id:
                        high_conf_domains.add(hit_id.split('_')[0])
            
            result['estimated_domains'] = len(high_conf_domains) if high_conf_domains else 1
            
        except Exception as e:
            print(f"Warning: Failed to parse {xml_path}: {e}")
        
        return result
    
    def categorize_proteins(self, proteins: List[Dict]) -> Dict[str, List[TestProteinCandidate]]:
        """Categorize proteins based on their characteristics"""
        
        categorized = defaultdict(list)
        
        for protein in proteins:
            # Analyze domain summary
            analysis = self.analyze_domain_summary(protein['domain_summary_path'])
            
            # Create candidate
            candidate = TestProteinCandidate(
                protein_id=f"{protein['pdb_id']}_{protein['chain_id']}",
                pdb_id=protein['pdb_id'],
                chain_id=protein['chain_id'],
                length=protein['length'],
                batch_id=protein['batch_id'],
                process_id=protein['process_id'],
                evidence_counts=analysis['evidence_counts'],
                tgroup_hits=analysis['tgroup_hits'],
                domain_count_estimate=analysis['estimated_domains'],
                has_discontinuous=analysis['has_discontinuous']
            )
            
            # Categorize by size
            if candidate.length < 150:
                categorized['small_proteins'].append(candidate)
            elif candidate.length > 800:
                categorized['large_proteins'].append(candidate)
            
            # Categorize by T-group hits
            if len(candidate.tgroup_hits) == 0:
                continue  # Skip if no T-group hits
            elif len(candidate.tgroup_hits) == 1:
                categorized['single_tgroup'].append(candidate)
            else:
                categorized['multi_tgroup'].append(candidate)
            
            # Check for specific T-groups of interest
            for tgroup in candidate.tgroup_hits:
                if tgroup.startswith("221."):  # Ubiquitin-like
                    categorized['ubiquitin_like'].append(candidate)
                elif tgroup.startswith("2484."):  # TIM barrel
                    categorized['tim_barrel'].append(candidate)
                elif tgroup.startswith("7579."):  # PBP-like
                    categorized['pbp_like'].append(candidate)
                elif tgroup.startswith("218."):  # Ankyrin repeat
                    categorized['repeat_proteins'].append(candidate)
            
            # Good decomposition candidates: multi-domain with chain BLAST hits
            if (candidate.domain_count_estimate > 1 and 
                candidate.evidence_counts['chain_blast'] > 0):
                categorized['decomposition_candidates'].append(candidate)
        
        return categorized
    
    def select_diverse_test_set(self, categorized: Dict[str, List[TestProteinCandidate]], 
                               target_count: int = 15) -> List[TestProteinCandidate]:
        """Select a diverse set of test proteins"""
        
        selected = []
        used_proteins = set()
        
        # Priority selection strategy
        selection_targets = [
            ('decomposition_candidates', 3),  # Like 8ovp_A
            ('single_tgroup', 3),             # Simple single domain
            ('multi_tgroup', 3),              # Multi-domain  
            ('small_proteins', 2),            # Size extremes
            ('large_proteins', 2),            # Size extremes
            ('ubiquitin_like', 1),            # Common fold
            ('repeat_proteins', 1),           # Complex architecture
        ]
        
        for category, target in selection_targets:
            candidates = categorized.get(category, [])
            
            # Sort by evidence quality
            sorted_candidates = sorted(
                candidates,
                key=lambda x: (
                    sum(x.evidence_counts.values()),  # Total evidence
                    len(x.tgroup_hits),                # T-group diversity
                    -x.length                          # Prefer shorter for testing
                ),
                reverse=True
            )
            
            # Select up to target, avoiding duplicates
            added = 0
            for candidate in sorted_candidates:
                if candidate.protein_id not in used_proteins and added < target:
                    candidate.category = category
                    candidate.description = f"{self.categories.get(category, category)}"
                    selected.append(candidate)
                    used_proteins.add(candidate.protein_id)
                    added += 1
            
            if len(selected) >= target_count:
                break
        
        return selected[:target_count]
    
    def export_test_proteins(self, candidates: List[TestProteinCandidate], output_file: str):
        """Export test proteins to CSV"""
        
        with open(output_file, 'w', newline='') as f:
            fieldnames = [
                'protein_id', 'pdb_id', 'chain_id', 'category', 'length',
                'domain_count_estimate', 'evidence_total', 'chain_blast',
                'domain_blast', 'hhsearch', 'tgroup_count', 'tgroup_hits',
                'has_discontinuous', 'description'
            ]
            
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            
            for candidate in candidates:
                writer.writerow({
                    'protein_id': candidate.protein_id,
                    'pdb_id': candidate.pdb_id,
                    'chain_id': candidate.chain_id,
                    'category': candidate.category,
                    'length': candidate.length,
                    'domain_count_estimate': candidate.domain_count_estimate,
                    'evidence_total': sum(candidate.evidence_counts.values()),
                    'chain_blast': candidate.evidence_counts['chain_blast'],
                    'domain_blast': candidate.evidence_counts['domain_blast'],
                    'hhsearch': candidate.evidence_counts['hhsearch'],
                    'tgroup_count': len(candidate.tgroup_hits),
                    'tgroup_hits': '|'.join(sorted(candidate.tgroup_hits)),
                    'has_discontinuous': candidate.has_discontinuous,
                    'description': candidate.description
                })
        
        print(f"\nExported {len(candidates)} test proteins to {output_file}")
    
    def analyze_tgroup_distribution(self, proteins: List[Dict]) -> Dict[str, int]:
        """Analyze T-group distribution across all proteins"""
        
        tgroup_counts = defaultdict(int)
        
        for protein in proteins:
            analysis = self.analyze_domain_summary(protein['domain_summary_path'])
            for tgroup in analysis['tgroup_hits']:
                tgroup_counts[tgroup] += 1
        
        return dict(sorted(tgroup_counts.items(), key=lambda x: x[1], reverse=True))


def main():
    parser = argparse.ArgumentParser(
        description='Find diverse test proteins from ecod_schema batches'
    )
    
    parser.add_argument('--config', required=True,
                       help='Database configuration file (YAML)')
    parser.add_argument('--batch-id', type=int, default=36,
                       help='Batch ID to analyze')
    parser.add_argument('--output', default='test_proteins_from_batch.csv',
                       help='Output CSV file')
    parser.add_argument('--target-count', type=int, default=15,
                       help='Target number of test proteins')
    parser.add_argument('--analyze-tgroups', action='store_true',
                       help='Analyze T-group distribution in batch')
    
    args = parser.parse_args()
    
    # Load database configuration
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    
    db_config = config['database']
    
    # Find test proteins
    finder = EcodSchemaTestFinder(db_config, args.batch_id)
    
    try:
        finder.connect()
        
        # Find all proteins with domain summaries
        proteins = finder.find_proteins_with_domain_summaries()
        
        if not proteins:
            print(f"No proteins with domain summaries found in batch {args.batch_id}")
            return
        
        # Analyze T-group distribution if requested
        if args.analyze_tgroups:
            print("\nAnalyzing T-group distribution...")
            tgroup_dist = finder.analyze_tgroup_distribution(proteins[:100])  # Sample first 100
            
            print("\nTop 20 T-groups in batch:")
            for tgroup, count in list(tgroup_dist.items())[:20]:
                ref_name = REFERENCE_TGROUPS.get(tgroup, "")
                print(f"  {tgroup}: {count} proteins {f'({ref_name})' if ref_name else ''}")
        
        # Categorize proteins
        print("\nCategorizing proteins...")
        categorized = finder.categorize_proteins(proteins)
        
        # Print category summary
        print("\nProteins by category:")
        for category, candidates in categorized.items():
            if candidates:
                print(f"  {category}: {len(candidates)} proteins")
        
        # Select diverse test set
        selected = finder.select_diverse_test_set(categorized, args.target_count)
        
        # Export results
        finder.export_test_proteins(selected, args.output)
        
        # Print selected proteins
        print("\nSELECTED TEST PROTEINS:")
        print("=" * 80)
        
        for candidate in selected:
            print(f"\n{candidate.protein_id}:")
            print(f"  Category: {candidate.category}")
            print(f"  Length: {candidate.length} residues")
            print(f"  Estimated domains: {candidate.domain_count_estimate}")
            print(f"  Evidence: CB={candidate.evidence_counts['chain_blast']}, "
                  f"DB={candidate.evidence_counts['domain_blast']}, "
                  f"HH={candidate.evidence_counts['hhsearch']}")
            print(f"  T-groups ({len(candidate.tgroup_hits)}): {', '.join(sorted(candidate.tgroup_hits))}")
        
    finally:
        if finder.conn:
            finder.conn.close()


if __name__ == "__main__":
    main()
