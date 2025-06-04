#!/usr/bin/env python3
"""
Collision Detection for Mini PyECOD Import

Analyzes potential data collisions between mini results and existing
partition_proteins/partition_domains records before importing.

Usage:
    python detect_collisions.py                    # Check all completed results
    python detect_collisions.py --batch-name xyz   # Check specific batch
    python detect_collisions.py --summary          # Summary report only
"""

import sys
import argparse
import yaml
import psycopg2
import psycopg2.extras
from pathlib import Path
from typing import Dict, List, Tuple
from datetime import datetime
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class CollisionDetector:
    """Detect potential data collisions for mini import"""
    
    def __init__(self, config_path: str = "config.local.yml"):
        self.config = self._load_config(config_path)
        self.db_conn = self._init_db_connection()
    
    def _load_config(self, config_path: str) -> Dict:
        """Load configuration"""
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    
    def _init_db_connection(self):
        """Initialize PostgreSQL connection"""
        try:
            db_config = self.config['database']
            conn = psycopg2.connect(**db_config)
            logger.info("✓ Connected to production database")
            return conn
        except Exception as e:
            logger.error(f"Database connection failed: {e}")
            raise
    
    def scan_mini_results(self, batch_name: str = None) -> List[Tuple[str, str, str]]:
        """Scan mini results and extract protein identifiers
        
        Returns:
            List of (pdb_id, chain_id, batch_name) tuples
        """
        batch_base = Path(self.config["paths"]["batch_base_dir"])
        mini_proteins = []
        
        batch_dirs = [batch_base / batch_name] if batch_name else batch_base.iterdir()
        
        for batch_dir in batch_dirs:
            if not batch_dir.is_dir():
                continue
                
            mini_domains_dir = batch_dir / "mini_domains"
            if not mini_domains_dir.exists():
                continue
            
            xml_files = list(mini_domains_dir.glob("*.mini.domains.xml"))
            
            for xml_file in xml_files:
                # Extract protein info from filename: "8ovp_A.mini.domains.xml"
                protein_id = xml_file.stem.replace('.mini.domains', '')
                parts = protein_id.split('_')
                pdb_id = parts[0]
                chain_id = parts[1] if len(parts) > 1 else 'A'
                
                mini_proteins.append((pdb_id, chain_id, batch_dir.name))
        
        logger.info(f"Found {len(mini_proteins)} mini results to check")
        return mini_proteins
    
    def check_database_collisions(self, mini_proteins: List[Tuple[str, str, str]]) -> Dict:
        """Check for existing records in partition_proteins that would collide"""
        
        collision_data = {
            'total_mini_proteins': len(mini_proteins),
            'existing_records': [],
            'collision_summary': {},
            'batch_mapping': {}
        }
        
        # Get batch mappings
        with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            cursor.execute("""
                SELECT id, batch_name FROM ecod_schema.batch 
                ORDER BY id
            """)
            batches = cursor.fetchall()
            collision_data['batch_mapping'] = {b['batch_name']: b['id'] for b in batches}
        
        # Check each mini protein against existing records
        for pdb_id, chain_id, batch_name in mini_proteins:
            batch_id = collision_data['batch_mapping'].get(batch_name)
            
            if not batch_id:
                logger.warning(f"No batch_id found for {batch_name}")
                continue
            
            with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
                cursor.execute("""
                    SELECT pp.id, pp.pdb_id, pp.chain_id, pp.batch_id, pp.timestamp,
                           pp.process_version, pp.reference_version, pp.is_classified,
                           pp.domains_with_evidence, pp.fully_classified_domains,
                           COUNT(pd.id) as domain_count
                    FROM pdb_analysis.partition_proteins pp
                    LEFT JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
                    WHERE pp.pdb_id = %s AND pp.chain_id = %s AND pp.batch_id = %s
                    GROUP BY pp.id, pp.pdb_id, pp.chain_id, pp.batch_id, pp.timestamp,
                             pp.process_version, pp.reference_version, pp.is_classified,
                             pp.domains_with_evidence, pp.fully_classified_domains
                    ORDER BY pp.timestamp DESC
                """, (pdb_id, chain_id, batch_id))
                
                existing = cursor.fetchall()
                
                if existing:
                    collision_info = {
                        'mini_protein': f"{pdb_id}_{chain_id}",
                        'batch_name': batch_name,
                        'batch_id': batch_id,
                        'existing_records': []
                    }
                    
                    for record in existing:
                        collision_info['existing_records'].append({
                            'protein_id': record['id'],
                            'process_version': record['process_version'],
                            'reference_version': record['reference_version'],
                            'is_classified': record['is_classified'],
                            'domains_with_evidence': record['domains_with_evidence'],
                            'fully_classified_domains': record['fully_classified_domains'],
                            'domain_count': record['domain_count'],
                            'timestamp': record['timestamp'].isoformat() if record['timestamp'] else None
                        })
                    
                    collision_data['existing_records'].append(collision_info)
        
        # Generate collision summary
        collision_summary = {
            'total_collisions': len(collision_data['existing_records']),
            'process_versions': {},
            'reference_versions': {},
            'classification_status': {'classified': 0, 'unclassified': 0},
            'domain_counts': []
        }
        
        for collision in collision_data['existing_records']:
            for record in collision['existing_records']:
                # Count process versions
                pv = record['process_version'] or 'unknown'
                collision_summary['process_versions'][pv] = collision_summary['process_versions'].get(pv, 0) + 1
                
                # Count reference versions  
                rv = record['reference_version'] or 'unknown'
                collision_summary['reference_versions'][rv] = collision_summary['reference_versions'].get(rv, 0) + 1
                
                # Count classification status
                if record['is_classified']:
                    collision_summary['classification_status']['classified'] += 1
                else:
                    collision_summary['classification_status']['unclassified'] += 1
                
                # Track domain counts
                collision_summary['domain_counts'].append(record['domain_count'])
        
        collision_data['collision_summary'] = collision_summary
        return collision_data
    
    def print_collision_report(self, collision_data: Dict):
        """Print a detailed collision report"""
        
        print("🚨 Mini PyECOD Collision Detection Report")
        print("=" * 60)
        print(f"Analysis Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        
        # Summary statistics
        summary = collision_data['collision_summary']
        total_mini = collision_data['total_mini_proteins']
        total_collisions = summary['total_collisions']
        
        print("📊 Summary:")
        print(f"  Mini proteins to import:    {total_mini:>8,}")
        print(f"  Potential collisions:       {total_collisions:>8,}")
        print(f"  Collision rate:             {total_collisions/max(1,total_mini)*100:>7.1f}%")
        print(f"  Safe imports:               {total_mini-total_collisions:>8,}")
        print()
        
        if total_collisions == 0:
            print("✅ No collisions detected! Safe to import all results.")
            return
        
        # Process version breakdown
        print("🔧 Existing Process Versions:")
        for pv, count in sorted(summary['process_versions'].items(), key=lambda x: -x[1]):
            print(f"  {pv:<25} {count:>8}")
        print()
        
        # Reference version breakdown
        print("📚 Existing Reference Versions:")
        for rv, count in sorted(summary['reference_versions'].items(), key=lambda x: -x[1]):
            print(f"  {rv:<25} {count:>8}")
        print()
        
        # Classification status
        classified = summary['classification_status']['classified']
        unclassified = summary['classification_status']['unclassified']
        print("📈 Existing Classification Status:")
        print(f"  Classified:             {classified:>8}")
        print(f"  Unclassified:           {unclassified:>8}")
        print(f"  Classification rate:    {classified/max(1,classified+unclassified)*100:>7.1f}%")
        print()
        
        # Domain count statistics
        if summary['domain_counts']:
            domain_counts = summary['domain_counts']
            avg_domains = sum(domain_counts) / len(domain_counts)
            max_domains = max(domain_counts)
            print("🧬 Existing Domain Statistics:")
            print(f"  Average domains/protein: {avg_domains:>7.1f}")
            print(f"  Maximum domains:         {max_domains:>8}")
            print(f"  Zero domains:            {domain_counts.count(0):>8}")
        print()
        
        # Recommendations
        print("💡 Recommendations:")
        
        if any('mini' in pv.lower() for pv in summary['process_versions']):
            print("  ⚠️  MINI RESULTS ALREADY EXIST - Check if previous import occurred")
        
        if total_collisions / total_mini > 0.5:
            print("  🔒 HIGH COLLISION RATE - Use 'separate' strategy to avoid overwrites")
        elif classified > unclassified:
            print("  📊 GOOD EXISTING DATA - Use 'skip' strategy to preserve existing results")
        else:
            print("  🔄 POOR EXISTING DATA - Consider 'update' strategy to improve results")
        
        print(f"  🎯 RECOMMENDED: Use --collision-strategy separate for safety")
        print()
        
        # Command suggestions
        print("🚀 Suggested Import Commands:")
        print("  # Safe import (creates separate records):")
        print("  python import_results.py --import-all --collision-strategy separate --limit 100")
        print()
        print("  # Skip existing records:")
        print("  python import_results.py --import-all --collision-strategy skip --limit 100")
        print()
        print("  # Just check collisions without importing:")
        print("  python import_results.py --check-collisions --limit 100")
    
    def print_summary_only(self, collision_data: Dict):
        """Print only the summary statistics"""
        
        summary = collision_data['collision_summary']
        total_mini = collision_data['total_mini_proteins']
        total_collisions = summary['total_collisions']
        
        print(f"Mini proteins: {total_mini:,} | Collisions: {total_collisions:,} | Rate: {total_collisions/max(1,total_mini)*100:.1f}%")
        
        if total_collisions > 0:
            print("Process versions:", ", ".join(f"{k}({v})" for k,v in summary['process_versions'].items()))


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Detect Potential Data Collisions for Mini PyECOD Import'
    )
    
    parser.add_argument('--batch-name', type=str,
                       help='Check specific batch only')
    parser.add_argument('--summary', action='store_true',
                       help='Show summary only (no detailed report)')
    parser.add_argument('--config', type=str, default='config.local.yml',
                       help='Config file path')
    
    args = parser.parse_args()
    
    # Initialize detector
    detector = CollisionDetector(args.config)
    
    # Scan mini results
    mini_proteins = detector.scan_mini_results(args.batch_name)
    
    if not mini_proteins:
        print("No mini results found to check.")
        return
    
    # Check for collisions
    collision_data = detector.check_database_collisions(mini_proteins)
    
    # Print report
    if args.summary:
        detector.print_summary_only(collision_data)
    else:
        detector.print_collision_report(collision_data)


if __name__ == "__main__":
    main()
