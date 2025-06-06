#!/usr/bin/env python3
"""
Reference Domain Recovery for Mini PyECOD Results

Recovers ECOD reference domain IDs for mini results by mapping
classifications back to known ECOD domains. Much cheaper than 
rerunning the entire pipeline.

Strategy:
1. Read existing mini XML results
2. Extract classification (t_group, h_group, etc.) for each domain
3. Query ECOD database for domains with matching classification  
4. Assign representative ECOD domain ID to each partition
5. Update results with reference domain links

Usage:
    python recover_reference_domains.py --batch-name batch_031
    python recover_reference_domains.py --scan-all
"""

import os
import sys
import argparse
import yaml
import xml.etree.ElementTree as ET
import psycopg2
import psycopg2.extras
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from collections import defaultdict
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class ClassificationLookup:
    """Represents a classification for domain lookup"""
    t_group: Optional[str]
    h_group: Optional[str] 
    x_group: Optional[str]
    a_group: Optional[str]
    
    def __hash__(self):
        return hash((self.t_group, self.h_group, self.x_group, self.a_group))
    
    def __eq__(self, other):
        return (self.t_group, self.h_group, self.x_group, self.a_group) == \
               (other.t_group, other.h_group, other.x_group, other.a_group)

@dataclass
class ReferenceDomain:
    """ECOD reference domain"""
    ecod_domain_id: str
    ecod_uid: int
    pdb_id: str
    chain_id: str
    classification: ClassificationLookup
    is_manual_rep: bool = False
    
class ReferenceDomainRecovery:
    """Recover reference domain IDs for mini results"""
    
    def __init__(self, config_path: str = "config.local.yml"):
        self.config = self._load_config(config_path)
        self.db_conn = self._init_db_connection()
        
        # Cache for classification -> reference domain mapping
        self.classification_cache = {}
        self.stats = {
            "processed": 0,
            "updated": 0,
            "no_classification": 0,
            "no_reference_found": 0,
            "cache_hits": 0
        }
    
    def _load_config(self, config_path: str) -> Dict:
        """Load configuration"""
        try:
            with open(config_path, 'r') as f:
                return yaml.safe_load(f)
        except FileNotFoundError:
            logger.error(f"Config file {config_path} not found")
            sys.exit(1)
    
    def _init_db_connection(self):
        """Initialize database connection to ECOD"""
        if not self.config.get('database'):
            logger.error("Database configuration required")
            sys.exit(1)
            
        try:
            return psycopg2.connect(**self.config['database'])
        except Exception as e:
            logger.error(f"Database connection failed: {e}")
            sys.exit(1)
    
    def build_classification_lookup(self) -> Dict[ClassificationLookup, List[ReferenceDomain]]:
        """Build lookup table: classification -> list of reference domains"""
        
        logger.info("Building classification lookup table from ECOD database...")
        
        with self.db_conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cursor:
            # Query from existing ECOD domain tables
            # Adjust table names based on your schema
            cursor.execute("""
                SELECT 
                    d.ecod_domain_id,
                    d.ecod_uid,
                    p.pdb_id,
                    p.chain_id,
                    eth.f_id as t_group,     -- T-group stored in f_id
                    eth.h_id as h_group,     -- H-group
                    eth.x_id as x_group,     -- X-group  
                    a.a_id as a_group,       -- A-group
                    d.is_manrep as is_manual_rep
                FROM 
                    public.domain d
                JOIN 
                    public.protein p ON d.ecod_source_id = p.source_id
                LEFT JOIN 
                    public.ecod_tmp_hier eth ON d.ecod_uid = eth.ecod_uid
                LEFT JOIN 
                    public.ecod_a_names a ON eth.a_id = a.a_id
                WHERE 
                    d.type = 'experimental structure'
                    AND eth.f_id IS NOT NULL  -- Must have classification
                ORDER BY 
                    d.is_manrep DESC,  -- Prefer manual representatives
                    d.ecod_uid ASC     -- Stable ordering
            """)
            
            results = cursor.fetchall()
            logger.info(f"Found {len(results)} classified ECOD domains")
        
        # Group by classification
        lookup = defaultdict(list)
        
        for row in results:
            classification = ClassificationLookup(
                t_group=row['t_group'],
                h_group=row['h_group'], 
                x_group=row['x_group'],
                a_group=row['a_group']
            )
            
            ref_domain = ReferenceDomain(
                ecod_domain_id=row['ecod_domain_id'],
                ecod_uid=row['ecod_uid'],
                pdb_id=row['pdb_id'],
                chain_id=row['chain_id'],
                classification=classification,
                is_manual_rep=row['is_manual_rep'] or False
            )
            
            lookup[classification].append(ref_domain)
        
        logger.info(f"Built lookup for {len(lookup)} unique classifications")
        return dict(lookup)
    
    def find_reference_domain(self, classification: ClassificationLookup) -> Optional[ReferenceDomain]:
        """Find best reference domain for a classification"""
        
        # Check cache first
        if classification in self.classification_cache:
            self.stats["cache_hits"] += 1
            return self.classification_cache[classification]
        
        # Query database if not cached
        if not hasattr(self, '_lookup_table'):
            self._lookup_table = self.build_classification_lookup()
        
        candidates = self._lookup_table.get(classification, [])
        
        if not candidates:
            self.classification_cache[classification] = None
            return None
        
        # Prefer manual representatives, then first by UID
        best_candidate = sorted(candidates, 
                               key=lambda x: (not x.is_manual_rep, x.ecod_uid))[0]
        
        self.classification_cache[classification] = best_candidate
        return best_candidate
    
    def process_mini_result(self, xml_file: Path) -> bool:
        """Process a single mini result file and add reference domains"""
        
        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()
            
            domain_elements = root.findall(".//domain")
            modified = False
            
            for domain_elem in domain_elements:
                self.stats["processed"] += 1
                
                # Extract current classification
                t_group = domain_elem.get('t_group')
                h_group = domain_elem.get('h_group')
                x_group = domain_elem.get('x_group')
                a_group = domain_elem.get('a_group')
                
                # Skip if no classification
                if not t_group:
                    self.stats["no_classification"] += 1
                    continue
                
                classification = ClassificationLookup(
                    t_group=t_group,
                    h_group=h_group,
                    x_group=x_group,
                    a_group=a_group
                )
                
                # Find reference domain
                ref_domain = self.find_reference_domain(classification)
                
                if not ref_domain:
                    self.stats["no_reference_found"] += 1
                    logger.warning(f"No reference domain found for classification: {classification}")
                    continue
                
                # Add reference domain information
                domain_elem.set('reference_ecod_domain_id', ref_domain.ecod_domain_id)
                domain_elem.set('reference_ecod_uid', str(ref_domain.ecod_uid))
                domain_elem.set('reference_pdb_id', ref_domain.pdb_id)
                domain_elem.set('reference_chain_id', ref_domain.chain_id)
                
                self.stats["updated"] += 1
                modified = True
            
            # Write back if modified
            if modified:
                # Create backup
                backup_path = xml_file.with_suffix('.xml.backup')
                if not backup_path.exists():
                    xml_file.rename(backup_path)
                    
                # Write updated XML
                ET.indent(tree, space="  ")
                tree.write(xml_file, encoding="utf-8", xml_declaration=True)
                
                return True
            
            return False
            
        except Exception as e:
            logger.error(f"Error processing {xml_file}: {e}")
            return False
    
    def process_batch(self, batch_name: str) -> Dict:
        """Process all mini results in a batch"""
        
        batch_base = Path(self.config["paths"]["batch_base_dir"])
        batch_dir = batch_base / batch_name
        mini_domains_dir = batch_dir / "mini_domains"
        
        if not mini_domains_dir.exists():
            logger.error(f"No mini results found in {batch_name}")
            return {"error": f"Batch {batch_name} not found"}
        
        xml_files = list(mini_domains_dir.glob("*.mini.domains.xml"))
        logger.info(f"Processing {len(xml_files)} results in {batch_name}")
        
        # Reset stats for this batch
        batch_stats = {
            "batch_name": batch_name,
            "total_files": len(xml_files),
            "files_modified": 0,
            "domains_processed": 0,
            "domains_updated": 0,
            "errors": []
        }
        
        initial_stats = self.stats.copy()
        
        for xml_file in xml_files:
            try:
                if self.process_mini_result(xml_file):
                    batch_stats["files_modified"] += 1
            except Exception as e:
                batch_stats["errors"].append(f"{xml_file.name}: {e}")
        
        # Calculate batch-specific stats
        batch_stats["domains_processed"] = self.stats["processed"] - initial_stats["processed"]
        batch_stats["domains_updated"] = self.stats["updated"] - initial_stats["updated"]
        
        return batch_stats
    
    def process_all_batches(self) -> List[Dict]:
        """Process all available mini results"""
        
        batch_base = Path(self.config["paths"]["batch_base_dir"])
        all_results = []
        
        for batch_dir in batch_base.iterdir():
            if not batch_dir.is_dir():
                continue
                
            mini_domains_dir = batch_dir / "mini_domains"
            if not mini_domains_dir.exists():
                continue
            
            batch_result = self.process_batch(batch_dir.name)
            all_results.append(batch_result)
            
            logger.info(f"Batch {batch_dir.name}: {batch_result['files_modified']} files modified")
        
        return all_results
    
    def print_summary(self, batch_results: List[Dict] = None):
        """Print recovery summary"""
        
        print("🔗 Reference Domain Recovery Summary")
        print("=" * 50)
        
        if batch_results:
            total_files = sum(r.get("total_files", 0) for r in batch_results)
            total_modified = sum(r.get("files_modified", 0) for r in batch_results)
            print(f"Batches processed:    {len(batch_results)}")
            print(f"Files processed:      {total_files:,}")
            print(f"Files modified:       {total_modified:,}")
        
        print(f"Domains processed:    {self.stats['processed']:,}")
        print(f"Domains updated:      {self.stats['updated']:,}")
        print(f"No classification:    {self.stats['no_classification']:,}")
        print(f"No reference found:   {self.stats['no_reference_found']:,}")
        print(f"Cache hit rate:       {self.stats['cache_hits'] / max(1, self.stats['processed']) * 100:.1f}%")
        
        success_rate = self.stats['updated'] / max(1, self.stats['processed']) * 100
        print(f"Success rate:         {success_rate:.1f}%")
        
        if success_rate > 95:
            print("✅ Reference recovery highly successful")
        elif success_rate > 80:
            print("⚠️  Reference recovery mostly successful")
        else:
            print("❌ Reference recovery needs investigation")


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Recover Reference Domain IDs for Mini PyECOD Results'
    )
    
    parser.add_argument('--scan-all', action='store_true',
                       help='Process all available results')
    parser.add_argument('--batch-name', type=str,
                       help='Process specific batch')
    parser.add_argument('--config', type=str, default='config.local.yml',
                       help='Config file path')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be done without modifying files')
    
    args = parser.parse_args()
    
    # Initialize recovery system
    recovery = ReferenceDomainRecovery(args.config)
    
    if args.dry_run:
        logger.info("DRY RUN MODE - no files will be modified")
        # TODO: Implement dry-run logic
        return
    
    # Process results
    if args.batch_name:
        result = recovery.process_batch(args.batch_name)
        recovery.print_summary([result])
    elif args.scan_all:
        results = recovery.process_all_batches()
        recovery.print_summary(results)
    else:
        parser.print_help()
        return


if __name__ == "__main__":
    main()
