#!/usr/bin/env python3
"""
Import Mini PyECOD Results to Production Database

Scans completed mini domain files and imports them to:
- pdb_analysis.partition_proteins (protein metadata)  
- pdb_analysis.partition_domains (individual domains)

Usage:
    python import_results.py --scan-completed                    # Find completed results
    python import_results.py --batch-name batch_036             # Import specific batch  
    python import_results.py --import-all --limit 100           # Import first 100 results
    python import_results.py --verify                           # Verify imported data
"""

import os
import sys
import argparse
import yaml
import psycopg2
import psycopg2.extras
import xml.etree.ElementTree as ET
import sqlite3
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
from datetime import datetime
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class MiniBatch:
    """Information about a batch with mini results"""
    batch_name: str
    batch_dir: Path
    mini_domains_dir: Path
    completed_files: List[Path]
    database_batch_id: Optional[int] = None

@dataclass
class MiniDomain:
    """Parsed domain from mini XML"""
    domain_id: str
    start_pos: int
    end_pos: int
    range_str: str
    source: str = "mini_pyecod"
    source_id: str = ""
    confidence: float = 0.8  # Default mini confidence
    t_group: Optional[str] = None
    h_group: Optional[str] = None
    x_group: Optional[str] = None
    a_group: Optional[str] = None
    is_manual_rep: bool = False
    is_f70: bool = False
    is_f40: bool = False
    is_f99: bool = False

@dataclass
class MiniProteinResult:
    """Complete mini result for a protein"""
    protein_id: str  # e.g., "8ovp_A"
    pdb_id: str      # e.g., "8ovp"
    chain_id: str    # e.g., "A"
    batch_name: str
    xml_file: Path
    is_classified: bool
    domains: List[MiniDomain]
    sequence_length: int = 0
    reference_version: str = "mini_develop291"

class MiniResultsImporter:
    """Import mini PyECOD results to production database"""
    
    def __init__(self, config_path: str = "config.local.yml"):
        self.config = self._load_config(config_path)
        self.db_conn = self._init_db_connection()
        self.tracking_db = self._init_tracking_db()
        
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
    
    def _init_tracking_db(self):
        """Initialize tracking database"""
        tracking_path = self.config.get('paths', {}).get('tracking_db', '/tmp/mini_production_status.db')
        return sqlite3.connect(tracking_path)
    
    def scan_completed_results(self) -> List[MiniBatch]:
        """Scan for completed mini results across all batches"""
        
        batch_base = Path(self.config["paths"]["batch_base_dir"])
        logger.info(f"🔍 Scanning for completed results in {batch_base}")
        
        mini_batches = []
        
        for batch_dir in batch_base.iterdir():
            if not batch_dir.is_dir():
                continue
                
            mini_domains_dir = batch_dir / "mini_domains"
            if not mini_domains_dir.exists():
                continue
            
            # Find completed XML files
            xml_files = list(mini_domains_dir.glob("*.mini.domains.xml"))
            
            if xml_files:
                mini_batch = MiniBatch(
                    batch_name=batch_dir.name,
                    batch_dir=batch_dir,
                    mini_domains_dir=mini_domains_dir,
                    completed_files=xml_files
                )
                mini_batches.append(mini_batch)
                logger.info(f"📁 {batch_dir.name}: {len(xml_files)} completed")
        
        total_files = sum(len(batch.completed_files) for batch in mini_batches)
        logger.info(f"✓ Found {total_files} completed results across {len(mini_batches)} batches")
        
        return mini_batches
    
    def parse_mini_xml(self, xml_file: Path) -> MiniProteinResult:
        """Parse a mini domain XML file"""
        
        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()
            
            # Extract protein info from filename: "8ovp_A.mini.domains.xml"
            protein_id = xml_file.stem.replace('.mini.domains', '')
            pdb_id = protein_id.split('_')[0]
            chain_id = protein_id.split('_')[1] if '_' in protein_id else 'A'
            
            # Get batch name from directory structure
            batch_name = xml_file.parent.parent.name
            
            # Parse domains
            domains = []
            domain_elements = root.findall(".//domain")
            
            for i, domain_elem in enumerate(domain_elements):
                # Extract domain attributes
                range_str = domain_elem.get('range', '')
                start_pos = int(domain_elem.get('start', '0'))
                end_pos = int(domain_elem.get('end', '0'))
                
                domain = MiniDomain(
                    domain_id=f"{protein_id}_d{i+1}",
                    start_pos=start_pos,
                    end_pos=end_pos,
                    range_str=range_str,
                    source="mini_pyecod",
                    source_id=domain_elem.get('source_id', ''),
                    confidence=float(domain_elem.get('confidence', '0.8')),
                    t_group=domain_elem.get('t_group'),
                    h_group=domain_elem.get('h_group'),
                    x_group=domain_elem.get('x_group'),
                    a_group=domain_elem.get('a_group'),
                    is_manual_rep=domain_elem.get('is_manual_rep', 'false').lower() == 'true',
                    is_f70=domain_elem.get('is_f70', 'false').lower() == 'true',
                    is_f40=domain_elem.get('is_f40', 'false').lower() == 'true',
                    is_f99=domain_elem.get('is_f99', 'false').lower() == 'true'
                )
                domains.append(domain)
            
            # Create result
            result = MiniProteinResult(
                protein_id=protein_id,
                pdb_id=pdb_id,
                chain_id=chain_id,
                batch_name=batch_name,
                xml_file=xml_file,
                is_classified=len(domains) > 0,
                domains=domains,
                sequence_length=0,  # Will get from database if needed
                reference_version="mini_develop291"
            )
            
            return result
            
        except Exception as e:
            logger.error(f"Error parsing {xml_file}: {e}")
            raise
    
    def get_or_create_batch_mapping(self, batch_name: str) -> Optional[int]:
        """Get database batch_id for a batch name, or None if not found"""
        
        try:
            with self.db_conn.cursor() as cursor:
                # Try to find existing batch in ecod_schema
                cursor.execute("""
                    SELECT id FROM ecod_schema.batch 
                    WHERE batch_name = %s
                    ORDER BY id DESC LIMIT 1
                """, (batch_name,))
                
                result = cursor.fetchone()
                if result:
                    return result[0]
                
                # If not found, return None - we'll handle this case
                logger.warning(f"No database batch_id found for batch: {batch_name}")
                return None
                
        except Exception as e:
            logger.error(f"Error getting batch mapping for {batch_name}: {e}")
            return None
    
    def import_protein_result(self, result: MiniProteinResult, collision_strategy: str = "separate") -> bool:
        """Import a single protein result to database with collision handling
        
        Args:
            result: Protein result to import
            collision_strategy: How to handle existing data
                - "separate": Use distinct process_version (recommended)
                - "skip": Skip if existing record found
                - "update": Update existing record (dangerous!)
                - "check": Check for conflicts and report only
        """
        
        try:
            with self.db_conn.cursor() as cursor:
                # Get batch_id
                batch_id = self.get_or_create_batch_mapping(result.batch_name)
                
                # Check for existing records first
                cursor.execute("""
                    SELECT id, process_version, is_classified, domains_with_evidence,
                           timestamp
                    FROM pdb_analysis.partition_proteins 
                    WHERE pdb_id = %s AND chain_id = %s AND batch_id = %s
                    ORDER BY timestamp DESC
                """, (result.pdb_id, result.chain_id, batch_id))
                
                existing_records = cursor.fetchall()
                
                # Handle collision strategy
                if existing_records and collision_strategy == "skip":
                    logger.info(f"⏭️  Skipping {result.protein_id}: existing record found")
                    return True
                
                elif existing_records and collision_strategy == "check":
                    for record in existing_records:
                        logger.warning(f"🔍 COLLISION: {result.protein_id} exists with process_version='{record[1]}', "
                                     f"classified={record[2]}, domains={record[3]}, timestamp={record[4]}")
                    return True
                
                elif existing_records and collision_strategy == "separate":
                    # Use distinct process_version to avoid collision
                    mini_process_version = "mini_pyecod_1.0"
                    logger.info(f"🔀 Separate import: {result.protein_id} (process_version={mini_process_version})")
                    
                elif existing_records and collision_strategy == "update":
                    logger.warning(f"⚠️  Updating existing record for {result.protein_id}")
                
                # Determine process version based on strategy
                if collision_strategy == "separate":
                    process_version = "mini_pyecod_1.0"
                    reference_version = "mini_develop291"
                else:
                    process_version = "mini_pyecod_1.0"  # Still distinguish mini results
                    reference_version = result.reference_version
                
                # Insert protein record
                if collision_strategy == "update" and existing_records:
                    # Update existing record
                    cursor.execute("""
                        UPDATE pdb_analysis.partition_proteins 
                        SET is_classified = %s,
                            domains_with_evidence = %s,
                            fully_classified_domains = %s,
                            process_version = %s,
                            reference_version = %s,
                            residues_assigned = %s
                        WHERE pdb_id = %s AND chain_id = %s AND batch_id = %s
                        RETURNING id
                    """, (
                        result.is_classified,
                        len(result.domains),
                        len([d for d in result.domains if d.t_group]),
                        process_version,
                        reference_version,
                        sum(d.end_pos - d.start_pos + 1 for d in result.domains),
                        result.pdb_id,
                        result.chain_id,
                        batch_id
                    ))
                    protein_db_id = cursor.fetchone()[0]
                    
                    # Delete existing domains for this protein
                    cursor.execute("""
                        DELETE FROM pdb_analysis.partition_domains 
                        WHERE protein_id = %s
                    """, (protein_db_id,))
                
                else:
                    # Insert new record (separate strategy creates new record even if collision exists)
                    cursor.execute("""
                        INSERT INTO pdb_analysis.partition_proteins 
                        (pdb_id, chain_id, batch_id, reference_version, is_classified, 
                         sequence_length, coverage, residues_assigned, domains_with_evidence, 
                         fully_classified_domains, process_version)
                        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                        RETURNING id
                    """, (
                        result.pdb_id,
                        result.chain_id,
                        batch_id,
                        reference_version,
                        result.is_classified,
                        result.sequence_length,
                        0.0,  # coverage - calculate later if needed
                        sum(d.end_pos - d.start_pos + 1 for d in result.domains),
                        len(result.domains),
                        len([d for d in result.domains if d.t_group]),
                        process_version
                    ))
                    protein_db_id = cursor.fetchone()[0]
                
                # Insert domain records
                for i, domain in enumerate(result.domains):
                    cursor.execute("""
                        INSERT INTO pdb_analysis.partition_domains
                        (protein_id, domain_number, domain_id, start_pos, end_pos, range,
                         source, source_id, confidence, t_group, h_group, x_group, a_group,
                         is_manual_rep, is_f70, is_f40, is_f99)
                        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                    """, (
                        protein_db_id,
                        i + 1,  # domain_number
                        domain.domain_id,
                        domain.start_pos,
                        domain.end_pos,
                        domain.range_str,
                        domain.source,
                        domain.source_id,
                        domain.confidence,
                        domain.t_group,
                        domain.h_group,
                        domain.x_group,
                        domain.a_group,
                        domain.is_manual_rep,
                        domain.is_f70,
                        domain.is_f40,
                        domain.is_f99
                    ))
                
                self.db_conn.commit()
                
                strategy_msg = f"({collision_strategy})" if existing_records else ""
                logger.info(f"✓ Imported {result.protein_id}: {len(result.domains)} domains {strategy_msg}")
                return True
                
        except Exception as e:
            logger.error(f"Error importing {result.protein_id}: {e}")
            self.db_conn.rollback()
            return False
    
    def import_batch_results(self, batch_name: str, limit: Optional[int] = None, 
                           collision_strategy: str = "separate") -> Dict[str, int]:
        """Import all results from a specific batch
        
        Args:
            batch_name: Name of batch to import
            limit: Maximum number of results to import
            collision_strategy: How to handle existing data (separate/skip/update/check)
        """
        
        logger.info(f"🚀 Importing results from batch: {batch_name} (strategy: {collision_strategy})")
        
        # Find the batch
        mini_batches = self.scan_completed_results()
        target_batch = None
        
        for batch in mini_batches:
            if batch.batch_name == batch_name:
                target_batch = batch
                break
        
        if not target_batch:
            logger.error(f"Batch not found: {batch_name}")
            return {'imported': 0, 'failed': 0, 'total': 0, 'skipped': 0, 'collisions': 0}
        
        # Process files
        files_to_process = target_batch.completed_files[:limit] if limit else target_batch.completed_files
        
        imported = 0
        failed = 0
        skipped = 0
        collisions = 0
        
        for xml_file in files_to_process:
            try:
                result = self.parse_mini_xml(xml_file)
                
                # Count collisions if checking
                if collision_strategy == "check":
                    # This will log collisions but not import
                    success = self.import_protein_result(result, collision_strategy)
                    if success:
                        collisions += 1
                else:
                    success = self.import_protein_result(result, collision_strategy)
                    if success:
                        imported += 1
                    else:
                        failed += 1
                        
            except Exception as e:
                logger.error(f"Failed to process {xml_file}: {e}")
                failed += 1
        
        stats = {
            'imported': imported,
            'failed': failed,
            'total': len(files_to_process),
            'skipped': skipped,
            'collisions': collisions if collision_strategy == "check" else 0
        }
        
        if collision_strategy == "check":
            logger.info(f"✓ Collision check {batch_name}: {collisions} potential collisions found")
        else:
            logger.info(f"✓ Batch {batch_name}: {imported}/{len(files_to_process)} imported")
        return stats
    
    def import_all_results(self, limit: Optional[int] = None, 
                         collision_strategy: str = "separate") -> Dict[str, int]:
        """Import all available mini results
        
        Args:
            limit: Maximum number of results to import across all batches
            collision_strategy: How to handle existing data (separate/skip/update/check)
        """
        
        logger.info(f"🚀 Importing all mini results (strategy: {collision_strategy})")
        
        mini_batches = self.scan_completed_results()
        
        total_imported = 0
        total_failed = 0
        total_processed = 0
        total_collisions = 0
        
        for batch in mini_batches:
            logger.info(f"Processing batch: {batch.batch_name}")
            
            # Apply limit across all batches
            remaining_limit = limit - total_processed if limit else None
            if remaining_limit is not None and remaining_limit <= 0:
                break
            
            batch_limit = min(len(batch.completed_files), remaining_limit) if remaining_limit else None
            
            stats = self.import_batch_results(batch.batch_name, batch_limit, collision_strategy)
            
            total_imported += stats['imported']
            total_failed += stats['failed']
            total_processed += stats['total']
            total_collisions += stats.get('collisions', 0)
        
        final_stats = {
            'imported': total_imported,
            'failed': total_failed,
            'total': total_processed,
            'batches_processed': len(mini_batches),
            'collisions': total_collisions
        }
        
        if collision_strategy == "check":
            logger.info(f"🔍 Collision check complete: {total_collisions} potential collisions across {len(mini_batches)} batches")
        else:
            logger.info(f"🎉 Import complete: {total_imported}/{total_processed} imported across {len(mini_batches)} batches")
        return final_stats
    
    def verify_import(self, sample_size: int = 10) -> Dict[str, any]:
        """Verify imported data quality"""
        
        logger.info(f"🔍 Verifying import quality (sample size: {sample_size})")
        
        with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            # Get mini import statistics
            cursor.execute("""
                SELECT 
                    COUNT(*) as total_proteins,
                    COUNT(*) FILTER (WHERE is_classified = true) as classified_proteins,
                    AVG(domains_with_evidence) as avg_domains_per_protein,
                    COUNT(DISTINCT reference_version) as versions_count
                FROM pdb_analysis.partition_proteins 
                WHERE process_version = 'mini_pyecod_1.0'
            """)
            
            protein_stats = cursor.fetchone()
            
            # Get domain statistics
            cursor.execute("""
                SELECT 
                    COUNT(*) as total_domains,
                    COUNT(*) FILTER (WHERE t_group IS NOT NULL) as classified_domains,
                    COUNT(DISTINCT t_group) as unique_t_groups,
                    AVG(confidence) as avg_confidence
                FROM pdb_analysis.partition_domains pd
                JOIN pdb_analysis.partition_proteins pp ON pd.protein_id = pp.id
                WHERE pp.process_version = 'mini_pyecod_1.0'
            """)
            
            domain_stats = cursor.fetchone()
            
            # Sample some results for detailed check
            cursor.execute("""
                SELECT pp.pdb_id, pp.chain_id, pp.domains_with_evidence,
                       pd.domain_id, pd.t_group, pd.confidence
                FROM pdb_analysis.partition_proteins pp
                JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
                WHERE pp.process_version = 'mini_pyecod_1.0'
                ORDER BY pp.id DESC 
                LIMIT %s
            """, (sample_size,))
            
            samples = cursor.fetchall()
        
        verification = {
            'protein_stats': dict(protein_stats),
            'domain_stats': dict(domain_stats),
            'samples': [dict(s) for s in samples],
            'quality_score': min(100, protein_stats['classified_proteins'] / max(1, protein_stats['total_proteins']) * 100)
        }
        
        logger.info(f"✓ Verification complete:")
        logger.info(f"  Proteins: {protein_stats['total_proteins']} ({protein_stats['classified_proteins']} classified)")
        logger.info(f"  Domains: {domain_stats['total_domains']} ({domain_stats['classified_domains']} classified)")
        logger.info(f"  Quality Score: {verification['quality_score']:.1f}%")
        
        return verification


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Import Mini PyECOD Results to Production Database'
    )
    
    parser.add_argument('--scan-completed', action='store_true',
                       help='Scan and show completed results')
    parser.add_argument('--batch-name', type=str,
                       help='Import specific batch')
    parser.add_argument('--import-all', action='store_true',
                       help='Import all available results')
    parser.add_argument('--limit', type=int,
                       help='Limit number of results to import')
    parser.add_argument('--verify', action='store_true',
                       help='Verify imported data')
    parser.add_argument('--check-collisions', action='store_true',
                       help='Check for potential data collisions without importing')
    parser.add_argument('--collision-strategy', type=str, 
                       choices=['separate', 'skip', 'update', 'check'], 
                       default='separate',
                       help='How to handle existing data (default: separate)')
    parser.add_argument('--config', type=str, default='config.local.yml',
                       help='Config file path')
    
    args = parser.parse_args()
    
    # Initialize importer
    importer = MiniResultsImporter(args.config)
    
    if args.scan_completed:
        mini_batches = importer.scan_completed_results()
        print(f"\n📊 Completed Results Summary:")
        print(f"{'Batch Name':<40} {'Files':<10}")
        print("-" * 55)
        for batch in mini_batches:
            print(f"{batch.batch_name:<40} {len(batch.completed_files):<10}")
        
        total = sum(len(b.completed_files) for b in mini_batches)
        print(f"\nTotal: {total} completed results across {len(mini_batches)} batches")
        return
    
    if args.verify:
        verification = importer.verify_import()
        return
    
    if args.check_collisions:
        args.collision_strategy = 'check'
        print("🔍 Checking for data collisions...")
    
    # Set collision strategy
    collision_strategy = args.collision_strategy
    
    if args.batch_name:
        stats = importer.import_batch_results(args.batch_name, args.limit, collision_strategy)
        print(f"\n✅ Import Results:")
        print(f"  Imported: {stats['imported']}")
        print(f"  Failed: {stats['failed']}")
        print(f"  Total: {stats['total']}")
        if collision_strategy == 'check':
            print(f"  Collisions: {stats['collisions']}")
        return
    
    if args.import_all or args.check_collisions:
        stats = importer.import_all_results(args.limit, collision_strategy)
        print(f"\n🎉 Final Results:")
        print(f"  Imported: {stats['imported']}")
        print(f"  Failed: {stats['failed']}")
        print(f"  Total: {stats['total']}")
        print(f"  Batches: {stats['batches_processed']}")
        if collision_strategy == 'check':
            print(f"  Collisions: {stats['collisions']}")
        return
    
    # Default: show help
    parser.print_help()


if __name__ == "__main__":
    main()
