#!/usr/bin/env python3
"""
Reconcile tracking database with actual filesystem results
Run this to fix the disconnect between tracking and reality
"""

import os
import sys
import yaml
import psycopg2
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List, Tuple

def load_config(config_path: str = "config/config.local.yml"):
    """Load database config"""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def count_domains_in_file(xml_path: str) -> int:
    """Count domains in a mini domains XML file"""
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        domains = root.findall(".//domain")
        return len(domains)
    except:
        return 0

def find_completed_results() -> List[Tuple[str, str, int]]:
    """Find all completed mini domain files and count domains"""
    
    batch_base = Path("/data/ecod/pdb_updates/batches")
    completed_results = []
    
    print("🔍 Scanning for completed results...")
    
    for batch_dir in batch_base.iterdir():
        if not batch_dir.is_dir():
            continue
            
        mini_domains_dir = batch_dir / "mini_domains"
        if not mini_domains_dir.exists():
            continue
            
        print(f"  Scanning {batch_dir.name}...")
        
        for result_file in mini_domains_dir.glob("*.mini.domains.xml"):
            protein_id = result_file.stem.replace('.mini.domains', '')
            domain_count = count_domains_in_file(str(result_file))
            completed_results.append((protein_id, str(result_file), domain_count))
    
    print(f"✓ Found {len(completed_results)} completed results")
    return completed_results

def reconcile_tracking_database():
    """Update tracking database to match filesystem reality"""
    
    config = load_config()
    
    # Get completed results from filesystem
    completed_results = find_completed_results()
    
    if not completed_results:
        print("No completed results found!")
        return
    
    # Connect to database
    conn = psycopg2.connect(**config['database'])
    conn.autocommit = True
    
    print(f"📊 Reconciling {len(completed_results)} completed results...")
    
    updated_count = 0
    inserted_count = 0
    
    with conn.cursor() as cursor:
        for protein_id, result_path, domain_count in completed_results:
            # Check if tracking record exists
            cursor.execute("""
                SELECT status FROM mini_production.tracking_status 
                WHERE protein_id = %s
            """, (protein_id,))
            
            existing = cursor.fetchone()
            
            if existing:
                if existing[0] != 'completed':
                    # Update existing record to completed
                    cursor.execute("""
                        UPDATE mini_production.tracking_status 
                        SET status = 'completed',
                            domains_found = %s,
                            completed_at = NOW()
                        WHERE protein_id = %s
                    """, (domain_count, protein_id))
                    updated_count += 1
            else:
                # Insert new record
                cursor.execute("""
                    INSERT INTO mini_production.tracking_status 
                    (protein_id, status, domains_found, submitted_at, completed_at)
                    VALUES (%s, 'completed', %s, NOW(), NOW())
                """, (protein_id, domain_count))
                inserted_count += 1
    
    print(f"✅ Reconciliation complete:")
    print(f"   Updated existing records: {updated_count}")
    print(f"   Inserted new records: {inserted_count}")
    print(f"   Total completed: {updated_count + inserted_count}")
    
    # Show new statistics
    with conn.cursor() as cursor:
        cursor.execute("""
            SELECT status, COUNT(*) 
            FROM mini_production.tracking_status 
            GROUP BY status
        """)
        
        print(f"\n📊 Updated tracking statistics:")
        for status, count in cursor.fetchall():
            print(f"   {status}: {count:,}")
    
    conn.close()

def cleanup_duplicate_submissions():
    """Clean up duplicate submissions that haven't started yet"""
    
    config = load_config()
    conn = psycopg2.connect(**config['database'])
    conn.autocommit = True
    
    print("🧹 Cleaning up duplicate submissions...")
    
    with conn.cursor() as cursor:
        # Find proteins with multiple submissions
        cursor.execute("""
            SELECT protein_id, COUNT(*) as submission_count
            FROM mini_production.tracking_status 
            GROUP BY protein_id
            HAVING COUNT(*) > 1
        """)
        
        duplicates = cursor.fetchall()
        print(f"Found {len(duplicates)} proteins with multiple submissions")
        
        if duplicates:
            # Keep only the most recent submission for each protein
            cursor.execute("""
                DELETE FROM mini_production.tracking_status 
                WHERE ctid NOT IN (
                    SELECT DISTINCT ON (protein_id) ctid
                    FROM mini_production.tracking_status 
                    ORDER BY protein_id, submitted_at DESC
                )
            """)
            
            deleted_count = cursor.rowcount
            print(f"✓ Deleted {deleted_count} duplicate submission records")
    
    conn.close()

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Reconcile tracking with filesystem results')
    parser.add_argument('--reconcile', action='store_true', help='Update tracking database')
    parser.add_argument('--cleanup-duplicates', action='store_true', help='Remove duplicate submissions')
    parser.add_argument('--show-disconnect', action='store_true', help='Show the disconnect only')
    
    args = parser.parse_args()
    
    if args.show_disconnect:
        config = load_config()
        conn = psycopg2.connect(**config['database'])
        
        with conn.cursor() as cursor:
            cursor.execute("SELECT COUNT(*) FROM mini_production.tracking_status WHERE status = 'completed'")
            tracking_completed = cursor.fetchone()[0]
            
            cursor.execute("SELECT COUNT(*) FROM mini_production.tracking_status")
            total_tracking = cursor.fetchone()[0]
        
        completed_results = find_completed_results()
        
        print(f"📊 Disconnect Analysis:")
        print(f"   Actual files on disk: {len(completed_results):,}")
        print(f"   Tracking says completed: {tracking_completed:,}")
        print(f"   Total tracking records: {total_tracking:,}")
        print(f"   Disconnect: {len(completed_results) - tracking_completed:,}")
        
        conn.close()
        
    elif args.reconcile:
        reconcile_tracking_database()
        
    elif args.cleanup_duplicates:
        cleanup_duplicate_submissions()
        
    else:
        print("Use --show-disconnect, --reconcile, or --cleanup-duplicates")
