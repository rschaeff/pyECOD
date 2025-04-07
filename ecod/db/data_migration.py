# ecod/db/data_migration.py
import logging
from typing import Dict, Any
import psycopg2
import psycopg2.extras

class DataMigration:
    def __init__(self, source_config: Dict[str, Any], target_config: Dict[str, Any]):
        self.source_config = source_config
        self.target_config = target_config
        self.logger = logging.getLogger("ecod.data_migration")
        
    def migrate_protein_data(self):
        """Migrate protein and sequence data from old to new schema"""
        source_conn = psycopg2.connect(**self.source_config)
        target_conn = psycopg2.connect(**self.target_config)
        
        try:
            # Get proteins from old schema
            self.logger.info("Fetching proteins from old schema")
            with source_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as source_cursor:
                source_cursor.execute("""
                    SELECT 
                        id, pdb_id, chain_id, source_id, length
                    FROM 
                        pdb_analysis.protein
                """)
                proteins = source_cursor.fetchall()
                
            self.logger.info(f"Found {len(proteins)} proteins")
            
            # Insert into new schema
            with target_conn.cursor() as target_cursor:
                for protein in proteins:
                    target_cursor.execute("""
                        INSERT INTO ecod_schema.protein
                        (id, pdb_id, chain_id, source_id, length)
                        VALUES (%s, %s, %s, %s, %s)
                        ON CONFLICT (source_id) DO NOTHING
                    """, (
                        protein['id'],
                        protein['pdb_id'],
                        protein['chain_id'],
                        protein['source_id'],
                        protein['length']
                    ))
                    
            target_conn.commit()
            self.logger.info("Proteins migrated successfully")
            
            # Migrate sequences
            self.logger.info("Fetching protein sequences")
            with source_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as source_cursor:
                source_cursor.execute("""
                    SELECT 
                        protein_id, sequence, md5(sequence) as md5_hash
                    FROM 
                        pdb_analysis.protein_sequence
                """)
                sequences = source_cursor.fetchall()
                
            self.logger.info(f"Found {len(sequences)} sequences")
            
            # Insert into new schema
            with target_conn.cursor() as target_cursor:
                for sequence in sequences:
                    target_cursor.execute("""
                        INSERT INTO ecod_schema.protein_sequence
                        (protein_id, sequence, md5_hash)
                        VALUES (%s, %s, %s)
                        ON CONFLICT (md5_hash) DO NOTHING
                    """, (
                        sequence['protein_id'],
                        sequence['sequence'],
                        sequence['md5_hash']
                    ))
                    
            target_conn.commit()
            self.logger.info("Sequences migrated successfully")
            
        finally:
            source_conn.close()
            target_conn.close()
            
    def migrate_reference_data(self):
        """Migrate ECOD version and reference data"""
        source_conn = psycopg2.connect(**self.source_config)
        target_conn = psycopg2.connect(**self.target_config)
        
        try:
            # Get versions from old schema
            self.logger.info("Fetching ECOD versions")
            with source_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as source_cursor:
                source_cursor.execute("""
                    SELECT 
                        id, version_name, release_date, is_current
                    FROM 
                        pdb_analysis.ecod_versions
                """)
                versions = source_cursor.fetchall()
                
            self.logger.info(f"Found {len(versions)} ECOD versions")
            
            # Insert into new schema
            with target_conn.cursor() as target_cursor:
                for version in versions:
                    target_cursor.execute("""
                        INSERT INTO ecod_schema.ecod_version
                        (id, version_name, release_date, is_current)
                        VALUES (%s, %s, %s, %s)
                        ON CONFLICT (version_name) DO NOTHING
                    """, (
                        version['id'],
                        version['version_name'],
                        version['release_date'],
                        version['is_current']
                    ))
                    
            target_conn.commit()
            self.logger.info("ECOD versions migrated successfully")
            
            # Get reference resources
            self.logger.info("Fetching reference resources")
            with source_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as source_cursor:
                source_cursor.execute("""
                    SELECT 
                        id, version_id, resource_type, resource_path
                    FROM 
                        pdb_analysis.ecod_reference_resources
                """)
                resources = source_cursor.fetchall()
                
            self.logger.info(f"Found {len(resources)} reference resources")
            
            # Insert into new schema
            with target_conn.cursor() as target_cursor:
                for resource in resources:
                    target_cursor.execute("""
                        INSERT INTO ecod_schema.reference_resource
                        (id, version_id, resource_type, resource_path)
                        VALUES (%s, %s, %s, %s)
                        ON CONFLICT (version_id, resource_type) DO NOTHING
                    """, (
                        resource['id'],
                        resource['version_id'],
                        resource['resource_type'],
                        resource['resource_path']
                    ))
                    
            target_conn.commit()
            self.logger.info("Reference resources migrated successfully")
            
        finally:
            source_conn.close()
            target_conn.close()