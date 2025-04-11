"""
Database management commands for the ECOD pipeline
"""

import argparse
import logging
import os
from typing import Dict, Any

from base_command import BaseCommand
from ecod.db.migration_manager import MigrationManager

logger = logging.getLogger("ecod.cli.db")

# Define commands in this group
COMMANDS = {
    'migrate': 'Apply database migrations',
    'import': 'Import data from external sources',
    'sample': 'Insert sample data for testing',
    'status': 'Check database status'
}

def setup_parser(parser: argparse.ArgumentParser) -> None:
    """Set up the argument parser for database commands"""
    subparsers = parser.add_subparsers(dest='command', help='Database command')
    
    # Migrate command
    migrate_parser = subparsers.add_parser('migrate', help=COMMANDS['migrate'])
    migrate_parser.add_argument('--migrations-dir', type=str, 
                             default='ecod/db/migrations',
                             help='Directory containing migration files')
    
    # Import command
    import_parser = subparsers.add_parser('import', help=COMMANDS['import'])
    import_parser.add_argument('--source', type=str, required=True,
                            help='Source type (pdb, uniprot, ecod)')
    import_parser.add_argument('--file', type=str,
                            help='File to import from')
    import_parser.add_argument('--limit', type=int,
                            help='Maximum items to import')
    
    # Sample command
    sample_parser = subparsers.add_parser('sample', help=COMMANDS['sample'])
    sample_parser.add_argument('--pdb-id', type=str, default='1abc',
                           help='PDB ID for sample protein')
    sample_parser.add_argument('--chain-id', type=str, default='A',
                           help='Chain ID for sample protein')
    sample_parser.add_argument('--sequence', type=str,
                           help='Amino acid sequence (uses default if not provided)')
    sample_parser.add_argument('--count', type=int, default=1,
                           help='Number of sample proteins to create')
    
    # Status command
    status_parser = subparsers.add_parser('status', help=COMMANDS['status'])
    status_parser.add_argument('--schema', type=str,
                            help='Schema to check (defaults to all)')

def run_command(args: argparse.Namespace) -> int:
    """Run the specified database command"""
    # Load configuration
    config_manager = ConfigManager(args.config)
    db_config = config_manager.get_db_config()
    
    # Handle different commands
    if args.command == 'migrate':
        return _run_migrations(args, db_config)
    elif args.command == 'import':
        return _import_data(args, db_config, config_manager)
    elif args.command == 'sample':
        return _insert_sample_data(args, db_config)
    elif args.command == 'status':
        return _check_status(args, db_config)
    else:
        logger.error(f"Unknown command: {args.command}")
        return 1

def _run_migrations(args: argparse.Namespace, db_config: Dict[str, Any]) -> int:
    """Run database migrations"""
    try:
        logger.info(f"Running migrations from {args.migrations_dir}")
        migration_manager = MigrationManager(db_config, args.migrations_dir)
        migration_manager.apply_migrations()
        logger.info("Migrations applied successfully")
        return 0
    except Exception as e:
        logger.error(f"Error applying migrations: {str(e)}")
        return 1

def _import_data(args: argparse.Namespace, db_config: Dict[str, Any], config_manager: ConfigManager) -> int:
    """Import data from external sources"""
    db = DBManager(db_config)
    
    if args.source.lower() == 'pdb':
        return _import_pdb_data(args, db, config_manager)
    elif args.source.lower() == 'uniprot':
        return _import_uniprot_data(args, db, config_manager)
    elif args.source.lower() == 'ecod':
        return _import_ecod_data(args, db, config_manager)
    else:
        logger.error(f"Unknown source type: {args.source}")
        return 1

def _import_pdb_data(args: argparse.Namespace, db: DBManager, config_manager: ConfigManager) -> int:
    """Import data from PDB"""
    # Import data migrator
    from ecod.db.data_migration import DataMigration
    
    # Configure source database
    source_config = db.config.copy()
    source_config['schema'] = 'pdb_analysis'
    
    # Create data migrator
    data_migrator = DataMigration(source_config, db.config)
    
    try:
        # Migrate protein data
        data_migrator.migrate_protein_data()
        logger.info("Successfully imported PDB protein data")
        return 0
    except Exception as e:
        logger.error(f"Error importing PDB data: {str(e)}")
        return 1

def _import_uniprot_data(args: argparse.Namespace, db: DBManager, config_manager: ConfigManager) -> int:
    """Import data from UniProt"""
    logger.error("UniProt import not implemented yet")
    return 1

def _import_ecod_data(args: argparse.Namespace, db: DBManager, config_manager: ConfigManager) -> int:
    """Import data from ECOD reference database"""
    # Import data migrator
    from ecod.db.data_migration import DataMigration
    
    # Configure source database
    source_config = db.config.copy()
    source_config['schema'] = 'public'
    
    # Create data migrator
    data_migrator = DataMigration(source_config, db.config)
    
    try:
        # Migrate reference data
        data_migrator.migrate_reference_data()
        logger.info("Successfully imported ECOD reference data")
        return 0
    except Exception as e:
        logger.error(f"Error importing ECOD data: {str(e)}")
        return 1

def _insert_sample_data(args: argparse.Namespace, db_config: Dict[str, Any]) -> int:
    """Insert sample data for testing"""
    db = DBManager(db_config)
    
    # Default sample sequence if not provided
    sequence = args.sequence
    if not sequence:
        sequence = "MKKRLTESQFQEGTARIMSSFGLSKKGVANITVAVTNGYSGKQLISQQEDITKSLLVTELGSSRTPETVRMVLSNMEKLSSADFVFLTADQVEEKILSKHQGVQVLESTLMTIPDTPASMLLKLTDVPAKLKEALDATLNDKEFLERKAQEKEAFLEALKTMTKDLPRAQKAAFEEQEKQWLEAAQVAPA"
    
    success = True
    for i in range(args.count):
        try:
            # For multiple proteins, append a number to the chain ID
            chain_id = args.chain_id if args.count == 1 else f"{args.chain_id}{i+1}"
            
            # Create the source_id (format used in the code is pdb_id_chain_id)
            source_id = f"{args.pdb_id}_{chain_id}"
            
            # Check if protein already exists
            check_query = "SELECT id FROM ecod_schema.protein WHERE source_id = %s"
            existing = db.execute_query(check_query, (source_id,))
            
            if existing:
                logger.info(f"Protein {source_id} already exists with ID {existing[0][0]}")
                continue
            
            # Insert the protein
            logger.info(f"Inserting new protein {source_id}")
            protein_id = db.insert(
                "ecod_schema.protein",
                {
                    "pdb_id": args.pdb_id,
                    "chain_id": chain_id,
                    "source_id": source_id,
                    "length": len(sequence)
                },
                "id"
            )
            
            # Calculate MD5 hash of the sequence
            import hashlib
            md5_hash = hashlib.md5(sequence.encode()).hexdigest()
            
            # Insert the sequence
            db.insert(
                "ecod_schema.protein_sequence",
                {
                    "protein_id": protein_id,
                    "sequence": sequence,
                    "md5_hash": md5_hash
                }
            )
            
            logger.info(f"Inserted protein {source_id} with ID {protein_id} and its sequence")
            
        except Exception as e:
            logger.error(f"Error inserting sample protein: {str(e)}")
            success = False
    
    return 0 if success else 1

def _check_status(args: argparse.Namespace, db_config: Dict[str, Any]) -> int:
    """Check database status"""
    db = DBManager(db_config)
    
    try:
        # Check schemas
        schemas = []
        if args.schema:
            schemas = [args.schema]
        else:
            query = """
            SELECT schema_name 
            FROM information_schema.schemata
            WHERE schema_name NOT LIKE 'pg_%' 
              AND schema_name != 'information_schema'
            ORDER BY schema_name
            """
            rows = db.execute_query(query)
            schemas = [row[0] for row in rows if row and len(row) > 0]
        
        print(f"Found {len(schemas)} schemas: {', '.join(schemas)}")
        
        # For each schema, check table counts
        for schema in schemas:
            print(f"\nSchema: {schema}")
            
            query = f"""
            SELECT 
                tablename, 
                (SELECT COUNT(*) FROM {schema}.\"{{tablename}}\" AS count
            FROM 
                pg_tables 
            WHERE 
                schemaname = %s
            ORDER BY 
                tablename
            """
            rows = db.execute_dict_query(query, (schema,))
            
            for row in rows:
                print(f"  - {row['tablename']}: {row['count']} rows")
        
        return 0
        
    except Exception as e:
        logger.error(f"Error checking database status: {str(e)}")
        return 1