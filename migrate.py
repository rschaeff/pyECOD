# migrate.py
#!/usr/bin/env python3

import os
import sys
import argparse
import logging
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from ecod.config import ConfigManager
from ecod.db.migration_manager import MigrationManager

def main():
    parser = argparse.ArgumentParser(description='ECOD Database Migration Tool')
    parser.add_argument('--config', type=str, help='Path to configuration file')
    parser.add_argument('--migrations-dir', type=str, default='ecod/db/migrations',
                       help='Directory containing migration files')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    )
    
    # Load configuration
    config_manager = ConfigManager(args.config)
    db_config = config_manager.get_db_config()
    
    if not db_config:
        print("Error: Database configuration not found.")
        return 1
    
    # Run migrations
    try:
        migration_manager = MigrationManager(db_config, args.migrations_dir)
        migration_manager.apply_migrations()
        print("Migrations applied successfully!")
        return 0
    except Exception as e:
        print(f"Error applying migrations: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())