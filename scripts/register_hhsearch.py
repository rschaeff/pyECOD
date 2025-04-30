#!/usr/bin/env python3
"""
register_hhsearch.py - Register HHR files and convert to XML format

This script:
1. Finds HHR files that exist but aren't tracked in the database
2. Registers them in the database
3. Converts HHR files to XML format for domain analysis
4. Registers the XML files in the database
5. Updates process status to reflect completed HHSearch step

Usage:
    python register_hhsearch.py --batch-id <batch_id> [--force] [--chains <pdb_id>_<chain_id> ...]
"""

import os
import sys
import logging
import argparse
from typing import List, Dict, Any, Optional, Tuple

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext
from ecod.db import DBManager
from ecod.pipelines.hhsearch.processor import HHRToXMLConverter
from ecod.utils.hhsearch_utils import HHRParser
from ecod.utils.path_utils import (
    get_standardized_paths,
    find_files_with_legacy_paths,
    migrate_file_to_standard_path,
    get_file_db_path
)


def setup_logging(verbose=False, log_file=None):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO

    handlers = [logging.StreamHandler()]
    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )


class HHRegistrar:
    """Register HHSearch results and convert to XML for domain analysis"""

    def __init__(self, context=None):
        """Initialize with application context"""
        self.context = context or ApplicationContext()
        self.db = self.context.db
        self.config = self.context.config_manager.config
        self.logger = logging.getLogger("ecod.hh_registrar")

        # Initialize parser and converter
        self.parser = HHRParser(self.logger)
        self.converter = HHRToXMLConverter(self.logger)

    def register_batch_results(self, batch_id: int, force_regenerate: bool = False) -> int:
        """
        Find, convert, and register HHSearch results for a batch

        Args:
            batch_id: Batch ID to process
            force_regenerate: Force regeneration of XML files

        Returns:
            Number of registered files
        """
        # Get batch information
        batch_info = self._get_batch_info(batch_id)
        if not batch_info:
            self.logger.error(f"Batch {batch_id} not found")
            return 0

        # Get chains with HHR files
        chains = self._find_chains_with_hhr(batch_id, batch_info)
        if not chains:
            self.logger.warning(f"No chains with HHR files found in batch {batch_id}")
            return 0

        # Register files
        registered_count = 0
        for chain in chains:
            if self._register_chain_results(
                chain['process_id'],
                chain['pdb_id'],
                chain['chain_id'],
                batch_info['ref_version'],
                batch_info['base_path'],
                force_regenerate
            ):
                registered_count += 1

        self.logger.info(f"Registered HHR results for {registered_count} out of {len(chains)} chains")
        return registered_count

    def register_specific_chains(self, batch_id: int, chain_ids: List[str],
                               force_regenerate: bool = False) -> int:
        """
        Find, convert, and register HHSearch results for specific chains

        Args:
            batch_id: Batch ID
            chain_ids: List of chain IDs to process (format: "pdbid_chainid")
            force_regenerate: Force regeneration of XML files

        Returns:
            Number of registered files
        """
        # Get batch information
        batch_info = self._get_batch_info(batch_id)
        if not batch_info:
            self.logger.error(f"Batch {batch_id} not found")
            return 0

        # Get specific chains
        chains = []
        for chain_id in chain_ids:
            try:
                pdb_id, chain_letter = chain_id.split('_')
                result = self._find_specific_chain(batch_id, pdb_id, chain_letter)
                if result:
                    chains.append(result)
            except ValueError:
                self.logger.warning(f"Invalid chain ID format: {chain_id}, expected pdbid_chainid")

        if not chains:
            self.logger.warning(f"No specified chains found in batch {batch_id}")
            return 0

        # Register files
        registered_count = 0
        for chain in chains:
            if self._register_chain_results(
                chain['process_id'],
                chain['pdb_id'],
                chain['chain_id'],
                batch_info['ref_version'],
                batch_info['base_path'],
                force_regenerate
            ):
                registered_count += 1

        self.logger.info(f"Registered HHR results for {registered_count} out of {len(chains)} chains")
        return registered_count

    def _get_batch_info(self, batch_id: int) -> Optional[Dict[str, Any]]:
        """Get batch information"""
        query = """
        SELECT
            id,
            batch_name,
            base_path,
            ref_version
        FROM
            ecod_schema.batch
        WHERE
            id = %s
        """

        try:
            rows = self.db.execute_dict_query(query, (batch_id,))
            if rows:
                return rows[0]
        except Exception as e:
            self.logger.error(f"Error retrieving batch information: {str(e)}")

        return None

    def _find_chains_with_hhr(self, batch_id: int, batch_info: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Find chains with HHR files in a batch"""
        # First attempt to find chains with registered HHR files
        query = """
        SELECT
            p.id as protein_id,
            p.pdb_id,
            p.chain_id,
            ps.id as process_id,
            ps.relative_path,
            pf.file_path as hhr_path
        FROM
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        JOIN
            ecod_schema.process_file pf ON ps.id = pf.process_id
        WHERE
            ps.batch_id = %s
            AND pf.file_type = 'hhr'
            AND pf.file_exists = TRUE
        """

        chains = self.db.execute_dict_query(query, (batch_id,))
        if chains:
            self.logger.info(f"Found {len(chains)} chains with registered HHR files")
            return chains

        # If no registered HHR files, search on filesystem
        self.logger.info("No registered HHR files found, searching filesystem...")

        query = """
        SELECT
            p.id as protein_id,
            p.pdb_id,
            p.chain_id,
            ps.id as process_id,
            ps.relative_path,
            b.base_path
        FROM
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        JOIN
            ecod_schema.batch b ON ps.batch_id = b.id
        WHERE
            ps.batch_id = %s
            AND ps.status = 'success'
        """

        result = []
        chains = self.db.execute_dict_query(query, (batch_id,))

        if not chains:
            return []

        # Check if HHR files exist in filesystem
        ref_version = batch_info['ref_version']
        base_path = batch_info['base_path']

        for chain in chains:
            pdb_id = chain['pdb_id']
            chain_id = chain['chain_id']

            # Find files using both standard and legacy paths
            file_info = find_files_with_legacy_paths(base_path, pdb_id, chain_id, ref_version)

            # If HHR file exists at any location, add to result
            if file_info['hhr']['exists_at']:
                chain['hhr_path'] = file_info['hhr']['exists_at']
                result.append(chain)

        self.logger.info(f"Found {len(result)} chains with unregistered HHR files")
        return result

    def _find_specific_chain(self, batch_id: int, pdb_id: str, chain_id: str) -> Optional[Dict[str, Any]]:
        """Find a specific chain in a batch and its associated HHsearch file"""

        query = """
        SELECT
            p.id as protein_id,
            p.pdb_id,
            p.chain_id,
            ps.id as process_id,
            ps.relative_path,
            b.base_path
        FROM
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        JOIN
            ecod_schema.batch b ON ps.batch_id = b.id
        WHERE
            ps.batch_id = %s
            AND p.pdb_id = %s
            AND p.chain_id = %s
        LIMIT 1
        """

        try:
            rows = self.db.execute_dict_query(query, (batch_id, pdb_id, chain_id))
            if not rows:
                self.logger.warning(f"Chain {pdb_id}_{chain_id} not found in batch {batch_id}")
                return None

            # Create a regular dictionary from the database result
            chain_data = dict(rows[0])

            # Get reference version
            ref_version = self._get_batch_ref_version(batch_id)
            base_path = chain_data['base_path']

            # Use path_utils to find standard and legacy paths
            file_info = find_files_with_legacy_paths(base_path, pdb_id, chain_id, ref_version)

            # Check if HHR file exists at any location
            if file_info['hhr']['exists_at']:
                hhr_path = file_info['hhr']['exists_at']
                self.logger.info(f"Found HHR file at: {hhr_path}")

                # Create a new dictionary with all chain data plus the hhr_path
                result = dict(chain_data)
                result['hhr_path'] = hhr_path
                return result

            # If no HHR file found with the utility, log the paths that were checked
            self.logger.warning(f"No HHR file found for {pdb_id}_{chain_id}")
            self.logger.debug(f"Standard path checked: {file_info['hhr']['standard']}")
            if file_info['hhr']['legacy']:
                self.logger.debug(f"Legacy path checked: {file_info['hhr']['legacy']}")

            return None

        except Exception as e:
            self.logger.error(f"Error finding chain {pdb_id}_{chain_id}: {str(e)}")
            self.logger.exception("Stack trace")
            return None

        def _get_batch_ref_version(self, batch_id: int) -> str:
            """Get reference version for a batch"""
            query = "SELECT ref_version FROM ecod_schema.batch WHERE id = %s"

            try:
                rows = self.db.execute_query(query, (batch_id,))
                if rows:
                    return rows[0][0]
            except Exception as e:
                self.logger.error(f"Error getting reference version: {str(e)}")

            # Default to configured reference version
            return self.config.get('reference', {}).get('current_version', 'develop291')

    def _register_chain_results(self, process_id: int, pdb_id: str, chain_id: str,
                              ref_version: str, base_path: str,
                              force_regenerate: bool = False) -> bool:
        """Register HHSearch results for a chain"""
        try:
            # Get standardized paths
            paths = get_standardized_paths(base_path, pdb_id, chain_id, ref_version)
            hhr_file = paths['hhr']
            xml_file = paths['hh_xml']

            # Check if XML file already exists (unless force regenerate)
            if os.path.exists(xml_file) and os.path.getsize(xml_file) > 0 and not force_regenerate:
                self.logger.info(f"XML file already exists for {pdb_id}_{chain_id}, using existing file")

                # Register XML file if needed
                self._check_and_register_file(process_id, 'hh_xml', xml_file, base_path)

                # Update process status
                self._update_process_status(process_id, "hhsearch_complete")

                return True

            # Check if HHR file exists at standard location
            if not os.path.exists(hhr_file) or os.path.getsize(hhr_file) == 0:
                self.logger.info(f"HHR file not found at standard location: {hhr_file}")

                # Find files using both standard and legacy paths
                file_info = find_files_with_legacy_paths(base_path, pdb_id, chain_id, ref_version)

                # If HHR file exists at legacy location, migrate to standard
                if file_info['hhr']['legacy'] and file_info['hhr']['legacy'] != hhr_file:
                    source_path = file_info['hhr']['legacy']
                    self.logger.info(f"Found HHR file at legacy location: {source_path}")

                    # Migrate to standard location
                    success = migrate_file_to_standard_path(source_path, hhr_file)
                    if not success:
                        self.logger.error(f"Failed to migrate HHR file to standard location")
                        return False
                else:
                    self.logger.error(f"No HHR file found for {pdb_id}_{chain_id}")
                    return False

            # Register HHR file
            self._check_and_register_file(process_id, 'hhr', hhr_file, base_path)

            # Convert HHR to XML
            self.logger.info(f"Converting {hhr_file} to XML...")
            hhr_data = self.parser.parse(hhr_file)
            if not hhr_data:
                self.logger.error(f"Failed to parse HHR file: {hhr_file}")
                return False

            # Convert to XML
            xml_string = self.converter.convert(hhr_data, pdb_id, chain_id, ref_version)
            if not xml_string:
                self.logger.error(f"Failed to convert HHR data to XML for {pdb_id}_{chain_id}")
                return False

            # Save XML file
            if not self.converter.save(xml_string, xml_file):
                self.logger.error(f"Failed to save XML file: {xml_file}")
                return False

            # Register XML file
            self._check_and_register_file(process_id, 'hh_xml', xml_file, base_path)

            # Update process status
            self._update_process_status(process_id, "hhsearch_complete")

            self.logger.info(f"Successfully registered HHSearch results for {pdb_id}_{chain_id}")
            return True

        except Exception as e:
            self.logger.error(f"Error registering results for {pdb_id}_{chain_id}: {str(e)}")
            return False

    def _check_and_register_file(self, process_id: int, file_type: str,
                                file_path: str, base_path: str) -> bool:
        """Check if file is registered in database and register if not"""
        try:
            # Check if file already registered
            query = """
            SELECT id FROM ecod_schema.process_file
            WHERE process_id = %s AND file_type = %s
            """

            existing = self.db.execute_query(query, (process_id, file_type))

            # Get relative path for database storage
            rel_path = get_file_db_path(base_path, file_path)

            if existing:
                # Update existing record
                self.db.update(
                    "ecod_schema.process_file",
                    {
                        "file_path": rel_path,
                        "file_exists": True,
                        "file_size": os.path.getsize(file_path),
                        "last_checked": "CURRENT_TIMESTAMP"
                    },
                    "id = %s",
                    (existing[0][0],)
                )
                self.logger.info(f"Updated {file_type} record for process {process_id}")
            else:
                # Create new record
                self.db.insert(
                    "ecod_schema.process_file",
                    {
                        "process_id": process_id,
                        "file_type": file_type,
                        "file_path": rel_path,
                        "file_exists": True,
                        "file_size": os.path.getsize(file_path),
                        "last_checked": "CURRENT_TIMESTAMP"
                    }
                )
                self.logger.info(f"Created new {file_type} record for process {process_id}")

            return True
        except Exception as e:
            self.logger.error(f"Error registering file {file_path}: {str(e)}")
            return False

    def _update_process_status(self, process_id: int, stage: str, error_message: str = None) -> bool:
        """Update process status in database"""
        try:
            status = "error" if error_message else "success"

            update_data = {
                "current_stage": stage,
                "status": status,
                "updated_at": "CURRENT_TIMESTAMP"
            }

            if error_message:
                update_data["error_message"] = error_message

            self.db.update(
                "ecod_schema.process_status",
                update_data,
                "id = %s",
                (process_id,)
            )

            return True
        except Exception as e:
            self.logger.error(f"Error updating process status for {process_id}: {str(e)}")
            return False


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Register and convert HHSearch results')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to process')
    parser.add_argument('--chains', nargs='+', default=None,
                      help='Specific chains to process (format: pdbid_chainid)')
    parser.add_argument('--force', action='store_true',
                      help='Force reprocessing of already processed results')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')

    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)

    logger = logging.getLogger("ecod.register_hhsearch")
    logger.info(f"Starting HHSearch registration and conversion for batch {args.batch_id}")

    # Use absolute path for config if provided
    if args.config and not os.path.isabs(args.config):
        args.config = os.path.abspath(args.config)

    # Initialize application context
    context = ApplicationContext(args.config)

    # Create registrar
    registrar = HHRegistrar(context)

    try:
        # Process batch
        if args.chains:
            result = registrar.register_specific_chains(
                args.batch_id,
                args.chains,
                args.force
            )
        else:
            result = registrar.register_batch_results(
                args.batch_id,
                args.force
            )

        logger.info(f"Successfully registered {result} HHR files and converted them to XML")

        if result == 0:
            return 1

        return 0
    except Exception as e:
        logger.error(f"Error processing batch: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
