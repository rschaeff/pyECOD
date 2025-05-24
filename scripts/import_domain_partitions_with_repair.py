#!/usr/bin/env python3
"""
Enhanced Import Domain Partitions Script with Complete Status Synchronization

This script reads domain partition XML files and imports them into the pdb_analysis partition tables
with comprehensive status tracking synchronization between ecod_schema and pdb_analysis.

Key enhancements:
- Complete status synchronization with ecod_schema.process_status
- Conflict resolution strategies (skip, update, replace, version)
- Repair mode for previously imported batches
- Process ID tracking and status updates
- Validation and verification modes
- Progress reporting with status verification

Usage:
    python import_domain_partitions_with_repair.py --config config.yml [options]

Status Synchronization:
    --update-process-status  Update ecod_schema.process_status after import (default: True)
    --verify-sync           Verify synchronization after import (default: True)
"""

import os
import sys
import logging
import argparse
import datetime
import xml.etree.ElementTree as ET
import psycopg2
from psycopg2.extras import execute_values, RealDictCursor
import yaml
from datetime import datetime
from enum import Enum
from typing import List, Dict, Any, Optional, Tuple


class ConflictStrategy(Enum):
    """Strategies for handling conflicts during import."""
    SKIP = "skip"           # INSERT ... ON CONFLICT DO NOTHING
    UPDATE = "update"       # INSERT ... ON CONFLICT DO UPDATE
    REPLACE = "replace"     # DELETE existing + INSERT new
    VERSION = "version"     # Keep existing, create new version


def setup_logging(log_file=None, verbose=False):
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    format_str = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

    handlers = [logging.StreamHandler()]
    if log_file:
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(level=level, format=format_str, handlers=handlers)
    return logging.getLogger(__name__)


def parse_config(config_path):
    """Parse configuration file with optional local overrides."""
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)

        # Check for local config in the same directory
        local_config_path = os.path.join(
            os.path.dirname(config_path),
            'config.local.yml'
        )

        # If local config exists, merge it with main config
        if os.path.exists(local_config_path):
            with open(local_config_path, 'r') as f:
                local_config = yaml.safe_load(f)

            if local_config:
                config = deep_merge(config, local_config)
                logging.debug(f"Merged local configuration from {local_config_path}")

        return config
    except Exception as e:
        logging.error(f"Error parsing config file: {str(e)}")
        raise


def deep_merge(dict1, dict2):
    """Deep merge two dictionaries."""
    result = dict1.copy()

    for key, value in dict2.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = value

    return result


def get_db_connection(config):
    """Create database connection from config with proper error handling."""
    db_config = config.get('database', {})

    try:
        conn_info = {
            'host': db_config.get('host', 'dione'),
            'port': db_config.get('port', 45000),
            'dbname': db_config.get('name', 'ecod_protein'),
            'user': db_config.get('user', 'ecod')
        }
        logging.debug(f"Connecting to database: {conn_info}")

        conn = psycopg2.connect(
            host=db_config.get('host', 'dione'),
            port=db_config.get('port', 45000),
            dbname=db_config.get('name', 'ecod_protein'),
            user=db_config.get('user', 'ecod'),
            password=db_config.get('password', '')
        )

        return conn
    except psycopg2.Error as e:
        logging.error(f"Database connection error: {e}")
        raise


def resolve_file_path(base_path, file_path):
    """Resolve a relative file path to its absolute location."""
    if os.path.isabs(file_path):
        logging.warning(f"Found unexpected absolute path: {file_path}")
        return file_path

    if file_path.startswith('domains/'):
        filename_part = file_path[len('domains/'):]
        domains_dir = os.path.join(base_path, 'domains')
        return os.path.join(domains_dir, filename_part)

    return os.path.join(base_path, file_path)


def get_partition_files_from_db(conn, batch_id=None, limit=None, rep_only=True, repair_mode=False):
    """Retrieve domain partition file paths from ecod_schema.process_file table with process ID tracking."""
    with conn.cursor() as cur:
        if repair_mode:
            # In repair mode, get all files regardless of previous import status
            query = """
                SELECT
                    pf.id, ps.batch_id, pf.file_type, pf.file_path,
                    p.pdb_id, p.chain_id, pf.last_checked as timestamp,
                    CASE WHEN pp.id IS NOT NULL THEN true ELSE false END as previously_imported,
                    ps.id as process_id
                FROM
                    ecod_schema.process_file pf
                JOIN
                    ecod_schema.process_status ps ON pf.process_id = ps.id
                JOIN
                    ecod_schema.protein p ON ps.protein_id = p.id
                LEFT JOIN
                    pdb_analysis.partition_proteins pp ON (
                        pp.pdb_id = p.pdb_id
                        AND pp.chain_id = p.chain_id
                        AND pp.batch_id = ps.batch_id
                    )
                WHERE
                    pf.file_type = 'domain_partition'
                    AND pf.file_exists = true
                """
        else:
            # Normal mode: only get files not previously imported
            query = """
                SELECT
                    pf.id, ps.batch_id, pf.file_type, pf.file_path,
                    p.pdb_id, p.chain_id, pf.last_checked as timestamp,
                    false as previously_imported,
                    ps.id as process_id
                FROM
                    ecod_schema.process_file pf
                JOIN
                    ecod_schema.process_status ps ON pf.process_id = ps.id
                JOIN
                    ecod_schema.protein p ON ps.protein_id = p.id
                LEFT JOIN
                    pdb_analysis.partition_proteins pp ON (
                        pp.pdb_id = p.pdb_id
                        AND pp.chain_id = p.chain_id
                        AND pp.batch_id = ps.batch_id
                    )
                WHERE
                    pf.file_type = 'domain_partition'
                    AND pf.file_exists = true
                    AND pp.id IS NULL  -- Only unimported files
                """

        params = []

        if rep_only:
            query += " AND ps.is_representative = true"

        if batch_id:
            query += " AND ps.batch_id = %s"
            params.append(batch_id)

        query += " ORDER BY pf.id"

        if limit:
            query += " LIMIT %s"
            params.append(limit)

        cur.execute(query, params)
        return cur.fetchall()


def update_import_status(conn, pdb_id, chain_id, batch_id, process_id=None, success=True, error_message=None):
    """
    Update process status after successful import to pdb_analysis.

    This is the critical synchronization step that was missing.
    """
    try:
        if process_id:
            # Direct update using process_id (faster and more reliable)
            if success:
                query = """
                UPDATE ecod_schema.process_status
                SET current_stage = 'pdb_analysis_imported',
                    status = 'success',
                    error_message = NULL,
                    updated_at = CURRENT_TIMESTAMP
                WHERE id = %s
                """
                params = (process_id,)
            else:
                query = """
                UPDATE ecod_schema.process_status
                SET current_stage = 'pdb_analysis_import_failed',
                    status = 'error',
                    error_message = %s,
                    updated_at = CURRENT_TIMESTAMP
                WHERE id = %s
                """
                params = (error_message[:500] if error_message else "Import failed", process_id)
        else:
            # Fallback: update by protein and batch (slower but works without process_id)
            if success:
                query = """
                UPDATE ecod_schema.process_status ps
                SET current_stage = 'pdb_analysis_imported',
                    status = 'success',
                    error_message = NULL,
                    updated_at = CURRENT_TIMESTAMP
                FROM ecod_schema.protein p
                WHERE ps.protein_id = p.id
                  AND p.pdb_id = %s
                  AND p.chain_id = %s
                  AND ps.batch_id = %s
                """
                params = (pdb_id, chain_id, batch_id)
            else:
                query = """
                UPDATE ecod_schema.process_status ps
                SET current_stage = 'pdb_analysis_import_failed',
                    status = 'error',
                    error_message = %s,
                    updated_at = CURRENT_TIMESTAMP
                FROM ecod_schema.protein p
                WHERE ps.protein_id = p.id
                  AND p.pdb_id = %s
                  AND p.chain_id = %s
                  AND ps.batch_id = %s
                """
                params = (error_message[:500] if error_message else "Import failed", pdb_id, chain_id, batch_id)

        with conn.cursor() as cur:
            cur.execute(query, params)
            updated_count = cur.rowcount

            if updated_count > 0:
                logging.debug(f"Updated import status for {pdb_id}_{chain_id} (process_id: {process_id})")
            else:
                logging.warning(f"No process status updated for {pdb_id}_{chain_id} (process_id: {process_id})")

            return updated_count > 0

    except Exception as e:
        logging.error(f"Error updating import status for {pdb_id}_{chain_id}: {str(e)}")
        return False


def verify_import_synchronization(conn, batch_id, logger):
    """Verify that imports are properly synchronized with process status."""

    try:
        # Query to check synchronization for this batch
        query = """
        SELECT
            COUNT(CASE WHEN ps.is_representative = true THEN 1 END) as rep_proteins,
            COUNT(pp.id) as imported_partitions,
            COUNT(CASE WHEN ps.is_representative = true
                       AND ps.current_stage = 'pdb_analysis_imported'
                       AND ps.status = 'success' THEN 1 END) as marked_imported,
            COUNT(CASE WHEN ps.is_representative = true
                       AND ps.current_stage = 'domain_partition_complete'
                       AND ps.status = 'success' THEN 1 END) as marked_complete
        FROM ecod_schema.batch b
        JOIN ecod_schema.process_status ps ON b.id = ps.batch_id
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        LEFT JOIN pdb_analysis.partition_proteins pp ON (
            pp.pdb_id = p.pdb_id AND pp.chain_id = p.chain_id AND pp.batch_id = ps.batch_id
        )
        WHERE b.id = %s
        """

        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(query, (batch_id,))
            result = cur.fetchone()

            logger.info(f"Import synchronization verification for batch {batch_id}:")
            logger.info(f"  Representative proteins: {result['rep_proteins']}")
            logger.info(f"  Imported partitions: {result['imported_partitions']}")
            logger.info(f"  Marked as imported: {result['marked_imported']}")
            logger.info(f"  Marked as complete: {result['marked_complete']}")

            # Check for synchronization
            if result['imported_partitions'] == result['marked_imported']:
                logger.info(f"✅ Batch {batch_id}: Import status fully synchronized")
                return True
            else:
                gap = result['imported_partitions'] - result['marked_imported']
                logger.warning(f"⚠️ Batch {batch_id}: Import gap of {gap} proteins")
                return False

    except Exception as e:
        logger.error(f"Error verifying import synchronization: {e}")
        return False


def update_batch_import_completion(conn, batch_id, logger):
    """Update batch completion status after import operations."""

    try:
        # Calculate import completion status
        query = """
        SELECT
            COUNT(CASE WHEN ps.is_representative = true THEN 1 END) as total_rep,
            COUNT(CASE WHEN ps.is_representative = true
                       AND ps.current_stage IN ('pdb_analysis_imported', 'domain_partition_complete')
                       AND ps.status = 'success' THEN 1 END) as completed_rep,
            COUNT(CASE WHEN ps.is_representative = true
                       AND ps.status = 'error' THEN 1 END) as error_rep
        FROM ecod_schema.process_status ps
        WHERE ps.batch_id = %s
        """

        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(query, (batch_id,))
            result = cur.fetchone()

            total = result['total_rep']
            completed = result['completed_rep']
            errors = result['error_rep']

            if total == 0:
                return

            # Determine batch status
            completion_rate = completed / total

            if completion_rate >= 0.95:  # 95% or more complete
                new_status = 'import_complete'
            elif completion_rate >= 0.8:  # 80% or more complete
                new_status = 'import_mostly_complete'
            elif completion_rate >= 0.5:  # 50% or more complete
                new_status = 'import_partial'
            else:
                new_status = 'import_incomplete'

            # Update batch status
            cur.execute(
                """
                UPDATE ecod_schema.batch
                SET status = %s,
                    completed_items = %s,
                    updated_at = CURRENT_TIMESTAMP
                WHERE id = %s
                """,
                (new_status, completed, batch_id)
            )

            logger.info(f"Updated batch {batch_id} import status to '{new_status}' "
                       f"({completed}/{total} completed, {errors} errors)")

    except Exception as e:
        logger.error(f"Error updating batch import completion: {e}")


def check_existing_partition(conn, pdb_id, chain_id, batch_id):
    """Check if partition already exists and return details."""
    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(
            """
            SELECT id, process_version, timestamp, is_classified,
                   sequence_length, coverage, residues_assigned
            FROM pdb_analysis.partition_proteins
            WHERE pdb_id = %s AND chain_id = %s AND batch_id = %s
            ORDER BY timestamp DESC
            LIMIT 1
            """,
            (pdb_id, chain_id, batch_id)
        )
        return cur.fetchone()


def delete_existing_partition(conn, pdb_id, chain_id, batch_id):
    """Delete existing partition and all associated data."""
    with conn.cursor() as cur:
        # Get partition protein ID
        cur.execute(
            """
            SELECT id FROM pdb_analysis.partition_proteins
            WHERE pdb_id = %s AND chain_id = %s AND batch_id = %s
            """,
            (pdb_id, chain_id, batch_id)
        )
        result = cur.fetchone()
        if not result:
            return 0

        protein_id = result[0]

        # Delete in reverse dependency order
        cur.execute(
            "DELETE FROM pdb_analysis.domain_evidence WHERE domain_id IN "
            "(SELECT id FROM pdb_analysis.partition_domains WHERE protein_id = %s)",
            (protein_id,)
        )

        cur.execute(
            "DELETE FROM pdb_analysis.partition_domains WHERE protein_id = %s",
            (protein_id,)
        )

        cur.execute(
            "DELETE FROM pdb_analysis.partition_proteins WHERE id = %s",
            (protein_id,)
        )

        return cur.rowcount


def parse_partition_xml(content, file_id=None, pdb_id=None, chain_id=None, batch_id=None, timestamp=None):
    """Parse domain partition XML content."""
    try:
        # Parse the XML
        root = ET.fromstring(content)

        # Extract basic metadata
        parsed_pdb_id = root.get('pdb_id', '')
        parsed_chain_id = root.get('chain_id', '')
        reference = root.get('reference', '')
        is_classified = root.get('is_classified', 'false').lower() == 'true'

        # Use provided values if XML doesn't have them
        if not parsed_pdb_id and pdb_id:
            parsed_pdb_id = pdb_id
        if not parsed_chain_id and chain_id:
            parsed_chain_id = chain_id

        # Get metadata
        metadata = {}
        metadata_elem = root.find('./metadata')
        if metadata_elem is not None:
            for child in metadata_elem:
                tag = child.tag
                value = child.text
                if value and value.strip():
                    try:
                        if tag in ('sequence_length', 'domain_count', 'residues_assigned',
                                  'domains_with_evidence', 'fully_classified_domains'):
                            metadata[tag] = int(value)
                        elif tag in ('coverage',):
                            metadata[tag] = float(value)
                        elif tag in ('timestamp',):
                            metadata[tag] = datetime.fromisoformat(value)
                        else:
                            metadata[tag] = value
                    except (ValueError, TypeError):
                        metadata[tag] = value

        # Extract each domain
        domains = []
        for idx, domain_elem in enumerate(root.findall('./domains/domain'), 1):
            domain = {
                'domain_number': idx,
                'domain_id': f"{parsed_pdb_id}_{parsed_chain_id}_d{idx}",
                'start_pos': int(domain_elem.get('start', 0)),
                'end_pos': int(domain_elem.get('end', 0)),
                'range': domain_elem.get('range', ''),
                'source': domain_elem.get('source', ''),
                'source_id': domain_elem.get('source_id', ''),
                'confidence': float(domain_elem.get('confidence', 0)),
                't_group': domain_elem.get('t_group'),
                'h_group': domain_elem.get('h_group'),
                'x_group': domain_elem.get('x_group'),
                'a_group': domain_elem.get('a_group'),
                'is_manual_rep': domain_elem.get('is_manual_rep', 'False').lower() == 'true',
                'is_f70': domain_elem.get('is_f70', 'False').lower() == 'true',
                'is_f40': domain_elem.get('is_f40', 'False').lower() == 'true',
                'is_f99': domain_elem.get('is_f99', 'False').lower() == 'true',
            }

            # Extract evidence items
            evidence_list = []
            evidence_elem = domain_elem.find('./evidence_list')
            if evidence_elem is not None:
                for ev in evidence_elem.findall('./evidence'):
                    evidence = {
                        'evidence_type': ev.get('type', ''),
                        'source_id': ev.get('source_id', ''),
                        'domain_ref_id': ev.get('domain_id', ''),
                        'hit_id': ev.get('hit_id', ''),
                        'pdb_id': ev.get('pdb_id', ''),
                        'chain_id': ev.get('chain_id', ''),
                        'confidence': float(ev.get('confidence', 0)) if ev.get('confidence') else None,
                        'probability': float(ev.get('probability', 0)) if ev.get('probability') else None,
                        'evalue': float(ev.get('evalue', 0)) if ev.get('evalue') else None,
                        'score': float(ev.get('score', 0)) if ev.get('score') else None,
                        'hsp_count': int(ev.get('hsp_count', 0)) if ev.get('hsp_count') else None,
                        'is_discontinuous': ev.get('discontinuous', 'False').lower() == 'true',
                        't_group': ev.get('t_group'),
                        'h_group': ev.get('h_group'),
                        'x_group': ev.get('x_group'),
                        'a_group': ev.get('a_group'),
                    }

                    # Get query and hit ranges
                    query_range = ev.find('./query_range')
                    hit_range = ev.find('./hit_range')

                    if query_range is not None and query_range.text:
                        evidence['query_range'] = query_range.text.strip()
                    if hit_range is not None and hit_range.text:
                        evidence['hit_range'] = hit_range.text.strip()

                    evidence_list.append(evidence)

            domain['evidence'] = evidence_list
            domains.append(domain)

        # Parse timestamp if available
        if 'timestamp' in metadata:
            partition_timestamp = metadata['timestamp']
        elif timestamp:
            partition_timestamp = timestamp
        else:
            partition_timestamp = datetime.now()

        return {
            'file_id': file_id,
            'pdb_id': parsed_pdb_id,
            'chain_id': parsed_chain_id,
            'batch_id': batch_id,
            'reference_version': reference,
            'is_classified': is_classified,
            'sequence_length': metadata.get('sequence_length', 0),
            'coverage': metadata.get('coverage', 0.0),
            'residues_assigned': metadata.get('residues_assigned', 0),
            'domains_with_evidence': metadata.get('domains_with_evidence', 0),
            'fully_classified_domains': metadata.get('fully_classified_domains', 0),
            'timestamp': partition_timestamp,
            'domains': domains
        }

    except Exception as e:
        logging.error(f"Error parsing XML content for file_id {file_id}: {str(e)}")
        return None


def insert_domain_partition_with_strategy(conn, partition_data, process_version, conflict_strategy,
                                         process_id=None, update_status=True, dry_run=False):
    """Insert domain partition using specified conflict resolution strategy with status tracking."""

    if dry_run:
        logging.info(f"DRY RUN - Would insert partition for {partition_data['pdb_id']}_{partition_data['chain_id']} "
                    f"with {len(partition_data['domains'])} domains using strategy {conflict_strategy.value}")
        return None, []

    pdb_id = partition_data['pdb_id']
    chain_id = partition_data['chain_id']
    batch_id = partition_data['batch_id']

    try:
        # Check if already exists
        existing = check_existing_partition(conn, pdb_id, chain_id, batch_id)

        if existing:
            logging.info(f"Found existing partition for {pdb_id}_{chain_id} "
                        f"(version {existing['process_version']}, timestamp {existing['timestamp']})")

            if conflict_strategy == ConflictStrategy.SKIP:
                logging.info(f"Skipping due to SKIP strategy")
                # Still update status to show it was processed
                if update_status:
                    update_import_status(conn, pdb_id, chain_id, batch_id, process_id, True)
                return existing['id'], []

            elif conflict_strategy == ConflictStrategy.REPLACE:
                logging.info(f"Deleting existing partition due to REPLACE strategy")
                delete_existing_partition(conn, pdb_id, chain_id, batch_id)

            elif conflict_strategy == ConflictStrategy.VERSION:
                # Increment version for new record
                old_version = existing.get('process_version', '1.0')
                try:
                    version_parts = old_version.split('.')
                    major, minor = int(version_parts[0]), int(version_parts[1]) if len(version_parts) > 1 else 0
                    process_version = f"{major}.{minor + 1}"
                    logging.info(f"Creating new version {process_version}")
                except:
                    process_version = f"{process_version}.1"

        protein_id = None
        domain_ids = []

        with conn.cursor() as cur:
            if conflict_strategy == ConflictStrategy.UPDATE and existing:
                # Update existing partition protein record
                cur.execute(
                    """
                    UPDATE pdb_analysis.partition_proteins SET
                        timestamp = %s, reference_version = %s, is_classified = %s,
                        sequence_length = %s, coverage = %s, residues_assigned = %s,
                        domains_with_evidence = %s, fully_classified_domains = %s,
                        source_file_id = %s, process_version = %s
                    WHERE pdb_id = %s AND chain_id = %s AND batch_id = %s
                    RETURNING id
                    """,
                    (
                        partition_data['timestamp'], partition_data['reference_version'],
                        partition_data['is_classified'], partition_data['sequence_length'],
                        partition_data['coverage'], partition_data['residues_assigned'],
                        partition_data['domains_with_evidence'], partition_data['fully_classified_domains'],
                        partition_data['file_id'], process_version,
                        pdb_id, chain_id, batch_id
                    )
                )
                protein_id = cur.fetchone()[0]

                # Delete existing domains and evidence for update
                cur.execute(
                    "DELETE FROM pdb_analysis.domain_evidence WHERE domain_id IN "
                    "(SELECT id FROM pdb_analysis.partition_domains WHERE protein_id = %s)",
                    (protein_id,)
                )
                cur.execute(
                    "DELETE FROM pdb_analysis.partition_domains WHERE protein_id = %s",
                    (protein_id,)
                )

            else:
                # Insert new partition protein record
                if conflict_strategy == ConflictStrategy.SKIP:
                    # Use ON CONFLICT DO NOTHING for protein level
                    cur.execute(
                        """
                        INSERT INTO pdb_analysis.partition_proteins (
                            pdb_id, chain_id, batch_id, timestamp, reference_version, is_classified,
                            sequence_length, coverage, residues_assigned, domains_with_evidence,
                            fully_classified_domains, source_file_id, process_version
                        ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                        ON CONFLICT (pdb_id, chain_id, batch_id) DO NOTHING
                        RETURNING id
                        """,
                        (
                            pdb_id, chain_id, batch_id,
                            partition_data['timestamp'], partition_data['reference_version'],
                            partition_data['is_classified'], partition_data['sequence_length'],
                            partition_data['coverage'], partition_data['residues_assigned'],
                            partition_data['domains_with_evidence'], partition_data['fully_classified_domains'],
                            partition_data['file_id'], process_version
                        )
                    )
                    result = cur.fetchone()
                    if result:
                        protein_id = result[0]
                    else:
                        # Record was skipped due to conflict
                        if update_status:
                            update_import_status(conn, pdb_id, chain_id, batch_id, process_id, True)
                        return existing['id'] if existing else None, []

                else:
                    # Regular insert for REPLACE and VERSION strategies
                    cur.execute(
                        """
                        INSERT INTO pdb_analysis.partition_proteins (
                            pdb_id, chain_id, batch_id, timestamp, reference_version, is_classified,
                            sequence_length, coverage, residues_assigned, domains_with_evidence,
                            fully_classified_domains, source_file_id, process_version
                        ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                        RETURNING id
                        """,
                        (
                            pdb_id, chain_id, batch_id,
                            partition_data['timestamp'], partition_data['reference_version'],
                            partition_data['is_classified'], partition_data['sequence_length'],
                            partition_data['coverage'], partition_data['residues_assigned'],
                            partition_data['domains_with_evidence'], partition_data['fully_classified_domains'],
                            partition_data['file_id'], process_version
                        )
                    )
                    protein_id = cur.fetchone()[0]

            # Insert domains (if we have a protein_id)
            if protein_id and partition_data['domains']:
                domain_values = []
                for domain in partition_data['domains']:
                    domain_values.append((
                        protein_id, domain['domain_number'], domain['domain_id'],
                        domain['start_pos'], domain['end_pos'], domain['range'],
                        domain['source'], domain['source_id'], domain['confidence'],
                        domain['t_group'], domain['h_group'], domain['x_group'], domain['a_group'],
                        domain['is_manual_rep'], domain['is_f70'], domain['is_f40'], domain['is_f99']
                    ))

                if domain_values:
                    domain_results = execute_values(
                        cur,
                        """
                        INSERT INTO pdb_analysis.partition_domains (
                            protein_id, domain_number, domain_id, start_pos, end_pos, range,
                            source, source_id, confidence, t_group, h_group, x_group, a_group,
                            is_manual_rep, is_f70, is_f40, is_f99
                        ) VALUES %s
                        RETURNING id
                        """,
                        domain_values,
                        fetch=True
                    )

                    # Insert evidence for each domain
                    evidence_values = []
                    for i, domain_id in enumerate(domain_results):
                        domain_ids.append(domain_id)
                        for evidence in partition_data['domains'][i].get('evidence', []):
                            evidence_values.append((
                                domain_id, evidence['evidence_type'], evidence.get('source_id'),
                                evidence.get('domain_ref_id'), evidence.get('hit_id'),
                                evidence.get('pdb_id'), evidence.get('chain_id'),
                                evidence.get('confidence'), evidence.get('probability'),
                                evidence.get('evalue'), evidence.get('score'),
                                evidence.get('hsp_count'), evidence.get('is_discontinuous'),
                                evidence.get('t_group'), evidence.get('h_group'),
                                evidence.get('x_group'), evidence.get('a_group'),
                                evidence.get('query_range'), evidence.get('hit_range')
                            ))

                    if evidence_values:
                        execute_values(
                            cur,
                            """
                            INSERT INTO pdb_analysis.domain_evidence (
                                domain_id, evidence_type, source_id, domain_ref_id, hit_id, pdb_id, chain_id,
                                confidence, probability, evalue, score, hsp_count, is_discontinuous,
                                t_group, h_group, x_group, a_group, query_range, hit_range
                            ) VALUES %s
                            """,
                            evidence_values
                        )

        # Update process status to show successful import
        if update_status and protein_id:
            success = update_import_status(conn, pdb_id, chain_id, batch_id, process_id, True)
            if success:
                logging.debug(f"Updated process status for successful import: {pdb_id}_{chain_id}")

        return protein_id, domain_ids

    except Exception as e:
        logging.error(f"Error inserting partition for {pdb_id}_{chain_id}: {str(e)}")

        # Update process status to show failed import
        if update_status:
            update_import_status(conn, pdb_id, chain_id, batch_id, process_id, False, str(e))

        raise


def validate_partition_data(partition_data):
    """Validate partition data before import."""
    issues = []

    if not partition_data['pdb_id'] or not partition_data['chain_id']:
        issues.append("Missing PDB ID or chain ID")

    if not partition_data['domains']:
        issues.append("No domains found")

    for i, domain in enumerate(partition_data['domains']):
        if domain['start_pos'] >= domain['end_pos']:
            issues.append(f"Domain {i+1}: Invalid range {domain['start_pos']}-{domain['end_pos']}")

        if not domain['evidence']:
            issues.append(f"Domain {i+1}: No evidence found")

    return issues


def read_file_content(file_path):
    """Read content from a file path."""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            return f.read()
    except Exception as e:
        logging.error(f"Error reading file {file_path}: {str(e)}")
        return None


def main():
    """Main function with comprehensive status synchronization."""
    parser = argparse.ArgumentParser(description='Import domain partitions with complete status synchronization')
    parser.add_argument('--config', required=True, help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, help='Process only a specific batch')
    parser.add_argument('--limit', type=int, help='Process only N files')
    parser.add_argument('--dry-run', action='store_true', help='Parse but don\'t write to database')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
    parser.add_argument('--log-file', help='Path to log file')
    parser.add_argument('--process-version', default='1.0', help='Version identifier for this processing run')
    parser.add_argument('--all-proteins', action='store_true', help='Process all proteins, not just representatives')

    # Repair mode options
    parser.add_argument('--repair-mode', action='store_true',
                       help='Enable repair mode (process previously imported files)')
    parser.add_argument('--conflict-strategy', type=str, default='skip',
                       choices=['skip', 'update', 'replace', 'version'],
                       help='Strategy for handling existing records')
    parser.add_argument('--validate-before-import', action='store_true',
                       help='Validate data before importing')

    # Status synchronization options
    parser.add_argument('--update-process-status', action='store_true', default=True,
                       help='Update ecod_schema.process_status after import (default: True)')
    parser.add_argument('--verify-sync', action='store_true', default=True,
                       help='Verify synchronization after import (default: True)')
    parser.add_argument('--update-batch-status', action='store_true',
                       help='Update batch completion status after import')

    args = parser.parse_args()

    # Set up logging
    logger = setup_logging(args.log_file, args.verbose)

    # Parse config
    try:
        config = parse_config(args.config)
    except Exception as e:
        logger.error(f"Error parsing config file: {str(e)}")
        sys.exit(1)

    # Convert strategy string to enum
    try:
        conflict_strategy = ConflictStrategy(args.conflict_strategy)
    except ValueError:
        logger.error(f"Invalid conflict strategy: {args.conflict_strategy}")
        sys.exit(1)

    logger.info(f"Starting import with conflict strategy: {conflict_strategy.value}")
    logger.info(f"Status tracking: {'enabled' if args.update_process_status else 'disabled'}")

    if args.repair_mode:
        logger.info("REPAIR MODE ENABLED - will process previously imported files")

    # Get database connection
    try:
        conn = get_db_connection(config)
    except Exception as e:
        logger.error(f"Error connecting to database: {str(e)}")
        sys.exit(1)

    try:
        conn.autocommit = False

        # Load batch information
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT id, batch_name, base_path
                FROM ecod_schema.batch
                """
            )
            batches = {row[0]: {"name": row[1], "base_path": row[2]} for row in cur.fetchall()}

        logger.info(f"Loaded information for {len(batches)} batches")

        # Get partition files with process ID tracking
        partition_files = get_partition_files_from_db(
            conn, args.batch_id, args.limit, not args.all_proteins, args.repair_mode
        )

        if args.repair_mode:
            previously_imported = sum(1 for f in partition_files if f[7])  # f[7] is previously_imported
            logger.info(f"Found {len(partition_files)} domain partition files "
                       f"({previously_imported} previously imported)")
        else:
            logger.info(f"Found {len(partition_files)} new domain partition files")

        if not partition_files:
            logger.info("No domain partition files to process. Exiting.")
            return

        # Process each file with enhanced status tracking
        stats = {
            'total_files': len(partition_files),
            'files_processed': 0,
            'files_skipped': 0,
            'files_not_found': 0,
            'files_invalid': 0,
            'partitions_imported': 0,
            'partitions_updated': 0,
            'partitions_skipped': 0,
            'domains_imported': 0,
            'evidence_imported': 0,
            'validation_errors': 0,
            'status_updates': 0
        }

        for file_data in partition_files:
            file_id, batch_id, file_type, file_path, pdb_id, chain_id, timestamp, previously_imported, process_id = file_data

            # Skip if we don't have batch information
            if batch_id not in batches:
                logger.warning(f"No batch information found for batch_id {batch_id}, skipping file_id {file_id}")
                stats['files_skipped'] += 1
                continue

            batch_info = batches[batch_id]

            # Resolve the file path
            full_path = resolve_file_path(batch_info["base_path"], file_path)

            # Check if file exists
            if not os.path.exists(full_path):
                logger.warning(f"File not found: {full_path}")
                stats['files_not_found'] += 1
                # Update status to show import failed due to missing file
                if args.update_process_status:
                    update_import_status(conn, pdb_id, chain_id, batch_id, process_id, False, "Domain file not found")
                continue

            # Read content from the file
            content = read_file_content(full_path)
            if not content:
                logger.warning(f"Failed to read content from file {full_path}")
                stats['files_skipped'] += 1
                if args.update_process_status:
                    update_import_status(conn, pdb_id, chain_id, batch_id, process_id, False, "Failed to read domain file")
                continue

            # Parse the XML content
            logger.debug(f"Processing file_id {file_id}: {pdb_id}_{chain_id} from batch {batch_id} "
                        f"({'previously imported' if previously_imported else 'new'}) process_id: {process_id}")

            partition = parse_partition_xml(content, file_id, pdb_id, chain_id, batch_id, timestamp)

            if not partition:
                logger.warning(f"Failed to parse partition for file_id {file_id}: {pdb_id}_{chain_id}")
                stats['files_invalid'] += 1
                if args.update_process_status:
                    update_import_status(conn, pdb_id, chain_id, batch_id, process_id, False, "Failed to parse domain XML")
                continue

            # Validate if requested
            if args.validate_before_import:
                validation_issues = validate_partition_data(partition)
                if validation_issues:
                    logger.warning(f"Validation issues for {pdb_id}_{chain_id}: {'; '.join(validation_issues)}")
                    stats['validation_errors'] += 1
                    if not args.dry_run:
                        if args.update_process_status:
                            update_import_status(conn, pdb_id, chain_id, batch_id, process_id, False,
                                               f"Validation failed: {'; '.join(validation_issues)}")
                        continue

            # Count evidence items
            evidence_count = sum(len(domain.get('evidence', [])) for domain in partition['domains'])

            # Import with conflict handling and status tracking
            try:
                protein_id, domain_ids = insert_domain_partition_with_strategy(
                    conn, partition, args.process_version, conflict_strategy,
                    process_id, args.update_process_status, args.dry_run
                )

                if not args.dry_run:
                    if protein_id:
                        if previously_imported and conflict_strategy != ConflictStrategy.SKIP:
                            stats['partitions_updated'] += 1
                        elif not previously_imported:
                            stats['partitions_imported'] += 1
                        else:
                            stats['partitions_skipped'] += 1

                        stats['domains_imported'] += len(domain_ids)
                        stats['evidence_imported'] += evidence_count

                        if args.update_process_status:
                            stats['status_updates'] += 1

                        logger.debug(f"Processed partition {protein_id} with {len(domain_ids)} domains")
                        conn.commit()
                    else:
                        stats['partitions_skipped'] += 1

                stats['files_processed'] += 1

                # Progress reporting
                if stats['files_processed'] % 50 == 0:
                    logger.info(f"Processed {stats['files_processed']} files, "
                               f"imported {stats['partitions_imported']} partitions, "
                               f"updated {stats['partitions_updated']} partitions, "
                               f"status updates: {stats['status_updates']}")

            except Exception as e:
                logger.error(f"Error importing partition for {pdb_id}_{chain_id}: {str(e)}")
                if args.verbose:
                    import traceback
                    logger.error(traceback.format_exc())
                stats['files_invalid'] += 1
                conn.rollback()
                continue

        # Verify synchronization if requested
        if args.verify_sync and args.batch_id:
            verify_import_synchronization(conn, args.batch_id, logger)

        # Update batch status if requested
        if args.update_batch_status and args.batch_id:
            update_batch_import_completion(conn, args.batch_id, logger)
            conn.commit()

        # Report final statistics
        logger.info("Enhanced Import Summary:")
        logger.info(f"  Total files: {stats['total_files']}")
        logger.info(f"  Files processed: {stats['files_processed']}")
        logger.info(f"  Files not found: {stats['files_not_found']}")
        logger.info(f"  Files invalid: {stats['files_invalid']}")
        logger.info(f"  Files skipped: {stats['files_skipped']}")
        logger.info(f"  Validation errors: {stats['validation_errors']}")
        logger.info(f"  Partitions imported: {stats['partitions_imported']}")
        logger.info(f"  Partitions updated: {stats['partitions_updated']}")
        logger.info(f"  Partitions skipped: {stats['partitions_skipped']}")
        logger.info(f"  Domains imported: {stats['domains_imported']}")
        logger.info(f"  Evidence items imported: {stats['evidence_imported']}")
        logger.info(f"  Process status updates: {stats['status_updates']}")

        if not args.dry_run:
            logger.info(f"Import completed successfully using {conflict_strategy.value} strategy with complete status synchronization")
        else:
            logger.info(f"DRY RUN completed successfully")

    except Exception as e:
        logger.error(f"Error during import: {str(e)}")
        if args.verbose:
            import traceback
            logger.error(traceback.format_exc())
        conn.rollback()
        sys.exit(1)

    finally:
        conn.close()


if __name__ == "__main__":
    main()
