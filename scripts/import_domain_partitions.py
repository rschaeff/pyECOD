        #!/usr/bin/env python3
        """
        Import Domain Partitions Script (Revised)

        This script reads domain partition XML files referenced in the ecod_schema.process_file table
        and imports them into the pdb_analysis partition tables.

        Usage:
            python import_domain_partitions.py --config config.yml [--batch-id BATCH_ID] [--limit LIMIT] [--dry-run]

        Options:
            --config CONFIG      Path to configuration file
            --batch-id BATCH_ID  Process only a specific batch
            --limit LIMIT        Process only N files (for testing)
            --dry-run            Parse and print domains but don't write to database
            --verbose            Enable verbose output
            --log-file LOG_FILE  Path to log file (default: stdout)
            --process-version VERSION  Version identifier for this processing run
        """

        import os
        import sys
        import logging
        import argparse
        import datetime
        import xml.etree.ElementTree as ET
        import psycopg2
        from psycopg2.extras import execute_values
        import yaml
        from datetime import datetime


        def setup_logging(log_file=None, verbose=False):
            """Set up logging configuration."""
            level = logging.DEBUG if verbose else logging.INFO
            format_str = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

            if log_file:
                logging.basicConfig(filename=log_file, level=level, format=format_str)
            else:
                logging.basicConfig(level=level, format=format_str)

            return logging.getLogger(__name__)


        def parse_config(config_path):
            """
            Parse configuration file with optional local overrides.

            Looks for a config.local.yml file in the same directory as the
            main config file to override sensitive settings.
            """
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

                    # Deep merge the configurations
                    if local_config:
                        config = deep_merge(config, local_config)
                        logging.debug(f"Merged local configuration from {local_config_path}")

                return config
            except Exception as e:
                logging.error(f"Error parsing config file: {str(e)}")
                raise


        def deep_merge(dict1, dict2):
            """
            Deep merge two dictionaries. Values in dict2 override values in dict1.
            For nested dictionaries, performs a recursive merge instead of simple override.
            """
            result = dict1.copy()

            for key, value in dict2.items():
                if key in result and isinstance(result[key], dict) and isinstance(value, dict):
                    # If both values are dicts, merge them recursively
                    result[key] = deep_merge(result[key], value)
                else:
                    # Otherwise, use the value from dict2
                    result[key] = value

            return result


        def get_db_connection(config):
            """Create database connection from config with proper error handling."""
            db_config = config.get('database', {})

            try:
                # Log connection attempt (without password)
                conn_info = {
                    'host': db_config.get('host', 'dione'),
                    'port': db_config.get('port', 45000),
                    'dbname': db_config.get('name', 'ecod_protein'),
                    'user': db_config.get('user', 'ecod')
                }
                logging.debug(f"Connecting to database: {conn_info}")

                # Establish connection
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
            """
            Resolve a relative file path to its absolute location.

            Args:
                base_path: Base path of the batch
                file_path: Relative path from database (standardized to start with 'domains/')

            Returns:
                Absolute file path
            """
            # Handle edge cases where path might still be absolute
            if os.path.isabs(file_path):
                logging.warning(f"Found unexpected absolute path: {file_path}")
                return file_path

            # For standardized 'domains/' paths, make sure we don't duplicate directories
            if file_path.startswith('domains/'):
                # Extract the filename part (after 'domains/')
                filename_part = file_path[len('domains/'):]

                # Combine with base path's domains directory
                domains_dir = os.path.join(base_path, 'domains')
                return os.path.join(domains_dir, filename_part)

            # For other relative paths, simply join with base path
            return os.path.join(base_path, file_path)


        def read_file_content(file_path):
            """Read content from a file path."""
            try:
                with open(file_path, 'r', encoding='utf-8') as f:
                    return f.read()
            except Exception as e:
                logging.error(f"Error reading file {file_path}: {str(e)}")
                return None


        def get_partition_files_from_db(conn, batch_id=None, limit=None, rep_only=True):
            """Retrieve domain partition file paths from ecod_schema.process_file table."""
            with conn.cursor() as cur:
                query = """
                    SELECT
                        pf.id, ps.batch_id, pf.file_type, pf.file_path,
                        p.pdb_id, p.chain_id, pf.last_checked as timestamp
                    FROM
                        ecod_schema.process_file pf
                    JOIN
                        ecod_schema.process_status ps ON pf.process_id = ps.id
                    JOIN
                        ecod_schema.protein p ON ps.protein_id = p.id
                    WHERE
                        pf.file_type = 'domain_partition'
                        AND pf.file_exists = true
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
                                elif tag in ('coverage'):
                                    metadata[tag] = float(value)
                                elif tag in ('timestamp'):
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


        def find_protein_id(conn, pdb_id, chain_id):
            """Find the protein ID in the database for a given PDB ID and chain ID."""
            with conn.cursor() as cur:
                cur.execute(
                    """
                    SELECT id
                    FROM pdb_analysis.protein
                    WHERE pdb_id = %s AND chain_id = %s
                    """,
                    (pdb_id, chain_id)
                )
                result = cur.fetchone()
                if result:
                    return result[0]
                return None


        def extract_domain_sequence(conn, protein_id, start_pos, end_pos):
            """Extract domain sequence from protein sequence."""
            with conn.cursor() as cur:
                cur.execute(
                    """
                    SELECT ps.sequence
                    FROM pdb_analysis.protein_sequence ps
                    WHERE ps.protein_id = %s
                    """,
                    (protein_id,)
                )
                result = cur.fetchone()
                if not result or not result[0]:
                    return None

                seq = result[0]
                # Adjust for 0-based indexing if necessary
                seq_start = max(0, start_pos - 1)  # Convert 1-based to 0-based
                seq_end = min(len(seq), end_pos)  # End is exclusive in Python slicing

                if seq_start >= seq_end or seq_start >= len(seq):
                    return None

                return seq[seq_start:seq_end]


        def insert_domain_partition(conn, partition_data, process_version, dry_run=False):
            """Insert domain partition and its domains into the database."""
            if dry_run:
                logging.info(f"DRY RUN - Would insert partition for {partition_data['pdb_id']}_{partition_data['chain_id']} with {len(partition_data['domains'])} domains")
                for domain in partition_data['domains']:
                    logging.info(f"  Domain {domain['domain_id']}: {domain['range']} ({domain.get('t_group', 'NONE')})")
                    for evidence in domain['evidence']:
                        logging.info(f"    Evidence: {evidence['evidence_type']} from {evidence.get('source_id', 'unknown')}")
                return None, []

            protein_id = None
            domain_ids = []

            # Find the corresponding protein_id
            pdb_protein_id = find_protein_id(conn, partition_data['pdb_id'], partition_data['chain_id'])
            if not pdb_protein_id:
                logging.warning(f"No protein found for {partition_data['pdb_id']}_{partition_data['chain_id']} in pdb_analysis.protein")

            with conn.cursor() as cur:
                # Insert the partition protein
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
                        partition_data['pdb_id'],
                        partition_data['chain_id'],
                        partition_data['batch_id'],
                        partition_data['timestamp'],
                        partition_data['reference_version'],
                        partition_data['is_classified'],
                        partition_data['sequence_length'],
                        partition_data['coverage'],
                        partition_data['residues_assigned'],
                        partition_data['domains_with_evidence'],
                        partition_data['fully_classified_domains'],
                        partition_data['file_id'],
                        process_version
                    )
                )
                protein_id = cur.fetchone()[0]

                # Insert domains
                domain_values = []
                for domain in partition_data['domains']:
                    domain_values.append((
                        protein_id,
                        domain['domain_number'],
                        domain['domain_id'],
                        domain['start_pos'],
                        domain['end_pos'],
                        domain['range'],
                        domain['source'],
                        domain['source_id'],
                        domain['confidence'],
                        domain['t_group'],
                        domain['h_group'],
                        domain['x_group'],
                        domain['a_group'],
                        domain['is_manual_rep'],
                        domain['is_f70'],
                        domain['is_f40'],
                        domain['is_f99']
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

                    # If we found a protein in the database, try to extract domain sequences
                    domain_sequences = {}
                    if pdb_protein_id:
                        for i, domain in enumerate(partition_data['domains']):
                            sequence = extract_domain_sequence(
                                conn, pdb_protein_id, domain['start_pos'], domain['end_pos']
                            )
                            if sequence:
                                domain_sequences[i] = sequence

                    # Insert evidence for each domain
                    evidence_values = []
                    for i, domain_id in enumerate(domain_results):
                        domain_ids.append(domain_id)
                        for evidence in partition_data['domains'][i].get('evidence', []):
                            evidence_values.append((
                                domain_id,
                                evidence['evidence_type'],
                                evidence.get('source_id'),
                                evidence.get('domain_ref_id'),
                                evidence.get('hit_id'),
                                evidence.get('pdb_id'),
                                evidence.get('chain_id'),
                                evidence.get('confidence'),
                                evidence.get('probability'),
                                evidence.get('evalue'),
                                evidence.get('score'),
                                evidence.get('hsp_count'),
                                evidence.get('is_discontinuous'),
                                evidence.get('t_group'),
                                evidence.get('h_group'),
                                evidence.get('x_group'),
                                evidence.get('a_group'),
                                evidence.get('query_range'),
                                evidence.get('hit_range')
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

            return protein_id, domain_ids


        def main():
            """Main function."""
            parser = argparse.ArgumentParser(description='Import domain partitions from process_file table')
            parser.add_argument('--config', required=True, help='Path to configuration file')
            parser.add_argument('--batch-id', type=int, help='Process only a specific batch')
            parser.add_argument('--limit', type=int, help='Process only N files')
            parser.add_argument('--dry-run', action='store_true', help='Parse but don\'t write to database')
            parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
            parser.add_argument('--log-file', help='Path to log file')
            parser.add_argument('--process-version', default='1.0', help='Version identifier for this processing run')
            parser.add_argument('--all-proteins', action='store_true', help='Process all proteins, not just representatives')

            args = parser.parse_args()

            # Set up logging
            logger = setup_logging(args.log_file, args.verbose)

            # Parse config
            try:
                config = parse_config(args.config)
            except Exception as e:
                logger.error(f"Error parsing config file: {str(e)}")
                sys.exit(1)

            # Get database connection
            try:
                conn = get_db_connection(config)
            except Exception as e:
                logger.error(f"Error connecting to database: {str(e)}")
                sys.exit(1)

            try:
                conn.autocommit = False

                # Load all batch information upfront for path resolution
                with conn.cursor() as cur:
                    cur.execute(
                        """
                        SELECT id, batch_name, base_path
                        FROM ecod_schema.batch
                        """
                    )
                    batches = {row[0]: {"name": row[1], "base_path": row[2]} for row in cur.fetchall()}

                logger.info(f"Loaded information for {len(batches)} batches")

                # Get partition files from the database - only representative proteins by default
                partition_files = get_partition_files_from_db(conn, args.batch_id, args.limit, not args.all_proteins)
                logger.info(f"Found {len(partition_files)} domain partition files")

                if not partition_files:
                    logger.info("No domain partition files found. Exiting.")
                    return

                # Process each file
                total_partitions = 0
                total_domains = 0
                total_evidence = 0
                files_not_found = 0

                for file_data in partition_files:
                    file_id, batch_id, file_type, file_path, pdb_id, chain_id, timestamp = file_data

                    # Skip if we don't have batch information
                    if batch_id not in batches:
                        logger.warning(f"No batch information found for batch_id {batch_id}, skipping file_id {file_id}")
                        continue

                    batch_info = batches[batch_id]

                    # Resolve the file path using batch information
                    # Use the batch's base_path directly without additional components
                    full_path = resolve_file_path(batch_info["base_path"], file_path)

                    # Check if file exists
                    if not os.path.exists(full_path):
                        logger.warning(f"File not found: {full_path}")
                        files_not_found += 1
                        continue

                    # Read content from the file
                    content = read_file_content(full_path)
                    if not content:
                        logger.warning(f"Failed to read content from file {full_path}")
                        continue

                    # Parse the XML content
                    logger.debug(f"Processing file_id {file_id}: {pdb_id}_{chain_id} from batch {batch_id} ({batch_info['name']})")
                    partition = parse_partition_xml(content, file_id, pdb_id, chain_id, batch_id, timestamp)

                    if not partition:
                        logger.warning(f"Failed to parse partition for file_id {file_id}: {pdb_id}_{chain_id}")
                        continue

                    if not partition.get('domains'):
                        logger.warning(f"No domains found in partition for file_id {file_id}: {pdb_id}_{chain_id}")
                        continue

                    # Count evidence items
                    evidence_count = sum(len(domain.get('evidence', [])) for domain in partition['domains'])

                    # Insert into database
                    protein_id, domain_ids = insert_domain_partition(
                        conn, partition, args.process_version, args.dry_run
                    )

                    if not args.dry_run and protein_id:
                        total_partitions += 1
                        total_domains += len(domain_ids)
                        total_evidence += evidence_count

                        logger.debug(f"Inserted partition {protein_id} with {len(domain_ids)} domains")
                        conn.commit()

                    if total_partitions % 50 == 0 and total_partitions > 0:
                        logger.info(f"Processed {total_partitions} partitions, inserted {total_domains} domains")

                # Report summary statistics
                logger.info(f"Import summary:")
                logger.info(f"  Total files: {len(partition_files)}")
                logger.info(f"  Files not found: {files_not_found}")
                logger.info(f"  Partitions imported: {total_partitions}")
                logger.info(f"  Domains imported: {total_domains}")
                logger.info(f"  Evidence items imported: {total_evidence}")

                if not args.dry_run:
                    logger.info(f"Import completed successfully")
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
