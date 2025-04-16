#!/usr/bin/env python3
"""
Diagnostic script for examining the 63 proteins stuck at domain_summary stage in batch 31.
This script checks file existence, integrity, and database synchronization issues.
"""

import os
import sys
import argparse
import logging
import xml.etree.ElementTree as ET
import psycopg2
from psycopg2.extras import DictCursor
import yaml
import json
from tabulate import tabulate

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("diagnose_batch31.log")
    ]
)
logger = logging.getLogger('fix_batch31')

class BatchFixer:
    """Tool for diagnosing domain processing issues"""
    
    def __init__(self, config):
        """Initialize the diagnostic tool"""
        self.config = config
        self.conn = None
        self.cursor = None
        self.batch_id = config['batch_id']
        
        # Initialize paths as None - will be loaded from database
        self.batch_dir = None
        self.data_dir = None
    
    def connect_db(self):
        """Connect to the database"""
        try:
            db_config = self.config.get('database', {})
            # Handle the case where 'database' is used instead of 'dbname'
            dbname = db_config.get('dbname', db_config.get('database', 'ecod_protein'))
            
            self.conn = psycopg2.connect(
                host=db_config.get('host', 'localhost'),
                port=db_config.get('port', 5432),
                dbname=dbname,
                user=db_config.get('user', 'ecod'),
                password=db_config.get('password', '')
            )
            self.cursor = self.conn.cursor(cursor_factory=DictCursor)
            logger.info("Connected to database")
            
            # Load batch information from database
            self._load_batch_info()
            
            return True
        except Exception as e:
            logger.error(f"Database connection error: {e}")
            return False
    
    def _load_batch_info(self):
        """Load batch information from database"""
        query = """
        SELECT base_path
        FROM ecod_schema.batch
        WHERE id = %s
        """
        try:
            self.cursor.execute(query, (self.batch_id,))
            result = self.cursor.fetchone()
            
            if result and result['base_path']:
                self.batch_dir = result['base_path']
                # Extract data_dir as the parent directory of base_path
                self.data_dir = os.path.dirname(os.path.dirname(self.batch_dir))
                logger.info(f"Loaded batch directory: {self.batch_dir}")
            else:
                logger.error(f"Batch {self.batch_id} not found in database or has no base_path")
                raise ValueError(f"Batch {self.batch_id} not found or has no base_path")
        except Exception as e:
            logger.error(f"Error loading batch information: {e}")
            raise
    
    def close_db(self):
        """Close the database connection"""
        if self.cursor:
            self.cursor.close()
        if self.conn:
            self.conn.close()
        logger.info("Database connection closed")
    
    def get_stuck_proteins(self):
        """Get proteins stuck at domain_summary stage"""
        query = """
        SELECT p.id, p.pdb_id, p.chain_id, p.length, ps.status, ps.current_stage, ps.id as process_id
        FROM ecod_schema.process_status ps
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        WHERE ps.batch_id = %s
        AND ps.current_stage = 'domain_summary'
        AND ps.status = 'success'
        ORDER BY p.pdb_id, p.chain_id
        """
        self.cursor.execute(query, (self.config['batch_id'],))
        return self.cursor.fetchall()
    
    def get_file_records(self, process_ids):
        """Get file records for processes"""
        if not process_ids:
            return []
        
        placeholders = ','.join(['%s'] * len(process_ids))
        query = f"""
        SELECT pf.process_id, p.pdb_id, p.chain_id, pf.file_type, pf.file_path
        FROM ecod_schema.process_file pf
        JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        WHERE ps.batch_id = %s
        AND pf.process_id IN ({placeholders})
        ORDER BY p.pdb_id, p.chain_id, pf.file_type
        """
        self.cursor.execute(query, (self.config['batch_id'],) + tuple(process_ids))
        return self.cursor.fetchall()
    
    def check_file_exists(self, file_path):
        """Check if a file exists in the filesystem"""
        if not file_path:
            return False
        
        # Handle both absolute and relative paths
        if os.path.isabs(file_path):
            full_path = file_path
        else:
            # Use data_dir from database
            full_path = os.path.join(self.data_dir, file_path)
        
        return os.path.exists(full_path)
    
    def validate_xml_file(self, file_path):
        """Validate XML file structure"""
        if not file_path:
            return False, "No file path provided"
        
        # Handle both absolute and relative paths
        if os.path.isabs(file_path):
            full_path = file_path
        else:
            # Use data_dir from database
            full_path = os.path.join(self.data_dir, file_path)
        
        if not os.path.exists(full_path):
            return False, f"File not found: {full_path}"
        
        try:
            tree = ET.parse(full_path)
            root = tree.getroot()
            return True, root.tag
        except ET.ParseError as e:
            return False, f"XML parsing error: {e}"
        except Exception as e:
            return False, f"Validation error: {e}"
    
    def check_expected_files(self, protein):
        """Check if all expected files exist for a protein"""
        pdb_id = protein['pdb_id']
        chain_id = protein['chain_id']
        
        # Define expected file paths
        expected_files = {
            'fasta': os.path.join(self.batch_dir, 'fastas', f"{pdb_id}_{chain_id}.fasta"),
            'chain_blast': os.path.join(self.batch_dir, 'blast', 'chain', f"{pdb_id}_{chain_id}.xml"),
            'domain_blast': os.path.join(self.batch_dir, 'blast', 'domain', f"{pdb_id}_{chain_id}.xml"),
            'domain_summary': os.path.join(self.batch_dir, 'domains', f"{pdb_id}_{chain_id}.domain_summary.xml")
        }
        
        results = {}
        for file_type, file_path in expected_files.items():
            exists = os.path.exists(file_path)
            
            if exists and file_path.endswith('.xml'):
                valid, message = self.validate_xml_file(file_path)
                results[file_type] = {
                    'path': file_path,
                    'exists': exists,
                    'valid': valid,
                    'message': message if not valid else None
                }
            else:
                results[file_type] = {
                    'path': file_path,
                    'exists': exists
                }
        
        return results
    
    def analyze_domain_summary(self, file_path):
        """Analyze domain summary file for content"""
        if not file_path or not os.path.exists(file_path):
            return {
                'status': 'missing',
                'message': 'File does not exist'
            }
        
        try:
            tree = ET.parse(file_path)
            root = tree.getroot()
            
            # Check root element
            if root.tag != 'blast_summ_doc':
                return {
                    'status': 'invalid',
                    'message': f'Unexpected root element: {root.tag}'
                }
            
            # Check for chain blast results
            chain_blast = root.find('.//chain_blast_run')
            chain_hits = []
            if chain_blast is not None:
                hits = chain_blast.findall('.//hit')
                for hit in hits:
                    hit_info = {
                        'id': hit.get('id'),
                        'query_regions': []
                    }
                    query_regions = hit.findall('.//query_reg')
                    for qr in query_regions:
                        range_str = qr.get('range') or qr.text
                        hit_info['query_regions'].append(range_str)
                    chain_hits.append(hit_info)
            
            # Check for domain blast results
            domain_blast = root.find('.//blast_run')
            domain_hits = []
            if domain_blast is not None:
                hits = domain_blast.findall('.//hit')
                for hit in hits:
                    hit_info = {
                        'id': hit.get('id'),
                        'query_regions': []
                    }
                    query_regions = hit.findall('.//query_reg')
                    for qr in query_regions:
                        range_str = qr.get('range') or qr.text
                        hit_info['query_regions'].append(range_str)
                    domain_hits.append(hit_info)
            
            return {
                'status': 'valid',
                'root_tag': root.tag,
                'chain_blast_hits': len(chain_hits),
                'domain_blast_hits': len(domain_hits),
                'chain_hits_details': chain_hits[:5],  # Limit to first 5 for brevity
                'domain_hits_details': domain_hits[:5]  # Limit to first 5 for brevity
            }
            
        except Exception as e:
            return {
                'status': 'error',
                'message': str(e)
            }
    
    def diagnose_protein(self, protein):
        """Run comprehensive diagnostic on a protein"""
        protein_id = protein['id']
        pdb_id = protein['pdb_id']
        chain_id = protein['chain_id']
        length = protein['length']
        process_id = protein['process_id']
        
        logger.info(f"Diagnosing protein {pdb_id}_{chain_id} (ID: {protein_id}, Process ID: {process_id}, Length: {length})")
        
        # Get file records from database
        self.cursor.execute("""
            SELECT file_type, file_path 
            FROM ecod_schema.process_file 
            WHERE process_id = %s
        """, (process_id,))
        
        file_records = {row['file_type']: row['file_path'] for row in self.cursor.fetchall()}
        
        # Check expected files
        expected_files = self.check_expected_files(protein)
        
        # Analyze domain summary if it exists
        domain_summary_path = expected_files.get('domain_summary', {}).get('path')
        domain_summary_analysis = None
        if domain_summary_path and os.path.exists(domain_summary_path):
            domain_summary_analysis = self.analyze_domain_summary(domain_summary_path)
        
        # Check database consistency
        db_consistent = True
        inconsistencies = []
        
        for file_type, details in expected_files.items():
            if details.get('exists') and file_type not in file_records:
                db_consistent = False
                inconsistencies.append(f"File exists but not in DB: {file_type}")
            elif not details.get('exists') and file_type in file_records:
                db_consistent = False
                inconsistencies.append(f"File in DB but doesn't exist: {file_type}")
        
        # Determine issue type
        issue_type = None
        issue_details = []
        
        # Check if this is a very small peptide
        if length is not None and length < 20:
            issue_type = "small_peptide"
            issue_details.append(f"Very small peptide: {length} residues")
        
        # Check if this is a very large protein
        if length is not None and length > 1500:
            issue_type = "large_protein"
            issue_details.append(f"Very large protein: {length} residues")
        
        # Check for missing domain summary
        if not expected_files.get('domain_summary', {}).get('exists'):
            issue_type = "missing_domain_summary"
            issue_details.append("Domain summary file is missing")
        
        # Check for invalid domain summary
        elif expected_files.get('domain_summary', {}).get('exists') and \
             not expected_files.get('domain_summary', {}).get('valid', True):
            issue_type = "invalid_domain_summary"
            issue_details.append(f"Domain summary XML is invalid: {expected_files['domain_summary'].get('message')}")
        
        # Check for empty BLAST results
        if expected_files.get('chain_blast', {}).get('exists') and \
           expected_files.get('chain_blast', {}).get('valid', True) and \
           domain_summary_analysis and \
           domain_summary_analysis.get('chain_blast_hits', 0) == 0:
            if not issue_type:
                issue_type = "empty_blast_results"
            issue_details.append("Chain BLAST results exist but contain no hits")
        
        # Check for database inconsistency
        if not db_consistent:
            if not issue_type:
                issue_type = "db_inconsistency"
            issue_details.append("Database and filesystem are inconsistent")
        
        # If no specific issue identified
        if not issue_type:
            issue_type = "unknown"
            issue_details.append("No specific issue identified")
        
        return {
            'protein_id': protein_id,
            'pdb_id': pdb_id,
            'chain_id': chain_id,
            'length': length,
            'process_id': process_id,
            'db_status': {
                'current_stage': protein['current_stage'],
                'status': protein['status']
            },
            'file_records': file_records,
            'expected_files': expected_files,
            'domain_summary_analysis': domain_summary_analysis,
            'db_consistent': db_consistent,
            'inconsistencies': inconsistencies,
            'issue_type': issue_type,
            'issue_details': issue_details
        }
    
    def run_diagnostics(self):
        """Run diagnostics on all stuck proteins"""
        # Get stuck proteins
        proteins = self.get_stuck_proteins()
        logger.info(f"Found {len(proteins)} proteins stuck at domain_summary stage")
        
        if not proteins:
            logger.info("No proteins to diagnose")
            return []
        
        # Run diagnostics on each protein
        results = []
        for protein in proteins:
            diagnostic = self.diagnose_protein(protein)
            results.append(diagnostic)
        
        return results
    
    def summarize_results(self, results):
        """Summarize diagnostic results"""
        if not results:
            return "No results to summarize"
        
        # Count issues by type
        issues_by_type = {}
        for result in results:
            issue_type = result['issue_type']
            if issue_type not in issues_by_type:
                issues_by_type[issue_type] = []
            issues_by_type[issue_type].append(result)
        
        # Generate summary
        summary = []
        summary.append(f"Total proteins diagnosed: {len(results)}")
        summary.append("\nIssue type breakdown:")
        
        for issue_type, proteins in issues_by_type.items():
            summary.append(f"  - {issue_type}: {len(proteins)} proteins")
        
        # Create a table of proteins with issues
        table_data = []
        for result in results:
            table_data.append([
                f"{result['pdb_id']}_{result['chain_id']}",
                result['length'],
                result['issue_type'],
                '; '.join(result['issue_details'])[:50] + ('...' if len('; '.join(result['issue_details'])) > 50 else '')
            ])
        
        summary.append("\nProtein issues summary:")
        summary.append(tabulate(table_data, headers=["Protein", "Length", "Issue Type", "Details"]))
        
        return '\n'.join(summary)
    
    def generate_recommendations(self, results):
        """Generate recommendations based on diagnostic results"""
        if not results:
            return "No results to generate recommendations from"
        
        recommendations = []
        recommendations.append("Recommendations based on diagnostic results:")
        
        # Group proteins by issue type
        issues_by_type = {}
        for result in results:
            issue_type = result['issue_type']
            if issue_type not in issues_by_type:
                issues_by_type[issue_type] = []
            issues_by_type[issue_type].append(result)
        
        # Generate recommendations for each issue type
        if 'small_peptide' in issues_by_type:
            count = len(issues_by_type['small_peptide'])
            recommendations.append(f"\n1. Handle {count} small peptides:")
            recommendations.append("   - Create single domain definitions for these proteins")
            recommendations.append("   - Consider lowering minimum domain length threshold")
            recommendations.append("   - Update database status to domain_partition_complete")
        
        if 'large_protein' in issues_by_type:
            count = len(issues_by_type['large_protein'])
            recommendations.append(f"\n2. Process {count} large proteins:")
            recommendations.append("   - Check for memory or timeout issues during processing")
            recommendations.append("   - Consider splitting processing into smaller segments")
            recommendations.append("   - Apply special handling for these exceptionally large proteins")
        
        if 'missing_domain_summary' in issues_by_type:
            count = len(issues_by_type['missing_domain_summary'])
            recommendations.append(f"\n3. Regenerate {count} missing domain summaries:")
            recommendations.append("   - Recreate domain summary files from BLAST results")
            recommendations.append("   - Update file paths in the database")
        
        if 'invalid_domain_summary' in issues_by_type:
            count = len(issues_by_type['invalid_domain_summary'])
            recommendations.append(f"\n4. Fix {count} invalid domain summaries:")
            recommendations.append("   - Apply XML structure corrections as noted in the Domain Summary Processing Improvements")
            recommendations.append("   - Regenerate these files with fixed XML structure")
        
        if 'empty_blast_results' in issues_by_type:
            count = len(issues_by_type['empty_blast_results'])
            recommendations.append(f"\n5. Handle {count} proteins with empty BLAST results:")
            recommendations.append("   - Create single domain definitions spanning the entire protein")
            recommendations.append("   - Consider re-running BLAST with different parameters")
        
        if 'db_inconsistency' in issues_by_type:
            count = len(issues_by_type['db_inconsistency'])
            recommendations.append(f"\n6. Fix {count} database inconsistencies:")
            recommendations.append("   - Align database records with filesystem state")
            recommendations.append("   - Add missing file records to the database")
            recommendations.append("   - Remove references to non-existent files")
        
        if 'unknown' in issues_by_type:
            count = len(issues_by_type['unknown'])
            recommendations.append(f"\n7. Investigate {count} proteins with unknown issues:")
            recommendations.append("   - Manually inspect domain summary files")
            recommendations.append("   - Check for any unusual characteristics in these proteins")
        
        # General recommendations
        recommendations.append("\nGeneral recommendations:")
        recommendations.append("1. Create a recovery script that:")
        recommendations.append("   - Processes each protein according to its issue type")
        recommendations.append("   - Updates database status to domain_partition_complete")
        recommendations.append("   - Creates domain partition files as needed")
        recommendations.append("2. Add validation checks to prevent similar issues in the future")
        recommendations.append("3. Consider adding automated recovery procedures for edge cases")
        
        return '\n'.join(recommendations)
    
    def export_results(self, results, output_dir):
        """Export diagnostic results to files"""
        os.makedirs(output_dir, exist_ok=True)
        
        # Export full results as JSON
        results_file = os.path.join(output_dir, "diagnostic_results.json")
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        # Export summary as text
        summary_file = os.path.join(output_dir, "diagnostic_summary.txt")
        with open(summary_file, 'w') as f:
            f.write(self.summarize_results(results))
        
        # Export recommendations as text
        recommendations_file = os.path.join(output_dir, "recommendations.txt")
        with open(recommendations_file, 'w') as f:
            f.write(self.generate_recommendations(results))
        
        # Export small peptides list
        small_peptides = [r for r in results if r['issue_type'] == 'small_peptide']
        small_peptides_file = os.path.join(output_dir, "small_peptides.txt")
        with open(small_peptides_file, 'w') as f:
            for protein in small_peptides:
                f.write(f"{protein['pdb_id']}_{protein['chain_id']}\t{protein['length']}\n")
        
        # Export large proteins list
        large_proteins = [r for r in results if r['issue_type'] == 'large_protein']
        large_proteins_file = os.path.join(output_dir, "large_proteins.txt")
        with open(large_proteins_file, 'w') as f:
            for protein in large_proteins:
                f.write(f"{protein['pdb_id']}_{protein['chain_id']}\t{protein['length']}\n")
        
        logger.info(f"Exported diagnostic results to {output_dir}")
        
        return {
            'results_file': results_file,
            'summary_file': summary_file,
            'recommendations_file': recommendations_file
        }

def load_config(config_path):
        """Load configuration from file"""
        try:
            # Load main config file
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)
            
            # Check for local config override
            config_dir = os.path.dirname(config_path)
            local_config_path = os.path.join(config_dir, 'config.local.yml')
            
            if os.path.exists(local_config_path):
                with open(local_config_path, 'r') as f:
                    local_config = yaml.safe_load(f)
                    
                # Merge configs, with local config taking precedence
                if local_config:
                    merge_configs(config, local_config)
            
            logger.info(f"Loaded configuration from {config_path} and local overrides")
            
            # Correction: the database name is 'database', not 'dbname'
            if 'database' in config.get('database', {}):
                config['database']['dbname'] = config['database'].get('database')
            
            return config
        except Exception as e:
            logger.error(f"Error loading configuration: {e}")
            raise

def merge_configs(config, local_config):
    """Recursively merge two configuration dictionaries"""
    for key, value in local_config.items():
        if key in config and isinstance(config[key], dict) and isinstance(value, dict):
            merge_configs(config[key], value)
        else:
            config[key] = value


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Diagnose proteins stuck at domain_summary stage")
    parser.add_argument("--batch-id", type=int, default=31, help="Batch ID to diagnose")
    parser.add_argument("--config", type=str, required=True, help="Configuration file path")
    parser.add_argument("--output-dir", type=str, default="diagnostic_results", help="Output directory for results")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main():
    """Main function"""
    args = parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    logger.info(f"Starting diagnostic for batch {args.batch_id}")
    
    # Load configuration
    config = load_config(args.config)
    config['batch_id'] = args.batch_id
    
    # Initialize diagnostic tool
    fixer = BatchFixer(config)
    
    # Connect to database
    if not fixer.connect_db():
        logger.error("Failed to connect to database")
        return 1
    
    try:
        # Run diagnostics
        results = fixer.run_diagnostics()
        
        # Generate and print summary
        summary = fixer.summarize_results(results)
        print("\n" + summary)
        
        # Generate and print recommendations
        recommendations = fixer.generate_recommendations(results)
        print("\n" + recommendations)
        
        # Export results
        fixer.export_results(results, args.output_dir)
        
    except Exception as e:
        logger.error(f"Error running diagnostics: {e}")
        return 1
    finally:
        # Close database connection
        fixer.close_db()
    
    logger.info("Diagnostics completed successfully")
    return 0


if __name__ == "__main__":
    sys.exit(main())