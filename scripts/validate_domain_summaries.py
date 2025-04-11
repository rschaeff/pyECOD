#!/usr/bin/env python3
"""
validate_domain_summaries.py - Validate domain summaries against source BLAST files
"""

import os
import sys
import logging
import argparse
import xml.etree.ElementTree as ET
import concurrent.futures
from typing import Dict, List, Any, Optional, Tuple
from datetime import datetime

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext

class SummaryValidator:
    """Validates domain summaries against their source BLAST files"""
    
    def __init__(self, context, batch_id, batch_path, remove_invalid=False, target_proteins=None, dry_run=False):
        """Initialize validator with context and batch info"""
        self.context = context
        self.batch_id = batch_id
        self.batch_path = batch_path
        self.logger = logging.getLogger("ecod.validate_summaries")
        self.remove_invalid = remove_invalid
        self.target_proteins = target_proteins
        self.dry_run = dry_run
        
        # Validation stats
        self.valid_count = 0
        self.issues_count = 0
        self.no_chain_blast_count = 0
        self.no_domain_blast_count = 0
        self.empty_chain_blast_count = 0
        self.empty_domain_blast_count = 0
        self.removed_count = 0
        
        # Issues log
        self.issues_log = []
    
    def validate_protein(self, protein: tuple) -> Tuple[bool, List[str]]:
        """Validate a single protein's domain summary"""

        protein_id, pdb_id, chain_id, summary_path, chain_blast_path, domain_blast_path = protein
        
        # Skip if we have target proteins and this isn't one of them
        if self.target_proteins and f"{pdb_id}_{chain_id}" not in self.target_proteins:
            return True, []
        
        # Prepare absolute paths if needed
        summary_path = os.path.join(self.batch_path, summary_path) if not os.path.isabs(summary_path) else summary_path
        chain_blast_path = os.path.join(self.batch_path, chain_blast_path) if chain_blast_path and not os.path.isabs(chain_blast_path) else chain_blast_path
        domain_blast_path = os.path.join(self.batch_path, domain_blast_path) if domain_blast_path and not os.path.isabs(domain_blast_path) else domain_blast_path
        
        issues = []
        
        self.logger.debug(f"Validating {pdb_id}_{chain_id}...")
        
        # Check files exist
        if not os.path.exists(summary_path):
            issues.append(f"Domain summary file not found: {summary_path}")
            return False, issues
        
        # Track BLAST file existence
        if not chain_blast_path:
            self.no_chain_blast_count += 1
        if not domain_blast_path:
            self.no_domain_blast_count += 1
            
        # Parse summary file
        try:
            summary_tree = ET.parse(summary_path)
            summary_root = summary_tree.getroot()
            
            # Check basic attributes
            blast_summ = summary_root.find(".//blast_summ")
            if blast_summ is None:
                issues.append(f"No blast_summ element found in {summary_path}")
                if self.remove_invalid:
                    self._remove_file(summary_path)
                return False, issues
                
            # Check PDB ID and chain ID
            if blast_summ.get("pdb") != pdb_id or blast_summ.get("chain") != chain_id:
                issues.append(f"PDB/chain mismatch in summary: {blast_summ.get('pdb')}_{blast_summ.get('chain')} vs {pdb_id}_{chain_id}")
                if self.remove_invalid:
                    self._remove_file(summary_path)
                return False, issues
                
            # Validate chain BLAST hits
            if chain_blast_path and os.path.exists(chain_blast_path):
                chain_issues = self.validate_chain_blast(summary_root, chain_blast_path)
                issues.extend(chain_issues)
                
                # Track empty BLAST files
                if "Chain BLAST file has no hits" in " ".join(chain_issues):
                    self.empty_chain_blast_count += 1
                    
            # Validate domain BLAST hits
            if domain_blast_path and os.path.exists(domain_blast_path):
                domain_issues = self.validate_domain_blast(summary_root, domain_blast_path)
                issues.extend(domain_issues)
                
                # Track empty BLAST files
                if "Domain BLAST file has no hits" in " ".join(domain_issues):
                    self.empty_domain_blast_count += 1
            
            if not issues:
                return True, []
            else:
                if self.remove_invalid:
                    self._remove_file(summary_path, self.dry_run)
                return False, issues
                
        except Exception as e:
            issues.append(f"Error parsing XML: {str(e)}")
            if self.remove_invalid:
                self._remove_file(summary_path)
            return False, issues

    def _remove_file(self, file_path, dry_run=False):
        """Remove a file and update the database status if needed"""
        try:
            if os.path.exists(file_path):
                # Get relative path if needed
                rel_path = file_path
                if file_path.startswith(self.batch_path):
                    rel_path = os.path.relpath(file_path, self.batch_path)
                
                if dry_run:
                    self.logger.info(f"[DRY RUN] Would remove invalid summary file: {file_path}")
                    self.logger.info(f"[DRY RUN] Would update database record for {rel_path}")
                    self.removed_count += 1
                    return
                    
                self.logger.info(f"Removing invalid summary file: {file_path}")
                
                # Update the database first - set file_exists to FALSE
                update_query = """
                    UPDATE ecod_schema.process_file
                    SET file_exists = FALSE, last_checked = CURRENT_TIMESTAMP
                    FROM ecod_schema.process_status ps
                    WHERE process_file.process_id = ps.id
                      AND ps.batch_id = %s
                      AND process_file.file_path = %s
                      AND process_file.file_type = 'domain_summary'
                """
                
                # Execute the update query
                result = self.context.db.execute_query(update_query, (self.batch_id, rel_path), is_update=True)
                
                if result:
                    self.logger.debug(f"Updated database record for {rel_path}")
                else:
                    self.logger.warning(f"Failed to update database record for {rel_path}")
                
                # Now remove the file
                os.remove(file_path)
                self.removed_count += 1
                
        except Exception as e:
            self.logger.error(f"Error handling file {file_path}: {str(e)}")

    def validate_chain_blast(self, summary_root: ET.Element, chain_blast_path: str) -> List[str]:
        """Validate chain BLAST hits in summary against source file"""
        issues = []
        
        try:
            # Parse source BLAST file
            chain_blast_tree = ET.parse(chain_blast_path)
            chain_blast_root = chain_blast_tree.getroot()
            
            # Get summary hits
            chain_blast_run = summary_root.find(".//chain_blast_run")
            if chain_blast_run is None:
                # Check if no_chain_blast flag is set
                if summary_root.find(".//blast_summ").get("no_chain_blast") == "true":
                    # This is expected
                    return []
                else:
                    issues.append("Missing chain_blast_run element in summary")
                    return issues
            
            summary_hits = chain_blast_run.findall("./hits/hit")
            
            # Get source hits
            source_hits = chain_blast_root.findall(".//Hit")
            
            # Check for empty BLAST file
            if len(source_hits) == 0:
                # Check if chain_blast_no_hits flag is set
                if summary_root.find(".//blast_summ").get("chain_blast_no_hits") == "true":
                    # This is expected
                    issues.append("Chain BLAST file has no hits (correctly flagged)")
                    return []
                else:
                    issues.append("Chain BLAST file has no hits but chain_blast_no_hits flag not set")
                    return issues
            
            # Skip further validation if summary has no hits
            if len(summary_hits) == 0:
                if len(source_hits) > 0:
                    issues.append(f"Source BLAST has {len(source_hits)} hits but summary has none")
                return issues
                
            # Validate hit details
            # We don't expect exact match due to filtering, but check first few hits
            checked_hit_count = min(5, len(summary_hits))
            
            for i in range(checked_hit_count):
                summary_hit = summary_hits[i]
                hit_num = summary_hit.get("num")
                pdb_id = summary_hit.get("pdb_id", "")
                chain_id = summary_hit.get("chain_id", "")
                
                # Try to find matching hit in source
                matching_source_hit = None
                
                for source_hit in source_hits:
                    source_hit_num = source_hit.findtext("Hit_num", "")
                    
                    if source_hit_num == hit_num:
                        matching_source_hit = source_hit
                        break
                
                if matching_source_hit is None:
                    issues.append(f"Hit {hit_num} in summary not found in source BLAST")
                    continue
                
                # Check hit definition - should contain PDB ID and chain ID
                hit_def = matching_source_hit.findtext("Hit_def", "")
                if pdb_id and chain_id and (pdb_id not in hit_def or chain_id not in hit_def):
                    issues.append(f"Hit {hit_num} definition mismatch: {hit_def} vs {pdb_id}_{chain_id}")
                
                # Check query/hit regions exist
                query_regions = summary_hit.findtext("query_reg", "")
                hit_regions = summary_hit.findtext("hit_reg", "")
                
                if not query_regions:
                    issues.append(f"Missing query regions for hit {hit_num}")
                
                if not hit_regions:
                    issues.append(f"Missing hit regions for hit {hit_num}")
                    
            # Check BLAST program and version matches
            summary_program = chain_blast_run.get("program", "")
            source_program = chain_blast_root.findtext(".//BlastOutput_program", "")
            
            if summary_program and source_program and summary_program != source_program:
                issues.append(f"BLAST program mismatch: {summary_program} vs {source_program}")
                
            # Check database matches
            summary_db = chain_blast_run.findtext("blast_db", "")
            source_db = chain_blast_root.findtext(".//BlastOutput_db", "")
            
            if summary_db and source_db and not (summary_db in source_db or source_db in summary_db):
                issues.append(f"BLAST database mismatch: {summary_db} vs {source_db}")
                
        except Exception as e:
            issues.append(f"Error validating chain BLAST: {str(e)}")
        
        return issues

    def validate_domain_blast(self, summary_root: ET.Element, domain_blast_path: str) -> List[str]:
        """Validate domain BLAST hits in summary against source file"""
        issues = []
        
        try:
            # Parse source BLAST file
            domain_blast_tree = ET.parse(domain_blast_path)
            domain_blast_root = domain_blast_tree.getroot()
            
            # Get summary hits
            blast_run = summary_root.find(".//blast_run")
            if blast_run is None:
                # Check if no_domain_blast flag is set
                if summary_root.find(".//blast_summ").get("no_domain_blast") == "true":
                    # This is expected
                    return []
                else:
                    issues.append("Missing blast_run element in summary")
                    return issues
            
            summary_hits = blast_run.findall("./hits/hit")
            
            # Get source hits
            source_hits = domain_blast_root.findall(".//Hit")
            
            # Check for empty BLAST file
            if len(source_hits) == 0:
                # Check if domain_blast_no_hits flag is set
                if summary_root.find(".//blast_summ").get("domain_blast_no_hits") == "true":
                    # This is expected
                    issues.append("Domain BLAST file has no hits (correctly flagged)")
                    return []
                else:
                    issues.append("Domain BLAST file has no hits but domain_blast_no_hits flag not set")
                    return issues
            
            # Skip further validation if summary has no hits
            if len(summary_hits) == 0:
                if len(source_hits) > 0:
                    issues.append(f"Source BLAST has {len(source_hits)} hits but summary has none")
                return issues
                
            # Validate hit details
            # We don't expect exact match due to filtering, but check first few hits
            checked_hit_count = min(5, len(summary_hits))
            
            for i in range(checked_hit_count):
                summary_hit = summary_hits[i]
                hit_num = summary_hit.get("num")
                domain_id = summary_hit.get("domain_id", "")
                
                # Try to find matching hit in source
                matching_source_hit = None
                
                for source_hit in source_hits:
                    source_hit_num = source_hit.findtext("Hit_num", "")
                    
                    if source_hit_num == hit_num:
                        matching_source_hit = source_hit
                        break
                
                if matching_source_hit is None:
                    issues.append(f"Hit {hit_num} in summary not found in source BLAST")
                    continue
                
                # Check hit definition - should contain domain ID
                hit_def = matching_source_hit.findtext("Hit_def", "")
                if domain_id and domain_id != "NA" and domain_id not in hit_def:
                    issues.append(f"Hit {hit_num} definition mismatch: {hit_def} vs {domain_id}")
                
                # Check query/hit regions exist
                query_regions = summary_hit.findtext("query_reg", "")
                hit_regions = summary_hit.findtext("hit_reg", "")
                
                if not query_regions:
                    issues.append(f"Missing query regions for hit {hit_num}")
                
                if not hit_regions:
                    issues.append(f"Missing hit regions for hit {hit_num}")
                    
                # Check evalue exists
                evalues = summary_hit.get("evalues", "")
                if not evalues:
                    issues.append(f"Missing evalues for hit {hit_num}")
                    
            # Check BLAST program and version matches
            summary_program = blast_run.get("program", "")
            source_program = domain_blast_root.findtext(".//BlastOutput_program", "")
            
            if summary_program and source_program and summary_program != source_program:
                issues.append(f"BLAST program mismatch: {summary_program} vs {source_program}")
                
            # Check database matches
            summary_db = blast_run.findtext("blast_db", "")
            source_db = domain_blast_root.findtext(".//BlastOutput_db", "")
            
            if summary_db and source_db and not (summary_db in source_db or source_db in summary_db):
                issues.append(f"BLAST database mismatch: {summary_db} vs {source_db}")
                
        except Exception as e:
            issues.append(f"Error validating domain BLAST: {str(e)}")
        
        return issues

def validate_batch_summaries(batch_id: int, config_path: str, verbose: bool = False, 
                           threads: int = 1, output_file: Optional[str] = None,
                           remove_invalid: bool = False, target_proteins: Optional[List[str]] = None,
                           dry_run: bool = False):
    """Validate all domain summaries in a batch against their source BLAST files"""
    context = ApplicationContext(config_path)
    logger = logging.getLogger("ecod.validate_summaries")
    
    # Get batch path
    batch_query = "SELECT base_path FROM ecod_schema.batch WHERE id = %s"
    batch_result = context.db.execute_query(batch_query, (batch_id,))
    
    if not batch_result:
        logger.error(f"Batch {batch_id} not found")
        return 1
        
    batch_path = batch_result[0][0]
    
    # Query for all proteins in batch with domain summaries
    query = """
    SELECT 
        p.id, p.pdb_id, p.chain_id, pf.file_path as summary_path,
        pfc.file_path as chain_blast_path, pfd.file_path as domain_blast_path
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    JOIN 
        ecod_schema.process_file pf ON ps.id = pf.process_id AND pf.file_type = 'domain_summary'
    LEFT JOIN
        ecod_schema.process_file pfc ON ps.id = pfc.process_id AND pfc.file_type = 'chain_blast_result'
    LEFT JOIN
        ecod_schema.process_file pfd ON ps.id = pfd.process_id AND pfd.file_type = 'domain_blast_result'
    WHERE 
        ps.batch_id = %s AND pf.file_exists = TRUE
    ORDER BY p.pdb_id, p.chain_id
    """
    
    proteins = context.db.execute_query(query, (batch_id,))
    
    if not proteins:
        logger.error(f"No domain summaries found for batch {batch_id}")
        return 1
    
    target_count = len(target_proteins) if target_proteins else len(proteins)
    logger.info(f"Found {len(proteins)} proteins with domain summaries, targeting {target_count} for validation")
    
    # Initialize validator
    validator = SummaryValidator(context, batch_id, batch_path, remove_invalid, target_proteins, dry_run)
    
    # Process proteins
    start_time = datetime.now()
    
    if threads > 1:
        logger.info(f"Validating with {threads} threads")
        results = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            future_to_protein = {executor.submit(validator.validate_protein, protein): protein for protein in proteins}
            
            for i, future in enumerate(concurrent.futures.as_completed(future_to_protein)):
                protein = future_to_protein[future]
                pdb_id = protein[1]
                chain_id = protein[2]
                
                # Skip logging for proteins we're not targeting
                if target_proteins and f"{pdb_id}_{chain_id}" not in target_proteins:
                    continue
                
                try:
                    valid, issues = future.result()
                    if valid:
                        validator.valid_count += 1
                        logger.debug(f"Validation successful for {pdb_id}_{chain_id}")
                    else:
                        validator.issues_count += 1
                        issue_str = f"Issues for {pdb_id}_{chain_id}: {', '.join(issues)}"
                        validator.issues_log.append(issue_str)
                        logger.warning(issue_str)
                    
                    # Progress update
                    if (i + 1) % 100 == 0 or (i + 1) == len(proteins):
                        logger.info(f"Progress: {i+1}/{len(proteins)} proteins processed")
                        
                except Exception as e:
                    validator.issues_count += 1
                    issue_str = f"Error validating {pdb_id}_{chain_id}: {str(e)}"
                    validator.issues_log.append(issue_str)
                    logger.error(issue_str, exc_info=verbose)
    else:
        logger.info("Validating proteins sequentially")
        for i, protein in enumerate(proteins):
            pdb_id = protein[1]
            chain_id = protein[2]
            
            # Skip processing for proteins we're not targeting
            if target_proteins and f"{pdb_id}_{chain_id}" not in target_proteins:
                continue
                
            try:
                valid, issues = validator.validate_protein(protein)
                if valid:
                    validator.valid_count += 1
                    logger.debug(f"Validation successful for {pdb_id}_{chain_id}")
                else:
                    validator.issues_count += 1
                    issue_str = f"Issues for {pdb_id}_{chain_id}: {', '.join(issues)}"
                    validator.issues_log.append(issue_str)
                    logger.warning(issue_str)
                
                # Progress update
                if (i + 1) % 100 == 0 or (i + 1) == len(proteins):
                    logger.info(f"Progress: {i+1}/{len(proteins)} proteins processed")
                    
            except Exception as e:
                validator.issues_count += 1
                issue_str = f"Error validating {pdb_id}_{chain_id}: {str(e)}"
                validator.issues_log.append(issue_str)
                logger.error(issue_str, exc_info=verbose)
    
    end_time = datetime.now()
    duration = (end_time - start_time).total_seconds()
    
    # Generate report
    logger.info(f"Validation complete: {validator.valid_count} valid, {validator.issues_count} with issues")

    if dry_run and remove_invalid:
        logger.info(f"[DRY RUN] Would have removed {validator.removed_count} invalid summary files")
    elif remove_invalid: 
        logger.info(f"Removed {validator.removed_count} invalid summary files")

    logger.info(f"BLAST file stats:")
    logger.info(f"  - No chain BLAST records: {validator.no_chain_blast_count}")
    logger.info(f"  - No domain BLAST records: {validator.no_domain_blast_count}")
    logger.info(f"  - Empty chain BLAST files: {validator.empty_chain_blast_count}")
    logger.info(f"  - Empty domain BLAST files: {validator.empty_domain_blast_count}")
    logger.info(f"Total validation time: {duration:.2f} seconds")
    
    # Write issues to file if requested
    if output_file and validator.issues_log:
        try:
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            with open(output_file, 'w') as f:
                f.write(f"Domain summary validation issues for batch {batch_id}\n")
                f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"Total proteins: {len(proteins)}\n")
                f.write(f"Valid: {validator.valid_count}\n")
                f.write(f"With issues: {validator.issues_count}\n")
                if remove_invalid:
                    f.write(f"Removed: {validator.removed_count}\n")
                f.write("\n--- Issues ---\n")
                for issue in validator.issues_log:
                    f.write(f"{issue}\n")
            logger.info(f"Wrote issues to {output_file}")
        except Exception as e:
            logger.error(f"Error writing issues file: {str(e)}")
    
    return 0 if validator.issues_count == 0 else 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Validate domain summaries against source BLAST files')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to validate')
    parser.add_argument('--threads', type=int, default=1,
                      help='Number of threads to use for parallel validation')
    parser.add_argument('--output-file', type=str,
                      help='Write issues to this file')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    parser.add_argument('--remove-invalid', action='store_true',
                      help='Remove domain summary files that fail validation')
    parser.add_argument('--target-proteins', type=str,
                      help='Comma-separated list of protein IDs to validate (format: pdb_chain, e.g., 8gh6_P,8gh7_D)')
    parser.add_argument('--dry-run', action='store_true',
                      help='Do not actually remove files or update database, just report what would happen')
    
    
    args = parser.parse_args()
    
    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Parse target proteins if provided
    target_proteins = None
    if args.target_proteins:
        target_proteins = [p.strip() for p in args.target_proteins.split(',')]
        
    sys.exit(validate_batch_summaries(
        args.batch_id, 
        args.config, 
        args.verbose, 
        args.threads, 
        args.output_file,
        args.remove_invalid,
        target_proteins,
        args.dry_run
    ))