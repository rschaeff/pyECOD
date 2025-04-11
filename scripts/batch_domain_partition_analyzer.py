#!/usr/bin/env python3
"""
batch_domain_partition_analyzer.py - Analyze protein domains for partition suitability

This script analyzes proteins in an ECOD batch to determine their suitability for 
blast-only partition vs. HHsearch processing. It parses domain summary XML files
to accurately assess BLAST hit distribution and generates detailed reports.
"""

import os
import sys
import logging
import argparse
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional, Set

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext

def setup_logging(verbose: bool = False, log_file: Optional[str] = None):
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
    
    return logging.getLogger("ecod.domain_partition")

class DomainPartitionAnalyzer:
    """Analyze proteins for partition suitability, focusing on BLAST hit coverage"""
    
    def __init__(self, config_path: str, peptide_threshold: int = 30):
        """Initialize analyzer with application context and settings
        
        Args:
            config_path: Path to ECOD configuration file
            peptide_threshold: Maximum residue length to consider as peptide
        """
        self.logger = logging.getLogger("ecod.domain_partition")
        self.config_path = config_path
        self.peptide_threshold = peptide_threshold
        
        # Initialize application context
        self.context = ApplicationContext(config_path)
        self.db = self.context.db_manager
        
        # Minimum BLAST hits required for confident partition
        self.min_blast_hits = 3
        self.hhsearch_probability_threshold = 90.0  # Minimum HHsearch probability
        
    def analyze_batch(self, batch_id: int, output_dir: Optional[str] = None):
        """Perform full analysis on a batch
        
        Args:
            batch_id: ECOD batch ID
            output_dir: Directory to save output files and visualizations
        """
        self.logger.info(f"Starting domain partition analysis for batch {batch_id}")
        
        # Get batch information
        batch_info = self._get_batch_info(batch_id)
        if not batch_info:
            self.logger.error(f"Batch {batch_id} not found")
            return
        
        # Get proteins and process status
        proteins = self._get_batch_proteins(batch_id)
        if not proteins:
            self.logger.error(f"No proteins found in batch {batch_id}")
            return
        
        self.logger.info(f"Found {len(proteins)} proteins in batch {batch_id}")
        
        # Get summary file paths
        summary_paths = self._get_domain_summary_paths(batch_id, proteins)
        self.logger.info(f"Found {len(summary_paths)} domain summary files")
        
        # Parse domain summary files
        self._parse_domain_summaries(proteins, summary_paths, batch_info['base_path'])
        
        # Analyze partition suitability
        results = self._analyze_partition_suitability(proteins)
        
        # Generate report
        self._print_analysis_report(batch_info, results)
        
        # Generate output files
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            self._generate_output_files(batch_id, results, output_dir)
            
        self.logger.info("Domain partition analysis completed")
        
        return results
    
    def _get_batch_info(self, batch_id: int) -> Optional[Dict[str, Any]]:
        """Get batch information from database
        
        Args:
            batch_id: Batch ID
            
        Returns:
            Dictionary with batch information or None if not found
        """
        query = """
        SELECT 
            id, batch_name, base_path, type, ref_version,
            total_items, completed_items, status,
            created_at, completed_at
        FROM 
            ecod_schema.batch
        WHERE 
            id = %s
        """
        
        try:
            result = self.db.execute_query(query, (batch_id,))
            
            if not result:
                return None
                
            return {
                'id': result[0][0],
                'name': result[0][1],
                'base_path': result[0][2],
                'type': result[0][3],
                'reference': result[0][4],
                'total_items': result[0][5],
                'completed_items': result[0][6],
                'status': result[0][7],
                'created_at': result[0][8],
                'completed_at': result[0][9]
            }
        except Exception as e:
            self.logger.error(f"Error retrieving batch information: {e}")
            return None
    
    def _get_batch_proteins(self, batch_id: int) -> List[Dict[str, Any]]:
        """Get protein information for the batch
        
        Args:
            batch_id: Batch ID
            
        Returns:
            List of protein dictionaries with metadata
        """
        query = """
        SELECT 
            p.id, p.pdb_id, p.chain_id, p.length,
            ps.id as process_id, ps.current_stage, ps.status
        FROM 
            ecod_schema.protein p
        JOIN 
            ecod_schema.process_status ps ON p.id = ps.protein_id
        WHERE 
            ps.batch_id = %s
        ORDER BY 
            p.pdb_id, p.chain_id
        """
        
        try:
            result = self.db.execute_query(query, (batch_id,))
            
            proteins = []
            for row in result:
                proteins.append({
                    'id': row[0],
                    'pdb_id': row[1],
                    'chain_id': row[2],
                    'length': row[3] or 0,
                    'process_id': row[4],
                    'current_stage': row[5],
                    'status': row[6],
                    'has_domain_blast_hits': False,
                    'has_chain_blast_hits': False,
                    'has_hhsearch_hits': False,
                    'domain_blast_hit_count': 0,
                    'chain_blast_hit_count': 0,
                    'hhsearch_hit_count': 0,
                    'best_domain_blast_evalue': 999.0,
                    'best_hhsearch_probability': 0.0,
                    'is_peptide': (row[3] or 0) < self.peptide_threshold if row[3] else False,
                    'partition_suitable': False,
                    'partition_confidence': 0.0,
                    'needs_hhsearch': False
                })
            
            return proteins
            
        except Exception as e:
            self.logger.error(f"Error retrieving batch proteins: {e}")
            return []
    
    def _get_domain_summary_paths(self, batch_id: int, proteins: List[Dict[str, Any]]) -> Dict[int, str]:
        """Get paths to domain summary files for proteins
        
        Args:
            batch_id: Batch ID
            proteins: List of protein dictionaries
            
        Returns:
            Dictionary mapping protein ID to summary file path
        """
        # Create lookup of process IDs to protein IDs
        process_to_protein = {p['process_id']: p['id'] for p in proteins}
        
        # Query for domain summary files
        query = """
        SELECT 
            pf.process_id, pf.file_path, pf.file_exists
        FROM 
            ecod_schema.process_file pf
        JOIN 
            ecod_schema.process_status ps ON pf.process_id = ps.id
        WHERE 
            ps.batch_id = %s AND pf.file_type = 'domain_summary'
            AND pf.file_exists = TRUE
        """
        
        try:
            result = self.db.execute_query(query, (batch_id,))
            
            summary_paths = {}
            for row in result:
                process_id = row[0]
                file_path = row[1]
                
                if process_id in process_to_protein:
                    protein_id = process_to_protein[process_id]
                    summary_paths[protein_id] = file_path
            
            return summary_paths
            
        except Exception as e:
            self.logger.error(f"Error retrieving domain summary files: {e}")
            return {}
    
    def _parse_domain_summaries(self, proteins: List[Dict[str, Any]], 
                               summary_paths: Dict[int, str],
                               base_path: str) -> None:
        """Parse domain summary XML files to extract hit information
        
        Args:
            proteins: List of protein dictionaries to update
            summary_paths: Dictionary mapping protein ID to summary file path
            base_path: Base directory path for resolving relative paths
        """
        self.logger.info("Parsing domain summary files...")
        
        # Create lookup table for faster protein access
        protein_lookup = {p['id']: p for p in proteins}
        
        # Track statistics
        processed = 0
        found_domain_hits = 0
        found_chain_hits = 0
        found_hhsearch_hits = 0
        errors = 0
        
        # Process each summary file
        for protein_id, file_path in summary_paths.items():
            try:
                # Resolve full path
                full_path = os.path.normpath(os.path.join(base_path, file_path))
                
                if not os.path.exists(full_path):
                    continue
                
                # Parse XML
                tree = ET.parse(full_path)
                root = tree.getroot()
                
                # Get blast_summ element
                blast_summ = root.find(".//blast_summ")
                if blast_summ is None:
                    continue
                
                # Track best scores/evalues
                min_evalue = 999.0
                max_probability = 0.0
                
                # Check for domain blast hits
                domain_blast = blast_summ.find("./blast_run")
                if domain_blast is not None:
                    hits = domain_blast.findall(".//hit")
                    if hits:
                        protein_lookup[protein_id]['has_domain_blast_hits'] = True
                        protein_lookup[protein_id]['domain_blast_hit_count'] = len(hits)
                        found_domain_hits += 1
                        
                        # Find best (lowest) e-value
                        for hit in hits:
                            evalues = hit.get("evalues", "").split(",")
                            if evalues and evalues[0]:
                                try:
                                    # Take the first e-value from each hit
                                    evalue = float(evalues[0])
                                    min_evalue = min(min_evalue, evalue)
                                except ValueError:
                                    pass
                
                # Check for chain blast hits
                chain_blast = blast_summ.find("./chain_blast_run")
                if chain_blast is not None:
                    hits = chain_blast.findall(".//hit")
                    if hits:
                        protein_lookup[protein_id]['has_chain_blast_hits'] = True
                        protein_lookup[protein_id]['chain_blast_hit_count'] = len(hits)
                        found_chain_hits += 1
                        
                        # Find best (lowest) e-value
                        for hit in hits:
                            evalues = hit.get("evalues", "").split(",")
                            if evalues and evalues[0]:
                                try:
                                    # Take the first e-value from each hit
                                    evalue = float(evalues[0])
                                    min_evalue = min(min_evalue, evalue)
                                except ValueError:
                                    pass
                
                # Check for HHsearch hits
                hhsearch = blast_summ.find("./hh_run")
                if hhsearch is not None:
                    hits = hhsearch.findall(".//hit")
                    if hits:
                        protein_lookup[protein_id]['has_hhsearch_hits'] = True
                        protein_lookup[protein_id]['hhsearch_hit_count'] = len(hits)
                        found_hhsearch_hits += 1
                        
                        # Find best (highest) probability
                        for hit in hits:
                            prob_str = hit.get("hh_prob", "0")
                            try:
                                probability = float(prob_str)
                                max_probability = max(max_probability, probability)
                            except ValueError:
                                pass
                
                # Update best scores
                protein_lookup[protein_id]['best_domain_blast_evalue'] = min_evalue
                protein_lookup[protein_id]['best_hhsearch_probability'] = max_probability
                
                processed += 1
                
                # Log progress periodically
                if processed % 1000 == 0:
                    self.logger.info(f"Processed {processed} summary files")
                    
            except Exception as e:
                self.logger.error(f"Error parsing summary file for protein {protein_id}: {e}")
                errors += 1
        
        self.logger.info(f"Parsed {processed} summary files")
        self.logger.info(f"Found: {found_domain_hits} with domain BLAST hits, "
                      f"{found_chain_hits} with chain BLAST hits, "
                      f"{found_hhsearch_hits} with HHSearch hits")
        
        if errors > 0:
            self.logger.warning(f"Encountered {errors} errors during parsing")
    
    def _analyze_partition_suitability(self, proteins: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze proteins for partition suitability
        
        Args:
            proteins: List of protein dictionaries with hit data
            
        Returns:
            Analysis results dictionary
        """
        # Initialize counters
        total_proteins = len(proteins)
        proteins_with_length = 0
        proteins_with_any_blast_hits = 0
        proteins_with_domain_blast_hits = 0
        proteins_with_chain_blast_hits = 0
        proteins_with_both_blast_hits = 0
        proteins_with_hhsearch_hits = 0
        
        peptides = 0
        peptides_with_blast_hits = 0
        non_peptides = 0
        non_peptides_with_blast_hits = 0
        non_peptides_without_blast_hits = 0
        
        suitable_for_partition = 0
        needs_hhsearch = 0
        
        # Calculate partition suitability for each protein
        for protein in proteins:
            # Check if length data is available
            if protein['length'] <= 0:
                continue
                
            proteins_with_length += 1
            
            # Calculate partition suitability
            has_blast_hits = protein['has_domain_blast_hits'] or protein['has_chain_blast_hits']
            protein['has_any_blast_hits'] = has_blast_hits
            
            if has_blast_hits:
                proteins_with_any_blast_hits += 1
            
            if protein['has_domain_blast_hits']:
                proteins_with_domain_blast_hits += 1
            
            if protein['has_chain_blast_hits']:
                proteins_with_chain_blast_hits += 1
            
            if protein['has_domain_blast_hits'] and protein['has_chain_blast_hits']:
                proteins_with_both_blast_hits += 1
            
            if protein['has_hhsearch_hits']:
                proteins_with_hhsearch_hits += 1
            
            # Check peptide status
            if protein['is_peptide']:
                peptides += 1
                if has_blast_hits:
                    peptides_with_blast_hits += 1
            else:
                non_peptides += 1
                if has_blast_hits:
                    non_peptides_with_blast_hits += 1
                else:
                    non_peptides_without_blast_hits += 1
            
            # Determine partition suitability
            if has_blast_hits and (
                protein['domain_blast_hit_count'] >= self.min_blast_hits or
                protein['chain_blast_hit_count'] >= self.min_blast_hits
            ):
                protein['partition_suitable'] = True
                suitable_for_partition += 1
                
                # Calculate confidence based on hit counts and e-values
                if protein['best_domain_blast_evalue'] < 1e-10:
                    confidence = 1.0  # Very high confidence
                elif protein['best_domain_blast_evalue'] < 1e-5:
                    confidence = 0.9  # High confidence
                elif protein['best_domain_blast_evalue'] < 1e-3:
                    confidence = 0.7  # Good confidence
                else:
                    confidence = 0.5  # Moderate confidence
                
                # Adjust confidence based on hit counts
                hit_count = max(protein['domain_blast_hit_count'], protein['chain_blast_hit_count'])
                if hit_count > 10:
                    confidence = min(1.0, confidence + 0.2)
                
                protein['partition_confidence'] = confidence
            else:
                # Not suitable for blast-only partition
                protein['partition_suitable'] = False
                protein['partition_confidence'] = 0.0
                
                # Check if protein needs HHsearch
                if not protein['is_peptide'] and not has_blast_hits:
                    protein['needs_hhsearch'] = True
                    needs_hhsearch += 1
        
        # Calculate percentages
        pct_with_any_blast = (proteins_with_any_blast_hits / proteins_with_length * 100) if proteins_with_length > 0 else 0
        pct_with_domain_blast = (proteins_with_domain_blast_hits / proteins_with_length * 100) if proteins_with_length > 0 else 0
        pct_with_chain_blast = (proteins_with_chain_blast_hits / proteins_with_length * 100) if proteins_with_length > 0 else 0
        pct_peptides = (peptides / proteins_with_length * 100) if proteins_with_length > 0 else 0
        
        pct_peptides_with_hits = (peptides_with_blast_hits / peptides * 100) if peptides > 0 else 0
        pct_non_peptides_with_hits = (non_peptides_with_blast_hits / non_peptides * 100) if non_peptides > 0 else 0
        pct_non_peptides_without_hits = (non_peptides_without_blast_hits / non_peptides * 100) if non_peptides > 0 else 0
        
        pct_suitable_for_partition = (suitable_for_partition / proteins_with_length * 100) if proteins_with_length > 0 else 0
        pct_needs_hhsearch = (needs_hhsearch / proteins_with_length * 100) if proteins_with_length > 0 else 0
        
        # Calculate average hit counts
        domain_blast_counts = [p['domain_blast_hit_count'] for p in proteins if p['has_domain_blast_hits']]
        chain_blast_counts = [p['chain_blast_hit_count'] for p in proteins if p['has_chain_blast_hits']]
        hhsearch_counts = [p['hhsearch_hit_count'] for p in proteins if p['has_hhsearch_hits']]
        
        avg_domain_blast_hits = sum(domain_blast_counts) / len(domain_blast_counts) if domain_blast_counts else 0
        avg_chain_blast_hits = sum(chain_blast_counts) / len(chain_blast_counts) if chain_blast_counts else 0
        avg_hhsearch_hits = sum(hhsearch_counts) / len(hhsearch_counts) if hhsearch_counts else 0
        
        # Get length statistics
        lengths = [p['length'] for p in proteins if p['length'] > 0]
        avg_length = sum(lengths) / len(lengths) if lengths else 0
        min_length = min(lengths) if lengths else 0
        max_length = max(lengths) if lengths else 0
        
        # Extract proteins that need different processing
        partition_suitable_proteins = [p for p in proteins if p['partition_suitable']]
        hhsearch_candidates = [p for p in proteins if p['needs_hhsearch']]
        
        # Collect results
        results = {
            'total_proteins': total_proteins,
            'proteins_with_length': proteins_with_length,
            'proteins_with_any_blast_hits': proteins_with_any_blast_hits,
            'proteins_with_domain_blast_hits': proteins_with_domain_blast_hits,
            'proteins_with_chain_blast_hits': proteins_with_chain_blast_hits,
            'proteins_with_both_blast_hits': proteins_with_both_blast_hits,
            'proteins_with_hhsearch_hits': proteins_with_hhsearch_hits,
            'peptides': peptides,
            'peptides_with_blast_hits': peptides_with_blast_hits,
            'non_peptides': non_peptides,
            'non_peptides_with_blast_hits': non_peptides_with_blast_hits,
            'non_peptides_without_blast_hits': non_peptides_without_blast_hits,
            'pct_with_any_blast': pct_with_any_blast,
            'pct_with_domain_blast': pct_with_domain_blast,
            'pct_with_chain_blast': pct_with_chain_blast,
            'pct_peptides': pct_peptides,
            'pct_peptides_with_hits': pct_peptides_with_hits,
            'pct_non_peptides_with_hits': pct_non_peptides_with_hits,
            'pct_non_peptides_without_hits': pct_non_peptides_without_hits,
            'avg_domain_blast_hits': avg_domain_blast_hits,
            'avg_chain_blast_hits': avg_chain_blast_hits,
            'avg_hhsearch_hits': avg_hhsearch_hits,
            'avg_length': avg_length,
            'min_length': min_length,
            'max_length': max_length,
            'suitable_for_partition': suitable_for_partition,
            'pct_suitable_for_partition': pct_suitable_for_partition,
            'needs_hhsearch': needs_hhsearch,
            'pct_needs_hhsearch': pct_needs_hhsearch,
            'partition_suitable_proteins': partition_suitable_proteins,
            'hhsearch_candidates': hhsearch_candidates,
            'all_proteins': proteins
        }
        
        return results
    
    def _print_analysis_report(self, batch_info: Dict[str, Any], results: Dict[str, Any]):
        """Print analysis report to console
        
        Args:
            batch_info: Batch information dictionary
            results: Analysis results dictionary
        """
        print("\n" + "="*80)
        print(f"BATCH DOMAIN PARTITION ANALYSIS REPORT - {batch_info['name']} (ID: {batch_info['id']})")
        print("="*80)
        
        print("\nBATCH INFORMATION:")
        print(f"  Name: {batch_info['name']}")
        print(f"  Type: {batch_info['type']}")
        print(f"  Reference: {batch_info['reference']}")
        print(f"  Total items: {batch_info['total_items']}")
        print(f"  Status: {batch_info['status']}")
        
        print("\nCHAIN STATISTICS:")
        print(f"  Total proteins in batch: {results['total_proteins']}")
        print(f"  Proteins with length data: {results['proteins_with_length']}")
        print(f"  Average chain length: {results['avg_length']:.1f} residues")
        print(f"  Length range: {results['min_length']} - {results['max_length']} residues")
        print(f"  Peptides (<{self.peptide_threshold} residues): {results['peptides']} ({results['pct_peptides']:.1f}%)")
        print(f"  Non-peptides (≥{self.peptide_threshold} residues): {results['non_peptides']} ({100-results['pct_peptides']:.1f}%)")
        
        print("\nBLAST HIT STATISTICS:")
        print(f"  Proteins with any BLAST hits: {results['proteins_with_any_blast_hits']} ({results['pct_with_any_blast']:.1f}%)")
        print(f"  Proteins with domain BLAST hits: {results['proteins_with_domain_blast_hits']} ({results['pct_with_domain_blast']:.1f}%)")
        print(f"  Proteins with chain BLAST hits: {results['proteins_with_chain_blast_hits']} ({results['pct_with_chain_blast']:.1f}%)")
        print(f"  Proteins with both types of hits: {results['proteins_with_both_blast_hits']}")
        
        if results['avg_domain_blast_hits'] > 0 or results['avg_chain_blast_hits'] > 0:
            print("\nBLAST HIT DETAILS:")
            if results['avg_domain_blast_hits'] > 0:
                print(f"  Average domain BLAST hits per protein: {results['avg_domain_blast_hits']:.1f}")
            if results['avg_chain_blast_hits'] > 0:
                print(f"  Average chain BLAST hits per protein: {results['avg_chain_blast_hits']:.1f}")
            if results['avg_hhsearch_hits'] > 0:
                print(f"  Average HHSearch hits per protein: {results['avg_hhsearch_hits']:.1f}")
        
        print("\nPEPTIDE ANALYSIS:")
        if results['peptides'] > 0:
            print(f"  Peptides with BLAST hits: {results['peptides_with_blast_hits']} ({results['pct_peptides_with_hits']:.1f}%)")
            print(f"  Peptides without BLAST hits: {results['peptides'] - results['peptides_with_blast_hits']} ({100-results['pct_peptides_with_hits']:.1f}%)")
        else:
            print("  No peptides found in this batch")
        
        print("\nNON-PEPTIDE ANALYSIS:")
        if results['non_peptides'] > 0:
            print(f"  Non-peptides with BLAST hits: {results['non_peptides_with_blast_hits']} ({results['pct_non_peptides_with_hits']:.1f}%)")
            print(f"  Non-peptides without BLAST hits: {results['non_peptides_without_blast_hits']} ({results['pct_non_peptides_without_hits']:.1f}%)")
        else:
            print("  No non-peptides found in this batch")
        
        print("\nPARTITION ANALYSIS:")
        print(f"  Proteins suitable for blast-only partition: {results['suitable_for_partition']} ({results['pct_suitable_for_partition']:.1f}%)")
        print(f"  Proteins requiring HHsearch analysis: {results['needs_hhsearch']} ({results['pct_needs_hhsearch']:.1f}%)")
        
        print("\nPARTITION RECOMMENDATION:")
        if results['pct_with_any_blast'] >= 50:
            print("  ✓ RECOMMENDED for --blast-only partition")
            print(f"    • More than 50% of chains ({results['pct_with_any_blast']:.1f}%) have BLAST hits")
        else:
            print("  ✗ NOT RECOMMENDED for --blast-only partition")
            print(f"    • Less than 50% of chains ({results['pct_with_any_blast']:.1f}%) have BLAST hits")
        
        if results['pct_needs_hhsearch'] >= 30:
            print("  ✓ RECOMMENDED for HHsearch pipeline")
            print(f"    • {results['pct_needs_hhsearch']:.1f}% of chains need HHsearch analysis")
            print(f"    • {results['needs_hhsearch']} non-peptide chains have no BLAST hits")
        else:
            print("  ✗ NOT necessary to run HHsearch pipeline")
            print(f"    • Only {results['pct_needs_hhsearch']:.1f}% of chains need HHsearch analysis")
        
        print("\n" + "="*80)
    
    def _generate_output_files(self, batch_id: int, results: Dict[str, Any], output_dir: str):
        """Generate output files and visualizations
        
        Args:
            batch_id: Batch ID
            results: Analysis results dictionary
            output_dir: Output directory path
        """
        self.logger.info(f"Generating output files in {output_dir}")
        
        # Make sure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # 1. Generate lists of proteins for different processing
        self._generate_protein_lists(batch_id, results, output_dir)
        
        # 2. Generate processing scripts
        self._generate_processing_scripts(batch_id, results, output_dir)
        
        # 3. Generate visualizations
        self._generate_visualizations(batch_id, results, output_dir)
        
        self.logger.info("Output files generated successfully")
    
    def _generate_protein_lists(self, batch_id: int, results: Dict[str, Any], output_dir: str):
        """Generate lists of proteins for different processing
        
        Args:
            batch_id: Batch ID
            results: Analysis results dictionary
            output_dir: Output directory path
        """
        # 1. List of proteins suitable for blast-only partition
        if results['partition_suitable_proteins']:
            # Sort by confidence (descending)
            partition_proteins = sorted(
                results['partition_suitable_proteins'],
                key=lambda p: p['partition_confidence'],
                reverse=True
            )
            
            with open(os.path.join(output_dir, f"batch_{batch_id}_blast_partition_proteins.txt"), 'w') as f:
                f.write(f"Proteins Suitable for Blast-Only Partition - Batch {batch_id}\n")
                f.write("=".ljust(80, "=") + "\n\n")
                f.write(f"Total proteins: {len(partition_proteins)}\n\n")
                
                f.write(f"{'PDB ID':<8} {'Chain':<6} {'Length':<8} {'BLAST Hits':<10} {'Confidence':<10}\n")
                f.write("-".ljust(50, "-") + "\n")
                
                for protein in partition_proteins:
                    hits = max(protein['domain_blast_hit_count'], protein['chain_blast_hit_count'])
                    f.write(f"{protein['pdb_id']:<8} {protein['chain_id']:<6} {protein['length']:<8} ")
                    f.write(f"{hits:<10} {protein['partition_confidence']:.2f}\n")
        
        # 2. List of proteins needing HHsearch
        if results['hhsearch_candidates']:
            # Sort by length (descending)
            hhsearch_proteins = sorted(
                results['hhsearch_candidates'],
                key=lambda p: p['length'],
                reverse=True
            )
            
            with open(os.path.join(output_dir, f"batch_{batch_id}_hhsearch_candidates.txt"), 'w') as f:
                f.write(f"Proteins Requiring HHsearch Analysis - Batch {batch_id}\n")
                f.write("=".ljust(80, "=") + "\n\n")
                f.write(f"Total candidates: {len(hhsearch_proteins)}\n\n")
                
                f.write(f"{'PDB ID':<8} {'Chain':<6} {'Length':<8} {'Current Stage':<20}\n")
                f.write("-".ljust(50, "-") + "\n")
                
                for protein in hhsearch_proteins:
                    f.write(f"{protein['pdb_id']:<8} {protein['chain_id']:<6} {protein['length']:<8} ")
                    f.write(f"{protein['current_stage']:<20}\n")
        
        # 3. List of peptides
        peptides = [p for p in results['all_proteins'] if p['is_peptide']]
        if peptides:
            # Sort by length (descending)
            peptides = sorted(peptides, key=lambda p: p['length'], reverse=True)
            
            with open(os.path.join(output_dir, f"batch_{batch_id}_peptides.txt"), 'w') as f:
                f.write(f"Peptides (Length < {self.peptide_threshold}) - Batch {batch_id}\n")
                f.write("=".ljust(80, "=") + "\n\n")
                f.write(f"Total peptides: {len(peptides)}\n\n")
                
                f.write(f"{'PDB ID':<8} {'Chain':<6} {'Length':<8} {'Has BLAST Hits':<15}\n")
                f.write("-".ljust(50, "-") + "\n")
                
                for protein in peptides:
                    has_hits = "Yes" if protein['has_any_blast_hits'] else "No"
                    f.write(f"{protein['pdb_id']:<8} {protein['chain_id']:<6} {protein['length']:<8} ")
                    f.write(f"{has_hits:<15}\n")
    
    def _generate_processing_scripts(self, batch_id: int, results: Dict[str, Any], output_dir: str):
        """Generate processing scripts for different protein groups
        
        Args:
            batch_id: Batch ID
            results: Analysis results dictionary
            output_dir: Output directory path
        """
        # 1. Script for running domain partition with blast-only
        with open(os.path.join(output_dir, f"batch_{batch_id}_run_blast_partition.sh"), 'w') as f:
            f.write("#!/bin/bash\n\n")
            f.write(f"# Script to run blast-only domain partition for Batch {batch_id}\n\n")
            
            f.write("# Set config path - modify if needed\n")
            f.write("CONFIG_PATH=\"config/config.yml\"\n\n")
            
            # Main command to run partition for the whole batch
            f.write("# Run domain analysis pipeline with blast-only\n")
            f.write(f"python run_domain_analysis.py --config $CONFIG_PATH \\\n")
            f.write(f"  --batch-id {batch_id} \\\n")
            f.write(f"  --blast-only \\\n")
            f.write(f"  --log-file logs/batch_{batch_id}_blast_partition.log\n\n")
            
            # Commands for individual proteins
            if results['partition_suitable_proteins']:
                f.write("# Alternatively, process individual proteins:\n")
                # Take first 5 proteins with highest confidence
                proteins = sorted(
                    results['partition_suitable_proteins'], 
                    key=lambda p: p['partition_confidence'],
                    reverse=True
                )[:5]
                
                for protein in proteins:
                    f.write(f"# python run_domain_analysis_single.py --config $CONFIG_PATH \\\n")
                    f.write(f"#   --protein-id {protein['id']} --batch-id {batch_id} --blast-only\n")
        
        # Make executable
        os.chmod(os.path.join(output_dir, f"batch_{batch_id}_run_blast_partition.sh"), 0o755)
        
        # 2. Script for running HHSearch analysis
        if results['hhsearch_candidates']:
            with open(os.path.join(output_dir, f"batch_{batch_id}_run_hhsearch.sh"), 'w') as f:
                f.write("#!/bin/bash\n\n")
                f.write(f"# Script to run HHsearch analysis for Batch {batch_id}\n\n")
                
                f.write("# Set config path - modify if needed\n")
                f.write("CONFIG_PATH=\"config/config.yml\"\n\n")
                
                # Main command for HHsearch pipeline
                f.write("# Run HHsearch pipeline for the batch\n")
                f.write(f"python run_hhsearch.py --config $CONFIG_PATH \\\n")
                f.write(f"  --batch-id {batch_id} \\\n")
                f.write(f"  --filter-no-blast-hits \\\n")
                f.write(f"  --min-length {self.peptide_threshold} \\\n")
                f.write(f"  --threads 8 \\\n")
                f.write(f"  --log-file logs/batch_{batch_id}_hhsearch.log\n\n")
                
                # Commands for individual proteins
                f.write("# Alternatively, process individual proteins:\n")
                # Take first 5 proteins with longest length
                proteins = sorted(
                    results['hhsearch_candidates'],
                    key=lambda p: p['length'],
                    reverse=True
                )[:5]
                
                for protein in proteins:
                    f.write(f"# python run_hhsearch_single.py --config $CONFIG_PATH \\\n")
                    f.write(f"#   --protein-id {protein['id']} --batch-id {batch_id}\n")
            
            # Make executable
            os.chmod(os.path.join(output_dir, f"batch_{batch_id}_run_hhsearch.sh"), 0o755)
    
    def _generate_visualizations(self, batch_id: int, results: Dict[str, Any], output_dir: str):
        """Generate visualizations of analysis results
        
        Args:
            batch_id: Batch ID
            results: Analysis results dictionary
            output_dir: Output directory path
        """
        # 1. Chain length distribution
        plt.figure(figsize=(10, 6))
        
        # Get lengths and create bins
        lengths = [p['length'] for p in results['all_proteins'] if p['length'] > 0]
        
        if lengths:
            # Create bins - use log scale for better visualization
            max_length = max(lengths)
            bins = np.logspace(np.log10(1), np.log10(max_length + 1), 50)
            
            # Split lengths by hit status
            lengths_with_hits = [p['length'] for p in results['all_proteins'] 
                               if p['length'] > 0 and p['has_any_blast_hits']]
            
            lengths_without_hits = [p['length'] for p in results['all_proteins'] 
                                  if p['length'] > 0 and not p['has_any_blast_hits']]
            
            # Plot histograms
            plt.hist(lengths, bins=bins, alpha=0.4, label='All proteins', color='blue')
            
            if lengths_with_hits:
                plt.hist(lengths_with_hits, bins=bins, alpha=0.5, label='With BLAST hits', color='green')
            
            if lengths_without_hits:
                plt.hist(lengths_without_hits, bins=bins, alpha=0.5, label='Without BLAST hits', color='red')
            
            # Add peptide threshold line
            plt.axvline(x=self.peptide_threshold, color='black', linestyle='--', 
                      label=f'Peptide threshold ({self.peptide_threshold})')
            
            plt.xscale('log')
            plt.xlabel('Chain Length (residues) - Log scale')
            plt.ylabel('Count')
            plt.title(f'Chain Length Distribution - Batch {batch_id}')
            plt.legend()
            plt.grid(alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'batch_{batch_id}_length_distribution.png'), dpi=300)
            plt.close()
        
        # 2. Partition category pie chart
        plt.figure(figsize=(10, 6))
        
        # Create data for pie chart
        labels = [
            'Non-peptides suitable for blast-only',
            'Non-peptides needing HHsearch',
            'Peptides with BLAST hits',
            'Peptides without BLAST hits'
        ]
        
        sizes = [
            results['suitable_for_partition'] - results['peptides_with_blast_hits'],
            results['needs_hhsearch'],
            results['peptides_with_blast_hits'],
            results['peptides'] - results['peptides_with_blast_hits']
        ]
        
        # Filter out empty categories
        non_zero_indices = [i for i, size in enumerate(sizes) if size > 0]
        filtered_labels = [labels[i] for i in non_zero_indices]
        filtered_sizes = [sizes[i] for i in non_zero_indices]
        
        if filtered_sizes:
            colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3']
            filtered_colors = [colors[i] for i in non_zero_indices]
            
            # Create explode array to highlight important categories
            explode = [0.1 if 'blast-only' in label or 'HHsearch' in label else 0 
                     for label in filtered_labels]
            
            plt.pie(filtered_sizes, explode=explode, labels=filtered_labels, colors=filtered_colors, 
                  autopct='%1.1f%%', shadow=True, startangle=90)
            plt.axis('equal')
            plt.title(f'Partition Categories - Batch {batch_id}')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'batch_{batch_id}_partition_categories.png'), dpi=300)
            plt.close()
        
        # 3. BLAST hit statistics
        plt.figure(figsize=(10, 6))
        
        categories = ['Any BLAST hits', 'Domain BLAST hits', 'Chain BLAST hits', 
                     'Suitable for partition', 'Needs HHsearch']
        
        counts = [
            results['proteins_with_any_blast_hits'],
            results['proteins_with_domain_blast_hits'],
            results['proteins_with_chain_blast_hits'],
            results['suitable_for_partition'],
            results['needs_hhsearch']
        ]
        
        # Calculate percentages
        percentages = [
            results['pct_with_any_blast'],
            results['pct_with_domain_blast'],
            results['pct_with_chain_blast'],
            results['pct_suitable_for_partition'],
            results['pct_needs_hhsearch']
        ]
        
        # Plot
        bars = plt.bar(categories, counts)
        
        # Add percentage labels
        for i, (bar, percentage) in enumerate(zip(bars, percentages)):
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                   f'{percentage:.1f}%',
                   ha='center', va='bottom', rotation=0)
        
        plt.xlabel('Category')
        plt.ylabel('Number of Proteins')
        plt.title(f'BLAST Hit and Partition Statistics - Batch {batch_id}')
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'batch_{batch_id}_hit_statistics.png'), dpi=300)
        plt.close()

def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='Analyze protein domains in a batch for partition suitability'
    )
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to analyze')
    parser.add_argument('--output-dir', type=str, 
                      help='Directory to save reports and visualizations')
    parser.add_argument('--peptide-threshold', type=int, default=30,
                      help='Maximum length to consider a chain as a peptide (default: 30)')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(args.verbose, args.log_file)
    
    # If no output directory specified, create one based on batch ID
    output_dir = args.output_dir
    if not output_dir:
        output_dir = f"batch_{args.batch_id}_analysis_{Path(__file__).stem}"
    
    logger.info(f"Starting domain partition analysis for batch {args.batch_id}")
    logger.info(f"Output will be saved to {output_dir}")
    
    # Initialize analyzer and run analysis
    analyzer = DomainPartitionAnalyzer(args.config, args.peptide_threshold)
    analyzer.analyze_batch(args.batch_id, output_dir)
    
    logger.info("Analysis completed successfully")

if __name__ == "__main__":
    main()