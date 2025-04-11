#!/usr/bin/env python3
"""
batch_domain_analyzer.py - Analyze proteins for blast-only partition suitability

This script analyzes a batch of proteins in the ECOD schema database to determine
which are suitable for blast-only partitioning and which require HHsearch. It examines
domain summary XML files for actual BLAST hits and generates reports and actionable
recommendations for domain processing.
"""

import os
import sys
import logging
import argparse
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple, Set
from datetime import datetime

# Import from ECOD framework
from ecod.core.context import ApplicationContext
from ecod.config import ConfigManager
from ecod.db.manager import DBManager
from ecod.exceptions import PipelineError

class BatchDomainAnalyzer:
    """Analyzer for determining blast-only partition suitability in ECOD batches"""
    
    def __init__(self, config_path: str, peptide_threshold: int = 30):
        """Initialize analyzer with application context
        
        Args:
            config_path: Path to ECOD configuration file
            peptide_threshold: Maximum residue length to consider a chain as peptide
        """
        self.config_path = config_path
        self.peptide_threshold = peptide_threshold
        
        # Initialize logger
        self.logger = logging.getLogger("ecod.domain_partition.analyzer")
        
        # Initialize application context
        self.context = ApplicationContext(config_path)
        self.db = self.context.db_manager
        
        # Initialize config
        self.config = ConfigManager(config_path).config
        
        # Set thresholds for analysis
        self.min_blast_hits = 3  # Minimum BLAST hits for confident partitioning
        self.min_hit_confidence = 70.0  # Minimum confidence (from HHsearch prob or BLAST e-value)
        
        # For tracking metrics on peptides, suitable proteins, etc.
        self.metrics = {}
        
    def analyze_batch(self, batch_id: int, output_dir: Optional[str] = None, 
                     include_hhsearch: bool = False) -> Dict[str, Any]:
        """Analyze a batch to determine domain partition strategy
        
        Args:
            batch_id: ID of the batch to analyze
            output_dir: Directory to save reports and visualizations
            include_hhsearch: Whether to include HHsearch results in analysis
            
        Returns:
            Analysis results dictionary
        """
        self.logger.info(f"Starting domain partition analysis for batch {batch_id}")
        
        # Get batch information
        batch_info = self._get_batch_info(batch_id)
        if not batch_info:
            self.logger.error(f"Batch {batch_id} not found in database")
            return {}
        
        # Get proteins in this batch
        proteins = self._get_batch_proteins(batch_id)
        if not proteins:
            self.logger.error(f"No proteins found in batch {batch_id}")
            return {}
            
        self.logger.info(f"Found {len(proteins)} proteins in batch {batch_id}")
        
        # Get domain summary file paths
        summary_files = self._get_domain_summary_files(batch_id, proteins)
        self.logger.info(f"Found {len(summary_files)} domain summary files")
        
        # Parse domain summaries to check for BLAST hits
        self._parse_domain_summaries(proteins, summary_files, batch_info['base_path'])
        
        # Perform partition suitability analysis
        analysis_results = self._analyze_partition_suitability(proteins, batch_info)
        
        # Generate reports and visualizations
        self._print_analysis_report(batch_info, analysis_results)
        
        if output_dir:
            self._generate_output_files(batch_id, analysis_results, output_dir)
        
        return analysis_results
    
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
            total_items, completed_items, status
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
                'status': result[0][7]
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
            ecod_schema.process_status ps ON ps.protein_id = p.id
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
    
    def _get_domain_summary_files(self, batch_id: int, proteins: List[Dict[str, Any]]) -> Dict[int, str]:
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
            pf.process_id, pf.file_path, pf.file_exists, p.pdb_id, p.chain_id
        FROM 
            ecod_schema.process_file pf
        JOIN 
            ecod_schema.process_status ps ON pf.process_id = ps.id
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        WHERE 
            ps.batch_id = %s 
            AND pf.file_type = 'domain_summary'
            AND pf.file_exists = TRUE
        """
        
        try:
            result = self.db.execute_query(query, (batch_id,))
            
            summary_files = {}
            for row in result:
                process_id = row[0]
                file_path = row[1]
                pdb_id = row[3]
                chain_id = row[4]
                
                if process_id in process_to_protein:
                    protein_id = process_to_protein[process_id]
                    summary_files[protein_id] = {
                        'path': file_path,
                        'pdb_id': pdb_id,
                        'chain_id': chain_id
                    }
            
            return summary_files
            
        except Exception as e:
            self.logger.error(f"Error retrieving domain summary files: {e}")
            return {}
    
    def _parse_domain_summaries(self, proteins: List[Dict[str, Any]], 
                               summary_files: Dict[int, Dict[str, str]], 
                               base_path: str) -> None:
        """Parse domain summary XML files to extract hit information
        
        Args:
            proteins: List of protein dictionaries to update
            summary_files: Dictionary mapping protein ID to summary file info
            base_path: Base directory path for resolving relative paths
        """
        self.logger.info("Parsing domain summary files for BLAST hits...")
        
        # Create protein ID lookup for faster access
        protein_lookup = {p['id']: p for p in proteins}
        
        # Counter for processed files
        processed_count = 0
        success_count = 0
        error_count = 0
        
        # Process each summary file
        for protein_id, file_info in summary_files.items():
            try:
                # Get file path
                file_path = file_info['path']
                
                # Construct full path
                full_path = os.path.normpath(os.path.join(base_path, file_path))
                
                if not os.path.exists(full_path):
                    self.logger.warning(f"Summary file not found: {full_path}")
                    continue
                
                # Parse XML
                tree = ET.parse(full_path)
                root = tree.getroot()
                
                # Get the blast_summ element
                blast_summ = root.find(".//blast_summ")
                if blast_summ is None:
                    continue
                
                # Check if it's marked as a peptide
                is_peptide = blast_summ.get("is_peptide", "false").lower() == "true"
                
                if is_peptide and protein_id in protein_lookup:
                    protein_lookup[protein_id]['is_peptide'] = True
                
                # Check for domain blast hits
                domain_blast_run = blast_summ.find("./blast_run")
                if domain_blast_run is not None:
                    hits = domain_blast_run.findall(".//hit")
                    hit_count = len(hits)
                    if hit_count > 0:
                        protein_lookup[protein_id]['has_domain_blast_hits'] = True
                        protein_lookup[protein_id]['domain_blast_hit_count'] = hit_count
                        
                        # Find best e-value
                        min_evalue = 999.0
                        for hit in hits:
                            evalues = hit.get("evalues", "").split(",")
                            if evalues and evalues[0]:
                                try:
                                    for evalue_str in evalues:
                                        if evalue_str and evalue_str.strip():
                                            evalue = float(evalue_str)
                                            if evalue < min_evalue:
                                                min_evalue = evalue
                                except ValueError:
                                    continue
                        
                        protein_lookup[protein_id]['best_domain_blast_evalue'] = min_evalue
                
                # Check for chain blast hits
                chain_blast_run = blast_summ.find("./chain_blast_run")
                if chain_blast_run is not None:
                    hits = chain_blast_run.findall(".//hit")
                    hit_count = len(hits)
                    if hit_count > 0:
                        protein_lookup[protein_id]['has_chain_blast_hits'] = True
                        protein_lookup[protein_id]['chain_blast_hit_count'] = hit_count
                
                # Check for HHSearch hits
                hh_run = blast_summ.find("./hh_run")
                if hh_run is not None:
                    hits = hh_run.findall(".//hit")
                    hit_count = len(hits)
                    if hit_count > 0:
                        protein_lookup[protein_id]['has_hhsearch_hits'] = True
                        protein_lookup[protein_id]['hhsearch_hit_count'] = hit_count
                        
                        # Find best probability
                        max_prob = 0.0
                        for hit in hits:
                            prob_str = hit.get("hh_prob", "0")
                            try:
                                prob = float(prob_str)
                                if prob > max_prob:
                                    max_prob = prob
                            except ValueError:
                                continue
                                
                        protein_lookup[protein_id]['best_hhsearch_probability'] = max_prob
                
                success_count += 1
            except Exception as e:
                self.logger.error(f"Error parsing summary file for protein {protein_id}: {e}")
                error_count += 1
            
            processed_count += 1
            
            # Log progress for large batches
            if processed_count % 1000 == 0:
                self.logger.info(f"Processed {processed_count} summary files")
        
        self.logger.info(f"Parsed {processed_count} summary files ({success_count} succeeded, {error_count} failed)")
    
    def _analyze_partition_suitability(self, proteins: List[Dict[str, Any]], 
                                     batch_info: Dict[str, Any]) -> Dict[str, Any]:
        """Analyze proteins for partition suitability
        
        Args:
            proteins: List of protein dictionaries
            batch_info: Batch information dictionary
            
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
        high_confidence_partition = 0
        medium_confidence_partition = 0
        low_confidence_partition = 0
        needs_hhsearch = 0
        
        # Classification criteria:
        # 1. Peptides (< peptide_threshold): Handle separately
        # 2. Non-peptides with sufficient BLAST hits: Suitable for blast-only partition
        # 3. Non-peptides without BLAST hits: Need HHsearch
        
        # Set classification for each protein
        for protein in proteins:
            # Skip if no length information
            if protein['length'] <= 0:
                continue
                
            proteins_with_length += 1
            
            # Check if protein has BLAST hits
            has_blast_hits = protein['has_domain_blast_hits'] or protein['has_chain_blast_hits']
            
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
                    
                    # Mark as needing HHsearch
                    protein['needs_hhsearch'] = True
                    needs_hhsearch += 1
            
            # Determine partition suitability based on BLAST hits
            hit_threshold = self.min_blast_hits
            
            if has_blast_hits and (
                protein['domain_blast_hit_count'] >= hit_threshold or
                protein['chain_blast_hit_count'] >= hit_threshold
            ):
                # Suitable for blast-only partition
                protein['partition_suitable'] = True
                suitable_for_partition += 1
                
                # Calculate confidence based on blast hit quality
                if protein['best_domain_blast_evalue'] < 1e-10:
                    confidence = 1.0  # High confidence
                    high_confidence_partition += 1
                elif protein['best_domain_blast_evalue'] < 1e-5:
                    confidence = 0.8  # Good confidence
                    high_confidence_partition += 1
                elif protein['best_domain_blast_evalue'] < 1e-3:
                    confidence = 0.6  # Medium confidence
                    medium_confidence_partition += 1
                else:
                    confidence = 0.4  # Low confidence
                    low_confidence_partition += 1
                
                # Adjust confidence based on hit counts
                hit_count = max(protein['domain_blast_hit_count'], protein['chain_blast_hit_count'])
                if hit_count > 10:
                    confidence = min(1.0, confidence + 0.2)
                
                protein['partition_confidence'] = confidence
        
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
        
        # Make partition recommendation
        recommended_for_blast_only = pct_with_any_blast >= 50
        recommended_for_hhsearch = pct_needs_hhsearch >= 30
        
        # Determine best approach based on analysis
        if recommended_for_blast_only and not recommended_for_hhsearch:
            recommended_approach = "blast_only"
        elif recommended_for_hhsearch and not recommended_for_blast_only:
            recommended_approach = "hhsearch"
        else:
            recommended_approach = "hybrid"
        
        # Split proteins into different categories
        partition_suitable_proteins = [p for p in proteins if p['partition_suitable']]
        hhsearch_candidates = [p for p in proteins if p['needs_hhsearch']]
        peptide_proteins = [p for p in proteins if p['is_peptide']]
        
        # Store overall metrics
        self.metrics = {
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
            'suitable_for_partition': suitable_for_partition,
            'pct_suitable_for_partition': pct_suitable_for_partition,
            'high_confidence_partition': high_confidence_partition,
            'medium_confidence_partition': medium_confidence_partition,
            'low_confidence_partition': low_confidence_partition,
            'needs_hhsearch': needs_hhsearch,
            'pct_needs_hhsearch': pct_needs_hhsearch,
            'recommended_for_blast_only': recommended_for_blast_only,
            'recommended_for_hhsearch': recommended_for_hhsearch,
            'recommended_approach': recommended_approach
        }
        
        # Collect full results
        results = {
            'batch_info': batch_info,
            'metrics': self.metrics,
            'partition_suitable_proteins': partition_suitable_proteins,
            'hhsearch_candidates': hhsearch_candidates,
            'peptide_proteins': peptide_proteins,
            'all_proteins': proteins
        }
        
        return results
    
    def _print_analysis_report(self, batch_info: Dict[str, Any], analysis: Dict[str, Any]):
        """Print analysis report to console
        
        Args:
            batch_info: Batch information dictionary
            analysis: Analysis results dictionary
        """
        metrics = analysis['metrics']
        
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
        print(f"  Total proteins in batch: {metrics['total_proteins']}")
        print(f"  Proteins with length data: {metrics['proteins_with_length']}")
        
        if metrics['proteins_with_length'] > 0:
            print(f"  Peptides (<{self.peptide_threshold} residues): {metrics['peptides']} ({metrics['pct_peptides']:.1f}%)")
            print(f"  Non-peptides (≥{self.peptide_threshold} residues): {metrics['non_peptides']} ({100-metrics['pct_peptides']:.1f}%)")
        
        print("\nBLAST HIT STATISTICS:")
        print(f"  Proteins with any BLAST hits: {metrics['proteins_with_any_blast_hits']} ({metrics['pct_with_any_blast']:.1f}%)")
        print(f"  Proteins with domain BLAST hits: {metrics['proteins_with_domain_blast_hits']} ({metrics['pct_with_domain_blast']:.1f}%)")
        print(f"  Proteins with chain BLAST hits: {metrics['proteins_with_chain_blast_hits']} ({metrics['pct_with_chain_blast']:.1f}%)")
        print(f"  Proteins with both types of hits: {metrics['proteins_with_both_blast_hits']}")
        
        print("\nPEPTIDE ANALYSIS:")
        if metrics['peptides'] > 0:
            print(f"  Peptides with BLAST hits: {metrics['peptides_with_blast_hits']} ({metrics['pct_peptides_with_hits']:.1f}%)")
            print(f"  Peptides without BLAST hits: {metrics['peptides'] - metrics['peptides_with_blast_hits']} "
                 f"({100-metrics['pct_peptides_with_hits']:.1f}%)")
        else:
            print("  No peptides found in this batch")
        
        print("\nNON-PEPTIDE ANALYSIS:")
        if metrics['non_peptides'] > 0:
            print(f"  Non-peptides with BLAST hits: {metrics['non_peptides_with_blast_hits']} ({metrics['pct_non_peptides_with_hits']:.1f}%)")
            print(f"  Non-peptides without BLAST hits: {metrics['non_peptides_without_blast_hits']} ({metrics['pct_non_peptides_without_hits']:.1f}%)")
        else:
            print("  No non-peptides found in this batch")
        
        print("\nPARTITION ANALYSIS:")
        print(f"  Proteins suitable for blast-only partition: {metrics['suitable_for_partition']} ({metrics['pct_suitable_for_partition']:.1f}%)")
        print(f"    • High confidence: {metrics['high_confidence_partition']}")
        print(f"    • Medium confidence: {metrics['medium_confidence_partition']}")
        print(f"    • Low confidence: {metrics['low_confidence_partition']}")
        print(f"  Proteins requiring HHsearch analysis: {metrics['needs_hhsearch']} ({metrics['pct_needs_hhsearch']:.1f}%)")
        
        print("\nPARTITION RECOMMENDATION:")
        
        if metrics['recommended_approach'] == "blast_only":
            print("  ✓ RECOMMENDED for --blast-only partition")
            print(f"    • More than 50% of chains ({metrics['pct_with_any_blast']:.1f}%) have BLAST hits")
            print("  ✗ NOT necessary to run HHsearch pipeline")
            print(f"    • Only {metrics['pct_needs_hhsearch']:.1f}% of chains need HHsearch analysis")
            
        elif metrics['recommended_approach'] == "hhsearch":
            print("  ✗ NOT RECOMMENDED for --blast-only partition")
            print(f"    • Less than 50% of chains ({metrics['pct_with_any_blast']:.1f}%) have BLAST hits")
            print("  ✓ RECOMMENDED for HHsearch pipeline")
            print(f"    • {metrics['pct_needs_hhsearch']:.1f}% of chains need HHsearch analysis")
            print(f"    • {metrics['needs_hhsearch']} chains have no BLAST hits and need more sensitive search")
            
        else:  # hybrid approach
            print("  ✓ RECOMMENDED for HYBRID approach (both blast-only and HHsearch)")
            print(f"    • {metrics['pct_suitable_for_partition']:.1f}% of chains suitable for blast-only partition")
            print(f"    • {metrics['pct_needs_hhsearch']:.1f}% of chains need HHsearch analysis")
            print(f"    • Consider running blast-only first, then HHsearch for remaining chains")
        
        print("\n" + "="*80)
    
    def _generate_output_files(self, batch_id: int, analysis: Dict[str, Any], output_dir: str):
        """Generate output files with analysis results
        
        Args:
            batch_id: Batch ID
            analysis: Analysis results dictionary
            output_dir: Output directory path
        """
        self.logger.info(f"Generating output files in {output_dir}")
        
        # Make sure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # 1. Generate lists of proteins for different processing
        self._generate_protein_lists(batch_id, analysis, output_dir)
        
        # 2. Generate processing scripts
        self._generate_processing_scripts(batch_id, analysis, output_dir)
        
        # 3. Generate visualizations
        self._generate_visualizations(batch_id, analysis, output_dir)
        
        self.logger.info("Output files generated successfully")
    
    def _generate_protein_lists(self, batch_id: int, analysis: Dict[str, Any], output_dir: str):
        """Generate lists of proteins for different processing
        
        Args:
            batch_id: Batch ID
            analysis: Analysis results dictionary
            output_dir: Output directory path
        """
        metrics = analysis['metrics']
        batch_info = analysis['batch_info']
        
        # 1. List of proteins suitable for blast-only partition
        partition_proteins = analysis['partition_suitable_proteins']
        if partition_proteins:
            # Sort by confidence (descending)
            partition_proteins = sorted(
                partition_proteins,
                key=lambda p: p['partition_confidence'],
                reverse=True
            )
            
            with open(os.path.join(output_dir, f"batch_{batch_id}_blast_partition_proteins.txt"), 'w') as f:
                f.write(f"Proteins Suitable for Blast-Only Partition - Batch {batch_id}\n")
                f.write("="*80 + "\n\n")
                f.write(f"Total proteins: {len(partition_proteins)}\n")
                f.write(f"Batch name: {batch_info['name']}\n")
                f.write(f"Reference: {batch_info['reference']}\n\n")
                
                f.write(f"{'PDB ID':<8} {'Chain':<6} {'Length':<8} {'BLAST Hits':<10} {'Confidence':<10}\n")
                f.write("-"*50 + "\n")
                
                for protein in partition_proteins:
                    hits = max(protein['domain_blast_hit_count'], protein['chain_blast_hit_count'])
                    f.write(f"{protein['pdb_id']:<8} {protein['chain_id']:<6} {protein['length']:<8} ")
                    f.write(f"{hits:<10} {protein['partition_confidence']:.2f}\n")
        
        # 2. List of proteins needing HHsearch
        hhsearch_proteins = analysis['hhsearch_candidates']
        if hhsearch_proteins:
            # Sort by length (descending)
            hhsearch_proteins = sorted(
                hhsearch_proteins,
                key=lambda p: p['length'],
                reverse=True
            )
            
            with open(os.path.join(output_dir, f"batch_{batch_id}_hhsearch_candidates.txt"), 'w') as f:
                f.write(f"Proteins Requiring HHsearch Analysis - Batch {batch_id}\n")
                f.write("="*80 + "\n\n")
                f.write(f"Total candidates: {len(hhsearch_proteins)}\n")
                f.write(f"Batch name: {batch_info['name']}\n")
                f.write(f"Reference: {batch_info['reference']}\n\n")
                
                f.write(f"{'PDB ID':<8} {'Chain':<6} {'Length':<8} {'Current Stage':<20}\n")
                f.write("-"*50 + "\n")
                
                for protein in hhsearch_proteins:
                    f.write(f"{protein['pdb_id']:<8} {protein['chain_id']:<6} {protein['length']:<8} ")
                    f.write(f"{protein['current_stage']:<20}\n")
        
        # 3. List of peptides
        peptide_proteins = analysis['peptide_proteins']
        if peptide_proteins:
            # Sort by length (descending)
            peptide_proteins = sorted(
                peptide_proteins,
                key=lambda p: p['length'],
                reverse=True
            )
            
            with open(os.path.join(output_dir, f"batch_{batch_id}_peptides.txt"), 'w') as f:
                f.write(f"Peptides (Length < {self.peptide_threshold}) - Batch {batch_id}\n")
                f.write("="*80 + "\n\n")
                f.write(f"Total peptides: {len(peptide_proteins)}\n")
                f.write(f"Batch name: {batch_info['name']}\n")
                f.write(f"Reference: {batch_info['reference']}\n\n")
                
                f.write(f"{'PDB ID':<8} {'Chain':<6} {'Length':<8} {'Has BLAST Hits':<15}\n")
                f.write("-"*50 + "\n")
                
                for protein in peptide_proteins:
                    has_hits = "Yes" if protein['has_domain_blast_hits'] or protein['has_chain_blast_hits'] else "No"
                    f.write(f"{protein['pdb_id']:<8} {protein['chain_id']:<6} {protein['length']:<8} ")
                    f.write(f"{has_hits:<15}\n")
    
    def _generate_processing_scripts(self, batch_id: int, analysis: Dict[str, Any], output_dir: str):
        """Generate processing scripts for different protein groups
        
        Args:
            batch_id: Batch ID
            analysis: Analysis results dictionary
            output_dir: Output directory path
        """
        metrics = analysis['metrics']
        batch_info = analysis['batch_info']
        
        # 1. Script for running domain partition with blast-only
        with open(os.path.join(output_dir, f"batch_{batch_id}_run_blast_partition.sh"), 'w') as f:
            f.write("#!/bin/bash\n\n")
            f.write(f"# Script to run blast-only domain partition for Batch {batch_id}\n")
            f.write(f"# Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("# Set configuration path - modify if needed\n")
            f.write("CONFIG_PATH=\"config/config.yml\"\n\n")
            
            # Main command to run partition for the whole batch
            f.write("# Run domain analysis pipeline with blast-only\n")
            f.write("python -m ecod.pipelines.domain_analysis.pipeline \\\n")
            f.write(f"  --config $CONFIG_PATH \\\n")
            f.write(f"  --batch-id {batch_id} \\\n")
            f.write(f"  --blast-only \\\n")
            f.write(f"  --log-file logs/batch_{batch_id}_blast_partition.log\n\n")
            
            # Alternative using the orchestrator
            f.write("# Alternative: use pipeline orchestrator\n")
            f.write("# python -m ecod.pipelines.orchestrator \\\n")
            f.write(f"#   --config $CONFIG_PATH \\\n")
            f.write(f"#   --run-domain-analysis \\\n")
            f.write(f"#   --batch-id {batch_id} \\\n")
            f.write(f"#   --blast-only\n\n")
            
            # Commands for individual proteins
            partition_proteins = analysis['partition_suitable_proteins']
            if partition_proteins:
                f.write("# Alternatively, process individual proteins:\n")
                # Take first 5 proteins with highest confidence
                proteins = sorted(
                    partition_proteins, 
                    key=lambda p: p['partition_confidence'],
                    reverse=True
                )[:5]
                
                for protein in proteins:
                    f.write(f"# python -m ecod.pipelines.domain_analysis.run_single \\\n")
                    f.write(f"#   --config $CONFIG_PATH \\\n")
                    f.write(f"#   --protein-id {protein['id']} \\\n")
                    f.write(f"#   --batch-id {batch_id} \\\n")
                    f.write(f"#   --blast-only\n")
        
        # Make executable
        os.chmod(os.path.join(output_dir, f"batch_{batch_id}_run_blast_partition.sh"), 0o755)
        
        # 2. Script for running HHSearch analysis
        hhsearch_proteins = analysis['hhsearch_candidates']
        if hhsearch_proteins:
            with open(os.path.join(output_dir, f"batch_{batch_id}_run_hhsearch.sh"), 'w') as f:
                f.write("#!/bin/bash\n\n")
                f.write(f"# Script to run HHsearch analysis for Batch {batch_id}\n")
                f.write(f"# Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                
                f.write("# Set configuration path - modify if needed\n")
                f.write("CONFIG_PATH=\"config/config.yml\"\n\n")
                
                # Main command for HHsearch pipeline
                f.write("# Run HHsearch pipeline for the batch\n")
                f.write("python -m ecod.pipelines.hhsearch_pipeline \\\n")
                f.write(f"  --config $CONFIG_PATH \\\n")
                f.write(f"  --batch-id {batch_id} \\\n")
                f.write(f"  --filter-no-blast-hits \\\n")
                f.write(f"  --min-length {self.peptide_threshold} \\\n")
                f.write(f"  --threads 8 \\\n")
                f.write(f"  --log-file logs/batch_{batch_id}_hhsearch.log\n\n")
                
                # Alternative using the orchestrator
                f.write("# Alternative: use pipeline orchestrator\n")
                f.write("# python -m ecod.pipelines.orchestrator \\\n")
                f.write(f"#   --config $CONFIG_PATH \\\n")
                f.write(f"#   --run-hhsearch \\\n")
                f.write(f"#   --batch-id {batch_id}\n\n")
                
                # Commands for individual proteins
                f.write("# Alternatively, process individual proteins:\n")
                # Take first 5 proteins with longest length
                proteins = sorted(
                    hhsearch_proteins,
                    key=lambda p: p['length'],
                    reverse=True
                )[:5]
                
                for protein in proteins:
                    f.write(f"# python -m ecod.pipelines.hhsearch_pipeline \\\n")
                    f.write(f"#   --config $CONFIG_PATH \\\n")
                    f.write(f"#   --protein-id {protein['id']} \\\n")
                    f.write(f"#   --batch-id {batch_id}\n")
            
            # Make executable
            os.chmod(os.path.join(output_dir, f"batch_{batch_id}_run_hhsearch.sh"), 0o755)
            
        # 3. Script for recommended hybrid approach
        if metrics['recommended_approach'] == "hybrid":
            with open(os.path.join(output_dir, f"batch_{batch_id}_run_hybrid.sh"), 'w') as f:
                f.write("#!/bin/bash\n\n")
                f.write(f"# Script to run hybrid approach (blast-only then HHsearch) for Batch {batch_id}\n")
                f.write(f"# Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                
                f.write("# Set configuration path - modify if needed\n")
                f.write("CONFIG_PATH=\"config/config.yml\"\n\n")
                
                f.write("# Step 1: Run blast-only domain partition\n")
                f.write("echo \"Running blast-only domain partition...\"\n")
                f.write("python -m ecod.pipelines.domain_analysis.pipeline \\\n")
                f.write(f"  --config $CONFIG_PATH \\\n")
                f.write(f"  --batch-id {batch_id} \\\n")
                f.write(f"  --blast-only \\\n")
                f.write(f"  --log-file logs/batch_{batch_id}_blast_partition.log\n\n")
                
                f.write("# Step 2: Get chains that still need HHsearch\n")
                f.write("echo \"Checking for chains that need HHsearch...\"\n")
                f.write("python -m ecod.pipelines.hhsearch_pipeline \\\n")
                f.write(f"  --config $CONFIG_PATH \\\n")
                f.write(f"  --batch-id {batch_id} \\\n")
                f.write(f"  --filter-no-blast-hits \\\n")
                f.write(f"  --min-length {self.peptide_threshold} \\\n")
                f.write(f"  --threads 8 \\\n")
                f.write(f"  --log-file logs/batch_{batch_id}_hhsearch.log\n\n")
                
                f.write("# Step 3: Run domain analysis on the complete batch\n")
                f.write("echo \"Running final domain analysis...\"\n")
                f.write("python -m ecod.pipelines.domain_analysis.pipeline \\\n")
                f.write(f"  --config $CONFIG_PATH \\\n")
                f.write(f"  --batch-id {batch_id} \\\n")
                f.write(f"  --log-file logs/batch_{batch_id}_domain_analysis.log\n\n")
                
                f.write("echo \"Hybrid approach completed\"\n")
                
            # Make executable
            os.chmod(os.path.join(output_dir, f"batch_{batch_id}_run_hybrid.sh"), 0o755)
    
    def _generate_visualizations(self, batch_id: int, analysis: Dict[str, Any], output_dir: str):
        """Generate visualizations of analysis results
        
        Args:
            batch_id: Batch ID
            analysis: Analysis results dictionary
            output_dir: Output directory path
        """
        metrics = analysis['metrics']
        batch_info = analysis['batch_info']
        proteins = analysis['all_proteins']
        
        # 1. Chain length distribution with blast hit status
        plt.figure(figsize=(10, 6))
        
        # Get lengths and hit status
        lengths = [p['length'] for p in proteins if p['length'] > 0]
        
        if lengths:
            # Create bins on log scale for better visualization
            max_length = max(lengths)
            bins = np.logspace(np.log10(1), np.log10(max_length + 1), 50)
            
            # Split lengths by hit status
            lengths_with_hits = [p['length'] for p in proteins 
                               if p['length'] > 0 and (p['has_domain_blast_hits'] or p['has_chain_blast_hits'])]
            
            lengths_without_hits = [p['length'] for p in proteins 
                                  if p['length'] > 0 and not (p['has_domain_blast_hits'] or p['has_chain_blast_hits'])]
            
            # Plot histograms
            if lengths_with_hits:
                plt.hist(lengths_with_hits, bins=bins, alpha=0.6, label='With BLAST hits', color='green')
            
            if lengths_without_hits:
                plt.hist(lengths_without_hits, bins=bins, alpha=0.6, label='Without BLAST hits', color='red')
            
            # Add peptide threshold line
            plt.axvline(x=self.peptide_threshold, color='black', linestyle='--', 
                      label=f'Peptide threshold ({self.peptide_threshold})')
            
            plt.xscale('log')
            plt.xlabel('Chain Length (residues) - Log scale')
            plt.ylabel('Count')
            plt.title(f'Chain Length Distribution by BLAST Hit Status - Batch {batch_id}')
            plt.grid(alpha=0.3)
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'batch_{batch_id}_length_distribution.png'), dpi=300)
            plt.close()
        
        # 2. Partition category pie chart
        plt.figure(figsize=(10, 6))
        
        # Get counts for each category
        suitable_count = len(analysis['partition_suitable_proteins'])
        hhsearch_count = len(analysis['hhsearch_candidates'])
        peptide_with_hits = sum(1 for p in analysis['peptide_proteins'] 
                              if p['has_domain_blast_hits'] or p['has_chain_blast_hits'])
        peptide_without_hits = len(analysis['peptide_proteins']) - peptide_with_hits
        
        # Create pie chart data
        labels = [
            'Non-peptides with BLAST hits',
            'Non-peptides needing HHsearch',
            'Peptides with BLAST hits',
            'Peptides without BLAST hits'
        ]
        
        sizes = [
            suitable_count,
            hhsearch_count,
            peptide_with_hits,
            peptide_without_hits
        ]
        
        # Remove empty categories
        non_zero_indices = [i for i, size in enumerate(sizes) if size > 0]
        filtered_labels = [labels[i] for i in non_zero_indices]
        filtered_sizes = [sizes[i] for i in non_zero_indices]
        
        if filtered_sizes:
            colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3']
            filtered_colors = [colors[i] for i in non_zero_indices]
            
            # Create explode array to highlight important segments
            explode = [0.1 if 'HHsearch' in label or 'with BLAST' in label else 0 
                     for label in filtered_labels]
            
            plt.pie(filtered_sizes, explode=explode, labels=filtered_labels, colors=filtered_colors, 
                  autopct='%1.1f%%', shadow=True, startangle=90)
            plt.axis('equal')
            plt.title(f'Partition Categories - Batch {batch_id}')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'batch_{batch_id}_partition_categories.png'), dpi=300)
            plt.close()
        
        # 3. Bar chart of partition recommendation
        plt.figure(figsize=(12, 6))
        
        categories = [
            'Any BLAST hits', 
            'Domain BLAST hits',
            'Chain BLAST hits',
            'Suitable for\nblast-only partition',
            'High confidence\nblast-only',
            'Needing HHsearch'
        ]
        
        values = [
            metrics['proteins_with_any_blast_hits'],
            metrics['proteins_with_domain_blast_hits'],
            metrics['proteins_with_chain_blast_hits'],
            metrics['suitable_for_partition'],
            metrics['high_confidence_partition'],
            metrics['needs_hhsearch']
        ]
        
        # Calculate percentages
        percentages = [
            metrics['pct_with_any_blast'],
            metrics['pct_with_domain_blast'],
            metrics['pct_with_chain_blast'],
            metrics['pct_suitable_for_partition'],
            (metrics['high_confidence_partition'] / metrics['proteins_with_length'] * 100) if metrics['proteins_with_length'] > 0 else 0,
            metrics['pct_needs_hhsearch']
        ]
        
        # Set different colors based on category
        colors = ['#729ECE', '#729ECE', '#729ECE', '#98df8a', '#2ca02c', '#ff9896']
        
        # Plot
        bars = plt.bar(categories, values, color=colors)
        
        # Add percentage labels
        for i, (bar, percentage) in enumerate(zip(bars, percentages)):
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                   f'{percentage:.1f}%',
                   ha='center', va='bottom', rotation=0)
        
        plt.xlabel('Category')
        plt.ylabel('Number of Proteins')
        plt.title(f'Domain Partition Analysis - Batch {batch_id}')
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'batch_{batch_id}_partition_analysis.png'), dpi=300)
        plt.close()
        
        # 4. Add a visual recommendation badge
        plt.figure(figsize=(10, 5))
        plt.axis('off')
        
        recommended_approach = metrics['recommended_approach']
        
        if recommended_approach == "blast_only":
            title = "RECOMMENDED: BLAST-ONLY PARTITION"
            color = "#2ca02c"  # Green
            details = [
                f"• {metrics['pct_with_any_blast']:.1f}% of chains have BLAST hits",
                f"• {metrics['suitable_for_partition']} chains suitable for blast-only partition",
                f"• Only {metrics['pct_needs_hhsearch']:.1f}% need HHsearch"
            ]
        elif recommended_approach == "hhsearch":
            title = "RECOMMENDED: HHSEARCH PIPELINE"
            color = "#d62728"  # Red
            details = [
                f"• Only {metrics['pct_with_any_blast']:.1f}% of chains have BLAST hits",
                f"• {metrics['needs_hhsearch']} chains need HHsearch",
                f"• {metrics['pct_needs_hhsearch']:.1f}% of chains have no BLAST hits"
            ]
        else:  # hybrid
            title = "RECOMMENDED: HYBRID APPROACH"
            color = "#ff7f0e"  # Orange
            details = [
                f"• {metrics['pct_suitable_for_partition']:.1f}% suitable for blast-only partition",
                f"• {metrics['pct_needs_hhsearch']:.1f}% need HHsearch",
                f"• Run blast-only first, then HHsearch for remaining chains"
            ]
        
        # Create recommendation box
        plt.text(0.5, 0.8, title, 
                size=20, weight="bold", ha="center", color="white",
                bbox=dict(facecolor=color, alpha=0.8, boxstyle="round,pad=0.5"))
        
        # Add details
        for i, detail in enumerate(details):
            plt.text(0.5, 0.6 - i*0.1, detail, size=14, ha="center")
        
        # Add batch info
        plt.text(0.5, 0.2, f"Batch {batch_id}: {batch_info['name']}", size=12, ha="center")
        plt.text(0.5, 0.15, f"Analysis date: {datetime.now().strftime('%Y-%m-%d')}", size=10, ha="center")
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'batch_{batch_id}_recommendation.png'), dpi=300)
        plt.close()

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

def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='Analyze proteins in a batch for blast-only partition suitability'
    )
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to analyze')
    parser.add_argument('--output-dir', type=str, 
                      help='Directory to save reports and visualizations')
    parser.add_argument('--peptide-threshold', type=int, default=30,
                      help='Maximum length to consider a chain as a peptide (default: 30)')
    parser.add_argument('--include-hhsearch', action='store_true',
                      help='Include HHsearch results in analysis')
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
        output_dir = f"batch_{args.batch_id}_domain_analysis_{datetime.now().strftime('%Y%m%d')}"
    
    logger.info(f"Starting domain partition analysis for batch {args.batch_id}")
    logger.info(f"Output will be saved to {output_dir}")
    
    # Initialize analyzer and run analysis
    analyzer = BatchDomainAnalyzer(args.config, args.peptide_threshold)
    results = analyzer.analyze_batch(args.batch_id, output_dir, args.include_hhsearch)
    
    if not results:
        logger.error("Analysis failed")
        return 1
    
    logger.info("Analysis completed successfully")
    return 0

if __name__ == "__main__":
    sys.exit(main())