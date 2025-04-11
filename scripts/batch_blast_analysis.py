#!/usr/bin/env python3
"""
batch_blast_analysis.py - Analyze batch proteins for blast-only partition suitability

This script analyzes the ECOD schema to find proteins suitable for the blast-only partitioning
approach by checking for BLAST hits in both domain and chain blast results. It also 
generates statistics on peptides and long chains with no hits (suitable for HHsearch).
"""

import os
import sys
import logging
import argparse
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from collections import Counter, defaultdict
from typing import Dict, List, Tuple, Any, Optional

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext
from ecod.config import ConfigManager

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

class BatchBlastAnalyzer:
    """Analyzer for batch proteins to determine blast-only partition suitability"""
    
    def __init__(self, config_path: str, peptide_threshold: int = 30):
        """Initialize analyzer with configuration
        
        Args:
            config_path: Path to configuration file
            peptide_threshold: Maximum length to consider a chain as a peptide
        """
        self.config_path = config_path
        self.peptide_threshold = peptide_threshold
        self.logger = logging.getLogger("ecod.batch_blast_analyzer")
        
        # Initialize application context
        self.context = ApplicationContext(config_path)
        self.db = self.context.db_manager
        
    def analyze_batch(self, batch_id: int, output_dir: Optional[str] = None):
        """Analyze a batch for blast-only partition suitability
        
        Args:
            batch_id: Batch ID to analyze
            output_dir: Directory to save reports and visualizations
        """
        self.logger.info(f"Analyzing batch {batch_id} for blast-only partition suitability")
        
        # Get batch info
        batch_info = self._get_batch_info(batch_id)
        if not batch_info:
            self.logger.error(f"Batch {batch_id} not found")
            return
            
        self.logger.info(f"Batch: {batch_info['name']} (ID: {batch_id})")
        
        # Get protein information for this batch
        proteins = self._get_batch_proteins(batch_id)
        if not proteins:
            self.logger.error(f"No proteins found in batch {batch_id}")
            return
            
        self.logger.info(f"Found {len(proteins)} proteins in batch {batch_id}")
        
        # Get BLAST hit information
        blast_results = self._get_blast_hit_info(batch_id, proteins)
        
        # Analyze proteins
        analysis_results = self._analyze_proteins(proteins, blast_results)
        
        # Print report
        self._print_analysis_report(batch_info, analysis_results)
        
        # Generate visualizations if output directory specified
        if output_dir:
            self._generate_visualizations(batch_id, analysis_results, output_dir)
    
    def _get_batch_info(self, batch_id: int) -> Optional[Dict[str, Any]]:
        """Get batch information
        
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
    
    def _get_batch_proteins(self, batch_id: int) -> List[Dict[str, Any]]:
        """Get all proteins in a batch
        
        Args:
            batch_id: Batch ID
            
        Returns:
            List of protein dictionaries
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
            p.id
        """
        
        result = self.db.execute_query(query, (batch_id,))
        
        proteins = []
        for row in result:
            proteins.append({
                'id': row[0],
                'pdb_id': row[1],
                'chain_id': row[2],
                'length': row[3] or 0,  # Handle NULL lengths
                'process_id': row[4],
                'current_stage': row[5],
                'status': row[6],
                'has_domain_blast_hits': False,
                'has_chain_blast_hits': False,
                'blast_hit_count': 0
            })
        
        return proteins
    
    def _get_blast_hit_info(self, batch_id: int, proteins: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
        """Get BLAST hit information for proteins in batch
        
        Args:
            batch_id: Batch ID
            proteins: List of protein dictionaries
            
        Returns:
            Dictionary mapping protein ID to blast hit info
        """
        # Create lookup of process IDs to protein IDs for faster processing
        process_to_protein = {p['process_id']: p['id'] for p in proteins}
        
        # Initialize results dictionary
        blast_results = {p['id']: {'domain_hits': 0, 'chain_hits': 0} for p in proteins}
        
        # Query for domain blast files
        domain_query = """
        SELECT 
            pf.process_id, pf.file_exists, pf.file_size
        FROM 
            ecod_schema.process_file pf
        JOIN 
            ecod_schema.process_status ps ON pf.process_id = ps.id
        WHERE 
            ps.batch_id = %s AND pf.file_type = 'domain_blast_result'
        """
        
        domain_result = self.db.execute_query(domain_query, (batch_id,))
        
        # Process domain blast results
        for row in domain_result:
            process_id = row[0]
            file_exists = row[1]
            file_size = row[2] or 0
            
            if process_id in process_to_protein:
                protein_id = process_to_protein[process_id]
                
                # A file that exists and has size > 1KB likely has hits
                # This is a heuristic - a more accurate approach would be to parse the XML
                if file_exists and file_size > 1024:
                    blast_results[protein_id]['domain_hits'] += 1
                    
                    # Update the protein dictionary as well
                    for p in proteins:
                        if p['id'] == protein_id:
                            p['has_domain_blast_hits'] = True
                            p['blast_hit_count'] += 1
                            break
        
        # Query for chain blast files
        chain_query = """
        SELECT 
            pf.process_id, pf.file_exists, pf.file_size
        FROM 
            ecod_schema.process_file pf
        JOIN 
            ecod_schema.process_status ps ON pf.process_id = ps.id
        WHERE 
            ps.batch_id = %s AND pf.file_type = 'chain_blast_result'
        """
        
        chain_result = self.db.execute_query(chain_query, (batch_id,))
        
        # Process chain blast results
        for row in chain_result:
            process_id = row[0]
            file_exists = row[1]
            file_size = row[2] or 0
            
            if process_id in process_to_protein:
                protein_id = process_to_protein[process_id]
                
                # A file that exists and has size > 1KB likely has hits
                if file_exists and file_size > 1024:
                    blast_results[protein_id]['chain_hits'] += 1
                    
                    # Update the protein dictionary as well
                    for p in proteins:
                        if p['id'] == protein_id:
                            p['has_chain_blast_hits'] = True
                            p['blast_hit_count'] += 1
                            break
        
        return blast_results
    
    def _analyze_proteins(self, proteins: List[Dict[str, Any]], 
                         blast_results: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze proteins for blast-only partition suitability
        
        Args:
            proteins: List of protein dictionaries
            blast_results: Dictionary mapping protein ID to blast hit info
            
        Returns:
            Dictionary with analysis results
        """
        # Initialize counters
        total_proteins = len(proteins)
        proteins_with_length = 0
        proteins_with_any_blast_hits = 0
        proteins_with_domain_blast_hits = 0
        proteins_with_chain_blast_hits = 0
        proteins_with_both_blast_hits = 0
        peptides = 0
        peptides_with_blast_hits = 0
        non_peptides = 0
        non_peptides_with_blast_hits = 0
        non_peptides_without_blast_hits = 0
        
        # Length statistics
        lengths = []
        length_with_hits = []
        length_without_hits = []
        
        # Process each protein
        for protein in proteins:
            length = protein['length']
            
            # Skip proteins without length information
            if length <= 0:
                continue
                
            proteins_with_length += 1
            lengths.append(length)
            
            # Check if protein has any blast hits
            has_blast_hits = protein['has_domain_blast_hits'] or protein['has_chain_blast_hits']
            if has_blast_hits:
                proteins_with_any_blast_hits += 1
                length_with_hits.append(length)
            else:
                length_without_hits.append(length)
            
            # Track specific blast hit types
            if protein['has_domain_blast_hits']:
                proteins_with_domain_blast_hits += 1
            
            if protein['has_chain_blast_hits']:
                proteins_with_chain_blast_hits += 1
            
            if protein['has_domain_blast_hits'] and protein['has_chain_blast_hits']:
                proteins_with_both_blast_hits += 1
            
            # Analyze peptides vs non-peptides
            is_peptide = length < self.peptide_threshold
            
            if is_peptide:
                peptides += 1
                if has_blast_hits:
                    peptides_with_blast_hits += 1
            else:
                non_peptides += 1
                if has_blast_hits:
                    non_peptides_with_blast_hits += 1
                else:
                    non_peptides_without_blast_hits += 1
        
        # Calculate percentages
        pct_with_any_blast = (proteins_with_any_blast_hits / proteins_with_length * 100) if proteins_with_length > 0 else 0
        pct_with_domain_blast = (proteins_with_domain_blast_hits / proteins_with_length * 100) if proteins_with_length > 0 else 0
        pct_with_chain_blast = (proteins_with_chain_blast_hits / proteins_with_length * 100) if proteins_with_length > 0 else 0
        pct_peptides = (peptides / proteins_with_length * 100) if proteins_with_length > 0 else 0
        pct_peptides_with_hits = (peptides_with_blast_hits / peptides * 100) if peptides > 0 else 0
        pct_non_peptides_with_hits = (non_peptides_with_blast_hits / non_peptides * 100) if non_peptides > 0 else 0
        pct_non_peptides_without_hits = (non_peptides_without_blast_hits / non_peptides * 100) if non_peptides > 0 else 0
        
        # Collect results
        results = {
            'total_proteins': total_proteins,
            'proteins_with_length': proteins_with_length,
            'proteins_with_any_blast_hits': proteins_with_any_blast_hits,
            'proteins_with_domain_blast_hits': proteins_with_domain_blast_hits,
            'proteins_with_chain_blast_hits': proteins_with_chain_blast_hits,
            'proteins_with_both_blast_hits': proteins_with_both_blast_hits,
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
            'lengths': lengths,
            'length_with_hits': length_with_hits,
            'length_without_hits': length_without_hits,
            'proteins': proteins  # Include the full protein list for detailed analysis
        }
        
        return results
    
    def _print_analysis_report(self, batch_info: Dict[str, Any], analysis: Dict[str, Any]):
        """Print analysis report
        
        Args:
            batch_info: Batch information dictionary
            analysis: Analysis results dictionary
        """
        print("\n" + "="*80)
        print(f"BATCH BLAST ANALYSIS REPORT - {batch_info['name']} (ID: {batch_info['id']})")
        print("="*80)
        
        print("\nBATCH INFORMATION:")
        print(f"  Name: {batch_info['name']}")
        print(f"  Type: {batch_info['type']}")
        print(f"  Reference: {batch_info['reference']}")
        print(f"  Total items: {batch_info['total_items']}")
        print(f"  Status: {batch_info['status']}")
        
        print("\nCHAIN STATISTICS:")
        print(f"  Total proteins in batch: {analysis['total_proteins']}")
        print(f"  Proteins with length data: {analysis['proteins_with_length']}")
        
        if analysis['proteins_with_length'] > 0:
            avg_length = sum(analysis['lengths']) / len(analysis['lengths']) if analysis['lengths'] else 0
            min_length = min(analysis['lengths']) if analysis['lengths'] else 0
            max_length = max(analysis['lengths']) if analysis['lengths'] else 0
            
            print(f"  Average chain length: {avg_length:.1f} residues")
            print(f"  Length range: {min_length} - {max_length} residues")
            print(f"  Peptides (<{self.peptide_threshold} residues): {analysis['peptides']} ({analysis['pct_peptides']:.1f}%)")
            print(f"  Non-peptides (≥{self.peptide_threshold} residues): {analysis['non_peptides']} ({100-analysis['pct_peptides']:.1f}%)")
        
        print("\nBLAST HIT STATISTICS:")
        print(f"  Proteins with any BLAST hits: {analysis['proteins_with_any_blast_hits']} ({analysis['pct_with_any_blast']:.1f}%)")
        print(f"  Proteins with domain BLAST hits: {analysis['proteins_with_domain_blast_hits']} ({analysis['pct_with_domain_blast']:.1f}%)")
        print(f"  Proteins with chain BLAST hits: {analysis['proteins_with_chain_blast_hits']} ({analysis['pct_with_chain_blast']:.1f}%)")
        print(f"  Proteins with both types of hits: {analysis['proteins_with_both_blast_hits']}")
        
        print("\nPEPTIDE ANALYSIS:")
        if analysis['peptides'] > 0:
            print(f"  Peptides with BLAST hits: {analysis['peptides_with_blast_hits']} ({analysis['pct_peptides_with_hits']:.1f}%)")
            print(f"  Peptides without BLAST hits: {analysis['peptides'] - analysis['peptides_with_blast_hits']} ({100-analysis['pct_peptides_with_hits']:.1f}%)")
        else:
            print("  No peptides found in this batch")
        
        print("\nNON-PEPTIDE ANALYSIS:")
        if analysis['non_peptides'] > 0:
            print(f"  Non-peptides with BLAST hits: {analysis['non_peptides_with_blast_hits']} ({analysis['pct_non_peptides_with_hits']:.1f}%)")
            print(f"  Non-peptides without BLAST hits: {analysis['non_peptides_without_blast_hits']} ({analysis['pct_non_peptides_without_hits']:.1f}%)")
        else:
            print("  No non-peptides found in this batch")
        
        print("\nPARTITION RECOMMENDATION:")
        if analysis['pct_with_any_blast'] >= 50:
            print("  ✓ RECOMMENDED for --blast-only partition")
            print(f"    • More than 50% of chains ({analysis['pct_with_any_blast']:.1f}%) have BLAST hits")
        else:
            print("  ✗ NOT RECOMMENDED for --blast-only partition")
            print(f"    • Less than 50% of chains ({analysis['pct_with_any_blast']:.1f}%) have BLAST hits")
        
        if analysis['pct_non_peptides_without_hits'] >= 30:
            print("  ✓ RECOMMENDED for HHsearch pipeline")
            print(f"    • {analysis['pct_non_peptides_without_hits']:.1f}% of non-peptide chains have no BLAST hits")
            print(f"    • These {analysis['non_peptides_without_blast_hits']} chains are suitable for HHsearch analysis")
        else:
            print("  ✗ NOT necessary to run HHsearch pipeline")
            print(f"    • Only {analysis['pct_non_peptides_without_hits']:.1f}% of non-peptide chains lack BLAST hits")
        
        print("\n" + "="*80)
    
    def _generate_visualizations(self, batch_id: int, analysis: Dict[str, Any], output_dir: str):
        """Generate visualizations of the analysis results
        
        Args:
            batch_id: Batch ID
            analysis: Analysis results dictionary
            output_dir: Directory to save visualizations
        """
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # 1. Chain length distribution
        if analysis['lengths']:
            plt.figure(figsize=(10, 6))
            plt.hist(analysis['lengths'], bins=50, alpha=0.7, label='All chains')
            plt.axvline(x=self.peptide_threshold, color='r', linestyle='--', 
                      label=f'Peptide threshold ({self.peptide_threshold})')
            
            if analysis['length_with_hits']:
                plt.hist(analysis['length_with_hits'], bins=50, alpha=0.5, color='g', label='Chains with BLAST hits')
            
            plt.xlabel('Chain Length (residues)')
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
            'Non-peptides with BLAST hits',
            'Non-peptides without BLAST hits (HHsearch candidates)',
            'Peptides with BLAST hits',
            'Peptides without BLAST hits'
        ]
        
        sizes = [
            analysis['non_peptides_with_blast_hits'],
            analysis['non_peptides_without_blast_hits'],
            analysis['peptides_with_blast_hits'],
            analysis['peptides'] - analysis['peptides_with_blast_hits']
        ]
        
        # Remove empty categories
        non_zero_indices = [i for i, size in enumerate(sizes) if size > 0]
        filtered_labels = [labels[i] for i in non_zero_indices]
        filtered_sizes = [sizes[i] for i in non_zero_indices]
        
        if filtered_sizes:
            colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3']
            filtered_colors = [colors[i] for i in non_zero_indices]
            
            # Create explode array to highlight HHsearch candidates
            explode = [0.1 if 'HHsearch' in label else 0 for label in filtered_labels]
            
            plt.pie(filtered_sizes, explode=explode, labels=filtered_labels, colors=filtered_colors, 
                  autopct='%1.1f%%', shadow=True, startangle=90)
            plt.axis('equal')
            plt.title(f'Chain Categories - Batch {batch_id}')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'batch_{batch_id}_categories.png'), dpi=300)
            plt.close()
        
        # 3. Blast hit distribution bar chart
        plt.figure(figsize=(10, 6))
        
        categories = ['Any BLAST hits', 'Domain BLAST hits', 'Chain BLAST hits', 'Both types']
        counts = [
            analysis['proteins_with_any_blast_hits'],
            analysis['proteins_with_domain_blast_hits'],
            analysis['proteins_with_chain_blast_hits'],
            analysis['proteins_with_both_blast_hits']
        ]
        
        # Calculate percentages for label display
        percentages = [
            analysis['pct_with_any_blast'],
            analysis['pct_with_domain_blast'],
            analysis['pct_with_chain_blast'],
            (analysis['proteins_with_both_blast_hits'] / analysis['proteins_with_length'] * 100) if analysis['proteins_with_length'] > 0 else 0
        ]
        
        # Plot
        bars = plt.bar(categories, counts, color=['#3274A1', '#E1812C', '#3A923A', '#C03D3E'])
        
        # Add percentage labels
        for i, (bar, percentage) in enumerate(zip(bars, percentages)):
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                   f'{percentage:.1f}%',
                   ha='center', va='bottom', rotation=0)
        
        plt.xlabel('BLAST Hit Type')
        plt.ylabel('Number of Proteins')
        plt.title(f'BLAST Hit Distribution - Batch {batch_id}')
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'batch_{batch_id}_blast_hits.png'), dpi=300)
        plt.close()
        
        # 4. Generate list of HHsearch candidates
        if analysis['non_peptides_without_blast_hits'] > 0:
            with open(os.path.join(output_dir, f'batch_{batch_id}_hhsearch_candidates.txt'), 'w') as f:
                f.write(f"HHsearch Candidates for Batch {batch_id}\n")
                f.write(f"================================\n\n")
                f.write(f"These proteins are non-peptides (length >= {self.peptide_threshold}) with no BLAST hits.\n")
                f.write(f"They are good candidates for more sensitive HHsearch analysis.\n\n")
                
                hhsearch_candidates = []
                for protein in analysis['proteins']:
                    if (protein['length'] >= self.peptide_threshold and 
                        not protein['has_domain_blast_hits'] and 
                        not protein['has_chain_blast_hits']):
                        hhsearch_candidates.append(protein)
                
                # Sort by length (descending)
                hhsearch_candidates.sort(key=lambda p: p['length'], reverse=True)
                
                f.write(f"{'PDB ID':<8} {'Chain':<6} {'Length':<10} {'Stage':<20} {'Status':<10}\n")
                f.write("-" * 60 + "\n")
                
                for protein in hhsearch_candidates:
                    f.write(f"{protein['pdb_id']:<8} {protein['chain_id']:<6} {protein['length']:<10} ")
                    f.write(f"{protein['current_stage']:<20} {protein['status']:<10}\n")
                
                f.write("\n\nTotal: {0} candidates".format(len(hhsearch_candidates)))
        
        # 5. Generate bash script to run HHsearch on candidates
        if analysis['non_peptides_without_blast_hits'] > 0:
            with open(os.path.join(output_dir, f'batch_{batch_id}_run_hhsearch.sh'), 'w') as f:
                f.write("#!/bin/bash\n\n")
                f.write(f"# Generated script to run HHsearch on candidates from batch {batch_id}\n")
                f.write("# These proteins are non-peptides with no BLAST hits\n\n")
                
                f.write("# Set config path - modify if needed\n")
                f.write("CONFIG_PATH=\"config/config.yml\"\n\n")
                
                f.write("# Run HHsearch pipeline\n")
                f.write(f"python ecod/scripts/run_hhsearch.py --config $CONFIG_PATH --batch-id {batch_id} \\\n")
                f.write("  --filter-no-blast-hits --min-length {0} --threads 8\n\n".format(self.peptide_threshold))
                
                f.write("# Alternatively, run for specific proteins:\n")
                
                hhsearch_candidates = []
                for protein in analysis['proteins']:
                    if (protein['length'] >= self.peptide_threshold and 
                        not protein['has_domain_blast_hits'] and 
                        not protein['has_chain_blast_hits']):
                        hhsearch_candidates.append(protein)
                
                # Write first 5 examples
                for i, protein in enumerate(hhsearch_candidates[:5]):
                    f.write(f"# python ecod/scripts/run_hhsearch_single.py --config $CONFIG_PATH ")
                    f.write(f"--protein-id {protein['id']} --batch-id {batch_id}\n")
                
                # Make executable
                os.chmod(os.path.join(output_dir, f'batch_{batch_id}_run_hhsearch.sh'), 0o755)
        
        self.logger.info(f"Generated visualizations in {output_dir}")

def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='Analyze ECOD batch proteins for blast-only partition suitability'
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
    setup_logging(args.verbose, args.log_file)
    
    # If no output directory specified, use current directory
    output_dir = args.output_dir
    if not output_dir:
        output_dir = f"batch_{args.batch_id}_analysis"
    
    # Initialize analyzer and run analysis
    analyzer = BatchBlastAnalyzer(args.config, args.peptide_threshold)
    analyzer.analyze_batch(args.batch_id, output_dir)

if __name__ == "__main__":
    main()