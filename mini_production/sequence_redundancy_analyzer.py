#!/usr/bin/env python3
"""
Sequence Redundancy Analysis for Unlimited Propagation

Analyzes the true potential of sequence-based propagation by examining
how many identical sequences exist for each mini classification.

This reveals the massive untapped propagation potential when the 
--limit-per-sequence constraint is removed.
"""

import yaml
import psycopg2
import psycopg2.extras
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import defaultdict
import argparse
from datetime import datetime

class SequenceRedundancyAnalyzer:
    """Analyze sequence redundancy and propagation potential"""
    
    def __init__(self, config_path: str = "config/config.local.yml"):
        self.config = self._load_config(config_path)
        self.db_conn = self._init_db_connection()
        
    def _load_config(self, config_path: str):
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    
    def _init_db_connection(self):
        try:
            return psycopg2.connect(**self.config['database'])
        except Exception as e:
            print(f"Database connection failed: {e}")
            raise
    
    def analyze_sequence_redundancy(self):
        """Analyze sequence redundancy for mini classifications"""
        
        print("🔍 Analyzing sequence redundancy for mini classifications...")
        
        with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            # Get all mini classifications with their sequence MD5s
            cursor.execute("""
                SELECT DISTINCT
                    pp.pdb_id,
                    pp.chain_id,
                    ps.sequence_md5,
                    pp.domains_with_evidence,
                    pp.process_version
                FROM pdb_analysis.partition_proteins pp
                JOIN pdb_analysis.protein p ON pp.pdb_id = p.pdb_id AND pp.chain_id = p.chain_id
                JOIN pdb_analysis.protein_sequence ps ON p.id = ps.protein_id
                WHERE pp.process_version = 'mini_pyecod_1.0'
                  AND ps.sequence_md5 IS NOT NULL
                ORDER BY pp.pdb_id, pp.chain_id
            """)
            
            mini_sequences = cursor.fetchall()
            print(f"✓ Found {len(mini_sequences)} mini classifications with sequence data")
            
            # For each mini sequence, find ALL other proteins with same MD5
            redundancy_data = []
            sequence_md5s = list(set(row['sequence_md5'] for row in mini_sequences))
            
            print(f"🔍 Analyzing redundancy for {len(sequence_md5s)} unique sequences...")
            
            for i, md5 in enumerate(sequence_md5s):
                if i % 100 == 0:
                    print(f"   Progress: {i}/{len(sequence_md5s)} sequences analyzed")
                
                # Find all proteins with this sequence MD5
                cursor.execute("""
                    SELECT 
                        p.id as protein_id,
                        p.pdb_id,
                        p.chain_id,
                        p.source_id,
                        ps.sequence_md5,
                        p.length,
                        -- Check if already classified by any system
                        CASE 
                            WHEN pp_any.id IS NOT NULL THEN true 
                            ELSE false 
                        END as already_classified,
                        -- Check classification type
                        CASE 
                            WHEN pp_mini.id IS NOT NULL THEN 'mini_direct'
                            WHEN pp_prop.id IS NOT NULL THEN 'mini_propagated'  
                            WHEN pp_main.id IS NOT NULL AND pp_main.is_classified THEN 'main_system'
                            ELSE 'unclassified'
                        END as classification_type
                    FROM pdb_analysis.protein p
                    JOIN pdb_analysis.protein_sequence ps ON p.id = ps.protein_id
                    -- Check for any existing classification
                    LEFT JOIN pdb_analysis.partition_proteins pp_any ON 
                        p.pdb_id = pp_any.pdb_id AND p.chain_id = pp_any.chain_id
                    -- Check for mini direct
                    LEFT JOIN pdb_analysis.partition_proteins pp_mini ON 
                        p.pdb_id = pp_mini.pdb_id AND p.chain_id = pp_mini.chain_id
                        AND pp_mini.process_version = 'mini_pyecod_1.0'
                    -- Check for mini propagated
                    LEFT JOIN pdb_analysis.partition_proteins pp_prop ON 
                        p.pdb_id = pp_prop.pdb_id AND p.chain_id = pp_prop.chain_id
                        AND pp_prop.process_version = 'mini_pyecod_propagated_1.0'
                    -- Check for main system
                    LEFT JOIN pdb_analysis.partition_proteins pp_main ON 
                        p.pdb_id = pp_main.pdb_id AND p.chain_id = pp_main.chain_id
                        AND (pp_main.process_version = '1.0' OR pp_main.process_version IS NULL)
                        AND pp_main.is_classified = true
                    WHERE ps.sequence_md5 = %s
                    ORDER BY p.pdb_id, p.chain_id
                """, (md5,))
                
                all_chains = cursor.fetchall()
                
                # Count classifications by type
                classification_counts = defaultdict(int)
                unclassified_candidates = []
                
                for chain in all_chains:
                    classification_counts[chain['classification_type']] += 1
                    if chain['classification_type'] == 'unclassified':
                        unclassified_candidates.append(chain)
                
                # Record redundancy data
                redundancy_data.append({
                    'sequence_md5': md5,
                    'total_chains': len(all_chains),
                    'mini_direct': classification_counts['mini_direct'],
                    'mini_propagated': classification_counts['mini_propagated'],
                    'main_system': classification_counts['main_system'],
                    'unclassified': classification_counts['unclassified'],
                    'propagation_potential': len(unclassified_candidates)
                })
        
        print(f"✓ Redundancy analysis complete for {len(redundancy_data)} sequences")
        return pd.DataFrame(redundancy_data)
    
    def calculate_propagation_potential(self, redundancy_df):
        """Calculate the true propagation potential"""
        
        print("📊 Calculating propagation potential...")
        
        # Current propagation (limited to 10 per sequence)
        current_propagated = redundancy_df['mini_propagated'].sum()
        
        # Unlimited propagation potential
        unlimited_potential = redundancy_df['propagation_potential'].sum()
        
        # Additional propagation possible
        additional_possible = unlimited_potential - current_propagated
        
        # High-impact sequences (those with many unclassified chains)
        high_impact = redundancy_df[redundancy_df['propagation_potential'] >= 50]
        
        # Sequence distribution analysis
        distribution = {
            'sequences_with_1_chain': len(redundancy_df[redundancy_df['total_chains'] == 1]),
            'sequences_with_2_10_chains': len(redundancy_df[(redundancy_df['total_chains'] >= 2) & (redundancy_df['total_chains'] <= 10)]),
            'sequences_with_11_50_chains': len(redundancy_df[(redundancy_df['total_chains'] >= 11) & (redundancy_df['total_chains'] <= 50)]),
            'sequences_with_51_100_chains': len(redundancy_df[(redundancy_df['total_chains'] >= 51) & (redundancy_df['total_chains'] <= 100)]),
            'sequences_with_100plus_chains': len(redundancy_df[redundancy_df['total_chains'] > 100])
        }
        
        potential_summary = {
            'current_propagated': current_propagated,
            'unlimited_potential': unlimited_potential,
            'additional_possible': additional_possible,
            'amplification_factor': unlimited_potential / max(1, redundancy_df['mini_direct'].sum()),
            'high_impact_sequences': len(high_impact),
            'high_impact_potential': high_impact['propagation_potential'].sum(),
            'distribution': distribution
        }
        
        return potential_summary, high_impact
    
    def create_redundancy_visualization(self, redundancy_df, output_path=None):
        """Create visualizations of sequence redundancy"""
        
        print("📊 Creating redundancy visualizations...")
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
        
        # 1. Distribution of total chains per sequence
        chain_counts = redundancy_df['total_chains']
        ax1.hist(chain_counts, bins=50, edgecolor='black', alpha=0.7, color='skyblue')
        ax1.set_xlabel('Number of Chains per Sequence')
        ax1.set_ylabel('Number of Sequences')
        ax1.set_title('Distribution of Sequence Redundancy\n(How many chains share each sequence)')
        ax1.set_yscale('log')
        
        # Add statistics
        mean_chains = chain_counts.mean()
        median_chains = chain_counts.median()
        max_chains = chain_counts.max()
        ax1.axvline(mean_chains, color='red', linestyle='--', label=f'Mean: {mean_chains:.1f}')
        ax1.axvline(median_chains, color='orange', linestyle='--', label=f'Median: {median_chains:.1f}')
        ax1.legend()
        ax1.text(0.7, 0.8, f'Max: {max_chains}', transform=ax1.transAxes, fontsize=10)
        
        # 2. Propagation potential distribution
        prop_potential = redundancy_df['propagation_potential']
        ax2.hist(prop_potential, bins=50, edgecolor='black', alpha=0.7, color='lightgreen')
        ax2.set_xlabel('Propagation Potential per Sequence')
        ax2.set_ylabel('Number of Sequences')
        ax2.set_title('Distribution of Propagation Potential\n(Unclassified chains per sequence)')
        ax2.set_yscale('log')
        
        # 3. High-impact sequences
        high_impact = redundancy_df[redundancy_df['propagation_potential'] >= 20]
        if len(high_impact) > 0:
            ax3.scatter(high_impact['total_chains'], high_impact['propagation_potential'], 
                       alpha=0.6, s=50, color='red')
            ax3.set_xlabel('Total Chains for Sequence')
            ax3.set_ylabel('Propagation Potential')
            ax3.set_title('High-Impact Sequences\n(≥20 propagation potential)')
            
            # Add trend line
            if len(high_impact) > 1:
                z = np.polyfit(high_impact['total_chains'], high_impact['propagation_potential'], 1)
                p = np.poly1d(z)
                ax3.plot(high_impact['total_chains'], p(high_impact['total_chains']), "r--", alpha=0.8)
        
        # 4. Classification breakdown
        classification_data = {
            'Mini Direct': redundancy_df['mini_direct'].sum(),
            'Mini Propagated': redundancy_df['mini_propagated'].sum(),
            'Main System': redundancy_df['main_system'].sum(),
            'Unclassified': redundancy_df['unclassified'].sum()
        }
        
        colors = ['#90EE90', '#32CD32', '#2E8B57', '#DC143C']
        wedges, texts, autotexts = ax4.pie(classification_data.values(), 
                                          labels=classification_data.keys(),
                                          autopct='%1.1f%%',
                                          colors=colors,
                                          startangle=90)
        ax4.set_title('Classification Status\n(All chains with mini sequence matches)')
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"✓ Saved visualization to: {output_path}")
        
        plt.show()
        return fig
    
    def print_redundancy_report(self, redundancy_df, potential_summary, high_impact):
        """Print comprehensive redundancy analysis report"""
        
        print("\n🔍 Sequence Redundancy Analysis Report")
        print("=" * 60)
        print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print()
        
        # Overall statistics
        total_sequences = len(redundancy_df)
        total_chains = redundancy_df['total_chains'].sum()
        avg_redundancy = redundancy_df['total_chains'].mean()
        max_redundancy = redundancy_df['total_chains'].max()
        
        print("📊 Overall Redundancy Statistics:")
        print(f"  Unique mini sequences:       {total_sequences:>8,}")
        print(f"  Total matching chains:       {total_chains:>8,}")
        print(f"  Average redundancy:          {avg_redundancy:>8.1f} chains/sequence")
        print(f"  Maximum redundancy:          {max_redundancy:>8,} chains/sequence")
        print()
        
        # Current vs potential propagation
        print("🚀 Propagation Analysis:")
        print(f"  Current propagated:          {potential_summary['current_propagated']:>8,}")
        print(f"  Unlimited potential:         {potential_summary['unlimited_potential']:>8,}")
        print(f"  Additional possible:         {potential_summary['additional_possible']:>8,}")
        print(f"  Total amplification:         {potential_summary['amplification_factor']:>8.1f}x")
        print()
        
        # Distribution breakdown
        print("📈 Redundancy Distribution:")
        dist = potential_summary['distribution']
        for category, count in dist.items():
            pct = count / total_sequences * 100
            print(f"  {category.replace('_', ' ').title():<25} {count:>6,} ({pct:>5.1f}%)")
        print()
        
        # High-impact sequences
        print("🎯 High-Impact Sequences (≥50 propagation potential):")
        print(f"  Count:                       {potential_summary['high_impact_sequences']:>8,}")
        print(f"  Total potential:             {potential_summary['high_impact_potential']:>8,}")
        
        if len(high_impact) > 0:
            print(f"\n  Top 10 highest impact sequences:")
            print(f"  {'Sequence MD5':<10} {'Total':<8} {'Unclass':<10} {'Potential':<12}")
            print("  " + "-" * 45)
            
            top_impact = high_impact.nlargest(10, 'propagation_potential')
            for _, row in top_impact.iterrows():
                md5_short = row['sequence_md5'][:8] + "..."
                print(f"  {md5_short:<10} {row['total_chains']:<8} {row['unclassified']:<10} {row['propagation_potential']:<12}")
        
        print()
        
        # Impact projection
        if potential_summary['additional_possible'] > 0:
            print("💡 Impact Projection:")
            current_classified = redundancy_df['mini_direct'].sum() + redundancy_df['mini_propagated'].sum()
            total_possible = current_classified + potential_summary['additional_possible']
            improvement = potential_summary['additional_possible'] / current_classified * 100
            
            print(f"  Current mini classified:     {current_classified:>8,}")
            print(f"  Potential total classified:  {total_possible:>8,}")
            print(f"  Improvement factor:          {improvement:>8.1f}% increase")
            
            # Database impact
            total_database_chains = 669925  # From previous analysis
            new_coverage = total_possible / total_database_chains * 100
            current_coverage = current_classified / total_database_chains * 100
            
            print(f"  Current database coverage:   {current_coverage:>8.1f}%")
            print(f"  Potential database coverage: {new_coverage:>8.1f}%")
    
    def find_unlimited_propagation_candidates(self, limit_threshold=50):
        """Find sequences with high propagation potential for unlimited propagation"""
        
        print(f"🎯 Finding sequences with ≥{limit_threshold} propagation potential...")
        
        with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            # Find mini sequences with many unclassified matches
            cursor.execute("""
                WITH mini_sequences AS (
                    SELECT DISTINCT ps.sequence_md5
                    FROM pdb_analysis.partition_proteins pp
                    JOIN pdb_analysis.protein p ON pp.pdb_id = p.pdb_id AND pp.chain_id = p.chain_id
                    JOIN pdb_analysis.protein_sequence ps ON p.id = ps.protein_id
                    WHERE pp.process_version = 'mini_pyecod_1.0'
                      AND ps.sequence_md5 IS NOT NULL
                ),
                sequence_analysis AS (
                    SELECT 
                        ms.sequence_md5,
                        COUNT(*) as total_chains,
                        COUNT(*) FILTER (WHERE pp.id IS NULL) as unclassified_chains
                    FROM mini_sequences ms
                    JOIN pdb_analysis.protein_sequence ps ON ms.sequence_md5 = ps.sequence_md5
                    JOIN pdb_analysis.protein p ON ps.protein_id = p.id
                    LEFT JOIN pdb_analysis.partition_proteins pp ON 
                        p.pdb_id = pp.pdb_id AND p.chain_id = pp.chain_id
                    GROUP BY ms.sequence_md5
                )
                SELECT *
                FROM sequence_analysis
                WHERE unclassified_chains >= %s
                ORDER BY unclassified_chains DESC
            """, (limit_threshold,))
            
            candidates = cursor.fetchall()
            
        print(f"✓ Found {len(candidates)} sequences with ≥{limit_threshold} propagation potential")
        return [dict(row) for row in candidates]


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Analyze sequence redundancy and propagation potential'
    )
    
    parser.add_argument('--config', type=str, default='config/config.local.yml',
                       help='Config file path')
    parser.add_argument('--output', type=str, 
                       help='Output file path for plots')
    parser.add_argument('--find-candidates', type=int, metavar='THRESHOLD',
                       help='Find sequences with ≥THRESHOLD propagation potential')
    parser.add_argument('--no-plot', action='store_true',
                       help='Skip visualization generation')
    
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = SequenceRedundancyAnalyzer(args.config)
    
    # Analyze redundancy
    redundancy_df = analyzer.analyze_sequence_redundancy()
    potential_summary, high_impact = analyzer.calculate_propagation_potential(redundancy_df)
    
    # Print report
    analyzer.print_redundancy_report(redundancy_df, potential_summary, high_impact)
    
    # Create visualization unless skipped
    if not args.no_plot:
        analyzer.create_redundancy_visualization(redundancy_df, args.output)
    
    # Find high-potential candidates if requested
    if args.find_candidates:
        candidates = analyzer.find_unlimited_propagation_candidates(args.find_candidates)
        print(f"\n🎯 High-Potential Propagation Candidates (≥{args.find_candidates}):")
        print(f"{'Sequence MD5':<35} {'Total Chains':<12} {'Unclassified':<12}")
        print("-" * 65)
        for candidate in candidates[:20]:  # Show top 20
            print(f"{candidate['sequence_md5']:<35} {candidate['total_chains']:<12} {candidate['unclassified_chains']:<12}")
    
    print("✅ Redundancy analysis complete!")


if __name__ == "__main__":
    main()
