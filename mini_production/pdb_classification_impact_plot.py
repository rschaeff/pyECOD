#!/usr/bin/env python3
"""
PDB Classification Impact Analysis

Creates a yearly breakdown showing:
- Existing PDB classifications (main system)
- New Mini PyECOD classifications  
- Remaining unclassified PDBs

This visualizes the impact of mini PyECOD on the classification crisis.
"""

import yaml
import psycopg2
import psycopg2.extras
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from datetime import datetime
import numpy as np
import argparse
from pathlib import Path

class PDBClassificationImpactAnalyzer:
    """Analyze and visualize PDB classification impact by year"""
    
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
    
    def get_pdb_classification_data(self):
        """Get comprehensive PDB classification data by year"""
        
        print("🔍 Analyzing PDB classification status by year...")
        
        with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            # Get comprehensive PDB classification data
            cursor.execute("""
                WITH pdb_years AS (
                    -- Get PDB entries with their deposition years
                    SELECT 
                        pe.pdb_id,
                        EXTRACT(YEAR FROM pe.deposition_date) as deposition_year,
                        pe.deposition_date
                    FROM pdb_analysis.pdb_entries pe
                    WHERE pe.deposition_date IS NOT NULL
                      AND EXTRACT(YEAR FROM pe.deposition_date) BETWEEN 2015 AND 2025
                ),
                chain_classifications AS (
                    -- Classify each protein chain by classification status
                    SELECT 
                        p.pdb_id,
                        p.chain_id,
                        p.source_id,
                        py.deposition_year,
                        py.deposition_date,
                        -- Check for existing (main system) classifications
                        CASE WHEN pp_main.id IS NOT NULL AND pp_main.is_classified = true 
                             THEN 'existing_classified'
                             ELSE NULL END as main_status,
                        -- Check for mini classifications  
                        CASE WHEN pp_mini.id IS NOT NULL 
                             THEN 'mini_classified'
                             ELSE NULL END as mini_status,
                        -- Check for mini propagated classifications
                        CASE WHEN pp_prop.id IS NOT NULL 
                             THEN 'mini_propagated'
                             ELSE NULL END as prop_status
                    FROM pdb_analysis.protein p
                    JOIN pdb_years py ON p.pdb_id = py.pdb_id
                    -- Left join with main system classifications
                    LEFT JOIN pdb_analysis.partition_proteins pp_main ON 
                        p.pdb_id = pp_main.pdb_id AND p.chain_id = pp_main.chain_id
                        AND (pp_main.process_version = '1.0' OR pp_main.process_version IS NULL)
                        AND pp_main.is_classified = true
                    -- Left join with mini direct classifications
                    LEFT JOIN pdb_analysis.partition_proteins pp_mini ON 
                        p.pdb_id = pp_mini.pdb_id AND p.chain_id = pp_mini.chain_id
                        AND pp_mini.process_version = 'mini_pyecod_1.0'
                    -- Left join with mini propagated classifications
                    LEFT JOIN pdb_analysis.partition_proteins pp_prop ON 
                        p.pdb_id = pp_prop.pdb_id AND p.chain_id = pp_prop.chain_id
                        AND pp_prop.process_version = 'mini_pyecod_propagated_1.0'
                ),
                final_classification AS (
                    -- Determine final classification status (priority: mini > main > unclassified)
                    SELECT 
                        pdb_id,
                        chain_id,
                        source_id,
                        deposition_year,
                        deposition_date,
                        CASE 
                            -- Mini classifications (direct or propagated) take priority
                            WHEN mini_status = 'mini_classified' OR prop_status = 'mini_propagated' 
                                THEN 'mini_classified'
                            -- Then existing main system classifications
                            WHEN main_status = 'existing_classified' 
                                THEN 'existing_classified'
                            -- Otherwise unclassified
                            ELSE 'unclassified'
                        END as classification_status
                    FROM chain_classifications
                )
                SELECT 
                    deposition_year,
                    classification_status,
                    COUNT(*) as chain_count,
                    COUNT(DISTINCT pdb_id) as structure_count
                FROM final_classification
                WHERE deposition_year IS NOT NULL
                GROUP BY deposition_year, classification_status
                ORDER BY deposition_year, classification_status
            """)
            
            results = cursor.fetchall()
            
        # Convert to DataFrame for easier manipulation
        df = pd.DataFrame([dict(row) for row in results])
        print(f"✓ Retrieved classification data for {len(df)} year/status combinations")
        
        return df
    
    def create_impact_visualization(self, df, output_path: str = None, 
                                  metric: str = 'chain_count'):
        """Create the impact visualization showing before/after mini PyECOD"""
        
        print(f"📊 Creating classification impact visualization (metric: {metric})...")
        
        # Pivot data for stacked bar chart
        pivot_df = df.pivot(index='deposition_year', 
                           columns='classification_status', 
                           values=metric).fillna(0)
        
        # Ensure all columns exist
        for col in ['existing_classified', 'mini_classified', 'unclassified']:
            if col not in pivot_df.columns:
                pivot_df[col] = 0
        
        # Reorder columns for proper stacking
        column_order = ['existing_classified', 'mini_classified', 'unclassified']
        pivot_df = pivot_df[column_order]
        
        # Create the plot
        fig, ax = plt.subplots(figsize=(14, 8))
        
        # Color scheme
        colors = {
            'existing_classified': '#2E8B57',     # Sea Green (darker green)
            'mini_classified': '#90EE90',         # Light Green (brighter)
            'unclassified': '#DC143C'             # Crimson Red
        }
        
        # Create stacked bar chart
        bottom = np.zeros(len(pivot_df))
        bars = []
        
        for column in column_order:
            if column in pivot_df.columns:
                bars.append(ax.bar(pivot_df.index, pivot_df[column], 
                                 bottom=bottom, 
                                 color=colors[column], 
                                 label=column.replace('_', ' ').title(),
                                 edgecolor='white', linewidth=0.5))
                bottom += pivot_df[column]
        
        # Customize the plot
        ax.set_xlabel('Year', fontsize=12, fontweight='bold')
        ylabel = 'Number of Protein Chains' if metric == 'chain_count' else 'Number of Structures'
        ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
        
        # Title with impact emphasis
        title = f'PDB Classification Impact: Before and After Mini PyECOD\n'
        title += f'Showing the transformation of unclassified proteins (metric: {ylabel.lower()})'
        ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
        
        # Format axes
        ax.set_xlim(pivot_df.index.min() - 0.5, pivot_df.index.max() + 0.5)
        ax.tick_params(axis='both', which='major', labelsize=10)
        
        # Add commas to y-axis labels
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x):,}'))
        
        # Custom legend
        legend_elements = [
            mpatches.Patch(color=colors['existing_classified'], 
                         label='Existing Classifications (Main System)'),
            mpatches.Patch(color=colors['mini_classified'], 
                         label='New Mini PyECOD Classifications'),
            mpatches.Patch(color=colors['unclassified'], 
                         label='Remaining Unclassified')
        ]
        
        ax.legend(handles=legend_elements, loc='upper left', fontsize=11)
        
        # Add grid for better readability
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_axisbelow(True)
        
        # Add annotation about the mini impact
        if 'mini_classified' in pivot_df.columns and pivot_df['mini_classified'].sum() > 0:
            total_mini = int(pivot_df['mini_classified'].sum())
            total_rescued = int(pivot_df['unclassified'].sum() + pivot_df['mini_classified'].sum())
            
            # Find the best year to place annotation
            annotation_year = pivot_df.index[-2] if len(pivot_df.index) > 1 else pivot_df.index[-1]
            annotation_height = pivot_df.loc[annotation_year].sum() * 0.8
            
            ax.annotate(f'Mini PyECOD Impact:\n{total_mini:,} new classifications',
                       xy=(annotation_year, annotation_height), 
                       xytext=(annotation_year - 1, annotation_height * 1.2),
                       fontsize=10, fontweight='bold',
                       bbox=dict(boxstyle="round,pad=0.3", facecolor='lightblue', alpha=0.8),
                       arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.1'))
        
        plt.tight_layout()
        
        # Save if output path provided
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"✓ Saved visualization to: {output_path}")
        
        plt.show()
        
        return fig, ax
    
    def generate_impact_statistics(self, df):
        """Generate detailed impact statistics"""
        
        print("📈 Generating impact statistics...")
        
        # Calculate totals by category
        totals = df.groupby('classification_status')[['chain_count', 'structure_count']].sum()
        
        # Calculate yearly trends
        yearly_totals = df.groupby('deposition_year')[['chain_count', 'structure_count']].sum()
        yearly_breakdown = df.pivot(index='deposition_year', 
                                  columns='classification_status', 
                                  values='chain_count').fillna(0)
        
        # Mini impact analysis
        mini_impact = {}
        if 'mini_classified' in totals.index:
            mini_chains = totals.loc['mini_classified', 'chain_count']
            mini_structures = totals.loc['mini_classified', 'structure_count']
            total_chains = totals['chain_count'].sum()
            total_structures = totals['structure_count'].sum()
            
            mini_impact = {
                'mini_chains': int(mini_chains),
                'mini_structures': int(mini_structures),
                'mini_chain_percentage': mini_chains / total_chains * 100,
                'mini_structure_percentage': mini_structures / total_structures * 100,
                'total_chains': int(total_chains),
                'total_structures': int(total_structures)
            }
        
        # Calculate classification rates by year
        classification_rates = yearly_breakdown.div(yearly_breakdown.sum(axis=1), axis=0) * 100
        
        return {
            'totals': totals,
            'yearly_totals': yearly_totals,
            'yearly_breakdown': yearly_breakdown,
            'classification_rates': classification_rates,
            'mini_impact': mini_impact
        }
    
    def print_impact_report(self, stats):
        """Print comprehensive impact report"""
        
        print("\n🎯 PDB Classification Impact Report")
        print("=" * 60)
        print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print()
        
        # Overall totals
        print("📊 Overall Classification Totals:")
        totals = stats['totals']
        for status in totals.index:
            chains = int(totals.loc[status, 'chain_count'])
            structures = int(totals.loc[status, 'structure_count'])
            print(f"  {status.replace('_', ' ').title():<25} {chains:>8,} chains, {structures:>8,} structures")
        
        total_chains = int(totals['chain_count'].sum())
        total_structures = int(totals['structure_count'].sum())
        print(f"  {'Total':<25} {total_chains:>8,} chains, {total_structures:>8,} structures")
        print()
        
        # Mini impact
        if stats['mini_impact']:
            impact = stats['mini_impact']
            print("🚀 Mini PyECOD Impact:")
            print(f"  New classifications:     {impact['mini_chains']:>8,} chains ({impact['mini_chain_percentage']:.1f}%)")
            print(f"  New structures:          {impact['mini_structures']:>8,} structures ({impact['mini_structure_percentage']:.1f}%)")
            print()
        
        # Yearly trends (recent years)
        print("📈 Recent Yearly Trends (2020-2025):")
        yearly = stats['yearly_breakdown']
        rates = stats['classification_rates']
        
        print(f"{'Year':<6} {'Total':<8} {'Existing':<10} {'Mini':<8} {'Unclass':<10} {'Class Rate':<12}")
        print("-" * 65)
        
        for year in sorted(yearly.index):
            if year >= 2020:
                total = int(yearly.loc[year].sum())
                existing = int(yearly.loc[year].get('existing_classified', 0))
                mini = int(yearly.loc[year].get('mini_classified', 0))
                unclass = int(yearly.loc[year].get('unclassified', 0))
                
                class_rate = (existing + mini) / total * 100 if total > 0 else 0
                
                print(f"{int(year):<6} {total:<8,} {existing:<10,} {mini:<8,} {unclass:<10,} {class_rate:<11.1f}%")
        
        print()
        
        # Key insights
        print("💡 Key Insights:")
        
        if stats['mini_impact']:
            impact = stats['mini_impact']
            print(f"  • Mini PyECOD classified {impact['mini_chains']:,} protein chains")
            print(f"  • This represents {impact['mini_chain_percentage']:.1f}% of all analyzed chains")
            
        # Calculate rescue rate
        if 'unclassified' in totals.index and 'mini_classified' in totals.index:
            unclassified = totals.loc['unclassified', 'chain_count']
            mini_classified = totals.loc['mini_classified', 'chain_count']
            total_needing_classification = unclassified + mini_classified
            rescue_rate = mini_classified / total_needing_classification * 100
            print(f"  • Rescue rate: {rescue_rate:.1f}% of previously unclassified proteins")
        
        print(f"  • Database contains {total_chains:,} total protein chains")
        print(f"  • Analysis covers {len(yearly.index)} years of PDB depositions")


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Analyze and visualize PDB classification impact'
    )
    
    parser.add_argument('--config', type=str, default='config/config.local.yml',
                       help='Config file path')
    parser.add_argument('--output', type=str, 
                       help='Output file path for plot (e.g., pdb_impact.png)')
    parser.add_argument('--metric', type=str, choices=['chain_count', 'structure_count'],
                       default='chain_count',
                       help='Metric to visualize (default: chain_count)')
    parser.add_argument('--stats-only', action='store_true',
                       help='Generate statistics only (no plot)')
    
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = PDBClassificationImpactAnalyzer(args.config)
    
    # Get data
    df = analyzer.get_pdb_classification_data()
    
    if df.empty:
        print("❌ No data found. Check database connection and data availability.")
        return
    
    # Generate statistics
    stats = analyzer.generate_impact_statistics(df)
    analyzer.print_impact_report(stats)
    
    # Create visualization unless stats-only
    if not args.stats_only:
        analyzer.create_impact_visualization(df, args.output, args.metric)
    
    print("✅ Analysis complete!")


if __name__ == "__main__":
    main()
