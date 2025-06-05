#!/usr/bin/env python3
"""
Example Protein Visualizations for Slide Deck

Generate detailed visualizations for specific protein examples to showcase
mini_production results in presentations.

Usage:
    python example_visualizations.py --protein 8ovp_A --comparison
    python example_visualizations.py --showcase-category excellent --count 3
    python example_visualizations.py --improvement-examples --top 5
"""

import os
import sys
import argparse
import yaml
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import xml.etree.ElementTree as ET
import logging

# Import our modules
sys.path.insert(0, str(Path(__file__).parent))
from assess_quality import QualityAssessment, ProteinQuality
from compare_results import ResultComparator, ComparisonResult

# Try to import PyMOL visualization
try:
    sys.path.insert(0, str(Path(__file__).parent / "../mini/core"))
    from visualization import PyMOLVisualizer
    PYMOL_AVAILABLE = True
except ImportError:
    PYMOL_AVAILABLE = False
    logging.warning("PyMOL visualization not available")

logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ExampleVisualizer:
    """Generate detailed visualizations for example proteins"""
    
    def __init__(self, config_path: str = "config.local.yml", output_dir: str = "/tmp/example_viz"):
        self.config_path = config_path
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Load config
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)
        
        # Initialize tools
        self.quality_assessor = QualityAssessment(config_path)
        self.comparator = ResultComparator(config_path)
        
        if PYMOL_AVAILABLE:
            self.pymol_viz = PyMOLVisualizer()
        
        # Set up plotting
        plt.style.use('seaborn-v0_8')
    
    def visualize_protein_example(self, protein_id: str, create_pymol: bool = True) -> Dict:
        """Create comprehensive visualization for a single protein example"""
        
        logger.info(f"Creating visualization for {protein_id}")
        
        # Find the mini result file
        mini_file = self._find_mini_result_file(protein_id)
        if not mini_file:
            logger.error(f"No mini result found for {protein_id}")
            return {}
        
        # Parse mini result
        quality_result = self.quality_assessor.parse_mini_result(mini_file)
        
        # Try to get comparison with main
        comparison_result = None
        pdb_id = protein_id.split('_')[0]
        chain_id = protein_id.split('_')[1] if '_' in protein_id else 'A'
        
        try:
            comparison_result = self.comparator.compare_protein(protein_id, mini_file)
        except Exception as e:
            logger.warning(f"Could not compare with main result: {e}")
        
        # Create domain visualization
        self._create_domain_visualization(protein_id, quality_result, comparison_result)
        
        # Create PyMOL script if available
        pymol_script = None
        if create_pymol and PYMOL_AVAILABLE:
            try:
                pymol_script = self._create_pymol_example(protein_id, mini_file.parent.parent)
            except Exception as e:
                logger.warning(f"Could not create PyMOL visualization: {e}")
        
        # Create summary card
        summary_data = self._create_protein_summary_card(protein_id, quality_result, comparison_result)
        
        return {
            "protein_id": protein_id,
            "quality_result": quality_result,
            "comparison_result": comparison_result,
            "pymol_script": pymol_script,
            "summary_data": summary_data,
            "visualization_files": [
                f"{protein_id}_domain_visualization.png",
                f"{protein_id}_summary_card.png"
            ]
        }
    
    def _find_mini_result_file(self, protein_id: str) -> Optional[Path]:
        """Find mini result file for a protein"""
        batch_base = Path(self.config["paths"]["batch_base_dir"])
        
        # Search across all batches
        for batch_dir in batch_base.iterdir():
            if not batch_dir.is_dir():
                continue
            
            mini_domains_dir = batch_dir / "mini_domains"
            if not mini_domains_dir.exists():
                continue
            
            result_file = mini_domains_dir / f"{protein_id}.mini.domains.xml"
            if result_file.exists():
                return result_file
        
        return None
    
    def _create_domain_visualization(self, protein_id: str, quality_result: ProteinQuality, 
                                   comparison_result: Optional[ComparisonResult]):
        """Create domain architecture visualization"""
        
        fig, ax = plt.subplots(figsize=(14, 8))
        
        # Draw sequence as a horizontal bar
        sequence_length = quality_result.sequence_length or 500  # Fallback
        
        # Main sequence bar
        sequence_bar = patches.Rectangle((0, 2), sequence_length, 1, 
                                       linewidth=2, edgecolor='black', 
                                       facecolor='lightgray', alpha=0.5)
        ax.add_patch(sequence_bar)
        
        # Domain colors
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FECA57', '#9B59B6', '#E67E22']
        
        # Draw mini domains
        y_mini = 2.5
        for i, domain in enumerate(quality_result.domains):
            color = colors[i % len(colors)]
            
            # Parse domain range - handle discontinuous domains
            try:
                from ecod.core.sequence_range import SequenceRange
                seq_range = SequenceRange.parse(domain.range)
                
                for start, end in seq_range.segments:
                    width = end - start + 1
                    domain_rect = patches.Rectangle((start-1, y_mini), width, 0.3,
                                                  linewidth=1, edgecolor='black',
                                                  facecolor=color, alpha=0.8)
                    ax.add_patch(domain_rect)
                    
                    # Add domain label
                    mid_point = start + width/2 - 1
                    ax.text(mid_point, y_mini + 0.15, f"D{i+1}", 
                           ha='center', va='center', fontsize=8, fontweight='bold')
            
            except Exception as e:
                logger.warning(f"Could not parse domain range '{domain.range}': {e}")
        
        # Draw main domains if comparison available
        if comparison_result:
            y_main = 1.5
            
            # Simplified main domains representation (would need actual main result parsing)
            main_domain_count = comparison_result.main_domains
            if main_domain_count > 0:
                domain_width = sequence_length / main_domain_count
                for i in range(main_domain_count):
                    start = i * domain_width
                    main_rect = patches.Rectangle((start, y_main), domain_width, 0.3,
                                                linewidth=1, edgecolor='gray',
                                                facecolor='lightblue', alpha=0.6)
                    ax.add_patch(main_rect)
                    
                    ax.text(start + domain_width/2, y_main + 0.15, f"M{i+1}",
                           ha='center', va='center', fontsize=8, fontweight='bold')
        
        # Formatting
        ax.set_xlim(0, sequence_length)
        ax.set_ylim(1, 4)
        ax.set_xlabel('Residue Position', fontsize=12)
        ax.set_title(f'Domain Architecture: {protein_id}', fontsize=16, fontweight='bold')
        
        # Legend
        legend_elements = [
            patches.Patch(color='lightgray', alpha=0.5, label='Protein Sequence'),
            patches.Patch(color=colors[0], alpha=0.8, label='Mini Domains'),
        ]
        
        if comparison_result:
            legend_elements.append(patches.Patch(color='lightblue', alpha=0.6, label='Main Domains'))
        
        ax.legend(handles=legend_elements, loc='upper right')
        
        # Add quality information
        info_text = f"Length: {sequence_length} residues\n"
        info_text += f"Mini: {quality_result.total_domains} domains ({quality_result.classified_domains} classified)\n"
        info_text += f"Coverage: {quality_result.coverage_percentage:.1f}%\n"
        info_text += f"Tier: {quality_result.tier}"
        
        if comparison_result:
            info_text += f"\nMain: {comparison_result.main_domains} domains\n"
            info_text += f"Winner: {comparison_result.winner}"
        
        ax.text(0.02, 0.98, info_text, transform=ax.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), fontsize=10)
        
        # Remove y-axis ticks and labels
        ax.set_yticks([])
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / f"{protein_id}_domain_visualization.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Domain visualization saved for {protein_id}")
    
    def _create_protein_summary_card(self, protein_id: str, quality_result: ProteinQuality,
                                   comparison_result: Optional[ComparisonResult]) -> Dict:
        """Create a summary card visualization"""
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.axis('off')
        
        # Title
        ax.text(0.5, 0.95, f"Protein Analysis: {protein_id}", 
               transform=ax.transAxes, ha='center', fontsize=20, fontweight='bold')
        
        # Quality metrics section
        y_pos = 0.8
        ax.text(0.05, y_pos, "MINI ALGORITHM RESULTS", transform=ax.transAxes,
               fontsize=14, fontweight='bold', color='darkblue')
        
        y_pos -= 0.1
        metrics = [
            f"Domains Found: {quality_result.total_domains}",
            f"Classified: {quality_result.classified_domains}",
            f"Coverage: {quality_result.coverage_percentage:.1f}%",
            f"Confidence: {quality_result.average_confidence:.2f}",
            f"Quality Tier: {quality_result.tier.upper()}"
        ]
        
        for metric in metrics:
            ax.text(0.1, y_pos, metric, transform=ax.transAxes, fontsize=12)
            y_pos -= 0.08
        
        # Comparison section if available
        if comparison_result:
            y_pos -= 0.05
            ax.text(0.05, y_pos, "COMPARISON WITH MAIN", transform=ax.transAxes,
                   fontsize=14, fontweight='bold', color='darkred')
            
            y_pos -= 0.1
            comparison_metrics = [
                f"Main Domains: {comparison_result.main_domains}",
                f"Winner: {comparison_result.winner.upper()}",
                f"Score Improvement: {comparison_result.score_difference:+.3f}",
                f"Improvement Type: {comparison_result.improvement_type}"
            ]
            
            for metric in comparison_metrics:
                ax.text(0.1, y_pos, metric, transform=ax.transAxes, fontsize=12)
                y_pos -= 0.08
        
        # Quality visualization (simple gauge)
        if quality_result.overall_score > 0:
            # Create a simple quality gauge on the right side
            gauge_x = 0.7
            gauge_y = 0.6
            gauge_radius = 0.15
            
            # Background circle
            circle = plt.Circle((gauge_x, gauge_y), gauge_radius, 
                              color='lightgray', alpha=0.3, transform=ax.transAxes)
            ax.add_patch(circle)
            
            # Score arc
            score = quality_result.overall_score
            angle = score * 180  # Half circle
            
            if score >= 0.8:
                color = 'green'
            elif score >= 0.6:
                color = 'orange'
            else:
                color = 'red'
            
            # Draw arc (simplified)
            ax.text(gauge_x, gauge_y - 0.25, f"Quality Score\n{score:.2f}", 
                   transform=ax.transAxes, ha='center', fontsize=12, fontweight='bold')
        
        # PDB info
        ax.text(0.05, 0.05, f"PDB: {quality_result.pdb_id} | Chain: {quality_result.chain_id} | "
                           f"Batch: {quality_result.batch_name}",
               transform=ax.transAxes, fontsize=10, style='italic')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / f"{protein_id}_summary_card.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        # Return data for JSON export
        summary_data = {
            "protein_id": protein_id,
            "pdb_id": quality_result.pdb_id,
            "mini_results": {
                "domains": quality_result.total_domains,
                "classified": quality_result.classified_domains,
                "coverage": quality_result.coverage_percentage,
                "confidence": quality_result.average_confidence,
                "tier": quality_result.tier,
                "score": quality_result.overall_score
            }
        }
        
        if comparison_result:
            summary_data["comparison"] = {
                "main_domains": comparison_result.main_domains,
                "winner": comparison_result.winner,
                "score_improvement": comparison_result.score_difference,
                "improvement_type": comparison_result.improvement_type
            }
        
        return summary_data
    
    def _create_pymol_example(self, protein_id: str, batch_dir: Path) -> Optional[str]:
        """Create PyMOL visualization script"""
        
        if not PYMOL_AVAILABLE:
            return None
        
        try:
            # Look for old domains file
            old_file = batch_dir / "domains" / f"{protein_id}.develop291.domains.xml"
            new_file = batch_dir / "mini_domains" / f"{protein_id}.mini.domains.xml"
            
            if not old_file.exists() or not new_file.exists():
                logger.warning(f"Missing domain files for PyMOL comparison: {protein_id}")
                return None
            
            script_path = self.pymol_viz.create_comparison(
                protein_id, str(old_file), str(new_file), str(self.output_dir)
            )
            
            logger.info(f"PyMOL script created: {script_path}")
            return script_path
            
        except Exception as e:
            logger.error(f"Failed to create PyMOL script for {protein_id}: {e}")
            return None
    
    def showcase_category_examples(self, category: str, count: int = 3) -> List[Dict]:
        """Generate examples from a specific quality category"""
        
        logger.info(f"Showcasing {count} examples from '{category}' category")
        
        # Get all quality results
        all_results = self.quality_assessor.assess_all_batches()
        
        # Filter by category
        category_results = [r for r in all_results if r.tier == category]
        
        if not category_results:
            logger.warning(f"No results found for category '{category}'")
            return []
        
        # Sort by score and take top examples
        category_results.sort(key=lambda x: x.overall_score, reverse=True)
        top_examples = category_results[:count]
        
        # Generate visualizations for each
        showcase_data = []
        for result in top_examples:
            viz_data = self.visualize_protein_example(result.protein_id, create_pymol=True)
            if viz_data:
                showcase_data.append(viz_data)
        
        # Create category summary
        self._create_category_summary(category, showcase_data)
        
        return showcase_data
    
    def showcase_improvement_examples(self, count: int = 5) -> List[Dict]:
        """Generate examples showing biggest improvements over main algorithm"""
        
        logger.info(f"Showcasing top {count} improvement examples")
        
        # Get comparison results
        comparison_results = self.comparator.sample_comparison(1000)  # Large sample
        
        # Filter to mini wins with significant improvements
        improvements = [c for c in comparison_results 
                       if c.winner == 'mini' and c.score_difference > 0.1]
        
        # Sort by score improvement
        improvements.sort(key=lambda x: x.score_difference, reverse=True)
        top_improvements = improvements[:count]
        
        # Generate visualizations
        showcase_data = []
        for comp in top_improvements:
            viz_data = self.visualize_protein_example(comp.protein_id, create_pymol=True)
            if viz_data:
                viz_data["comparison_highlight"] = {
                    "score_improvement": comp.score_difference,
                    "improvement_type": comp.improvement_type,
                    "mini_domains": comp.mini_domains,
                    "main_domains": comp.main_domains
                }
                showcase_data.append(viz_data)
        
        # Create improvement summary
        self._create_improvement_summary(showcase_data)
        
        return showcase_data
    
    def _create_category_summary(self, category: str, showcase_data: List[Dict]):
        """Create summary visualization for category showcase"""
        
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.axis('off')
        
        # Title
        ax.text(0.5, 0.95, f"{category.upper()} QUALITY EXAMPLES", 
               transform=ax.transAxes, ha='center', fontsize=18, fontweight='bold')
        
        # Examples grid
        y_start = 0.8
        for i, viz_data in enumerate(showcase_data):
            y_pos = y_start - (i * 0.25)
            
            protein_id = viz_data["protein_id"]
            quality = viz_data["quality_result"]
            
            # Protein name
            ax.text(0.05, y_pos, f"{i+1}. {protein_id}", transform=ax.transAxes,
                   fontsize=14, fontweight='bold')
            
            # Metrics
            y_pos -= 0.05
            metrics_text = (f"Domains: {quality.total_domains} | "
                          f"Classified: {quality.classified_domains} | "
                          f"Coverage: {quality.coverage_percentage:.1f}% | "
                          f"Score: {quality.overall_score:.3f}")
            
            ax.text(0.1, y_pos, metrics_text, transform=ax.transAxes, fontsize=11)
            
            # Comparison if available
            if viz_data.get("comparison_result"):
                comp = viz_data["comparison_result"]
                y_pos -= 0.04
                comp_text = f"vs Main: {comp.mini_domains} vs {comp.main_domains} domains, {comp.winner} wins"
                ax.text(0.1, y_pos, comp_text, transform=ax.transAxes, 
                       fontsize=10, style='italic', color='darkblue')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / f"{category}_examples_summary.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Category summary created for {category}")
    
    def _create_improvement_summary(self, showcase_data: List[Dict]):
        """Create summary visualization for improvement examples"""
        
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.axis('off')
        
        # Title
        ax.text(0.5, 0.95, "TOP ALGORITHM IMPROVEMENTS", 
               transform=ax.transAxes, ha='center', fontsize=18, fontweight='bold')
        
        # Examples grid
        y_start = 0.8
        for i, viz_data in enumerate(showcase_data):
            y_pos = y_start - (i * 0.15)
            
            protein_id = viz_data["protein_id"]
            highlight = viz_data.get("comparison_highlight", {})
            
            # Protein name and improvement
            improvement = highlight.get("score_improvement", 0)
            ax.text(0.05, y_pos, f"{i+1}. {protein_id} (+{improvement:.3f})", 
                   transform=ax.transAxes, fontsize=14, fontweight='bold')
            
            # Details
            y_pos -= 0.05
            mini_d = highlight.get("mini_domains", 0)
            main_d = highlight.get("main_domains", 0)
            imp_type = highlight.get("improvement_type", "unknown")
            
            details_text = f"Mini: {mini_d} domains vs Main: {main_d} domains | Improvement: {imp_type}"
            ax.text(0.1, y_pos, details_text, transform=ax.transAxes, fontsize=11)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / "improvement_examples_summary.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info("Improvement summary created")


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Generate Example Protein Visualizations for Slide Deck'
    )
    
    parser.add_argument('--protein', type=str,
                       help='Generate visualization for specific protein')
    parser.add_argument('--comparison', action='store_true',
                       help='Include comparison with main algorithm')
    parser.add_argument('--showcase-category', type=str,
                       choices=['excellent', 'good', 'acceptable', 'poor'],
                       help='Showcase examples from quality category')
    parser.add_argument('--improvement-examples', action='store_true',
                       help='Showcase top improvement examples')
    parser.add_argument('--count', type=int, default=3,
                       help='Number of examples to showcase')
    parser.add_argument('--top', type=int, default=5,
                       help='Number of top improvements to show')
    
    parser.add_argument('--output-dir', type=str, default='/tmp/example_viz',
                       help='Output directory for visualizations')
    parser.add_argument('--config', type=str, default='config.local.yml',
                       help='Config file path')
    
    args = parser.parse_args()
    
    # Initialize visualizer
    visualizer = ExampleVisualizer(args.config, args.output_dir)
    
    if args.protein:
        # Single protein visualization
        viz_data = visualizer.visualize_protein_example(args.protein, args.comparison)
        if viz_data:
            print(f"✅ Visualization created for {args.protein}")
            print(f"📁 Output directory: {args.output_dir}")
            for file in viz_data["visualization_files"]:
                print(f"📊 {file}")
        else:
            print(f"❌ Failed to create visualization for {args.protein}")
    
    elif args.showcase_category:
        # Category showcase
        showcase_data = visualizer.showcase_category_examples(args.showcase_category, args.count)
        print(f"✅ Created showcase for {len(showcase_data)} {args.showcase_category} examples")
        print(f"📁 Output directory: {args.output_dir}")
        
    elif args.improvement_examples:
        # Improvement showcase
        showcase_data = visualizer.showcase_improvement_examples(args.top)
        print(f"✅ Created showcase for {len(showcase_data)} improvement examples")
        print(f"📁 Output directory: {args.output_dir}")
        
    else:
        parser.print_help()
        print(f"\n💡 Examples:")
        print(f"  python example_visualizations.py --protein 8ovp_A --comparison")
        print(f"  python example_visualizations.py --showcase-category excellent --count 3")
        print(f"  python example_visualizations.py --improvement-examples --top 5")


if __name__ == "__main__":
    main()
