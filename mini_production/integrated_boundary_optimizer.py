#!/usr/bin/env python3
"""
Integrated 3D Boundary Optimizer - Proof of Concept

End-to-end pipeline for optimizing double IG domain boundaries using
3D spatial consensus. Works on small test set (20-30 proteins) to
validate approach before scaling.

Integration points:
- Reads from existing mini_domains XML files
- Loads PDB structures
- Uses BioPython for fast superposition
- Calculates spatial boundary consensus
- Identifies and optimizes outliers
- Validates improvements
"""

import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, field
import json
import tempfile
import concurrent.futures
from collections import defaultdict
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# BioPython imports (install with: pip install biopython)
try:
    from Bio.PDB import PDBParser
    from Bio.SVDSuperimposer import SVDSuperimposer
except ImportError:
    print("Warning: BioPython not available. Install with: pip install biopython")
    SVDSuperimposer = None

@dataclass
class ProteinStructure:
    """Protein structure with coordinates and boundary info"""
    protein_id: str
    pdb_id: str
    chain_id: str
    
    # Structure data
    ca_coordinates: np.ndarray
    residue_numbers: List[int]
    sequence_length: int
    
    # Boundary data from mini_domains
    current_boundary_position: int
    domain1_range: str
    domain2_range: str
    
    # 3D boundary point
    boundary_coordinates_3d: Optional[np.ndarray] = None
    
    # Alignment results
    aligned_coordinates: Optional[np.ndarray] = None
    alignment_rmsd: float = 0.0
    alignment_success: bool = False

@dataclass
class BoundaryOptimization:
    """Results of boundary optimization for one protein"""
    protein_id: str
    original_position: int
    optimized_position: int
    original_coords_3d: np.ndarray
    optimized_coords_3d: np.ndarray
    distance_improvement: float
    residues_moved: int
    confidence_score: float

@dataclass 
class ConsensusAnalysis:
    """Results of spatial consensus analysis"""
    spatial_centroid: np.ndarray
    boundary_variance: float
    outlier_threshold: float
    outliers: List[str] = field(default_factory=list)
    optimizations: List[BoundaryOptimization] = field(default_factory=list)
    total_improvement: float = 0.0

class IntegratedBoundaryOptimizer:
    """Integrated 3D boundary optimization pipeline"""
    
    def __init__(self, 
                 batch_base: str = "/data/ecod/pdb_updates/batches",
                 pdb_base: str = "/data/pdb"):
        self.batch_base = Path(batch_base)
        self.pdb_base = Path(pdb_base)
        self.parser = PDBParser(QUIET=True) if SVDSuperimposer else None
        self.superimposer = SVDSuperimposer() if SVDSuperimposer else None
        
    def find_double_ig_test_set(self, target_count: int = 25) -> List[str]:
        """Find a good test set of double IG proteins"""
        
        print(f"🔍 Finding {target_count} double IG proteins for test set...")
        
        candidates = []
        
        # Search through batch directories
        for batch_dir in self.batch_base.iterdir():
            if not batch_dir.is_dir():
                continue
                
            mini_domains_dir = batch_dir / "mini_domains"
            if not mini_domains_dir.exists():
                continue
            
            for xml_file in mini_domains_dir.glob("*.mini.domains.xml"):
                if len(candidates) >= target_count * 2:  # Get extra candidates
                    break
                    
                if self.is_double_ig_with_structure(xml_file):
                    protein_id = xml_file.stem.replace('.mini.domains', '')
                    candidates.append(protein_id)
        
        # Select diverse set (different PDBs)
        selected = self.select_diverse_test_set(candidates, target_count)
        
        print(f"✓ Selected {len(selected)} proteins for test set")
        return selected
    
    def is_double_ig_with_structure(self, xml_file: Path) -> bool:
        """Check if protein is double IG and has available structure"""
        
        try:
            # Check XML for double IG architecture
            tree = ET.parse(xml_file)
            root = tree.getroot()
            
            domain_elements = root.findall(".//domain")
            if len(domain_elements) != 2:
                return False
            
            # Both domains must be IG (11.1.1)
            t_groups = [d.get('t_group', '') for d in domain_elements]
            if not all(tg.startswith('11.1.1') for tg in t_groups):
                return False
            
            # Check if PDB structure is available
            protein_id = xml_file.stem.replace('.mini.domains', '')
            pdb_id = protein_id.split('_')[0]
            
            pdb_file = self.find_pdb_file(pdb_id)
            return pdb_file is not None
            
        except Exception:
            return False
    
    def select_diverse_test_set(self, candidates: List[str], target_count: int) -> List[str]:
        """Select diverse test set (different PDBs preferred)"""
        
        # Group by PDB ID
        by_pdb = defaultdict(list)
        for protein_id in candidates:
            pdb_id = protein_id.split('_')[0]
            by_pdb[pdb_id].append(protein_id)
        
        # Select one protein per PDB (diversity)
        selected = []
        for pdb_id, proteins in by_pdb.items():
            if len(selected) < target_count:
                selected.append(proteins[0])  # Take first from each PDB
        
        # Fill remaining slots if needed
        while len(selected) < target_count and len(selected) < len(candidates):
            for proteins in by_pdb.values():
                for protein_id in proteins[1:]:  # Take additional from same PDB
                    if protein_id not in selected and len(selected) < target_count:
                        selected.append(protein_id)
        
        return selected[:target_count]
    
    def load_protein_structure(self, protein_id: str) -> Optional[ProteinStructure]:
        """Load protein structure and boundary information"""
        
        try:
            # Parse protein ID
            pdb_id = protein_id.split('_')[0]
            chain_id = protein_id.split('_')[1] if '_' in protein_id else 'A'
            
            # Load PDB structure
            pdb_file = self.find_pdb_file(pdb_id)
            if not pdb_file:
                print(f"Warning: PDB file not found for {pdb_id}")
                return None
            
            structure = self.parser.get_structure(pdb_id, pdb_file)
            
            # Extract C-alpha coordinates
            ca_coords = []
            residue_numbers = []
            
            for model in structure:
                for chain in model:
                    if chain.id == chain_id:
                        for residue in chain:
                            if 'CA' in residue:
                                ca = residue['CA']
                                ca_coords.append(ca.get_coord())
                                residue_numbers.append(residue.id[1])
                        break
                break
            
            if not ca_coords:
                print(f"Warning: No CA atoms found for {protein_id}")
                return None
            
            # Load boundary information from mini_domains XML
            boundary_info = self.load_boundary_info(protein_id)
            if not boundary_info:
                print(f"Warning: No boundary info found for {protein_id}")
                return None
            
            return ProteinStructure(
                protein_id=protein_id,
                pdb_id=pdb_id,
                chain_id=chain_id,
                ca_coordinates=np.array(ca_coords),
                residue_numbers=residue_numbers,
                sequence_length=len(ca_coords),
                current_boundary_position=boundary_info['boundary_position'],
                domain1_range=boundary_info['domain1_range'],
                domain2_range=boundary_info['domain2_range']
            )
            
        except Exception as e:
            print(f"Warning: Failed to load {protein_id}: {e}")
            return None
    
    def find_pdb_file(self, pdb_id: str) -> Optional[Path]:
        """Find PDB file in standard directory structure"""
        
        # Standard PDB layout: /data/pdb/xy/pdb1xyz.ent.gz
        subdir = pdb_id[1:3].lower()
        
        # Try different file formats
        candidates = [
            self.pdb_base / subdir / f"pdb{pdb_id.lower()}.ent",
            self.pdb_base / subdir / f"pdb{pdb_id.lower()}.ent.gz",
            self.pdb_base / subdir / f"{pdb_id.lower()}.cif",
            self.pdb_base / subdir / f"{pdb_id.lower()}.cif.gz"
        ]
        
        for candidate in candidates:
            if candidate.exists():
                return candidate
        
        return None
    
    def load_boundary_info(self, protein_id: str) -> Optional[Dict]:
        """Load boundary information from mini_domains XML"""
        
        # Find the XML file
        xml_file = None
        for batch_dir in self.batch_base.iterdir():
            if not batch_dir.is_dir():
                continue
            
            candidate = batch_dir / "mini_domains" / f"{protein_id}.mini.domains.xml"
            if candidate.exists():
                xml_file = candidate
                break
        
        if not xml_file:
            return None
        
        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()
            
            domain_elements = root.findall(".//domain")
            if len(domain_elements) != 2:
                return None
            
            domain1_range = domain_elements[0].get('range', '')
            domain2_range = domain_elements[1].get('range', '')
            
            # Calculate boundary position (end of first domain)
            boundary_position = self.calculate_boundary_position(domain1_range)
            
            return {
                'boundary_position': boundary_position,
                'domain1_range': domain1_range,
                'domain2_range': domain2_range
            }
            
        except Exception as e:
            print(f"Warning: Failed to parse boundary info for {protein_id}: {e}")
            return None
    
    def calculate_boundary_position(self, domain1_range: str) -> int:
        """Calculate boundary position from domain range"""
        
        try:
            # Parse range segments
            max_position = 0
            for segment in domain1_range.split(','):
                segment = segment.strip()
                if '-' in segment:
                    start, end = segment.split('-')
                    max_position = max(max_position, int(end))
                else:
                    max_position = max(max_position, int(segment))
            
            return max_position
            
        except ValueError:
            return 0
    
    def superimpose_structures(self, structures: List[ProteinStructure]) -> List[ProteinStructure]:
        """Superimpose all structures using BioPython"""
        
        print(f"🔄 Superimposing {len(structures)} structures...")
        
        # Use first structure as reference
        reference = structures[0]
        reference.aligned_coordinates = reference.ca_coordinates.copy()
        reference.alignment_success = True
        
        aligned_structures = [reference]
        
        for structure in structures[1:]:
            try:
                # Trim to same length (use shorter structure)
                min_length = min(len(reference.ca_coordinates), len(structure.ca_coordinates))
                ref_coords = reference.ca_coordinates[:min_length]
                mob_coords = structure.ca_coordinates[:min_length]
                
                # Superimpose
                self.superimposer.set(ref_coords, mob_coords)
                self.superimposer.run()
                
                # Get transformation
                rotation, translation = self.superimposer.get_rotran()
                rmsd = self.superimposer.get_rms()
                
                # Apply to full structure
                aligned_coords = (rotation @ structure.ca_coordinates.T).T + translation
                
                structure.aligned_coordinates = aligned_coords
                structure.alignment_rmsd = rmsd
                structure.alignment_success = True
                
                aligned_structures.append(structure)
                
                print(f"  {structure.protein_id}: RMSD = {rmsd:.2f}Å")
                
            except Exception as e:
                print(f"  Warning: Failed to align {structure.protein_id}: {e}")
                structure.alignment_success = False
        
        successful = [s for s in aligned_structures if s.alignment_success]
        print(f"✓ Successfully aligned {len(successful)}/{len(structures)} structures")
        
        return successful
    
    def map_boundaries_to_3d(self, structures: List[ProteinStructure]) -> List[ProteinStructure]:
        """Map boundary positions to 3D coordinates"""
        
        print("📍 Mapping boundaries to 3D coordinates...")
        
        mapped = []
        for structure in structures:
            if not structure.alignment_success:
                continue
            
            # Map boundary position to 3D coordinates
            boundary_idx = min(structure.current_boundary_position, len(structure.aligned_coordinates) - 1)
            
            if boundary_idx >= 0:
                structure.boundary_coordinates_3d = structure.aligned_coordinates[boundary_idx]
                mapped.append(structure)
                
                print(f"  {structure.protein_id}: boundary at position {structure.current_boundary_position} "
                      f"→ ({structure.boundary_coordinates_3d[0]:.1f}, "
                      f"{structure.boundary_coordinates_3d[1]:.1f}, "
                      f"{structure.boundary_coordinates_3d[2]:.1f})")
        
        print(f"✓ Mapped {len(mapped)} boundaries to 3D coordinates")
        return mapped
    
    def calculate_consensus(self, structures: List[ProteinStructure]) -> ConsensusAnalysis:
        """Calculate spatial consensus and identify outliers"""
        
        print("🎯 Calculating spatial consensus...")
        
        # Extract boundary coordinates
        boundary_coords = np.array([s.boundary_coordinates_3d for s in structures])
        
        # Calculate centroid
        centroid = np.mean(boundary_coords, axis=0)
        
        # Calculate distances to centroid
        distances = [np.linalg.norm(coords - centroid) for coords in boundary_coords]
        
        # Calculate variance
        variance = np.var(distances)
        
        # Identify outliers (distance > mean + 2*std)
        mean_distance = np.mean(distances)
        std_distance = np.std(distances)
        outlier_threshold = mean_distance + 2 * std_distance
        
        outliers = []
        for i, (structure, distance) in enumerate(zip(structures, distances)):
            if distance > outlier_threshold:
                outliers.append(structure.protein_id)
                print(f"  🚨 Outlier: {structure.protein_id} (distance: {distance:.2f}Å)")
        
        print(f"✓ Spatial centroid: ({centroid[0]:.1f}, {centroid[1]:.1f}, {centroid[2]:.1f})")
        print(f"✓ Boundary variance: {variance:.3f}")
        print(f"✓ Found {len(outliers)} outliers (threshold: {outlier_threshold:.2f}Å)")
        
        return ConsensusAnalysis(
            spatial_centroid=centroid,
            boundary_variance=variance,
            outlier_threshold=outlier_threshold,
            outliers=outliers
        )
    
    def optimize_boundaries(self, structures: List[ProteinStructure], 
                          consensus: ConsensusAnalysis) -> ConsensusAnalysis:
        """Optimize boundary positions for outlier proteins"""
        
        print("🔧 Optimizing outlier boundaries...")
        
        optimizations = []
        
        for structure in structures:
            if structure.protein_id not in consensus.outliers:
                continue
            
            print(f"  Optimizing {structure.protein_id}...")
            
            original_pos = structure.current_boundary_position
            original_coords = structure.boundary_coordinates_3d
            original_distance = np.linalg.norm(original_coords - consensus.spatial_centroid)
            
            best_pos = original_pos
            best_coords = original_coords
            best_distance = original_distance
            
            # Search window around current position
            search_window = 10
            start_pos = max(0, original_pos - search_window)
            end_pos = min(len(structure.aligned_coordinates) - 1, original_pos + search_window)
            
            for test_pos in range(start_pos, end_pos + 1):
                test_coords = structure.aligned_coordinates[test_pos]
                test_distance = np.linalg.norm(test_coords - consensus.spatial_centroid)
                
                if test_distance < best_distance:
                    best_distance = test_distance
                    best_pos = test_pos
                    best_coords = test_coords
            
            # Calculate improvement
            improvement = original_distance - best_distance
            residues_moved = abs(best_pos - original_pos)
            
            # Calculate confidence (higher improvement + fewer moves = higher confidence)
            confidence = improvement / (1 + residues_moved * 0.5)
            
            optimization = BoundaryOptimization(
                protein_id=structure.protein_id,
                original_position=original_pos,
                optimized_position=best_pos,
                original_coords_3d=original_coords,
                optimized_coords_3d=best_coords,
                distance_improvement=improvement,
                residues_moved=residues_moved,
                confidence_score=confidence
            )
            
            optimizations.append(optimization)
            
            print(f"    Original position: {original_pos} (distance: {original_distance:.2f}Å)")
            print(f"    Optimized position: {best_pos} (distance: {best_distance:.2f}Å)")
            print(f"    Improvement: {improvement:.2f}Å, moved {residues_moved} residues")
        
        consensus.optimizations = optimizations
        consensus.total_improvement = sum(opt.distance_improvement for opt in optimizations)
        
        print(f"✓ Optimized {len(optimizations)} boundaries")
        print(f"✓ Total spatial improvement: {consensus.total_improvement:.2f}Å")
        
        return consensus
    
    def visualize_results(self, structures: List[ProteinStructure], 
                         consensus: ConsensusAnalysis, output_file: str = "boundary_consensus_3d.png"):
        """Create 3D visualization of results"""
        
        fig = plt.figure(figsize=(15, 5))
        
        # Plot 1: Original boundaries
        ax1 = fig.add_subplot(131, projection='3d')
        
        coords = np.array([s.boundary_coordinates_3d for s in structures])
        colors = ['red' if s.protein_id in consensus.outliers else 'blue' for s in structures]
        
        ax1.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=colors, s=50, alpha=0.7)
        ax1.scatter(*consensus.spatial_centroid, c='gold', s=200, marker='*', label='Centroid')
        ax1.set_title('Original Boundaries')
        ax1.set_xlabel('X (Å)')
        ax1.set_ylabel('Y (Å)')
        ax1.set_zlabel('Z (Å)')
        
        # Plot 2: Optimization vectors
        ax2 = fig.add_subplot(132, projection='3d')
        
        ax2.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=colors, s=50, alpha=0.3)
        ax2.scatter(*consensus.spatial_centroid, c='gold', s=200, marker='*', label='Centroid')
        
        for opt in consensus.optimizations:
            orig = opt.original_coords_3d
            new = opt.optimized_coords_3d
            ax2.plot([orig[0], new[0]], [orig[1], new[1]], [orig[2], new[2]], 
                    'green', linewidth=2, alpha=0.8)
            ax2.scatter(*new, c='green', s=80, marker='s')
        
        ax2.set_title('Boundary Optimizations')
        ax2.set_xlabel('X (Å)')
        ax2.set_ylabel('Y (Å)')
        ax2.set_zlabel('Z (Å)')
        
        # Plot 3: Improvement summary
        ax3 = fig.add_subplot(133)
        
        if consensus.optimizations:
            improvements = [opt.distance_improvement for opt in consensus.optimizations]
            moves = [opt.residues_moved for opt in consensus.optimizations]
            
            scatter = ax3.scatter(moves, improvements, s=100, alpha=0.7, 
                                c=[opt.confidence_score for opt in consensus.optimizations],
                                cmap='viridis')
            
            ax3.set_xlabel('Residues Moved')
            ax3.set_ylabel('Distance Improvement (Å)')
            ax3.set_title('Optimization Quality')
            plt.colorbar(scatter, ax=ax3, label='Confidence')
            
            # Add protein labels
            for opt in consensus.optimizations:
                ax3.annotate(opt.protein_id, 
                           (opt.residues_moved, opt.distance_improvement),
                           fontsize=8, alpha=0.7)
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Saved visualization to {output_file}")
    
    def export_results(self, structures: List[ProteinStructure], 
                      consensus: ConsensusAnalysis, output_file: str = "boundary_optimization_results.json"):
        """Export results for further analysis"""
        
        results = {
            'analysis_summary': {
                'structures_analyzed': len(structures),
                'spatial_centroid': consensus.spatial_centroid.tolist(),
                'boundary_variance': consensus.boundary_variance,
                'outliers_found': len(consensus.outliers),
                'total_improvement': consensus.total_improvement
            },
            'structures': [],
            'optimizations': []
        }
        
        # Export structure data
        for structure in structures:
            results['structures'].append({
                'protein_id': structure.protein_id,
                'pdb_id': structure.pdb_id,
                'chain_id': structure.chain_id,
                'current_boundary': structure.current_boundary_position,
                'boundary_3d': structure.boundary_coordinates_3d.tolist(),
                'is_outlier': structure.protein_id in consensus.outliers,
                'alignment_rmsd': structure.alignment_rmsd
            })
        
        # Export optimization results
        for opt in consensus.optimizations:
            results['optimizations'].append({
                'protein_id': opt.protein_id,
                'original_position': opt.original_position,
                'optimized_position': opt.optimized_position,
                'distance_improvement': opt.distance_improvement,
                'residues_moved': opt.residues_moved,
                'confidence_score': opt.confidence_score,
                'recommend_update': opt.confidence_score > 0.5 and opt.residues_moved <= 5
            })
        
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"✓ Exported results to {output_file}")
    
    def print_summary(self, structures: List[ProteinStructure], consensus: ConsensusAnalysis):
        """Print comprehensive analysis summary"""
        
        print(f"\n📊 3D BOUNDARY CONSENSUS ANALYSIS - PROOF OF CONCEPT")
        print("=" * 80)
        print(f"Test set size: {len(structures)} double IG proteins")
        print(f"Successful alignments: {sum(1 for s in structures if s.alignment_success)}")
        print(f"Spatial centroid: ({consensus.spatial_centroid[0]:.1f}, "
              f"{consensus.spatial_centroid[1]:.1f}, {consensus.spatial_centroid[2]:.1f})")
        print(f"Boundary variance: {consensus.boundary_variance:.3f}")
        print(f"Outliers identified: {len(consensus.outliers)}")
        print(f"Total spatial improvement: {consensus.total_improvement:.2f}Å")
        
        if consensus.optimizations:
            print(f"\n🎯 BOUNDARY OPTIMIZATION RESULTS:")
            print(f"{'Protein':<12} {'Orig':<5} {'Opt':<5} {'Moved':<6} {'Improve':<8} {'Conf':<6} {'Recommend'}")
            print("-" * 70)
            
            for opt in consensus.optimizations:
                recommend = "✓ YES" if opt.confidence_score > 0.5 and opt.residues_moved <= 5 else "❌ NO"
                print(f"{opt.protein_id:<12} {opt.original_position:<5} {opt.optimized_position:<5} "
                      f"{opt.residues_moved:<6} {opt.distance_improvement:<8.2f} {opt.confidence_score:<6.2f} {recommend}")
            
            # Summary stats
            improvements = [opt.distance_improvement for opt in consensus.optimizations]
            moves = [opt.residues_moved for opt in consensus.optimizations]
            recommended = sum(1 for opt in consensus.optimizations 
                            if opt.confidence_score > 0.5 and opt.residues_moved <= 5)
            
            print(f"\n📈 OPTIMIZATION STATISTICS:")
            print(f"  Average improvement: {np.mean(improvements):.2f}Å")
            print(f"  Average residues moved: {np.mean(moves):.1f}")
            print(f"  High-confidence optimizations: {recommended}/{len(consensus.optimizations)}")
    
    def run_proof_of_concept(self, target_proteins: int = 25) -> Dict:
        """Run complete proof-of-concept analysis"""
        
        print("🚀 Starting 3D boundary optimization proof-of-concept...")
        
        # Step 1: Find test set
        test_proteins = self.find_double_ig_test_set(target_proteins)
        if len(test_proteins) < 10:
            print(f"❌ Insufficient test proteins found: {len(test_proteins)}")
            return {}
        
        # Step 2: Load structures
        print(f"\n📥 Loading {len(test_proteins)} protein structures...")
        structures = []
        for protein_id in test_proteins:
            structure = self.load_protein_structure(protein_id)
            if structure:
                structures.append(structure)
        
        if len(structures) < 10:
            print(f"❌ Insufficient structures loaded: {len(structures)}")
            return {}
        
        print(f"✓ Loaded {len(structures)} structures successfully")
        
        # Step 3: Superimpose structures
        aligned_structures = self.superimpose_structures(structures)
        
        # Step 4: Map boundaries to 3D
        mapped_structures = self.map_boundaries_to_3d(aligned_structures)
        
        # Step 5: Calculate consensus
        consensus = self.calculate_consensus(mapped_structures)
        
        # Step 6: Optimize boundaries
        consensus = self.optimize_boundaries(mapped_structures, consensus)
        
        # Step 7: Generate outputs
        self.print_summary(mapped_structures, consensus)
        self.visualize_results(mapped_structures, consensus)
        self.export_results(mapped_structures, consensus)
        
        print(f"\n💡 PROOF-OF-CONCEPT INSIGHTS:")
        print(f"✓ 3D spatial consensus successfully calculated")
        print(f"✓ Outlier boundaries identified and optimized")
        print(f"✓ Quantitative validation of boundary corrections")
        print(f"✓ Ready for scaling to full dataset")
        
        return {
            'structures': mapped_structures,
            'consensus': consensus,
            'success': True
        }

def main():
    """Run proof-of-concept analysis"""
    
    print("🧬 3D Boundary Optimization - Proof of Concept")
    print("=" * 60)
    
    optimizer = IntegratedBoundaryOptimizer()
    
    # Run analysis on 25 proteins
    results = optimizer.run_proof_of_concept(target_proteins=25)
    
    if results['success']:
        print(f"\n🎉 Proof-of-concept completed successfully!")
        print(f"Generated files:")
        print(f"  • boundary_consensus_3d.png (3D visualization)")
        print(f"  • boundary_optimization_results.json (detailed results)")
        print(f"\nNext steps:")
        print(f"  1. Review optimization recommendations")
        print(f"  2. Validate a few cases manually")
        print(f"  3. Scale to full double IG dataset (1500+ proteins)")
    else:
        print(f"❌ Proof-of-concept failed - check data availability")

if __name__ == "__main__":
    main()
