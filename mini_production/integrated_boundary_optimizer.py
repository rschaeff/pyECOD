#!/usr/bin/env python3
"""
Integrated 3D Boundary Optimizer - Uses Existing Architecture Analysis

Integrates with existing analyze_domain_architectures.py to find double IG proteins
and performs 3D spatial consensus optimization on their boundaries.

Usage:
    python integrated_boundary_optimizer.py --architecture "11.1.1 → 11.1.1" --max-proteins 25
    python integrated_boundary_optimizer.py --auto-find-double-ig
"""

import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import json
import argparse
import tempfile
import concurrent.futures
from collections import defaultdict
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Import the existing architecture analyzer
import sys
sys.path.append(str(Path(__file__).parent))
from analyze_domain_architectures import DomainArchitectureAnalyzer, DomainArchitecture, ProteinResult

# BioPython imports
try:
    from Bio.PDB import MMCIFParser, PDBParser
    from Bio.SVDSuperimposer import SVDSuperimposer
    BIOPYTHON_AVAILABLE = True
except ImportError:
    print("Warning: BioPython not available. Install with: pip install biopython")
    BIOPYTHON_AVAILABLE = False

@dataclass
class Protein3DStructure:
    """3D structure with boundary information from existing analysis"""
    protein_result: ProteinResult  # From existing analyzer
    ca_coordinates: Optional[np.ndarray] = None
    aligned_coordinates: Optional[np.ndarray] = None
    boundary_position_seq: Optional[int] = None  # Sequence position
    boundary_coordinates_3d: Optional[np.ndarray] = None  # 3D coordinates
    alignment_rmsd: float = 0.0
    structure_loaded: bool = False
    alignment_success: bool = False

@dataclass
class BoundaryOptimization:
    """Results of 3D boundary optimization"""
    protein_id: str
    original_position: int
    optimized_position: int
    original_coords_3d: np.ndarray
    optimized_coords_3d: np.ndarray
    distance_improvement: float
    residues_moved: int
    confidence_score: float

@dataclass
class SpatialConsensusResults:
    """Results of spatial consensus analysis"""
    architecture: str
    total_proteins: int
    structures_loaded: int
    aligned_successfully: int
    spatial_centroid: Optional[np.ndarray] = None
    boundary_variance: float = 0.0
    outlier_threshold: float = 0.0
    outliers: List[str] = None
    optimizations: List[BoundaryOptimization] = None
    total_improvement: float = 0.0

    def __post_init__(self):
        if self.outliers is None:
            self.outliers = []
        if self.optimizations is None:
            self.optimizations = []

class IntegratedBoundaryOptimizer:
    """3D boundary optimization using existing architecture analysis"""

    def __init__(self, config_path: str = "config/config.local.yml", pdb_base: str = "/usr2/pdb/data"):
        self.architecture_analyzer = DomainArchitectureAnalyzer(config_path)
        self.pdb_base = Path(pdb_base)

        if BIOPYTHON_AVAILABLE:
            self.mmcif_parser = MMCIFParser(QUIET=True)
            self.pdb_parser = PDBParser(QUIET=True)  # Fallback
            self.superimposer = SVDSuperimposer()
        else:
            self.mmcif_parser = None
            self.pdb_parser = None
            self.superimposer = None
            print("Warning: BioPython not available - structure loading disabled")

    def find_double_ig_architecture(self) -> Optional[DomainArchitecture]:
        """Find the double IG architecture using existing analyzer"""

        print("🔍 Scanning for double IG proteins using existing architecture analyzer...")

        # Load all results using existing analyzer
        results = self.architecture_analyzer.scan_mini_results()

        if not results:
            print("❌ No mini results found")
            return None

        # Get two-domain architectures from excellent tier
        two_domain_architectures = self.architecture_analyzer.get_two_domain_architectures(results, "excellent")

        # Look for double IG architecture (11.1.1 → 11.1.1)
        double_ig_candidates = []
        for arch_str, architecture in two_domain_architectures.items():
            t_groups = architecture.t_group_sequence
            if len(t_groups) == 2 and all(tg and tg.startswith('11.1.1') for tg in t_groups):
                double_ig_candidates.append((arch_str, architecture))

        if not double_ig_candidates:
            print("❌ No double IG architectures found in excellent tier")
            print("Available two-domain architectures:")
            for arch_str in list(two_domain_architectures.keys())[:10]:
                print(f"  {arch_str}")
            return None

        # Take the most common double IG architecture
        double_ig_candidates.sort(key=lambda x: len(x[1].examples), reverse=True)
        arch_str, architecture = double_ig_candidates[0]

        print(f"✓ Found double IG architecture: {arch_str}")
        print(f"✓ Found {len(architecture.examples)} examples in excellent tier")

        return architecture

    def find_specific_architecture(self, architecture_string: str) -> Optional[DomainArchitecture]:
        """Find specific architecture using existing analyzer"""

        print(f"🔍 Looking for architecture: {architecture_string}")

        # Load all results
        results = self.architecture_analyzer.scan_mini_results()

        if not results:
            print("❌ No mini results found")
            return None

        # Get all excellent tier architectures
        excellent_architectures = self.architecture_analyzer.analyze_tier_architectures(results, "excellent")

        if architecture_string in excellent_architectures:
            architecture = excellent_architectures[architecture_string]
            print(f"✓ Found architecture: {architecture_string}")
            print(f"✓ Found {len(architecture.examples)} examples in excellent tier")
            return architecture
        else:
            print(f"❌ Architecture '{architecture_string}' not found in excellent tier")
            print("Available architectures:")
            for arch_str in list(excellent_architectures.keys())[:10]:
                print(f"  {arch_str}")
            return None

    def calculate_boundary_position(self, protein_result: ProteinResult) -> Optional[int]:
        """Calculate boundary position between first and second domain"""

        if len(protein_result.domains) < 2:
            return None

        # Get the range of the first domain
        first_domain = protein_result.domains[0]

        if not first_domain.segments:
            return None

        # Boundary is at the end of the first domain
        boundary_position = max(end for start, end in first_domain.segments)
        return boundary_position

    def find_pdb_structure(self, pdb_id: str) -> Optional[Tuple[Path, str]]:
        """Find PDB structure file - prefer mmCIF format

        Returns:
            Tuple of (file_path, format) where format is 'cif' or 'pdb'
        """

        pdb_id = pdb_id.lower()

        # Prefer mmCIF format - it's the modern standard and handles edge cases better
        possible_paths = [
            # mmCIF format (preferred)
            (self.pdb_base / "structures" / "divided" / "mmCIF" / pdb_id[1:3] / f"{pdb_id}.cif.gz", 'cif'),
            (self.pdb_base / "structures" / "divided" / "mmCIF" / pdb_id[1:3] / f"{pdb_id}.cif", 'cif'),
            (self.pdb_base / f"{pdb_id}.cif.gz", 'cif'),
            (self.pdb_base / f"{pdb_id}.cif", 'cif'),
            # PDB format (fallback)
            (self.pdb_base / "structures" / "divided" / "pdb" / pdb_id[1:3] / f"pdb{pdb_id}.ent.gz", 'pdb'),
            (self.pdb_base / "structures" / "divided" / "pdb" / pdb_id[1:3] / f"pdb{pdb_id}.ent", 'pdb'),
            (self.pdb_base / f"pdb{pdb_id}.ent.gz", 'pdb'),
            (self.pdb_base / f"pdb{pdb_id}.ent", 'pdb'),
        ]

        for path, format_type in possible_paths:
            if path.exists():
                return path, format_type

        return None

    def load_structure_coordinates(self, protein_result: ProteinResult) -> Optional[Protein3DStructure]:
        """Load 3D coordinates using modern mmCIF format with fallback to PDB"""

        if not BIOPYTHON_AVAILABLE:
            return None

        try:
            structure_info = self.find_pdb_structure(protein_result.pdb_id)
            if not structure_info:
                print(f"  ❌ {protein_result.protein_id}: No structure file found")
                return None

            pdb_file, file_format = structure_info

            # Log which format we're using
            print(f"  📁 {protein_result.protein_id}: Loading {file_format.upper()} format from {pdb_file.name}")

            # Choose appropriate parser
            if file_format == 'cif':
                parser = self.mmcif_parser
            else:
                parser = self.pdb_parser

            # Handle compressed files
            if pdb_file.suffix == '.gz':
                import gzip
                with gzip.open(pdb_file, 'rt') as f:
                    content = f.read()

                # Save to temporary file
                if file_format == 'cif':
                    suffix = '.cif'
                else:
                    suffix = '.pdb'

                with tempfile.NamedTemporaryFile(mode='w', suffix=suffix, delete=False) as tmp:
                    tmp.write(content)
                    tmp_path = tmp.name

                try:
                    structure = parser.get_structure(protein_result.pdb_id, tmp_path)
                finally:
                    Path(tmp_path).unlink()  # Clean up
            else:
                structure = parser.get_structure(protein_result.pdb_id, pdb_file)

            # Extract C-alpha coordinates
            ca_coords = []
            residue_numbers = []

            for model in structure:
                for chain in model:
                    if chain.id == protein_result.chain_id:
                        for residue in chain:
                            try:
                                # Skip non-standard residues (HETATMs)
                                if residue.id[0] != ' ':
                                    continue

                                if 'CA' not in residue:
                                    continue

                                ca = residue['CA']
                                coords = ca.get_coord()

                                # Get residue number - mmCIF handles this more robustly
                                res_num = residue.id[1]
                                if isinstance(res_num, int) and res_num > 0:
                                    ca_coords.append(coords)
                                    residue_numbers.append(res_num)

                            except (ValueError, KeyError, AttributeError):
                                # Skip problematic residues - mmCIF should minimize these
                                continue
                        break
                break

            if len(ca_coords) < 50:  # Need reasonable number of residues
                print(f"  ❌ {protein_result.protein_id}: Only {len(ca_coords)} residues found (need ≥50)")
                return None

            # Calculate boundary position
            boundary_position = self.calculate_boundary_position(protein_result)

            # Log coordinate range for debugging
            coords_array = np.array(ca_coords)
            coord_ranges = {
                'x': (coords_array[:, 0].min(), coords_array[:, 0].max()),
                'y': (coords_array[:, 1].min(), coords_array[:, 1].max()),
                'z': (coords_array[:, 2].min(), coords_array[:, 2].max())
            }

            print(f"  📊 {protein_result.protein_id}: {len(ca_coords)} residues, "
                  f"coord ranges: X({coord_ranges['x'][0]:.1f}-{coord_ranges['x'][1]:.1f}), "
                  f"Y({coord_ranges['y'][0]:.1f}-{coord_ranges['y'][1]:.1f}), "
                  f"Z({coord_ranges['z'][0]:.1f}-{coord_ranges['z'][1]:.1f})")

            protein_3d = Protein3DStructure(
                protein_result=protein_result,
                ca_coordinates=coords_array,
                boundary_position_seq=boundary_position,
                structure_loaded=True
            )

            return protein_3d

        except Exception as e:
            print(f"  ❌ {protein_result.protein_id}: Failed to load - {e}")
            return None

    def load_structures_parallel(self, protein_results: List[ProteinResult], max_proteins: int = 25) -> List[Protein3DStructure]:
        """Load structures in parallel using modern mmCIF format"""

        # Take the first max_proteins * 1.5 to allow for some failures
        proteins_to_load = protein_results[:int(max_proteins * 1.5)]

        print(f"📥 Loading {len(proteins_to_load)} protein structures...")

        structures = []
        format_counts = {'cif': 0, 'pdb': 0}

        # Use parallel loading with reasonable batch size
        batch_size = 10

        for i in range(0, len(proteins_to_load), batch_size):
            if len(structures) >= max_proteins:
                break

            batch = proteins_to_load[i:i+batch_size]

            # Load sequentially within batch to see format logging clearly
            for protein in batch:
                if len(structures) >= max_proteins:
                    break

                structure = self.load_structure_coordinates(protein)
                if structure:
                    structures.append(structure)

                    # Track format usage
                    structure_info = self.find_pdb_structure(protein.pdb_id)
                    if structure_info:
                        format_counts[structure_info[1]] += 1

        print(f"\n📊 FORMAT USAGE SUMMARY:")
        print(f"  mmCIF files loaded: {format_counts['cif']}")
        print(f"  PDB files loaded: {format_counts['pdb']}")
        print(f"  Total structures loaded: {len(structures)}")

        if format_counts['pdb'] > format_counts['cif']:
            print(f"  ⚠️  WARNING: More PDB than mmCIF files used - may indicate coordinate issues")

        if len(structures) < 5:
            print(f"⚠️  Only {len(structures)} structures loaded - trying more proteins...")

            # Try a larger pool if we don't have enough
            if len(proteins_to_load) < len(protein_results):
                additional_proteins = protein_results[len(proteins_to_load):len(proteins_to_load) + max_proteins]
                additional_structures = self.load_structures_parallel(additional_proteins, max_proteins - len(structures))
                structures.extend(additional_structures)

        return structures

    def superimpose_structures(self, structures: List[Protein3DStructure]) -> List[Protein3DStructure]:
        """Superimpose all structures using BioPython"""

        if not BIOPYTHON_AVAILABLE or len(structures) < 2:
            return structures

        print(f"🔄 Superimposing {len(structures)} structures...")

        # Check coordinate systems before superposition
        print(f"📊 PRE-SUPERPOSITION COORDINATE ANALYSIS:")
        for structure in structures[:3]:  # Check first 3
            coords = structure.ca_coordinates
            print(f"  {structure.protein_result.protein_id}: "
                  f"X({coords[:, 0].min():.1f}-{coords[:, 0].max():.1f}), "
                  f"Y({coords[:, 1].min():.1f}-{coords[:, 1].max():.1f}), "
                  f"Z({coords[:, 2].min():.1f}-{coords[:, 2].max():.1f})")

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

                print(f"  ✓ {structure.protein_result.protein_id}: RMSD = {rmsd:.2f}Å")

            except Exception as e:
                print(f"  ❌ Failed to align {structure.protein_result.protein_id}: {e}")
                structure.alignment_success = False

        successful = [s for s in aligned_structures if s.alignment_success]

        # Check coordinate systems after superposition
        print(f"📊 POST-SUPERPOSITION COORDINATE ANALYSIS:")
        for structure in successful[:3]:  # Check first 3
            coords = structure.aligned_coordinates
            print(f"  {structure.protein_result.protein_id}: "
                  f"X({coords[:, 0].min():.1f}-{coords[:, 0].max():.1f}), "
                  f"Y({coords[:, 1].min():.1f}-{coords[:, 1].max():.1f}), "
                  f"Z({coords[:, 2].min():.1f}-{coords[:, 2].max():.1f})")

        print(f"✓ Successfully aligned {len(successful)}/{len(structures)} structures")

        return successful

    def map_boundaries_to_3d(self, structures: List[Protein3DStructure]) -> List[Protein3DStructure]:
        """Map sequence boundaries to 3D coordinates"""

        print("📍 Mapping boundaries to 3D coordinates...")

        mapped = []
        boundary_coords_list = []

        for structure in structures:
            if not structure.alignment_success or structure.boundary_position_seq is None:
                continue

            # Map boundary position to 3D coordinates
            boundary_idx = min(structure.boundary_position_seq - 1,  # Convert to 0-based
                             len(structure.aligned_coordinates) - 1)
            boundary_idx = max(0, boundary_idx)  # Ensure non-negative

            structure.boundary_coordinates_3d = structure.aligned_coordinates[boundary_idx]
            boundary_coords_list.append(structure.boundary_coordinates_3d)
            mapped.append(structure)

            print(f"  ✓ {structure.protein_result.protein_id}: position {structure.boundary_position_seq} "
                  f"→ ({structure.boundary_coordinates_3d[0]:.1f}, "
                  f"{structure.boundary_coordinates_3d[1]:.1f}, "
                  f"{structure.boundary_coordinates_3d[2]:.1f})")

        # Analyze boundary coordinate distribution
        if boundary_coords_list:
            boundary_array = np.array(boundary_coords_list)
            print(f"\n📊 BOUNDARY COORDINATE ANALYSIS:")
            print(f"  X range: {boundary_array[:, 0].min():.1f} - {boundary_array[:, 0].max():.1f} "
                  f"(spread: {boundary_array[:, 0].max() - boundary_array[:, 0].min():.1f}Å)")
            print(f"  Y range: {boundary_array[:, 1].min():.1f} - {boundary_array[:, 1].max():.1f} "
                  f"(spread: {boundary_array[:, 1].max() - boundary_array[:, 1].min():.1f}Å)")
            print(f"  Z range: {boundary_array[:, 2].min():.1f} - {boundary_array[:, 2].max():.1f} "
                  f"(spread: {boundary_array[:, 2].max() - boundary_array[:, 2].min():.1f}Å)")

            # Check for obviously wrong coordinates
            max_spread = max(
                boundary_array[:, 0].max() - boundary_array[:, 0].min(),
                boundary_array[:, 1].max() - boundary_array[:, 1].min(),
                boundary_array[:, 2].max() - boundary_array[:, 2].min()
            )

            if max_spread > 300:
                print(f"  ⚠️  WARNING: Large coordinate spread ({max_spread:.1f}Å) suggests alignment/coordinate issues")

        print(f"✓ Mapped {len(mapped)} boundaries to 3D coordinates")
        return mapped

    def calculate_spatial_consensus(self, structures: List[Protein3DStructure]) -> SpatialConsensusResults:
        """Calculate spatial consensus of boundary positions"""

        print("🎯 Calculating spatial consensus...")

        if not structures:
            return SpatialConsensusResults(
                architecture="unknown",
                total_proteins=0,
                structures_loaded=0,
                aligned_successfully=0
            )

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
        for structure, distance in zip(structures, distances):
            if distance > outlier_threshold:
                outliers.append(structure.protein_result.protein_id)
                print(f"  🚨 Outlier: {structure.protein_result.protein_id} (distance: {distance:.2f}Å)")

        architecture_str = structures[0].protein_result.architecture_string

        print(f"✓ Spatial centroid: ({centroid[0]:.1f}, {centroid[1]:.1f}, {centroid[2]:.1f})")
        print(f"✓ Boundary variance: {variance:.3f}")
        print(f"✓ Found {len(outliers)} outliers (threshold: {outlier_threshold:.2f}Å)")

        return SpatialConsensusResults(
            architecture=architecture_str,
            total_proteins=len(structures),
            structures_loaded=len(structures),
            aligned_successfully=len(structures),
            spatial_centroid=centroid,
            boundary_variance=variance,
            outlier_threshold=outlier_threshold,
            outliers=outliers
        )

    def optimize_outlier_boundaries(self, structures: List[Protein3DStructure],
                                  consensus: SpatialConsensusResults) -> SpatialConsensusResults:
        """Optimize boundary positions for outlier proteins"""

        print("🔧 Optimizing outlier boundaries...")

        optimizations = []

        for structure in structures:
            if structure.protein_result.protein_id not in consensus.outliers:
                continue

            print(f"  Optimizing {structure.protein_result.protein_id}...")

            original_pos = structure.boundary_position_seq
            original_coords = structure.boundary_coordinates_3d
            original_distance = np.linalg.norm(original_coords - consensus.spatial_centroid)

            best_pos = original_pos
            best_coords = original_coords
            best_distance = original_distance

            # Search window around current position
            search_window = 10
            start_pos = max(1, original_pos - search_window)
            end_pos = min(len(structure.aligned_coordinates), original_pos + search_window)

            for test_pos in range(start_pos, end_pos + 1):
                test_idx = min(test_pos - 1, len(structure.aligned_coordinates) - 1)
                test_idx = max(0, test_idx)
                test_coords = structure.aligned_coordinates[test_idx]
                test_distance = np.linalg.norm(test_coords - consensus.spatial_centroid)

                if test_distance < best_distance:
                    best_distance = test_distance
                    best_pos = test_pos
                    best_coords = test_coords

            # Calculate improvement metrics
            improvement = original_distance - best_distance
            residues_moved = abs(best_pos - original_pos)
            confidence = improvement / (1 + residues_moved * 0.5)

            optimization = BoundaryOptimization(
                protein_id=structure.protein_result.protein_id,
                original_position=original_pos,
                optimized_position=best_pos,
                original_coords_3d=original_coords,
                optimized_coords_3d=best_coords,
                distance_improvement=improvement,
                residues_moved=residues_moved,
                confidence_score=confidence
            )

            optimizations.append(optimization)

            print(f"    Original: position {original_pos} (distance: {original_distance:.2f}Å)")
            print(f"    Optimized: position {best_pos} (distance: {best_distance:.2f}Å)")
            print(f"    Improvement: {improvement:.2f}Å, moved {residues_moved} residues")

        consensus.optimizations = optimizations
        consensus.total_improvement = sum(opt.distance_improvement for opt in optimizations)

        print(f"✓ Optimized {len(optimizations)} boundaries")
        print(f"✓ Total spatial improvement: {consensus.total_improvement:.2f}Å")

        return consensus

    def visualize_results(self, structures: List[Protein3DStructure],
                         consensus: SpatialConsensusResults,
                         output_file: str = "boundary_consensus_3d.png"):
        """Create 3D visualization of results"""

        fig = plt.figure(figsize=(15, 5))

        # Plot 1: Original boundaries
        ax1 = fig.add_subplot(131, projection='3d')

        coords = np.array([s.boundary_coordinates_3d for s in structures])
        colors = ['red' if s.protein_result.protein_id in consensus.outliers else 'blue'
                 for s in structures]

        ax1.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=colors, s=50, alpha=0.7)
        ax1.scatter(*consensus.spatial_centroid, c='gold', s=200, marker='*', label='Centroid')
        ax1.set_title(f'Original Boundaries\n{consensus.architecture}')
        ax1.set_xlabel('X (Å)')
        ax1.set_ylabel('Y (Å)')
        ax1.set_zlabel('Z (Å)')

        # Plot 2: Optimization vectors
        ax2 = fig.add_subplot(132, projection='3d')

        ax2.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=colors, s=50, alpha=0.3)
        ax2.scatter(*consensus.spatial_centroid, c='gold', s=200, marker='*')

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
                protein_short = opt.protein_id.split('_')[0]  # Just PDB ID
                ax3.annotate(protein_short,
                           (opt.residues_moved, opt.distance_improvement),
                           fontsize=8, alpha=0.7)

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Saved visualization to {output_file}")

    def export_results(self, consensus: SpatialConsensusResults,
                      output_file: str = "boundary_optimization_results.json"):
        """Export results for further analysis"""

        # Convert numpy types to native Python types for JSON serialization
        def convert_numpy(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, (np.float32, np.float64)):
                return float(obj)
            elif isinstance(obj, (np.int32, np.int64)):
                return int(obj)
            return obj

        results = {
            'analysis_summary': {
                'architecture': consensus.architecture,
                'total_proteins': int(consensus.total_proteins),
                'structures_loaded': int(consensus.structures_loaded),
                'aligned_successfully': int(consensus.aligned_successfully),
                'spatial_centroid': [float(x) for x in consensus.spatial_centroid] if consensus.spatial_centroid is not None else None,
                'boundary_variance': float(consensus.boundary_variance),
                'outliers_found': len(consensus.outliers),
                'total_improvement': float(consensus.total_improvement)
            },
            'optimizations': []
        }

        # Export optimization results with robust type conversion
        for opt in consensus.optimizations:
            # Convert all numpy types to native Python
            results['optimizations'].append({
                'protein_id': str(opt.protein_id),
                'original_position': int(opt.original_position),
                'optimized_position': int(opt.optimized_position),
                'distance_improvement': float(opt.distance_improvement),
                'residues_moved': int(opt.residues_moved),
                'confidence_score': float(opt.confidence_score),
                'recommend_update': bool(opt.confidence_score > 0.5 and opt.residues_moved <= 5)
            })

        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)

        print(f"✓ Exported results to {output_file}")

    def export_detailed_diagnostics(self, structures: List[Protein3DStructure],
                                   consensus: SpatialConsensusResults,
                                   output_file: str = "boundary_diagnostics.json"):
        """Export detailed diagnostics for analysis"""

        diagnostics = {
            'analysis_info': {
                'architecture': consensus.architecture,
                'total_structures': len(structures),
                'spatial_centroid': [float(x) for x in consensus.spatial_centroid] if consensus.spatial_centroid is not None else None,
                'boundary_variance': float(consensus.boundary_variance),
                'outlier_threshold': float(consensus.outlier_threshold)
            },
            'structure_details': [],
            'boundary_coordinates': [],
            'distance_from_centroid': []
        }

        # Export detailed structure information
        for structure in structures:
            details = {
                'protein_id': str(structure.protein_result.protein_id),
                'pdb_id': str(structure.protein_result.pdb_id),
                'chain_id': str(structure.protein_result.chain_id),
                'sequence_length': int(len(structure.ca_coordinates)),
                'boundary_position_seq': int(structure.boundary_position_seq) if structure.boundary_position_seq else None,
                'boundary_coordinates_3d': [float(x) for x in structure.boundary_coordinates_3d] if structure.boundary_coordinates_3d is not None else None,
                'alignment_rmsd': float(structure.alignment_rmsd),
                'domain1_range': str(structure.protein_result.domains[0].range_str) if structure.protein_result.domains else None,
                'domain2_range': str(structure.protein_result.domains[1].range_str) if len(structure.protein_result.domains) > 1 else None,
                'is_outlier': str(structure.protein_result.protein_id) in consensus.outliers
            }

            # Calculate distance from centroid
            if structure.boundary_coordinates_3d is not None and consensus.spatial_centroid is not None:
                distance = float(np.linalg.norm(structure.boundary_coordinates_3d - consensus.spatial_centroid))
                details['distance_from_centroid'] = distance
            else:
                details['distance_from_centroid'] = None

            diagnostics['structure_details'].append(details)

        with open(output_file, 'w') as f:
            json.dump(diagnostics, f, indent=2)

        print(f"✓ Exported detailed diagnostics to {output_file}")

        # Print some immediate insights
        print(f"\n🔍 DIAGNOSTIC INSIGHTS:")

        # RMSD distribution
        rmsds = [s.alignment_rmsd for s in structures]
        print(f"RMSD range: {min(rmsds):.1f} - {max(rmsds):.1f}Å (average: {np.mean(rmsds):.1f}Å)")

        # High RMSD structures (likely superposition problems)
        high_rmsd = [s for s in structures if s.alignment_rmsd > 10.0]
        print(f"High RMSD structures (>10Å): {len(high_rmsd)}/{len(structures)}")
        for s in high_rmsd[:5]:  # Show first 5
            print(f"  {s.protein_result.protein_id}: {s.alignment_rmsd:.1f}Å")

        # Boundary position distribution
        boundary_positions = [s.boundary_position_seq for s in structures if s.boundary_position_seq]
        if boundary_positions:
            print(f"Boundary positions: {min(boundary_positions)} - {max(boundary_positions)} (average: {np.mean(boundary_positions):.1f})")

        # Distance distribution
        distances = [details['distance_from_centroid'] for details in diagnostics['structure_details']
                    if details['distance_from_centroid'] is not None]
        if distances:
            print(f"Distances from centroid: {min(distances):.1f} - {max(distances):.1f}Å (average: {np.mean(distances):.1f}Å)")

        return diagnostics

    def print_summary(self, consensus: SpatialConsensusResults):
        """Print comprehensive analysis summary"""

        print(f"\n📊 3D BOUNDARY CONSENSUS ANALYSIS")
        print("=" * 80)
        print(f"Architecture: {consensus.architecture}")
        print(f"Total proteins: {consensus.total_proteins}")
        print(f"Structures loaded: {consensus.structures_loaded}")
        print(f"Successfully aligned: {consensus.aligned_successfully}")

        if consensus.spatial_centroid is not None:
            print(f"Spatial centroid: ({consensus.spatial_centroid[0]:.1f}, "
                  f"{consensus.spatial_centroid[1]:.1f}, {consensus.spatial_centroid[2]:.1f})")
            print(f"Boundary variance: {consensus.boundary_variance:.3f}")
            print(f"Outliers identified: {len(consensus.outliers)}")
            print(f"Total spatial improvement: {consensus.total_improvement:.2f}Å")

        if consensus.optimizations:
            print(f"\n🎯 BOUNDARY OPTIMIZATION RESULTS:")
            print(f"{'Protein':<15} {'Orig':<5} {'Opt':<5} {'Moved':<6} {'Improve':<8} {'Conf':<6} {'Recommend'}")
            print("-" * 75)

            for opt in consensus.optimizations:
                recommend = "✓ YES" if opt.confidence_score > 0.5 and opt.residues_moved <= 5 else "❌ NO"
                print(f"{opt.protein_id:<15} {opt.original_position:<5} {opt.optimized_position:<5} "
                      f"{opt.residues_moved:<6} {opt.distance_improvement:<8.2f} {opt.confidence_score:<6.2f} {recommend}")

    def run_analysis(self, architecture_string: Optional[str] = None, max_proteins: int = 25) -> SpatialConsensusResults:
        """Run complete 3D boundary optimization analysis"""

        print("🚀 Starting 3D boundary optimization analysis...")

        # Step 1: Find target architecture
        if architecture_string:
            architecture = self.find_specific_architecture(architecture_string)
        else:
            architecture = self.find_double_ig_architecture()

        if not architecture:
            return SpatialConsensusResults(
                architecture="not_found",
                total_proteins=0,
                structures_loaded=0,
                aligned_successfully=0
            )

        # Step 2: Load 3D structures
        structures = self.load_structures_parallel(architecture.examples, max_proteins)

        if len(structures) < 5:
            print(f"❌ Insufficient structures loaded: {len(structures)} (need at least 5)")
            print(f"💡 Attempting to load more structures...")

            # Try a larger pool if needed
            if len(architecture.examples) > max_proteins:
                additional_structures = self.load_structures_parallel(
                    architecture.examples[max_proteins:max_proteins*2],
                    max_proteins - len(structures)
                )
                structures.extend(additional_structures)

            if len(structures) < 5:
                print(f"❌ Still insufficient structures: {len(structures)}")
                return SpatialConsensusResults(
                    architecture=architecture.architecture_string,
                    total_proteins=len(architecture.examples),
                    structures_loaded=len(structures),
                    aligned_successfully=0
                )

        # Step 3: Superimpose structures
        aligned_structures = self.superimpose_structures(structures)

        # Step 3.5: Filter out poor alignments
        rmsd_threshold = 8.0  # Angstroms
        good_alignments = [s for s in aligned_structures if s.alignment_rmsd < rmsd_threshold]

        print(f"🔍 RMSD Filtering:")
        print(f"  Structures before filtering: {len(aligned_structures)}")
        print(f"  RMSD threshold: {rmsd_threshold}Å")
        print(f"  Structures after filtering: {len(good_alignments)}")

        if len(good_alignments) < 5:
            print(f"❌ Insufficient well-aligned structures: {len(good_alignments)} (need at least 5)")
            print(f"💡 Poor alignments (RMSD > {rmsd_threshold}Å):")
            poor_alignments = [s for s in aligned_structures if s.alignment_rmsd >= rmsd_threshold]
            for s in poor_alignments[:10]:  # Show first 10
                print(f"    {s.protein_result.protein_id}: {s.alignment_rmsd:.1f}Å")

            return SpatialConsensusResults(
                architecture=architecture.architecture_string,
                total_proteins=len(architecture.examples),
                structures_loaded=len(structures),
                aligned_successfully=len(good_alignments)
            )

        # Step 4: Map boundaries to 3D (using only well-aligned structures)
        mapped_structures = self.map_boundaries_to_3d(good_alignments)

        if len(mapped_structures) < 5:
            print(f"❌ Insufficient boundaries mapped: {len(mapped_structures)}")
            return SpatialConsensusResults(
                architecture=architecture.architecture_string,
                total_proteins=len(architecture.examples),
                structures_loaded=len(structures),
                aligned_successfully=len(aligned_structures)
            )

        # Step 5: Calculate spatial consensus
        consensus = self.calculate_spatial_consensus(mapped_structures)

        # Step 6: Optimize outlier boundaries
        consensus = self.optimize_outlier_boundaries(mapped_structures, consensus)

        # Step 7: Generate outputs
        self.print_summary(consensus)
        self.visualize_results(mapped_structures, consensus)
        self.export_results(consensus)
        self.export_detailed_diagnostics(mapped_structures, consensus, good_alignments, aligned_structures)

        return consensus


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='3D Boundary Optimization using Existing Architecture Analysis'
    )

    parser.add_argument('--architecture', type=str,
                       help='Specific architecture to analyze (e.g., "11.1.1 → 11.1.1")')
    parser.add_argument('--auto-find-double-ig', action='store_true',
                       help='Automatically find and analyze double IG architecture')
    parser.add_argument('--max-proteins', type=int, default=50,
                       help='Maximum proteins to analyze (will try more if needed)')
    parser.add_argument('--config', type=str, default='config/config.local.yml',
                       help='Config file path')
    parser.add_argument('--pdb-base', type=str, default='/usr2/pdb/data',
                       help='PDB repository base path')

    args = parser.parse_args()

    print("🧬 Integrated 3D Boundary Optimization")
    print("=" * 60)

    if not BIOPYTHON_AVAILABLE:
        print("❌ BioPython not available - cannot perform 3D analysis")
        print("Install with: pip install biopython")
        return

    # Initialize optimizer
    optimizer = IntegratedBoundaryOptimizer(args.config, args.pdb_base)

    # Run analysis
    if args.auto_find_double_ig:
        results = optimizer.run_analysis(None, args.max_proteins)
    elif args.architecture:
        results = optimizer.run_analysis(args.architecture, args.max_proteins)
    else:
        print("❌ Must specify --architecture or --auto-find-double-ig")
        return

    # Print final summary
    if results.aligned_successfully >= 5:
        print(f"\n🎉 Analysis completed successfully!")
        print(f"Architecture: {results.architecture}")
        print(f"Structures analyzed: {results.aligned_successfully}")

        if results.optimizations:
            recommended = sum(1 for opt in results.optimizations
                            if opt.confidence_score > 0.5 and opt.residues_moved <= 5)
            print(f"High-confidence optimizations: {recommended}/{len(results.optimizations)}")

        print(f"\nGenerated files:")
        print(f"  • boundary_consensus_3d.png (visualization)")
        print(f"  • boundary_optimization_results.json (detailed results)")
    else:
        print(f"❌ Analysis failed - insufficient data")
        print(f"Architecture: {results.architecture}")
        print(f"Structures loaded: {results.structures_loaded}")
        print(f"Successfully aligned: {results.aligned_successfully}")


if __name__ == "__main__":
    main()
