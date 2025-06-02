#!/usr/bin/env python3
"""
PyMOL Domain Comparison Visualization

Creates side-by-side comparison of domain assignments:
- Left: Previous algorithm results
- Right: New mini_pyecod results

Usage:
    python pymol_comparison.py 8ovp_A
    python pymol_comparison.py 8ovp_A --batch-mode
"""

import sys
import os
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass

@dataclass
class Domain:
    """Domain representation for visualization"""
    id: str
    range_str: str
    family: str
    color: str
    segments: List[Tuple[int, int]]

class PDBRepository:
    """Handle local PDB repository access"""
    
    def __init__(self, pdb_root: str = "/usr2/pdb/data"):
        self.pdb_root = pdb_root
    
    def find_structure(self, pdb_id: str) -> Optional[str]:
        """Find mmCIF file for PDB ID in local repository"""
        pdb_id = pdb_id.lower()
        
        # Local repository structure: /usr2/pdb/data/structures/divided/mmCIF/ov/8ovp.cif.gz
        possible_paths = [
            f"{self.pdb_root}/structures/divided/mmCIF/{pdb_id[1:3]}/{pdb_id}.cif.gz",
            f"{self.pdb_root}/structures/divided/mmCIF/{pdb_id[1:3]}/{pdb_id}.cif",
            # Fallback paths
            f"{self.pdb_root}/mmCIF/{pdb_id}.cif.gz",
            f"{self.pdb_root}/mmCIF/{pdb_id}.cif",
            f"{self.pdb_root}/{pdb_id}.cif.gz",
            f"{self.pdb_root}/{pdb_id}.cif",
        ]
        
        for path in possible_paths:
            if os.path.exists(path):
                return path
        
        return None

class DomainParser:
    """Parse domain assignments from XML files"""
    
    @staticmethod
    def parse_new_domains(xml_path: str) -> List[Domain]:
        """Parse domains from our mini_pyecod output XML"""
        domains = []
        
        try:
            tree = ET.parse(xml_path)
            root = tree.getroot()
            
            # Same color scheme as old domains for consistency
            base_colors = ["red", "blue", "green", "cyan", "magenta", "orange"]
            variant_colors = ["salmon", "lightblue", "lightgreen", "lightcyan", "pink", "gold"]
            all_colors = base_colors + variant_colors
            
            for i, domain_elem in enumerate(root.findall(".//domain")):
                domain_id = domain_elem.get("id", f"d{i+1}")
                range_str = domain_elem.get("range", "")
                family = domain_elem.get("family", "unknown")
                
                # Color rotation
                color = all_colors[i % len(all_colors)]
                
                # Parse range into segments
                segments = DomainParser._parse_range_segments(range_str)
                
                domains.append(Domain(
                    id=domain_id,
                    range_str=range_str,
                    family=family,
                    color=color,
                    segments=segments
                ))
                
        except Exception as e:
            print(f"Error parsing new domains from {xml_path}: {e}")
        
        return domains
    
    @staticmethod
    def parse_old_domains(xml_path: str) -> List[Domain]:
        """Parse domains from previous algorithm output XML (develop291 format)"""
        domains = []
        
        try:
            tree = ET.parse(xml_path)
            root = tree.getroot()
            
            # Extended color scheme: 6 base + 6 variations = 12 colors, then rotate
            base_colors = ["red", "blue", "green", "cyan", "magenta", "orange"]
            variant_colors = ["salmon", "lightblue", "lightgreen", "lightcyan", "pink", "gold"]
            all_colors = base_colors + variant_colors
            
            # Parse old format: <domain start="5" end="245" range="5-245" ...>
            for i, domain_elem in enumerate(root.findall(".//domain")):
                range_str = domain_elem.get("range", "")
                source = domain_elem.get("source", "unknown")
                source_id = domain_elem.get("source_id", "")
                t_group = domain_elem.get("t_group", "")
                
                # Create domain ID from available info
                domain_id = f"old_d{i+1}"
                if source_id:
                    domain_id = source_id
                elif t_group:
                    domain_id = t_group
                
                # Family name from source_id or t_group
                family = source_id or t_group or source or "unknown"
                
                # Color rotation
                color = all_colors[i % len(all_colors)]
                
                # Parse range into segments
                segments = DomainParser._parse_range_segments(range_str)
                
                domains.append(Domain(
                    id=domain_id,
                    range_str=range_str,
                    family=family,
                    color=color,
                    segments=segments
                ))
                
        except Exception as e:
            print(f"Error parsing old domains from {xml_path}: {e}")
        
        return domains
    
    @staticmethod
    def _parse_range_segments(range_str: str) -> List[Tuple[int, int]]:
        """Parse range string into list of (start, end) tuples"""
        segments = []
        
        if not range_str:
            return segments
        
        try:
            for segment in range_str.split(','):
                segment = segment.strip()
                if '-' in segment:
                    start, end = segment.split('-')
                    segments.append((int(start), int(end)))
                else:
                    # Single residue
                    pos = int(segment)
                    segments.append((pos, pos))
        except ValueError as e:
            print(f"Warning: Could not parse range '{range_str}': {e}")
        
        return segments

class PyMOLVisualizer:
    """Generate PyMOL commands for domain visualization"""
    
    def __init__(self, pdb_repo: PDBRepository):
        self.pdb_repo = pdb_repo
        
    def create_comparison_script(self, protein_id: str, 
                               old_domains: List[Domain],
                               new_domains: List[Domain],
                               output_dir: str = "/tmp/pymol_comparison") -> str:
        """Create PyMOL script for side-by-side comparison"""
        
        pdb_id = protein_id.split('_')[0]
        chain_id = protein_id.split('_')[1] if '_' in protein_id else 'A'
        
        # Find structure file
        structure_path = self.pdb_repo.find_structure(pdb_id)
        if not structure_path:
            raise FileNotFoundError(f"Structure not found for {pdb_id}")
        
        os.makedirs(output_dir, exist_ok=True)
        script_path = os.path.join(output_dir, f"{protein_id}_comparison.pml")
        
        with open(script_path, 'w') as f:
            f.write(self._generate_pymol_commands(
                protein_id, chain_id, structure_path, old_domains, new_domains
            ))
        
        return script_path
    
    def _generate_pymol_commands(self, protein_id: str, chain_id: str,
                               structure_path: str,
                               old_domains: List[Domain],
                               new_domains: List[Domain]) -> str:
        """Generate the actual PyMOL command script with grid mode"""
        
        commands = []
        pdb_id = protein_id.split('_')[0]
        
        # Header
        commands.append(f"# PyMOL comparison script for {protein_id}")
        commands.append(f"# Generated by mini_pyecod visualization system")
        commands.append(f"# OLD: {len(old_domains)} domains vs NEW: {len(new_domains)} domains")
        commands.append("")
        
        # Load structure
        commands.append("# Load structure")
        commands.append(f"load {structure_path}, {pdb_id}")
        commands.append("")
        
        # Create two copies for comparison
        commands.append("# Create two copies for comparison")
        commands.append(f"create old_algorithm, {pdb_id}")
        commands.append(f"create new_algorithm, {pdb_id}")
        commands.append(f"delete {pdb_id}")
        commands.append("")
        
        # Hide other chains for clarity
        commands.append("# Hide other chains for clarity")
        commands.append(f"hide everything, not chain {chain_id}")
        commands.append("")
        
        # Basic styling
        commands.append("# Basic styling")
        commands.append("hide everything")
        commands.append("show cartoon")
        commands.append("set cartoon_side_chain_helper, on")
        commands.append("set cartoon_fancy_helices, 1")
        commands.append("color gray90, all")
        commands.append("")
        
        # Color old domains
        commands.append("# Color OLD algorithm domains")
        commands.append(f"# Found {len(old_domains)} domains")
        for i, domain in enumerate(old_domains):
            commands.append(f"# Domain {i+1}: {domain.family} ({domain.range_str})")
            for start, end in domain.segments:
                selection = f"old_algorithm and chain {chain_id} and resi {start}-{end}"
                commands.append(f"color {domain.color}, {selection}")
        commands.append("")
        
        # Color new domains  
        commands.append("# Color NEW algorithm domains")
        commands.append(f"# Found {len(new_domains)} domains")
        for i, domain in enumerate(new_domains):
            commands.append(f"# Domain {i+1}: {domain.family} ({domain.range_str})")
            for start, end in domain.segments:
                selection = f"new_algorithm and chain {chain_id} and resi {start}-{end}"
                commands.append(f"color {domain.color}, {selection}")
        commands.append("")
        
        # Set up grid mode for side-by-side comparison
        commands.append("# Set up grid mode for side-by-side comparison")
        commands.append("set grid_mode, 1")
        commands.append("set grid_slot, 1, old_algorithm")
        commands.append("set grid_slot, 2, new_algorithm")
        commands.append("")
        
        # Orient both structures the same way
        commands.append("# Orient structures")
        commands.append("orient old_algorithm")
        commands.append("orient new_algorithm")
        commands.append("")
        
        # Add text labels for clarity
        commands.append("# Add informative labels")
        commands.append("set label_size, 20")
        commands.append("set label_color, black")
        commands.append("")
        commands.append("# Create label objects")
        commands.append("pseudoatom old_label")
        commands.append("pseudoatom new_label")
        commands.append(f"label old_label, 'Previous Algorithm ({len(old_domains)} domains)'")
        commands.append(f"label new_label, 'Mini PyECOD ({len(new_domains)} domains)'")
        commands.append("set grid_slot, 1, old_label")
        commands.append("set grid_slot, 2, new_label")
        commands.append("")
        
        # Final view settings
        commands.append("# Final view settings")
        commands.append("zoom all")
        commands.append("set depth_cue, 0")
        commands.append("set ray_shadows, 0")
        commands.append("bg_color white")
        commands.append("")
        
        # Save session
        commands.append("# Save session")
        session_path = f"/tmp/pymol_comparison/{protein_id}_comparison.pse"
        commands.append(f"save {session_path}")
        commands.append("")
        
        # Create domain legend as comments
        commands.append("# DOMAIN COMPARISON SUMMARY")
        commands.append(f"# Protein: {protein_id}")
        commands.append("#")
        commands.append(f"# OLD ALGORITHM ({len(old_domains)} domains):")
        for i, domain in enumerate(old_domains, 1):
            commands.append(f"#   {i}. {domain.family}: {domain.range_str} ({domain.color})")
        commands.append("#")
        commands.append(f"# NEW ALGORITHM ({len(new_domains)} domains):")
        for i, domain in enumerate(new_domains, 1):
            commands.append(f"#   {i}. {domain.family}: {domain.range_str} ({domain.color})")
        commands.append("#")
        commands.append("# Expected improvement: More biologically meaningful domains")
        commands.append("# with better coverage and fewer spurious assignments")
        
        return "\n".join(commands)

def create_comparison_visualization(protein_id: str,
                                  batch_dir: str,
                                  pdb_repo_path: str = "/usr2/pdb/data",
                                  output_dir: str = "/tmp/pymol_comparison") -> str:
    """Main function to create comparison visualization"""
    
    print(f"Creating comparison visualization for {protein_id}")
    
    # Initialize components
    pdb_repo = PDBRepository(pdb_repo_path)
    parser = DomainParser()
    visualizer = PyMOLVisualizer(pdb_repo)
    
    # File paths
    new_domains_file = f"/tmp/{protein_id}_range_cache.domains.xml"  # Our output
    old_domains_file = os.path.join(batch_dir, "domains", f"{protein_id}.develop291.domains.xml")  # Previous algorithm
    
    # Check files exist
    if not os.path.exists(new_domains_file):
        raise FileNotFoundError(f"New domains file not found: {new_domains_file}")
    
    if not os.path.exists(old_domains_file):
        raise FileNotFoundError(f"Old domains file not found: {old_domains_file}")
    
    # Parse domains
    print("Parsing domain assignments...")
    old_domains = parser.parse_old_domains(old_domains_file)
    new_domains = parser.parse_new_domains(new_domains_file)
    
    print(f"Old algorithm: {len(old_domains)} domains")
    for i, domain in enumerate(old_domains, 1):
        print(f"  {i}. {domain.family}: {domain.range_str}")
    
    print(f"New algorithm: {len(new_domains)} domains")
    for i, domain in enumerate(new_domains, 1):
        print(f"  {i}. {domain.family}: {domain.range_str}")
    
    # Create PyMOL script
    script_path = visualizer.create_comparison_script(
        protein_id, old_domains, new_domains, output_dir
    )
    
    print(f"\nPyMOL script created: {script_path}")
    print(f"Session will be saved to: /tmp/pymol_comparison/{protein_id}_comparison.pse")
    print(f"\nTo view: pymol {script_path}")
    
    return script_path

def main():
    """Main CLI interface"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Create PyMOL domain comparison visualization')
    parser.add_argument('protein_id', help='Protein ID (e.g., 8ovp_A)')
    parser.add_argument('--batch-dir', 
                        default='/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424',
                        help='Batch directory containing old domain files')
    parser.add_argument('--pdb-repo', default='/data/pdb',
                        help='PDB repository path')
    parser.add_argument('--old-suffix', default='domains.xml',
                        help='Suffix for old domain files')
    parser.add_argument('--output-dir', default='/tmp/pymol_comparison',
                        help='Output directory for PyMOL scripts')
    
    args = parser.parse_args()
    
    try:
        script_path = create_comparison_visualization(
            args.protein_id,
            args.batch_dir,
            args.pdb_repo,
            args.old_suffix,
            args.output_dir
        )
        
        print(f"\n✅ Visualization ready!")
        print(f"Run: pymol {script_path}")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
