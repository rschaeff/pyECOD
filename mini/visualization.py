#!/usr/bin/env python3
"""
Clean PyMOL domain comparison visualization for mini_pyecod

Simplified interface for creating side-by-side domain comparisons.
"""

import os
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List, Tuple, Optional
from dataclasses import dataclass

@dataclass
class Domain:
    """Domain representation for visualization"""
    id: str
    range_str: str
    family: str
    segments: List[Tuple[int, int]]

class PyMOLVisualizer:
    """Clean PyMOL visualization for domain comparisons"""
    
    def __init__(self, pdb_repo_path: str = "/usr2/pdb/data"):
        self.pdb_repo_path = pdb_repo_path
        self.colors = [
            "red", "blue", "green", "cyan", "magenta", "orange",
            "salmon", "lightblue", "lightgreen", "lightcyan", "pink", "gold"
        ]
    
    def create_comparison(self, protein_id: str, 
                         old_domains_file: str,
                         new_domains_file: str,
                         output_dir: str = "/tmp/pymol_comparison") -> str:
        """
        Create PyMOL comparison between old and new domain assignments.
        
        Args:
            protein_id: Protein ID (e.g., "8ovp_A")
            old_domains_file: Path to old domain XML file
            new_domains_file: Path to new domain XML file  
            output_dir: Output directory for PyMOL script
            
        Returns:
            Path to generated PyMOL script
        """
        # Parse domains from both files
        old_domains = self._parse_domain_file(old_domains_file, is_old_format=True)
        new_domains = self._parse_domain_file(new_domains_file, is_old_format=False)
        
        # Find structure file
        pdb_id = protein_id.split('_')[0]
        structure_path = self._find_structure(pdb_id)
        if not structure_path:
            raise FileNotFoundError(f"Structure not found for {pdb_id}")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Generate PyMOL script
        script_path = os.path.join(output_dir, f"{protein_id}_comparison.pml")
        self._write_pymol_script(protein_id, structure_path, old_domains, 
                                new_domains, script_path, output_dir)
        
        return script_path
    
    def _parse_domain_file(self, xml_path: str, is_old_format: bool = False) -> List[Domain]:
        """Parse domain XML file"""
        domains = []
        
        try:
            tree = ET.parse(xml_path)
            root = tree.getroot()
            
            for i, domain_elem in enumerate(root.findall(".//domain")):
                if is_old_format:
                    # Old format: <domain range="5-245" source_id="e4q0cB4" ...>
                    range_str = domain_elem.get("range", "")
                    family = domain_elem.get("source_id", "") or domain_elem.get("t_group", "unknown")
                    domain_id = f"old_d{i+1}"
                else:
                    # New format: <domain id="d1" range="2-248" family="2ia4" ...>
                    domain_id = domain_elem.get("id", f"d{i+1}")
                    range_str = domain_elem.get("range", "")
                    family = domain_elem.get("family", "unknown")
                
                # Parse range into segments
                segments = self._parse_range_segments(range_str)
                
                domains.append(Domain(
                    id=domain_id,
                    range_str=range_str,
                    family=family,
                    segments=segments
                ))
                
        except Exception as e:
            print(f"Warning: Error parsing {xml_path}: {e}")
        
        return domains
    
    def _parse_range_segments(self, range_str: str) -> List[Tuple[int, int]]:
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
    
    def _find_structure(self, pdb_id: str) -> Optional[str]:
        """Find structure file in PDB repository"""
        pdb_id = pdb_id.lower()
        
        possible_paths = [
            f"{self.pdb_repo_path}/structures/divided/mmCIF/{pdb_id[1:3]}/{pdb_id}.cif.gz",
            f"{self.pdb_repo_path}/structures/divided/mmCIF/{pdb_id[1:3]}/{pdb_id}.cif",
            f"{self.pdb_repo_path}/mmCIF/{pdb_id}.cif.gz",
            f"{self.pdb_repo_path}/mmCIF/{pdb_id}.cif",
        ]
        
        for path in possible_paths:
            if os.path.exists(path):
                return path
        
        return None
    
    def _write_pymol_script(self, protein_id: str, structure_path: str,
                           old_domains: List[Domain], new_domains: List[Domain],
                           script_path: str, output_dir: str):
        """Write PyMOL comparison script"""
        
        pdb_id = protein_id.split('_')[0]
        chain_id = protein_id.split('_')[1] if '_' in protein_id else 'A'
        
        lines = [
            f"# PyMOL domain comparison for {protein_id}",
            f"# Old algorithm: {len(old_domains)} domains",
            f"# New algorithm: {len(new_domains)} domains",
            "",
            "# Load structure and create copies",
            f"load {structure_path}, {pdb_id}",
            f"create old_algorithm, {pdb_id}",
            f"create new_algorithm, {pdb_id}",
            f"delete {pdb_id}",
            "",
            "# Basic styling",
            "hide everything",
            "show cartoon",
            "set cartoon_fancy_helices, 1",
            "color gray90, all",
            "",
            "# Color old algorithm domains",
            f"# {len(old_domains)} domains from previous algorithm"
        ]
        
        # Color old domains
        for i, domain in enumerate(old_domains):
            color = self.colors[i % len(self.colors)]
            lines.append(f"# Domain {i+1}: {domain.family} ({domain.range_str})")
            for start, end in domain.segments:
                selection = f"old_algorithm and chain {chain_id} and resi {start}-{end}"
                lines.append(f"color {color}, {selection}")
        
        lines.extend([
            "",
            "# Color new algorithm domains",
            f"# {len(new_domains)} domains from mini_pyecod"
        ])
        
        # Color new domains
        for i, domain in enumerate(new_domains):
            color = self.colors[i % len(self.colors)]
            lines.append(f"# Domain {i+1}: {domain.family} ({domain.range_str})")
            for start, end in domain.segments:
                selection = f"new_algorithm and chain {chain_id} and resi {start}-{end}"
                lines.append(f"color {color}, {selection}")
        
        # Grid layout and final setup
        lines.extend([
            "",
            "# Grid layout for side-by-side comparison",
            "set grid_mode, 1",
            "set grid_slot, 1, old_algorithm",
            "set grid_slot, 2, new_algorithm",
            "",
            "# Final view",
            "orient",
            "zoom all",
            "bg_color white",
            "",
            f"# Save session",
            f"save {output_dir}/{protein_id}_comparison.pse"
        ])
        
        # Write script
        with open(script_path, 'w') as f:
            f.write('\n'.join(lines))
    
    def quick_comparison(self, protein_id: str, batch_dir: str) -> str:
        """
        Quick comparison using standard file locations.
        
        Args:
            protein_id: Protein ID (e.g., "8ovp_A")
            batch_dir: Batch directory containing old domain files
            
        Returns:
            Path to PyMOL script
        """
        old_file = f"{batch_dir}/domains/{protein_id}.develop291.domains.xml"
        new_file = f"/tmp/{protein_id}_mini.domains.xml"
        
        if not os.path.exists(old_file):
            raise FileNotFoundError(f"Old domains file not found: {old_file}")
        if not os.path.exists(new_file):
            raise FileNotFoundError(f"New domains file not found: {new_file}")
        
        return self.create_comparison(protein_id, old_file, new_file)

# Convenience functions
def create_comparison(protein_id: str, old_domains_file: str, new_domains_file: str,
                     output_dir: str = "/tmp/pymol_comparison") -> str:
    """Create PyMOL comparison with default settings"""
    visualizer = PyMOLVisualizer()
    return visualizer.create_comparison(protein_id, old_domains_file, 
                                      new_domains_file, output_dir)

def quick_comparison(protein_id: str, 
                    batch_dir: str = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424") -> str:
    """Quick comparison using standard file locations"""
    visualizer = PyMOLVisualizer()
    return visualizer.quick_comparison(protein_id, batch_dir)

if __name__ == "__main__":
    # Command line interface
    import argparse
    
    parser = argparse.ArgumentParser(description='Create PyMOL domain comparison')
    parser.add_argument('protein_id', help='Protein ID (e.g., 8ovp_A)')
    parser.add_argument('--old-file', help='Old domains XML file')
    parser.add_argument('--new-file', help='New domains XML file')
    parser.add_argument('--batch-dir', 
                        default='/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424',
                        help='Batch directory (for auto-detection)')
    parser.add_argument('--output-dir', default='/tmp/pymol_comparison',
                        help='Output directory')
    
    args = parser.parse_args()
    
    try:
        if args.old_file and args.new_file:
            script_path = create_comparison(args.protein_id, args.old_file, 
                                          args.new_file, args.output_dir)
        else:
            script_path = quick_comparison(args.protein_id, args.batch_dir)
        
        print(f"✅ PyMOL script created: {script_path}")
        print(f"Run: pymol {script_path}")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        exit(1)
