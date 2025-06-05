#!/usr/bin/env python3
"""
Fixed PyMOL domain comparison visualization with proper coordinate translation

Critical fix: Translates seqid ranges to PDB residue numbers using pdbx_poly_seq_scheme
"""

import os
import gzip
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List, Tuple, Optional, Dict
from dataclasses import dataclass
import re

@dataclass
class Domain:
    """Domain representation for visualization"""
    id: str
    range_str: str
    family: str
    segments: List[Tuple[int, int]]  # seqid coordinates
    pdb_segments: List[Tuple[str, str]] = None  # PDB coordinates

@dataclass
class ResidueMapping:
    """Mapping between seqid and PDB coordinates"""
    seqid: int
    pdb_num: str  # Can include insertion codes like "27A"
    pdb_ins_code: str = ""

class CoordinateTranslator:
    """Handles seqid to PDB coordinate translation using mmCIF pdbx_poly_seq_scheme"""

    def __init__(self, structure_path: str, chain_id: str):
        self.structure_path = structure_path
        self.chain_id = chain_id
        self.seqid_to_pdb = {}  # seqid -> PDB residue number
        self.pdb_to_seqid = {}  # PDB residue number -> seqid
        self._parse_coordinate_mapping()

    def _parse_coordinate_mapping(self):
        """Parse pdbx_poly_seq_scheme from mmCIF file"""

        if self.structure_path.endswith('.gz'):
            opener = gzip.open
            mode = 'rt'
        else:
            opener = open
            mode = 'r'

        try:
            with opener(self.structure_path, mode) as f:
                lines = f.readlines()

                in_poly_seq_scheme = False
                header_processed = False

                for i, line in enumerate(lines):
                    line = line.strip()

                    # Look for the start of pdbx_poly_seq_scheme section
                    if line.startswith('loop_'):
                        # Check if next few lines contain pdbx_poly_seq_scheme
                        for j in range(i+1, min(i+20, len(lines))):
                            if 'pdbx_poly_seq_scheme' in lines[j]:
                                in_poly_seq_scheme = True
                                break
                        continue

                    # Skip header lines in the poly_seq_scheme section
                    if in_poly_seq_scheme and line.startswith('_pdbx_poly_seq_scheme.'):
                        header_processed = True
                        continue

                    # Parse data lines
                    if in_poly_seq_scheme and header_processed and line:
                        # End of section?
                        if line.startswith('loop_') or line.startswith('#') or line.startswith('_'):
                            break

                        # Parse the data line
                        parts = line.split()
                        if len(parts) >= 6:
                            try:
                                asym_id = parts[0]
                                entity_id = parts[1]
                                seq_id = int(parts[2])  # seqid position
                                pdb_strand_id = parts[3]  # chain ID
                                auth_seq_num = parts[4]  # PDB residue number
                                pdb_ins_code = parts[5] if len(parts) > 5 else ""

                                # Only process our chain
                                if pdb_strand_id == self.chain_id:
                                    # Handle insertion codes
                                    if pdb_ins_code and pdb_ins_code != '.' and pdb_ins_code != '?':
                                        pdb_res_id = f"{auth_seq_num}{pdb_ins_code}"
                                    else:
                                        pdb_res_id = auth_seq_num

                                    # Store mapping
                                    self.seqid_to_pdb[seq_id] = pdb_res_id
                                    self.pdb_to_seqid[pdb_res_id] = seq_id

                            except (ValueError, IndexError):
                                continue

        except Exception as e:
            print(f"Warning: Could not parse coordinate mapping from {self.structure_path}: {e}")
            # Fall back to 1:1 mapping
            self._create_fallback_mapping()

    def _create_fallback_mapping(self):
        """Create 1:1 fallback mapping when mmCIF parsing fails"""
        print(f"Using fallback 1:1 coordinate mapping for chain {self.chain_id}")
        # This assumes seqid == PDB residue number (not always true!)
        for i in range(1, 1000):  # Reasonable upper bound
            self.seqid_to_pdb[i] = str(i)
            self.pdb_to_seqid[str(i)] = i

    def translate_seqid_range(self, seqid_start: int, seqid_end: int) -> Tuple[str, str]:
        """Translate seqid range to PDB residue range"""

        pdb_start = self.seqid_to_pdb.get(seqid_start)
        pdb_end = self.seqid_to_pdb.get(seqid_end)

        if pdb_start is None or pdb_end is None:
            print(f"Warning: Could not translate seqid range {seqid_start}-{seqid_end} "
                  f"to PDB coordinates for chain {self.chain_id}")
            # Fallback to original numbers
            return str(seqid_start), str(seqid_end)

        return pdb_start, pdb_end

    def get_mapping_info(self) -> Dict:
        """Get mapping statistics for debugging"""
        return {
            'total_residues_mapped': len(self.seqid_to_pdb),
            'seqid_range': f"{min(self.seqid_to_pdb.keys())}-{max(self.seqid_to_pdb.keys())}" if self.seqid_to_pdb else "None",
            'pdb_range': f"{min(self.pdb_to_seqid.keys(), key=lambda x: int(re.sub(r'[A-Z]', '', x)))}-{max(self.pdb_to_seqid.keys(), key=lambda x: int(re.sub(r'[A-Z]', '', x)))}" if self.pdb_to_seqid else "None",
            'has_insertion_codes': any(not x.isdigit() for x in self.pdb_to_seqid.keys())
        }

class PyMOLVisualizer:
    """PyMOL visualization with proper coordinate translation"""

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
        Now with proper coordinate translation!
        """
        # Parse domains from both files
        old_domains = self._parse_domain_file(old_domains_file, is_old_format=True)
        new_domains = self._parse_domain_file(new_domains_file, is_old_format=False)

        # Find structure file
        pdb_id = protein_id.split('_')[0]
        chain_id = protein_id.split('_')[1] if '_' in protein_id else 'A'
        structure_path = self._find_structure(pdb_id)
        if not structure_path:
            raise FileNotFoundError(f"Structure not found for {pdb_id}")

        # Create coordinate translator
        translator = CoordinateTranslator(structure_path, chain_id)
        mapping_info = translator.get_mapping_info()
        print(f"Coordinate mapping for {protein_id}: {mapping_info}")

        # Translate domain coordinates
        old_domains = self._translate_domain_coordinates(old_domains, translator)
        new_domains = self._translate_domain_coordinates(new_domains, translator)

        # Create output directory
        os.makedirs(output_dir, exist_ok=True)

        # Generate PyMOL script
        script_path = os.path.join(output_dir, f"{protein_id}_comparison.pml")
        self._write_pymol_script(protein_id, structure_path, old_domains,
                                new_domains, script_path, output_dir, mapping_info)

        return script_path

    def _translate_domain_coordinates(self, domains: List[Domain],
                                    translator: CoordinateTranslator) -> List[Domain]:
        """Translate all domain coordinates from seqid to PDB"""

        for domain in domains:
            pdb_segments = []
            for seqid_start, seqid_end in domain.segments:
                pdb_start, pdb_end = translator.translate_seqid_range(seqid_start, seqid_end)
                pdb_segments.append((pdb_start, pdb_end))
            domain.pdb_segments = pdb_segments

        return domains

    def _parse_domain_file(self, xml_path: str, is_old_format: bool = False) -> List[Domain]:
        """Parse domain XML file (unchanged from original)"""
        domains = []

        try:
            tree = ET.parse(xml_path)
            root = tree.getroot()

            for i, domain_elem in enumerate(root.findall(".//domain")):
                if is_old_format:
                    range_str = domain_elem.get("range", "")
                    family = domain_elem.get("source_id", "") or domain_elem.get("t_group", "unknown")
                    domain_id = f"old_d{i+1}"
                else:
                    domain_id = domain_elem.get("id", f"d{i+1}")
                    range_str = domain_elem.get("range", "")
                    family = domain_elem.get("family", "unknown")

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
        """Parse range string into list of (start, end) tuples (unchanged)"""
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
                    pos = int(segment)
                    segments.append((pos, pos))
        except ValueError as e:
            print(f"Warning: Could not parse range '{range_str}': {e}")

        return segments

    def _find_structure(self, pdb_id: str) -> Optional[str]:
        """Find structure file in PDB repository (unchanged)"""
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
                           script_path: str, output_dir: str, mapping_info: Dict):
        """Write PyMOL comparison script with proper PDB coordinates"""

        pdb_id = protein_id.split('_')[0]
        chain_id = protein_id.split('_')[1] if '_' in protein_id else 'A'

        lines = [
            f"# PyMOL domain comparison for {protein_id}",
            f"# Old algorithm: {len(old_domains)} domains",
            f"# New algorithm: {len(new_domains)} domains",
            f"# Coordinate mapping: {mapping_info['total_residues_mapped']} residues",
            f"# SeqID range: {mapping_info['seqid_range']}, PDB range: {mapping_info['pdb_range']}",
            f"# Has insertion codes: {mapping_info['has_insertion_codes']}",
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
            "# Color old algorithm domains (using PDB coordinates)",
            f"# {len(old_domains)} domains from previous algorithm"
        ]

        # Color old domains using PDB coordinates
        for i, domain in enumerate(old_domains):
            color = self.colors[i % len(self.colors)]
            lines.append(f"# Domain {i+1}: {domain.family} (seqid: {domain.range_str})")

            if domain.pdb_segments:
                for pdb_start, pdb_end in domain.pdb_segments:
                    # Handle both numeric and alphanumeric residue numbers
                    if pdb_start == pdb_end:
                        selection = f"old_algorithm and chain {chain_id} and resi {pdb_start}"
                    else:
                        selection = f"old_algorithm and chain {chain_id} and resi {pdb_start}-{pdb_end}"
                    lines.append(f"color {color}, {selection}")
                    lines.append(f"# PDB coordinates: {pdb_start}-{pdb_end}")
            else:
                lines.append(f"# Warning: No PDB coordinates available for this domain")

        lines.extend([
            "",
            "# Color new algorithm domains (using PDB coordinates)",
            f"# {len(new_domains)} domains from mini_pyecod"
        ])

        # Color new domains using PDB coordinates
        for i, domain in enumerate(new_domains):
            color = self.colors[i % len(self.colors)]
            lines.append(f"# Domain {i+1}: {domain.family} (seqid: {domain.range_str})")

            if domain.pdb_segments:
                for pdb_start, pdb_end in domain.pdb_segments:
                    if pdb_start == pdb_end:
                        selection = f"new_algorithm and chain {chain_id} and resi {pdb_start}"
                    else:
                        selection = f"new_algorithm and chain {chain_id} and resi {pdb_start}-{pdb_end}"
                    lines.append(f"color {color}, {selection}")
                    lines.append(f"# PDB coordinates: {pdb_start}-{pdb_end}")
            else:
                lines.append(f"# Warning: No PDB coordinates available for this domain")

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
            f"save {output_dir}/{protein_id}_comparison.pse",
            "",
            "# Debug info - show coordinate mapping quality",
            f"print 'Coordinate translation summary:'",
            f"print 'Total residues mapped: {mapping_info['total_residues_mapped']}'",
            f"print 'SeqID range: {mapping_info['seqid_range']}'",
            f"print 'PDB range: {mapping_info['pdb_range']}'",
            f"print 'Has insertion codes: {mapping_info['has_insertion_codes']}'"
        ])

        # Write script
        with open(script_path, 'w') as f:
            f.write('\n'.join(lines))

    def quick_comparison(self, protein_id: str, batch_dir: str) -> str:
        """Quick comparison using standard file locations (unchanged interface)"""
        old_file = f"{batch_dir}/domains/{protein_id}.develop291.domains.xml"
        new_file = f"/tmp/{protein_id}_mini.domains.xml"

        if not os.path.exists(old_file):
            raise FileNotFoundError(f"Old domains file not found: {old_file}")
        if not os.path.exists(new_file):
            raise FileNotFoundError(f"New domains file not found: {new_file}")

        return self.create_comparison(protein_id, old_file, new_file)

# Convenience functions (unchanged interface)
def create_comparison(protein_id: str, old_domains_file: str, new_domains_file: str,
                     output_dir: str = "/tmp/pymol_comparison") -> str:
    """Create PyMOL comparison with coordinate translation"""
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

    parser = argparse.ArgumentParser(description='Create PyMOL domain comparison with coordinate translation')
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
