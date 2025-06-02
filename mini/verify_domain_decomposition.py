#!/usr/bin/env python3
"""
Verify that domain decomposition setup is correct for 2ia4
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.decomposer import load_domain_definitions
from mini.parser import load_protein_lengths
from mini.blast_parser import load_chain_blast_alignments

def verify_setup():
    """Verify that decomposition setup is correct"""
    
    print("Verifying Domain Decomposition Setup")
    print("=" * 40)
    
    # Check domain definitions
    domain_def_file = Path(__file__).parent / "test_data" / "domain_definitions.csv"
    if not domain_def_file.exists():
        print(f"❌ Domain definitions file not found: {domain_def_file}")
        return False
    
    print(f"✓ Found domain definitions: {domain_def_file}")
    
    # Load domain definitions
    domain_definitions = load_domain_definitions(str(domain_def_file), verbose=True)
    
    if not domain_definitions:
        print("❌ No domain definitions loaded")
        return False
    
    print(f"✓ Loaded domain definitions for {len(domain_definitions)} protein chains")
    
    # Check specifically for 2ia4
    target_keys = [
        ('2ia4', 'A'),   # Standard format
        ('2ia4', 'a'),   # Lowercase chain
    ]
    
    found_2ia4 = False
    for key in target_keys:
        if key in domain_definitions:
            found_2ia4 = True
            domains = domain_definitions[key]
            print(f"\n✓ Found 2ia4 domains ({len(domains)} domains):")
            for domain in domains:
                print(f"  - {domain.domain_id}: {domain.range} (length: {domain.length})")
            break
    
    if not found_2ia4:
        print("\n❌ 2ia4 not found in domain definitions")
        print("Available proteins:")
        for (pdb, chain), domains in list(domain_definitions.items())[:10]:
            print(f"  {pdb}_{chain}: {len(domains)} domains")
        return False
    
    # Check protein lengths
    protein_lengths_file = Path(__file__).parent / "test_data" / "protein_lengths.csv"
    if protein_lengths_file.exists():
        protein_lengths = load_protein_lengths(str(protein_lengths_file))
        print(f"\n✓ Loaded protein lengths for {len(protein_lengths)} proteins")
        
        # Check for 2ia4
        found_2ia4_length = False
        for key in [('2ia4', 'A'), ('2ia4', 'a')]:
            if key in protein_lengths:
                print(f"✓ 2ia4_A protein length: {protein_lengths[key]} residues")
                found_2ia4_length = True
                break
        
        if not found_2ia4_length:
            print("⚠️  2ia4 protein length not found (may affect chain BLAST confidence)")
    
    # Check BLAST data availability
    batch_dir = Path("/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424")
    blast_dir = batch_dir / "blast" / "chain"
    blast_file = blast_dir / "8ovp_A.develop291.xml"
    
    if blast_file.exists():
        print(f"\n✓ Found BLAST XML: {blast_file}")
        
        # Try to load alignments
        alignments = load_chain_blast_alignments(str(blast_dir), "8ovp", "A", verbose=True)
        print(f"✓ Loaded {len(alignments)} BLAST alignments")
        
        # Check for 2ia4 alignment
        if ('2ia4', 'A') in alignments or ('2ia4', 'a') in alignments:
            print("✓ Found 2ia4 alignment - decomposition should work!")
        else:
            print("⚠️  No 2ia4 alignment found - check BLAST results")
            print("Available alignments:")
            for (pdb, chain) in list(alignments.keys())[:5]:
                print(f"  {pdb}_{chain}")
    else:
        print(f"\n❌ BLAST XML not found: {blast_file}")
        return False
    
    print(f"\n{'='*40}")
    print("Decomposition setup verification complete!")
    print("If all checks passed, 2ia4 should decompose properly.")
    
    return True

def show_expected_decomposition():
    """Show what we expect to happen with 2ia4 decomposition"""
    
    print(f"\n{'='*40}")
    print("Expected Decomposition Process")
    print(f"{'='*40}")
    
    print("Current situation:")
    print("  Input: 2ia4 @ 2-248,491-517 (discontinuous, 274 residues)")
    print("  Problem: This is treated as ONE domain")
    
    print("\nWith proper decomposition:")
    print("  1. Load 2ia4 domain definitions from CSV")
    print("  2. Load BLAST alignment query<->hit mapping")
    print("  3. Map query regions to 2ia4 reference domains")
    print("  4. Replace single domain with multiple domains")
    
    print("\nExpected output:")
    print("  Domain 1: e2ia4A1 @ 2-248 (247 residues)")
    print("  Domain 2: e2ia4A2 @ 491-517 (27 residues)")
    print("  Plus any other domains from other evidence")

if __name__ == "__main__":
    if verify_setup():
        show_expected_decomposition()
        print(f"\n✓ Ready to test with: ./pyecod_mini 8ovp_A --verbose")
    else:
        print(f"\n❌ Setup issues found - fix before testing")
