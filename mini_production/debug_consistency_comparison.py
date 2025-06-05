#!/usr/bin/env python3
"""
Fixed consistency comparison debug script
"""

import sys
import yaml
import psycopg2
import psycopg2.extras
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict, Counter
import numpy as np
from dataclasses import dataclass
import logging

# Import from the fixed consistency comparison
from consistency_comparison import ConsistencyComparator

def debug_consistency_comparison_fixed(protein_id: str = "8ovp_A"):
    """Debug consistency-based comparison for specific protein with FIXED implementation"""

    print(f"🔍 FIXED CONSISTENCY-BASED Analysis for {protein_id}")
    print("=" * 70)

    comparator = ConsistencyComparator()

    # Find files
    batch_base = Path("../data/ecod/pdb_updates/batches")
    mini_file = None
    main_file = None

    for batch_dir in batch_base.iterdir():
        if not batch_dir.is_dir():
            continue

        # Look for mini file
        mini_domains_dir = batch_dir / "mini_domains"
        if mini_domains_dir.exists():
            test_file = mini_domains_dir / f"{protein_id}.mini.domains.xml"
            if test_file.exists():
                mini_file = test_file

        # Look for main file
        domains_dir = batch_dir / "domains"
        if domains_dir.exists():
            test_file = domains_dir / f"{protein_id}.develop291.domains.xml"
            if test_file.exists():
                main_file = test_file

    if not mini_file:
        print("❌ No mini file found!")
        return

    # Parse domains with FIXED parser
    mini_domains = comparator.parse_partition_domains(mini_file)
    main_domains = comparator.parse_partition_domains(main_file) if main_file and main_file.exists() else []

    print(f"📊 DOMAIN ANNOTATIONS (FIXED PARSING):")
    print(f"  Mini: {len(mini_domains)} domains")
    for i, domain in enumerate(mini_domains):
        print(f"    D{i+1}: {domain.t_group} | {domain.domain_range} | length={domain.length}")

    print(f"  Main: {len(main_domains)} domains")
    for i, domain in enumerate(main_domains):
        print(f"    D{i+1}: {domain.t_group} | {domain.domain_range} | length={domain.length}")

    # Get family patterns
    all_families = list(set(d.t_group for d in mini_domains + main_domains if d.t_group))
    family_patterns = comparator.get_family_patterns(all_families)

    print(f"\n📚 FAMILY PATTERNS (FIXED):")
    for family_id, pattern in family_patterns.items():
        print(f"  {family_id}:")
        print(f"    Typical length: {pattern.typical_length_range[0]}-{pattern.typical_length_range[1]}")
        print(f"    Split frequency: {pattern.split_frequency:.3f}")
        print(f"    Examples: {', '.join(pattern.representative_examples[:3])}")

    # Calculate consistency scores with FIXED implementation
    mini_scores = comparator.calculate_overall_consistency_score(mini_domains)
    main_scores = comparator.calculate_overall_consistency_score(main_domains)

    print(f"\n🎯 FIXED CONSISTENCY SCORES:")
    print(f"  Mini Algorithm:")
    print(f"    Overall: {mini_scores['overall']:.3f}")
    print(f"    Family consistency: {mini_scores['family_consistency']:.3f}")
    print(f"    Boundary quality: {mini_scores['boundary_quality']:.3f}")
    print(f"    Coverage: {mini_scores['coverage']:.3f}")

    print(f"  Main Algorithm:")
    print(f"    Overall: {main_scores['overall']:.3f}")
    print(f"    Family consistency: {main_scores['family_consistency']:.3f}")
    print(f"    Boundary quality: {main_scores['boundary_quality']:.3f}")
    print(f"    Coverage: {main_scores['coverage']:.3f}")

    # Winner
    score_diff = mini_scores['overall'] - main_scores['overall']
    if score_diff > 0.05:
        winner = "MINI"
    elif score_diff < -0.05:
        winner = "MAIN"
    else:
        winner = "TIE"

    print(f"\n🏆 CONSISTENCY WINNER: {winner}")
    print(f"  Score difference: {score_diff:+.3f}")

    # Analysis
    print(f"\n💡 CONSISTENCY ANALYSIS:")

    for domain in mini_domains:
        if domain.t_group in family_patterns:
            pattern = family_patterns[domain.t_group]
            in_range = (domain.length > 0 and
                       pattern.typical_length_range[0] <= domain.length <= pattern.typical_length_range[1])
            print(f"  Mini {domain.t_group}: length {domain.length} "
                  f"{'✅' if in_range else '⚠️'} (typical: {pattern.typical_length_range[0]}-{pattern.typical_length_range[1]})")

    for domain in main_domains:
        if domain.t_group in family_patterns:
            pattern = family_patterns[domain.t_group]
            in_range = (domain.length > 0 and
                       pattern.typical_length_range[0] <= domain.length <= pattern.typical_length_range[1])
            print(f"  Main {domain.t_group}: length {domain.length} "
                  f"{'✅' if in_range else '⚠️'} (typical: {pattern.typical_length_range[0]}-{pattern.typical_length_range[1]})")

    # Biological interpretation
    print(f"\n🧬 BIOLOGICAL INTERPRETATION:")

    # Count family occurrences
    mini_families = Counter(d.t_group for d in mini_domains)
    main_families = Counter(d.t_group for d in main_domains)

    print(f"  Family distribution:")
    print(f"    Mini: {dict(mini_families)}")
    print(f"    Main: {dict(main_families)}")

    # Specific analysis for known families
    if "7523.1.1" in family_patterns:  # GFP-like
        gfp_pattern = family_patterns["7523.1.1"]
        mini_gfp = mini_families.get("7523.1.1", 0)
        main_gfp = main_families.get("7523.1.1", 0)

        print(f"\n  GFP family (7523.1.1) analysis:")
        print(f"    Split frequency in database: {gfp_pattern.split_frequency:.3f}")
        print(f"    Mini splits GFP into: {mini_gfp} pieces")
        print(f"    Main splits GFP into: {main_gfp} pieces")

        if gfp_pattern.split_frequency < 0.1:  # Rarely split
            if mini_gfp == 1:
                print(f"    ✅ Mini preserves GFP integrity (consistent with database)")
            else:
                print(f"    ⚠️ Mini splits GFP (inconsistent with database)")

            if main_gfp == 1:
                print(f"    ✅ Main preserves GFP integrity (consistent with database)")
            else:
                print(f"    ⚠️ Main splits GFP (inconsistent with database)")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Fixed consistency-based domain comparison')
    parser.add_argument('--protein', default='8ovp_A', help='Protein to analyze')

    args = parser.parse_args()

    debug_consistency_comparison_fixed(args.protein)
