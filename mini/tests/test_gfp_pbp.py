# mini_pyecod/test_gfp_pbp.py
"""Test the GFP+PBP case"""

from parser import parse_domain_summary
from partitioner import partition_domains
from writer import write_domain_partition

def test_gfp_pbp():
    # Parse evidence
    evidence = parse_domain_summary("test_data/8ovp_A.domain_summary.xml")
    
    print(f"Found {len(evidence)} evidence items")
    
    # Show what we found
    print("\nKey evidence:")
    for ev in evidence[:5]:  # First 5
        print(f"  {ev.type}: {ev.source_pdb} @ {ev.query_range}")
    
    # Partition (NO MERGING across families!)
    domains = partition_domains(evidence, sequence_length=518)
    
    print(f"\nFound {len(domains)} domains:")
    for domain in domains:
        print(f"  Domain {domain.id}:")
        print(f"    Family: {domain.family}")  
        print(f"    Range: {domain.range}")
        print(f"    Discontinuous: {domain.range.is_discontinuous}")
        print(f"    Evidence: {domain.evidence_count} hits")
    
    # Write output
    write_domain_partition(domains, "8ovp", "A", "output/8ovp_A.domains.xml")
    
    # Verify we got separate GFP and PBP domains
    gfp_domains = [d for d in domains if d.family == "6dgv"]
    pbp_domains = [d for d in domains if d.family == "2ia4"]
    
    assert len(gfp_domains) >= 1, "Should find GFP domain"
    assert len(pbp_domains) >= 1, "Should find PBP domain"
    
    print("\n✅ TEST PASSED: Found separate GFP and PBP domains!")

if __name__ == "__main__":
    test_gfp_pbp()
