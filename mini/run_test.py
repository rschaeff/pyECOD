# mini/run_test.py
def test_8ovp():
    """Test with residue blocking"""
    xml_path = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424/domains/8ovp_A.develop291.domain_summary.xml"

    evidence = parse_domain_summary(xml_path)
    print(f"Found {len(evidence)} evidence items")

    # Filter to high-quality evidence first
    good_evidence = [e for e in evidence if e.confidence > 0.6 or (e.evalue and e.evalue < 1e-5)]
    print(f"Filtered to {len(good_evidence)} high-quality evidence items")

    # Partition with blocking
    domains = partition_domains(good_evidence, sequence_length=518)

    print(f"\nFound {len(domains)} domains:")
    for domain in domains:
        print(f"  Domain {domain.id}:")
        print(f"    Family: {domain.family}")
        print(f"    Range: {domain.range}")
        print(f"    Size: {domain.range.size} residues")

    # Should be close to 3!
    if len(domains) == 3:
        print("\n✅ SUCCESS: Found expected 3 domains!")
    else:
        print(f"\n⚠️  Found {len(domains)} domains (expected 3)")
