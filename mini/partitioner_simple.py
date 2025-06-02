# mini_pyecod/partitioner_simple.py
"""Ultra-simple: Each evidence hit = one domain. Period."""

def partition_domains_simple(evidence_list: List[Evidence]) -> List[Domain]:
    """
    Simplest possible algorithm:
    - Each evidence hit becomes a domain
    - Group identical hits together
    - NO MERGING EVER
    """
    # Group evidence by exact position and source
    evidence_groups = defaultdict(list)
    
    for evidence in evidence_list:
        key = (evidence.source_pdb, str(evidence.query_range))
        evidence_groups[key].append(evidence)
    
    # Convert to domains
    domains = []
    domain_num = 1
    
    for (source_pdb, range_str), group in evidence_groups.items():
        # All evidence in group has same range, just take first
        domain = Domain(
            id=f"d{domain_num}",
            range=group[0].query_range,
            family=source_pdb,
            evidence_count=len(group),
            source=group[0].type,
            evidence_items=group
        )
        domains.append(domain)
        domain_num += 1
    
    return domains
