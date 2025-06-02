# mini/partitioner.py
"""Ultra-simple partitioning - each unique evidence location is a domain"""

from collections import defaultdict
from typing import List, Dict, Tuple
from .models import Evidence, Domain

def partition_domains(evidence_list: List[Evidence], sequence_length: int) -> List[Domain]:
    """
    Simplest correct algorithm:
    - Group evidence by (source_family, range)
    - Each unique combination becomes a domain
    - NO MERGING
    """
    # Group evidence by exact family AND position
    evidence_groups: Dict[Tuple[str, str], List[Evidence]] = defaultdict(list)

    for evidence in evidence_list:
        # Determine family identifier
        if evidence.t_group:
            family = evidence.t_group  # Prefer classification
        elif evidence.source_pdb:
            family = evidence.source_pdb
        else:
            family = "unknown"

        # Create unique key
        range_key = str(evidence.query_range)
        key = (family, range_key)
        evidence_groups[key].append(evidence)

    # Convert to domains
    domains = []
    domain_num = 1

    for (family, range_str), group in sorted(evidence_groups.items()):
        # Pick best evidence from group
        best_evidence = max(group, key=lambda e: (e.confidence, -e.evalue if e.evalue else 0))

        domain = Domain(
            id=f"d{domain_num}",
            range=best_evidence.query_range,
            family=family,
            evidence_count=len(group),
            source=best_evidence.type,
            evidence_items=group
        )
        domains.append(domain)
        domain_num += 1

    return domains
