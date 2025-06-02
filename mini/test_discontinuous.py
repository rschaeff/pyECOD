#!/usr/bin/env python3
"""Test discontinuous chain BLAST decomposition"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.models import Evidence
from mini.decomposer import decompose_chain_blast_discontinuous
from ecod.core.sequence_range import SequenceRange

def test_discontinuous_decomposition():
    """Test that discontinuous hits are properly decomposed"""
    
    # Create a mock discontinuous chain BLAST hit like 2ia4
    evidence = Evidence(
        type="chain_blast",
        source_pdb="2ia4",
        query_range=SequenceRange.parse("2-248,491-517"),
        confidence=0.95,
        evalue=1e-100,
        domain_id="2ia4_A",
        alignment_coverage=0.95
    )
    
    print(f"Original evidence: {evidence.source_pdb} @ {evidence.query_range}")
    print(f"Is discontinuous: {evidence.query_range.is_discontinuous}")
    
    # Decompose
    decomposed = decompose_chain_blast_discontinuous(evidence)
    
    print(f"\nDecomposed into {len(decomposed)} pieces:")
    for i, dec in enumerate(decomposed):
        print(f"  {i+1}. {dec.domain_id}: {dec.query_range} ({dec.query_range.total_length} residues)")

if __name__ == "__main__":
    test_discontinuous_decomposition()
