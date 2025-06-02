# mini/__init__.py
"""
Mini pyECOD - Minimal domain partitioning

RULES:
1. NO imports from pipelines/
2. NO database access
3. NO external dependencies beyond stdlib + lxml
4. ONLY import from ecod.core (shared utilities)
5. Must work standalone
"""

# Enforce isolation
__all__ = [
    'partition_domains',
    'parse_domain_summary',
    'load_domain_definitions',
    'load_reference_lengths',
    'load_protein_lengths',
    'write_domain_partition',
    'load_chain_blast_alignments'
]

# Guard against bad imports
import sys
if 'ecod.pipelines' in sys.modules:
    raise ImportError("mini_pyecod must not import from pipelines!")

# Import main functions for convenience
from .parser import parse_domain_summary, load_reference_lengths, load_protein_lengths
from .partitioner import partition_domains
from .decomposer import load_domain_definitions
from .writer import write_domain_partition
from .blast_parser import load_chain_blast_alignments
