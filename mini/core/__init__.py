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
# 1. UPDATE mini/core/__init__.py - Add enhanced writer functions
__all__ = [
    'partition_domains',
    'parse_domain_summary',
    'load_domain_definitions',
    'load_reference_lengths',
    'load_protein_lengths',
    'write_domain_partition',
    'write_domain_partition_from_layout',  # NEW: Enhanced writer
    'create_metadata_from_batch',          # NEW: Metadata creation
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
from .writer import (
    write_domain_partition,
    write_domain_partition_from_layout,
    create_metadata_from_batch
)
from .evidence_utils import (
    calculate_evidence_confidence,
    populate_evidence_provenance,
    validate_evidence_provenance
)
from .domain_utils import (
    create_domain_with_provenance,
    validate_domain_provenance,
    get_domain_coverage_stats
)
from .blast_parser import load_chain_blast_alignments
