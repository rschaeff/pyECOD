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
__all__ = ['partition_domains', 'parse_domain_summary']

# Guard against bad imports
import sys
if 'ecod.pipelines' in sys.modules:
    raise ImportError("mini_pyecod must not import from pipelines!")
