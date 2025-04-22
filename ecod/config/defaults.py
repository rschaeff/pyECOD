#!/usr/bin/env python3
"""
Default configuration values for the ECOD pipeline
"""

DEFAULT_CONFIG = {
    'database': {
        'database': 'ecod_protein',
        'host': 'dione',
        'port': 45000,
        'user': 'ecod',
    },
    'paths': {
        'output_dir': './output',
    },
    'tools': {
        'blast_path': 'blastp',
        'hhblits_path': 'hhblits',
        'hhsearch_path': 'hhsearch',
        'hhmake_path': 'hhmake',
    },
    'reference': {
        'current_version': 'develop291',
    },
    'logging': {
        'level': 'INFO',
        'format': '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    },
    'pipeline': {
        'batch_size': 100,
        'threads': 8,
        'memory': '8G',
        'timeout': 3600,
        'check_interval': 300,
        'force_overwrite': False
    }
}