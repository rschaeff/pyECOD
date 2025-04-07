#!/usr/bin/env python3
"""
ECOD (Evolutionary Classification Of protein Domains) Pipeline

Command-line interface for protein domain classification and analysis.
"""

import sys
import os

# Add the parent directory to the path so we can import the ecod package
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

# Import the main entry point
from ecod.cli.main import main

if __name__ == "__main__":
    sys.exit(main())