"""
HHSearch-related commands for the ECOD pipeline 
"""

import argparse
import logging
from typing import Dict, Any

from ecod.config import ConfigManager
from ecod.pipelines.hhsearch_pipeline import HHSearchPipeline
from ecod.db.manager import DBManager
