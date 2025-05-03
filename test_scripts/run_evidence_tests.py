#!/usr/bin/env python3
# run_evidence_tests.py

import os
import sys
import unittest

# Add the project root to the path
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

# Import the test case
from ecod.tests.test_utils.test_evidence_bridge import TestEvidenceIntegration

if __name__ == "__main__":
    # Create a test suite with just this test case
    suite = unittest.TestLoader().loadTestsFromTestCase(TestEvidenceIntegration)
    
    # Run the tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Exit with non-zero code if tests failed
    sys.exit(not result.wasSuccessful())
