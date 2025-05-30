#!/bin/bash
# scripts/run_simple_e2e_tests.sh

set -e

echo "Running simplified E2E integration tests..."

# Set up environment
export PYTHONPATH="$(pwd):$PYTHONPATH"

# Run the simplified tests first
echo "Running script-only tests..."
python -m pytest ecod/tests/integration/test_domain_partition_e2e_simple.py -v --tb=short

echo "Running existing integration tests (non-E2E)..."
python -m pytest ecod/tests/integration/test_simple_integration.py -v --tb=short

echo "Running domain processing debug tests..."
python -m pytest ecod/tests/integration/test_domain_processing_debug.py -v --tb=short

echo "Integration tests complete!"
