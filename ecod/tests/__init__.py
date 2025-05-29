[tool:pytest]
# Main pytest configuration for ECOD test suite

# Test discovery
testpaths = .
python_files = test_*.py *_test.py
python_classes = Test* *Tests
python_functions = test_*

# Output and reporting
addopts =
    -v
    --tb=short
    --strict-markers
    --strict-config
    --disable-warnings
    --color=yes

# Custom markers
markers =
    unit: Unit tests
    integration: Integration tests
    slow: Slow running tests (may be skipped)
    performance: Performance tests
    golden: Golden dataset tests
    evidence: Evidence processing tests
    weights: Evidence weight tests
    service: Service integration tests
    regression: Regression detection tests

# Filtering
filterwarnings =
    ignore::DeprecationWarning
    ignore::PendingDeprecationWarning
    ignore::UserWarning:psycopg2

# Minimum version requirement
minversion = 6.0

# Test timeout (in seconds)
timeout = 300

# Parallel execution settings (if using pytest-xdist)
# addopts = -n auto
