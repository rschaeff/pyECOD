#!/bin/bash
# Complete workflow for PyMOL domain comparison

# Configuration
PROTEIN_ID="${1:-8ovp_A}"
BATCH_DIR="/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424"
REFERENCE_FILE="test_data_range_cache/domain_lengths.csv"

echo "============================================================"
echo "Complete PyMOL Domain Comparison Workflow"
echo "============================================================"
echo "Protein: $PROTEIN_ID"
echo "Batch dir: $BATCH_DIR"
echo ""

# Step 1: Check if reference files exist
echo "Step 1: Checking reference files..."
if [ ! -f "$REFERENCE_FILE" ]; then
    echo "❌ Reference file missing: $REFERENCE_FILE"
    echo ""
    echo "Generating reference files from range cache..."
    
    # Check if range cache exists
    CACHE_FILE="/data/ecod/database_versions/v291/ecod.develop291.range_cache.txt"
    if [ ! -f "$CACHE_FILE" ]; then
        echo "❌ Range cache file missing: $CACHE_FILE"
        echo "Please ensure the range cache file is available"
        exit 1
    fi
    
    # Generate reference files
    python range_cache_parser.py --cache-file "$CACHE_FILE" --output-dir test_data_range_cache
    
    if [ $? -ne 0 ]; then
        echo "❌ Failed to generate reference files"
        exit 1
    fi
    
    echo "✅ Reference files generated"
else
    echo "✅ Reference files exist"
fi

echo ""

# Step 2: Run domain partitioning if output doesn't exist
echo "Step 2: Checking new domain partitioning results..."
NEW_DOMAINS="/tmp/${PROTEIN_ID}_mini.domains.xml"

if [ ! -f "$NEW_DOMAINS" ]; then
    echo "❌ New domains file missing: $NEW_DOMAINS"
    echo ""
    echo "Running mini_pyecod domain partitioning..."
    
    python scripts/quick_test.py "$PROTEIN_ID" \
        --batch-dir "$BATCH_DIR" \
        --reference-lengths "$REFERENCE_FILE"
    
    if [ $? -ne 0 ]; then
        echo "❌ Domain partitioning failed"
        exit 1
    fi
    
    echo "✅ Domain partitioning completed"
else
    echo "✅ New domains file exists"
fi

echo ""

# Step 3: Check old domains exist
echo "Step 3: Checking old domain partitioning results..."
OLD_DOMAINS="${BATCH_DIR}/domains/${PROTEIN_ID}.develop291.domains.xml"

if [ ! -f "$OLD_DOMAINS" ]; then
    echo "❌ Old domains file missing: $OLD_DOMAINS"
    echo "Cannot proceed without original algorithm results"
    exit 1
else
    echo "✅ Old domains file exists"
fi

echo ""

# Step 4: Generate PyMOL comparison
echo "Step 4: Generating PyMOL comparison..."
python scripts/run_pymol_comparison.py "$PROTEIN_ID" \
    --batch-dir "$BATCH_DIR" \
    --output-dir "/tmp"

if [ $? -ne 0 ]; then
    echo "❌ PyMOL comparison generation failed"
    exit 1
fi

echo ""
echo "============================================================"
echo "✅ WORKFLOW COMPLETED SUCCESSFULLY"
echo "============================================================"
echo ""
echo "To view the comparison:"
echo "  pymol /tmp/pymol_comparison/${PROTEIN_ID}_comparison.pml"
echo ""
echo "The session will be saved to:"
echo "  /tmp/pymol_comparison/${PROTEIN_ID}_comparison.pse"
echo ""
echo "Use 'Ctrl+S' in PyMOL to save the session if you make changes."
