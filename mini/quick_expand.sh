#!/bin/bash
"""
Quick Test Suite Expansion Script

This script automates the most common workflow for expanding your test suite.

Usage:
    ./quick_expand.sh                    # Full workflow  
    ./quick_expand.sh --setup-only       # Just setup validation
    ./quick_expand.sh --test-only        # Just run tests
    ./quick_expand.sh --curate-only      # Just generate curation templates
"""

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Output directory
OUTPUT_DIR="/tmp/test_suite_expansion"

# Functions
print_header() {
    echo -e "\n${BLUE}============================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}============================================${NC}\n"
}

print_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

print_error() {
    echo -e "${RED}❌ $1${NC}"
}

check_prerequisites() {
    print_header "Checking Prerequisites"
    
    # Check we're in the right directory
    if [ ! -f "pyecod_mini" ]; then
        print_error "pyecod_mini not found. Run this script from the mini/ directory."
        exit 1
    fi
    print_success "Found pyecod_mini executable"
    
    # Check Python scripts exist
    if [ ! -f "expand_test_suite.py" ]; then
        print_error "expand_test_suite.py not found"
        exit 1
    fi
    print_success "Found expand_test_suite.py"
    
    if [ ! -f "curation_validator.py" ]; then
        print_error "curation_validator.py not found"
        exit 1
    fi
    print_success "Found curation_validator.py"
    
    # Check batch_test_proteins.py
    if [ ! -f "batch_test_proteins.py" ]; then
        print_error "batch_test_proteins.py not found"
        exit 1
    fi
    
    PROTEIN_COUNT=$(python -c "from batch_test_proteins import BATCH_TEST_PROTEINS; print(len(BATCH_TEST_PROTEINS))")
    print_success "Found $PROTEIN_COUNT test protein candidates"
}

validate_setup() {
    print_header "Validating Setup"
    
    echo "Checking pyecod_mini configuration..."
    if python expand_test_suite.py --validate-setup; then
        print_success "Setup validation passed"
        return 0
    else
        print_error "Setup validation failed"
        echo
        echo "Common fixes:"
        echo "  - Run: python pyecod_mini --setup-references"
        echo "  - Check batch directory permissions"
        echo "  - Verify reference data files exist"
        return 1
    fi
}

run_tests() {
    print_header "Running Algorithm Tests"
    
    echo "Testing all candidate proteins..."
    echo "This may take 5-10 minutes depending on protein count and batch size."
    echo
    
    # Start with a smaller subset to validate quickly
    echo "Phase 1: Testing simple cases first..."
    if python expand_test_suite.py --run-all --categories validated single_domain_small single_domain_medium --verbose; then
        print_success "Phase 1 tests completed"
    else
        print_warning "Some Phase 1 tests failed - continuing with remaining categories"
    fi
    
    echo
    echo "Phase 2: Testing remaining categories..."
    if python expand_test_suite.py --run-all --categories multi_domain_clear chain_blast_multi large_complex minimal_evidence --verbose; then
        print_success "Phase 2 tests completed"
    else
        print_warning "Some Phase 2 tests failed - check individual protein issues"
    fi
    
    echo
    echo "Getting test summary..."
    python expand_test_suite.py --summary
}

generate_curation() {
    print_header "Generating Curation Templates"
    
    echo "Creating curation templates for successful proteins..."
    if python expand_test_suite.py --curate; then
        print_success "Curation templates generated"
        
        CURATION_DIR="$OUTPUT_DIR/curation"
        if [ -d "$CURATION_DIR" ]; then
            CURATION_COUNT=$(find "$CURATION_DIR" -name "*_curation.md" | wc -l)
            print_success "Created $CURATION_COUNT curation files in $CURATION_DIR"
            
            echo
            echo "Next steps for manual curation:"
            echo "  1. Open: $CURATION_DIR/master_curation.md"
            echo "  2. Edit individual curation files to fill in expected domain boundaries"
            echo "  3. Mark proteins as test-ready when curation is complete"
            echo "  4. Run: python curation_validator.py --check-progress"
        fi
    else
        print_error "Failed to generate curation templates"
        return 1
    fi
}

show_progress() {
    print_header "Current Progress"
    
    echo "Checking curation progress..."
    if python curation_validator.py --check-progress 2>/dev/null; then
        echo
    else
        print_warning "No curation progress yet - run curation generation first"
    fi
    
    echo "Checking test results..."
    if python expand_test_suite.py --summary 2>/dev/null; then
        echo
    else
        print_warning "No test results yet - run algorithm tests first"
    fi
}

show_next_steps() {
    print_header "Next Steps"
    
    echo "Your test suite expansion is ready for manual curation!"
    echo
    echo "Manual Curation Workflow:"
    echo "  1. cd $OUTPUT_DIR/curation"
    echo "  2. Open master_curation.md to see progress"
    echo "  3. Edit individual *_curation.md files:"
    echo "     - Fill in expected domain boundaries under 'Manual Curation'"
    echo "     - Check validation boxes when ready"
    echo "  4. python curation_validator.py --check-progress"
    echo "  5. python curation_validator.py --validate-boundaries"
    echo "  6. python curation_validator.py --generate-tests"
    echo
    echo "Curation Tips:"
    echo "  - Use PyMOL to visualize protein structures"
    echo "  - Consult ECOD/Pfam databases for domain annotations"
    echo "  - Focus on clear structural boundaries"
    echo "  - Start with simple single-domain proteins"
    echo
    echo "Files created:"
    echo "  📊 Test results: $OUTPUT_DIR/test_results.json"
    echo "  📝 Curation files: $OUTPUT_DIR/curation/*.md"
    echo "  📋 Master tracking: $OUTPUT_DIR/curation/master_curation.md"
}

# Main workflow
main() {
    print_header "Mini PyECOD Test Suite Expansion"
    echo "This script will expand your test suite from 1 to 10-15 proteins"
    echo "with systematic testing and curation template generation."
    echo
    
    # Parse arguments
    SETUP_ONLY=false
    TEST_ONLY=false  
    CURATE_ONLY=false
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            --setup-only)
                SETUP_ONLY=true
                shift
                ;;
            --test-only)
                TEST_ONLY=true
                shift
                ;;
            --curate-only)
                CURATE_ONLY=true
                shift
                ;;
            --help)
                echo "Usage: $0 [--setup-only|--test-only|--curate-only]"
                echo
                echo "Options:"
                echo "  --setup-only    Only validate setup and exit"
                echo "  --test-only     Only run algorithm tests"
                echo "  --curate-only   Only generate curation templates"
                echo "  --help          Show this help"
                exit 0
                ;;
            *)
                print_error "Unknown option: $1"
                exit 1
                ;;
        esac
    done
    
    # Always check prerequisites
    check_prerequisites
    
    # Validate setup
    if ! validate_setup; then
        exit 1
    fi
    
    if $SETUP_ONLY; then
        print_success "Setup validation complete!"
        exit 0
    fi
    
    # Run tests
    if $TEST_ONLY || ! $CURATE_ONLY; then
        if ! run_tests; then
            print_warning "Some tests failed, but continuing..."
        fi
    fi
    
    if $TEST_ONLY; then
        show_progress
        exit 0
    fi
    
    # Generate curation templates
    if $CURATE_ONLY || ! $TEST_ONLY; then
        if ! generate_curation; then
            exit 1
        fi
    fi
    
    if $CURATE_ONLY; then
        show_progress
        exit 0
    fi
    
    # Show final status and next steps
    show_progress
    show_next_steps
    
    print_success "Test suite expansion setup complete!"
    echo
    echo "Time to start manual curation of domain boundaries."
    echo "See the workflow guide above for detailed instructions."
}

# Run main function
main "$@"
