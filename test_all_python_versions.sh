#!/bin/bash
# Multi-Python Version Compatibility Test Script
# Tests iPA across Python 3.8, 3.9, 3.10, and 3.11

set -e  # Exit on error

echo "======================================================================"
echo "  iPA Multi-Python Version Compatibility Test"
echo "======================================================================"
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Python versions to test
PYTHON_VERSIONS=("3.8" "3.9" "3.10" "3.11")

# Results tracking
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

# Function to test a specific Python version
test_python_version() {
    local py_version=$1
    local env_name="ipa_test_py${py_version//./}"
    
    echo ""
    echo "----------------------------------------------------------------------"
    echo -e "${YELLOW}Testing Python ${py_version}${NC}"
    echo "----------------------------------------------------------------------"
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    # Check if conda is available
    if ! command -v conda &> /dev/null; then
        echo -e "${RED}❌ Conda is not installed. Please install Miniconda or Anaconda.${NC}"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        return 1
    fi
    
    # Create conda environment
    echo "Creating conda environment: $env_name"
    if ! conda create -y -n "$env_name" python="$py_version" > /dev/null 2>&1; then
        echo -e "${RED}❌ Failed to create conda environment${NC}"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        return 1
    fi
    
    # Install iPA
    echo "Installing iPA..."
    if ! conda run -n "$env_name" pip install -e . > /dev/null 2>&1; then
        echo -e "${RED}❌ Failed to install iPA${NC}"
        conda env remove -y -n "$env_name" > /dev/null 2>&1
        FAILED_TESTS=$((FAILED_TESTS + 1))
        return 1
    fi
    
    # Run compatibility test
    echo "Running compatibility tests..."
    if conda run -n "$env_name" python test_python_compatibility.py; then
        echo -e "${GREEN}✅ Python ${py_version} tests PASSED${NC}"
        PASSED_TESTS=$((PASSED_TESTS + 1))
    else
        echo -e "${RED}❌ Python ${py_version} tests FAILED${NC}"
        FAILED_TESTS=$((FAILED_TESTS + 1))
    fi
    
    # Clean up
    echo "Removing test environment..."
    conda env remove -y -n "$env_name" > /dev/null 2>&1
    
    return 0
}

# Main execution
echo "This script will test iPA compatibility with multiple Python versions."
echo "Each test will:"
echo "  1. Create a fresh conda environment"
echo "  2. Install iPA and dependencies"
echo "  3. Run compatibility tests"
echo "  4. Remove the test environment"
echo ""
echo "Python versions to test: ${PYTHON_VERSIONS[*]}"
echo ""
read -p "Continue? (y/n) " -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Aborted."
    exit 0
fi

# Test each Python version
for version in "${PYTHON_VERSIONS[@]}"; do
    test_python_version "$version" || true
done

# Print summary
echo ""
echo "======================================================================"
echo "  Test Summary"
echo "======================================================================"
echo "Total Versions Tested:  $TOTAL_TESTS"
echo -e "Passed:               ${GREEN}${PASSED_TESTS}${NC}"
echo -e "Failed:               ${RED}${FAILED_TESTS}${NC}"
echo ""

if [ $FAILED_TESTS -eq 0 ]; then
    echo -e "${GREEN}🎉 All Python version tests passed!${NC}"
    exit 0
else
    echo -e "${RED}❌ Some tests failed. Please check the output above.${NC}"
    exit 1
fi
