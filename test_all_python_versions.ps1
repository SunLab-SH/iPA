# Multi-Python Version Compatibility Test Script (PowerShell)
# Tests iPA across Python 3.8, 3.9, 3.10, and 3.11

Write-Host "======================================================================" -ForegroundColor Cyan
Write-Host "  iPA Multi-Python Version Compatibility Test (Windows)" -ForegroundColor Cyan
Write-Host "======================================================================" -ForegroundColor Cyan
Write-Host ""

# Python versions to test
$PYTHON_VERSIONS = @("3.8", "3.9", "3.10", "3.11")

# Results tracking
$TOTAL_TESTS = 0
$PASSED_TESTS = 0
$FAILED_TESTS = 0

# Function to test a specific Python version
function Test-PythonVersion {
    param(
        [string]$py_version
    )
    
    $env_name = "ipa_test_py" + $py_version.Replace(".", "")
    
    Write-Host ""
    Write-Host "----------------------------------------------------------------------" -ForegroundColor Yellow
    Write-Host "Testing Python $py_version" -ForegroundColor Yellow
    Write-Host "----------------------------------------------------------------------" -ForegroundColor Yellow
    
    $script:TOTAL_TESTS++
    
    # Check if conda is available
    if (-not (Get-Command conda -ErrorAction SilentlyContinue)) {
        Write-Host "❌ Conda is not installed. Please install Miniconda or Anaconda." -ForegroundColor Red
        $script:FAILED_TESTS++
        return $false
    }
    
    try {
        # Create conda environment
        Write-Host "Creating conda environment: $env_name"
        conda create -y -n $env_name python=$py_version | Out-Null
        
        # Install iPA
        Write-Host "Installing iPA..."
        conda run -n $env_name pip install -e . | Out-Null
        
        # Run compatibility test
        Write-Host "Running compatibility tests..."
        $test_result = conda run -n $env_name python test_python_compatibility.py
        
        if ($LASTEXITCODE -eq 0) {
            Write-Host "✅ Python $py_version tests PASSED" -ForegroundColor Green
            $script:PASSED_TESTS++
            $success = $true
        } else {
            Write-Host "❌ Python $py_version tests FAILED" -ForegroundColor Red
            $script:FAILED_TESTS++
            $success = $false
        }
        
        # Clean up
        Write-Host "Removing test environment..."
        conda env remove -y -n $env_name | Out-Null
        
        return $success
        
    } catch {
        Write-Host "❌ Error during testing: $_" -ForegroundColor Red
        $script:FAILED_TESTS++
        
        # Try to clean up even on error
        try {
            conda env remove -y -n $env_name 2>$null
        } catch {}
        
        return $false
    }
}

# Main execution
Write-Host "This script will test iPA compatibility with multiple Python versions." -ForegroundColor White
Write-Host "Each test will:" -ForegroundColor White
Write-Host "  1. Create a fresh conda environment" -ForegroundColor White
Write-Host "  2. Install iPA and dependencies" -ForegroundColor White
Write-Host "  3. Run compatibility tests" -ForegroundColor White
Write-Host "  4. Remove the test environment" -ForegroundColor White
Write-Host ""
Write-Host "Python versions to test: $($PYTHON_VERSIONS -join ', ')" -ForegroundColor White
Write-Host ""

$confirmation = Read-Host "Continue? (y/n)"
if ($confirmation -notmatch '^[Yy]') {
    Write-Host "Aborted." -ForegroundColor Yellow
    exit 0
}

# Test each Python version
foreach ($version in $PYTHON_VERSIONS) {
    Test-PythonVersion -py_version $version
}

# Print summary
Write-Host ""
Write-Host "======================================================================" -ForegroundColor Cyan
Write-Host "  Test Summary" -ForegroundColor Cyan
Write-Host "======================================================================" -ForegroundColor Cyan
Write-Host "Total Versions Tested:  $TOTAL_TESTS" -ForegroundColor White
Write-Host "Passed:               $PASSED_TESTS" -ForegroundColor Green
Write-Host "Failed:               $FAILED_TESTS" -ForegroundColor Red
Write-Host ""

if ($FAILED_TESTS -eq 0) {
    Write-Host "🎉 All Python version tests passed!" -ForegroundColor Green
    exit 0
} else {
    Write-Host "❌ Some tests failed. Please check the output above." -ForegroundColor Red
    exit 1
}
