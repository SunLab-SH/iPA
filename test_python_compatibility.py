#!/usr/bin/env python
"""
Python Version Compatibility Test Script

This script tests iPA installation and basic functionality across different Python versions.
It verifies that all dependencies are compatible and core modules can be imported.

Usage:
    python test_python_compatibility.py [--python-version VERSION]

Example:
    # Test current Python version
    python test_python_compatibility.py
    
    # Test specific version (requires conda)
    python test_python_compatibility.py --python-version 3.11
"""

import sys
import subprocess
import platform
from pathlib import Path


def print_header(title):
    """Print formatted header"""
    print("\n" + "=" * 70)
    print(f"  {title}")
    print("=" * 70)


def print_section(title):
    """Print formatted section"""
    print(f"\n--- {title} ---")


def get_python_info():
    """Get current Python environment information"""
    return {
        'version': platform.python_version(),
        'implementation': platform.python_implementation(),
        'platform': platform.platform(),
        'executable': sys.executable,
    }


def test_import(module_name, optional=False):
    """Test if a module can be imported"""
    try:
        __import__(module_name)
        status = "✅ PASS"
        error = None
    except ImportError as e:
        if optional:
            status = "⚠️  OPTIONAL (not installed)"
        else:
            status = "❌ FAIL"
        error = str(e)
    
    return status, error


def check_package_version(package_name):
    """Get installed package version"""
    try:
        import importlib
        module = importlib.import_module(package_name.replace('-', '_'))
        version = getattr(module, '__version__', 'unknown')
        return version
    except Exception:
        return "not installed"


def run_tests():
    """Run all compatibility tests"""
    results = {
        'passed': 0,
        'failed': 0,
        'warnings': 0,
    }
    
    # Get Python info
    py_info = get_python_info()
    print_header("Python Environment Information")
    print(f"Python Version:    {py_info['version']}")
    print(f"Implementation:    {py_info['implementation']}")
    print(f"Platform:          {py_info['platform']}")
    print(f"Executable:        {py_info['executable']}")
    
    # Check Python version compatibility
    print_section("Python Version Check")
    major, minor = sys.version_info[:2]
    
    if major == 3 and 8 <= minor <= 11:
        print(f"✅ Python {major}.{minor} is supported (3.8-3.11)")
        results['passed'] += 1
    elif major == 3 and minor == 12:
        print(f"⚠️  Python {major}.{minor} is experimental (not officially supported yet)")
        results['warnings'] += 1
    else:
        print(f"❌ Python {major}.{minor} is NOT supported")
        print("   Supported versions: 3.8, 3.9, 3.10, 3.11")
        results['failed'] += 1
        return results
    
    # Test core dependencies
    print_section("Core Dependencies")
    core_packages = [
        ('numpy', False),
        ('scipy', False),
        ('pandas', False),
        ('matplotlib', False),
        ('skimage', False),  # scikit-image
        ('tifffile', False),
        ('PIL', False),  # Pillow
        ('cv2', True),  # opencv-python (optional for some features)
        ('mrcfile', False),
        ('plotly', False),
        ('seaborn', False),
        ('torch', False),
        ('torchvision', False),
        ('tqdm', False),
        ('h5py', False),
        ('sklearn', False),  # scikit-learn
    ]
    
    for pkg, optional in core_packages:
        status, error = test_import(pkg, optional)
        version = check_package_version(pkg)
        
        if "PASS" in status:
            results['passed'] += 1
            print(f"{status} | {pkg:20s} | v{version}")
        elif "OPTIONAL" in status:
            results['warnings'] += 1
            print(f"{status} | {pkg:20s} | {version}")
        else:
            results['failed'] += 1
            print(f"{status} | {pkg:20s} | Error: {error}")
    
    # Test file format support
    print_section("File Format Support")
    format_packages = [
        ('readlif', True),
        ('czifile', True),
        ('nd2reader', True),
        ('aicsimageio', True),
        ('imageio', False),
    ]
    
    for pkg, optional in format_packages:
        status, error = test_import(pkg, optional)
        version = check_package_version(pkg)
        
        if "PASS" in status:
            results['passed'] += 1
            print(f"{status} | {pkg:20s} | v{version}")
        elif "OPTIONAL" in status:
            results['warnings'] += 1
            print(f"{status} | {pkg:20s} | {version}")
        else:
            results['failed'] += 1
            print(f"{status} | {pkg:20s} | Error: {error}")
    
    # Test iPA modules
    print_section("iPA Module Imports")
    ipa_modules = [
        'ipa',
        'ipa.data_loader',
        'ipa.processing',
        'ipa.analysis',
    ]
    
    for module in ipa_modules:
        status, error = test_import(module, False)
        if "PASS" in status:
            results['passed'] += 1
            print(f"{status} | {module}")
        else:
            results['failed'] += 1
            print(f"{status} | {module} | Error: {error}")
    
    # Test key classes and functions
    print_section("Key API Components")
    api_tests = [
        ("UniversalDataLoader", "from ipa.data_loader import UniversalDataLoader"),
        ("QuickLogger", "from ipa.data_loader import QuickLogger"),
        ("Partitioning", "from ipa.processing.partitioning import Partitioning"),
        ("N2V", "from ipa.processing.denoising import N2V"),
        ("actin_to_vesicle_analysis", "from ipa.analysis import actin_to_vesicle_analysis"),
    ]
    
    for name, import_stmt in api_tests:
        try:
            exec(import_stmt)
            print(f"✅ PASS | {name}")
            results['passed'] += 1
        except Exception as e:
            print(f"❌ FAIL | {name} | Error: {e}")
            results['failed'] += 1
    
    # Print summary
    print_header("Test Summary")
    total = results['passed'] + results['failed'] + results['warnings']
    print(f"Total Tests:    {total}")
    print(f"Passed:         {results['passed']} ✅")
    print(f"Failed:         {results['failed']} ❌")
    print(f"Warnings:       {results['warnings']} ⚠️")
    
    if results['failed'] == 0:
        print("\n🎉 All critical tests passed!")
        if results['warnings'] > 0:
            print(f"⚠️  {results['warnings']} optional packages not installed (this is OK)")
        return True
    else:
        print(f"\n❌ {results['failed']} test(s) failed. Please check the errors above.")
        return False


def test_with_conda_env(python_version):
    """Test with a specific Python version using conda"""
    env_name = f"ipa_test_py{python_version.replace('.', '')}"
    
    print_header(f"Testing with Python {python_version}")
    print(f"Creating conda environment: {env_name}")
    
    try:
        # Create conda environment
        subprocess.run([
            'conda', 'create', '-y', '-n', env_name, 
            f'python={python_version}'
        ], check=True)
        
        # Activate and install iPA
        print(f"Installing iPA in {env_name}...")
        subprocess.run([
            'conda', 'run', '-n', env_name,
            'pip', 'install', '-e', '.'
        ], check=True, cwd=Path(__file__).parent)
        
        # Run tests
        print(f"Running tests in {env_name}...")
        result = subprocess.run([
            'conda', 'run', '-n', env_name,
            'python', 'test_python_compatibility.py'
        ], cwd=Path(__file__).parent)
        
        # Clean up
        print(f"Removing test environment: {env_name}")
        subprocess.run(['conda', 'env', 'remove', '-y', '-n', env_name], check=True)
        
        return result.returncode == 0
        
    except subprocess.CalledProcessError as e:
        print(f"❌ Error during testing: {e}")
        return False


def main():
    """Main entry point"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Test iPA Python version compatibility'
    )
    parser.add_argument(
        '--python-version',
        type=str,
        default=None,
        help='Test with specific Python version (requires conda). E.g., 3.11'
    )
    
    args = parser.parse_args()
    
    if args.python_version:
        # Test with specific Python version using conda
        success = test_with_conda_env(args.python_version)
        sys.exit(0 if success else 1)
    else:
        # Test current Python version
        print_header("iPA Python Version Compatibility Test")
        success = run_tests()
        sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
