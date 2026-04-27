"""
iPA - Integrated Processing and Analysis Toolkit

A comprehensive toolkit for multi-modal cellular and subcellular imaging data analysis.
"""

__version__ = "0.1.0"
__author__ = "liad"
__email__ = "liad@shanghaitech.edu.cn"

# Import main modules from the new structure
from . import data_loader
from . import processing
from . import analysis



__all__ = [
    'data_loader',
    'processing',
    'analysis',
    '__version__',
]
