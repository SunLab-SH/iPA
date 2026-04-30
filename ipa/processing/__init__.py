"""
Processing Module

Image processing module including segmentation, denoising, and partitioning.
"""

# Import submodules
from . import segmentation
from . import denoising
from . import partitioning

# Export commonly used classes/functions for convenient access
from .partitioning import Partitioning, visualize_partitions, plot_partition_features
from .denoising import N2V, N2N, predict_3d
from .segmentation import create_segmenter

__all__ = [
    # Submodules
    'segmentation',
    'denoising',
    'partitioning',
    
    # Commonly used classes/functions (Unified API)
    'Partitioning',
    'visualize_partitions',
    'plot_partition_features',
    'N2V',
    'N2N',
    'predict_3d',
    'create_segmenter',
]