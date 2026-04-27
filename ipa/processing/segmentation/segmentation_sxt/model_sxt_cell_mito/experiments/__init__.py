"""
SXT Model Experiments Module

This module provides segmentation experiments for SXT (Soft X-ray Tomography) data,
including membrane/nucleus segmentation and mitochondria segmentation.
"""

import os
import sys
from pathlib import Path



# Import membrane and nucleus segmentation
from .mem_nu_seg_with_val import run_cell_segmentation


from.mito_seg_with_val import run_mito_segmentation
from .parser import args, args_nu, args_mito

