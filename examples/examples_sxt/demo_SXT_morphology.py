"""
SXT Morphology Analysis Demo

This script demonstrates morphological analysis of organelles (e.g., ISG) 
using SXT data. It extracts features such as volume, surface area, and sphericity.
"""

import os
import sys
import numpy as np
import tifffile
import pandas as pd

# Add project root to path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
sys.path.insert(0, PROJECT_ROOT)

from ipa.analysis import extract_isg_features
from ipa.data_loader import QuickLogger

def main():
    # Initialize logger
    output_dir = os.path.join(PROJECT_ROOT, 'results', 'sxt_morphology_demo')
    os.makedirs(output_dir, exist_ok=True)
    
    logger = QuickLogger("sxt_morphology", log_dir=os.path.join(PROJECT_ROOT, 'logs'))
    logger.step("Starting SXT Morphology Analysis Demo")

    # --- Configuration ---
    data_id = '784_5'
    data_root = os.path.join(PROJECT_ROOT, 'data', 'sxt')
    isg_mask_path = os.path.join(data_root, f'{data_id}_isg_label.tiff')
    
    # Morphology Analysis: ISG Features
    logger.step("Running Morphology Analysis (ISG Features)...")
    if os.path.exists(isg_mask_path):
        isg_mask = tifffile.imread(isg_mask_path)
        isg_features_df, _ = extract_isg_features(
            isg_mask, 
            min_distance=5, 
            min_size=20, 
            save_csv=True, 
            output_path=os.path.join(output_dir, f'{data_id}_isg_morphology.csv')
        )
        
        # Convert volume column to numeric
        isg_features_df['Vol (pix)'] = pd.to_numeric(isg_features_df['Vol (pix)'], errors='coerce')
        valid_vols = isg_features_df['Vol (pix)'].dropna()
        
        logger.step(f"   Found {len(valid_vols)} ISGs.")
        logger.step(f"   Average Volume: {valid_vols.mean():.2f} voxels")
        logger.step(f"   Max Volume: {valid_vols.max():.2f} voxels")
    else:
        logger.step(f"   Skipping Morphology: {isg_mask_path} not found.")

    logger.step(f"Analysis Completed! Check '{output_dir}' folder for outputs.")

if __name__ == '__main__':
    main()
