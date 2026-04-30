"""
SIM/WFM Interaction Analysis Demo

This script demonstrates how to analyze the spatial relationship between 
ISG (Insulin Secretory Granules) and the Plasma Membrane (PM) in SIM/WFM data.

It performs the following steps:
1. Loads real ISG and PM segmentation masks from SIM data.
2. Extracts ISG granule centers and volumes.
3. Calculates the distance from each granule to the PM.
4. Classifies granules into 'Docked' vs 'Reserve' pools using a fixed threshold.
"""

import os
import sys
import numpy as np
from skimage import measure

# --- Path Configuration ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
sys.path.insert(0, PROJECT_ROOT)

from ipa.analysis import analyze_docked_granules
from ipa.data_loader import QuickLogger

def extract_isg_features_from_mask(isg_mask):
    """Helper function to extract center coordinates and volume from a 3D binary mask."""
    labeled_mask = measure.label(isg_mask > 0)
    regions = measure.regionprops(labeled_mask)
    
    # Prepare data in the format: [Index, ..., CX, CY, CZ, ..., Vol]
    data_rows = []
    for i, region in enumerate(regions):
        if region.area < 5: continue # Filter noise
        
        z, y, x = region.centroid
        row = [0] * 20
        row[0] = i
        row[4], row[5], row[6] = x, y, z # CX, CY, CZ
        row[19] = region.area # Volume in voxels
        data_rows.append(row)
        
    return np.array(data_rows)

def main():
    logger = QuickLogger("sim_interaction_demo", log_dir=os.path.join(PROJECT_ROOT, 'logs'))
    logger.step("Starting SIM/WFM Interaction Analysis Demo")

    # --- 1. Data Preparation ---
    data_root = os.path.join(PROJECT_ROOT, 'data', 'sim')
    
    isg_path = os.path.join(data_root, '20220909_30-2-1-SIM_ISG_seg.tif')
    pm_path = os.path.join(data_root, '20220909_30-2-1-SIM_PM.mrc')
    
    if not os.path.exists(isg_path) or not os.path.exists(pm_path):
        logger.step(f"Error: Real SIM data files not found at {data_root}")
        return

    import tifffile
    import mrcfile
    
    isg_mask = tifffile.imread(isg_path)
    
    with mrcfile.open(pm_path, permissive=True) as mrc:
        pm_mask = mrc.data > 0 # Convert to binary

    logger.step(f"Loaded ISG mask shape: {isg_mask.shape}")
    logger.step(f"Loaded PM mask shape: {pm_mask.shape}")

    # Extract features from the mask
    isg_data = extract_isg_features_from_mask(isg_mask)
    logger.step(f"Extracted {len(isg_data)} ISG granules")

    # --- 2. Docking Analysis ---
    logger.step("Analyzing ISG-PM docking...")
    
    results = analyze_docked_granules(
        img_isg=isg_mask,
        mem_mask=pm_mask.astype(np.uint8),
        isg_data=isg_data,
        distance_threshold=5, # Fixed threshold in pixels
        show_visualization=False
    )
    
    logger.step("✅ Analysis Complete!")
    logger.step(f"   Total Granules: {results['total_granules']}")
    logger.step(f"   Docked Granules (RRP): {results['docked_granules']}")
    logger.step(f"   Docking Ratio: {results['docked_ratio']:.2%}")

    # Save results
    output_dir = os.path.join(PROJECT_ROOT, 'results', 'sim_interaction_demo')
    os.makedirs(output_dir, exist_ok=True)
    
    import json
    output_path = os.path.join(output_dir, 'docking_results.json')
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=4)
    
    logger.file_out(output_path)
    logger.step("Example finished successfully.")

if __name__ == '__main__':
    main()
