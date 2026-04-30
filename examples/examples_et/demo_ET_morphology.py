"""
Cryo-ET Morphology Analysis Demo

This script demonstrates morphological analysis of tubular structures 
(e.g., Actin filaments) from Cryo-ET data. It calculates total length 
and individual filament lengths.
"""

import os
import sys
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Add project root to path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
sys.path.insert(0, PROJECT_ROOT)

from ipa.analysis import calculate_tubular_length
from ipa.data_loader import QuickLogger

def main():
    # Initialize logger
    output_dir = os.path.join(PROJECT_ROOT, 'results', 'et_morphology_demo')
    os.makedirs(output_dir, exist_ok=True)

    logger = QuickLogger("et_morphology", log_dir=os.path.join(PROJECT_ROOT, 'logs'))
    logger.step("Starting Cryo-ET Morphology Analysis Demo")

    # --- 1. Data Preparation ---
    data_root = os.path.join(PROJECT_ROOT, 'data', 'cryoET')
    json_path = os.path.join(data_root, '20210517_30_009_actin_filled_points.json')
    
    if not os.path.exists(json_path):
        logger.step(f"Error: Input file not found: {json_path}")
        return

    logger.step(f"Loading actin coordinates from: {os.path.basename(json_path)}")
    with open(json_path, 'r') as f:
        actin_data = json.load(f)

    # Cryo-ET voxel size: 1.34 nm (from Manuscript)
    voxel_size = (1.34, 1.34, 1.34) 
    logger.step(f"Using voxel size: {voxel_size} nm")

    # --- 2. Skeleton Mask Reconstruction ---
    logger.step("Reconstructing 3D skeleton mask...")
    
    all_coords = []
    if isinstance(actin_data, dict):
        for points in actin_data.values():
            all_coords.extend(points)
    else:
        all_coords = actin_data

    if not all_coords:
        logger.step("Error: No coordinates found.")
        return

    coords_array = np.array(all_coords)
    
    min_vals = np.min(coords_array, axis=0).astype(int)
    max_vals = np.max(coords_array, axis=0).astype(int)
    
    shape = (max_vals[0] - min_vals[0] + 5, 
             max_vals[1] - min_vals[1] + 5, 
             max_vals[2] - min_vals[2] + 5)
    
    skeleton_mask = np.zeros(shape, dtype=np.uint8)
    
    for coord in coords_array:
        x, y, z = int(coord[0]) - min_vals[0], int(coord[1]) - min_vals[1], int(coord[2]) - min_vals[2]
        if 0 <= x < shape[0] and 0 <= y < shape[1] and 0 <= z < shape[2]:
            skeleton_mask[x, y, z] = 1

    logger.step(f"   Mask shape: {skeleton_mask.shape}")

    # --- 3. Length Calculation ---
    logger.step("Calculating tubular lengths...")
    
    total_length, individual_lengths = calculate_tubular_length(
        skeleton_mask, 
        voxel_size=voxel_size, 
        return_individual=True
    )
    
    logger.step("✅ Analysis Complete!")
    logger.step(f"   Total Length: {total_length:.2f} nm")
    logger.step(f"   Filament Count: {len(individual_lengths)}")

    # --- 4. Visualization ---
    viz_path = os.path.join(output_dir, 'actin_skeleton_example.png')
    
    mid_z = shape[2] // 2
    plt.figure(figsize=(8, 6))
    plt.imshow(skeleton_mask[:, :, mid_z], cmap='hot', interpolation='nearest')
    plt.title(f"Actin Skeleton Slice (Z={mid_z})\nTotal Length: {total_length:.1f} nm")
    plt.axis('off')
    plt.savefig(viz_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    logger.file_out(viz_path)
    logger.step("Demo finished successfully.")

if __name__ == '__main__':
    main()
