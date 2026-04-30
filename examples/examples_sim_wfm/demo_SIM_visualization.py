"""
SIM Data 3D Visualization Demo
"""

import sys
import os
import numpy as np
import tifffile
import mrcfile
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --- Path Configuration ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
sys.path.insert(0, PROJECT_ROOT)

from ipa.processing.partitioning import visualize_complete_scene
from ipa.data_loader import UniversalDataLoader
from ipa.data_loader import QuickLogger

def main():
    """Main function to run the complete SIM visualization pipeline."""
    # Initialize logger
    log_dir = f'{PROJECT_ROOT}/logs'
    logger = QuickLogger("sim_visualization", log_dir=log_dir)
    logger.step("Starting SIM Data 3D Visualization Demo")
    
    # Configuration
    rootpath = f'{PROJECT_ROOT}/data/sim/'
    dataid = "20220909_30-2-1-SIM"
    
    # File paths
    isg_file_path = f"{rootpath}/{dataid}_ISG_seg.tif"
    pm_mask_file_path = f'{rootpath}/{dataid[:-4]}-SIM_PM.mrc'
    ne_mask_file_path = f'{rootpath}/{dataid[:-4]}-SIM_N.mrc'
    
    logger.step(f"Working with dataset: {dataid}")
    logger.file_in(isg_file_path)
    logger.file_in(pm_mask_file_path)
    logger.file_in(ne_mask_file_path)
    
    # Load data
    logger.step("Loading data files")
    img_isg = UniversalDataLoader.load_data(isg_file_path)
    mem_mask = UniversalDataLoader.load_data(pm_mask_file_path)
    ne_mask = UniversalDataLoader.load_data(ne_mask_file_path)

    logger.step(f"ISG data loaded: {img_isg.shape}")
    logger.step(f"Plasma membrane mask loaded: {mem_mask.shape}")
    logger.step(f"Nuclear envelope mask loaded: {ne_mask.shape}")

    # Preprocess data
    logger.step("Preprocessing data")
    nonzero_coords = np.where(mem_mask != 0)
    x0 = (nonzero_coords[1].min() // 100) * 100
    y0 = (nonzero_coords[2].min() // 100) * 100
    x1 = (nonzero_coords[1].max() // 100) * 100 + 100
    y1 = (nonzero_coords[2].max() // 100) * 100 + 100
    
    logger.step(f"Cropping region: x[{x0}:{x1}], y[{y0}:{y1}]")
    
    # Crop and process masks
    mem_mask_cropped = mem_mask[:, x0:x1, y0:y1]
    img_isg = np.flip(img_isg, axis=0)
    img_isg_cropped = img_isg[:, x0:x1, y0:y1]
    ne_mask_cropped = ne_mask[:, x0:x1, y0:y1]
    
    # Create binary masks
    mem_mask_binary = (mem_mask_cropped > 0).astype(np.uint8)
    ne_mask_binary = (ne_mask_cropped > 0).astype(np.uint8)
    
    isg_mask_binary = img_isg_cropped.copy()
    isg_mask_binary[isg_mask_binary > 0] = 1
    isg_mask_binary[mem_mask_binary < 1] = 0  # Remove ISG outside membrane
    isg_mask_binary[ne_mask_binary > 0] = 0  # Remove ISG inside nucleus
    
    logger.step(f"Processed shapes: Membrane {mem_mask_binary.shape}, ISG {isg_mask_binary.shape}")
    
    # Quick alignment check
    logger.step("Creating alignment check visualization")
    mid_z = isg_mask_binary.shape[0] // 2
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    axes[0].imshow(isg_mask_binary[mid_z], cmap='Blues')
    axes[0].set_title(f'ISG (Z={mid_z})')
    axes[1].imshow(mem_mask_binary[mid_z], cmap='Reds')
    axes[1].set_title(f'Membrane (Z={mid_z})')
    
    overlay = np.zeros((*mem_mask_binary[mid_z].shape, 3))
    overlay[:, :, 0] = mem_mask_binary[mid_z] * 0.6
    overlay[:, :, 2] = isg_mask_binary[mid_z] * 0.8
    overlay[:, :, 1] = ne_mask_binary[mid_z] * 0.4
    axes[2].imshow(overlay)
    axes[2].set_title(f'Overlay (Z={mid_z})')
    
    for ax in axes:
        ax.axis('off')
    plt.tight_layout()
    viz_path = f'{PROJECT_ROOT}/results/sim_visualization_2d.png'
    os.makedirs(os.path.dirname(viz_path), exist_ok=True)
    plt.savefig(viz_path, dpi=150, bbox_inches='tight')
    logger.step(f"2D visualization saved to: {viz_path}")
    plt.close()
    
    # 3D Visualization
    logger.step("Creating 3D visualization")
    masks = [mem_mask_binary, ne_mask_binary, isg_mask_binary]
    mask_names = ['Cell Membrane', 'Nuclear Envelope', 'ISG Data']
    
    # Create visualizations with different z-scales
    logger.step("Generating visualization with z-scale=1")
    visualize_complete_scene(
        masks=masks,
        mask_names=mask_names,
        title=f"Complete_Scene_{dataid}",
        z_scale=1,
        save=True
    )
    
    logger.step("Generating visualization with z-scale=5")
    visualize_complete_scene(
        masks=masks,
        mask_names=mask_names,
        title=f"ZScale5_{dataid}",
        z_scale=5,
        save=True
    )
    
    logger.step("Visualization completed successfully")
    

if __name__ == "__main__":
    main()
