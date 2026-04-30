"""
SIM Spatial Partitioning Analysis Demo

This script demonstrates how to perform cellular compartment partitioning analysis
on SIM (Structured Illumination Microscopy) imaging data. It divides the cellular 
space between nucleus and plasma membrane into radial partitions for downstream 
spatial analysis (e.g., RDF calculation).

Usage:
    python demo_SIM_partitioning.py

Requirements:
    - Input: Plasma membrane (PM) and nuclear envelope (NE) masks (.mrc format)
    - Output: Partition masks, coordinate files (.xvg), and visualizations
    
The script performs:
1. Boundary extraction from PM and NE masks
2. Radial partition creation using NE-PM pairs method
3. Coordinate extraction and XVG file generation
4. Visualization of partitions
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use Agg backend for non-interactive mode
import matplotlib.pyplot as plt
import mrcfile

# --- Path Configuration ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
sys.path.insert(0, PROJECT_ROOT)

from ipa.processing.partitioning import Partitioning, visualize_partitions
from ipa.data_loader import QuickLogger, UniversalDataLoader

def main():
    """Main function to run SIM partitioning analysis."""
    
    # Initialize logger
    log_dir = f'{PROJECT_ROOT}/logs'
    logger = QuickLogger("sim_partitioning", log_dir=log_dir)
    logger.step("Starting SIM Partitioning Analysis Demo")
    
    # Configuration
    root_dir = f'{PROJECT_ROOT}/data/sim/'
    dataid = '20220909_30-2-2-SIM'
    
    logger.step(f"Working with dataset: {dataid}")
    
    # Create results directory
    results_dir = f'{root_dir}results/'
    os.makedirs(results_dir, exist_ok=True)

    # File paths for PM and NE masks (MRC format for SIM, stored in sim/)
    mask_pm_path = f'{root_dir}{dataid}_PM.mrc'
    mask_ne_path = f'{root_dir}{dataid}_N.mrc'

    logger.step("Loading mask data")
    logger.file_in(mask_pm_path)
    logger.file_in(mask_ne_path)

    # Check if files exist
    if not os.path.exists(mask_pm_path) or not os.path.exists(mask_ne_path):
        logger.step(f"Warning: Mask files not found. Using synthetic data for demo.")
        # Generate synthetic masks for demonstration
        logger.step("Generating synthetic PM and NE masks...")
        shape = (64, 128, 128)
        mask_data_pm = np.zeros(shape, dtype=np.float32)
        mask_data_ne = np.zeros(shape, dtype=np.float32)
        
        # Create spherical masks
        z, y, x = np.ogrid[:shape[0], :shape[1], :shape[2]]
        center_z, center_y, center_x = shape[0]//2, shape[1]//2, shape[2]//2
        
        # NE: inner sphere (radius=20)
        ne_mask = ((z - center_z)**2 + (y - center_y)**2 + (x - center_x)**2) < 20**2
        mask_data_ne[ne_mask] = 1
        
        # PM: outer sphere (radius=40)
        pm_mask = ((z - center_z)**2 + (y - center_y)**2 + (x - center_x)**2) < 40**2
        mask_data_pm[pm_mask] = 1
    else:
        # Load real data (MRC format)
        mask_data_pm = mrcfile.open(mask_pm_path, mode='r', permissive=True).data
        mask_data_ne = mrcfile.open(mask_ne_path, mode='r', permissive=True).data
        logger.step(f"PM mask loaded: {mask_data_pm.shape}")
        logger.step(f"NE mask loaded: {mask_data_ne.shape}")

    # Visualize masks
    middle_slice = mask_data_pm.shape[0] // 2
    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    plt.imshow(mask_data_pm[middle_slice], cmap='gray')
    plt.title('Plasma Membrane (PM)')
    plt.axis('off')
    
    plt.subplot(1, 2, 2)
    plt.imshow(mask_data_ne[middle_slice], cmap='gray')
    plt.title('Nuclear Envelope (NE)')
    plt.axis('off')
    
    mask_viz_path = f'{results_dir}sim_mask_visualization.png'
    plt.savefig(mask_viz_path, dpi=150, bbox_inches='tight')
    logger.step(f"Mask visualization saved to: {mask_viz_path}")
    plt.close()

    # Initialize partitioning processor
    logger.step("Initializing partitioning processor")
    partitioner = Partitioning(root_dir)

    # Extract boundaries for analysis
    logger.step("Extracting NE and PM boundaries")
    center, ne_edge, pm_edge = partitioner.extract_ne_pm_edges(mask_data_pm, mask_data_ne)
    
    logger.step(f"NE boundary points: {len(ne_edge)}")
    logger.step(f"PM boundary points: {len(pm_edge)}")

    # Create radial partitions
    logger.step("Creating radial partitions (8 shells)")
    n_slices = 8
    partition_mask = partitioner.create_nepm_radial_partitions(
        ne_edge, pm_edge,
        shape=mask_data_pm.shape,
        n_slices=n_slices,
        pm_mask=mask_data_pm,
        ne_mask=mask_data_ne
    )
    
    logger.step(f"Partition mask created: {partition_mask.shape}")
    logger.step(f"Unique partition IDs: {np.unique(partition_mask)}")

    # Save partition mask
    partition_mask_path = f'{results_dir}{dataid}_partition_mask.mrc'
    with mrcfile.new(partition_mask_path, overwrite=True) as mrc:
        mrc.set_data(partition_mask.astype(np.uint16))
    logger.file_out(partition_mask_path)
    logger.step(f"Partition mask saved to: {partition_mask_path}")

    # Visualize partitions
    logger.step("Visualizing partitions")
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    middle_slice = partition_mask.shape[0] // 2
    
    # Original PM mask
    axes[0].imshow(mask_data_pm[middle_slice], cmap='gray')
    axes[0].set_title('Plasma Membrane')
    axes[0].axis('off')
    
    # Partition mask
    im = axes[1].imshow(partition_mask[middle_slice], cmap='viridis')
    axes[1].set_title(f'Radial Partitions ({n_slices} shells)')
    axes[1].axis('off')
    plt.colorbar(im, ax=axes[1])
    
    # Overlay
    overlay = mask_data_pm[middle_slice] * 0.3 + partition_mask[middle_slice] * 0.7
    axes[2].imshow(overlay, cmap='viridis')
    axes[2].set_title('Overlay: Partitions on PM')
    axes[2].axis('off')
    
    plt.tight_layout()
    comparison_path = f"{results_dir}{dataid}_partition_comparison.png"
    plt.savefig(comparison_path, dpi=300, bbox_inches='tight')
    logger.step(f"Partition comparison saved to: {comparison_path}")
    plt.close()

    # Save partition coordinates to XVG (for RDF analysis)
    logger.step("Saving partition coordinates to XVG format")
    coords_output = f'{results_dir}{dataid}_partition_coords.xvg'
    
    # Extract coordinates for each shell
    with open(coords_output, 'w') as f:
        f.write("# SIM Partition Coordinates\n")
        f.write(f"# Data ID: {dataid}\n")
        f.write(f"# Number of shells: {n_slices}\n")
        f.write("# Shell_ID  X  Y  Z\n")
        
        for shell_id in range(1, n_slices + 1):
            coords = np.argwhere(partition_mask == shell_id)
            for z, y, x in coords[::10]:  # Sample every 10th point to reduce file size
                f.write(f"{shell_id}  {x:.2f}  {y:.2f}  {z:.2f}\n")
    
    logger.file_out(coords_output)
    logger.step(f"Partition coordinates saved to: {coords_output}")

    logger.step("SIM Partitioning Analysis completed successfully!")
    print(f"\n✅ Results saved to: {results_dir}")
    print(f"   - Partition mask: {dataid}_partition_mask.mrc")
    print(f"   - Coordinates: {dataid}_partition_coords.xvg")
    print(f"   - Visualizations: {dataid}_partition_comparison.png")

if __name__ == '__main__':
    main()
