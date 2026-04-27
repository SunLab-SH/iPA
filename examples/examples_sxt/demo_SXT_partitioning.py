#%%

# examples/basic_usage/demo_partitioning.py
import numpy as np
from ipa.processing.partitioning import Partitioning
from ipa.processing.partitioning import visualize_partitions, plot_partition_features
import matplotlib
matplotlib.use('TkAgg')  # Use TkAgg backend to support graphical display
import matplotlib.pyplot as plt
import tifffile
import os

import sys
from parsers import get_args
from ipa.data_loader import QuickLogger

def main():
    mainpath = get_args().main_path
    
    # Initialize logger
    log_dir = f'{mainpath}/logs'
    logger = QuickLogger("sxt_partitioning", log_dir=log_dir)
    logger.step("Starting SXT Partitioning Demo")
    
    logger.step(f"Main path: {mainpath}")
    root_dir = f'{mainpath}/data/sxt_images/'  # Output directory, recommended to create new
    dataid = '784_5'
    logger.step(f"Processing dataset: {dataid}")
    
    # Create results directory
    results_dir = f'{root_dir}/results/'
    os.makedirs(results_dir, exist_ok=True)
    logger.step(f"Results directory: {results_dir}")

    mask_pm_path = f'{mainpath}/data/sxt_images/{dataid}_wholecell_label.tiff'
    mask_ne_path = f'{mainpath}/data/sxt_images/{dataid}_NC_label.tiff'

    logger.file_in(mask_pm_path)
    logger.file_in(mask_ne_path)
    logger.step("Loading mask data...")
    
    mask_data_pm = tifffile.imread(mask_pm_path).astype(int)
    mask_data_ne = tifffile.imread(mask_ne_path).astype(int)

    logger.step(f"PM mask shape: {mask_data_pm.shape}, dtype: {mask_data_pm.dtype}")
    logger.step(f"NE mask shape: {mask_data_ne.shape}, dtype: {mask_data_ne.dtype}")
    
    # Add data validation
    logger.step(f"PM mask unique values: {np.unique(mask_data_pm)}")
    logger.step(f"NE mask unique values: {np.unique(mask_data_ne)}")
    logger.step(f"PM mask non-zero voxels: {np.sum(mask_data_pm > 0)}")
    logger.step(f"NE mask non-zero voxels: {np.sum(mask_data_ne > 0)}")

    # Convert boolean masks to int for proper addition
    plt.figure(figsize=(12, 4))
    
    plt.subplot(1, 3, 1)
    plt.imshow(mask_data_pm[255], cmap='Reds')
    plt.title(f"PM Mask - Slice 255")
    plt.axis('off')
    
    plt.subplot(1, 3, 2)
    plt.imshow(mask_data_ne[255], cmap='Blues')
    plt.title(f"NE Mask - Slice 255")
    plt.axis('off')
    
    plt.subplot(1, 3, 3)
    plt.imshow(mask_data_pm[255] + mask_data_ne[255], cmap='gray')
    plt.title(f"Combined Masks - Data ID: {dataid}")
    plt.axis('off')
    
    plt.tight_layout()
    plt.show()

    # 2. Initialize partitioning processor
    logger.step("Initializing partitioning processor...")
    partitioner = Partitioning(root_dir)

    # 3. Extract boundaries for analysis
    logger.step("Extracting boundaries for analysis...")
    center, ne_edge, pm_edge = partitioner.extract_ne_pm_edges(mask_data_pm, mask_data_ne)
    
    logger.step(f"Center shape: {center.shape if hasattr(center, 'shape') else type(center)}")
    logger.step(f"NE edge shape: {ne_edge.shape if hasattr(ne_edge, 'shape') else type(ne_edge)}")
    logger.step(f"PM edge shape: {pm_edge.shape if hasattr(pm_edge, 'shape') else type(pm_edge)}")
    
 

    logger.step("Creating radial partitions using NE-PM pairs...")
    partition_masks = partitioner.create_nepm_radial_partitions(
        ne_edge, pm_edge, mask_data_pm.shape,
        n_slices=8, pm_mask=mask_data_pm, ne_mask=mask_data_ne
    )
    method_name = "shell_based_radial_partitioning"

    
    logger.step(f"Partition masks shape: {partition_masks.shape}")
    logger.step(f"Unique partition values: {np.unique(partition_masks)}")
    logger.step(f"Method used: {method_name}")

    # Extract coordinate points from partition masks and save as XVG format
    logger.step("Extracting partition coordinates for XVG output...")
    partition_coords = partitioner.extract_partition_coordinates(partition_masks, sampling_density=1)
    xvg_path = partitioner.save_partition_coords_to_xvg(partition_coords, f"{dataid}_{method_name}", results_dir)
    logger.file_out(xvg_path)

    # 5. Visualize partitions
    logger.step("Generating visualizations...")
    
    # Add direct visualization of partition masks
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Display original masks
    middle_slice = mask_data_pm.shape[0] // 2

    unique_partitions = np.unique(partition_masks)
    n_partitions = len(unique_partitions[unique_partitions != 0])  # Exclude background (0)
    logger.step(f"Number of partitions: {n_partitions}")

    axes[0].imshow(mask_data_pm[middle_slice] + mask_data_ne[middle_slice], cmap='gray')
    axes[0].set_title('Original Masks (PM + NE)')
    axes[0].axis('off')
    
    # Display partition masks
    axes[1].imshow(partition_masks[middle_slice], cmap='tab10', vmin=0, vmax=n_partitions)
    axes[1].set_title(f'Partition Masks ({method_name}, Slice {middle_slice})')
    axes[1].axis('off')
    
    # Display overlay image
    overlay = mask_data_pm[middle_slice] * 0.3 + partition_masks[middle_slice] * 0.7
    axes[2].imshow(overlay, cmap='viridis')
    axes[2].set_title(f'Overlay: {method_name.upper()} Partitions on Cell')
    axes[2].axis('off')
    
    plt.tight_layout()
    comparison_path = f"{results_dir}{dataid}_{method_name}_partition_comparison.png"
    plt.savefig(comparison_path, dpi=300, bbox_inches='tight')
    logger.file_out(comparison_path)
    plt.show()
    
    # Check voxel count for each partition
    for partition_id in unique_partitions:
        if partition_id != 0:
            count = np.sum(partition_masks == partition_id)
            logger.step(f"Partition {partition_id}: {count} voxels")
    
    # 2D visualization at middle slice
    logger.step(f"Visualizing slice {middle_slice} of {mask_data_pm.shape[0]} slices")
    
    # Check the slice data before visualization
    slice_data = partition_masks[middle_slice, :, :]
    logger.step(f"Slice data shape: {slice_data.shape}")
    logger.step(f"Slice data unique values: {np.unique(slice_data)}")
    
    # Convert partition_masks to list format expected by visualize_partitions
    unique_partitions = np.unique(partition_masks)
    partition_list = []
    for partition_id in unique_partitions:
        if partition_id != 0:  # Skip background
            mask = (partition_masks == partition_id)
            partition_list.append(mask)
    
    logger.step(f"Created {len(partition_list)} partition masks for visualization")
    
    # 2D visualization
    viz_2d_path = f"{results_dir}{dataid}_{method_name}_partitions_2d.png"
    visualize_partitions(partition_list, slice_idx=middle_slice, save_path=viz_2d_path)
    logger.file_out(viz_2d_path)

    logger.step(f"Processing completed for dataset: {dataid} using {method_name} method")

if __name__ == '__main__':
    main()

