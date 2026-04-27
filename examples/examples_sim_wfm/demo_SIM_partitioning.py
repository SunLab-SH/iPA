#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
SIM Cellular Partitioning Analysis Demo

This script demonstrates how to perform cellular compartment partitioning analysis
on SIM imaging data. It divides the cellular space between nucleus and plasma 
membrane into radial partitions for downstream spatial analysis.

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
import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # Use TkAgg backend to support graphical display
import matplotlib.pyplot as plt

from ipa.processing.partitioning import Partitioning, visualize_partitions
from ipa.data_loader import UniversalDataLoader, QuickLogger
from parsers import get_args


def main():
    """
    Main function to perform cellular partitioning analysis on SIM data.
    
    This function:
    1. Loads PM and NE mask data
    2. Extracts boundaries
    3. Creates radial partitions
    4. Saves partition coordinates
    5. Generates visualizations
    """
    mainpath = get_args().main_path
    
    # Initialize logger
    log_dir = f'{mainpath}/logs'
    logger = QuickLogger("sim_partitioning_analysis", log_dir=log_dir)
    logger.step("=" * 60)
    logger.step("SIM Cellular Partitioning Analysis Demo")
    logger.step("=" * 60)
    
    root_dir = f'{mainpath}/data/sim_images/'  # Output directory, recommended to create new
    dataid = '20220909_30-2-1-SIM'
    
    logger.step(f"Working with dataset: {dataid}")
    
    # Create results directory
    results_dir = f'{root_dir}/results/'
    os.makedirs(results_dir, exist_ok=True)

    mask_pm_path = f'{mainpath}/data/sim_images/{dataid}_PM.mrc'
    mask_ne_path = f'{mainpath}/data/sim_images/{dataid}_N.mrc'

    logger.step("Loading mask data")
    logger.file_in(mask_pm_path)
    logger.file_in(mask_ne_path)

    mask_data_pm = UniversalDataLoader.load_data(mask_pm_path)
    mask_data_ne = UniversalDataLoader.load_data(mask_ne_path)

    logger.step(f"PM mask loaded: {mask_data_pm.shape}, dtype: {mask_data_pm.dtype}")
    logger.step(f"NE mask loaded: {mask_data_ne.shape}, dtype: {mask_data_ne.dtype}")
    
    # Add data validation
    logger.step(f"PM mask unique values: {np.unique(mask_data_pm)}")
    logger.step(f"NE mask unique values: {np.unique(mask_data_ne)}")
    logger.step(f"PM mask non-zero voxels: {np.sum(mask_data_pm > 0)}")
    logger.step(f"NE mask non-zero voxels: {np.sum(mask_data_ne > 0)}")
    
    print(mask_data_pm.shape, mask_data_ne.shape)
    print(f"PM mask dtype: {mask_data_pm.dtype}")
    print(f"NE mask dtype: {mask_data_ne.dtype}")

    # # 1. Load mask data
    # mask_data = np.load(mask_path)  # Shape: (T, Z, 2, X, Y)
    
    # mask_data_pm = mask_data[tp, :, 0, :, :]  # Extract mask data from frame 0
    # mask_data_ne = mask_data[tp, :, 1, :, :]  # Extract mask data from frame 0
    # mask_data_pm = np.transpose(mask_data_pm, (0, 1, 2))
    # mask_data_ne = np.transpose(mask_data_ne, (0, 1, 2))
    tp = 0  # Process only frame 0, use loop for multiple frames

    # Convert boolean masks to int for proper addition
    plt.figure(figsize=(12, 4))
    
    plt.subplot(1, 3, 1)
    plt.imshow(mask_data_pm[10], cmap='Reds')
    plt.title(f"PM Mask - Slice 10")
    plt.axis('off')
    
    plt.subplot(1, 3, 2)
    plt.imshow(mask_data_ne[10], cmap='Blues')
    plt.title(f"NE Mask - Slice 10")
    plt.axis('off')
    
    plt.subplot(1, 3, 3)
    plt.imshow(mask_data_pm[10] + mask_data_ne[10], cmap='gray')
    plt.title(f"Mask Data - Data ID: {dataid}")
    plt.axis('off')
    
    plt.tight_layout()
    plt.show()

    # 2. Initialize partitioning processor
    logger.step("Initializing partitioning processor")
    partitioner = Partitioning(root_dir)

    # 3. Extract boundaries for analysis
    logger.step("Extracting NE and PM boundaries")
    center, ne_edge, pm_edge = partitioner.extract_ne_pm_edges(mask_data_pm, mask_data_ne)
    
    logger.step(f"Center shape: {center.shape if hasattr(center, 'shape') else type(center)}")
    logger.step(f"NE edge shape: {ne_edge.shape if hasattr(ne_edge, 'shape') else type(ne_edge)}")
    logger.step(f"PM edge shape: {pm_edge.shape if hasattr(pm_edge, 'shape') else type(pm_edge)}")
    print(f"Center shape: {center.shape if hasattr(center, 'shape') else type(center)}")
    print(f"NE edge shape: {ne_edge.shape if hasattr(ne_edge, 'shape') else type(ne_edge)}")
    print(f"PM edge shape: {pm_edge.shape if hasattr(pm_edge, 'shape') else type(pm_edge)}")
    print(f"Input mask shape: {mask_data_pm.shape}")
    
    # 4. Create radial partitions using NE-PM pairs method
    logger.step("Creating radial partitions using NE-PM pairs...")
    print("Creating radial partitions using NE-PM pairs...")
    partition_masks = partitioner.create_nepm_radial_partitions_pure_pairs(
        ne_edge, pm_edge, mask_data_pm.shape, 
        n_slices=8, 
        pm_mask=mask_data_pm, 
        ne_mask=mask_data_ne
    )
    
    logger.step(f"Partition masks created: {partition_masks.shape}, unique values: {np.unique(partition_masks)}")
    print(partition_masks.shape, np.unique(partition_masks))

    # Check partition distribution before saving
    unique_partitions = np.unique(partition_masks)
    logger.step("Partition voxel distribution:")
    for partition_id in unique_partitions:
        if partition_id != 0:
            count = np.sum(partition_masks == partition_id)
            logger.step(f"Partition {partition_id}: {count} voxels")
            print(f"Partition {partition_id}: {count} voxels")

    # Extract coordinate points from partition masks and save as XVG format
    logger.step("Extracting partition coordinates for XVG output...")
    print("Extracting partition coordinates for XVG output...")
    partition_coords = partitioner.extract_partition_coordinates(partition_masks, sampling_density=0.05)
    partitioner.save_partition_coords_to_xvg(partition_coords, dataid, results_dir)

    # 5. Visualize partitions
    logger.step("Generating visualizations...")
    print("Generating visualizations...")
    
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
    axes[1].set_title(f'Partition Masks (Slice {middle_slice})')
    axes[1].axis('off')
    
    # Display overlay image
    overlay = mask_data_pm[middle_slice] * 0.3 + partition_masks[middle_slice] * 0.7
    axes[2].imshow(overlay, cmap='viridis')
    axes[2].set_title('Overlay: Partitions on Cell')
    axes[2].axis('off')
    
    plt.tight_layout()
    plt.savefig(f"{results_dir}{dataid}_partition_comparison.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    # Get the number of unique partitions for proper vmax
    unique_partitions = np.unique(partition_masks)
    n_partitions = len(unique_partitions[unique_partitions != 0])  # Exclude background (0)
    
    logger.step(f"Number of partitions: {n_partitions}")
    logger.step(f"Unique partition values: {unique_partitions}")
    print(f"Number of partitions: {n_partitions}")
    print(f"Unique partition values: {unique_partitions}")
    
    # Check voxel count for each partition
    for partition_id in unique_partitions:
        if partition_id != 0:
            count = np.sum(partition_masks == partition_id)
            print(f"Partition {partition_id}: {count} voxels")
    
    # Save output files
    comparison_plot_path = f"{results_dir}{dataid}_partition_comparison.png"
    plt.tight_layout()
    plt.savefig(comparison_plot_path, dpi=300, bbox_inches='tight')
    logger.file_out(comparison_plot_path)
    plt.show()
    
    # 2D visualization at middle slice
    middle_slice = mask_data_pm.shape[0] // 2
    print(f"Visualizing slice {middle_slice} of {mask_data_pm.shape[0]} slices")
    
    # Check the slice data before visualization
    slice_data = partition_masks[middle_slice, :, :]
    print(f"Slice data shape: {slice_data.shape}")
    print(f"Slice data unique values: {np.unique(slice_data)}")
    
    # Convert partition_masks to list format expected by visualize_partitions
    unique_partitions = np.unique(partition_masks)
    partition_list = []
    for partition_id in unique_partitions:
        if partition_id != 0:  # Skip background
            mask = (partition_masks == partition_id)
            partition_list.append(mask)
    
    logger.step(f"Created {len(partition_list)} partition masks for visualization")
    
    # 2D visualization
    partition_2d_path = f"{results_dir}{dataid}_partitions_2d.png"
    visualize_partitions(partition_list, slice_idx=middle_slice, 
                        save_path=partition_2d_path)
    logger.file_out(partition_2d_path)

    logger.step("=" * 60)
    logger.step("Partitioning analysis completed successfully!")
    logger.step(f"Generated {len(partition_list)} partitions for dataset: {dataid}")
    logger.step("=" * 60)


if __name__ == '__main__':
    main()

