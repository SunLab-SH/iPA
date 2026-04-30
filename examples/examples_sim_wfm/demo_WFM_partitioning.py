#%%

# examples/basic_usage/demo_partitioning.py
import numpy as np
from ipa.processing.partitioning import Partitioning
from ipa.processing.partitioning import visualize_partitions, plot_partition_features
from ipa.data_loader import QuickLogger
import matplotlib
matplotlib.use('Agg')  # Use Agg backend for non-interactive mode
import matplotlib.pyplot as plt
import mrcfile
import pickle
import json
import os
import sys
from pathlib import Path

# --- Path Configuration ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
sys.path.insert(0, PROJECT_ROOT)

def main():
    mainpath = PROJECT_ROOT
    
    # Initialize logger
    log_dir = f'{mainpath}/logs'
    logger = QuickLogger("wfm_partitioning_analysis", log_dir=log_dir)
    logger.step("Starting WFM Partitioning Analysis Demo")
    
    print(mainpath)
    root_dir = f'{mainpath}/data/wfm/'  # Output directory, recommended to create new
    dataid = '2.8-30_P3S3-1'
    
    logger.step(f"Working with dataset: {dataid}")
    
    # Create results directory
    results_dir = f'{root_dir}results/'
    os.makedirs(results_dir, exist_ok=True)

    mask_pm_path = f'{mainpath}/data/wfm/{dataid}_Memb.mrc'
    mask_ne_path = f'{mainpath}/data/wfm/{dataid}_Nuc.mrc'

    logger.step("Loading mask data")
    logger.file_in(mask_pm_path)
    logger.file_in(mask_ne_path)

    mask_data_pm = mrcfile.open(mask_pm_path, mode='r', permissive=True).data
    mask_data_ne = mrcfile.open(mask_ne_path, mode='r', permissive=True).data

    logger.step(f"PM mask loaded: {mask_data_pm.shape}")
    logger.step(f"NE mask loaded: {mask_data_ne.shape}")
    print(mask_data_pm.shape, mask_data_ne.shape)

    # # 1. Load mask data
    # mask_data = np.load(mask_path)  # Shape: (T, Z, 2, X, Y)
    
    # mask_data_pm = mask_data[tp, :, 0, :, :]  # Extract mask data from frame 0
    # mask_data_ne = mask_data[tp, :, 1, :, :]  # Extract mask data from frame 0
    # mask_data_pm = np.transpose(mask_data_pm, (0, 1, 2))
    # mask_data_ne = np.transpose(mask_data_ne, (0, 1, 2))
    tp = 0  # Process only frame 0, use loop for multiple frames


    plt.imshow(mask_data_pm[15] + mask_data_ne[15], cmap='gray')
    plt.title(f"Mask Data - Data ID: {dataid}, Time Point: {tp}")
    plt.axis('off')
    mask_viz_path = f'{results_dir}mask_visualization.png'
    plt.savefig(mask_viz_path, dpi=150, bbox_inches='tight')
    logger.step(f"Mask visualization saved to: {mask_viz_path}")
    plt.close()

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
    
    # 4. Use create_radial_partitions to generate continuous partition masks, pass original masks
    logger.step("Creating continuous radial partitions using distance transform...")
    print("Creating continuous radial partitions using distance transform...")
    partition_masks = partitioner.create_nepm_radial_partitions(
        ne_edge, pm_edge, mask_data_pm.shape, 
        n_slices=8, 
        pm_mask=mask_data_pm, 
        ne_mask=mask_data_ne
    )
    
    logger.step(f"Partition masks created: {partition_masks.shape}, unique values: {np.unique(partition_masks)}")
    print(partition_masks.shape, np.unique(partition_masks))

    # No need to save the full partition masks - XVG coordinates are sufficient for RDF analysis
    # partition_save_path = f"{results_dir}{dataid}_partition_masks.npy"
    # np.save(partition_save_path, partition_masks)
    # print(f"Partition masks saved to: {partition_save_path}")

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
    comparison_path = f"{results_dir}{dataid}_partition_comparison.png"
    plt.savefig(comparison_path, dpi=300, bbox_inches='tight')
    logger.step(f"Partition comparison saved to: {comparison_path}")
    plt.close()
    
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
            logger.step(f"Partition {partition_id}: {count} voxels")
            print(f"Partition {partition_id}: {count} voxels")
    
    # Save output files
    comparison_plot_path = f"{results_dir}{dataid}_partition_comparison.png"
    plt.tight_layout()
    plt.savefig(comparison_plot_path, dpi=300, bbox_inches='tight')
    logger.file_out(comparison_plot_path)
    logger.step(f"Final comparison plot saved to: {comparison_plot_path}")
    plt.close()
    
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
    
    print(f"Created {len(partition_list)} partition masks for visualization")
    
    # 2D visualization
    partition_2d_path = f"{results_dir}{dataid}_partitions_2d.png"
    visualize_partitions(partition_list, slice_idx=middle_slice, 
                        save_path=partition_2d_path)
    logger.file_out(partition_2d_path)
    logger.step("Partitioning analysis completed successfully")




if __name__ == '__main__':
    main()
