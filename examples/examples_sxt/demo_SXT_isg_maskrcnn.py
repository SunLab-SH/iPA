"""
SXT ISG Mask R-CNN Instance Segmentation Demo

This script demonstrates instance segmentation using the trained Mask R-CNN model.
Note: This is an alternative to the blob_fit post-processing method.
For production use, we recommend demo_SXT_isg_instance.py (blob_fit method).

Model Performance:
- mAP @ [0.50:0.95]: 0.2276
- mAP @ 0.50: 0.6137
- Recall @ 100: 0.367
"""

import os
import sys
import numpy as np
import mrcfile
import tifffile
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from ipa.processing.segmentation.unified import create_segmenter


def main():
    # Configuration
    TEST_CELL_ID = "784_5"
    IMG_DIR = "data/sxt/sxt_isg_training/for_24_datasets"
    
    # Find the mrc file
    img_path = None
    for f in os.listdir(IMG_DIR):
        if f.endswith('.mrc') and TEST_CELL_ID in f:
            img_path = os.path.join(IMG_DIR, f)
            break
    
    if not img_path:
        print(f"Error: Could not find mrc file for cell {TEST_CELL_ID}")
        return
    
    print("=" * 60)
    print("SXT ISG Mask R-CNN Instance Segmentation Demo")
    print("=" * 60)
    print(f"Test Cell: {TEST_CELL_ID}")
    print(f"Image Path: {img_path}")
    print()
    print("Note: This uses Mask R-CNN (mAP=0.23)")
    print("For better results, use demo_SXT_isg_instance.py (blob_fit)")
    print()
    
    # 1. Load data
    print("[1/4] Loading 3D volume...")
    with mrcfile.open(img_path, permissive=True) as mrc:
        vol = mrc.data.astype(np.float32)
    print(f"  Volume shape: {vol.shape}")
    
    # 2. Create segmenter and load model
    print("\n[2/4] Loading Mask R-CNN model...")
    segmenter = create_segmenter(modality='sxt', task='isg_maskrcnn')
    segmenter.load_model()
    print("  Model loaded successfully!")
    
    # 3. Predict
    print("\n[3/4] Running Mask R-CNN prediction...")
    instance_mask = segmenter.predict(vol, threshold=0.5)
    
    # 4. Statistics
    print("\n[4/4] Analyzing instances...")
    unique_ids = np.unique(instance_mask)
    num_instances = len(unique_ids[unique_ids > 0])
    
    print(f"  Total instances detected: {num_instances}")
    
    if num_instances > 0:
        volumes = []
        for uid in unique_ids:
            if uid > 0:
                volumes.append(np.sum(instance_mask == uid))
        
        volumes = np.array(volumes)
        print(f"  Average volume: {np.mean(volumes):.2f} voxels")
        print(f"  Median volume:  {np.median(volumes):.2f} voxels")
        print(f"  Max volume:     {np.max(volumes)} voxels")
        print(f"  Min volume:     {np.min(volumes)} voxels")
    
    # 5. Visualization
    print("\nGenerating visualization...")
    z_idx = vol.shape[0] // 2
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # Original
    axes[0].imshow(vol[z_idx], cmap='gray')
    axes[0].set_title("Original Slice", fontsize=12)
    axes[0].axis('off')
    
    # Binary Mask
    binary_mask = (instance_mask > 0).astype(np.uint8)
    axes[1].imshow(binary_mask[z_idx], cmap='viridis')
    axes[1].set_title("Binary Mask", fontsize=12)
    axes[1].axis('off')
    
    # Instance Mask (Labeled)
    if num_instances > 0:
        # Create a random color map for instances
        np.random.seed(42)
        colors = np.random.rand(num_instances, 3)
        # Add black for background
        colors = np.vstack([[0, 0, 0], colors])
        cmap = ListedColormap(colors)
        
        axes[2].imshow(instance_mask[z_idx], cmap=cmap)
        axes[2].set_title(f"Instance Mask ({num_instances} ISGs)", fontsize=12)
    else:
        axes[2].text(0.5, 0.5, 'No instances found', ha='center', va='center', transform=axes[2].transAxes)
        axes[2].set_title("Instance Mask", fontsize=12)
    axes[2].axis('off')
    
    plt.tight_layout()
    output_path = f"sxt_isg_maskrcnn_result_{TEST_CELL_ID}.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Visualization saved to: {output_path}")
    
    # Save the instance mask
    mask_output = f"sxt_isg_maskrcnn_{TEST_CELL_ID}.tiff"
    tifffile.imwrite(mask_output, instance_mask)
    print(f"Instance mask saved to: {mask_output}")
    
    print("\n" + "=" * 60)
    print("Demo completed!")
    print("=" * 60)
    print("\nComparison:")
    print("- Mask R-CNN: Fast inference, mAP=0.23")
    print("- Blob Fit (demo_SXT_isg_instance.py): More accurate, professional algorithm")
    print("=" * 60)


if __name__ == "__main__":
    main()
