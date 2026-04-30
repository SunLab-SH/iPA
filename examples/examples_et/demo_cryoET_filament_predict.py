"""
Cryo-ET Filament Segmentation - Prediction Example

This script demonstrates filament segmentation and skeletonization on Cryo-ET data
using traditional image processing methods (no deep learning model required).

Pipeline:
1. Threshold-based segmentation
2. Morphological cleaning
3. 3D skeletonization to extract filament structure

Data source: Real Cryo-ET denoised volume (20210517_30_009_denoised.mrc)
Note: This analysis is not used in the main paper; provided for API completeness.

Usage:
    python demo_cryoET_filament_predict.py
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add project root to path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
sys.path.insert(0, PROJECT_ROOT)

from ipa.processing.segmentation import create_segmenter


def load_real_cryoet_data():
    """Load real Cryo-ET denoised volume from MRC file.
    
    Note: Uses compressed center subvolume (1/8 of original) for faster processing.
    Full volume available as: 20210517_30_009_denoised.mrc.bak
    """
    import mrcfile
    
    data_dir = os.path.join(PROJECT_ROOT, 'data', 'cryoET')
    # Use compressed center subvolume for demo (100 MB vs 2.4 GB)
    mrc_path = Path(data_dir) / '20210517_30_009_denoised_center_subvolume.mrc'
    
    if not mrc_path.exists():
        raise FileNotFoundError(f"Cryo-ET data not found: {mrc_path}")
    
    print(f"Loading Cryo-ET volume from: {mrc_path}")
    print(f"Note: Using compressed center subvolume (1/8 of original, 95.8% compression)")
    with mrcfile.open(mrc_path, mode='r') as mrc:
        volume = mrc.data.astype(np.float32)
    
    print(f"Volume shape: {volume.shape}, dtype: {volume.dtype}")
    print(f"Value range: [{volume.min():.2f}, {volume.max():.2f}]")
    
    return volume


def visualize_3d_results(volume, mask, skeleton, output_path: Path):
    """Visualize 3D segmentation results with multiple slices."""
    slice_z = volume.shape[0] // 2
    slice_y = volume.shape[1] // 2
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # Z-slice
    axes[0, 0].imshow(volume[slice_z], cmap='gray')
    axes[0, 0].set_title(f'Input Volume (Z={slice_z})')
    axes[0, 0].axis('off')
    
    axes[0, 1].imshow(mask[slice_z], cmap='jet')
    axes[0, 1].set_title('Segmentation Mask')
    axes[0, 1].axis('off')
    
    axes[0, 2].imshow(skeleton[slice_z], cmap='Reds')
    axes[0, 2].set_title('Skeleton')
    axes[0, 2].axis('off')
    
    # Y-slice
    axes[1, 0].imshow(volume[:, slice_y, :], cmap='gray')
    axes[1, 0].set_title(f'Input Volume (Y={slice_y})')
    axes[1, 0].axis('off')
    
    axes[1, 1].imshow(mask[:, slice_y, :], cmap='jet')
    axes[1, 1].set_title('Filament Mask')
    axes[1, 1].axis('off')
    
    axes[1, 2].imshow(skeleton[:, slice_y, :], cmap='Reds')
    axes[1, 2].set_title('Skeleton')
    axes[1, 2].axis('off')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Visualization saved to: {output_path}")
    plt.close()


def main():
    """Main prediction pipeline using real Cryo-ET data."""
    
    print("="*70)
    print("Cryo-ET Filament Segmentation - Real Data")
    print("="*70)
    
    # Setup output directory
    output_dir = Path(PROJECT_ROOT) / 'results' / 'et_filament_prediction_demo'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Load real Cryo-ET data
    print("\n[Step 1] Loading Real Cryo-ET Data")
    print("-" * 70)
    
    volume_full = load_real_cryoet_data()
    
    # Extract center 1/2 of each dimension (volume becomes 1/8)
    # ⚠️  NOTE: This sub-volume extraction is for DEMONSTRATION ONLY to reduce processing time.
    #          For actual analysis, use the full volume or adjust region based on your needs.
    print(f"\nExtracting center sub-volume (1/2 of each dimension)...")
    print(f"⚠️  Note: Sub-volume extraction is for demonstration only. Use full volume for real analysis.")
    z_size, y_size, x_size = volume_full.shape
    z_start, z_end = z_size // 4, 3 * z_size // 4
    y_start, y_end = y_size // 4, 3 * y_size // 4
    x_start, x_end = x_size // 4, 3 * x_size // 4
    
    volume = volume_full[z_start:z_end, y_start:y_end, x_start:x_end].copy()
    print(f"Original shape: {volume_full.shape}")
    print(f"Sub-volume shape: {volume.shape} (1/8 of original)")
    print(f"Region: Z[{z_start}:{z_end}], Y[{y_start}:{y_end}], X[{x_start}:{x_end}]")
    
    # Save input
    np.save(output_dir / 'input_volume.npy', volume)
    print(f"Input saved to: {output_dir / 'input_volume.npy'}")
    
    # Step 2: Create segmenter
    print("\n[Step 2] Creating Segmenter")
    print("-" * 70)
    print("Using unified API: create_segmenter('et', 'filament')")
    segmenter = create_segmenter(modality='et', task='filament')
    print(f"Segmenter info: {segmenter.get_info()}")
    print("Note: Filament segmentation uses traditional methods (no model needed)")
    
    # Step 3: Run prediction with skeletonization
    print("\n[Step 3] Running Segmentation and Skeletonization")
    print("-" * 70)
    
    try:
        # Run segmentation (returns dict with 'mask' and 'skeleton')
        result = segmenter.predict(volume, skeletonize=True, threshold_multiplier=0.5)
        mask = result['mask']
        skeleton = result['skeleton']
        
        print(f"✅ Prediction successful!")
        print(f"   Mask shape: {mask.shape}")
        print(f"   Skeleton shape: {skeleton.shape}")
        print(f"   Filament voxels: {mask.sum()}")
        print(f"   Skeleton voxels: {skeleton.sum()}")
        
        # Save results
        np.save(output_dir / 'filament_mask.npy', mask)
        np.save(output_dir / 'filament_skeleton.npy', skeleton)
        print(f"   Saved to: {output_dir}")
    except Exception as e:
        print(f"❌ Prediction failed: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # Step 4: Visualize
    print("\n[Step 4] Visualizing Results")
    print("-" * 70)
    visualize_3d_results(volume, mask, skeleton, output_dir / 'results.png')
    
    # Summary
    print("\n" + "="*70)
    print("Summary")
    print("="*70)
    print(f"\n✅ Cryo-ET filament segmentation completed!")
    print(f"\nOutput files:")
    print(f"  - Input volume:       {output_dir / 'input_volume.npy'}")
    print(f"  - Filament mask:      {output_dir / 'filament_mask.npy'}")
    print(f"  - Filament skeleton:  {output_dir / 'filament_skeleton.npy'}")
    print(f"  - Visualization:      {output_dir / 'results.png'}")
    
    print(f"\nPipeline: Thresholding → Morphology → 3D Skeletonization")
    print(f"Note: Cryo-ET filament segmentation is not used in the main analysis.")
    
    print("\n" + "="*70)
    print("Demo completed!")
    print("="*70)


if __name__ == "__main__":
    main()
