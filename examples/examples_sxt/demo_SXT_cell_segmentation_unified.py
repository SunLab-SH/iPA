"""
Demo: SXT Cell Segmentation using Unified API

This example demonstrates the new unified segmentation API.
Note: This is a simplified demo with synthetic data.
For real usage, you need trained models and actual SXT data.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Add iPA module path
current_dir = os.path.dirname(os.path.abspath(__file__))
ipa_root = os.path.dirname(os.path.dirname(current_dir))
sys.path.insert(0, ipa_root)

from ipa.processing.segmentation import create_segmenter


def generate_synthetic_sxt_data(shape=(64, 256, 256)):
    """
    Generate synthetic SXT-like 3D data for demonstration.
    
    Args:
        shape: Volume shape (depth, height, width)
        
    Returns:
        tuple: (volume, cell_mask, nucleus_mask)
    """
    print("Generating synthetic SXT data...")
    
    z, y, x = np.meshgrid(
        np.arange(shape[0]),
        np.arange(shape[1]),
        np.arange(shape[2]),
        indexing='ij'
    )
    
    # Create cell mask (ellipsoid)
    center_cell = (shape[0]//2, shape[1]//2, shape[2]//2)
    radii_cell = (shape[0]//2 - 2, shape[1]//2 - 10, shape[2]//2 - 10)
    
    cell_mask = ((z - center_cell[0])**2 / radii_cell[0]**2 +
                 (y - center_cell[1])**2 / radii_cell[1]**2 +
                 (x - center_cell[2])**2 / radii_cell[2]**2) <= 1
    
    # Create nucleus mask (smaller ellipsoid inside cell)
    center_nuc = center_cell
    radii_nuc = (shape[0]//4, shape[1]//4 - 5, shape[2]//4 - 5)
    
    nucleus_mask = ((z - center_nuc[0])**2 / radii_nuc[0]**2 +
                    (y - center_nuc[1])**2 / radii_nuc[1]**2 +
                    (x - center_nuc[2])**2 / radii_nuc[2]**2) <= 1
    
    # Create intensity volume (cell brighter than background)
    volume = np.zeros(shape, dtype=np.float32)
    volume[cell_mask] = 0.7 + np.random.normal(0, 0.1, cell_mask.sum())
    volume[nucleus_mask] = 0.9 + np.random.normal(0, 0.05, nucleus_mask.sum())
    
    # Add noise
    volume += np.random.normal(0, 0.05, shape)
    volume = np.clip(volume, 0, 1)
    
    print(f"Volume shape: {volume.shape}")
    print(f"Cell mask shape: {cell_mask.shape}")
    print(f"Nucleus mask shape: {nucleus_mask.shape}")
    
    return volume.astype(np.float32), cell_mask.astype(np.uint8), nucleus_mask.astype(np.uint8)


def visualize_results(volume, pred_mask, gt_mask, slice_idx=None):
    """
    Visualize segmentation results.
    
    Args:
        volume: Input volume
        pred_mask: Predicted mask
        gt_mask: Ground truth mask
        slice_idx: Slice index to visualize
    """
    if slice_idx is None:
        slice_idx = volume.shape[0] // 2
    
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    
    # Volume
    axes[0].imshow(volume[slice_idx], cmap='gray')
    axes[0].set_title('Input Volume')
    axes[0].axis('off')
    
    # Ground truth
    axes[1].imshow(gt_mask[slice_idx], cmap='jet')
    axes[1].set_title('Ground Truth')
    axes[1].axis('off')
    
    # Prediction
    axes[2].imshow(pred_mask[slice_idx], cmap='jet')
    axes[2].set_title('Prediction')
    axes[2].axis('off')
    
    # Overlay
    overlay = volume[slice_idx].copy()
    overlay = np.stack([overlay] * 3, axis=-1)  # RGB
    overlay[pred_mask[slice_idx] > 0] = [1, 0, 0]  # Red for prediction
    axes[3].imshow(overlay)
    axes[3].set_title('Overlay (Red=Prediction)')
    axes[3].axis('off')
    
    plt.tight_layout()
    plt.savefig('sxt_segmentation_demo.png', dpi=150, bbox_inches='tight')
    print("Results saved to: sxt_segmentation_demo.png")
    plt.show()


def main():
    """Main demonstration."""
    
    print("=" * 70)
    print("SXT Cell Segmentation Demo (Unified API)")
    print("=" * 70)
    
    # 1. Generate synthetic data
    print("\n[Step 1] Generating synthetic data...")
    volume, cell_mask_gt, nucleus_mask_gt = generate_synthetic_sxt_data(
        shape=(32, 128, 128)
    )
    
    # 2. Create segmenter using unified API
    print("\n[Step 2] Creating segmenter...")
    segmenter = create_segmenter(modality='sxt', task='cell')
    print(f"Segmenter info: {segmenter.get_info()}")
    
    # 3. Note about model loading
    print("\n[Step 3] Model Loading")
    print("Note: For real usage, you need to:")
    print("  1. Train a model or download pretrained weights")
    print("  2. Call: segmenter.load_model('path/to/model.pth')")
    print("\nFor this demo, we'll show the API structure.")
    print("The actual segmentation requires trained models.")
    
    # 4. Show how to use (commented out since we don't have real model)
    print("\n[Step 4] Usage Example (pseudo-code)")
    print("""
# Load model (when you have one)
segmenter.load_model('sxt_cell_model.pth')

# Predict
result = segmenter.predict(volume)

# Access masks
cell_mask = result['cell_mask']
nucleus_mask = result['nucleus_mask']
    """)
    
    # 5. Show legacy API still works
    print("\n[Step 5] Backward Compatibility")
    print("Legacy API still works:")
    print("  from ipa.processing.segmentation import run_cell_segmentation")
    print("  result = run_cell_segmentation(dataid_list=['784_5'])")
    
    # 6. Visualize synthetic data
    print("\n[Step 6] Visualizing synthetic data...")
    visualize_results(
        volume=volume,
        pred_mask=cell_mask_gt,  # Using GT as "prediction" for demo
        gt_mask=cell_mask_gt,
        slice_idx=volume.shape[0] // 2
    )
    
    print("\n" + "=" * 70)
    print("Demo completed!")
    print("=" * 70)
    print("\nNext steps:")
    print("1. Prepare your SXT data")
    print("2. Train or obtain a segmentation model")
    print("3. Use the unified API for prediction")
    print("\nSee documentations/SEGMENTATION_QUICK_REFERENCE.md for more details")


if __name__ == "__main__":
    main()
