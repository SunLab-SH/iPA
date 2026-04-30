"""
SIM Mitochondria Segmentation - Prediction Example

This script demonstrates mitochondria segmentation on fluorescence microscopy data
using a thresholding-based pipeline inspired by Lefebvre et al. (2021).

Pipeline steps:
1. Background correction (Gaussian blur + subtraction)
2. Global intensity thresholding (Otsu's method)
3. Morphological filtering (opening/closing for denoising)
4. Connected-component labeling with size filtering

Data source: Synthetic data (for demonstration only)
Note: This analysis is not used in the main paper; provided for API completeness.

Usage:
    python demo_SIM_Mito_predict.py
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add iPA module path
current_dir = os.path.dirname(os.path.abspath(__file__))
ipa_root = os.path.dirname(os.path.dirname(current_dir))
sys.path.insert(0, ipa_root)

from ipa.processing.segmentation import create_segmenter


def generate_synthetic_mito_data(shape=(512, 512), num_mitos=30):
    """Generate synthetic mitochondria-like structures."""
    print("Generating synthetic mitochondria-like data...")
    
    # Create background
    image = np.random.normal(0.1, 0.05, shape).astype(np.float32)
    gt_mask = np.zeros(shape, dtype=np.uint8)
    
    # Create elongated mito structures
    np.random.seed(42)
    
    for _ in range(num_mitos):
        # Random starting point
        y_start = np.random.randint(50, shape[0] - 50)
        x_start = np.random.randint(50, shape[1] - 50)
        
        # Create curved elongated structure
        length = np.random.randint(30, 80)
        angle = np.random.uniform(0, 2 * np.pi)
        
        y_coords = []
        x_coords = []
        
        for i in range(length):
            angle += np.random.normal(0, 0.15)
            y = int(y_start + i * np.cos(angle))
            x = int(x_start + i * np.sin(angle))
            
            if 0 <= y < shape[0] and 0 <= x < shape[1]:
                y_coords.append(y)
                x_coords.append(x)
        
        # Draw mito (with thickness)
        if len(y_coords) > 0:
            for dy in range(-2, 3):
                for dx in range(-2, 3):
                    for y, x in zip(y_coords, x_coords):
                        yy, xx = y + dy, x + dx
                        if 0 <= yy < shape[0] and 0 <= xx < shape[1]:
                            gt_mask[yy, xx] = 1
                            image[yy, xx] = 0.6 + np.random.normal(0, 0.1)
    
    print(f"Generated {num_mitos} synthetic mitochondria")
    return image, gt_mask


def main():
    """Main prediction pipeline using synthetic data."""
    
    print("="*70)
    print("Mitochondria Segmentation - API Demo (Synthetic Data)")
    print("="*70)
    
    # Setup output directory
    output_dir = Path('sim_mito_output')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Generate synthetic data
    print("\n[Step 1] Generating Synthetic Mito Data")
    print("-" * 70)
    print("Note: This demo uses synthetic data for API demonstration only.")
    print("      Mito segmentation is not used in the main analysis.")
    
    image, gt_mask = generate_synthetic_mito_data(shape=(512, 512), num_mitos=30)
    
    # Save input
    np.save(output_dir / 'input_image.npy', image)
    print(f"Input saved to: {output_dir / 'input_image.npy'}")
    
    # Step 2: Create segmenter
    print("\n[Step 2] Creating Segmenter")
    print("-" * 70)
    print("Using unified API: create_segmenter('sim', 'mito')")
    segmenter = create_segmenter(modality='sim', task='mito')
    print(f"Segmenter info: {segmenter.get_info()}")
    print("Note: Mito segmentation uses thresholding-based approach (no deep learning model)")
    
    # Step 3: Run prediction
    print("\n[Step 3] Running Segmentation")
    print("-" * 70)
    
    try:
        mito_mask = segmenter.predict(image, threshold=None, min_size=100)
        print(f"✅ Prediction successful!")
        print(f"   Mask shape: {mito_mask.shape}")
        print(f"   Mito voxels: {mito_mask.sum()}")
        
        # Save prediction
        np.save(output_dir / 'mito_mask.npy', mito_mask)
        print(f"   Saved to: {output_dir / 'mito_mask.npy'}")
    except Exception as e:
        print(f"❌ Prediction failed: {e}")
        print("   This may require implementing the thresholding-based mito segmentation algorithm")
        return
    
    # Step 4: Visualize
    print("\n[Step 4] Visualizing Results")
    print("-" * 70)
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    axes[0].imshow(image, cmap='gray')
    axes[0].set_title('Input Image')
    axes[0].axis('off')
    
    axes[1].imshow(mito_mask, cmap='jet')
    axes[1].set_title('Mito Segmentation')
    axes[1].axis('off')
    
    overlay = np.stack([image] * 3, axis=-1) if image.ndim == 2 else image[..., :3]
    overlay[mito_mask > 0] = [1, 0, 0]
    axes[2].imshow(overlay)
    axes[2].set_title('Overlay (Red=Mito)')
    axes[2].axis('off')
    
    plt.tight_layout()
    viz_path = output_dir / 'results.png'
    plt.savefig(viz_path, dpi=150, bbox_inches='tight')
    print(f"Visualization saved to: {viz_path}")
    plt.close()
    
    # Summary
    print("\n" + "="*70)
    print("Summary")
    print("="*70)
    print(f"\n✅ Mito segmentation completed!")
    print(f"\nOutput files:")
    print(f"  - Input image:      {output_dir / 'input_image.npy'}")
    print(f"  - Mito mask:        {output_dir / 'mito_mask.npy'}")
    print(f"  - Visualization:    {viz_path}")
    
    print(f"\n⚠️  Note: This demo uses synthetic data for API demonstration.")
    print(f"   Mito segmentation is not used in the main analysis.")
    print(f"   Pipeline follows Lefebvre et al. (2021) thresholding-based approach.")
    
    print(f"\nAlgorithm: Thresholding-based pipeline")
    print(f"  - Background correction: Gaussian blur (sigma=20) + subtraction")
    print(f"  - Thresholding: Otsu's method (automatic)")
    print(f"  - Morphology: Opening (disk/ball r=1) + Closing (r=2)")
    print(f"  - Size filtering: min_size=100 pixels")
    
    print("\n" + "="*70)
    print("Demo completed!")
    print("="*70)


if __name__ == "__main__":
    main()
