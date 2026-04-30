"""
Demo: SIM ISG Instance Segmentation

This example demonstrates how to perform instance segmentation for Insulin Secretory Granules (ISGs)
in SIM (Structured Illumination Microscopy) images using a lightweight intensity-based approach.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from skimage import measure, color

# Add iPA module path
current_dir = os.path.dirname(os.path.abspath(__file__))
ipa_root = os.path.dirname(os.path.dirname(current_dir))
sys.path.insert(0, ipa_root)

from ipa.processing.segmentation import create_segmenter


def main():
    print("=" * 70)
    print("SIM ISG Instance Segmentation Demo")
    print("=" * 70)

    # 1. Load Real SIM Data
    data_path = os.path.join(ipa_root, 'data', 'sim', '20220909_30-2-1-SIM_raw_ISG.tif')
    
    if not os.path.exists(data_path):
        print(f"Error: SIM data not found at {data_path}")
        return

    print(f"\n[Step 1] Loading data from {os.path.basename(data_path)}...")
    from tifffile import imread
    vol = imread(data_path).astype(np.float32)
    
    # Normalize
    vol = (vol - vol.min()) / (vol.max() - vol.min())
    print(f"Volume shape: {vol.shape}")

    # 2. Create Segmenter
    print("\n[Step 2] Creating SIM ISG segmenter...")
    seg = create_segmenter('sim', 'isg')

    # 3. Predict (Instance Segmentation)
    print("\n[Step 3] Running instance segmentation on the first slice...")
    # sigma controls the expected spot size, min_size filters out noise
    labeled_mask = seg.predict(vol[0], sigma=3.0, min_size=15)
    
    num_isgs = labeled_mask.max()
    print(f"Detected {num_isgs} individual ISGs.")

    # 4. Visualize Results
    print("\n[Step 4] Visualizing results...")
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # Original Image
    axes[0].imshow(vol[0], cmap='gray')
    axes[0].set_title('Original SIM Image')
    axes[0].axis('off')
    
    # Labeled Mask (Colored)
    colored_mask = color.label2rgb(labeled_mask, bg_label=0)
    axes[1].imshow(colored_mask)
    axes[1].set_title(f'Instance Masks ({num_isgs} ISGs)')
    axes[1].axis('off')
    
    # Overlay
    axes[2].imshow(vol[0], cmap='gray')
    axes[2].imshow(color.label2rgb(labeled_mask, bg_label=0, alpha=0.5))
    axes[2].set_title('Overlay')
    axes[2].axis('off')
    
    plt.tight_layout()
    out_path = os.path.join(ipa_root, 'results', 'demo_sim_isg_segmentation.png')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"Results saved to: {out_path}")
    plt.show()

    print("\n" + "=" * 70)
    print("Demo completed successfully!")
    print("=" * 70)


if __name__ == "__main__":
    main()
