"""
Demo: N2N Denoising Prediction on Real SXT Data

This example demonstrates how to use a trained N2N model to denoise real SXT data.
"""

import os
import sys
import numpy as np
import mrcfile
import matplotlib.pyplot as plt

# Add iPA module path
current_dir = os.path.dirname(os.path.abspath(__file__))
ipa_root = os.path.dirname(os.path.dirname(current_dir))
sys.path.insert(0, ipa_root)

from ipa.processing.denoising import N2N


def main():
    print("=" * 70)
    print("N2N Denoising Prediction (Real SXT Data)")
    print("=" * 70)

    # 1. Load Real Data
    data_path = os.path.join(ipa_root, 'data', 'sxt', 'Stevens_pancreatic_INS_1E_784_5_pre_rec.mrc')
    if not os.path.exists(data_path):
        print(f"Error: Real SXT data not found at {data_path}")
        return

    print(f"\n[Step 1] Loading data from {os.path.basename(data_path)}...")
    with mrcfile.open(data_path, permissive=True) as mrc:
        clean_data = mrc.data.astype(np.float32)
    
    # Normalize to [0, 1]
    min_val, max_val = clean_data.min(), clean_data.max()
    clean_data = (clean_data - min_val) / (max_val - min_val)
    
    # Use a center crop
    z, y, x = clean_data.shape
    cz, cy, cx = z//2, y//2, x//2
    crop_size = 64
    clean_crop = clean_data[cz-crop_size:cz+crop_size, cy-crop_size:cy+crop_size, cx-crop_size:cx+crop_size]
    
    # Add noise to simulate input
    noisy_data = np.clip(clean_crop + np.random.normal(0, 0.15, clean_crop.shape), 0, 1)

    # 2. Load Model
    # Use the trained model from the unified models directory
    model_dir = os.path.join(PROJECT_ROOT, 'ipa', 'processing', 'denoising', 'models', 'n2n', 'gaussian-0012')
    if not os.path.exists(model_dir):
        print(f"Error: Trained model directory not found at {model_dir}. Please run demo_SXT_n2n_train.py first.")
        return

    # Find the best model (lowest loss)
    import glob
    model_files = glob.glob(os.path.join(model_dir, '*.pt'))
    if not model_files:
        print("Error: No .pt model files found in the directory.")
        return
    
    # Sort by filename to get the last epoch or use a specific one
    model_path = sorted(model_files)[-1] 

    print(f"\n[Step 2] Loading model from {model_path}...")
    n2n = N2N(n_channels=1, n_filters=32)
    n2n.load_model(model_path)

    # 3. Predict
    print("\n[Step 3] Running prediction...")
    denoised_data = n2n.predict(noisy_data, batch_size=4)

    # 4. Visualize
    print("\n[Step 4] Visualizing results...")
    slice_idx = noisy_data.shape[0] // 2
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    axes[0].imshow(noisy_data[slice_idx], cmap='gray'); axes[0].set_title('Noisy Input')
    axes[1].imshow(denoised_data[slice_idx], cmap='gray'); axes[1].set_title('Denoised (N2N)')
    axes[2].imshow(np.abs(noisy_data[slice_idx] - denoised_data[slice_idx]), cmap='hot'); axes[2].set_title('Difference')
    for ax in axes: ax.axis('off')
    
    out_path = os.path.join(ipa_root, 'results', 'demo_n2n_predict_sxt.png')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"Results saved to: {out_path}")
    plt.show()

    print("\n" + "=" * 70)
    print("Prediction completed successfully!")
    print("=" * 70)


if __name__ == "__main__":
    main()
