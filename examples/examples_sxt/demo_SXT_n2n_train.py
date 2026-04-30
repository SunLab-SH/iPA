"""
Demo: N2N Denoising Training on Real SXT Data

This example demonstrates how to train a Noise2Noise (N2N) model using real SXT data.
It creates noisy pairs by adding Gaussian noise to the original volume.
"""

import os
import sys
import numpy as np
import mrcfile

# Add iPA module path
current_dir = os.path.dirname(os.path.abspath(__file__))
ipa_root = os.path.dirname(os.path.dirname(current_dir))
sys.path.insert(0, ipa_root)

from ipa.processing.denoising import N2N


def main():
    print("=" * 70)
    print("N2N Denoising Training (Real SXT Data)")
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
    
    # Use a center crop for faster demo training
    z, y, x = clean_data.shape
    cz, cy, cx = z//2, y//2, x//2
    crop_size = 64
    clean_crop = clean_data[cz-crop_size:cz+crop_size, cy-crop_size:cy+crop_size, cx-crop_size:cx+crop_size]

    # Create Noisy Pairs (N2N requirement)
    print("\n[Step 2] Generating noisy pairs...")
    noisy_1 = np.clip(clean_crop + np.random.normal(0, 0.15, clean_crop.shape), 0, 1)
    noisy_2 = np.clip(clean_crop + np.random.normal(0, 0.15, clean_crop.shape), 0, 1)
    
    # Split into train/val
    split_idx = int(clean_crop.shape[0] * 0.8)
    train_n1, val_n1 = noisy_1[:split_idx], noisy_1[split_idx:]
    train_n2, val_n2 = noisy_2[:split_idx], noisy_2[split_idx:]

    # 2. Initialize and Train
    print("\n[Step 3] Training N2N model...")
    n2n = N2N(n_channels=1, n_filters=32)
    
    n2n.train(
        noisy_data_1=train_n1,
        noisy_data_2=train_n2,
        val_data_1=val_n1,
        val_data_2=val_n2,
        epochs=20,
        batch_size=4,
        lr=1e-3
    )

    # 3. Save Model
    model_path = os.path.join(current_dir, 'n2n_model_demo.pth')
    n2n.save_model(model_path)
    print(f"\nModel saved to: {model_path}")

    print("\n" + "=" * 70)
    print("Training completed successfully!")
    print("=" * 70)


if __name__ == "__main__":
    main()
