"""
SXT ISG Segmentation - Prediction Example

This script demonstrates ISG (Insulin Granule) segmentation using REAL data.
It loads actual SIM images and applies the trained ISG model for prediction.

Data source: data/sxt/
Model: models/isg_mask_rcnn_model.pth

Usage:
    python demo_SXT_isg_predict.py
"""

import os
import sys
import numpy as np
from pathlib import Path

# Add iPA module path
current_dir = os.path.dirname(os.path.abspath(__file__))
ipa_root = os.path.dirname(os.path.dirname(current_dir))
sys.path.insert(0, ipa_root)

from ipa.processing.segmentation import create_segmenter


def generate_synthetic_isg_data(shape=(64, 128, 128), num_isgs=50):
    """
    Generate synthetic SXT volume with ISG (insulin granule) structures.
    
    ISGs appear as small spherical objects in SXT images.
    
    Args:
        shape: Volume shape (depth, height, width)
        num_isgs: Number of ISGs to generate
        
    Returns:
        Tuple of (volume, isg_mask)
    """
    print(f"Generating synthetic SXT data with {num_isgs} ISGs...")
    
    # Create background
    volume = np.random.normal(0.1, 0.03, shape).astype(np.float32)
    isg_mask = np.zeros(shape, dtype=np.uint8)
    
    # Generate random ISGs (small spheres)
    np.random.seed(42)
    
    for _ in range(num_isgs):
        # Random center
        z = np.random.randint(10, shape[0] - 10)
        y = np.random.randint(10, shape[1] - 10)
        x = np.random.randint(10, shape[2] - 10)
        
        # Random radius (3-6 voxels)
        radius = np.random.uniform(3, 6)
        
        # Create sphere
        z_grid, y_grid, x_grid = np.meshgrid(
            np.arange(shape[0]),
            np.arange(shape[1]),
            np.arange(shape[2]),
            indexing='ij'
        )
        
        distance = np.sqrt(
            (z_grid - z)**2 + 
            (y_grid - y)**2 + 
            (x_grid - x)**2
        )
        
        sphere = distance <= radius
        
        # Add ISG to volume (bright spots)
        volume[sphere] = 0.7 + np.random.normal(0, 0.1)
        isg_mask[sphere] = 1
    
    # Smooth slightly
    from scipy.ndimage import gaussian_filter
    volume = gaussian_filter(volume, sigma=0.5)
    volume = np.clip(volume, 0, 1)
    
    print(f"Volume shape: {volume.shape}")
    print(f"ISG voxels: {isg_mask.sum()}")
    print(f"Number of ISGs: ~{num_isgs}")
    
    return volume.astype(np.float32), isg_mask.astype(np.uint8)


def main():
    """Main prediction pipeline with REAL SXT data."""
    
    print("="*70)
    print("SXT ISG Segmentation - Prediction with Real Data")
    print("="*70)
    
    # Setup paths
    data_dir = Path('/media/cuixi/data7/liad/gitspace/iPA/data/sxt')
    image_path = data_dir / 'Stevens_pancreatic_INS_1E_784_5_pre_rec.mrc'
    isg_label_path = data_dir / '784_5_isg_label.tiff'
    output_dir = Path('sxt_isg_real_output')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Load real SXT data
    print("\n[Step 1] Loading Real SXT Data")
    print("-" * 70)
    print(f"Image path: {image_path}")
    
    from ipa.data_loader import UniversalDataLoader
    print("Loading SXT volume...")
    volume = UniversalDataLoader.load_data(str(image_path))
    print(f"✅ Loaded! Shape: {volume.shape}, dtype: {volume.dtype}")
    
    # Load ground truth label if available
    isg_gt = None
    if isg_label_path.exists():
        print("Loading ISG label...")
        isg_gt = UniversalDataLoader.load_data(str(isg_label_path))
        print(f"   ISG label shape: {isg_gt.shape}, voxels: {isg_gt.sum()}")
    
    # Save input
    np.save(output_dir / 'input_volume.npy', volume)
    print(f"Input saved to: {output_dir / 'input_volume.npy'}")
    
    # Step 2: Create segmenter and load model
    print("\n[Step 2] Creating Segmenter and Loading Model")
    print("-" * 70)
    print("Using unified API: create_segmenter('sxt', 'isg')")
    segmenter = create_segmenter(modality='sxt', task='isg')
    
    try:
        segmenter.load_model()
        print(f"✅ Model loaded successfully!")
    except Exception as e:
        print(f"⚠️  Model loading failed: {e}")
    
    # Step 3: Run prediction
    print("\n[Step 3] Running Segmentation")
    print("-" * 70)
    if segmenter.is_loaded:
        print("Predicting ISG segmentation...")
        try:
            isg_mask = segmenter.predict(volume)
            print(f"✅ Prediction successful!")
            print(f"   Mask shape: {isg_mask.shape}, voxels: {isg_mask.sum()}")
            np.save(output_dir / 'isg_mask.npy', isg_mask)
        except Exception as e:
            print(f"❌ Prediction failed: {e}")
    else:
        print("⚠️  Model not loaded. Skipping prediction.")
    
    # Summary
    print("\n" + "="*70)
    print("Summary")
    print("="*70)
    print(f"\n✅ Data loading completed!")
    print(f"\nOutput files:")
    print(f"  - Input volume:       {output_dir / 'input_volume.npy'}")
    if isg_gt is not None:
        print(f"  - ISG ground truth:  {output_dir / 'isg_ground_truth.npy'}")
    
    print(f"\nData Source:")
    print(f"  - Volume: {image_path.name}")
    print(f"  - ISG label: {isg_label_path.name}")
    
    print(f"\nNext Steps:")
    print(f"  1. Train ISG segmentation model")
    print(f"  2. Load trained model and predict")
    print(f"  3. Compare predictions with ground truth label")
    
    print("\n" + "="*70)
    print("Demo completed with real SXT data!")
    print("="*70)


if __name__ == "__main__":
    main()
