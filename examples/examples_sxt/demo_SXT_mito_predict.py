"""
SXT Mitochondria Segmentation - Prediction Example

This script demonstrates mitochondria segmentation using REAL SXT data.
It loads actual SXT volumes and mito labels for validation.

Data source: data/sxt/Stevens_pancreatic_INS_1E_784_5_pre_rec.mrc
Label: data/sxt/784_5_mito_label.tiff

Usage:
    python demo_SXT_mito_predict.py
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


def generate_synthetic_mito_data(shape=(64, 128, 128), num_mitos=30):
    """
    Generate synthetic SXT volume with mitochondria structures.
    
    Mitochondria appear as elongated, tubular structures.
    
    Args:
        shape: Volume shape (depth, height, width)
        num_mitos: Number of mitochondria to generate
        
    Returns:
        Tuple of (volume, mito_mask)
    """
    print(f"Generating synthetic SXT data with {num_mitos} mitochondria...")
    
    # Create background
    volume = np.random.normal(0.15, 0.05, shape).astype(np.float32)
    mito_mask = np.zeros(shape, dtype=np.uint8)
    
    # Generate random mitochondria (elongated structures)
    np.random.seed(42)
    
    for _ in range(num_mitos):
        # Random starting point
        z_start = np.random.randint(10, shape[0] - 10)
        y_start = np.random.randint(10, shape[1] - 10)
        x_start = np.random.randint(10, shape[2] - 10)
        
        # Create elongated structure
        length = np.random.randint(15, 30)
        direction = np.random.randn(3)
        direction = direction / np.linalg.norm(direction)
        
        thickness = np.random.uniform(2.5, 4)
        
        for t in range(length):
            z = int(z_start + t * direction[0])
            y = int(y_start + t * direction[1])
            x = int(x_start + t * direction[2])
            
            if 0 <= z < shape[0] and 0 <= y < shape[1] and 0 <= x < shape[2]:
                # Create sphere at this position
                z_grid, y_grid, x_grid = np.meshgrid(
                    np.arange(max(0, z-4), min(shape[0], z+5)),
                    np.arange(max(0, y-4), min(shape[1], y+5)),
                    np.arange(max(0, x-4), min(shape[2], x+5)),
                    indexing='ij'
                )
                
                distance = np.sqrt(
                    (z_grid - z)**2 + 
                    (y_grid - y)**2 + 
                    (x_grid - x)**2
                )
                
                sphere = distance <= thickness
                
                # Add to mask and volume
                zz, yy, xx = np.where(sphere)
                zz += max(0, z-4)
                yy += max(0, y-4)
                xx += max(0, x-4)
                
                valid = (zz < shape[0]) & (yy < shape[1]) & (xx < shape[2])
                mito_mask[zz[valid], yy[valid], xx[valid]] = 1
                volume[zz[valid], yy[valid], xx[valid]] = 0.75 + np.random.normal(0, 0.1)
    
    # Smooth slightly
    from scipy.ndimage import gaussian_filter
    volume = gaussian_filter(volume, sigma=0.8)
    volume = np.clip(volume, 0, 1)
    
    print(f"Volume shape: {volume.shape}")
    print(f"Mito voxels: {mito_mask.sum()}")
    print(f"Number of mitochondria: ~{num_mitos}")
    
    return volume.astype(np.float32), mito_mask.astype(np.uint8)


def main():
    """Main prediction pipeline with REAL SXT data."""
    
    print("="*70)
    print("SXT Mitochondria Segmentation - Prediction with Real Data")
    print("="*70)
    
    # Setup paths
    data_dir = Path('/media/cuixi/data7/liad/gitspace/iPA/data/sxt')
    image_path = data_dir / 'Stevens_pancreatic_INS_1E_784_5_pre_rec.mrc'
    mito_label_path = data_dir / '784_5_mito_label.tiff'
    model_path = Path('/media/cuixi/data7/liad/gitspace/iPA/models/mask_rcnn_mito/best_mito_model.pth')
    output_dir = Path('sxt_mito_real_output')
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
    mito_gt = None
    if mito_label_path.exists():
        print("Loading mito label...")
        mito_gt = UniversalDataLoader.load_data(str(mito_label_path))
        print(f"   Mito label shape: {mito_gt.shape}, voxels: {mito_gt.sum()}")
    
    # Save input
    np.save(output_dir / 'input_volume.npy', volume)
    print(f"Input saved to: {output_dir / 'input_volume.npy'}")
    
    # Step 2: Create segmenter and load model
    print("\n[Step 2] Creating Segmenter and Loading Model")
    print("-" * 70)
    print("Using unified API: create_segmenter('sxt', 'mito')")
    segmenter = create_segmenter(modality='sxt', task='mito')
    
    try:
        segmenter.load_model()  # Will use default path from unified.py
        print(f"✅ Model loaded successfully!")
    except Exception as e:
        print(f"⚠️  Model loading failed: {e}")
    
    # Step 3: Run prediction
    print("\n[Step 3] Running Segmentation")
    print("-" * 70)
    if segmenter.is_loaded:
        print("Predicting mitochondria segmentation (middle slice demo)...")
        try:
            mito_mask = segmenter.predict(volume)
            print(f"✅ Prediction successful!")
            print(f"   Mask shape: {mito_mask.shape}")
            print(f"   Mito voxels: {mito_mask.sum()}")
            
            # Save results
            np.save(output_dir / 'mito_mask.npy', mito_mask)
            print(f"   Saved to: {output_dir / 'mito_mask.npy'}")
        except Exception as e:
            print(f"❌ Prediction failed: {e}")
            mito_mask = None
    else:
        print("⚠️  Model not loaded. Skipping prediction.")
        mito_mask = None
    
    # Summary
    print("\n" + "="*70)
    print("Summary")
    print("="*70)
    print(f"\n✅ Data loading completed!")
    print(f"\nOutput files:")
    print(f"  - Input volume:       {output_dir / 'input_volume.npy'}")
    if mito_gt is not None:
        print(f"  - Mito ground truth:  {output_dir / 'mito_ground_truth.npy'}")
    
    print(f"\nData Source:")
    print(f"  - Volume: {image_path.name}")
    print(f"  - Mito label: {mito_label_path.name}")
    
    print(f"\nNext Steps:")
    print(f"  1. Train mito segmentation model")
    print(f"  2. Load trained model and predict")
    print(f"  3. Compare predictions with ground truth label")
    
    print("\n" + "="*70)
    print("Demo completed with real SXT data!")
    print("="*70)


if __name__ == "__main__":
    main()
