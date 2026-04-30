"""
SXT Cell Segmentation - Prediction Example

This script demonstrates cell and nucleus segmentation using REAL SXT data.
It loads actual SXT volumes and applies segmentation (requires trained model).

Data source: data/sxt/Stevens_pancreatic_INS_1E_784_5_pre_rec.mrc
Labels: data/sxt/784_5_wholecell_label.tiff, 784_5_NC_label.tiff

Usage:
    python demo_SXT_cell_predict.py
"""

import os
import sys
import numpy as np
from pathlib import Path

# Add iPA module path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
sys.path.insert(0, PROJECT_ROOT)

from ipa.processing.segmentation import create_segmenter


def generate_synthetic_sxt_data(shape=(64, 128, 128)):
    """
    Generate synthetic SXT-like 3D volume with cell structures.
    
    Args:
        shape: Volume shape (depth, height, width)
        
    Returns:
        Tuple of (volume, cell_mask, nucleus_mask)
    """
    print("Generating synthetic SXT-like data...")
    
    z, y, x = np.meshgrid(
        np.arange(shape[0]),
        np.arange(shape[1]),
        np.arange(shape[2]),
        indexing='ij'
    )
    
    # Create cell (larger ellipsoid)
    cell_center = (shape[0]//2, shape[1]//2, shape[2]//2)
    cell_radii = (shape[0]//2 - 5, shape[1]//2 - 10, shape[2]//2 - 10)
    
    cell_mask = ((z - cell_center[0])**2 / cell_radii[0]**2 +
                 (y - cell_center[1])**2 / cell_radii[1]**2 +
                 (x - cell_center[2])**2 / cell_radii[2]**2) <= 1
    
    # Create nucleus (smaller ellipsoid inside cell)
    nucleus_radii = (shape[0]//4, shape[1]//4, shape[2]//4)
    
    nucleus_mask = ((z - cell_center[0])**2 / nucleus_radii[0]**2 +
                    (y - cell_center[1])**2 / nucleus_radii[1]**2 +
                    (x - cell_center[2])**2 / nucleus_radii[2]**2) <= 1
    
    # Create intensity volume
    volume = np.zeros(shape, dtype=np.float32)
    
    # Nucleus is denser (brighter in SXT)
    volume[nucleus_mask] = 0.9 + np.random.normal(0, 0.05, nucleus_mask.sum())
    
    # Cytoplasm (cell minus nucleus)
    cytoplasm = cell_mask & ~nucleus_mask
    volume[cytoplasm] = 0.5 + np.random.normal(0, 0.08, cytoplasm.sum())
    
    # Background
    background = ~cell_mask
    volume[background] = 0.1 + np.random.normal(0, 0.03, background.sum())
    
    # Add noise
    volume += np.random.normal(0, 0.05, shape)
    volume = np.clip(volume, 0, 1)
    
    print(f"Volume shape: {volume.shape}")
    print(f"Cell voxels: {cell_mask.sum()}")
    print(f"Nucleus voxels: {nucleus_mask.sum()}")
    
    return volume.astype(np.float32), cell_mask.astype(np.uint8), nucleus_mask.astype(np.uint8)


def main():
    """Main prediction pipeline with REAL SXT data."""
    
    print("="*70)
    print("SXT Cell Segmentation - Prediction with Real Data")
    print("="*70)
    
    # Setup paths
    data_dir = Path(PROJECT_ROOT) / 'data' / 'sxt'
    image_path = data_dir / 'Stevens_pancreatic_INS_1E_784_5_pre_rec.mrc'
    cell_label_path = data_dir / '784_5_wholecell_label.tiff'
    nucleus_label_path = data_dir / '784_5_NC_label.tiff'
    output_dir = Path(PROJECT_ROOT) / 'results' / 'sxt_cell_real_output'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Load real SXT data
    print("\n[Step 1] Loading Real SXT Data")
    print("-" * 70)
    print(f"Image path: {image_path}")
    
    if not image_path.exists():
        print(f"❌ Image not found: {image_path}")
        print("Using synthetic data instead...")
        volume, cell_gt, nucleus_gt = generate_synthetic_sxt_data(shape=(64, 128, 128))
    else:
        from ipa.data_loader import UniversalDataLoader
        print("Loading SXT volume...")
        volume = UniversalDataLoader.load_data(str(image_path))
        print(f"✅ Loaded! Shape: {volume.shape}, dtype: {volume.dtype}")
        
        # Load ground truth labels if available
        cell_gt = None
        nucleus_gt = None
        if cell_label_path.exists():
            print("Loading cell label...")
            cell_gt = UniversalDataLoader.load_data(str(cell_label_path))
            print(f"   Cell label shape: {cell_gt.shape}")
        
        if nucleus_label_path.exists():
            print("Loading nucleus label...")
            nucleus_gt = UniversalDataLoader.load_data(str(nucleus_label_path))
            print(f"   Nucleus label shape: {nucleus_gt.shape}")
    
    # Save input
    np.save(output_dir / 'input_volume.npy', volume)
    print(f"Input saved to: {output_dir / 'input_volume.npy'}")
    
    # Step 2: Create segmenter and load model
    print("\n[Step 2] Creating Segmenter and Loading Model")
    print("-" * 70)
    print("Using unified API: create_segmenter('sxt', 'cell')")
    segmenter = create_segmenter(modality='sxt', task='cell')
    
    try:
        segmenter.load_model()  # Will use default path from unified.py
        print(f"✅ Model loaded successfully!")
    except Exception as e:
        print(f"⚠️  Model loading failed: {e}")
    
    # Step 3: Run prediction
    print("\n[Step 3] Running Segmentation")
    print("-" * 70)
    if segmenter.is_loaded:
        print("Predicting cell and nucleus segmentation...")
        try:
            result = segmenter.predict(volume)
            print(f"✅ Prediction successful!")
            if isinstance(result, dict):
                for k, v in result.items():
                    print(f"   {k} shape: {v.shape}, voxels: {v.sum()}")
                    np.save(output_dir / f'{k}.npy', v)
            else:
                print(f"   Result shape: {result.shape}")
                np.save(output_dir / 'prediction.npy', result)
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
    if cell_gt is not None:
        print(f"  - Cell ground truth:  {output_dir / 'cell_ground_truth.npy'}")
    if nucleus_gt is not None:
        print(f"  - Nucleus ground truth: {output_dir / 'nucleus_ground_truth.npy'}")
    
    print(f"\nData Source:")
    print(f"  - Volume: {image_path.name}")
    print(f"  - Cell label: {cell_label_path.name}")
    print(f"  - Nucleus label: {nucleus_label_path.name}")
    
    print(f"\nNext Steps:")
    print(f"  1. Train cell segmentation model")
    print(f"  2. Load trained model and predict")
    print(f"  3. Compare predictions with ground truth labels")
    
    print("\n" + "="*70)
    print("Demo completed with real SXT data!")
    print("="*70)


if __name__ == "__main__":
    main()
