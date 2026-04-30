"""
SXT Composite Mask Generation - Combining multiple organelle masks into one.

This script demonstrates how to generate a composite mask that includes:
- Background (0)
- Cytoplasm/Cell Membrane (1)
- Nucleus (2)
- Mitochondria (3)
- Insulin Secretory Granules (ISGs) (4)
"""

import numpy as np
from pathlib import Path
from ipa.processing.segmentation import create_segmenter
from ipa.data_loader import UniversalDataLoader

def main():
    print("="*70)
    print("SXT Composite Mask Generation (Full Organelle Set)")
    print("="*70)
    
    # Setup paths
    data_dir = Path('/media/cuixi/data7/liad/gitspace/iPA/data/sxt')
    image_path = data_dir / 'Stevens_pancreatic_INS_1E_784_5_pre_rec.mrc'
    output_dir = Path('sxt_composite_output')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Load Data
    print("\n[Step 1] Loading Real SXT Data...")
    volume = UniversalDataLoader.load_data(str(image_path))
    print(f"✅ Loaded! Shape: {volume.shape}")
    
    # Step 2: Create Segmenters
    print("\n[Step 2] Initializing Segmenters...")
    cell_seg = create_segmenter('sxt', 'cell')
    mito_seg = create_segmenter('sxt', 'mito')
    isg_seg = create_segmenter('sxt', 'isg')
    
    # Load models
    try:
        cell_seg.load_model()
        mito_seg.load_model()
        isg_seg.load_model()
        print("✅ All models loaded successfully!")
    except Exception as e:
        print(f"⚠️ Model loading error: {e}")
        return

    # Step 3: Run Segmentation / Load Labels
    print("\n[Step 3] Processing Center Slice...")
    mid_idx = volume.shape[0] // 2
    slice_vol = volume[mid_idx:mid_idx+1] 
    
    print(f"   Processing slice {mid_idx}...")
    
    # 1. Cell & Nucleus
    cell_result = cell_seg.predict(slice_vol)
    cell_mask = cell_result['cell_mask'][0]
    nucleus_mask = cell_result['nucleus_mask'][0]
    
    # 2. Mito (Using real label for demo)
    mito_label_path = data_dir / '784_5_mito_label.tiff'
    if mito_label_path.exists():
        full_mito_vol = UniversalDataLoader.load_data(str(mito_label_path))
        mito_mask = (full_mito_vol[mid_idx] > 0).astype(np.uint8)
    else:
        mito_mask = mito_seg.predict(slice_vol)[0]
    
    # 3. ISG (Using real label for demo)
    isg_label_path = data_dir / '784_5_isg_label.tiff'
    if isg_label_path.exists():
        full_isg_vol = UniversalDataLoader.load_data(str(isg_label_path))
        isg_mask = (full_isg_vol[mid_idx] > 0).astype(np.uint8)
    else:
        isg_mask = isg_seg.predict(slice_vol)
    
    print(f"   Cell voxels: {cell_mask.sum()}")
    print(f"   Nucleus voxels: {nucleus_mask.sum()}")
    print(f"   Mito voxels: {mito_mask.sum()}")
    print(f"   ISG voxels: {isg_mask.sum()}")
    
    # Step 4: Composite the masks
    print("\n[Step 4] Generating Composite Mask...")
    # Priority: ISG (4) > Mito (3) > Nucleus (2) > Cytoplasm (1) > Background (0)
    
    composite = np.zeros_like(cell_mask, dtype=np.uint8)
    
    # 1. Fill Cell area (Cytoplasm + others)
    composite[cell_mask > 0] = 1
    
    # 2. Fill Nucleus
    composite[nucleus_mask > 0] = 2
    
    # 3. Fill Mito (in cytoplasm)
    composite[(mito_mask > 0) & (nucleus_mask == 0)] = 3
    
    # 4. Fill ISG (in cytoplasm, highest priority among organelles)
    composite[(isg_mask > 0) & (nucleus_mask == 0)] = 4
    
    print(f"   Composite mask stats:")
    print(f"   - Background (0):     {(composite == 0).sum()}")
    print(f"   - Cytoplasm (1):      {(composite == 1).sum()}")
    print(f"   - Nucleus (2):        {(composite == 2).sum()}")
    print(f"   - Mitochondria (3):   {(composite == 3).sum()}")
    print(f"   - ISGs (4):           {(composite == 4).sum()}")
    
    # Save results
    np.save(output_dir / f'slice_{mid_idx}_composite.npy', composite)
    print(f"\n✅ Composite mask saved to: {output_dir / f'slice_{mid_idx}_composite.npy'}")
    
    print("\n" + "="*70)
    print("Demo completed! You can now use this composite mask for spatial analysis.")
    print("="*70)

if __name__ == "__main__":
    main()
