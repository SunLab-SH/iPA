import sys
import os
# Add project root to path
sys.path.append(r'd:\Gitspace\ipa_full\iPA')

# Create necessary __init__.py files if missing (I already did this, but good to be safe in thought)
# Now import
try:
    from ipa.processing.denoising.n2n_package.src.datasets import NoisyDataset
    from argparse import Namespace
    from torch.utils.data import DataLoader
    import torch

    params = Namespace(
        noise_type='gaussian',
        noise_param=25,
        crop_size=64,
        clean_targets=False,
        seed=None,
        batch_size=4
    )

    root_dir = r'd:\Gitspace\ipa_full\data2\data\sxt_images'
    print(f"Testing dataset with root: {root_dir}")
    dataset = NoisyDataset(root_dir, redux=0, crop_size=64, clean_targets=False, noise_dist=('gaussian', 25.0))
    print(f"Dataset length: {len(dataset)}")
    
    if len(dataset) == 0:
        print("Dataset is empty!")
        sys.exit(1)

    loader = DataLoader(dataset, batch_size=4)
    print("Attempting to load one batch...")
    for source, target in loader:
        print(f"Source shape: {source.shape}, Target shape: {target.shape}")
        print(f"Source dtype: {source.dtype}, Target dtype: {target.dtype}")
        break
    print("Success!")

except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
