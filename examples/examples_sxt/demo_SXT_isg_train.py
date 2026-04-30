"""
SXT ISG Segmentation - Training Example (Demo)

This script demonstrates the training workflow for SXT ISG (Insulin Granule) segmentation.
ISGs are small spherical organelles that require specialized detection methods.

Note: This is a simplified demo with synthetic data to show the API.
For real training, you need actual SXT datasets with labeled ISG masks.

Usage:
    python demo_SXT_isg_train.py
"""

import os
import sys
import numpy as np
from pathlib import Path

# Add iPA module path
current_dir = os.path.dirname(os.path.abspath(__file__))
ipa_root = os.path.dirname(os.path.dirname(current_dir))
sys.path.insert(0, ipa_root)


def generate_sxt_isg_dataset(num_samples=10, shape=(32, 64, 64), num_isgs=40):
    """
    Generate synthetic SXT training dataset with ISG structures.
    
    ISGs appear as small spherical objects in SXT images.
    
    Args:
        num_samples: Number of training samples
        shape: Volume shape (depth, height, width)
        num_isgs: Number of ISGs per sample
        
    Returns:
        Tuple of (volumes, isg_masks) lists
    """
    print(f"Generating {num_samples} synthetic SXT ISG samples...")
    
    volumes = []
    isg_masks = []
    
    for i in range(num_samples):
        # Create background
        volume = np.random.normal(0.1, 0.03, shape).astype(np.float32)
        isg_mask = np.zeros(shape, dtype=np.uint8)
        
        # Generate random ISGs (small spheres)
        np.random.seed(i + 42)
        
        for _ in range(num_isgs):
            # Random center
            z = np.random.randint(5, shape[0] - 5)
            y = np.random.randint(5, shape[1] - 5)
            x = np.random.randint(5, shape[2] - 5)
            
            # Random radius (2-5 voxels)
            radius = np.random.uniform(2, 5)
            
            # Create sphere
            z_grid, y_grid, x_grid = np.meshgrid(
                np.arange(max(0, int(z-radius)-1), min(shape[0], int(z+radius)+2)),
                np.arange(max(0, int(y-radius)-1), min(shape[1], int(y+radius)+2)),
                np.arange(max(0, int(x-radius)-1), min(shape[2], int(x+radius)+2)),
                indexing='ij'
            )
            
            distance = np.sqrt(
                (z_grid - z)**2 + 
                (y_grid - y)**2 + 
                (x_grid - x)**2
            )
            
            sphere = distance <= radius
            
            # Add to mask and volume
            zz, yy, xx = np.where(sphere)
            zz += max(0, int(z-radius)-1)
            yy += max(0, int(y-radius)-1)
            xx += max(0, int(x-radius)-1)
            
            valid = (zz < shape[0]) & (yy < shape[1]) & (xx < shape[2])
            isg_mask[zz[valid], yy[valid], xx[valid]] = 1
            volume[zz[valid], yy[valid], xx[valid]] = 0.7 + np.random.normal(0, 0.1)
        
        # Smooth slightly
        from scipy.ndimage import gaussian_filter
        volume = gaussian_filter(volume, sigma=0.5)
        volume = np.clip(volume, 0, 1)
        
        volumes.append(volume)
        isg_masks.append(isg_mask)
        
        if (i + 1) % 3 == 0:
            print(f"  Generated {i + 1}/{num_samples} samples")
    
    print(f"✅ Dataset generated!")
    print(f"   Volumes: {len(volumes)} samples, shape {volumes[0].shape}")
    print(f"   ISG masks: {len(isg_masks)} samples")
    print(f"   Avg ISG voxels per sample: {np.mean([m.sum() for m in isg_masks]):.0f}")
    
    return volumes, isg_masks


def demonstrate_training_workflow():
    """Demonstrate the training workflow (pseudo-code for real implementation)."""
    
    print("\n" + "="*70)
    print("SXT ISG Segmentation - Training Workflow Demo")
    print("="*70)
    
    # Setup output directory
    output_dir = Path('sxt_isg_training_demo')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Prepare dataset
    print("\n[Step 1] Preparing Training Dataset")
    print("-" * 70)
    train_volumes, train_isg_masks = generate_sxt_isg_dataset(
        num_samples=10, shape=(32, 64, 64), num_isgs=40
    )
    
    val_volumes, val_isg_masks = generate_sxt_isg_dataset(
        num_samples=3, shape=(32, 64, 64), num_isgs=30
    )
    
    # Save dataset
    np.savez(output_dir / 'train_data.npz',
             volumes=np.array(train_volumes),
             isg_masks=np.array(train_isg_masks))
    
    np.savez(output_dir / 'val_data.npz',
             volumes=np.array(val_volumes),
             isg_masks=np.array(val_isg_masks))
    
    print(f"Dataset saved to: {output_dir}")
    
    # Step 2: Show training code structure
    print("\n[Step 2] Training Code Structure")
    print("-" * 70)
    print("""
# Real training would use this pattern:

import torch
from torch.utils.data import DataLoader, Dataset
from ipa.processing.segmentation import create_segmenter

class SXTISGDataset(Dataset):
    def __init__(self, volumes, isg_masks):
        self.volumes = volumes
        self.isg_masks = isg_masks
    
    def __len__(self):
        return len(self.volumes)
    
    def __getitem__(self, idx):
        volume = torch.from_numpy(self.volumes[idx]).unsqueeze(0).float()
        isg_mask = torch.from_numpy(self.isg_masks[idx]).unsqueeze(0).float()
        return volume, isg_mask

# Create segmenter
segmenter = create_segmenter('sxt', 'isg')

# Initialize model (sphere-like organelle detection network)
model = segmenter.model

# Setup training
train_dataset = SXTISGDataset(train_volumes, train_isg_masks)
train_loader = DataLoader(train_dataset, batch_size=2, shuffle=True)

optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
criterion = torch.nn.BCEWithLogitsLoss()

# Training loop
for epoch in range(50):
    model.train()
    total_loss = 0
    
    for volumes, isg_masks in train_loader:
        optimizer.zero_grad()
        
        # Forward pass
        isg_pred = model(volumes)
        
        # Calculate loss
        loss = criterion(isg_pred, isg_masks)
        
        # Backward pass
        loss.backward()
        optimizer.step()
        
        total_loss += loss.item()
    
    avg_loss = total_loss / len(train_loader)
    print(f"Epoch {epoch+1}/50, Loss: {avg_loss:.4f}")

# Save model
segmenter.save_model('sxt_isg_trained.pth')
    """)
    
    # Step 3: Explain key points
    print("\n[Step 3] Key Training Points")
    print("-" * 70)
    print("""
Important considerations for SXT ISG segmentation:

1. Data Characteristics:
   - ISGs are small spherical organelles (insulin granules)
   - Size: typically 200-400 nm diameter (2-5 voxels in SXT)
   - Appear as bright spots in SXT images
   - Can be numerous (tens to hundreds per cell)

2. Model Architecture:
   - 3D U-Net or specialized sphere detection network
   - Input: 1 channel (grayscale intensity)
   - Output: 1 channel (binary ISG mask)
   - Consider using multi-scale features for different ISG sizes

3. Training Parameters:
   - Learning rate: 1e-3 (Adam optimizer)
   - Batch size: 2-4 (small objects, can use larger batches)
   - Epochs: 50-100
   - Loss: BCEWithLogitsLoss + Focal loss (for class imbalance)
   - Weight decay: 1e-5

4. Challenges:
   - Small object detection (easy to miss)
   - Class imbalance (few ISG voxels vs many background voxels)
   - Dense packing (ISGs can be close together)
   - Need good resolution to detect small granules

5. Data Augmentation:
   - Random rotation (3D)
   - Random flip
   - Intensity scaling
   - Gaussian noise
   - Random cropping (focus on regions with ISGs)

6. Evaluation Metrics:
   - Dice coefficient
   - Precision and Recall (important for small objects)
   - F1 score
   - Object-level detection rate (count accuracy)

7. Post-processing:
   - Connected component analysis
   - Size filtering (remove artifacts)
   - Watershed separation for touching ISGs

8. Hardware Requirements:
   - GPU with ≥8GB VRAM
   - RAM: 16GB+
   - Storage: 30GB+ for datasets
    """)
    
    # Step 4: Next steps
    print("\n[Step 4] Next Steps for Real Training")
    print("-" * 70)
    print("""
To train on real SXT ISG data:

1. Prepare your dataset:
   - Collect SXT volumes with insulin granules
   - Manually annotate ISG masks (challenging due to small size!)
   - Split into train/val/test sets
   - Verify annotation quality carefully

2. Update the training script:
   - Replace synthetic data with real data loader
   - Use focal loss or weighted BCE for class imbalance
   - Implement multi-scale training
   - Add validation during training

3. Train:
   python demo_SXT_isg_train.py --data_dir /path/to/data --epochs 100

4. Evaluate:
   - Calculate Dice score on test set
   - Count ISGs and compare with ground truth
   - Check for false positives (artifacts)
   - Analyze missed detections (small ISGs)

5. Post-processing:
   - Apply size filters
   - Separate touching ISGs
   - Validate count accuracy

6. Use trained model:
   segmenter = create_segmenter('sxt', 'isg')
   segmenter.load_model('sxt_isg_trained.pth')
   isg_mask = segmenter.predict(new_volume)
   
   # Analyze ISG properties
   from skimage.measure import label, regionprops
   labeled = label(isg_mask)
   props = regionprops(labeled)
   print(f"Number of ISGs: {len(props)}")
    """)
    
    # Summary
    print("\n" + "="*70)
    print("Training Demo Completed")
    print("="*70)
    print(f"\nOutput:")
    print(f"  - Training dataset: {output_dir / 'train_data.npz'}")
    print(f"  - Validation set: {output_dir / 'val_data.npz'}")
    print(f"  - Training template: See code above")
    
    print(f"\nNote:")
    print(f"  This is a demonstration with synthetic data.")
    print(f"  For real training, replace with actual SXT datasets.")
    print(f"  The unified API (create_segmenter) remains the same.")
    
    print("\n" + "="*70)


if __name__ == "__main__":
    demonstrate_training_workflow()
