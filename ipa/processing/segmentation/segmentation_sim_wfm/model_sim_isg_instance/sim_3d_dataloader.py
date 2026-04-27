#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Custom 3D SIM Data Loader for PM+NE Segmentation
Dynamically extracts 2D slices from 3D volumes during training
"""

import os
import glob
import numpy as np
import torch
from torch.utils.data import Dataset
from PIL import Image
import random
from ipa.data_loader import UniversalDataLoader


class SIM3DDataset(Dataset):
    """
    Dataset for 3D SIM data that dynamically extracts 2D slices during training
    """
    
    def __init__(self, img_dir, mask_dir, ids_list=None, img_scale=1.0, mask_suffix='_label'):
        """
        Initialize SIM 3D Dataset
        
        Args:
            img_dir: Directory containing 3D image files
            mask_dir: Directory containing 3D mask files  
            ids_list: List of IDs to use (if None, use all available)
            img_scale: Scale factor for images
            mask_suffix: Suffix for mask files
        """
        self.img_dir = img_dir
        self.mask_dir = mask_dir
        self.img_scale = img_scale
        self.mask_suffix = mask_suffix
        
        # Find all available 3D files
        img_files = glob.glob(os.path.join(img_dir, '*.tiff'))
        
        self.data_pairs = []
        
        for img_file in img_files:
            filename = os.path.basename(img_file)
            if filename.endswith('.tiff'):
                # Extract ID from filename
                file_id = filename.replace('.tiff', '')
                
                # Check if this ID should be included
                if ids_list is not None:
                    id_without_date = file_id.replace('20220909_', '')
                    if f'20220909_{id_without_date}' not in ids_list:
                        continue
                
                # Find corresponding mask file
                mask_file = os.path.join(mask_dir, f"{file_id}{mask_suffix}.tiff")
                
                if os.path.exists(mask_file):
                    # Load 3D data to get slice information
                    try:
                        img_3d = UniversalDataLoader.load_data(img_file)
                        mask_3d = UniversalDataLoader.load_data(mask_file)
                        
                        # Find slices with mask data
                        z_with_data = np.any(mask_3d > 0, axis=(1, 2))
                        valid_z_indices = np.where(z_with_data)[0]
                        
                        # Add each valid slice as a training sample
                        for z_idx in valid_z_indices:
                            self.data_pairs.append({
                                'img_file': img_file,
                                'mask_file': mask_file,
                                'z_index': z_idx,
                                'id': f"{file_id}_z{z_idx:03d}"
                            })
                            
                    except Exception as e:
                        print(f"Error loading {file_id}: {e}")
                        continue
        
        print(f"SIM3DDataset initialized with {len(self.data_pairs)} 2D slices from 3D volumes")
    
    def __len__(self):
        return len(self.data_pairs)
    
    def __getitem__(self, idx):
        """
        Get a 2D slice from 3D volume
        """
        data_info = self.data_pairs[idx]
        
        # Load 3D volumes
        img_3d = UniversalDataLoader.load_data(data_info['img_file'])
        mask_3d = UniversalDataLoader.load_data(data_info['mask_file'])
        
        # Extract 2D slice
        z_idx = data_info['z_index']
        img_2d = img_3d[z_idx]  # Shape: (H, W, C)
        mask_2d = mask_3d[z_idx]  # Shape: (H, W)
        
        # Convert to PIL Images for transforms
        img_pil = Image.fromarray(img_2d.astype(np.uint8))
        mask_pil = Image.fromarray(mask_2d.astype(np.uint8))
        
        # Apply scaling if needed
        if self.img_scale != 1.0:
            new_width = int(img_2d.shape[1] * self.img_scale)
            new_height = int(img_2d.shape[0] * self.img_scale)
            img_pil = img_pil.resize((new_width, new_height), Image.BILINEAR)
            mask_pil = mask_pil.resize((new_width, new_height), Image.NEAREST)
        
        # Convert back to numpy
        img_array = np.array(img_pil)
        mask_array = np.array(mask_pil)
        
        # Convert to torch tensors
        # Image: (H, W, C) -> (C, H, W) and normalize to [0, 1]
        img_tensor = torch.from_numpy(img_array).permute(2, 0, 1).float() / 255.0
        
        # Mask: (H, W) -> (H, W) as long tensor
        mask_tensor = torch.from_numpy(mask_array).long()
        
        return {
            'image': img_tensor,
            'mask': mask_tensor,
            'id': data_info['id']
        }


def create_sim_data_loaders(img_dir, mask_dir, allowed_ids=None, batch_size=4, 
                           val_percent=0.2, img_scale=1.0, num_workers=4, data_usage_ratio=1.0):
    """
    Create train and validation data loaders for SIM 3D data
    
    Args:
        img_dir: Directory containing 3D image files
        mask_dir: Directory containing 3D mask files
        allowed_ids: List of allowed IDs (if None, use all)
        batch_size: Batch size for training
        val_percent: Percentage of data to use for validation
        img_scale: Scale factor for images
        num_workers: Number of worker processes
        data_usage_ratio: Fraction of data to use (0.1 = 10%, 1.0 = 100%)
    
    Returns:
        train_loader, val_loader
    """
    # Create full dataset
    full_dataset = SIM3DDataset(
        img_dir=img_dir,
        mask_dir=mask_dir,
        ids_list=allowed_ids,
        img_scale=img_scale
    )
    
    # Apply data usage ratio (randomly sample subset)
    if data_usage_ratio < 1.0:
        import random
        total_samples = len(full_dataset)
        n_samples_to_use = int(total_samples * data_usage_ratio)
        
        # Randomly sample indices
        random.seed(42)  # For reproducibility
        indices_to_use = random.sample(range(total_samples), n_samples_to_use)
        
        # Create subset dataset
        from torch.utils.data import Subset
        full_dataset = Subset(full_dataset, indices_to_use)
        
        print(f"Using {data_usage_ratio*100:.1f}% of data: {n_samples_to_use}/{total_samples} samples")
    
    # Split into train and validation
    n_val = int(len(full_dataset) * val_percent)
    n_train = len(full_dataset) - n_val
    
    train_dataset, val_dataset = torch.utils.data.random_split(
        full_dataset, [n_train, n_val],
        generator=torch.Generator().manual_seed(42)
    )
    
    # Create data loaders
    train_loader = torch.utils.data.DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=num_workers,
        pin_memory=True,
        drop_last=True
    )
    
    val_loader = torch.utils.data.DataLoader(
        val_dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=num_workers,
        pin_memory=True,
        drop_last=False
    )
    
    return train_loader, val_loader


def create_isg_data_loaders(img_dir, mask_dir, allowed_ids=None, batch_size=4, 
                           val_percent=0.2, img_scale=1.0, num_workers=4, data_usage_ratio=1.0):
    """
    Create train and validation data loaders for ISG 3D data
    This is an alias for create_sim_data_loaders for ISG training compatibility
    
    Args:
        img_dir: Directory containing 3D image files
        mask_dir: Directory containing 3D mask files
        allowed_ids: List of allowed IDs (if None, use all)
        batch_size: Batch size for training
        val_percent: Percentage of data to use for validation
        img_scale: Scale factor for images
        num_workers: Number of worker processes
        data_usage_ratio: Fraction of data to use (0.1 = 10%, 1.0 = 100%)
    
    Returns:
        train_loader, val_loader
    """
    # Just call the existing function - same logic works for ISG data
    return create_sim_data_loaders(
        img_dir=img_dir,
        mask_dir=mask_dir,
        allowed_ids=allowed_ids,
        batch_size=batch_size,
        val_percent=val_percent,
        img_scale=img_scale,
        num_workers=num_workers,
        data_usage_ratio=data_usage_ratio
    )