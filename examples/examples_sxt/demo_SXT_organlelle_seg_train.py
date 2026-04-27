#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
SXT Training Examples for different segmentation models
"""

import os
import sys
import time
import torch
import glob
import warnings
import numpy as np
from pathlib import Path

# Filter MRC file warnings at the very beginning
warnings.filterwarnings("ignore", category=RuntimeWarning, module="mrcfile")
warnings.filterwarnings("ignore", message=".*Unrecognised machine stamp.*")
warnings.filterwarnings("ignore", message=".*Map ID string not found.*")
warnings.filterwarnings("ignore", message=".*data block cannot be read.*")

from ipa.data_loader import QuickLogger
from parsers import get_args






def validate_data_pairs(img_dir, mask_dir, all_ids, model_type, logger):
    """Validate image-mask pairs and return valid training IDs"""
    
    # Log available files and check for ID matches
    img_files = glob.glob(os.path.join(img_dir, '*.mrc'))
    mask_files = glob.glob(os.path.join(mask_dir, '*.tiff'))
    
    logger.step(f"Found {len(img_files)} image files and {len(mask_files)} mask files")
    
    # Check which IDs are actually available and validate file pairs
    available_ids = set()
    valid_pairs = {}
    missing_details = {}
    
    for img_file in img_files:
        filename = os.path.basename(img_file)
        # Extract ID from filename with improved logic
        parts = filename.split('_')
        file_id = None
        
        # Try multiple extraction methods
        for i in range(len(parts) - 1):
            if parts[i].isdigit() and parts[i+1].isdigit():
                file_id = f"{parts[i]}_{parts[i+1]}"
                break
        
        # Additional extraction method for edge cases
        if file_id is None:
            for i in range(len(parts)):
                if parts[i].isdigit():
                    for j in range(i+1, len(parts)):
                        if parts[j].isdigit():
                            file_id = f"{parts[i]}_{parts[j]}"
                            break
                    if file_id:
                        break
        
        if file_id:
            # Check if corresponding mask exists with multiple patterns
            if model_type == 'isg':
                mask_patterns = [
                    f"{file_id}_isg_label.tiff",
                    f"{file_id}_isg_label.tif",
                    f"{file_id}.tiff",
                    f"{file_id}.tif",
                    f"{file_id}_label.tiff",
                    f"{file_id}_label.tif",
                    f"label_{file_id}.tiff",
                    f"label_{file_id}.tif"
                ]
            elif model_type == 'mito':
                mask_patterns = [
                    f"{file_id}_mito_label.tiff",
                    f"{file_id}_mito_label.tif",
                    f"{file_id}_mito.tiff",
                    f"{file_id}_mito.tif",
                    f"{file_id}.tiff",
                    f"{file_id}.tif"
                ]
            elif model_type == 'cell_nucleus':
                mask_patterns = [
                    f"{file_id}_nu_label.tiff",
                    f"{file_id}_nu_label.tif",
                    f"{file_id}_nu.tiff",
                    f"{file_id}_nu.tif",
                    f"{file_id}_cell_nucleus.tiff",
                    f"{file_id}_cell_nucleus.tif",
                    f"{file_id}.tiff",
                    f"{file_id}.tif"
                ]
            elif model_type == 'pm_ne':
                mask_patterns = [
                    f"{file_id}_pm_ne_label.tiff",
                    f"{file_id}_pm_ne_label.tif",
                    f"{file_id}_pm_ne.tiff",
                    f"{file_id}_pm_ne.tif",
                    f"{file_id}.tiff",
                    f"{file_id}.tif"
                ]
            else:
                # Default patterns
                mask_patterns = [
                    f"{file_id}.tiff",
                    f"{file_id}.tif",
                    f"{file_id}_label.tiff",
                    f"{file_id}_label.tif"
                ]
            
            mask_found = False
            attempted_paths = []
            
            for pattern in mask_patterns:
                mask_path = os.path.join(mask_dir, pattern)
                attempted_paths.append(pattern)
                if os.path.exists(mask_path):
                    available_ids.add(file_id)
                    valid_pairs[file_id] = {
                        'img': img_file,
                        'mask': mask_path
                    }
                    mask_found = True
                    logger.step(f"✓ Found pair: {file_id} -> {pattern}")
                    break
            
            if not mask_found:
                # Try partial matching
                for mask_file in mask_files:
                    mask_name = os.path.basename(mask_file)
                    if file_id in mask_name and (model_type in mask_name.lower() or 'label' in mask_name.lower()):
                        available_ids.add(file_id)
                        valid_pairs[file_id] = {
                            'img': img_file,
                            'mask': mask_file
                        }
                        mask_found = True
                        logger.step(f"✓ Found pair (partial match): {file_id} -> {mask_name}")
                        break
            
            if not mask_found:
                missing_details[file_id] = {
                    'img_file': filename,
                    'attempted_patterns': attempted_paths
                }
                logger.step(f"✗ Missing mask for: {file_id}")
        else:
            logger.step(f"✗ Could not extract ID from: {filename}")
    
    # Report detailed missing information
    if missing_details:
        logger.step(f"Detailed missing mask information:")
        for file_id, details in missing_details.items():
            logger.step(f"  ID: {file_id}")
            logger.step(f"    Source image: {details['img_file']}")
            logger.step(f"    Tried patterns: {details['attempted_patterns']}")
            
            # Show actual mask files for comparison
            similar_masks = [f for f in mask_files if file_id.split('_')[0] in os.path.basename(f) or file_id.split('_')[1] in os.path.basename(f)]
            if similar_masks:
                logger.step(f"    Similar mask files found: {[os.path.basename(f) for f in similar_masks[:3]]}")
    
    logger.step(f"Available valid pairs: {len(valid_pairs)}")
    logger.step(f"Available IDs: {sorted(available_ids)}")
    
    # Check for missing IDs
    missing_ids = set(all_ids) - available_ids
    if missing_ids:
        logger.step(f"WARNING: Missing IDs (no valid pairs): {sorted(missing_ids)}")
    
    # Use only available IDs for training
    training_ids = [id for id in all_ids if id in available_ids]
    
    return training_ids, valid_pairs

def train_isg_model():
    """Train ISG segmentation model using 3D data"""
    
    training_start = time.time()
    mainpath = get_args().main_path

    # Define all IDs to use for training (mixed together)
    all_ids = ['784_5', '766_8', '842_17', '783_5', '766_5', '842_12', '766_2', '766_7', 
               '766_10', '766_11', '769_5', '769_7', '783_6', '783_12', '784_4', '784_6', 
               '784_7', '785_7', '822_4', '822_6', '822_7', '842_13', '931_9', '931_14']
    
    # Initialize logger
    log_dir = f'{mainpath}/logs'
    logger = QuickLogger("sxt_isg_training", log_dir=log_dir)
    logger.step("Starting SXT ISG Model Training")

    logger.step(f"Using all {len(all_ids)} IDs for training (mixed together):")
    logger.step(f"IDs: {all_ids}")

    # Data paths
    training_img_dir = r'D:\Gitspace\ipa_full\data\SXT\for_24_datasets'
    training_mask_dir = r'D:\Gitspace\ipa_full\data\SXT\tiff_24_datasets_isg'

    logger.step(f"Image directory: {training_img_dir}")
    logger.step(f"Mask directory: {training_mask_dir}")
    
    # Verify data directories
    if not os.path.exists(training_img_dir) or not os.path.exists(training_mask_dir):
        logger.step("ERROR: Training data directories not found!")
        return

    # Validate data pairs
    training_ids, valid_pairs = validate_data_pairs(
        training_img_dir, training_mask_dir, all_ids, 'isg', logger
    )
    
    logger.step(f"Using {len(training_ids)} available IDs for training: {training_ids}")
    
    if len(training_ids) == 0:
        logger.step("ERROR: No valid training data found!")
        return

    try:
        # Import training modules
        from ipa.processing.segmentation.segmentation_sxt.model_sxt_isg.train import train_net
        from ipa.processing.segmentation.segmentation_sxt.model_sxt_isg.unet import UNet
        
        logger.step("Starting ISG model training...")
        
        # Setup device and model
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        logger.step(f"Using device: {device}")
        
        # Create UNet model for ISG (2 classes: background=0, ISG=1)
        net = UNet(n_channels=3, n_classes=2, bilinear=False)
        net.to(device)
        logger.step("ISG Model created with 2 classes (background, ISG)")
        
        # Training configuration with optimized settings
        training_config = {
            'epochs': 25,
            'batch_size': 16,
            'learning_rate': 2e-4,
            'val_percent': 0.15,
            'img_scale': 0.5,
            'amp': True,
            'custom_img_dir': training_img_dir,
            'custom_mask_dir': training_mask_dir,
            'use_wandb': False,
            'allowed_ids': training_ids,
            'sample_ratio': 0.8,
            'organelle': 'isg'
        }
        
        logger.step(f"Training config: {training_config}")
        logger.step("Starting ISG training...")
        
        # Start training with original train_net function
        train_net(
            net=net,
            device=device,
            external_logger=logger,  # Pass logger to train_net
            **training_config
        )
        
        training_time = time.time() - training_start
        logger.step(f"ISG training completed successfully in {training_time:.2f}s")
        
    except Exception as e:
        training_time = time.time() - training_start
        logger.step(f"Training failed after {training_time:.2f}s: {str(e)}")
        import traceback
        logger.step(f"Error details: {traceback.format_exc()}")

def train_mito_model():
    """Train mitochondria segmentation model using 3D data"""
    
    training_start = time.time()
    mainpath = get_args().main_path

    # Define all IDs to use for training (mixed together)
    all_ids = ['784_5', '766_8', '842_17', '783_5', '766_5', '842_12', '766_2', '766_7', 
               '766_10', '766_11', '769_5', '769_7', '783_6', '783_12', '784_4', '784_6', 
               '784_7', '785_7', '822_4', '822_6', '822_7', '842_13', '931_9', '931_14']
    
    # Initialize logger
    log_dir = f'{mainpath}/logs'
    logger = QuickLogger("sxt_mito_training", log_dir=log_dir)
    logger.step("Starting SXT Mitochondria Model Training")

    logger.step(f"Using all {len(all_ids)} IDs for training (mixed together):")
    logger.step(f"IDs: {all_ids}")

    # Data paths
    training_img_dir = r'D:\Gitspace\ipa_full\data\SXT\for_24_datasets'
    training_mask_dir = r'D:\Gitspace\ipa_full\data\SXT\tiff_24_datasets_mito'

    logger.step(f"Image directory: {training_img_dir}")
    logger.step(f"Mask directory: {training_mask_dir}")
    
    # Verify data directories
    if not os.path.exists(training_img_dir) or not os.path.exists(training_mask_dir):
        logger.step("ERROR: Training data directories not found!")
        return

    # Validate data pairs
    training_ids, valid_pairs = validate_data_pairs(
        training_img_dir, training_mask_dir, all_ids, 'mito', logger
    )
    
    logger.step(f"Using {len(training_ids)} available IDs for training: {training_ids}")
    
    if len(training_ids) == 0:
        logger.step("ERROR: No valid training data found!")
        return

    try:
        # Import training modules
        from ipa.processing.segmentation.segmentation_sxt.model_sxt_isg.train import train_net
        from ipa.processing.segmentation.segmentation_sxt.model_sxt_isg.unet import UNet
        
        logger.step("Starting Mitochondria model training...")
        
        # Setup device and model
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        logger.step(f"Using device: {device}")
        
        # Create UNet model for Mito (2 classes: background=0, mito=1)
        net = UNet(n_channels=3, n_classes=2, bilinear=False)
        net.to(device)
        logger.step("Mito Model created with 2 classes (background, mitochondria)")
        
        # Training configuration with optimized settings
        training_config = {
            'epochs': 25,
            'batch_size': 12,
            'learning_rate': 2e-4,
            'val_percent': 0.15,
            'img_scale': 0.5,
            'amp': True,
            'custom_img_dir': training_img_dir,
            'custom_mask_dir': training_mask_dir,
            'use_wandb': False,
            'allowed_ids': training_ids,
            'sample_ratio': 0.8,
            'organelle': 'mito'
        }
        
        logger.step(f"Training config: {training_config}")
        logger.step("Starting Mitochondria training...")
        
        # Start training
        train_net(
            net=net,
            device=device,
            external_logger=logger,  # Pass logger to train_net
            **training_config
        )
        
        training_time = time.time() - training_start
        logger.step(f"Mito training completed successfully in {training_time:.2f}s")
        
    except Exception as e:
        training_time = time.time() - training_start
        logger.step(f"Training failed after {training_time:.2f}s: {str(e)}")
        import traceback
        logger.step(f"Error details: {traceback.format_exc()}")

def train_pm_ne_model():
    """Train PM (cytosol) and NE (nucleus) segmentation model using 3D data"""
    
    training_start = time.time()
    mainpath = get_args().main_path

    # Define all IDs to use for training (mixed together)
    all_ids = ['784_5', '766_8', '842_17', '783_5', '766_5', '842_12', '766_2', '766_7', 
               '766_10', '766_11', '769_5', '769_7', '783_6', '783_12', '784_4', '784_6', 
               '784_7', '785_7', '822_4', '822_6', '822_7', '842_13', '931_9', '931_14']

    # Initialize logger
    log_dir = f'{mainpath}/logs'
    logger = QuickLogger("sxt_pm_ne_training", log_dir=log_dir)
    logger.step("Starting SXT PM+NE Model Training")

    logger.step(f"Using all {len(all_ids)} IDs for training (mixed together):")
    logger.step(f"IDs: {all_ids}")

    # Data paths - use pre-generated combined masks
    training_img_dir = r'D:\Gitspace\ipa_full\data\SXT\for_24_datasets'
    training_mask_dir = r'D:\Gitspace\ipa_full\data\SXT\tiff_24_datasets_pm_ne'  # Pre-generated combined masks

    logger.step(f"Image directory: {training_img_dir}")
    logger.step(f"Combined PM+NE mask directory: {training_mask_dir}")
    
    # Verify data directories
    if not os.path.exists(training_img_dir) or not os.path.exists(training_mask_dir):
        logger.step("ERROR: Training data directories not found!")
        logger.step("Please run 'generate_pm_ne_masks.py' first to generate combined masks!")
        return

    # Validate data pairs using existing function (add pm_ne to the validation patterns)
    training_ids, valid_pairs = validate_data_pairs(
        training_img_dir, training_mask_dir, all_ids, 'pm_ne', logger
    )
    
    logger.step(f"Using {len(training_ids)} available IDs for training: {training_ids}")
    
    if len(training_ids) == 0:
        logger.step("ERROR: No valid training data found!")
        logger.step("Please run 'generate_pm_ne_masks.py' first to generate combined masks!")
        return

    try:
        # Import training modules
        from ipa.processing.segmentation.segmentation_sxt.model_sxt_isg.train import train_net
        from ipa.processing.segmentation.segmentation_sxt.model_sxt_isg.unet import UNet
        
        logger.step("Starting PM+NE model training...")
        
        # Setup device and model
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        logger.step(f"Using device: {device}")
        
        # Create UNet model for PM+NE (3 classes: background=0, cytosol=1, nucleus=2)
        net = UNet(n_channels=3, n_classes=3, bilinear=False)
        net.to(device)
        logger.step("PM+NE Model created with 3 classes (background, cytosol, nucleus)")
        
        # Training configuration with optimized settings
        training_config = {
            'epochs': 2,
            'batch_size': 8,   # Smaller batch size due to 3 classes
            'learning_rate': 2e-4,
            'val_percent': 0.15,
            'img_scale': 0.5,
            'amp': True,
            'custom_img_dir': training_img_dir,
            'custom_mask_dir': training_mask_dir,
            'use_wandb': False,
            'allowed_ids': training_ids,
            'sample_ratio': 0.04,
            'organelle': 'pm_ne'
        }
        
        logger.step(f"Training config: {training_config}")
        logger.step("Starting PM+NE training...")
        
        # Start training
        train_net(
            net=net,
            device=device,
            external_logger=logger,
            **training_config
        )
        
        training_time = time.time() - training_start
        logger.step(f"PM+NE training completed successfully in {training_time:.2f}s")
        
    except Exception as e:
        training_time = time.time() - training_start
        logger.step(f"Training failed after {training_time:.2f}s: {str(e)}")
        import traceback
        logger.step(f"Error details: {traceback.format_exc()}")


def main():
    """Main training script with model selection"""
    print("SXT Training Demo")
    print("Available models:")
    print("1. ISG (Insulin Secretory Granules) - 2 classes")
    print("2. Mitochondria segmentation - 2 classes")
    print("3. PM + NE segmentation - 3 classes (cytosol=1, nucleus=2)")
    print()
    
    # Add simple command line model selection
    import sys

    # Train PM+NE model
    print("\nTraining PM+NE model...")
    train_pm_ne_model()
    # train_isg_model()
    # train_mito_model()


if __name__ == "__main__":
    main()
