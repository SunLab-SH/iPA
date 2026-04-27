#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
SIM Cell Segmentation Training - PM and NE Segmentation Model Training

This script provides a production-ready interface for training deep learning models
for cell segmentation on SIM (Structured Illumination Microscopy) imaging data.
It includes data preparation, model training, validation, and checkpoint management.

Supported operations:
- Data preparation: Combine PM and NE channels from raw SIM data
- Model training: Train UNet model for 3-class segmentation (background, PM, NE)
- Validation: Track metrics including loss, accuracy, and AUC
- Checkpoint saving: Save model checkpoints and training curves

Usage:
    # Run data preparation and training
    python demo_SIM_cell_seg_train.py --main_path /path/to/ipa --run_data_prep --run_training
    
    # Training only (if data already prepared)
    python demo_SIM_cell_seg_train.py --main_path /path/to/ipa --run_training --epochs 30 --batch_size 16
    
    # Full pipeline with custom parameters
    python demo_SIM_cell_seg_train.py --main_path /path/to/ipa --run_data_prep --run_training --epochs 20 --learning_rate 0.0001 --visualization
"""

# ═══════════════════════════════════════════════════════════════════
# Imports
# ═══════════════════════════════════════════════════════════════════
import os
import sys
import json
import time
import glob
import argparse
import warnings
import numpy as np
import torch
import matplotlib
matplotlib.use('Agg')  # Server-friendly backend

from pathlib import Path

# Filter MRC file warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, module="mrcfile")
warnings.filterwarnings("ignore", message=".*Unrecognised machine stamp.*")
warnings.filterwarnings("ignore", message=".*Map ID string not found.*")
warnings.filterwarnings("ignore", message=".*data block cannot be read.*")

from ipa.data_loader import QuickLogger


# ═══════════════════════════════════════════════════════════════════
# Argument Parsing
# ═══════════════════════════════════════════════════════════════════
def parse_args():
    """Parse command line arguments for SIM cell segmentation training"""
    p = argparse.ArgumentParser(
        description="SIM Cell Segmentation Training",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    p.add_argument("--main_path", type=str, required=True,
                   help="Project main path (containing iPA directory)")
    
    # Optional arguments
    p.add_argument("--raw_img_dir", type=str, 
                   default=r'F:\Salilab\Projects\IPA-toolbox\raw_image\SIM_data\Singlecolor',
                   help="Directory containing raw SIM images")
    p.add_argument("--mask_dir", type=str,
                   default=r'F:\Salilab\Projects\IPA-toolbox\raw_image\SIM_data\from_bing\SIM\output_mask',
                   help="Directory containing mask files")
    p.add_argument("--output_dir", type=str, default=None,
                   help="Output directory (default: {main_path}/outputs/sim_training)")
    p.add_argument("--seed", type=int, default=42,
                   help="Random seed for reproducibility")
    
    # Analysis selection flags
    p.add_argument("--run_data_prep", action="store_true",
                   help="Run data preparation (combine PM and NE channels)")
    p.add_argument("--run_training", action="store_true",
                   help="Run model training")
    p.add_argument("--force_reprocess", action="store_true",
                   help="Force reprocessing of existing data")
    
    # Training parameters
    p.add_argument("--epochs", type=int, default=20,
                   help="Number of training epochs")
    p.add_argument("--batch_size", type=int, default=8,
                   help="Batch size for training")
    p.add_argument("--learning_rate", type=float, default=1e-4,
                   help="Learning rate")
    p.add_argument("--val_percent", type=float, default=0.2,
                   help="Validation set percentage")
    p.add_argument("--img_scale", type=float, default=0.5,
                   help="Image scale factor")
    p.add_argument("--data_usage_ratio", type=float, default=0.8,
                   help="Ratio of data to use for training")
    
    # Visualization and performance
    p.add_argument("--visualization", action="store_true",
                   help="Enable visualization output")
    p.add_argument("--num_workers", type=int, default=2,
                   help="Number of data loader workers")
    
    return p.parse_args()


# ═══════════════════════════════════════════════════════════════════
# Helper Functions
# ═══════════════════════════════════════════════════════════════════

def set_seed(seed: int):
    """Set random seed for reproducibility"""
    import random
    import numpy as np
    
    random.seed(seed)
    np.random.seed(seed)
    
    try:
        import torch
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    except ImportError:
        pass


def ensure_outdir(base: str) -> str:
    """Create timestamped output directory"""
    os.makedirs(base, exist_ok=True)
    run_id = time.strftime("%Y%m%d-%H%M%S")
    outdir = os.path.join(base, run_id)
    os.makedirs(outdir, exist_ok=True)
    return outdir


def save_json(outdir: str, name: str, obj: dict) -> str:
    """Save dict as JSON file"""
    fp = os.path.join(outdir, f"{name}.json")
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)
    return fp


def check_exists(logger: QuickLogger, **kwargs):
    """Check if required files exist"""
    for k, v in kwargs.items():
        if v and not os.path.exists(v):
            logger.step(f"Warning: Missing {k}: {v}")


# ═══════════════════════════════════════════════════════════════════
# Data Preparation Functions
# ═══════════════════════════════════════════════════════════════════

def plot_sim_data_visualization(pm_img, ne_img, combined_img, combined_mask, file_id, output_dir, logger):
    """Plot visualization of SIM data processing"""
    import matplotlib.pyplot as plt
    
    try:
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        
        # Enhance and plot PM channel
        pm_display = (pm_img - pm_img.min()) / (pm_img.max() - pm_img.min()) if pm_img.max() > pm_img.min() else pm_img
        pm_display = np.power(pm_display, 0.5)
        axes[0, 0].imshow(pm_display, cmap='viridis')
        axes[0, 0].set_title(f'PM Channel\n{file_id}')
        axes[0, 0].axis('off')
        
        # Enhance and plot NE channel
        ne_display = (ne_img - ne_img.min()) / (ne_img.max() - ne_img.min()) if ne_img.max() > ne_img.min() else ne_img
        ne_display = np.power(ne_display, 0.5)
        axes[0, 1].imshow(ne_display, cmap='magma')
        axes[0, 1].set_title(f'NE Channel\n{file_id}')
        axes[0, 1].axis('off')
        
        # Plot combined RGB
        rgb_display = combined_img.astype(np.float32) / 255.0 if combined_img.max() > 1 else combined_img
        for i in range(3):
            if rgb_display[:,:,i].max() > rgb_display[:,:,i].min():
                channel = rgb_display[:,:,i]
                channel = (channel - channel.min()) / (channel.max() - channel.min())
                rgb_display[:,:,i] = np.power(channel, 0.5)
        axes[0, 2].imshow(rgb_display)
        axes[0, 2].set_title('Combined RGB (PM=Red, NE=Green)')
        axes[0, 2].axis('off')
        
        # Plot overlays
        rgb_for_overlay = combined_img / 255.0 if combined_img.max() > 1 else combined_img
        
        axes[1, 0].imshow(rgb_for_overlay)
        pm_overlay = np.zeros((*combined_mask.shape, 4))
        pm_overlay[combined_mask == 1] = [1, 0, 0, 0.5]
        axes[1, 0].imshow(pm_overlay)
        axes[1, 0].set_title('PM Mask Overlay')
        axes[1, 0].axis('off')
        
        axes[1, 1].imshow(rgb_for_overlay)
        ne_overlay = np.zeros((*combined_mask.shape, 4))
        ne_overlay[combined_mask == 2] = [0, 0, 1, 0.5]
        axes[1, 1].imshow(ne_overlay)
        axes[1, 1].set_title('NE Mask Overlay')
        axes[1, 1].axis('off')
        
        # Plot combined mask
        mask_display = np.zeros((*combined_mask.shape, 3), dtype=np.uint8)
        mask_display[combined_mask == 1] = [255, 0, 0]
        mask_display[combined_mask == 2] = [0, 0, 255]
        axes[1, 2].imshow(mask_display)
        axes[1, 2].set_title('Combined Mask (PM=Red, NE=Blue)')
        axes[1, 2].axis('off')
        
        # Add statistics
        pm_pixels = np.sum(combined_mask == 1)
        ne_pixels = np.sum(combined_mask == 2)
        total_pixels = combined_mask.size
        stats_text = f"""PM: {pm_pixels} ({pm_pixels/total_pixels*100:.1f}%), NE: {ne_pixels} ({ne_pixels/total_pixels*100:.1f}%)
PM range: {pm_img.min():.3f}-{pm_img.max():.3f}, NE range: {ne_img.min():.3f}-{ne_img.max():.3f}"""
        
        fig.text(0.02, 0.02, stats_text, fontsize=10, verticalalignment='bottom',
                bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
        
        plt.tight_layout()
        viz_path = os.path.join(output_dir, f'visualization_{file_id.strip("_")}.png')
        plt.savefig(viz_path, dpi=150, bbox_inches='tight')
        plt.close()
        
        logger.step(f"Visualization saved: {viz_path}")
        
    except Exception as e:
        logger.step(f"Warning: Visualization failed for {file_id}: {str(e)}")


def prepare_sim_data(base_date, file_ids, raw_dir, mask_dir, output_dir, force_reprocess, logger):
    """Prepare SIM training data by combining PM and NE channels"""
    from ipa.data_loader import UniversalDataLoader
    import tifffile
    
    combined_img_dir = os.path.join(output_dir, 'combined_images')
    combined_mask_dir = os.path.join(output_dir, 'combined_masks')
    
    # Check if data already processed
    if not force_reprocess and os.path.exists(combined_img_dir):
        existing_imgs = glob.glob(os.path.join(combined_img_dir, '*.tiff'))
        existing_masks = glob.glob(os.path.join(combined_mask_dir, '*_label.tiff'))
        
        if len(existing_imgs) > 0 and len(existing_masks) > 0:
            logger.step(f"Found {len(existing_imgs)} processed images, {len(existing_masks)} masks - skipping preparation")
            processed_ids = [os.path.basename(f).replace(f'{base_date}_', '').replace('.tiff', '') 
                           for f in existing_imgs if f'{base_date}_' in os.path.basename(f)]
            return processed_ids, combined_img_dir, combined_mask_dir
    
    os.makedirs(combined_img_dir, exist_ok=True)
    os.makedirs(combined_mask_dir, exist_ok=True)
    
    processed_ids = []
    
    for file_id in file_ids:
        try:
            # Build file paths
            pm_img_path = os.path.join(raw_dir, f"{base_date}{file_id}_Actin.tif")
            ne_img_path = os.path.join(raw_dir, f"{base_date}{file_id}_N.tif")
            pm_mask_path = os.path.join(mask_dir, f"{base_date}{file_id}_PM.mrc")
            ne_mask_path = os.path.join(mask_dir, f"{base_date}{file_id}_N.mrc")
            
            # Check all files exist
            if not all(os.path.exists(p) for p in [pm_img_path, ne_img_path, pm_mask_path, ne_mask_path]):
                missing = [os.path.basename(p) for p in [pm_img_path, ne_img_path, pm_mask_path, ne_mask_path] 
                          if not os.path.exists(p)]
                logger.step(f"Skip {file_id} - missing: {missing}")
                continue
            
            # Load data
            pm_img = UniversalDataLoader.load_data(pm_img_path)
            ne_img = UniversalDataLoader.load_data(ne_img_path)
            pm_mask = UniversalDataLoader.load_data(pm_mask_path)
            ne_mask = UniversalDataLoader.load_data(ne_mask_path)
            
            # Validate shapes
            if pm_img.shape != ne_img.shape or pm_mask.shape != ne_mask.shape:
                logger.step(f"Skip {file_id} - shape mismatch")
                continue
            
            # Normalize images to 0-1
            pm_img_norm = (pm_img - pm_img.min()) / (pm_img.max() - pm_img.min()) if pm_img.max() > pm_img.min() else np.zeros_like(pm_img)
            ne_img_norm = (ne_img - ne_img.min()) / (ne_img.max() - ne_img.min()) if ne_img.max() > ne_img.min() else np.zeros_like(ne_img)
            
            # Combine as 3-channel RGB-like image
            combined_img = np.stack([
                (pm_img_norm * 255).astype(np.uint8),
                (ne_img_norm * 255).astype(np.uint8),
                np.zeros_like(pm_img, dtype=np.uint8)
            ], axis=-1)
            
            # Combine masks: PM=1, NE=2
            combined_mask = np.zeros_like(pm_mask, dtype=np.uint8)
            combined_mask[pm_mask > 0] = 1
            combined_mask[ne_mask > 0] = 2
            
            # Check for mask data
            z_with_data = np.any(combined_mask > 0, axis=(1,2))
            if not np.any(z_with_data):
                logger.step(f"Skip {file_id} - no mask data")
                continue
            
            # Save 3D volumes
            base_id = f"{base_date}{file_id}"
            img_output_path = os.path.join(combined_img_dir, f"{base_id}.tiff")
            mask_output_path = os.path.join(combined_mask_dir, f"{base_id}_label.tiff")
            
            tifffile.imwrite(img_output_path, combined_img)
            tifffile.imwrite(mask_output_path, combined_mask)
            
            processed_ids.append(f"{file_id.strip('_')}")
            
            # Create visualization for middle slice
            middle_slice = pm_img.shape[0] // 2 if len(pm_img.shape) > 2 else 0
            pm_vis = pm_img[middle_slice] if len(pm_img.shape) > 2 else pm_img
            ne_vis = ne_img[middle_slice] if len(ne_img.shape) > 2 else ne_img
            pm_mask_vis = pm_mask[middle_slice] if len(pm_mask.shape) > 2 else pm_mask
            ne_mask_vis = ne_mask[middle_slice] if len(ne_mask.shape) > 2 else ne_mask
            
            pm_vis_norm = (pm_vis - pm_vis.min()) / (pm_vis.max() - pm_vis.min()) if pm_vis.max() > pm_vis.min() else np.zeros_like(pm_vis)
            ne_vis_norm = (ne_vis - ne_vis.min()) / (ne_vis.max() - ne_vis.min()) if ne_vis.max() > ne_vis.min() else np.zeros_like(ne_vis)
            
            combined_img_vis = np.stack([
                (pm_vis_norm * 255).astype(np.uint8),
                (ne_vis_norm * 255).astype(np.uint8),
                np.zeros_like(pm_vis, dtype=np.uint8)
            ], axis=-1)
            
            combined_mask_vis = np.zeros_like(pm_mask_vis, dtype=np.uint8)
            combined_mask_vis[pm_mask_vis > 0] = 1
            combined_mask_vis[ne_mask_vis > 0] = 2
            
            plot_sim_data_visualization(pm_vis_norm, ne_vis_norm, combined_img_vis, 
                                       combined_mask_vis, file_id, output_dir, logger)
            
        except Exception as e:
            logger.step(f"Error processing {file_id}: {str(e)}")
            continue
    
    logger.step(f"Processed {len(processed_ids)} data pairs")
    return processed_ids, combined_img_dir, combined_mask_dir


def validate_sim_data(img_dir, mask_dir, logger):
    """Validate SIM image-mask pairs"""
    img_files = glob.glob(os.path.join(img_dir, '*.tiff'))
    mask_files = glob.glob(os.path.join(mask_dir, '*_label.tiff'))
    
    valid_pairs = {}
    available_ids = set()
    
    for img_file in img_files:
        filename = os.path.basename(img_file)
        if filename.endswith('.tiff') and '20220909_' in filename:
            file_id = filename.replace('20220909_', '').replace('.tiff', '')
            mask_path = os.path.join(mask_dir, f"20220909_{file_id}_label.tiff")
            
            if os.path.exists(mask_path):
                available_ids.add(file_id)
                valid_pairs[file_id] = {'img': img_file, 'mask': mask_path}
    
    logger.step(f"Found {len(valid_pairs)} valid image-mask pairs")
    return list(available_ids), valid_pairs


def train_model(img_dir, mask_dir, training_ids, config, outdir, logger):
    """Train PM+NE segmentation model"""
    from ipa.processing.segmentation import UNet_sim as UNet
    from ipa.processing.segmentation import create_sim_data_loaders
    import torch.nn as nn
    import torch.optim as optim
    from torch.cuda.amp import autocast, GradScaler
    
    # Setup device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logger.step(f"Using device: {device}")
    
    # Create model
    net = UNet(n_channels=3, n_classes=3, bilinear=False)
    net.to(device)
    
    # Create data loaders
    train_loader, val_loader = create_sim_data_loaders(
        img_dir=img_dir,
        mask_dir=mask_dir,
        allowed_ids=[f"20220909_{id}" for id in training_ids],
        batch_size=config['batch_size'],
        val_percent=config['val_percent'],
        img_scale=config['img_scale'],
        data_usage_ratio=config['data_usage_ratio'],
        num_workers=config['num_workers']
    )
    
    logger.step(f"Data loaders: {len(train_loader)} train batches, {len(val_loader)} val batches")
    
    # Setup training
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(net.parameters(), lr=config['learning_rate'])
    scaler = GradScaler() if config['amp'] else None
    
    checkpoint_dir = os.path.join(outdir, 'checkpoints')
    os.makedirs(checkpoint_dir, exist_ok=True)
    
    best_val_loss = float('inf')
    train_losses, val_losses = [], []
    
    # Training loop
    for epoch in range(config['epochs']):
        epoch_start = time.time()
        
        # Training phase
        net.train()
        epoch_loss = 0
        for batch_idx, batch in enumerate(train_loader):
            images = batch['image'].to(device)
            masks = batch['mask'].to(device)
            
            optimizer.zero_grad()
            
            if config['amp'] and scaler:
                with autocast():
                    outputs = net(images)
                    loss = criterion(outputs, masks)
                scaler.scale(loss).backward()
                scaler.step(optimizer)
                scaler.update()
            else:
                outputs = net(images)
                loss = criterion(outputs, masks)
                loss.backward()
                optimizer.step()
            
            epoch_loss += loss.item()
            
            if batch_idx % 10 == 0:
                logger.step(f"Epoch {epoch+1}/{config['epochs']}, Batch {batch_idx}/{len(train_loader)}, Loss: {loss.item():.4f}")
        
        avg_train_loss = epoch_loss / len(train_loader)
        train_losses.append(avg_train_loss)
        
        # Validation phase
        net.eval()
        val_loss = 0
        correct = 0
        total = 0
        
        all_predictions = []
        all_targets = []
        
        with torch.no_grad():
            for batch in val_loader:
                images = batch['image'].to(device)
                masks = batch['mask'].to(device)
                
                outputs = net(images)
                loss = criterion(outputs, masks)
                val_loss += loss.item()
                
                predictions = torch.argmax(outputs, dim=1)
                correct += (predictions == masks).sum().item()
                total += masks.numel()
                
                probs = torch.softmax(outputs, dim=1)
                all_predictions.append(probs.cpu().numpy().reshape(-1, 3))
                all_targets.append(masks.cpu().numpy().reshape(-1))
        
        avg_val_loss = val_loss / len(val_loader)
        val_losses.append(avg_val_loss)
        accuracy = correct / total
        
        # Calculate AUC
        try:
            from sklearn.metrics import roc_auc_score
            all_preds = np.concatenate(all_predictions, axis=0)
            all_targs = np.concatenate(all_targets, axis=0)
            
            auc_scores = {}
            for class_idx, class_name in enumerate(['Background', 'PM', 'NE']):
                binary_targets = (all_targs == class_idx).astype(int)
                class_probs = all_preds[:, class_idx]
                if len(np.unique(binary_targets)) > 1:
                    auc_scores[class_name] = roc_auc_score(binary_targets, class_probs)
                else:
                    auc_scores[class_name] = 0.0
            macro_auc = np.mean(list(auc_scores.values()))
        except:
            auc_scores = {'Background': 0.0, 'PM': 0.0, 'NE': 0.0}
            macro_auc = 0.0
        
        epoch_time = time.time() - epoch_start
        logger.step(f"Epoch {epoch+1}/{config['epochs']}: Train Loss={avg_train_loss:.4f}, Val Loss={avg_val_loss:.4f}, Acc={accuracy:.4f}, AUC={macro_auc:.4f} ({epoch_time:.2f}s)")
        
        # Save checkpoint
        checkpoint_path = os.path.join(checkpoint_dir, f'epoch_{epoch+1}.pth')
        torch.save({
            'epoch': epoch + 1,
            'model_state_dict': net.state_dict(),
            'optimizer_state_dict': optimizer.state_dict(),
            'train_loss': avg_train_loss,
            'val_loss': avg_val_loss,
            'accuracy': accuracy,
            'macro_auc': macro_auc,
            'class_aucs': auc_scores,
        }, checkpoint_path)
        
        # Save best model
        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            best_model_path = os.path.join(outdir, 'best_model.pth')
            torch.save({
                'epoch': epoch + 1,
                'model_state_dict': net.state_dict(),
                'train_loss': avg_train_loss,
                'val_loss': avg_val_loss,
                'accuracy': accuracy,
                'macro_auc': macro_auc,
            }, best_model_path)
            logger.step(f"Best model saved (val_loss={avg_val_loss:.4f})")
    
    # Save final model
    final_model_path = os.path.join(outdir, 'final_model.pth')
    torch.save({
        'model_state_dict': net.state_dict(),
        'train_losses': train_losses,
        'val_losses': val_losses,
    }, final_model_path)
    
    # Plot training curves
    if config.get('visualization', False):
        try:
            import matplotlib.pyplot as plt
            plt.figure(figsize=(12, 4))
            plt.subplot(1, 2, 1)
            plt.plot(range(1, len(train_losses) + 1), train_losses, 'b-', label='Train')
            plt.plot(range(1, len(val_losses) + 1), val_losses, 'r-', label='Val')
            plt.xlabel('Epoch')
            plt.ylabel('Loss')
            plt.title('Training and Validation Loss')
            plt.legend()
            plt.grid(True)
            
            plt.subplot(1, 2, 2)
            plt.plot(range(1, len(val_losses) + 1), val_losses, 'r-')
            plt.xlabel('Epoch')
            plt.ylabel('Loss')
            plt.title('Validation Loss')
            plt.grid(True)
            
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, 'training_curves.png'), dpi=150)
            plt.close()
        except:
            pass
    
    return {
        'final_train_loss': train_losses[-1],
        'final_val_loss': val_losses[-1],
        'best_val_loss': best_val_loss,
        'train_losses': train_losses,
        'val_losses': val_losses
    }


# ═══════════════════════════════════════════════════════════════════
# Main Function
# ═══════════════════════════════════════════════════════════════════

def main():
    """Main entry point"""
    
    # Step 1: Initialize and parse arguments
    args = parse_args()
    set_seed(args.seed)
    
    if args.output_dir is None:
        args.output_dir = os.path.join(args.main_path, "outputs", "sim_training")
    
    outdir = ensure_outdir(args.output_dir)
    
    # Step 2: Initialize logger
    logger = QuickLogger(name="sim_cell_seg_training", log_dir=outdir)
    logger.step("SIM Cell Segmentation Training")
    
    # Step 3: Build paths and config
    base_date = "20220909"
    file_ids = [
        "_0-2-1-SIM", "_0-2-2-SIM", "_0-3-1-SIM", "_0-3-1-1-SIM", "_0-3-3-SIM",
        "_0-3-5-SIM", "_0-4-1-1-SIM", "_0-4-1-SIM", "_0-4-2-SIM",
        "_5-2-1-SIM", "_5-2-2-SIM", "_5-2-3-1-SIM", "_5-2-3-SIM",
        "_5-3-2-SIM", "_5-3-3-SIM", "_5-3-4-1-SIM", "_5-3-4-SIM",
        "_30-2-1-SIM", "_30-2-2-SIM", "_30-2-3-SIM", "_30-2-4-SIM",
        "_30-3-2-SIM", "_30-3-3-SIM", "_30-4-1-SIM", "_30-4-3-SIM",
        "_30-4-4-SIM", "_30-4-5-SIM"
    ]
    
    logger.step(f"Main path: {args.main_path}")
    logger.step(f"Output: {outdir}")
    logger.step(f"Processing {len(file_ids)} file IDs")
    
    training_config = {
        'epochs': args.epochs,
        'batch_size': args.batch_size,
        'learning_rate': args.learning_rate,
        'val_percent': args.val_percent,
        'img_scale': args.img_scale,
        'amp': True,
        'data_usage_ratio': args.data_usage_ratio,
        'num_workers': args.num_workers,
        'visualization': args.visualization
    }
    
    # Step 4: Save configuration
    config_data = {
        'training_config': training_config,
        'file_ids': file_ids,
        'base_date': base_date,
        'raw_img_dir': args.raw_img_dir,
        'mask_dir': args.mask_dir,
        'seed': args.seed,
        'timestamp': time.strftime("%Y-%m-%d %H:%M:%S")
    }
    save_json(outdir, "config", config_data)
    
    # Step 5: Validate input directories
    check_exists(logger, raw_img_dir=args.raw_img_dir, mask_dir=args.mask_dir)
    
    # Step 6: Define execution wrapper
    def timed_run(name: str, fn, **kwargs):
        t0 = time.time()
        try:
            logger.step(f"Starting {name}...")
            result = fn(**kwargs)
            dt = time.time() - t0
            logger.step(f"{name} completed in {dt:.2f}s")
            return result
        except FileNotFoundError as e:
            logger.step(f"{name} failed: file not found: {e}")
            return None
        except MemoryError:
            logger.step(f"{name} failed: out of memory")
            return None
        except Exception as e:
            logger.step(f"{name} failed: {type(e).__name__}: {e}")
            return None
    
    # Step 7: Conditional execution
    combined_img_dir = None
    combined_mask_dir = None
    
    # Data preparation
    if args.run_data_prep:
        result = timed_run(
            "data_preparation",
            prepare_sim_data,
            base_date=base_date,
            file_ids=file_ids,
            raw_dir=args.raw_img_dir,
            mask_dir=args.mask_dir,
            output_dir=outdir,
            force_reprocess=args.force_reprocess,
            logger=logger
        )
        if result:
            processed_ids, combined_img_dir, combined_mask_dir = result
    
    # Training
    if args.run_training:
        # Determine data directories
        if combined_img_dir is None:
            combined_img_dir = os.path.join(outdir, 'combined_images')
            combined_mask_dir = os.path.join(outdir, 'combined_masks')
        
        if not os.path.exists(combined_img_dir) or not os.path.exists(combined_mask_dir):
            logger.step("Error: Training data not found. Run with --run_data_prep first")
        else:
            # Validate data
            training_ids, _ = validate_sim_data(combined_img_dir, combined_mask_dir, logger)
            
            if len(training_ids) > 0:
                result = timed_run(
                    "model_training",
                    train_model,
                    img_dir=combined_img_dir,
                    mask_dir=combined_mask_dir,
                    training_ids=training_ids,
                    config=training_config,
                    outdir=outdir,
                    logger=logger
                )
                if result:
                    save_json(outdir, "training_results", result)
            else:
                logger.step("Error: No valid training data found")
    
    # Step 8: Completion
    logger.step(f"All operations completed. Results: {outdir}")
    print(f"\nOutputs saved to: {outdir}")


# ═══════════════════════════════════════════════════════════════════
# Entry Point
# ═══════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    main()
