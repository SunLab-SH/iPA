#%%
import argparse
import logging
import sys
from pathlib import Path
import os
import glob
import json
from datetime import datetime
import time
import warnings

import torch
import torch.nn as nn
import torch.nn.functional as F
import wandb
from torch import optim
from torch.utils.data import DataLoader, random_split, Dataset
from tqdm import tqdm
import numpy as np
import tifffile
import mrcfile
from PIL import Image

# Filter MRC file warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, module="mrcfile")
warnings.filterwarnings("ignore", message=".*Unrecognised machine stamp.*")
warnings.filterwarnings("ignore", message=".*Map ID string not found.*")
warnings.filterwarnings("ignore", message=".*data block cannot be read.*")

from .utils.data_loading import BasicDataset, CarvanaDataset
from .utils.dice_score import dice_loss
from .evaluate import evaluate
from .unet import UNet

curpath = os.path.dirname(os.path.abspath(__file__))
os.chdir(curpath)
os.path.pardir

sys.path.append(os.getcwd())
parpath = os.getcwd()
# print(parpath)

# Updated paths for ISG segmentation training data
dir_img = Path(f'./data/isg_training/imgs/')
dir_mask = Path(f'./data/isg_training/masks/')
dir_checkpoint = Path(f'./checkpoints/isg_models/')

class ISG3DDataset(Dataset):
    """Dataset class for handling 3D data and extracting 2D slices for training"""
    
    def __init__(self, img_dir, mask_dir, scale=0.5, allowed_ids=None, sample_ratio=1.0):
        dataset_load_start = time.time()
        logging.info("Starting ISG dataset loading...")
        
        self.img_dir = Path(img_dir)
        self.mask_dir = Path(mask_dir)
        self.scale = scale
        self.slice_info = []
        self.file_slice_count = {}
        self.allowed_ids = allowed_ids
        self.sample_ratio = sample_ratio
        
        # Find all 3D image files
        img_files = list(self.img_dir.glob('*.tif')) + list(self.img_dir.glob('*.tiff')) + list(self.img_dir.glob('*.mrc'))
        mask_files = list(self.mask_dir.glob('*.tif')) + list(self.mask_dir.glob('*.tiff')) + list(self.mask_dir.glob('*.mrc'))
        
        logging.info(f"Found {len(img_files)} image files and {len(mask_files)} mask files")
        
        # Filter by allowed IDs if provided
        if allowed_ids:
            logging.info(f"Filtering files by allowed IDs: {allowed_ids}")
            filtered_img_files = []
            for img_file in img_files:
                img_identifier = self._extract_identifier_from_img(img_file.name)
                if img_identifier in allowed_ids:
                    filtered_img_files.append(img_file)
            img_files = filtered_img_files
            logging.info(f"After filtering: {len(img_files)} image files")
        
        successful_matches = 0
        
        for img_file in img_files:
            # Extract identifier from image filename
            img_identifier = self._extract_identifier_from_img(img_file.name)
            mask_file = self._find_matching_mask(img_identifier, mask_files)
            
            if mask_file is None:
                logging.warning(f"No mask found for {img_file.name} (identifier: {img_identifier})")
                continue
            
            # Load and validate 3D data
            try:
                img_3d, mask_3d = self._load_and_validate_data(img_file, mask_file)
                if img_3d is None or mask_3d is None:
                    logging.warning(f"Skipping {img_file.name} due to data loading failure")
                    continue
                
                # Find valid mask range (frames with actual mask data)
                mask_sums = np.sum(mask_3d, axis=(1, 2))  # Sum over height and width for each frame
                valid_mask_frames = np.where(mask_sums > 0)[0]  # Frames with non-zero masks
                
                if len(valid_mask_frames) == 0:
                    logging.warning(f"No valid mask frames found in {mask_file.name}")
                    continue
                
                # Expand sampling range by ±5 frames around valid mask frames
                min_valid_frame = max(0, int(valid_mask_frames.min()) - 5)
                max_valid_frame = min(img_3d.shape[0] - 1, int(valid_mask_frames.max()) + 5)
                
                total_available_slices = max_valid_frame - min_valid_frame + 1
                
                logging.info(f"Valid mask frames for {img_file.name}: {int(valid_mask_frames.min())}-{int(valid_mask_frames.max())}, "
                           f"Expanded sampling range: {min_valid_frame}-{max_valid_frame} ({total_available_slices} slices)")
                
                # Store slice information
                self.file_slice_count[img_file.name] = total_available_slices
                
                # Apply sampling ratio - only use a portion of slices from the expanded range
                if sample_ratio < 1.0:
                    num_sampled_slices = max(1, int(total_available_slices * sample_ratio))
                    # Sample evenly across the expanded range
                    slice_indices = np.linspace(min_valid_frame, max_valid_frame, num_sampled_slices, dtype=int)
                    logging.info(f"Sampling {num_sampled_slices}/{total_available_slices} slices from expanded range")
                else:
                    slice_indices = range(min_valid_frame, max_valid_frame + 1)
                
                for slice_idx in slice_indices:
                    self.slice_info.append({
                        'img_file': str(img_file),
                        'mask_file': str(mask_file),
                        'slice_idx': int(slice_idx),  # Convert to Python int
                        'source_file': img_file.name,
                        'identifier': img_identifier,
                        'global_idx': len(self.slice_info),
                        'is_valid_mask_frame': int(slice_idx) in valid_mask_frames
                    })
                
                logging.info(f"✓ {img_file.name}: {len(slice_indices)}/{total_available_slices} slices (range: {min_valid_frame}-{max_valid_frame}) -> {mask_file.name}")
                successful_matches += 1
                
            except Exception as e:
                logging.error(f"Error loading {img_file}: {e}")
                continue
        
        dataset_load_time = time.time() - dataset_load_start
        logging.info(f"ISG dataset loading completed in {dataset_load_time:.2f}s")
        logging.info(f"Successfully matched {successful_matches}/{len(img_files)} files")
        logging.info(f"Total slices: {len(self.slice_info)}")
        
        if len(self.slice_info) == 0:
            raise ValueError("No valid ISG training data found. Check file naming and paths.")
    
    def _extract_identifier_from_img(self, filename):
        """Extract identifier for matching with improved logic"""
        base = filename.rsplit('.', 1)[0]
        parts = base.split('_')
        
        # Look for two consecutive numeric parts
        for i in range(len(parts) - 1):
            if parts[i].isdigit() and parts[i+1].isdigit():
                return f"{parts[i]}_{parts[i+1]}"
        
        # Fallback: look for last two numeric parts before 'pre'
        if 'pre' in parts:
            pre_idx = parts.index('pre')
            if pre_idx >= 2 and parts[pre_idx-2].isdigit() and parts[pre_idx-1].isdigit():
                return f"{parts[pre_idx-2]}_{parts[pre_idx-1]}"
        
        # Additional fallback: look for pattern like "3_842"
        for i in range(len(parts)):
            if parts[i].isdigit():
                for j in range(i+1, len(parts)):
                    if parts[j].isdigit():
                        return f"{parts[i]}_{parts[j]}"
        
        return None
    
    def _find_matching_mask(self, img_identifier, mask_files):
        """Find matching mask file with improved matching"""
        if img_identifier is None:
            return None
        
        # Try multiple matching patterns
        patterns_to_try = [
            f"{img_identifier}_isg_label.tiff",
            f"{img_identifier}_isg_label.tif", 
            f"{img_identifier}.tiff",
            f"{img_identifier}.tif",
            f"{img_identifier}_label.tiff",
            f"{img_identifier}_label.tif",
            f"label_{img_identifier}.tiff",
            f"label_{img_identifier}.tif"
        ]
        
        for mask_file in mask_files:
            mask_name = mask_file.name
            
            # Direct pattern matching
            for pattern in patterns_to_try:
                if mask_name == pattern:
                    return mask_file
            
            # Partial matching - if identifier is contained in mask filename
            if img_identifier in mask_name and ('isg' in mask_name.lower() or 'label' in mask_name.lower()):
                return mask_file
        
        return None
    
    def _load_and_validate_data(self, img_file, mask_file):
        """Load and validate image and mask data with better error handling"""
        try:
            # Load image with error handling
            img_3d = None
            if str(img_file).endswith(('.tif', '.tiff')):
                try:
                    img_3d = tifffile.imread(str(img_file))
                except Exception as e:
                    logging.error(f"Failed to load TIFF image {img_file}: {e}")
                    return None, None
            else:
                try:
                    with mrcfile.open(str(img_file), permissive=True) as mrc:
                        if mrc.data is not None:
                            img_3d = mrc.data.copy()
                        else:
                            logging.error(f"MRC file {img_file} has no data")
                            return None, None
                except Exception as e:
                    logging.error(f"Failed to load MRC image {img_file}: {e}")
                    return None, None
            
            # Load mask with error handling
            mask_3d = None
            if str(mask_file).endswith(('.tif', '.tiff')):
                try:
                    mask_3d = tifffile.imread(str(mask_file))
                except Exception as e:
                    logging.error(f"Failed to load TIFF mask {mask_file}: {e}")
                    return None, None
            else:
                try:
                    with mrcfile.open(str(mask_file), permissive=True) as mrc:
                        if mrc.data is not None:
                            mask_3d = mrc.data.copy()
                        else:
                            logging.error(f"MRC file {mask_file} has no data")
                            return None, None
                except Exception as e:
                    logging.error(f"Failed to load MRC mask {mask_file}: {e}")
                    return None, None
            
            # Validate loaded data
            if img_3d is None or mask_3d is None:
                logging.error(f"One or both files failed to load: {img_file}, {mask_file}")
                return None, None
            
            # Handle 2D files
            if len(img_3d.shape) == 2:
                img_3d = img_3d[np.newaxis, ...]
            if len(mask_3d.shape) == 2:
                mask_3d = mask_3d[np.newaxis, ...]
            
            # Check dimension match
            if img_3d.shape[0] != mask_3d.shape[0]:
                logging.error(f"Slice count mismatch: {img_3d.shape[0]} vs {mask_3d.shape[0]}")
                return None, None
            
            return img_3d, mask_3d
            
        except Exception as e:
            logging.error(f"Unexpected error loading {img_file} or {mask_file}: {e}")
            return None, None
    
    def __len__(self):
        return len(self.slice_info)
    
    def __getitem__(self, idx):
        slice_data = self.slice_info[idx]
        
        try:
            # Load 3D data with better error handling and correct format detection
            img_3d = None
            mask_3d = None
            
            # Check image file format
            if slice_data['img_file'].endswith(('.tif', '.tiff')):
                img_3d = tifffile.imread(slice_data['img_file'])
            else:
                # MRC file
                try:
                    with mrcfile.open(slice_data['img_file'], permissive=True) as mrc:
                        if mrc.data is not None:
                            img_3d = mrc.data.copy()
                        else:
                            raise ValueError("MRC image data is None")
                except Exception as e:
                    logging.error(f"Failed to load MRC image {slice_data['img_file']}: {e}")
                    raise
            
            # Check mask file format - fix the bug here
            if slice_data['mask_file'].endswith(('.tif', '.tiff')):
                # TIFF mask file
                try:
                    mask_3d = tifffile.imread(slice_data['mask_file'])
                except Exception as e:
                    logging.error(f"Failed to load TIFF mask {slice_data['mask_file']}: {e}")
                    raise
            else:
                # MRC mask file
                try:
                    with mrcfile.open(slice_data['mask_file'], permissive=True) as mrc:
                        if mrc.data is not None:
                            mask_3d = mrc.data.copy()
                        else:
                            raise ValueError("MRC mask data is None")
                except Exception as e:
                    logging.error(f"Failed to load MRC mask {slice_data['mask_file']}: {e}")
                    raise
            
            # Validate data
            if img_3d is None or mask_3d is None:
                raise ValueError(f"Failed to load data for {slice_data['source_file']}")
            
            # Handle dimensions
            if len(img_3d.shape) == 2:
                img_3d = img_3d[np.newaxis, ...]
            if len(mask_3d.shape) == 2:
                mask_3d = mask_3d[np.newaxis, ...]
            
            # Extract slice
            slice_idx = slice_data['slice_idx']
            
            # Check bounds
            if slice_idx >= img_3d.shape[0]:
                raise IndexError(f"Image slice index {slice_idx} out of bounds (max: {img_3d.shape[0]-1})")
            
            img_slice = img_3d[slice_idx]
            
            # Handle mask slice - if outside mask range, create empty mask
            if slice_idx >= mask_3d.shape[0]:
                # Create empty mask with same spatial dimensions as image
                mask_slice = np.zeros_like(img_slice, dtype=np.uint8)
                logging.debug(f"Using empty mask for slice {slice_idx} (beyond mask range)")
            else:
                mask_slice = mask_3d[slice_idx]
            
            # Convert using ISG-specific processing
            img_slice = self._convert_img_slice(img_slice)
            mask_slice = self._convert_mask_slice(mask_slice)
            
            # Convert to PIL and then tensor
            img_rgb = np.stack([img_slice]*3, axis=-1).astype(np.uint8)
            pil_img = Image.fromarray(img_rgb)
            pil_mask = Image.fromarray(mask_slice.astype(np.uint8))
            
            # Create writable copies to avoid PyTorch warning
            img_tensor = torch.from_numpy(BasicDataset.preprocess(pil_img, self.scale, is_mask=False).copy())
            mask_tensor = torch.from_numpy(BasicDataset.preprocess(pil_mask, self.scale, is_mask=True).copy())
            
            return {
                'image': img_tensor,
                'mask': mask_tensor.long().squeeze()
            }
            
        except Exception as e:
            logging.error(f"Error in __getitem__ for {slice_data['source_file']}, slice {slice_data['slice_idx']}: {e}")
            # Return a dummy sample to avoid crashing the training
            dummy_img = torch.zeros((3, 256, 256))
            dummy_mask = torch.zeros((256, 256), dtype=torch.long)
            logging.warning(f"Returning dummy data for failed sample: {slice_data['source_file']}")
            return {
                'image': dummy_img,
                'mask': dummy_mask
            }
    
    def _convert_img_slice(self, img_slice):
        """Convert image slice using ISG processing (matching isg_mask_pred.py) with overflow protection"""
        # Convert to float64 to prevent overflow
        img_slice = img_slice.astype(np.float64)
        
        lac_factor = 33.33
        img_slice = img_slice * lac_factor
        
        # Clip values to prevent overflow
        img_slice = np.clip(img_slice, 0, 1)
        
        # Convert to 0-255 range with safe handling
        min_val = np.min(img_slice)
        max_val = 1.0
        
        # Avoid division by zero
        if max_val > min_val:
            img_slice = (img_slice - min_val) / (max_val - min_val) * 255
        else:
            img_slice = np.zeros_like(img_slice)
        
        # Safe conversion to uint8
        img_slice = np.clip(img_slice, 0, 255)
        return img_slice.astype(np.uint8)
    
    def _convert_mask_slice(self, mask_slice):
        """Convert mask slice to binary format"""
        return (mask_slice > 0).astype(np.uint8)

def train_net(net,
              device,
              epochs: int = 5,
              batch_size: int = 1,
              learning_rate: float = 1e-5,
              val_percent: float = 0.1,
              save_checkpoint: bool = True,
              img_scale: float = 0.5,
              amp: bool = False,
              custom_img_dir=None,
              custom_mask_dir=None,
              use_wandb: bool = False,
              allowed_ids=None,
              sample_ratio: float = 1.0,
              organelle: str = "isg",
              external_logger=None):  # Add external_logger parameter
    
    # Override paths if custom directories provided
    global dir_img, dir_mask
    if custom_img_dir:
        dir_img = Path(custom_img_dir)
    if custom_mask_dir:
        dir_mask = Path(custom_mask_dir)
    
    # Create organelle-specific checkpoint directory
    organelle_checkpoint_dir = Path(dir_checkpoint) / f'sxt_{organelle}_checkpoints'
    organelle_checkpoint_dir.mkdir(parents=True, exist_ok=True)
    
    # Setup training logs directory - organelle specific
    training_logs_dir = organelle_checkpoint_dir / 'training_logs'
    training_logs_dir.mkdir(parents=True, exist_ok=True)
    
    # Setup detailed logging
    log_file = training_logs_dir / f'training_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logging.getLogger().addHandler(file_handler)
    
    logging.info(f"Training log file: {log_file}")
    
    # Training metadata
    training_metadata = {
        'start_time': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'device': str(device),
        'epochs': epochs,
        'batch_size': batch_size,
        'learning_rate': learning_rate,
        'val_percent': val_percent,
        'img_scale': img_scale,
        'amp': amp,
        'img_dir': str(dir_img),
        'mask_dir': str(dir_mask),
        'checkpoint_dir': str(dir_checkpoint)
    }
    
    # Save training metadata
    metadata_file = training_logs_dir / f'metadata_{datetime.now().strftime("%Y%m%d_%H%M%S")}.json'
    with open(metadata_file, 'w') as f:
        json.dump(training_metadata, f, indent=2)
    
    logging.info(f"Training metadata saved to: {metadata_file}")

    # 1. Create dataset with ID filtering and sampling
    try:
        dataset = ISG3DDataset(dir_img, dir_mask, img_scale, allowed_ids=allowed_ids, sample_ratio=sample_ratio)
        logging.info(f"Using ISG3DDataset with {len(dataset)} slices (sample ratio: {sample_ratio})")
        
        # Log which IDs were used
        if allowed_ids:
            used_ids = set()
            for slice_info in dataset.slice_info:
                used_ids.add(slice_info['identifier'])
            logging.info(f"Dataset uses {len(used_ids)} IDs: {sorted(used_ids)}")
        
        # Log dataset statistics
        dataset_stats = {
            'total_slices': len(dataset),
            'total_files': len(dataset.file_slice_count),
            'files_info': {k: int(v) for k, v in dataset.file_slice_count.items()}  # Convert to Python int
        }
        
        dataset_stats_file = training_logs_dir / f'dataset_stats_{datetime.now().strftime("%Y%m%d_%H%M%S")}.json'
        with open(dataset_stats_file, 'w') as f:
            json.dump(dataset_stats, f, indent=2)
        
        logging.info(f"Dataset statistics saved to: {dataset_stats_file}")
        
    except Exception as e:
        logging.error(f"Failed to create ISG3DDataset: {e}")
        logging.error("ISG3DDataset creation failed. Please check your data paths and file structure.")
        raise RuntimeError(f"Dataset creation failed: {e}")

    # 2. Split into train / validation partitions
    n_val = int(len(dataset) * val_percent)
    n_train = len(dataset) - n_val
    train_set, val_set = random_split(dataset, [n_train, n_val], generator=torch.Generator().manual_seed(0))

    # Analyze split by data ID
    train_ids = set()
    val_ids = set()
    train_slice_count = {}
    val_slice_count = {}
    
    # Count slices per ID in training set
    for idx in train_set.indices:
        slice_info = dataset.slice_info[idx]
        data_id = slice_info['identifier']
        train_ids.add(data_id)
        train_slice_count[data_id] = train_slice_count.get(data_id, 0) + 1
    
    # Count slices per ID in validation set
    for idx in val_set.indices:
        slice_info = dataset.slice_info[idx]
        data_id = slice_info['identifier']
        val_ids.add(data_id)
        val_slice_count[data_id] = val_slice_count.get(data_id, 0) + 1
    
    # Prepare split information
    split_info = {
        'split_method': 'random_slice_split',
        'total_samples': len(dataset),
        'train_samples': n_train,
        'val_samples': n_val,
        'train_percentage': (n_train / len(dataset)) * 100,
        'val_percentage': (n_val / len(dataset)) * 100,
        'train_ids': sorted(list(train_ids)),
        'val_ids': sorted(list(val_ids)),
        'overlapping_ids': sorted(list(train_ids & val_ids)),
        'train_slice_count': dict(sorted(train_slice_count.items())),
        'val_slice_count': dict(sorted(val_slice_count.items()))
    }
    
    # Print and log split information
    split_summary = (f"Data Split Summary:\n"
                    f"  Split Method: Random slice-level split\n"
                    f"  Total 2D Slices: {len(dataset)}\n"
                    f"  Train Slices: {n_train} ({(n_train / len(dataset)) * 100:.1f}%)\n"
                    f"  Val Slices: {n_val} ({(n_val / len(dataset)) * 100:.1f}%)\n"
                    f"  Train IDs ({len(train_ids)}): {sorted(list(train_ids))}\n"
                    f"  Val IDs ({len(val_ids)}): {sorted(list(val_ids))}\n"
                    f"  Overlapping IDs ({len(train_ids & val_ids)}): {sorted(list(train_ids & val_ids))}")
    
    print(f"[DATA SPLIT] {split_summary}")
    logging.info(f"Data Split: {split_summary}")
    
    if external_logger:
        external_logger.step(f"[DATA SPLIT] {split_summary}")
    
    # Log detailed slice counts per ID
    train_detail = "Train slice counts: " + ", ".join([f"{id}:{count}" for id, count in sorted(train_slice_count.items())])
    val_detail = "Val slice counts: " + ", ".join([f"{id}:{count}" for id, count in sorted(val_slice_count.items())])
    
    print(f"[TRAIN DETAILS] {train_detail}")
    print(f"[VAL DETAILS] {val_detail}")
    logging.info(train_detail)
    logging.info(val_detail)
    
    if external_logger:
        external_logger.step(f"[TRAIN DETAILS] {train_detail}")
        external_logger.step(f"[VAL DETAILS] {val_detail}")

    # Log dataset split
    split_file = training_logs_dir / f'dataset_split_{datetime.now().strftime("%Y%m%d_%H%M%S")}.json'
    with open(split_file, 'w') as f:
        json.dump(split_info, f, indent=2)
    
    logging.info(f"Dataset split info: Train: {n_train}, Val: {n_val}")
    logging.info(f"Dataset split saved to: {split_file}")
    
    if external_logger:
        external_logger.step(f"Dataset split details saved to: {split_file}")

    # 3. Create data loaders - optimize for speed
    loader_args = dict(batch_size=batch_size, num_workers=4, pin_memory=True, prefetch_factor=2)
    train_loader = DataLoader(train_set, shuffle=True, **loader_args)
    val_loader = DataLoader(val_set, shuffle=False, drop_last=True, **loader_args)

    # Initialize wandb with better error handling
    experiment = None
    if use_wandb:
        try:
            # Disable wandb completely to avoid network issues
            os.environ["WANDB_MODE"] = "disabled"
            logging.info("W&B disabled to avoid network issues")
        except Exception as e:
            logging.warning(f"W&B setup failed: {e}")
    else:
        logging.info("W&B logging disabled by user")

    # Training progress tracking
    training_progress = {
        'epochs': [],
        'train_losses': [],
        'val_scores': [],
        'learning_rates': [],
        'epoch_times': []
    }

    logging.info(f'''Starting ISG training:
        Epochs:          {epochs}
        Batch size:      {batch_size}
        Learning rate:   {learning_rate}
        Training size:   {n_train}
        Validation size: {n_val}
        Sample ratio:    {sample_ratio}
        Checkpoints:     {save_checkpoint}
        Device:          {device.type}
        Images scaling:  {img_scale}
        Mixed Precision: {amp}
        Log directory:   {training_logs_dir}
        Data loading:    {loader_args['num_workers']} workers with prefetch
    ''')

    # 4. Set up the optimizer, the loss, the learning rate scheduler and the loss scaling for AMP
    optimizer = optim.RMSprop(net.parameters(), lr=learning_rate, weight_decay=1e-8, momentum=0.9)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'max', patience=2)  # goal: maximize Dice score
    grad_scaler = torch.cuda.amp.GradScaler(enabled=amp)
    criterion = nn.CrossEntropyLoss()
    global_step = 0

    # 5. Begin training
    best_val_score = 0.0
   
    for epoch in range(1, epochs + 1):
        epoch_start_time = time.time()
        net.train()
        epoch_loss = 0
        batch_losses = []

        logging.info(f"Starting epoch {epoch}/{epochs}")

        total_batches = len(train_loader)
        mid_batch = total_batches // 2 

        with tqdm(total=n_train, desc=f'Epoch {epoch}/{epochs}', unit='img') as pbar:
            for batch_idx, batch in enumerate(train_loader):
                images = batch['image']
                true_masks = batch['mask']
                assert images.shape[1] == net.n_channels, (
                    f'Network has been defined with {net.n_channels} input channels, '
                    f'but loaded images have {images.shape[1]} channels. Please check that '
                    'the images are loaded correctly.'
                )
                images = images.to(device=device, dtype=torch.float32)
                true_masks = true_masks.to(device=device, dtype=torch.long)
                with torch.cuda.amp.autocast(enabled=amp):
                    masks_pred = net(images)
                    loss = criterion(masks_pred, true_masks) + dice_loss(
                        F.softmax(masks_pred, dim=1).float(),
                        F.one_hot(true_masks, net.n_classes).permute(0, 3, 1, 2).float(),
                        multiclass=True
                    )
                optimizer.zero_grad(set_to_none=True)
                grad_scaler.scale(loss).backward()
                grad_scaler.step(optimizer)
                grad_scaler.update()
                pbar.update(images.shape[0])
                global_step += 1
                epoch_loss += loss.item()
                batch_losses.append(loss.item())

                pbar.set_postfix(**{'loss (batch)': loss.item()})

                # ---- Only validate at mid-epoch ----
                if batch_idx == mid_batch:
                    val_score = evaluate(net, val_loader, device)
                    scheduler.step(val_score)
                    logging.info(f'Mid-epoch validation Dice score: {val_score:.4f}')
                    val_score_float = float(val_score)
                    if val_score_float > best_val_score:
                        best_val_score = val_score_float
                        best_model_path = organelle_checkpoint_dir / f'best_{organelle}_unet.pth'
                        torch.save(net.state_dict(), str(best_model_path))
                        print(f"[BEST MODEL] New best model saved with Dice score: {val_score_float:.4f}")
                        logging.info(f'New best model saved with Dice score: {val_score_float:.4f}')
                        if external_logger:
                            external_logger.step(f"[BEST MODEL] New best model saved: {best_model_path} (Dice: {val_score_float:.4f})")

            # ---- Epoch-end validation ----
            val_score = evaluate(net, val_loader, device)
            scheduler.step(val_score)
            logging.info(f'Epoch-end validation Dice score: {val_score:.4f}')
            val_score_float = float(val_score)
            if val_score_float > best_val_score:
                best_val_score = val_score_float
                best_model_path = organelle_checkpoint_dir / f'best_{organelle}_unet.pth'
                torch.save(net.state_dict(), str(best_model_path))
                print(f"[BEST MODEL] Final best model saved with Dice score: {val_score_float:.4f}")
                logging.info(f'New best model saved with Dice score: {val_score_float:.4f}')
                if external_logger:
                    external_logger.step(f"[BEST MODEL] Final best model saved: {best_model_path} (Dice: {val_score_float:.4f})")

        # Calculate epoch statistics
        epoch_time = time.time() - epoch_start_time
        avg_epoch_loss = epoch_loss / len(train_loader)
        current_lr = optimizer.param_groups[0]['lr']
        
        # Final validation for this epoch
        final_val_score = evaluate(net, val_loader, device)
        final_val_score_float = float(final_val_score)  # Convert tensor to float
        
        # Calculate additional metrics (AUC and Accuracy)
        net.eval()
        total_correct = 0
        total_pixels = 0
        all_predictions = []
        all_targets = []
        
        with torch.no_grad():
            for batch in val_loader:
                images = batch['image'].to(device=device, dtype=torch.float32)
                true_masks = batch['mask'].to(device=device, dtype=torch.long)
                
                masks_pred = net(images)
                probs = F.softmax(masks_pred, dim=1)
                predictions = torch.argmax(masks_pred, dim=1)
                
                # Calculate accuracy
                total_correct += (predictions == true_masks).sum().item()
                total_pixels += true_masks.numel()
                
                # Store for AUC calculation (for binary/multi-class)
                if net.n_classes == 2:  # Binary classification
                    all_predictions.extend(probs[:, 1].cpu().numpy().flatten())  # Positive class probability
                    all_targets.extend(true_masks.cpu().numpy().flatten())
                else:  # Multi-class - use macro average
                    # For multi-class, we'll calculate per-class and then average
                    for c in range(net.n_classes):
                        class_probs = probs[:, c].cpu().numpy().flatten()
                        class_targets = (true_masks == c).cpu().numpy().flatten().astype(int)
                        all_predictions.append(class_probs)
                        all_targets.append(class_targets)
        
        # Calculate accuracy
        accuracy = total_correct / total_pixels if total_pixels > 0 else 0.0
        
        # Calculate AUC
        auc = 0.0
        try:
            from sklearn.metrics import roc_auc_score
            import numpy as np
            
            if net.n_classes == 2:  # Binary classification
                if len(np.unique(all_targets)) > 1:  # Check if both classes exist
                    auc = roc_auc_score(all_targets, all_predictions)
            else:  # Multi-class classification
                auc_scores = []
                for c in range(net.n_classes):
                    if len(np.unique(all_targets[c])) > 1:  # Check if both classes exist
                        auc_c = roc_auc_score(all_targets[c], all_predictions[c])
                        auc_scores.append(auc_c)
                auc = np.mean(auc_scores) if auc_scores else 0.0
                
        except ImportError:
            logging.warning("scikit-learn not available for AUC calculation")
        except Exception as e:
            logging.warning(f"AUC calculation failed: {e}")

        # Print and log epoch results with all metrics
        epoch_summary = (f"Epoch {epoch}/{epochs} completed in {epoch_time:.2f}s - "
                        f"Avg Loss: {avg_epoch_loss:.4f}, Val Dice: {final_val_score_float:.4f}, "
                        f"Accuracy: {accuracy:.4f}, AUC: {auc:.4f}, "
                        f"LR: {current_lr:.6f}, Best Val: {best_val_score:.4f}")
        
        print(f"[EPOCH SUMMARY] {epoch_summary}")
        logging.info(epoch_summary)
        
        # Also log to external logger if provided
        if external_logger:
            external_logger.step(f"[EPOCH SUMMARY] {epoch_summary}")

        # Log epoch results
        epoch_info = {
            'epoch': epoch,
            'avg_loss': avg_epoch_loss,
            'val_score': final_val_score,
            'accuracy': accuracy,
            'auc': auc,
            'learning_rate': current_lr,
            'epoch_time': epoch_time,
            'best_val_score': best_val_score
        }
        
        # Convert tensors to float before storing in training_progress
        training_progress['epochs'].append(epoch)
        training_progress['train_losses'].append(float(avg_epoch_loss))  # Convert to float
        training_progress['val_scores'].append(float(final_val_score))   # Convert to float
        training_progress['learning_rates'].append(float(current_lr))    # Convert to float
        training_progress['epoch_times'].append(float(epoch_time))       # Convert to float
        
        # Print additional detailed metrics
        detailed_metrics = f"[EPOCH METRICS] Epoch {epoch}: Time={epoch_time:.2f}s, Train_Loss={avg_epoch_loss:.4f}, Val_Dice={final_val_score_float:.4f}, Accuracy={accuracy:.4f}, AUC={auc:.4f}, LR={current_lr:.6f}"
        print(detailed_metrics)
        
        logging.info(f'Epoch {epoch} completed in {epoch_time:.2f}s - '
                    f'Avg Loss: {avg_epoch_loss:.4f}, Val Dice: {final_val_score:.4f}, '
                    f'Accuracy: {accuracy:.4f}, AUC: {auc:.4f}, '
                    f'LR: {current_lr:.6f}')
        
        # Also log to external logger if provided
        if external_logger:
            external_logger.step(detailed_metrics)

        # Save checkpoint every 2 epochs (or on the last epoch)
        if save_checkpoint and (epoch % 2 == 0 or epoch == epochs):
            checkpoint_name = f'{organelle}_unet_epoch{epoch}.pth'
            checkpoint_path = organelle_checkpoint_dir / checkpoint_name
            torch.save({
                'epoch': epoch,
                'model_state_dict': net.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'train_loss': avg_epoch_loss,
                'val_loss': final_val_score_float,
                'accuracy': accuracy,
                'auc': auc,
                'learning_rate': current_lr,
                'best_val_score': best_val_score,
                'training_config': {
                    'epochs': epochs,
                    'batch_size': batch_size,
                    'learning_rate': learning_rate,
                    'organelle': organelle
                }
            }, str(checkpoint_path))
            print(f"[CHECKPOINT] Saved: {checkpoint_name}")
            logging.info(f'Checkpoint saved: {checkpoint_name}')
            if external_logger:
                external_logger.step(f"[CHECKPOINT] Saved: {checkpoint_path}")

    # Save final training progress
    progress_file = training_logs_dir / f'training_progress_{datetime.now().strftime("%Y%m%d_%H%M%S")}.json'
    with open(progress_file, 'w') as f:
        json.dump(training_progress, f, indent=2)
    
    # Training summary - ensure all values are JSON serializable
    training_summary = {
        'organelle': organelle,
        'total_epochs': int(epochs),
        'final_val_score': float(training_progress['val_scores'][-1]) if training_progress['val_scores'] else 0.0,
        'best_val_score': float(best_val_score),
        'final_loss': float(training_progress['train_losses'][-1]) if training_progress['train_losses'] else 0.0,
        'total_training_time': float(sum(training_progress['epoch_times'])) if training_progress['epoch_times'] else 0.0,
        'average_epoch_time': float(sum(training_progress['epoch_times']) / len(training_progress['epoch_times'])) if training_progress['epoch_times'] else 0.0,
        'end_time': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'checkpoint_dir': str(organelle_checkpoint_dir)
    }
    
    summary_file = training_logs_dir / f'{organelle}_training_summary_{datetime.now().strftime("%Y%m%d_%H%M%S")}.json'
    with open(summary_file, 'w') as f:
        json.dump(training_summary, f, indent=2)
    
    # Save final model with organelle-specific name
    final_model_path = organelle_checkpoint_dir / f'final_{organelle}_unet.pth'
    torch.save({
        'model_state_dict': net.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'training_summary': training_summary,
        'training_progress': training_progress,
        'final_epoch': epochs
    }, str(final_model_path))
    
    logging.info("="*60)
    logging.info(f"TRAINING COMPLETED - {organelle.upper()}")
    logging.info("="*60)
    logging.info(f"Best validation Dice score: {best_val_score:.4f}")
    logging.info(f"Final validation Dice score: {training_progress['val_scores'][-1]:.4f}")
    logging.info(f"Total training time: {training_summary['total_training_time']:.2f}s")
    logging.info(f"Checkpoints saved to: {organelle_checkpoint_dir}")
    logging.info(f"Final model saved to: {final_model_path}")
    logging.info(f"Training logs saved to: {training_logs_dir}")
    logging.info(f"Training summary saved to: {summary_file}")
    logging.info("="*60)
    
    if external_logger:
        external_logger.step(f"[TRAINING COMPLETE] {organelle.upper()} model training finished")
        external_logger.step(f"[FINAL RESULTS] Best Dice: {best_val_score:.4f}, Final Dice: {training_progress['val_scores'][-1]:.4f}")
        external_logger.step(f"[FINAL RESULTS] Total time: {training_summary['total_training_time']:.2f}s")
        external_logger.step(f"[SAVED] Final model: {final_model_path}")
        external_logger.step(f"[SAVED] Checkpoints: {organelle_checkpoint_dir}")
        external_logger.step(f"[SAVED] Training summary: {summary_file}")

def get_args():
    parser = argparse.ArgumentParser(description='Train the UNet on images and target masks')
    parser.add_argument('--epochs', '-e', metavar='E', type=int, default=30, help='Number of epochs')
    parser.add_argument('--batch-size', '-b', dest='batch_size', metavar='B', type=int, default=4, help='Batch size')
    parser.add_argument('--learning-rate', '-l', metavar='LR', type=float, default=1e-5,
                        help='Learning rate', dest='lr')
    parser.add_argument('--load', '-f', type=str, default=False, help='Load model from a .pth file')
    parser.add_argument('--scale', '-s', type=float, default=0.5, help='Downscaling factor of the images')
    parser.add_argument('--validation', '-v', dest='val', type=float, default=10.0,
                        help='Percent of the data that is used as validation (0-100)')
    parser.add_argument('--amp', action='store_true', default=False, help='Use mixed precision')
    parser.add_argument('--bilinear', action='store_true', default=False, help='Use bilinear upsampling')
    parser.add_argument('--classes', '-c', type=int, default=2, help='Number of classes')

    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()

    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logging.info(f'Using device {device}')

    # Change here to adapt to your data
    # n_channels=3 for RGB images, n_channels=1 for gray images  
    # n_classes is the number of probabilities you want to get per pixel
    net = UNet(n_channels= 3, n_classes=args.classes, bilinear=args.bilinear)

    logging.info(f'Network:\n'
                 f'\t{net.n_channels} input channels\n'
                 f'\t{net.n_classes} output channels (classes)\n'
                 f'\t{"Bilinear" if net.bilinear else "Transposed conv"} upscaling')

    if args.load:
        net.load_state_dict(torch.load(args.load, map_location=device))
        logging.info(f'Model loaded from {args.load}')

    net.to(device=device)
    try:
        train_net(net=net,
                  epochs=args.epochs,
                  batch_size=args.batch_size,
                  learning_rate=args.lr,
                  device=device,
                  img_scale=args.scale,
                  val_percent=args.val / 100,
                  amp=args.amp)
    except KeyboardInterrupt:
        torch.save(net.state_dict(), 'INTERRUPTED.pth')
        logging.info('Saved interrupt')
        raise

# %%
