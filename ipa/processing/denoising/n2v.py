"""
Noise2Void (N2V) denoising implementation.

N2V is a self-supervised denoising method that trains on single noisy images
without requiring clean ground truth or noise pairs.
"""

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Dataset
from typing import Optional

from .base import BaseDenoiser
from .noise2void.model import UNet
from .noise2void.training import train_loop


class N2VDataset(Dataset):
    """Dataset for N2V training with blind-spot masking."""
    
    def __init__(self, data: np.ndarray, mask_ratio: float = 0.195, 
                 patch_size: tuple = (64, 64)):
        """
        Args:
            data: Noisy input data, shape (D, H, W) or (D, H, W, C)
            mask_ratio: Ratio of pixels to mask (default: 0.195 from original N2V)
            patch_size: Size of patches to extract (height, width)
        """
        self.data = self._prepare_data(data)
        self.mask_ratio = mask_ratio
        self.patch_size = patch_size
        
    def _prepare_data(self, data: np.ndarray) -> np.ndarray:
        """Prepare and normalize data."""
        data = data.astype(np.float32)
        if data.max() > 1.0:
            data = data / 255.0
        
        # Ensure 4D: (Samples, H, W, C)
        if data.ndim == 3:
            data = data[..., np.newaxis]
        
        return data
    
    def __len__(self) -> int:
        return len(self.data)
    
    def __getitem__(self, idx: int):
        """
        Get a sample with blind-spot masking applied.
        
        Returns:
            tuple: (masked_input, target, mask)
        """
        img = self.data[idx]
        
        # Extract random patch
        h, w = img.shape[:2]
        ph, pw = self.patch_size
        
        if h > ph:
            y = np.random.randint(0, h - ph)
        else:
            y = 0
            ph = h
            
        if w > pw:
            x = np.random.randint(0, w - pw)
        else:
            x = 0
            pw = w
        
        patch = img[y:y+ph, x:x+pw, :].copy()
        
        # Apply blind-spot masking
        masked_patch, mask = self._apply_masking(patch)
        
        # Convert to tensor: (C, H, W)
        masked_tensor = torch.from_numpy(masked_patch.transpose(2, 0, 1)).float()
        target_tensor = torch.from_numpy(patch.transpose(2, 0, 1)).float()
        mask_tensor = torch.from_numpy(mask.transpose(2, 0, 1)).float()
        
        return masked_tensor, target_tensor, mask_tensor
    
    def _apply_masking(self, patch: np.ndarray) -> tuple:
        """
        Apply blind-spot masking to patch.
        
        Args:
            patch: Input patch (H, W, C)
            
        Returns:
            tuple: (masked_patch, mask)
        """
        masked = patch.copy()
        mask = np.ones_like(patch)
        
        h, w, c = patch.shape
        n_pixels = h * w
        n_mask = int(n_pixels * self.mask_ratio)
        
        # Randomly select pixels to mask
        indices = np.random.choice(n_pixels, n_mask, replace=False)
        
        for idx in indices:
            y, x = divmod(idx, w)
            # Replace with neighbor value (simple strategy)
            ny, nx = y, max(0, x - 1)  # Use left neighbor
            masked[y, x, :] = patch[ny, nx, :]
            mask[y, x, :] = 0.0
        
        return masked, mask


class N2V(BaseDenoiser):
    """
    Noise2Void denoiser.
    
    Example:
        >>> from ipa.processing.denoising import N2V
        >>> 
        >>> # Create denoiser
        >>> n2v = N2V(n_filters=64)
        >>> 
        >>> # Train on noisy data
        >>> n2v.train(noisy_data, epochs=50, batch_size=4)
        >>> 
        >>> # Predict
        >>> denoised = n2v.predict(noisy_data)
        >>> 
        >>> # Save/Load model
        >>> n2v.save_model('n2v_model.pth')
        >>> n2v.load_model('n2v_model.pth')
    """
    
    def __init__(self, n_channels: int = 1, n_filters: int = 64, 
                 device: Optional[str] = None):
        """
        Initialize N2V denoiser.
        
        Args:
            n_channels: Number of input channels (default: 1)
            n_filters: Number of convolutional filters (default: 64)
            device: Device to run on ('cuda' or 'cpu')
        """
        super().__init__(n_channels=n_channels, n_filters=n_filters, device=device)
        
        # Initialize UNet model
        self.model = UNet(
            nch_in=n_channels,
            nch_out=n_channels,
            nch_ker=n_filters,
            norm='bnorm'
        )
        self.model.to(self.device)
    
    def train(self, train_data: np.ndarray, val_data: Optional[np.ndarray] = None,
              epochs: int = 50, batch_size: int = 4, lr: float = 1e-3,
              mask_ratio: float = 0.195, patch_size: tuple = (64, 64)):
        """
        Train the N2V model.
        
        Args:
            train_data: Training data, shape (D, H, W) or (D, H, W, C)
            val_data: Optional validation data
            epochs: Number of training epochs (default: 50)
            batch_size: Batch size (default: 4)
            lr: Learning rate (default: 1e-3)
            mask_ratio: Ratio of pixels to mask (default: 0.195)
            patch_size: Patch size for training (default: (64, 64))
        """
        print(f"Training N2V model on {self.device}")
        print(f"Training data shape: {train_data.shape}")
        print(f"Epochs: {epochs}, Batch size: {batch_size}, LR: {lr}")
        
        # Create datasets
        train_dataset = N2VDataset(train_data, mask_ratio=mask_ratio, patch_size=patch_size)
        train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, num_workers=0)
        
        val_loader = None
        if val_data is not None:
            val_dataset = N2VDataset(val_data, mask_ratio=mask_ratio, patch_size=patch_size)
            val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False, num_workers=0)
        
        # Setup optimizer and loss
        optimizer = torch.optim.Adam(self.model.parameters(), lr=lr)
        criterion = nn.MSELoss()
        
        # Train
        for epoch in range(epochs):
            self.model.train()
            total_loss = 0.0
            n_batches = 0
            
            for masked_input, target, mask in train_loader:
                masked_input = masked_input.to(self.device)
                target = target.to(self.device)
                mask = mask.to(self.device)
                
                # Forward pass
                output = self.model(masked_input)
                
                # Compute loss only on masked pixels
                loss = criterion(output * (1 - mask), target * (1 - mask))
                
                # Backward pass
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                
                total_loss += loss.item()
                n_batches += 1
            
            avg_loss = total_loss / n_batches
            
            # Validation
            val_loss = None
            if val_loader is not None:
                self.model.eval()
                val_total_loss = 0.0
                val_n_batches = 0
                
                with torch.no_grad():
                    for masked_input, target, mask in val_loader:
                        masked_input = masked_input.to(self.device)
                        target = target.to(self.device)
                        mask = mask.to(self.device)
                        
                        output = self.model(masked_input)
                        loss = criterion(output * (1 - mask), target * (1 - mask))
                        
                        val_total_loss += loss.item()
                        val_n_batches += 1
                
                val_loss = val_total_loss / val_n_batches
            
            # Print progress
            if (epoch + 1) % 10 == 0 or epoch == 0:
                if val_loss is not None:
                    print(f"Epoch [{epoch+1}/{epochs}] Loss: {avg_loss:.6f}, Val Loss: {val_loss:.6f}")
                else:
                    print(f"Epoch [{epoch+1}/{epochs}] Loss: {avg_loss:.6f}")
        
        self.is_trained = True
        print("Training completed!")
    
    def predict(self, data: np.ndarray, batch_size: int = 4) -> np.ndarray:
        """
        Predict denoised output.
        
        Args:
            data: Noisy input data, shape (D, H, W) or (D, H, W, C)
            batch_size: Batch size for prediction (default: 4)
            
        Returns:
            Denoised data, same shape as input
        """
        if not self.is_trained:
            raise RuntimeError("Model has not been trained! Call train() first.")
        
        print(f"Predicting on {self.device}")
        
        # Prepare data
        data_normalized = self._normalize_data(data)
        data_4d = self._ensure_4d(data_normalized)
        
        # Create dataset (no masking for prediction)
        dataset = SimpleDataset(data_4d)
        loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=0)
        
        # Predict
        self.model.eval()
        predictions = []
        
        with torch.no_grad():
            for batch in loader:
                batch = batch.to(self.device)
                output = self.model(batch)
                predictions.append(output.cpu().numpy())
        
        # Concatenate and reshape
        result = np.concatenate(predictions, axis=0)  # (D, C, H, W)
        
        # Transpose to (D, H, W, C)
        result = result.transpose(0, 2, 3, 1)
        
        # Remove channel dimension if input was 3D
        if data.ndim == 3:
            result = result[..., 0]  # (D, H, W)
        
        print(f"Prediction completed! Output shape: {result.shape}")
        return result


class SimpleDataset(Dataset):
    """Simple dataset for prediction (no masking)."""
    
    def __init__(self, data: np.ndarray):
        self.data = data
    
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        img = self.data[idx]
        # Convert to tensor: (C, H, W)
        return torch.from_numpy(img.transpose(2, 0, 1)).float()
