"""
Noise2Noise (N2N) denoising implementation.

N2N trains on pairs of noisy images without requiring clean ground truth.
"""

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Dataset
from typing import Optional, Tuple

from .base import BaseDenoiser
from .noise2void.model import UNet


class N2NDataset(Dataset):
    """Dataset for N2N training with noisy image pairs."""
    
    def __init__(self, noisy_data_1: np.ndarray, noisy_data_2: np.ndarray,
                 patch_size: tuple = (64, 64)):
        """
        Args:
            noisy_data_1: First set of noisy images, shape (D, H, W) or (D, H, W, C)
            noisy_data_2: Second set of noisy images (same content, different noise)
            patch_size: Size of patches to extract (height, width)
        """
        self.data_1 = self._prepare_data(noisy_data_1)
        self.data_2 = self._prepare_data(noisy_data_2)
        self.patch_size = patch_size
        
        assert len(self.data_1) == len(self.data_2), "Both datasets must have same length"
        
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
        return len(self.data_1)
    
    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Get a pair of noisy patches.
        
        Returns:
            tuple: (noisy_input, noisy_target)
        """
        img1 = self.data_1[idx]
        img2 = self.data_2[idx]
        
        # Extract random patch
        h, w = img1.shape[:2]
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
        
        patch1 = img1[y:y+ph, x:x+pw, :]
        patch2 = img2[y:y+ph, x:x+pw, :]
        
        # Convert to tensor: (C, H, W)
        tensor1 = torch.from_numpy(patch1.transpose(2, 0, 1)).float()
        tensor2 = torch.from_numpy(patch2.transpose(2, 0, 1)).float()
        
        return tensor1, tensor2


class N2N(BaseDenoiser):
    """
    Noise2Noise denoiser.
    
    Example:
        >>> from ipa.processing.denoising import N2N
        >>> 
        >>> # Create denoiser
        >>> n2n = N2N(n_filters=64)
        >>> 
        >>> # Train on pairs of noisy images
        >>> n2n.train(noisy_data_1, noisy_data_2, epochs=50, batch_size=4)
        >>> 
        >>> # Predict
        >>> denoised = n2n.predict(noisy_data)
        >>> 
        >>> # Save/Load model
        >>> n2n.save_model('n2n_model.pth')
        >>> n2n.load_model('n2n_model.pth')
    """
    
    def __init__(self, n_channels: int = 1, n_filters: int = 64,
                 device: Optional[str] = None):
        """
        Initialize N2N denoiser.
        
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
    
    def train(self, noisy_data_1: np.ndarray, noisy_data_2: np.ndarray,
              val_data_1: Optional[np.ndarray] = None,
              val_data_2: Optional[np.ndarray] = None,
              epochs: int = 50, batch_size: int = 4, lr: float = 1e-3,
              patch_size: tuple = (64, 64), loss_type: str = 'l1'):
        """
        Train the N2N model.
        
        Args:
            noisy_data_1: First set of noisy images
            noisy_data_2: Second set of noisy images (paired with data_1)
            val_data_1: Optional validation data (first set)
            val_data_2: Optional validation data (second set)
            epochs: Number of training epochs (default: 50)
            batch_size: Batch size (default: 4)
            lr: Learning rate (default: 1e-3)
            patch_size: Patch size for training (default: (64, 64))
            loss_type: Loss function type, 'l1' or 'mse' (default: 'l1')
        """
        print(f"Training N2N model on {self.device}")
        print(f"Training data shapes: {noisy_data_1.shape}, {noisy_data_2.shape}")
        print(f"Epochs: {epochs}, Batch size: {batch_size}, LR: {lr}")
        
        # Create datasets
        train_dataset = N2NDataset(noisy_data_1, noisy_data_2, patch_size=patch_size)
        train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, num_workers=0)
        
        val_loader = None
        if val_data_1 is not None and val_data_2 is not None:
            val_dataset = N2NDataset(val_data_1, val_data_2, patch_size=patch_size)
            val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False, num_workers=0)
        
        # Setup optimizer and loss
        optimizer = torch.optim.Adam(self.model.parameters(), lr=lr)
        
        if loss_type == 'l1':
            criterion = nn.L1Loss()
        elif loss_type == 'mse':
            criterion = nn.MSELoss()
        else:
            raise ValueError(f"Unknown loss type: {loss_type}. Use 'l1' or 'mse'.")
        
        # Train
        for epoch in range(epochs):
            self.model.train()
            total_loss = 0.0
            n_batches = 0
            
            for noisy_input, noisy_target in train_loader:
                noisy_input = noisy_input.to(self.device)
                noisy_target = noisy_target.to(self.device)
                
                # Forward pass
                output = self.model(noisy_input)
                
                # Compute loss
                loss = criterion(output, noisy_target)
                
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
                    for noisy_input, noisy_target in val_loader:
                        noisy_input = noisy_input.to(self.device)
                        noisy_target = noisy_target.to(self.device)
                        
                        output = self.model(noisy_input)
                        loss = criterion(output, noisy_target)
                        
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
        
        # Create dataset
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
    """Simple dataset for prediction."""
    
    def __init__(self, data: np.ndarray):
        self.data = data
    
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        img = self.data[idx]
        # Convert to tensor: (C, H, W)
        return torch.from_numpy(img.transpose(2, 0, 1)).float()
