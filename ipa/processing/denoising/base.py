"""
Base class for denoising methods.

Provides a unified interface for different denoising algorithms
including Noise2Void (N2V) and Noise2Noise (N2N).
"""

from abc import ABC, abstractmethod
import numpy as np
import torch
from typing import Optional, Union


class BaseDenoiser(ABC):
    """Abstract base class for all denoising methods."""
    
    def __init__(self, n_channels: int = 1, n_filters: int = 64, device: Optional[str] = None):
        """
        Initialize the denoiser.
        
        Args:
            n_channels: Number of input channels (default: 1 for grayscale)
            n_filters: Number of convolutional filters (default: 64)
            device: Device to run on ('cuda' or 'cpu'). Auto-detected if None.
        """
        self.n_channels = n_channels
        self.n_filters = n_filters
        
        # Auto-detect device
        if device is None:
            self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        else:
            self.device = torch.device(device)
            
        self.model = None
        self.is_trained = False
        
    @abstractmethod
    def train(self, *args, **kwargs):
        """
        Train the denoising model.
        
        Must be implemented by subclasses.
        """
        pass
    
    @abstractmethod
    def predict(self, data: np.ndarray, **kwargs) -> np.ndarray:
        """
        Predict denoised output.
        
        Args:
            data: Input noisy data (numpy array)
            
        Returns:
            Denoised data (numpy array)
        """
        pass
    
    def save_model(self, path: str):
        """
        Save trained model to file.
        
        Args:
            path: Path to save the model
        """
        if not self.is_trained:
            raise RuntimeError("Model has not been trained yet!")
        
        torch.save({
            'model_state_dict': self.model.state_dict(),
            'n_channels': self.n_channels,
            'n_filters': self.n_filters,
        }, path)
        print(f"Model saved to: {path}")
    
    def load_model(self, path: str):
        """
        Load trained model from file.
        
        Args:
            path: Path to the saved model
        """
        checkpoint = torch.load(path, map_location=self.device)
        
        # Verify compatibility
        if checkpoint['n_channels'] != self.n_channels:
            raise ValueError(f"Channel mismatch: expected {self.n_channels}, got {checkpoint['n_channels']}")
        if checkpoint['n_filters'] != self.n_filters:
            raise ValueError(f"Filter mismatch: expected {self.n_filters}, got {checkpoint['n_filters']}")
        
        self.model.load_state_dict(checkpoint['model_state_dict'])
        self.is_trained = True
        print(f"Model loaded from: {path}")
    
    def _normalize_data(self, data: np.ndarray) -> np.ndarray:
        """
        Normalize data to [0, 1] range.
        
        Args:
            data: Input data
            
        Returns:
            Normalized data
        """
        data = data.astype(np.float32)
        if data.max() > 1.0:
            data = data / 255.0
        return data
    
    def _ensure_4d(self, data: np.ndarray) -> np.ndarray:
        """
        Ensure data is 4D: (Samples, Height, Width, Channels).
        
        Args:
            data: Input data (2D, 3D, or 4D)
            
        Returns:
            4D data array
        """
        if data.ndim == 2:
            # (H, W) -> (1, H, W, 1)
            return data[np.newaxis, ..., np.newaxis]
        elif data.ndim == 3:
            # (D, H, W) -> (D, H, W, 1)
            return data[..., np.newaxis]
        elif data.ndim == 4:
            # Already 4D
            return data
        else:
            raise ValueError(f"Unsupported data dimensions: {data.ndim}D")
