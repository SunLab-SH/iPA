"""
Base class for segmentation methods.

Provides a unified interface for different segmentation tasks
across multiple imaging modalities (SXT, SIM, WFM, Cryo-ET).
"""

from abc import ABC, abstractmethod
import numpy as np
import torch
from typing import Optional, Dict, Any


class BaseSegmenter(ABC):
    """Abstract base class for all segmentation methods."""
    
    def __init__(self, modality: str, task: str, device: Optional[str] = None):
        """
        Initialize the segmenter.
        
        Args:
            modality: Imaging modality ('sxt', 'sim', 'wfm', 'et')
            task: Segmentation task ('cell', 'mito', 'er', 'nucleus', 'filament', etc.)
            device: Device to run on ('cuda' or 'cpu'). Auto-detected if None.
        """
        self.modality = modality.lower()
        self.task = task.lower()
        
        # Auto-detect device
        if device is None:
            self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        else:
            self.device = torch.device(device)
            
        self.model = None
        self.is_loaded = False
        
    @abstractmethod
    def predict(self, data: np.ndarray, **kwargs) -> np.ndarray:
        """
        Predict segmentation mask.
        
        Args:
            data: Input image/volume data
            
        Returns:
            Segmentation mask (numpy array)
        """
        pass
    
    def load_model(self, path: str = None, **kwargs):
        """
        Load trained model from file.
        
        Args:
            path: Path to the saved model. If None, uses default path if available.
            **kwargs: Additional parameters for model loading
        """
        if not hasattr(self, '_load_model_impl'):
            raise NotImplementedError("Model loading not implemented for this segmenter")
        
        self._load_model_impl(path, **kwargs)
        self.is_loaded = True
        if path:
            print(f"Model loaded from: {path}")
        else:
            print(f"Model loaded using default path.")
    
    def save_model(self, path: str):
        """
        Save trained model to file.
        
        Args:
            path: Path to save the model
        """
        if not self.is_loaded:
            raise RuntimeError("No model loaded to save!")
        
        if not hasattr(self, '_save_model_impl'):
            raise NotImplementedError("Model saving not implemented for this segmenter")
        
        self._save_model_impl(path)
        print(f"Model saved to: {path}")
    
    def train(self, train_data: Any, val_data: Optional[Any] = None, **kwargs):
        """
        Train the segmentation model (optional).
        
        Args:
            train_data: Training data
            val_data: Optional validation data
            **kwargs: Training parameters
        """
        if not hasattr(self, '_train_impl'):
            raise NotImplementedError("Training not implemented for this segmenter")
        
        self._train_impl(train_data, val_data, **kwargs)
        self.is_loaded = True
        print("Training completed!")
    
    def _normalize_data(self, data: np.ndarray) -> np.ndarray:
        """
        Normalize data to [0, 1] range.
        
        Args:
            data: Input data
            
        Returns:
            Normalized data
        """
        data = data.astype(np.float32)
        min_val = data.min()
        max_val = data.max()
        
        if max_val > min_val:
            data = (data - min_val) / (max_val - min_val)
        
        return data
    
    def _ensure_channel_dim(self, data: np.ndarray, n_channels: int = 1) -> np.ndarray:
        """
        Ensure data has channel dimension.
        
        Args:
            data: Input data
            n_channels: Number of channels
            
        Returns:
            Data with channel dimension
        """
        if data.ndim == 2:
            # (H, W) -> (1, H, W) or (H, W, 1)
            return data[np.newaxis, ...] if n_channels == 1 else data[..., np.newaxis]
        elif data.ndim == 3:
            # Could be (D, H, W) or (H, W, C)
            if n_channels == 1:
                return data[np.newaxis, ...]  # (1, D, H, W)
            else:
                return data  # Assume already (H, W, C)
        else:
            return data
    
    def get_info(self) -> Dict[str, str]:
        """Get segmenter information."""
        return {
            'modality': self.modality,
            'task': self.task,
            'device': str(self.device),
            'is_loaded': str(self.is_loaded)
        }
