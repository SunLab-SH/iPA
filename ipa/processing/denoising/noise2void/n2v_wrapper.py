import os
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Dataset
import copy

from .model import UNet
from .training import train_loop

class MemoryDataset(Dataset):
    """
    Dataset for N2V that works with in-memory numpy arrays.
    It expects 3D data (Depth, Height, Width) or 4D data (Depth, Height, Width, Channel).
    If 3D, it adds a channel dimension.
    
    IMPORTANT: For N2V training, input data should ALREADY be noisy images.
    N2V is self-supervised and learns from the noise itself.
    """
    def __init__(self, data, sgm=0, ratio=0.9, size_window=(5, 5), add_noise=False):
        """
        Args:
            data: Numpy array. Shape (D, H, W) or (D, H, W, C). Should be noisy data for N2V.
            sgm: Sigma for Gaussian noise. Only used if add_noise=True (for data augmentation).
            ratio: Ratio of pixels to keep (1 - masking_ratio).
            size_window: Window size for neighbor replacement.
            add_noise: Whether to add additional Gaussian noise. Default False for N2V.
        """
        self.data = data
        self.sgm = sgm
        self.ratio = ratio
        self.size_window = size_window
        self.add_noise = add_noise
        
        # Ensure 4D shape: (Samples, Height, Width, Channels)
        if self.data.ndim == 3:
            # Assume (D, H, W) -> (D, H, W, 1)
            self.data = self.data[..., np.newaxis]
        elif self.data.ndim == 2:
            # Assume (H, W) -> (1, H, W, 1)
            self.data = self.data[np.newaxis, ..., np.newaxis]
            
        # Normalize to 0-1 if likely in 0-255 range
        if self.data.max() > 1.0:
            self.data = self.data.astype(np.float32) / 255.0
        else:
            self.data = self.data.astype(np.float32)
            
        # Only add noise if explicitly requested (for augmentation)
        # For standard N2V, input data is already noisy, so we don't add more noise
        if self.add_noise and self.sgm > 0:
            self.noise = self.sgm / 255.0 * np.random.randn(*self.data.shape)
        else:
            self.noise = np.zeros_like(self.data)

    def __len__(self):
        return self.data.shape[0]

    def __getitem__(self, index):
        # Get single slice: (H, W, C)
        clean_slice = self.data[index]
        noise_slice = self.noise[index]
        
        # For N2V: input is already noisy, so label = input (no extra noise added)
        # If add_noise=True, we add small augmentation noise
        label = clean_slice + noise_slice
        input_img, mask = self.generate_mask(copy.deepcopy(label))
        
        # Transpose to PyTorch format: (C, H, W)
        input_tensor = torch.from_numpy(input_img.transpose((2, 0, 1)).astype(np.float32))
        label_tensor = torch.from_numpy(label.transpose((2, 0, 1)).astype(np.float32))
        mask_tensor = torch.from_numpy(mask.transpose((2, 0, 1)).astype(np.float32))
        
        return {'input': input_tensor, 'label': label_tensor, 'mask': mask_tensor}

    def generate_mask(self, input_img):
        """
        Generates N2V blind-spot mask.
        Input shape: (H, W, C)
        """
        ratio = self.ratio
        size_window = self.size_window
        size_data = input_img.shape
        num_sample = int(size_data[0] * size_data[1] * (1 - ratio))

        mask = np.ones(size_data, dtype=np.float32)
        output = input_img.copy()

        for ich in range(size_data[2]):
            idy_msk = np.random.randint(0, size_data[0], num_sample)
            idx_msk = np.random.randint(0, size_data[1], num_sample)

            idy_neigh = np.random.randint(-size_window[0] // 2 + size_window[0] % 2, 
                                          size_window[0] // 2 + size_window[0] % 2, num_sample)
            idx_neigh = np.random.randint(-size_window[1] // 2 + size_window[1] % 2, 
                                          size_window[1] // 2 + size_window[1] % 2, num_sample)

            idy_msk_neigh = idy_msk + idy_neigh
            idx_msk_neigh = idx_msk + idx_neigh

            # Boundary checking (wrapping)
            idy_msk_neigh = idy_msk_neigh + (idy_msk_neigh < 0) * size_data[0] - (idy_msk_neigh >= size_data[0]) * size_data[0]
            idx_msk_neigh = idx_msk_neigh + (idx_msk_neigh < 0) * size_data[1] - (idx_msk_neigh >= size_data[1]) * size_data[1]

            id_msk = (idy_msk, idx_msk, ich)
            id_msk_neigh = (idy_msk_neigh, idx_msk_neigh, ich)

            output[id_msk] = input_img[id_msk_neigh]
            mask[id_msk] = 0.0

        return output, mask

class N2V:
    """
    N2V Wrapper for PyTorch implementation.
    """
    def __init__(self, model_name='n2v_unet', checkpoint_dir='./checkpoints', n_channels=1, n_filters=64):
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.model_name = model_name
        self.checkpoint_dir = checkpoint_dir
        self.n_channels = n_channels
        
        # Initialize UNet
        # nch_in=1 (masked input), nch_out=1 (prediction), nch_ker=64 (filters)
        self.model = UNet(nch_in=n_channels, nch_out=n_channels, nch_ker=n_filters, norm='bnorm')
        self.model.to(self.device)
        
    def train(self, train_data, val_data=None, batch_size=4, epochs=20, lr=1e-3):
        """
        Train the N2V model.
        
        Args:
            train_data: Numpy array (D, H, W) or (D, H, W, C).
            val_data: Optional validation data.
            batch_size: Batch size.
            epochs: Number of epochs.
            lr: Learning rate.
        """
        # Prepare datasets
        train_dataset = MemoryDataset(train_data)
        loader_train = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, num_workers=0)
        
        loader_val = None
        if val_data is not None:
            val_dataset = MemoryDataset(val_data)
            loader_val = DataLoader(val_dataset, batch_size=batch_size, shuffle=False, num_workers=0)
            
        # Setup Optimizer and Loss
        optimizer = torch.optim.Adam(self.model.parameters(), lr=lr)
        criterion = nn.MSELoss().to(self.device) # N2V uses MSE or L1. Using MSE here.
        
        # Train
        self.model, history = train_loop(
            self.model, 
            loader_train, 
            loader_val, 
            optimizer, 
            criterion, 
            epochs, 
            self.device
        )
        
        # Save model
        if not os.path.exists(self.checkpoint_dir):
            os.makedirs(self.checkpoint_dir)
        save_path = os.path.join(self.checkpoint_dir, f'{self.model_name}_final.pth')
        torch.save(self.model.state_dict(), save_path)
        print(f"Model saved to {save_path}")
        
        return history

    def predict(self, data, batch_size=4):
        """
        Denoise the data.
        
        Args:
            data: Numpy array (D, H, W) or (D, H, W, C).
            batch_size: Batch size for inference.
            
        Returns:
            Denoised numpy array same shape as input.
        """
        self.model.eval()
        
        # Preprocess input
        original_ndim = data.ndim
        if original_ndim == 3: # (D, H, W)
             # (D, H, W) -> (D, H, W, 1)
            data_proc = data[..., np.newaxis]
        elif original_ndim == 2: # (H, W)
            data_proc = data[np.newaxis, ..., np.newaxis]
        else:
            data_proc = data
            
        # Normalize
        is_uint8 = data_proc.max() > 1.0
        if is_uint8:
            data_proc = data_proc.astype(np.float32) / 255.0
        else:
            data_proc = data_proc.astype(np.float32)
            
        # Inference Loop
        output_list = []
        num_samples = data_proc.shape[0]
        
        with torch.no_grad():
            for i in range(0, num_samples, batch_size):
                batch_data = data_proc[i:min(i+batch_size, num_samples)]
                # (B, H, W, C) -> (B, C, H, W)
                batch_tensor = torch.from_numpy(batch_data.transpose((0, 3, 1, 2))).to(self.device)
                
                pred = self.model(batch_tensor)
                
                # (B, C, H, W) -> (B, H, W, C)
                pred_np = pred.cpu().numpy().transpose((0, 2, 3, 1))
                output_list.append(pred_np)
                
        result = np.concatenate(output_list, axis=0)
        
        # Restore original scaling if needed
        if is_uint8:
            result = np.clip(result * 255.0, 0, 255).astype(np.uint8)
            
        # Restore original shape
        if original_ndim == 3:
            result = result[..., 0]
        elif original_ndim == 2:
            result = result[0, ..., 0]
            
        return result

    def load_model(self, model_path):
        """Load weights from a .pth file."""
        self.model.load_state_dict(torch.load(model_path, map_location=self.device))
        print(f"Loaded model from {model_path}")


# Convenience functions for easy API usage (matching documentation)

def predict_3d(noisy_volume, model_path=None, axes='ZYX', patch_size=None, **kwargs):
    """
    Denoise 3D volume using pre-trained N2V model.
    
    Args:
        noisy_volume: Noisy 3D numpy array (Z, Y, X) or (D, H, W)
        model_path: Path to trained model (.pth file). If None, uses default model.
        axes: Axis order string (e.g., 'ZYX', 'XYZ')
        patch_size: Optional patch size for large volumes
        **kwargs: Additional arguments for N2V class
    
    Returns:
        Denoised numpy array with same shape as input
    
    Example:
        >>> from ipa.processing.denoising import predict_3d
        >>> denoised = predict_3d(
        ...     noisy_volume=noisy_data,
        ...     model_path="models/n2v_3D/trained_model.pth",
        ...     axes='ZYX'
        ... )
    """
    # Initialize N2V
    n2v = N2V(**kwargs)
    
    # Load model if path provided
    if model_path and os.path.exists(model_path):
        n2v.load_model(model_path)
        print(f"Using model: {model_path}")
    else:
        print("Warning: No model loaded. Using randomly initialized model.")
        print("Please train a model first or provide model_path.")
    
    # Predict
    denoised = n2v.predict(noisy_volume)
    
    return denoised


def train_and_predict_3d(train_data, test_data=None, epochs=20, model_path=None, **kwargs):
    """
    Train N2V model and optionally predict on test data.
    
    Args:
        train_data: Noisy training data (3D numpy array)
        test_data: Optional test data for prediction
        epochs: Number of training epochs
        model_path: Path to save/load model
        **kwargs: Additional arguments
    
    Returns:
        If test_data provided: denoised result
        Otherwise: trained N2V instance
    """
    # Initialize and train
    n2v = N2V(**kwargs)
    history = n2v.train(train_data, epochs=epochs)
    
    # Save model
    if model_path:
        os.makedirs(os.path.dirname(model_path), exist_ok=True)
        torch.save(n2v.model.state_dict(), model_path)
        print(f"Model saved to {model_path}")
    
    # Predict if test data provided
    if test_data is not None:
        denoised = n2v.predict(test_data)
        return denoised
    
    return n2v

