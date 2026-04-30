"""
Denoising Module

Self-supervised image denoising using Noise2Void (N2V) and Noise2Noise (N2N).
Supports 3D volume data processing with GPU acceleration.

Example:
    >>> from ipa.processing.denoising import N2V, N2N
    >>> 
    >>> # N2V: Train on single noisy images
    >>> n2v = N2V(n_filters=64)
    >>> n2v.train(noisy_data, epochs=50)
    >>> denoised = n2v.predict(noisy_data)
    >>> 
    >>> # N2N: Train on pairs of noisy images
    >>> n2n = N2N(n_filters=64)
    >>> n2n.train(noisy_data_1, noisy_data_2, epochs=50)
    >>> denoised = n2n.predict(noisy_data)
"""

from .n2v import N2V
from .n2n import N2N
import numpy as np
from typing import Optional

def predict_3d(data: np.ndarray, method: str = 'n2v', model_path: Optional[str] = None, **kwargs) -> np.ndarray:
    """
    Convenient function to denoise 3D volume data.
    
    Args:
        data: Input noisy 3D volume.
        method: Denoising method, 'n2v' or 'n2n'.
        model_path: Path to pre-trained model (optional).
        **kwargs: Additional arguments for the denoiser class.
        
    Returns:
        Denoised 3D volume.
    """
    if method.lower() == 'n2v':
        denoiser = N2V(**kwargs)
    elif method.lower() == 'n2n':
        # Note: N2N usually requires pairs for training, but for prediction we can use one instance
        denoiser = N2N(**kwargs)
    else:
        raise ValueError(f"Unsupported denoising method: {method}")
        
    if model_path:
        denoiser.load_model(model_path)
        
    return denoiser.predict(data)

__all__ = ['N2V', 'N2N', 'predict_3d']
