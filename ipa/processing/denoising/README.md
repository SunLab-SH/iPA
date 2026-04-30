# Denoising Module

Self-supervised image denoising using Noise2Void (N2V) and Noise2Noise (N2N).

## Overview

This module provides two state-of-the-art self-supervised denoising methods:

- **Noise2Void (N2V)**: Trains on single noisy images without requiring clean ground truth
- **Noise2Noise (N2N)**: Trains on pairs of noisy images with different noise realizations

Both methods support 3D volume data and GPU acceleration.

## Quick Start

### N2V: Single Image Denoising

```python
from ipa.processing.denoising import N2V
import numpy as np

# Load your noisy 3D data (shape: D x H x W)
noisy_data = np.load('noisy_volume.npy')

# Create denoiser
n2v = N2V(n_channels=1, n_filters=64)

# Train
n2v.train(noisy_data, epochs=50, batch_size=4)

# Predict
denoised = n2v.predict(noisy_data)

# Save model
n2v.save_model('n2v_model.pth')
```

### N2N: Paired Image Denoising

```python
from ipa.processing.denoising import N2N
import numpy as np

# Load pairs of noisy images
noisy_data_1 = np.load('noisy_volumes_1.npy')
noisy_data_2 = np.load('noisy_volumes_2.npy')  # Same content, different noise

# Create denoiser
n2n = N2N(n_channels=1, n_filters=64)

# Train
n2n.train(noisy_data_1, noisy_data_2, epochs=50, batch_size=4)

# Predict
denoised = n2n.predict(noisy_data_1)

# Save model
n2n.save_model('n2n_model.pth')
```

## Examples

See the example scripts in `examples/examples_sxt/`:

- `demo_SXT_n2v_train.py` - Train N2V model
- `demo_SXT_n2v_predict.py` - Use trained N2V model
- `demo_SXT_n2n_train.py` - Train N2N model  
- `demo_SXT_n2n_predict.py` - Use trained N2N model

Run examples:

```bash
cd examples/examples_sxt
python demo_SXT_n2v_train.py
python demo_SXT_n2v_predict.py
```

## API Reference

### N2V Class

#### `__init__(n_channels=1, n_filters=64, device=None)`

Initialize N2V denoiser.

**Parameters:**
- `n_channels` (int): Number of input channels (default: 1)
- `n_filters` (int): Number of convolutional filters (default: 64)
- `device` (str): Device to run on ('cuda' or 'cpu', auto-detected if None)

#### `train(train_data, val_data=None, epochs=50, batch_size=4, lr=1e-3, mask_ratio=0.195, patch_size=(64, 64))`

Train the N2V model.

**Parameters:**
- `train_data` (np.ndarray): Training data, shape (D, H, W) or (D, H, W, C)
- `val_data` (np.ndarray, optional): Validation data
- `epochs` (int): Number of training epochs (default: 50)
- `batch_size` (int): Batch size (default: 4)
- `lr` (float): Learning rate (default: 1e-3)
- `mask_ratio` (float): Ratio of pixels to mask (default: 0.195)
- `patch_size` (tuple): Patch size for training (default: (64, 64))

#### `predict(data, batch_size=4)`

Predict denoised output.

**Parameters:**
- `data` (np.ndarray): Noisy input data
- `batch_size` (int): Batch size for prediction (default: 4)

**Returns:**
- np.ndarray: Denoised data

#### `save_model(path)`

Save trained model to file.

#### `load_model(path)`

Load trained model from file.

### N2N Class

Similar API to N2V, but `train()` requires two datasets:

```python
n2n.train(noisy_data_1, noisy_data_2, val_data_1=None, val_data_2=None, 
          epochs=50, batch_size=4, lr=1e-3, patch_size=(64, 64), loss_type='l1')
```

## Tips

1. **Data Preparation**: Normalize your data to [0, 1] range before training
2. **Patch Size**: Use smaller patches (e.g., 64x64) for faster training
3. **Epochs**: Start with 50 epochs, increase if needed
4. **Batch Size**: Adjust based on your GPU memory (smaller = less memory)
5. **Validation**: Always use validation data to monitor overfitting

## References

- Noise2Void: [Krull et al., CVPR 2019](https://openaccess.thecvf.com/content_CVPR_2019/html/Krull_Noise2Void_Learning_Denoising_From_Single_Noisy_Images_CVPR_2019_paper.html)
- Noise2Noise: [Lehtinen et al., ICML 2018](http://proceedings.mlr.press/v80/lehtinen18a.html)
