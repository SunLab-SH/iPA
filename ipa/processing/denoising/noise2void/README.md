# Noise2Void (PyTorch Implementation)

## Overview
Noise2Void (N2V) is a self-supervised training scheme that allows learning denoising from single noisy images. This implementation is based on PyTorch and integrated into iPA.

## Usage

### Python API

```python
from ipa.processing.denoising.noise2void.n2v_wrapper import N2V

# Initialize model
n2v = N2V(model_name='my_n2v_model', checkpoint_dir='./checkpoints')

# Train
# train_data can be a numpy array (e.g., [D, H, W] or [N, H, W])
history = n2v.train(train_data, batch_size=4, epochs=20, lr=1e-3)

# Predict
denoised_data = n2v.predict(noisy_data, batch_size=4)
```

## References
- [Noise2Void - Learning Denoising from Single Noisy Images](https://arxiv.org/abs/1811.10980)
