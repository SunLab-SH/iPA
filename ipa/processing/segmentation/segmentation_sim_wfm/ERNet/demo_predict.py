#!/usr/bin/env python3
"""
Demo script for ERNet prediction (inference only) with external data and parameters
This script is for pure prediction without ground truth labels
"""

import os
import glob
import numpy as np
from types import SimpleNamespace
from PIL import Image
import torch
import torch.nn as nn

# Import from local modules
try:
    from .models import GetModel
    from .datahandler import toTensor, toPIL
except ImportError:
    import sys
    import os
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from models import GetModel
    from datahandler import toTensor, toPIL


def create_pred_config(
    model_path,
    image_size=1000,
    model='rcan',
    nch_in=1,
    nch_out=2,
    n_resgroups=5,
    n_resblocks=10,
    use_cpu=False,
    undomulti=False
):
    """
    Create configuration object for ERNet prediction
    
    Args:
        model_path: Path to model weights file
        image_size: Size for image processing (default: 1000)
        model: Model type (default: 'rcan')
        nch_in: Number of input channels (default: 1)
        nch_out: Number of output channels (default: 2)
        n_resgroups: Number of residual groups (default: 5)
        n_resblocks: Number of residual blocks (default: 10)
        use_cpu: Use CPU instead of GPU (default: False)
        undomulti: Remove DataParallel wrapper (default: False)
    
    Returns:
        SimpleNamespace: Configuration object
    """
    
    opt = SimpleNamespace()
    
    # Core parameters
    opt.weights = model_path
    opt.imageSize = image_size
    opt.model = model
    opt.nch_in = nch_in
    opt.nch_out = nch_out
    opt.n_resgroups = n_resgroups
    opt.n_resblocks = n_resblocks
    opt.cpu = use_cpu
    opt.undomulti = undomulti
    
    # Additional required parameters
    opt.n_feats = 64
    opt.reduction = 16
    opt.narch = 0
    opt.multigpu = False
    
    return opt


def remove_dataparallel_wrapper(state_dict):
    """Remove DataParallel wrapper from state dict"""
    from collections import OrderedDict
    
    new_state_dict = OrderedDict()
    for k, vl in state_dict.items():
        name = k[7:] if k.startswith('module.') else k  # remove 'module.' prefix
        new_state_dict[name] = vl
    
    return new_state_dict


def load_model(config):
    """
    Load the trained model
    
    Args:
        config: Configuration object
    
    Returns:
        torch.nn.Module: Loaded model
    """
    
    # Check if weights file exists
    if not os.path.exists(config.weights):
        raise FileNotFoundError(f"Model weights file not found: {config.weights}")
    
    # Get model
    net = GetModel(config)
    
    # Load checkpoint
    print(f"Loading model from: {config.weights}")
    checkpoint = torch.load(config.weights, map_location='cpu' if config.cpu else None)
    
    # Handle DataParallel wrapper if needed
    if config.undomulti:
        checkpoint['state_dict'] = remove_dataparallel_wrapper(checkpoint['state_dict'])
    
    # Load state dict
    net.load_state_dict(checkpoint['state_dict'])
    
    # Set to evaluation mode
    net.eval()
    
    # Move to device
    if config.cpu:
        net = net.cpu()
    else:
        net = net.cuda()
    
    print("Model loaded successfully!")
    return net


def predict_single_image(model, img, config, save_path=None):
    """
    Predict segmentation for a single image
    
    Args:
        model: Loaded PyTorch model
        img: numpy array of the image (already loaded and normalized)
        config: Configuration object
        save_path: Optional path to save result
    
    Returns:
        numpy.ndarray: Prediction result
    """
    
    # Handle color images
    if len(img.shape) > 2:
        print("Converting color image to grayscale")
        img = img[:, :, 0]
    
    print(f"Processing image, shape: {img.shape}")
    
    h, w = img.shape
    imageSize = config.imageSize
    
    # Auto-adjust image size if needed
    if imageSize == 0:
        imageSize = 250
        while imageSize + 250 < h and imageSize + 250 < w:
            imageSize += 250
        print(f'Auto-adjusted imageSize to: {imageSize}')
    
    # Create image patches for processing
    images = []
    images.append(img[:imageSize, :imageSize])
    images.append(img[h-imageSize:, :imageSize])
    images.append(img[:imageSize, w-imageSize:])
    images.append(img[h-imageSize:, w-imageSize:])
    
    # Process each patch
    proc_images = []
    for idx, sub_img in enumerate(images):
        print(f'\rProcessing patch [{idx+1}/4]', end='', flush=True)
        
        # Convert to PIL and tensor
        pil_sub_img = Image.fromarray((sub_img * 255).astype('uint8'))
        sub_tensor = toTensor(pil_sub_img)
        sub_tensor = sub_tensor.unsqueeze(0)
        
        # Predict
        with torch.no_grad():
            if config.cpu:
                sr = model(sub_tensor)
            else:
                sr = model(sub_tensor.cuda())
            sr = sr.cpu()
            
            # Apply softmax and get prediction
            m = nn.LogSoftmax(dim=1)
            sr = m(sr)
            sr = sr.argmax(dim=1, keepdim=True)
            
            # Remove batch dimension and convert to PIL
            sr_squeezed = sr.squeeze(0)  # Remove batch dimension: [1, 1, H, W] -> [1, H, W]
            pil_sr_img = toPIL(sr_squeezed.float() / (config.nch_out - 1))
            proc_images.append(pil_sr_img)
    
    print()  # New line after processing
    
    # Stitch patches back together
    img1, img2, img3, img4 = proc_images
    
    woffset = (2 * imageSize - w) // 2
    hoffset = (2 * imageSize - h) // 2
    
    img1 = np.array(img1)[:imageSize-hoffset, :imageSize-woffset]
    img3 = np.array(img3)[:imageSize-hoffset, woffset:]
    top = np.concatenate((img1, img3), axis=1)
    
    img2 = np.array(img2)[hoffset:, :imageSize-woffset]
    img4 = np.array(img4)[hoffset:, woffset:]
    bot = np.concatenate((img2, img4), axis=1)
    
    result = np.concatenate((top, bot), axis=0)
    
    # Clean up borders
    result[:10, :] = 0
    result[-10:, :] = 0
    result[:, :10] = 0
    result[:, -10:] = 0
    
    # Save result if path provided
    if save_path:
        Image.fromarray(result.astype('uint8')).save(save_path)
        print(f"Result saved to: {save_path}")
    
    return result





def demo_prediction():
    """
    Demo function for ERNet prediction with single image
    """
    
    # Get script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Define paths
    weights_path = os.path.join(script_dir, 'pretrained_model', 'final.pth')
    image_path = os.path.join(script_dir, 'testdata', 'input', '190624_COS7_3_LAMP1-mCherry_ER-YFP_optokin_2s_40x_1pt5_01_2_MMStack_Pos0-10000.png')
    output_path = os.path.join(script_dir, 'prediction_output', 'result.png')
    
    # Check paths
    if not os.path.exists(weights_path):
        print(f"Error: Model weights not found at: {weights_path}")
        return
    
    if not os.path.exists(image_path):
        print(f"Error: Input image not found at: {image_path}")
        return
    
    # Create output directory
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Create configuration
    config = create_pred_config(
        model_path=weights_path,
        image_size=1000,
        model='rcan',
        nch_in=1,
        nch_out=2,
        n_resgroups=5,
        n_resblocks=10,
        use_cpu=False,  # Set to True if no GPU
        undomulti=False
    )
    
    print("=" * 60)
    print("ERNet Single Image Prediction Demo")
    print("=" * 60)
    print(f"Model weights: {config.weights}")
    print(f"Input image: {image_path}")
    print(f"Output path: {output_path}")
    print(f"Image size: {config.imageSize}")
    print(f"Model: {config.model}")
    print(f"Using CPU: {config.cpu}")
    print("=" * 60)
    
    try:
        # Load model
        model = load_model(config)
        
        # Load and preprocess image
        img = np.array(Image.open(image_path)) / 255.0
        print(f"Loaded image with shape: {img.shape}")
        
        # Run prediction
        result = predict_single_image(model, img, config, output_path)
        
        print("\nPrediction completed successfully!")
        print(f"Result saved to: {output_path}")
        
    except Exception as e:
        print(f"Error: {str(e)}")
        raise


def predict_custom_image(model_path, img_array, output_path, **kwargs):
    """
    Convenience function for custom prediction with external image array
    
    Args:
        model_path: Path to model weights
        img_array: numpy array of the image (normalized 0-1)
        output_path: Path to save result
        **kwargs: Additional model parameters
    
    Example:
        img = np.array(Image.open('image.png')) / 255.0
        predict_custom_image(
            model_path='model.pth',
            img_array=img,
            output_path='result.png',
            image_size=1000,
            use_cpu=False
        )
    """
    
    # Create config
    config = create_pred_config(model_path, **kwargs)
    
    # Load model
    model = load_model(config)
    
    # Predict
    result = predict_single_image(model, img_array, config, output_path)
    
    return result


if __name__ == '__main__':
    # Run demo
    demo_prediction()
    
    # Example usage:
    # img = np.array(Image.open('my_image.png')) / 255.0
    # predict_custom_image(
    #     model_path='pretrained_model/final.pth',
    #     img_array=img,
    #     output_path='my_result.png',
    #     image_size=1000,
    #     use_cpu=False
    # )