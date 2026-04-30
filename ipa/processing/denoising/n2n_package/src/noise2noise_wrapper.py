#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import torch
import torch.nn as nn
import os
import sys
from argparse import Namespace

# Define base path relative to this file for portability
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

from .datasets import load_dataset
from .noise2noise import Noise2Noise


def create_default_params(**kwargs):
    """Create default parameter configuration using relative paths"""
    
    default_params = {
        # Data parameters (relative to project root or current working directory)
        'train_dir': os.path.join(BASE_DIR, '..', '..', '..', 'data', 'denoised', 'train'),
        'valid_dir': os.path.join(BASE_DIR, '..', '..', '..', 'data', 'denoised', 'valid'),
        'data_dir': os.path.join(BASE_DIR, '..', '..', '..', 'data', 'denoised', 'test'),
        'ckpt_save_path': os.path.join(BASE_DIR, '..', 'models', 'n2n'),
        'load_ckpt': os.path.join(BASE_DIR, '..', 'models', 'n2n', 'best_model.pt'),
        'ckpt_overwrite': False,
        'report_interval': 100,
        'train_size': None,
        'valid_size': None,
        
        # Training hyperparameters
        'learning_rate': 0.001,
        'adam': [0.9, 0.99, 1e-8],
        'batch_size': 4,
        'nb_epochs': 50,
        'loss': 'l1',
        'cuda': torch.cuda.is_available(),
        'plot_stats': True,
        
        # Noise parameters
        'noise_type': 'gaussian',
        'noise_param': 25.0,
        'seed': None,
        'crop_size': 128,
        'clean_targets': False,
        
        # Prediction parameters
        'show_output': 1,
        'redux': False
    }
    
    # Update with user provided parameters
    default_params.update(kwargs)
    
    # Convert to Namespace object
    return Namespace(**default_params)


def train_noise2noise(train_dir=None, valid_dir=None, **kwargs):
    """Train Noise2Noise model
    
    Args:
        train_dir (str): Training dataset path
        valid_dir (str): Validation dataset path
        **kwargs: Other training parameters
        
    Returns:
        Noise2Noise: Trained model instance
    """
    
    print("=" * 50)
    print("Start Training Noise2Noise Model")
    print("=" * 50)
    
    # Create parameter configuration
    params = create_default_params(**kwargs)
    if train_dir:
        params.train_dir = train_dir
    if valid_dir:
        params.valid_dir = valid_dir
    
    # Create save directory
    os.makedirs(params.ckpt_save_path, exist_ok=True)
    
    print(f"Training Data Dir: {params.train_dir}")
    print(f"Validation Data Dir: {params.valid_dir}")
    print(f"Model Save Path: {params.ckpt_save_path}")
    print(f"Device: {'CUDA' if params.cuda and torch.cuda.is_available() else 'CPU'}")
    print(f"Noise Type: {params.noise_type}, Param: {params.noise_param}")
    print(f"Epochs: {params.nb_epochs}, Batch Size: {params.batch_size}")
    print("-" * 50)
    
    try:
        # Load datasets
        print("Loading training dataset...")
        train_loader = load_dataset(params.train_dir, params.train_size, params, shuffled=True)
        print("Loading validation dataset...")
        valid_loader = load_dataset(params.valid_dir, params.valid_size, params, shuffled=False)
        
        # Initialize and train model
        print("Initializing model...")
        n2n = Noise2Noise(params, trainable=True)
        
        print("Starting training...")
        n2n.train(train_loader, valid_loader)
        
        print("=" * 50)
        print("Model training completed!")
        print("=" * 50)
        
        return n2n
        
    except Exception as e:
        print(f"Error during training: {e}")
        raise


def predict_noise2noise(data_dir=None, model_path=None, **kwargs):
    """Predict using trained model
    
    Args:
        data_dir (str): Test data directory
        model_path (str): Model checkpoint path
        **kwargs: Other prediction parameters
        
    Returns:
        None
    """
    
    print("=" * 50)
    print("Start Noise2Noise Prediction")
    print("=" * 50)
    
    # Create parameter configuration
    params = create_default_params(**kwargs)
    if data_dir:
        params.data_dir = data_dir
    if model_path:
        params.load_ckpt = model_path
    
    # Set prediction specific parameters
    params.redux = False
    params.clean_targets = True
    
    print(f"Test Data Dir: {params.data_dir}")
    print(f"Model Path: {params.load_ckpt}")
    print(f"Device: {'CUDA' if params.cuda and torch.cuda.is_available() else 'CPU'}")
    print(f"Noise Type: {params.noise_type}, Param: {params.noise_param}")
    print("-" * 50)
    
    try:
        # Check if model file exists
        if not os.path.exists(params.load_ckpt):
            raise FileNotFoundError(f"Model file not found: {params.load_ckpt}")
        
        # Initialize model
        print("Initializing model...")
        n2n = Noise2Noise(params, trainable=False)
        
        # Load test data
        print("Loading test data...")
        test_loader = load_dataset(params.data_dir, 0, params, shuffled=False, single=True)
        
        # Load model weights
        print("Loading model weights...")
        n2n.load_model(params.load_ckpt)
        
        # Predict
        print("Starting prediction...")
        n2n.test(test_loader, show=params.show_output)
        
        print("=" * 50)
        print("Prediction completed!")
        print("=" * 50)
        
    except Exception as e:
        print(f"Error during prediction: {e}")
        raise


def train_and_predict(train_dir=None, valid_dir=None, test_dir=None, **kwargs):
    """Train model and then predict
    
    Args:
        train_dir (str): Training dataset path
        valid_dir (str): Validation dataset path
        test_dir (str): Test dataset path
        **kwargs: Other parameters
        
    Returns:
        Noise2Noise: Trained model instance
    """
    
    print("=" * 50)
    print("Noise2Noise Full Pipeline: Train + Predict")
    print("=" * 50)
    
    # Train model
    n2n = train_noise2noise(train_dir=train_dir, valid_dir=valid_dir, **kwargs)
    
    # Set prediction to use the trained model
    params = create_default_params(**kwargs)
    model_path = os.path.join(params.ckpt_save_path, 'best_model.pt')
    
    if os.path.exists(model_path):
        # Predict
        predict_noise2noise(data_dir=test_dir or params.data_dir, model_path=model_path, **kwargs)
    else:
        print("Warning: Model file not found after training, skipping prediction")
        print(f"Expected model path: {model_path}")
    
    return n2n


def get_available_devices():
    """Get available compute devices"""
    
    devices = {
        'cpu': True,
        'cuda': torch.cuda.is_available(),
        'cuda_devices': []
    }
    
    if devices['cuda']:
        devices['cuda_devices'] = [torch.cuda.get_device_name(i) for i in range(torch.cuda.device_count())]
    
    return devices


if __name__ == '__main__':
    """Test function"""
    
    print("Noise2Noise Wrapper Test")
    print("Available devices:", get_available_devices())
    
    # Add simple test call here
    # e.g.: train_noise2noise(nb_epochs=1, batch_size=2)
