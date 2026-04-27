#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple usage example for SIM Mitochondria Segmentation
"""

import numpy as np
from pathlib import Path
import tifffile
from quick_sim_segmentation import mito_sim_segmentation

def simple_example():
    """Simple example with minimal code"""
    
    # 1. Load your image data
    image_path = r"D:\Gitspace\ipa_full\iPA\data\sim_images\20220909_30-2-1-SIM_raw_Actin.tif"
    
    print("Loading image...")
    with tifffile.TiffFile(image_path) as tif:
        image_data = tif.asarray()
    print(f"Image shape: {image_data.shape}")
    
    # 2. Run segmentation with default parameters
    print("\nRunning segmentation...")
    results = mito_sim_segmentation(
        image_data=image_data,
        output_dir="simple_results"
    )
    
    # 3. Access results
    print(f"\nResults:")
    print(f"  Detected objects: {results['num_objects']}")
    print(f"  Processing time: {results['processing_time']:.1f} seconds")
    print(f"  Mask shape: {results['mask'].shape}")
    print(f"  Output directory: {results['output_dir']}")
    
    # 4. The main result is the binary mask
    final_mask = results['mask']  # 0-1 binary mask with original image size
    
    print(f"\nMain result saved as: {results['output_dir']}/mitochondria_mask.tif") 
    return results

def custom_parameters_example():
    """Example with custom parameters"""
    
    # Load image
    image_path = r"D:\Gitspace\ipa_full\iPA\data\sim_images\20220909_30-2-1-SIM_raw_Actin.tif"
    
    if not Path(image_path).exists():
        print(f"Image not found: {image_path}")
        return None
        
    print("Loading image...")
    with tifffile.TiffFile(image_path) as tif:
        image_data = tif.asarray()
    
    # Custom parameters for different needs
    print("\nRunning with custom parameters...")
    results = mito_sim_segmentation(
        image_data=image_data,
        output_dir="custom_results",
        
        # Processing speed control
        downsample_factor=4,        # 4x downsampling (higher = faster, lower quality)
        max_z_layers=30,            # Process max 30 Z layers (None = all layers) 
        max_crop_size=1024,         # Crop to 1024x1024 for XY (larger = slower)
        
        # Enhancement parameters
        gaussian_sigma1=1.0,        # Fine structure sigma
        gaussian_sigma2=2.0,        # Coarse structure sigma  
        laplacian_weight=0.5,       # Laplacian enhancement weight
        
        # Segmentation parameters
        min_object_size_3d=50,      # Minimum object size for 3D
        min_object_size_2d=20,      # Minimum object size for 2D
        morphology_ball_size=1,     # Morphological operations size
        
        # Output control
        verbose=True                # Print progress information
    )
    
    print(f"\nCustom segmentation completed!")
    print(f"  Objects: {results['num_objects']}")
    print(f"  Time: {results['processing_time']:.1f}s")
    
    return results

# Recommended parameter sets for different scenarios
PARAMETER_PRESETS = {
    'fast': {
        'downsample_factor': 8,
        'max_z_layers': 15,
        'max_crop_size': 512,
        'min_object_size_3d': 20,
        'gaussian_sigma1': 0.8,
        'gaussian_sigma2': 1.5,
        'description': 'Very fast processing, good for quick preview'
    },
    
    'balanced': {
        'downsample_factor': 4,
        'max_z_layers': 30,
        'max_crop_size': 1024,
        'min_object_size_3d': 50,
        'gaussian_sigma1': 1.0,
        'gaussian_sigma2': 2.0,
        'description': 'Good balance of speed and quality (default)'
    },
    
    'high_quality': {
        'downsample_factor': 2,
        'max_z_layers': None,  # All layers
        'max_crop_size': 2048,
        'min_object_size_3d': 100,
        'gaussian_sigma1': 1.2,
        'gaussian_sigma2': 2.5,
        'laplacian_weight': 0.3,
        'morphology_ball_size': 2,
        'description': 'High quality results, slower processing'
    },
    
    'fine_structures': {
        'downsample_factor': 3,
        'gaussian_sigma1': 0.7,
        'gaussian_sigma2': 1.3,
        'laplacian_weight': 0.7,
        'min_object_size_3d': 20,
        'min_object_size_2d': 10,
        'description': 'Better for detecting fine mitochondrial structures'
    }
}

def preset_example():
    """Example using parameter presets"""
    
    # Load image
    image_path = r"D:\Gitspace\ipa_full\iPA\data\sim_images\20220909_30-2-1-SIM_raw_Actin.tif"
    
    if not Path(image_path).exists():
        print(f"Image not found: {image_path}")
        return None
        
    print("Loading image...")
    with tifffile.TiffFile(image_path) as tif:
        image_data = tif.asarray()
    
    # Choose a preset
    preset_name = 'balanced'  # or 'fast', 'high_quality', 'fine_structures'
    preset = PARAMETER_PRESETS[preset_name]
    
    print(f"\nUsing '{preset_name}' preset: {preset['description']}")
    
    # Run with preset parameters
    results = mito_sim_segmentation(
        image_data=image_data,
        output_dir=f"preset_results_{preset_name}",
        **{k: v for k, v in preset.items() if k != 'description'}
    )
    
    print(f"\nPreset segmentation completed!")
    print(f"  Objects: {results['num_objects']}")
    print(f"  Time: {results['processing_time']:.1f}s")
    
    return results

if __name__ == "__main__":
    print("🔬 SIM Mitochondria Segmentation - Simple Examples")
    print("=" * 55)
    
    # Run examples
    try:
        print("\n1. Simple example with default parameters:")
        simple_example()
        
        print("\n" + "="*55)
        print("\n2. Custom parameters example:")
        custom_parameters_example()
        
        print("\n" + "="*55)
        print("\n3. Preset parameters example:")
        preset_example()
        
        print("\n" + "="*55)
        print("✅ All examples completed!")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n📋 Available parameter presets:")
    for name, preset in PARAMETER_PRESETS.items():
        print(f"   - '{name}': {preset['description']}")
    
    input("\nPress any key to exit...")