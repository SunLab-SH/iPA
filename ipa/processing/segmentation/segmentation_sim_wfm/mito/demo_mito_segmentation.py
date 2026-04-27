#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Demo script for SIM Mitochondria Segmentation
Shows different parameter configurations and usage examples
"""

import numpy as np
from pathlib import Path
import tifffile
import time
from quick_sim_segmentation import mito_sim_segmentation

def demo_basic_usage():
    """Basic usage example"""
    print("="*60)
    print("🔬 Demo 1: Basic Usage")
    print("="*60)
    
    # Load test image
    sim_path = r"D:\Gitspace\ipa_full\iPA\data\sim_images\20220909_30-2-1-SIM_raw_Actin.tif"
    
    if not Path(sim_path).exists():
        print(f"❌ Test image not found: {sim_path}")
        return None
    
    print("📂 Loading image data...")
    with tifffile.TiffFile(sim_path) as tif:
        image_data = tif.asarray()
    
    print(f"   Loaded image shape: {image_data.shape}")
    print(f"   Data type: {image_data.dtype}")
    print(f"   Memory usage: {image_data.nbytes / 1024**3:.2f} GB")
    
    # Run with default parameters
    results = mito_sim_segmentation(
        image_data=image_data,
        output_dir="demo_results/basic",
        verbose=True
    )
    
    print(f"\n✅ Basic segmentation completed!")
    print(f"   Objects detected: {results['num_objects']}")
    print(f"   Processing time: {results['processing_time']:.1f}s")
    print(f"   Output mask shape: {results['mask'].shape}")
    
    return results

def demo_fast_processing():
    """Fast processing with aggressive downsampling"""
    print("\n" + "="*60)
    print("⚡ Demo 2: Fast Processing (High Downsampling)")
    print("="*60)
    
    # Load test image
    sim_path = r"D:\Gitspace\ipa_full\iPA\data\sim_images\20220909_30-2-1-SIM_raw_Actin.tif"
    
    if not Path(sim_path).exists():
        print(f"❌ Test image not found: {sim_path}")
        return None
    
    print("📂 Loading image data...")
    with tifffile.TiffFile(sim_path) as tif:
        image_data = tif.asarray()
    
    # Fast processing parameters
    results = mito_sim_segmentation(
        image_data=image_data,
        output_dir="demo_results/fast",
        downsample_factor=8,           # Aggressive downsampling
        max_z_layers=20,               # Process only 20 Z layers
        max_crop_size=512,             # Smaller crop size
        min_object_size_3d=20,         # Smaller minimum object size
        gaussian_sigma1=0.8,           # Slightly smaller sigma
        gaussian_sigma2=1.5,
        verbose=True
    )
    
    print(f"\n⚡ Fast segmentation completed!")
    print(f"   Objects detected: {results['num_objects']}")
    print(f"   Processing time: {results['processing_time']:.1f}s")
    
    return results

def demo_high_quality():
    """High quality processing with conservative parameters"""
    print("\n" + "="*60)
    print("🎯 Demo 3: High Quality Processing")
    print("="*60)
    
    # Load test image
    sim_path = r"D:\Gitspace\ipa_full\iPA\data\sim_images\20220909_30-2-1-SIM_raw_Actin.tif"
    
    if not Path(sim_path).exists():
        print(f"❌ Test image not found: {sim_path}")
        return None
    
    print("📂 Loading image data...")
    with tifffile.TiffFile(sim_path) as tif:
        image_data = tif.asarray()
    
    # High quality parameters
    results = mito_sim_segmentation(
        image_data=image_data,
        output_dir="demo_results/high_quality",
        downsample_factor=2,           # Minimal downsampling
        max_z_layers=None,             # Process all Z layers
        max_crop_size=2048,            # Larger crop size
        min_object_size_3d=100,        # Larger minimum object size
        min_object_size_2d=50,
        gaussian_sigma1=1.2,           # Slightly larger sigma for smoother results
        gaussian_sigma2=2.5,
        laplacian_weight=0.3,          # Reduced Laplacian weight
        morphology_ball_size=2,        # Larger morphological operations
        verbose=True
    )
    
    print(f"\n🎯 High quality segmentation completed!")
    print(f"   Objects detected: {results['num_objects']}")
    print(f"   Processing time: {results['processing_time']:.1f}s")
    
    return results

def demo_parameter_comparison():
    """Compare different parameter settings"""
    print("\n" + "="*60)
    print("📊 Demo 4: Parameter Comparison")
    print("="*60)
    
    # Load test image
    sim_path = r"D:\Gitspace\ipa_full\iPA\data\sim_images\20220909_30-2-1-SIM_raw_Actin.tif"
    
    if not Path(sim_path).exists():
        print(f"❌ Test image not found: {sim_path}")
        return None
    
    print("📂 Loading image data...")
    with tifffile.TiffFile(sim_path) as tif:
        image_data = tif.asarray()
    
    # Test different gaussian sigma combinations
    sigma_configs = [
        (0.8, 1.5, "Fine details"),
        (1.0, 2.0, "Default"),
        (1.5, 3.0, "Smooth structures")
    ]
    
    results_comparison = []
    
    for sigma1, sigma2, description in sigma_configs:
        print(f"\n🔍 Testing {description} (σ1={sigma1}, σ2={sigma2})...")
        
        start_time = time.time()
        results = mito_sim_segmentation(
            image_data=image_data,
            output_dir=f"demo_results/comparison_sigma{sigma1}_{sigma2}",
            downsample_factor=6,  # Medium downsampling for speed
            gaussian_sigma1=sigma1,
            gaussian_sigma2=sigma2,
            max_z_layers=15,      # Limit Z layers for speed
            verbose=False         # Reduce output
        )
        
        results_comparison.append({
            'config': description,
            'sigma1': sigma1,
            'sigma2': sigma2,
            'objects': results['num_objects'],
            'time': results['processing_time']
        })
        
        print(f"   → Objects: {results['num_objects']}, Time: {results['processing_time']:.1f}s")
    
    # Print comparison summary
    print(f"\n📊 Parameter Comparison Summary:")
    print("-" * 60)
    for result in results_comparison:
        print(f"   {result['config']:<15} | Objects: {result['objects']:<3} | Time: {result['time']:>5.1f}s")
    
    return results_comparison

def demo_batch_processing():
    """Demo batch processing multiple images or regions"""
    print("\n" + "="*60)
    print("🔄 Demo 5: Batch Processing Simulation")
    print("="*60)
    
    # Load test image
    sim_path = r"D:\Gitspace\ipa_full\iPA\data\sim_images\20220909_30-2-1-SIM_raw_Actin.tif"
    
    if not Path(sim_path).exists():
        print(f"❌ Test image not found: {sim_path}")
        return None
    
    print("📂 Loading image data...")
    with tifffile.TiffFile(sim_path) as tif:
        full_image = tif.asarray()
    
    # Simulate batch processing by splitting image into regions
    if len(full_image.shape) == 3:
        z_chunks = 3  # Split into 3 Z regions
        z_size = full_image.shape[0] // z_chunks
        
        batch_results = []
        total_objects = 0
        total_time = 0
        
        for i in range(z_chunks):
            z_start = i * z_size
            z_end = (i + 1) * z_size if i < z_chunks - 1 else full_image.shape[0]
            
            chunk_data = full_image[z_start:z_end]
            print(f"\n🔍 Processing chunk {i+1}/{z_chunks} (Z: {z_start}-{z_end})...")
            
            results = mito_sim_segmentation(
                image_data=chunk_data,
                output_dir=f"demo_results/batch_chunk_{i+1}",
                downsample_factor=4,
                verbose=False
            )
            
            batch_results.append({
                'chunk': i+1,
                'z_range': f"{z_start}-{z_end}",
                'objects': results['num_objects'],
                'time': results['processing_time']
            })
            
            total_objects += results['num_objects']
            total_time += results['processing_time']
            
            print(f"   → Chunk {i+1}: {results['num_objects']} objects, {results['processing_time']:.1f}s")
        
        print(f"\n🔄 Batch Processing Summary:")
        print("-" * 50)
        for result in batch_results:
            print(f"   Chunk {result['chunk']} (Z: {result['z_range']:<10}) | Objects: {result['objects']:<3} | Time: {result['time']:>5.1f}s")
        print("-" * 50)
        print(f"   Total: {total_objects} objects, {total_time:.1f}s")
        
        return batch_results
    
    else:
        print("   Image is 2D, skipping batch processing demo")
        return None

def main():
    """Run all demos"""
    print("🚀 SIM Mitochondria Segmentation - Demo Suite")
    print("=" * 60)
    
    # Create output directory
    Path("demo_results").mkdir(exist_ok=True)
    
    try:
        # Run demos
        demo_basic_usage()
        demo_fast_processing()
        demo_high_quality()
        demo_parameter_comparison()
        demo_batch_processing()
        
        print("\n" + "="*60)
        print("🎉 All demos completed successfully!")
        print("="*60)
        print("\n📁 Check the following directories for results:")
        print("   - demo_results/basic/")
        print("   - demo_results/fast/")
        print("   - demo_results/high_quality/")
        print("   - demo_results/comparison_*/")
        print("   - demo_results/batch_chunk_*/")
        
    except Exception as e:
        print(f"\n❌ Demo failed with error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
    input("\nPress any key to exit...")