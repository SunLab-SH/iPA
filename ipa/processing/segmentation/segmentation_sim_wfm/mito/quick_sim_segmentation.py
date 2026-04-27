#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SIM Image Mitochondria Segmentation - Fast Version
Optimized for large image processing with high speed
"""

import numpy as np
from pathlib import Path
import tifffile
from scipy import ndimage as ndi
from skimage import transform, filters, morphology, measure
import time
import warnings
warnings.filterwarnings('ignore')

def mito_sim_segmentation(im_path, 
                         output_dir=None, 
                         downsample_factor=4,
                         crop_size=1024,
                         z_range_ratio=1.0,
                         sigma1=1.0,
                         sigma2=2.0,
                         laplacian_weight=0.5,
                         min_object_size_3d=50,
                         min_object_size_2d=20,
                         morphology_ball_size=1,
                         morphology_disk_size=1):
    """
    Fast SIM image segmentation (preserving original image size output)
    
    Parameters:
    -----------
    im_path : str
        Path to input SIM image
    output_dir : str, optional
        Output directory path (default: same as input with suffix)
    downsample_factor : int, default=4
        Downsampling factor for processing speed (1=no downsampling, 4=4x faster)
    crop_size : int, default=1024
        Size for central crop region (larger = more coverage, slower)
    z_range_ratio : float, default=1.0
        Ratio of Z slices to process (1.0=all slices, 0.5=middle half)
    sigma1 : float, default=1.0
        First Gaussian filter sigma for vessel enhancement
    sigma2 : float, default=2.0
        Second Gaussian filter sigma for vessel enhancement
    laplacian_weight : float, default=0.5
        Weight for Laplacian enhancement component
    min_object_size_3d : int, default=50
        Minimum object size for 3D morphological cleanup
    min_object_size_2d : int, default=20
        Minimum object size for 2D morphological cleanup
    morphology_ball_size : int, default=1
        Ball size for 3D morphological opening
    morphology_disk_size : int, default=1
        Disk size for 2D morphological opening
        
    Returns:
    --------
    bool : True if successful
    """
    
    print("🚀 Fast SIM Mitochondria Segmentation")
    print("=" * 40)
    
    start_time = time.time()
    
    # Setup output directory
    if output_dir is None:
        output_dir = Path(im_path).parent / "quick_mitochondria_results"
    Path(output_dir).mkdir(exist_ok=True)
    
    # Step 1: Fast image loading and preprocessing
    print("\n📂 Step 1: Fast Image Loading")
    with tifffile.TiffFile(im_path) as tif:
        original_image = tif.asarray()
    
    print(f"   Original size: {original_image.shape}")
    print(f"   Original data size: {original_image.nbytes / 1024**3:.2f} GB")
    
    # Save original shape and crop information
    original_shape = original_image.shape
    crop_info = None
    
    # Create processing image copy
    image_array = original_image.copy()
    
    # Select central region and downsample
    if len(image_array.shape) == 3:
        # Select Z slices based on z_range_ratio
        z_center = image_array.shape[0] // 2
        z_range = max(1, int(image_array.shape[0] * z_range_ratio))
        z_start = max(0, z_center - z_range // 2)
        z_end = min(image_array.shape[0], z_start + z_range)
        image_array = image_array[z_start:z_end]
        print(f"   Using {z_range}/{original_shape[0]} Z slices (ratio: {z_range_ratio})")
        
        # Record Z cropping info
        z_crop_info = (z_start, z_end, original_shape[0])
        
        # Crop central region
        h, w = image_array.shape[1:]
        if h > 1024 or w > 1024:  # Increase crop size
            center_h, center_w = h // 2, w // 2
            crop_size = 1024
            start_h = max(0, center_h - crop_size // 2)
            end_h = min(h, start_h + crop_size)
            start_w = max(0, center_w - crop_size // 2)
            end_w = min(w, start_w + crop_size)
            
            image_array = image_array[:, start_h:end_h, start_w:end_w]
            crop_info = (start_h, end_h, start_w, end_w, h, w)
            print(f"   Crop region: H[{start_h}:{end_h}], W[{start_w}:{end_w}]")
        else:
            crop_info = None
            z_crop_info = (z_start, z_end, original_shape[0])
        
        print(f"   Shape after cropping: {image_array.shape}")
        
        # Downsample
        if downsample_factor > 1:
            downsampled_shape = (image_array.shape[0],
                               image_array.shape[1] // downsample_factor,
                               image_array.shape[2] // downsample_factor)
            image_array = transform.resize(image_array, downsampled_shape, 
                                         preserve_range=True, 
                                         anti_aliasing=True).astype('uint16')
    else:
        crop_info = None
        z_crop_info = None
    
    print(f"   Final processing size: {image_array.shape}")
    print(f"   Processing data size: {image_array.nbytes / 1024**2:.1f} MB")
    
    # Step 2: Fast vessel enhancement (simplified Frangi)
    print("\n🔍 Step 2: Fast Vessel Enhancement")
    
    if len(image_array.shape) == 3:
        # Process 3D
        enhanced = np.zeros_like(image_array, dtype='float32')
        
        for i, slice_2d in enumerate(image_array):
            if i % 5 == 0:
                print(f"   Processing slice {i+1}/{len(image_array)}")
            
            # Simplified tubular structure enhancement
            gaussian1 = ndi.gaussian_filter(slice_2d.astype('float32'), sigma=sigma1)
            gaussian2 = ndi.gaussian_filter(slice_2d.astype('float32'), sigma=sigma2)
            
            # Calculate Laplacian response
            laplacian = ndi.laplace(gaussian1)
            
            # Combine different scales
            enhanced_slice = gaussian1 - gaussian2 + np.abs(laplacian) * laplacian_weight
            enhanced_slice[enhanced_slice < 0] = 0
            
            enhanced[i] = enhanced_slice
    else:
        # Process 2D
        gaussian1 = ndi.gaussian_filter(image_array.astype('float32'), sigma=sigma1)
        gaussian2 = ndi.gaussian_filter(image_array.astype('float32'), sigma=sigma2)
        laplacian = ndi.laplace(gaussian1)
        enhanced = gaussian1 - gaussian2 + np.abs(laplacian) * laplacian_weight
        enhanced[enhanced < 0] = 0
    
    print("   ✅ Vessel enhancement completed")
    
    # Step 3: Fast segmentation
    print("\n🎯 Step 3: Fast Segmentation")
    
    # Adaptive threshold
    threshold = filters.threshold_otsu(enhanced[enhanced > 0])
    print(f"   Using threshold: {threshold:.3f}")
    
    # Binarization
    binary = enhanced > threshold
    
    # Morphological cleanup
    if len(binary.shape) == 3:
        # 3D morphological operations
        binary = morphology.remove_small_objects(binary, min_size=min_object_size_3d)
        binary = ndi.binary_fill_holes(binary)
        binary = morphology.binary_opening(binary, morphology.ball(morphology_ball_size))
    else:
        # 2D morphological operations
        binary = morphology.remove_small_objects(binary, min_size=min_object_size_2d)
        binary = ndi.binary_fill_holes(binary)
        binary = morphology.binary_opening(binary, morphology.disk(morphology_disk_size))
    
    # Connected component labeling
    labels = measure.label(binary)
    
    # Statistics
    num_objects = len(np.unique(labels)) - 1  # Exclude background
    print(f"   Detected {num_objects} mitochondria structures")
    
    # Step 4: Map results back to original image size
    print("\n🔄 Step 4: Mapping Results to Original Size")
    
    # Create result arrays with original image size
    final_mask = np.zeros(original_shape, dtype='uint8')
    
    if len(original_shape) == 3:
        # Upsample to pre-crop size
        if downsample_factor > 1:
            # Calculate upsampled shape
            upsampled_shape = (image_array.shape[0],
                             image_array.shape[1] * downsample_factor,
                             image_array.shape[2] * downsample_factor)
            
            print(f"   Upsampling to size: {upsampled_shape}")
            
            # Upsample binary mask
            binary_upsampled = transform.resize(binary.astype('float32'), upsampled_shape, 
                                              preserve_range=True, 
                                              anti_aliasing=False, 
                                              order=0) > 0.5
        else:
            binary_upsampled = binary
        
        # Map back to original image region
        if crop_info is not None and z_crop_info is not None:
            start_h, end_h, start_w, end_w, orig_h, orig_w = crop_info
            z_start, z_end, orig_z = z_crop_info
            
            # Ensure size matching
            actual_h = min(binary_upsampled.shape[1], end_h - start_h)
            actual_w = min(binary_upsampled.shape[2], end_w - start_w)
            actual_z = min(binary_upsampled.shape[0], z_end - z_start)
            
            final_mask[z_start:z_start+actual_z, 
                      start_h:start_h+actual_h, 
                      start_w:start_w+actual_w] = binary_upsampled[:actual_z, :actual_h, :actual_w]
            
            print(f"   Mapped to original region: Z[{z_start}:{z_start+actual_z}], H[{start_h}:{start_h+actual_h}], W[{start_w}:{start_w+actual_w}]")
        elif z_crop_info is not None:
            # No cropping, directly place in Z range
            z_start, z_end, orig_z = z_crop_info
            actual_z = min(binary_upsampled.shape[0], z_end - z_start)
            
            final_mask[z_start:z_start+actual_z] = binary_upsampled[:actual_z]
            
            print(f"   Mapped to original Z layers: [{z_start}:{z_start+actual_z}]")
    
    # Save final binary mask (0-1 values)
    tifffile.imwrite(Path(output_dir) / "mitochondria_mask.tif", final_mask)
    
    print(f"   ✅ Final mask saved with original size: {original_shape}")
    
    # Count final objects
    final_num_objects = len(np.unique(measure.label(final_mask))) - 1
    print(f"   Final detected objects: {final_num_objects}")
    
    # Step 5: Generate middle slice visualization
    print("\n🖼️  Step 5: Generate Middle Slice Visualization")
    
    # Get middle slice from original and final mask
    if len(original_shape) == 3:
        mid_z = original_shape[0] // 2
        original_slice = original_image[mid_z]
        mask_slice = final_mask[mid_z]
        
        print(f"   Middle slice: Z={mid_z}")
    else:
        original_slice = original_image
        mask_slice = final_mask
    
    # Save middle slice visualization
    tifffile.imwrite(Path(output_dir) / "middle_slice_original.tif", original_slice)
    tifffile.imwrite(Path(output_dir) / "middle_slice_mask.tif", mask_slice * 255)  # 0-255 for visualization
    
    print(f"   ✅ Middle slice visualization saved")
    
    total_time = time.time() - start_time
    print(f"\n🎉 Fast segmentation completed!")
    print(f"📁 Results saved in: {output_dir}")
    print(f"⏱️  Total time: {total_time:.1f} seconds")
    print(f"📊 Detected {final_num_objects} mitochondria structures")
    
    return True

def main():
    """Main function"""
    print("⚡ SIM Mitochondria Fast Segmentation Tool")
    
    # SIM image path
    sim_path = r"D:\Gitspace\ipa_full\iPA\data\sim_images\20220909_30-2-1-SIM_raw_Actin.tif"
    
    if not Path(sim_path).exists():
        print(f"❌ Image file does not exist: {sim_path}")
        return
    
    # Run fast segmentation with default parameters
    success = mito_sim_segmentation(
        im_path=sim_path,
        downsample_factor=4  # 4x downsampling, balance speed and accuracy
    )
    
    if success:
        print("\n✅ Fast segmentation completed successfully!")
        print("\n📋 Output files:")
        print("   - mitochondria_mask.tif: Final binary mask (0-1 values, original size) 🎯")
        print("   - middle_slice_original.tif: Original middle slice")
        print("   - middle_slice_mask.tif: Mask middle slice (0-255 for visualization)")
        print("\n🎯 Main result: mitochondria_mask.tif (binary mask with original image size)")
        print("\n💡 For more parameter options, run demo_sim_segmentation.py")
    else:
        print("\n❌ Error occurred during segmentation")
    

if __name__ == "__main__":
    main()
    input("\nPress any key to exit...")