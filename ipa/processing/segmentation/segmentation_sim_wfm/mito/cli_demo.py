#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple CLI demo for SIM Mitochondria Segmentation
Usage: python cli_demo.py --input path/to/image.tif [options]
"""

import argparse
import numpy as np
from pathlib import Path
import tifffile
import time
from scipy import ndimage as ndi
from skimage import transform, filters, morphology, measure
import time
import warnings
warnings.filterwarnings('ignore')




def mito_sim_segmentation(image_data, 
                         output_dir=None,
                         downsample_factor=4,
                         gaussian_sigma1=1.0,
                         gaussian_sigma2=2.0,
                         laplacian_weight=0.5,
                         min_object_size_3d=50,
                         min_object_size_2d=20,
                         morphology_ball_size=1,
                         max_z_layers=None,
                         max_crop_size=1024,
                         verbose=True):
    """
    Fast SIM image mitochondria segmentation with adjustable parameters
    
    Parameters:
    -----------
    image_data : numpy.ndarray
        Input image data (2D or 3D)
    output_dir : str or Path, optional
        Output directory for results
    downsample_factor : int, default=4
        Downsampling factor for processing speed
    gaussian_sigma1 : float, default=1.0
        First Gaussian filter sigma for enhancement
    gaussian_sigma2 : float, default=2.0
        Second Gaussian filter sigma for enhancement
    laplacian_weight : float, default=0.5
        Weight for Laplacian enhancement
    min_object_size_3d : int, default=50
        Minimum object size for 3D morphological cleanup
    min_object_size_2d : int, default=20
        Minimum object size for 2D morphological cleanup
    morphology_ball_size : int, default=1
        Size of morphological ball/disk for opening operation
    max_z_layers : int, optional
        Maximum Z layers to process (None = all layers)
    max_crop_size : int, default=1024
        Maximum crop size for X-Y dimensions
    verbose : bool, default=True
        Print processing information
        
    Returns:
    --------
    dict : Dictionary containing results
        - 'mask': Final binary mask (0-1 values, original size)
        - 'middle_slice_original': Original middle slice
        - 'middle_slice_mask': Mask middle slice
        - 'num_objects': Number of detected objects
        - 'processing_time': Total processing time
    """
    
    if verbose:
        print("🚀 Fast SIM Mitochondria Segmentation")
        print("=" * 40)
    
    start_time = time.time()
    
    # Setup output directory
    if output_dir is None:
        output_dir = Path.cwd() / "mitochondria_results"
    else:
        output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Step 1: Image preprocessing
    if verbose:
        print("\n📂 Step 1: Image Preprocessing")
    
    original_image = image_data.copy()
    
    if verbose:
        print(f"   Original size: {original_image.shape}")
        print(f"   Original data size: {original_image.nbytes / 1024**3:.2f} GB")
    
    # Save original shape and crop information
    original_shape = original_image.shape
    crop_info = None
    
    # Create processing image copy
    image_array = original_image.copy()
    
    # Select central region and downsample
    if len(image_array.shape) == 3:
        # Select Z slices based on max_z_layers parameter
        if max_z_layers is not None and max_z_layers < image_array.shape[0]:
            z_center = image_array.shape[0] // 2
            z_range = min(max_z_layers, image_array.shape[0])
            z_start = max(0, z_center - z_range // 2)
            z_end = min(image_array.shape[0], z_start + z_range)
            image_array = image_array[z_start:z_end]
            if verbose:
                print(f"   Selected Z layers: {z_start}-{z_end} (center region)")
        else:
            # Use all Z layers
            z_start = 0
            z_end = image_array.shape[0]
            if verbose:
                print(f"   Using all Z layers: 0-{z_end}")
        
        # Record Z cropping info
        z_crop_info = (z_start, z_end, original_shape[0])
        
        # Crop central region
        h, w = image_array.shape[1:]
        if h > max_crop_size or w > max_crop_size:
            center_h, center_w = h // 2, w // 2
            crop_size = max_crop_size
            start_h = max(0, center_h - crop_size // 2)
            end_h = min(h, start_h + crop_size)
            start_w = max(0, center_w - crop_size // 2)
            end_w = min(w, start_w + crop_size)
            
            image_array = image_array[:, start_h:end_h, start_w:end_w]
            crop_info = (start_h, end_h, start_w, end_w, h, w)
            if verbose:
                print(f"   Crop region: H[{start_h}:{end_h}], W[{start_w}:{end_w}]")
        else:
            crop_info = None
            z_crop_info = (z_start, z_end, original_shape[0])
        
        if verbose:
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
            gaussian1 = ndi.gaussian_filter(slice_2d.astype('float32'), sigma=gaussian_sigma1)
            gaussian2 = ndi.gaussian_filter(slice_2d.astype('float32'), sigma=gaussian_sigma2)
            
            # Calculate Laplacian response
            laplacian = ndi.laplace(gaussian1)
            
            # Combine different scales
            enhanced_slice = gaussian1 - gaussian2 + np.abs(laplacian) * laplacian_weight
            enhanced_slice[enhanced_slice < 0] = 0
            
            enhanced[i] = enhanced_slice
    else:
        # Process 2D
        gaussian1 = ndi.gaussian_filter(image_array.astype('float32'), sigma=gaussian_sigma1)
        gaussian2 = ndi.gaussian_filter(image_array.astype('float32'), sigma=gaussian_sigma2)
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
        binary = morphology.binary_opening(binary, morphology.disk(morphology_ball_size))
    
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
    
    if verbose:
        print(f"\n🎉 Fast segmentation completed!")
        print(f"📁 Results saved in: {output_dir}")
        print(f"⏱️  Total time: {total_time:.1f} seconds")
        print(f"📊 Detected {final_num_objects} mitochondria structures")
    
    # Return results dictionary
    results = {
        'mask': final_mask,
        'middle_slice_original': original_slice,
        'middle_slice_mask': mask_slice,
        'num_objects': final_num_objects,
        'processing_time': total_time,
        'output_dir': str(output_dir)
    }
    
    return results




def main():
    """Main function with command line argument parsing"""
    
    parser = argparse.ArgumentParser(
        description='SIM Mitochondria Segmentation Demo',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('--input', '-i', 
                       type=str, 
                       required=True,
                       help='Input image file path')
    
    parser.add_argument('--output', '-o',
                       type=str,
                       default='segmentation_results',
                       help='Output directory')
    
    # Processing parameters
    parser.add_argument('--downsample', 
                       type=int, 
                       default=4,
                       help='Downsampling factor (higher = faster, lower quality)')
    
    parser.add_argument('--max_z_layers',
                       type=int,
                       default=30,
                       help='Maximum Z layers to process (0 = all layers)')
    
    parser.add_argument('--max_crop_size',
                       type=int,
                       default=1024,
                       help='Maximum crop size for XY dimensions')
    
    # Enhancement parameters
    parser.add_argument('--sigma1',
                       type=float,
                       default=1.0,
                       help='Fine structure Gaussian sigma')
    
    parser.add_argument('--sigma2',
                       type=float,
                       default=2.0,
                       help='Coarse structure Gaussian sigma')
    
    parser.add_argument('--laplacian_weight',
                       type=float,
                       default=0.5,
                       help='Laplacian enhancement weight')
    
    # Segmentation parameters
    parser.add_argument('--min_size_3d',
                       type=int,
                       default=50,
                       help='Minimum object size for 3D')
    
    parser.add_argument('--min_size_2d',
                       type=int,
                       default=20,
                       help='Minimum object size for 2D')
    
    parser.add_argument('--morph_size',
                       type=int,
                       default=1,
                       help='Morphological operations size')
    
    # Output control
    parser.add_argument('--quiet', '-q',
                       action='store_true',
                       help='Suppress progress output')
    
    args = parser.parse_args()
    
    # Validate input file
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"❌ Error: Input file does not exist: {input_path}")
        return 1
    
    # Load image data
    print(f"📂 Loading image: {input_path}")
    try:
        with tifffile.TiffFile(input_path) as tif:
            image_data = tif.asarray()
        print(f"   Image shape: {image_data.shape}")
        print(f"   Data type: {image_data.dtype}")
        print(f"   Memory usage: {image_data.nbytes / 1024**3:.2f} GB")
    except Exception as e:
        print(f"❌ Error loading image: {e}")
        return 1
    
    # Prepare parameters
    max_z = None if args.max_z_layers == 0 else args.max_z_layers
    
    # Print parameters
    if not args.quiet:
        print(f"\n🔧 Processing parameters:")
        print(f"   Downsample factor: {args.downsample}")
        print(f"   Max Z layers: {'All' if max_z is None else max_z}")
        print(f"   Max crop size: {args.max_crop_size}")
        print(f"   Gaussian sigmas: {args.sigma1}, {args.sigma2}")
        print(f"   Laplacian weight: {args.laplacian_weight}")
        print(f"   Min object sizes: 3D={args.min_size_3d}, 2D={args.min_size_2d}")
        print(f"   Morphology size: {args.morph_size}")
        print(f"   Output directory: {args.output}")
    
    # Run segmentation
    try:
        results = mito_sim_segmentation(
            image_data=image_data,
            output_dir=args.output,
            downsample_factor=args.downsample,
            max_z_layers=max_z,
            max_crop_size=args.max_crop_size,
            gaussian_sigma1=args.sigma1,
            gaussian_sigma2=args.sigma2,
            laplacian_weight=args.laplacian_weight,
            min_object_size_3d=args.min_size_3d,
            min_object_size_2d=args.min_size_2d,
            morphology_ball_size=args.morph_size,
            verbose=not args.quiet
        )
        
        # Print results
        print(f"\n🎉 Segmentation completed successfully!")
        print(f"   Detected objects: {results['num_objects']}")
        print(f"   Processing time: {results['processing_time']:.1f} seconds")
        print(f"   Output mask shape: {results['mask'].shape}")
        print(f"   Results saved in: {results['output_dir']}")
        print(f"\n🎯 Main result: {results['output_dir']}/mitochondria_mask.tif")
        
        return 0
        
    except Exception as e:
        print(f"❌ Error during segmentation: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    exit(main())