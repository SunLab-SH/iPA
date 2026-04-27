#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SIM Mitochondria Segmentation - Vessel Enhancement Based Segmentation

This script provides a production-ready interface for performing mitochondria
segmentation on SIM imaging data using vessel enhancement filters and
morphological operations.

Supported operations:
- Mitochondria segmentation: Enhanced vessel detection with multi-scale filtering
- 3D morphological operations: Noise reduction and object size filtering
- Parameter tuning: Adjustable quality parameters

Usage:
    # Run mitochondria segmentation with default parameters
    python demo_SIM_mito_segmentation.py --main_path /path/to/ipa --input_image /path/to/image.tif --run_segmentation
    
    # Run with high-quality parameters
    python demo_SIM_mito_segmentation.py --main_path /path/to/ipa --input_image /path/to/image.tif --run_segmentation --downsample_factor 2 --z_range_ratio 1.0
    
    # Run with custom output directory
    python demo_SIM_mito_segmentation.py --main_path /path/to/ipa --input_image /path/to/image.tif --run_segmentation --output_dir /custom/output
"""

# ═══════════════════════════════════════════════════════════════════
# Imports
# ═══════════════════════════════════════════════════════════════════
import os
import sys
import json
import time
import argparse
import matplotlib
matplotlib.use('Agg')

from pathlib import Path
from ipa.processing.segmentation import mito_sim_segmentation
from ipa.data_loader import QuickLogger


# ═══════════════════════════════════════════════════════════════════
# Argument Parsing
# ═══════════════════════════════════════════════════════════════════
def parse_args():
    """Parse command line arguments for SIM mitochondria segmentation"""
    p = argparse.ArgumentParser(
        description="SIM Mitochondria Segmentation",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    p.add_argument("--main_path", type=str, required=True,
                   help="Project main path (containing data/ directory)")
    
    # Optional arguments
    p.add_argument("--input_image", type=str, default=None,
                   help="Path to input SIM image (default: {main_path}/data/sim_images/20220909_30-2-1-SIM_raw_Actin.tif)")
    p.add_argument("--output_dir", type=str, default=None,
                   help="Output directory (default: {main_path}/outputs/mito_segmentation)")
    p.add_argument("--seed", type=int, default=42,
                   help="Random seed for reproducibility")
    
    # Analysis selection flags
    p.add_argument("--run_segmentation", action="store_true",
                   help="Run mitochondria segmentation")
    
    # Performance parameters
    p.add_argument("--downsample_factor", type=int, default=2,
                   help="Downsampling factor for processing")
    p.add_argument("--z_range_ratio", type=float, default=1.0,
                   help="Z range ratio (0-1) to process")
    p.add_argument("--crop_size", type=int, default=1536,
                   help="Crop size for processing")
    
    # Domain-specific parameters
    p.add_argument("--sigma1", type=float, default=0.8,
                   help="Sigma1 for vessel enhancement (fine scale)")
    p.add_argument("--sigma2", type=float, default=2.5,
                   help="Sigma2 for vessel enhancement (coarse scale)")
    p.add_argument("--laplacian_weight", type=float, default=0.7,
                   help="Laplacian weight for edge enhancement")
    p.add_argument("--min_object_size_3d", type=int, default=80,
                   help="Minimum 3D object size")
    p.add_argument("--min_object_size_2d", type=int, default=30,
                   help="Minimum 2D object size")
    p.add_argument("--morphology_ball_size", type=int, default=2,
                   help="Morphological ball size")
    p.add_argument("--morphology_disk_size", type=int, default=2,
                   help="Morphological disk size")
    
    return p.parse_args()


# ═══════════════════════════════════════════════════════════════════
# Helper Functions
# ═══════════════════════════════════════════════════════════════════

def set_seed(seed: int):
    """Set random seed for reproducibility"""
    import random
    import numpy as np
    
    random.seed(seed)
    np.random.seed(seed)
    
    try:
        import torch
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    except ImportError:
        pass


def ensure_outdir(base: str) -> str:
    """Create timestamped output directory"""
    os.makedirs(base, exist_ok=True)
    run_id = time.strftime("%Y%m%d-%H%M%S")
    outdir = os.path.join(base, run_id)
    os.makedirs(outdir, exist_ok=True)
    return outdir


def save_json(outdir: str, name: str, obj: dict) -> str:
    """Save dict as JSON file"""
    fp = os.path.join(outdir, f"{name}.json")
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)
    return fp


def check_exists(logger: QuickLogger, **kwargs):
    """Check if required files exist"""
    for k, v in kwargs.items():
        if v and not os.path.exists(v):
            logger.step(f"Warning: Missing {k}: {v}")


# ═══════════════════════════════════════════════════════════════════
# Main Function
# ═══════════════════════════════════════════════════════════════════

def main():
    """Main entry point"""
    
    # Step 1: Initialize and parse arguments
    args = parse_args()
    set_seed(args.seed)
    
    if args.input_image is None:
        args.input_image = os.path.join(args.main_path, "data", "sim_images", "20220909_30-2-1-SIM_raw_Actin.tif")
    
    if args.output_dir is None:
        args.output_dir = os.path.join(args.main_path, "outputs", "mito_segmentation")
    
    outdir = ensure_outdir(args.output_dir)
    
    # Step 2: Initialize logger
    logger = QuickLogger(name="sim_mito_segmentation", log_dir=outdir)
    logger.step("SIM Mitochondria Segmentation")
    
    # Step 3: Build paths and config
    logger.step(f"Main path: {args.main_path}")
    logger.step(f"Output: {outdir}")
    logger.step(f"Input image: {args.input_image}")
    
    config = {
        'downsample_factor': args.downsample_factor,
        'z_range_ratio': args.z_range_ratio,
        'crop_size': args.crop_size,
        'sigma1': args.sigma1,
        'sigma2': args.sigma2,
        'laplacian_weight': args.laplacian_weight,
        'min_object_size_3d': args.min_object_size_3d,
        'min_object_size_2d': args.min_object_size_2d,
        'morphology_ball_size': args.morphology_ball_size,
        'morphology_disk_size': args.morphology_disk_size
    }
    
    logger.step(f"Configuration: downsample={args.downsample_factor}x, z_ratio={args.z_range_ratio}, crop={args.crop_size}")
    
    # Step 4: Save configuration
    config_data = {
        'config': config,
        'input_image': args.input_image,
        'seed': args.seed,
        'timestamp': time.strftime("%Y-%m-%d %H:%M:%S")
    }
    save_json(outdir, "config", config_data)
    
    # Step 5: Validate input files
    check_exists(logger, input_image=args.input_image)
    
    # Step 6: Define execution wrapper
    def timed_run(name: str, fn, **kwargs):
        t0 = time.time()
        try:
            logger.step(f"Starting {name}...")
            result = fn(**kwargs)
            dt = time.time() - t0
            logger.step(f"{name} completed in {dt:.2f}s")
            return result
        except FileNotFoundError as e:
            logger.step(f"{name} failed: file not found: {e}")
            return None
        except MemoryError:
            logger.step(f"{name} failed: out of memory. Try increasing downsample_factor or reducing crop_size")
            return None
        except Exception as e:
            logger.step(f"{name} failed: {type(e).__name__}: {e}")
            return None
    
    # Step 7: Conditional execution
    if args.run_segmentation and os.path.exists(args.input_image):
        result = timed_run(
            "mitochondria_segmentation",
            mito_sim_segmentation,
            im_path=args.input_image,
            output_dir=outdir,
            downsample_factor=config['downsample_factor'],
            z_range_ratio=config['z_range_ratio'],
            crop_size=config['crop_size'],
            sigma1=config['sigma1'],
            sigma2=config['sigma2'],
            laplacian_weight=config['laplacian_weight'],
            min_object_size_3d=config['min_object_size_3d'],
            min_object_size_2d=config['min_object_size_2d'],
            morphology_ball_size=config['morphology_ball_size'],
            morphology_disk_size=config['morphology_disk_size']
        )
        
        if result is not None:
            save_json(outdir, "segmentation_results", {
                'status': 'completed',
                'parameters': config
            })
    
    # Step 8: Completion
    logger.step(f"All operations completed. Results: {outdir}")
    print(f"\nOutputs saved to: {outdir}")


# ═══════════════════════════════════════════════════════════════════
# Entry Point
# ═══════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    main()
