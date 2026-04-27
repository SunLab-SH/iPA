"""
CryoET Actin Filament Segmentation and Skeletonization

This script provides a production-ready interface for performing actin filament 
skeletonization on cryoET denoised tomograms with proper logging, error handling, 
and reproducibility features.

Supported Analysis:
- Actin Filament Skeletonization: 3D skeleton extraction from denoised cryoET data
- Filament Branch Detection: Extract and save skeleton branch coordinates

Usage:
    # Basic usage
    python demo_cryoET_segmentation.py --main_path /path/to/data --data_id 20210517_30_009
    
"""

# ═══════════════════════════════════════════════════════════════════
# Import Modules
# ═══════════════════════════════════════════════════════════════════
import os
import sys
import json
import time
import argparse
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')  # Headless mode for server compatibility
import matplotlib.pyplot as plt

from ipa.data_loader import UniversalDataLoader, QuickLogger
from ipa.processing.segmentation import (
    skeletonization_et_segmentation, 
    save_filament_branches_json
)
import tifffile


# ═══════════════════════════════════════════════════════════════════
# Argument Parsing
# ═══════════════════════════════════════════════════════════════════

def parse_args():
    """Parse command line arguments for cryoET segmentation"""
    p = argparse.ArgumentParser(
        description="CryoET Actin Filament Segmentation and Skeletonization",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # ─────────────────────────────────────────────────────────────
    # Required arguments
    # ─────────────────────────────────────────────────────────────
    p.add_argument(
        "--main_path", 
        type=str, 
        required=True,
        help="Project main path (containing data/cryoET/)"
    )
    p.add_argument(
        "--data_id", 
        type=str, 
        required=True,
        help="Data ID (e.g., 20210517_30_009)"
    )
    
    # ─────────────────────────────────────────────────────────────
    # Optional arguments
    # ─────────────────────────────────────────────────────────────
    p.add_argument(
        "--output_dir", 
        type=str, 
        default=None,
        help="Output directory (default: {main_path}/outputs/demo_cryoET_segmentation)"
    )
    p.add_argument(
        "--seed", 
        type=int, 
        default=42,
        help="Random seed for reproducibility"
    )
    
    # ─────────────────────────────────────────────────────────────
    # Segmentation parameters
    # ─────────────────────────────────────────────────────────────
    p.add_argument(
        "--gaussian_sigma", 
        type=float, 
        default=1.5,
        help="Gaussian smoothing sigma for noise reduction"
    )
    p.add_argument(
        "--threshold_multiplier", 
        type=float, 
        default=1.0,
        help="Threshold multiplier for segmentation"
    )
    p.add_argument(
        "--erosion_radius", 
        type=int, 
        default=1,
        help="Erosion radius for morphological operations"
    )
    p.add_argument(
        "--min_object_size", 
        type=int, 
        default=100,
        help="Minimum object size to keep (voxels)"
    )
    p.add_argument(
        "--dilation_radius", 
        type=int, 
        default=1,
        help="Dilation radius for morphological operations"
    )
    p.add_argument(
        "--min_skeleton_size", 
        type=int, 
        default=20,
        help="Minimum skeleton fragment size"
    )
    p.add_argument(
        "--max_connect_distance", 
        type=int, 
        default=8,
        help="Maximum distance to connect skeleton fragments"
    )
    p.add_argument(
        "--final_min_size", 
        type=int, 
        default=10,
        help="Final minimum skeleton size"
    )
    
    # ─────────────────────────────────────────────────────────────
    # Output control
    # ─────────────────────────────────────────────────────────────
    p.add_argument(
        "--visualization", 
        action="store_true",
        help="Enable visualization output"
    )
    p.add_argument(
        "--interpolate_branches", 
        action="store_true",
        help="Save interpolated skeleton branches as JSON"
    )
    p.add_argument(
        "--points_per_branch", 
        type=int, 
        default=100,
        help="Number of points per branch for interpolation"
    )
    
    return p.parse_args()


# ═══════════════════════════════════════════════════════════════════
# Helper Functions
# ═══════════════════════════════════════════════════════════════════

def set_seed(seed: int):
    """
    Set random seeds for reproducibility
    
    Args:
        seed: Random seed value
    """
    import random
    random.seed(seed)
    np.random.seed(seed)
    
    try:
        import torch
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    except ImportError:
        pass


def ensure_outdir(base: str) -> str:
    """
    Create timestamped output directory
    
    Args:
        base: Base output directory path
        
    Returns:
        Full path to timestamped output directory
    """
    os.makedirs(base, exist_ok=True)
    run_id = time.strftime("%Y%m%d-%H%M%S")
    outdir = os.path.join(base, run_id)
    os.makedirs(outdir, exist_ok=True)
    return outdir


def build_paths(main_path: str, data_id: str) -> dict:
    """
    Build file paths for all data files
    
    Args:
        main_path: Project main path
        data_id: Dataset ID
        
    Returns:
        Dictionary containing all file paths
    """
    datapath = os.path.join(main_path, "data", "cryoET")
    
    paths = {
        "denoised_file": os.path.join(datapath, f"{data_id}_denoised.mrc"),
        "vesicle_file": os.path.join(datapath, f"{data_id}_vesicle.mrc"),
        "datapath": datapath,
        "data_id": data_id,
    }
    
    return paths


def check_exists(logger: QuickLogger, **kwargs):
    """
    Check if required files exist and log warnings for missing files
    
    Args:
        logger: QuickLogger instance
        **kwargs: File name and path key-value pairs
    """
    for k, v in kwargs.items():
        if v and not os.path.exists(v):
            logger.error(f"Missing file for {k}: {v}")


def save_json(outdir: str, name: str, obj: dict) -> str:
    """
    Save analysis results to JSON file
    
    Args:
        outdir: Output directory
        name: File name (without extension)
        obj: Dictionary object to save
        
    Returns:
        Full path to saved file
    """
    fp = os.path.join(outdir, f"{name}.json")
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)
    return fp


# ═══════════════════════════════════════════════════════════════════
# Core Analysis Functions
# ═══════════════════════════════════════════════════════════════════

def run_actin_skeletonization(paths: dict, config: dict, logger: QuickLogger) -> dict:
    """
    Perform actin filament skeletonization on cryoET denoised data
    
    Args:
        paths: Dictionary with file paths
        config: Configuration dictionary with segmentation parameters
        logger: QuickLogger instance
        
    Returns:
        Dictionary with skeleton results and statistics
    """
    denoised_file = paths["denoised_file"]
    vesicle_file = paths["vesicle_file"]
    data_id = paths["data_id"]
    
    # Load denoised cryoET data
    logger.file_in(denoised_file)
    logger.step(f"Loading cryoET denoised data from: {denoised_file}")
    image_data = UniversalDataLoader.load_data(denoised_file)
    
    logger.step(f"Image shape: {image_data.shape}")
    logger.step(f"Data type: {image_data.dtype}")
    logger.step(f"Data range: {image_data.min():.4f} to {image_data.max():.4f}")
    
    # Optional: Load vesicle mask for context
    vesicle_mask = None
    if os.path.exists(vesicle_file):
        logger.file_in(vesicle_file)
        vesicle_mask = UniversalDataLoader.load_data(vesicle_file)
        logger.step(f"Vesicle mask loaded: {vesicle_mask.shape}")
    
    # Perform actin filament skeletonization
    logger.step("Performing 3D actin filament skeletonization...")
    skeleton_result = skeletonization_et_segmentation(
        image_data,
        gaussian_sigma=config["gaussian_sigma"],
        threshold_multiplier=config["threshold_multiplier"],
        erosion_radius=config["erosion_radius"],
        min_object_size=config["min_object_size"],
        dilation_radius=config["dilation_radius"],
        min_skeleton_size=config["min_skeleton_size"],
        max_connect_distance=config["max_connect_distance"],
        final_min_size=config["final_min_size"]
    )
    
    logger.step(f"Skeleton result shape: {skeleton_result.shape}")
    total_skeleton_points = int(np.sum(skeleton_result > 0))
    logger.step(f"Skeleton points detected: {total_skeleton_points}")
    
    # Visualize results if enabled
    if config.get("visualization", False):
        logger.step("Visualizing cryoET actin skeleton results")
        visualize_cryoet_skeleton_results(
            image_data, skeleton_result, vesicle_mask, 
            f"CryoET Actin Filaments - {data_id}",
            config["output_dir"]
        )
    
    # Save results
    save_skeleton_results(skeleton_result, paths, config, logger)
    
    # Prepare results dictionary
    skeleton_coords = np.where(skeleton_result > 0)
    results = {
        "data_id": data_id,
        "total_skeleton_points": total_skeleton_points,
        "skeleton_density": float(total_skeleton_points / skeleton_result.size),
        "image_shape": list(image_data.shape),
        "segmentation_parameters": {
            "gaussian_sigma": config["gaussian_sigma"],
            "threshold_multiplier": config["threshold_multiplier"],
            "min_object_size": config["min_object_size"],
            "min_skeleton_size": config["min_skeleton_size"],
            "max_connect_distance": config["max_connect_distance"],
            "final_min_size": config["final_min_size"],
        }
    }
    
    logger.step("CryoET actin skeletonization completed successfully.")
    return results


def visualize_cryoet_skeleton_results(original, skeleton, vesicle_mask, title, output_dir):
    """
    Visualize cryoET skeleton segmentation results with optional vesicle context
    
    Args:
        original: Original image data
        skeleton: Skeleton segmentation result
        vesicle_mask: Optional vesicle mask for context
        title: Plot title
        output_dir: Directory to save visualization
    """
    fig, axes = plt.subplots(3, 3, figsize=(18, 15))
    fig.suptitle(f'{title} - Skeletonization Results', fontsize=16)
    
    mid_z = original.shape[0] // 2
    mid_y = original.shape[1] // 2
    mid_x = original.shape[2] // 2
    
    # Original image views
    axes[0, 0].imshow(original[mid_z], cmap='gray')
    axes[0, 0].set_title(f'Original - Z slice {mid_z}')
    axes[0, 0].axis('off')
    
    axes[0, 1].imshow(original[:, mid_y], cmap='gray')
    axes[0, 1].set_title('Original - Y projection')
    axes[0, 1].axis('off')
    
    axes[0, 2].imshow(original[:, :, mid_x], cmap='gray')
    axes[0, 2].set_title('Original - X projection')
    axes[0, 2].axis('off')
    
    # Skeleton views
    axes[1, 0].imshow(skeleton[mid_z], cmap='hot')
    axes[1, 0].set_title(f'Skeleton - Z slice {mid_z}')
    axes[1, 0].axis('off')
    
    axes[1, 1].imshow(skeleton[:, mid_y], cmap='hot')
    axes[1, 1].set_title('Skeleton - Y projection')
    axes[1, 1].axis('off')
    
    axes[1, 2].imshow(skeleton[:, :, mid_x], cmap='hot')
    axes[1, 2].set_title('Skeleton - X projection')
    axes[1, 2].axis('off')
    
    # Overlay views (skeleton on original)
    # Normalize original data for overlay
    original_norm = (original - original.min()) / (original.max() - original.min())
    
    overlay_z = original_norm[mid_z].copy()
    overlay_z[skeleton[mid_z] > 0] = 1
    axes[2, 0].imshow(overlay_z, cmap='gray')
    axes[2, 0].set_title(f'Overlay - Z slice {mid_z}')
    axes[2, 0].axis('off')
    
    overlay_y = original_norm[:, mid_y].copy()
    overlay_y[skeleton[:, mid_y] > 0] = 1
    axes[2, 1].imshow(overlay_y, cmap='gray')
    axes[2, 1].set_title('Overlay - Y projection')
    axes[2, 1].axis('off')
    
    overlay_x = original_norm[:, :, mid_x].copy()
    overlay_x[skeleton[:, :, mid_x] > 0] = 1
    axes[2, 2].imshow(overlay_x, cmap='gray')
    axes[2, 2].set_title('Overlay - X projection')
    axes[2, 2].axis('off')
    
    # Add vesicle context if available
    if vesicle_mask is not None:
        # Add cyan contours for vesicle boundaries
        from scipy import ndimage
        
        # Get vesicle contours
        vesicle_edges_z = ndimage.binary_dilation(vesicle_mask[mid_z]) ^ vesicle_mask[mid_z]
        vesicle_edges_y = ndimage.binary_dilation(vesicle_mask[:, mid_y]) ^ vesicle_mask[:, mid_y]
        vesicle_edges_x = ndimage.binary_dilation(vesicle_mask[:, :, mid_x]) ^ vesicle_mask[:, :, mid_x]
        
        # Overlay vesicle contours on skeleton views
        axes[1, 0].contour(vesicle_edges_z, colors='cyan', linewidths=1)
        axes[1, 1].contour(vesicle_edges_y, colors='cyan', linewidths=1)
        axes[1, 2].contour(vesicle_edges_x, colors='cyan', linewidths=1)
    
    plt.tight_layout()
    
    # Save figure instead of showing
    vis_output = os.path.join(output_dir, "actin_skeleton_visualization.png")
    plt.savefig(vis_output, dpi=150, bbox_inches='tight')
    plt.close(fig)
    
    return vis_output


def save_skeleton_results(skeleton_result, paths: dict, config: dict, logger: QuickLogger):
    """
    Save cryoET actin segmentation results
    
    Args:
        skeleton_result: Skeleton segmentation array
        paths: Dictionary with file paths
        config: Configuration dictionary
        logger: QuickLogger instance
    """
    output_dir = config["output_dir"]
    data_id = paths["data_id"]
    
    # Save skeleton as TIFF
    skeleton_output_path = os.path.join(output_dir, f"{data_id}_actin_skeleton.tif")
    tifffile.imwrite(skeleton_output_path, skeleton_result.astype(np.uint8) * 255)
    
    logger.file_out(skeleton_output_path)
    logger.step(f"Skeleton TIFF saved to: {skeleton_output_path}")
    
    # Save skeleton coordinates as text file for further analysis
    skeleton_coords = np.where(skeleton_result > 0)
    coords_output_path = os.path.join(output_dir, f"{data_id}_actin_skeleton_coords.txt")
    
    with open(coords_output_path, 'w') as f:
        f.write("# Actin skeleton coordinates (Z, Y, X)\n")
        for i in range(len(skeleton_coords[0])):
            f.write(f"{skeleton_coords[0][i]}\t{skeleton_coords[1][i]}\t{skeleton_coords[2][i]}\n")
    
    logger.file_out(coords_output_path)
    logger.step(f"Skeleton coordinates saved to: {coords_output_path}")
    
    # Print summary statistics
    total_skeleton_points = int(np.sum(skeleton_result > 0))
    logger.step(f"Total skeleton points: {total_skeleton_points}")
    logger.step(f"Skeleton density: {total_skeleton_points / skeleton_result.size:.6f}")
    
    # Save filament branches if requested
    if config.get("interpolate_branches", False):
        logger.step("Saving interpolated filament branches...")
        save_filament_branches_json(
            skeleton_result, 
            data_id, 
            output_dir, 
            logger, 
            interpolate=True, 
            points_per_branch=config["points_per_branch"]
        )


# ═══════════════════════════════════════════════════════════════════
# Main Function
# ═══════════════════════════════════════════════════════════════════

def main():
    """Main entry point for cryoET actin segmentation"""
    
    # ─────────────────────────────────────────────────────────────
    # Step 1: Parse arguments and set seed
    # ─────────────────────────────────────────────────────────────
    args = parse_args()
    set_seed(args.seed)
    
    # Set default output directory
    if args.output_dir is None:
        args.output_dir = os.path.join(
            args.main_path, 
            "outputs", 
            "demo_cryoET_segmentation"
        )
    
    outdir = ensure_outdir(args.output_dir)
    
    # ─────────────────────────────────────────────────────────────
    # Step 2: Initialize logger
    # ─────────────────────────────────────────────────────────────
    logger = QuickLogger(
        name="cryoET_segmentation",
        log_dir=outdir
    )
    
    logger.step("=" * 70)
    logger.step("CryoET Actin Filament Segmentation Pipeline")
    logger.step("=" * 70)
    
    # ─────────────────────────────────────────────────────────────
    # Step 3: Build file paths and configuration
    # ─────────────────────────────────────────────────────────────
    paths = build_paths(args.main_path, args.data_id)
    
    logger.step(f"Main path: {args.main_path}")
    logger.step(f"Data ID: {args.data_id}")
    logger.step(f"Output directory: {outdir}")
    
    # Build configuration
    config = {
        "gaussian_sigma": args.gaussian_sigma,
        "threshold_multiplier": args.threshold_multiplier,
        "erosion_radius": args.erosion_radius,
        "min_object_size": args.min_object_size,
        "dilation_radius": args.dilation_radius,
        "min_skeleton_size": args.min_skeleton_size,
        "max_connect_distance": args.max_connect_distance,
        "final_min_size": args.final_min_size,
        "visualization": args.visualization,
        "interpolate_branches": args.interpolate_branches,
        "points_per_branch": args.points_per_branch,
        "output_dir": outdir,
    }
    
    # ─────────────────────────────────────────────────────────────
    # Step 4: Save configuration
    # ─────────────────────────────────────────────────────────────
    config_data = {
        "config": config,
        "paths": paths,
        "seed": args.seed,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "command_line_args": vars(args)
    }
    save_json(outdir, "config", config_data)
    logger.step(f"Configuration saved to {outdir}/config.json")
    
    # ─────────────────────────────────────────────────────────────
    # Step 5: Check file existence
    # ─────────────────────────────────────────────────────────────
    check_exists(
        logger,
        denoised_file=paths["denoised_file"]
    )
    
    # ─────────────────────────────────────────────────────────────
    # Step 6: Run actin skeletonization
    # ─────────────────────────────────────────────────────────────
    if not os.path.exists(paths["denoised_file"]):
        logger.error(f"Denoised data not found at {paths['denoised_file']}")
        logger.step("Available files in cryoET directory:")
        if os.path.exists(paths["datapath"]):
            for file in os.listdir(paths["datapath"]):
                if args.data_id in file:
                    logger.step(f"  {file}")
        logger.step("CryoET actin segmentation pipeline failed")
        return
    
    # Run skeletonization with timing and error handling
    t0 = time.time()
    try:
        logger.step("Starting actin filament skeletonization...")
        results = run_actin_skeletonization(paths, config, logger)
        dt = time.time() - t0
        
        # Save results
        fp = save_json(outdir, "actin_segmentation_results", results)
        logger.step(f"Results saved in {dt:.2f}s -> {fp}")
        
    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
    except MemoryError:
        logger.error(
            "Out of memory. Try reducing image size, cropping data, "
            "or adjusting segmentation parameters."
        )
    except Exception as e:
        logger.error(f"{type(e).__name__}: {e}")
        import traceback
        logger.error(traceback.format_exc())
    
    # ─────────────────────────────────────────────────────────────
    # Step 7: Complete
    # ─────────────────────────────────────────────────────────────
    logger.step("=" * 70)
    logger.step(f"Pipeline completed. Results saved to: {outdir}")
    logger.step("=" * 70)
    print(f"\nOutputs saved to: {outdir}")


# ═══════════════════════════════════════════════════════════════════
# Program Entry Point
# ═══════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    main()
