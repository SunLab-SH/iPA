"""
SIM Organelle Segmentation - Multi-Organelle Segmentation Suite

This script provides a production-ready interface for segmenting various
organelles from SIM imaging data including actin filaments, ISG spheres,
cell shapes, and nuclei.

Supported operations:
- Actin skeletonization: 3D skeleton extraction from actin filaments
- ISG segmentation: Sphere-like organelle detection and labeling
- Cell shape segmentation: Cell membrane/boundary segmentation
- Nucleus segmentation: Nuclear envelope detection and labeling

Usage:
    # Run all segmentation types
    python demo_SIM_organelle_segmentation.py --main_path /path/to/ipa --run_all
    
    # Run specific segmentation
    python demo_SIM_organelle_segmentation.py --main_path /path/to/ipa --run_actin --run_nucleus
    
    # Run with visualization
    python demo_SIM_organelle_segmentation.py --main_path /path/to/ipa --run_all --visualization
"""

# ═══════════════════════════════════════════════════════════════════
# Imports
# ═══════════════════════════════════════════════════════════════════
import os
import sys
import json
import time
import argparse
import numpy as np
import tifffile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pathlib import Path

from ipa.data_loader import UniversalDataLoader, QuickLogger
from ipa.processing.segmentation import (
    skeletonize_organelle, 
    segment_sphere_like_organelle, 
    segment_cell_shape, 
    segment_nucleus
)


# ═══════════════════════════════════════════════════════════════════
# Argument Parsing
# ═══════════════════════════════════════════════════════════════════
def parse_args():
    """Parse command line arguments for SIM organelle segmentation"""
    p = argparse.ArgumentParser(
        description="SIM Organelle Segmentation Suite",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    p.add_argument("--main_path", type=str, required=True,
                   help="Project main path (containing data/ directory)")
    
    # Optional arguments
    p.add_argument("--output_dir", type=str, default=None,
                   help="Output directory (default: {main_path}/outputs/organelle_segmentation)")
    p.add_argument("--seed", type=int, default=42,
                   help="Random seed for reproducibility")
    
    # Analysis selection flags
    p.add_argument("--run_actin", action="store_true",
                   help="Run actin filament skeletonization")
    p.add_argument("--run_isg", action="store_true",
                   help="Run ISG sphere segmentation")
    p.add_argument("--run_cell_shape", action="store_true",
                   help="Run cell shape segmentation")
    p.add_argument("--run_nucleus", action="store_true",
                   help="Run nucleus segmentation")
    p.add_argument("--run_all", action="store_true",
                   help="Run all segmentation types")
    
    # Visualization parameters
    p.add_argument("--visualization", action="store_true",
                   help="Enable visualization output")
    p.add_argument("--auto_save", action="store_true", default=True,
                   help="Auto-save results without prompting")
    
    # Domain-specific parameters
    p.add_argument("--actin_file", type=str, default="20220909_30-2-1-SIM_raw_Actin.tif",
                   help="Actin image filename")
    p.add_argument("--isg_file", type=str, default="20220909_30-2-1-SIM_raw_ISG.tif",
                   help="ISG image filename")
    p.add_argument("--nucleus_file", type=str, default="20220909_30-2-1-SIM_raw_N.tif",
                   help="Nucleus image filename")
    
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


def build_paths(main_path: str, actin_file: str, isg_file: str, nucleus_file: str) -> dict:
    """Build file paths for organelle segmentation"""
    img_dir = os.path.join(main_path, "data", "sim_images")
    
    paths = {
        "actin_path": os.path.join(img_dir, actin_file),
        "isg_path": os.path.join(img_dir, isg_file),
        "nucleus_path": os.path.join(img_dir, nucleus_file),
        "img_dir": img_dir
    }
    
    return paths


def check_exists(logger: QuickLogger, **kwargs):
    """Check if required files exist"""
    for k, v in kwargs.items():
        if v and not os.path.exists(v):
            logger.step(f"Warning: Missing {k}: {v}")


def save_results(result_data, filename, outdir, logger):
    """Save segmentation results to TIFF file"""
    output_path = os.path.join(outdir, f"{filename}_result.tif")
    
    if result_data.dtype == bool:
        tifffile.imwrite(output_path, result_data.astype(np.uint8) * 255)
    else:
        tifffile.imwrite(output_path, result_data.astype(np.uint16))
    
    logger.step(f"Saved: {output_path}")
    return output_path


def visualize_and_save(original, result, title, output_path, is_skeleton=False, is_labeled=False, num_objects=0):
    """Create and save visualization"""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(title, fontsize=16)
    
    mid_z = original.shape[0] // 2
    mid_y = original.shape[1] // 2
    mid_x = original.shape[2] // 2
    
    # Original views
    axes[0, 0].imshow(original[mid_z], cmap='gray')
    axes[0, 0].set_title(f'Original - Z slice {mid_z}')
    axes[0, 0].axis('off')
    
    axes[0, 1].imshow(original[:, mid_y], cmap='gray')
    axes[0, 1].set_title('Original - Y projection')
    axes[0, 1].axis('off')
    
    axes[0, 2].imshow(original[:, :, mid_x], cmap='gray')
    axes[0, 2].set_title('Original - X projection')
    axes[0, 2].axis('off')
    
    # Result views
    cmap = 'hot' if is_skeleton else ('tab20' if is_labeled else 'Reds')
    axes[1, 0].imshow(result[mid_z], cmap=cmap)
    axes[1, 0].set_title(f'Result - Z slice {mid_z}')
    axes[1, 0].axis('off')
    
    axes[1, 1].imshow(result[:, mid_y], cmap=cmap)
    axes[1, 1].set_title('Result - Y projection')
    axes[1, 1].axis('off')
    
    axes[1, 2].imshow(result[:, :, mid_x], cmap=cmap)
    axes[1, 2].set_title('Result - X projection')
    axes[1, 2].axis('off')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


# ═══════════════════════════════════════════════════════════════════
# Segmentation Functions
# ═══════════════════════════════════════════════════════════════════

def run_actin_segmentation(paths, outdir, visualization, logger):
    """Run actin filament skeletonization"""
    if not os.path.exists(paths["actin_path"]):
        logger.step(f"Skip actin - file not found: {paths['actin_path']}")
        return None
    
    image_data = UniversalDataLoader.load_data(paths["actin_path"])
    logger.step(f"Loaded actin data: shape={image_data.shape}")
    
    skeleton_result = skeletonize_organelle(image_data)
    
    result_path = save_results(skeleton_result, "actin_skeleton", outdir, logger)
    
    if visualization:
        viz_path = os.path.join(outdir, "actin_skeleton_viz.png")
        visualize_and_save(image_data, skeleton_result, "Actin Filament Skeletonization", 
                          viz_path, is_skeleton=True)
        logger.step(f"Visualization saved: {viz_path}")
    
    return {'result_path': result_path, 'shape': skeleton_result.shape}


def run_isg_segmentation(paths, outdir, visualization, logger):
    """Run ISG sphere segmentation"""
    if not os.path.exists(paths["isg_path"]):
        logger.step(f"Skip ISG - file not found: {paths['isg_path']}")
        return None
    
    image_data = UniversalDataLoader.load_data(paths["isg_path"])
    logger.step(f"Loaded ISG data: shape={image_data.shape}")
    
    labeled_spheres, num_spheres = segment_sphere_like_organelle(
        image_data, 
        threshold=None,
        min_size=100
    )
    
    logger.step(f"Detected {num_spheres} ISG spheres")
    
    result_path = save_results(labeled_spheres, "isg_spheres", outdir, logger)
    
    if visualization:
        viz_path = os.path.join(outdir, "isg_spheres_viz.png")
        visualize_and_save(image_data, labeled_spheres, f"ISG Spheres ({num_spheres} detected)", 
                          viz_path, is_labeled=True, num_objects=num_spheres)
        logger.step(f"Visualization saved: {viz_path}")
    
    return {'result_path': result_path, 'num_spheres': num_spheres, 'shape': labeled_spheres.shape}


def run_cell_shape_segmentation(paths, outdir, visualization, logger):
    """Run cell shape segmentation"""
    if not os.path.exists(paths["actin_path"]):
        logger.step(f"Skip cell shape - file not found: {paths['actin_path']}")
        return None
    
    image_data = UniversalDataLoader.load_data(paths["actin_path"])
    logger.step(f"Loaded cell shape data: shape={image_data.shape}")
    
    labeled_cell = segment_cell_shape(
        image_data,
        threshold=None,
        min_size=1000
    )
    
    result_path = save_results(labeled_cell, "cell_shape", outdir, logger)
    
    if visualization:
        viz_path = os.path.join(outdir, "cell_shape_viz.png")
        visualize_and_save(image_data, labeled_cell > 0, "Cell Shape Segmentation", 
                          viz_path, is_labeled=False)
        logger.step(f"Visualization saved: {viz_path}")
    
    return {'result_path': result_path, 'shape': labeled_cell.shape}


def run_nucleus_segmentation(paths, outdir, visualization, logger):
    """Run nucleus segmentation"""
    if not os.path.exists(paths["nucleus_path"]):
        logger.step(f"Skip nucleus - file not found: {paths['nucleus_path']}")
        return None
    
    image_data = UniversalDataLoader.load_data(paths["nucleus_path"])
    logger.step(f"Loaded nucleus data: shape={image_data.shape}")
    
    labeled_nuclei = segment_nucleus(image_data)
    num_nuclei = np.max(labeled_nuclei)
    
    logger.step(f"Detected {num_nuclei} nuclei")
    
    result_path = save_results(labeled_nuclei, "nuclei", outdir, logger)
    
    if visualization:
        viz_path = os.path.join(outdir, "nuclei_viz.png")
        visualize_and_save(image_data, labeled_nuclei, f"Nucleus Segmentation ({num_nuclei} detected)", 
                          viz_path, is_labeled=True, num_objects=num_nuclei)
        logger.step(f"Visualization saved: {viz_path}")
    
    return {'result_path': result_path, 'num_nuclei': num_nuclei, 'shape': labeled_nuclei.shape}


# ═══════════════════════════════════════════════════════════════════
# Main Function
# ═══════════════════════════════════════════════════════════════════

def main():
    """Main entry point"""
    
    # Step 1: Initialize and parse arguments
    args = parse_args()
    set_seed(args.seed)
    
    if args.output_dir is None:
        args.output_dir = os.path.join(args.main_path, "outputs", "organelle_segmentation")
    
    outdir = ensure_outdir(args.output_dir)
    
    # Step 2: Initialize logger
    logger = QuickLogger(name="sim_organelle_segmentation", log_dir=outdir)
    logger.step("SIM Organelle Segmentation Suite")
    
    # Step 3: Build paths and config
    paths = build_paths(args.main_path, args.actin_file, args.isg_file, args.nucleus_file)
    
    logger.step(f"Main path: {args.main_path}")
    logger.step(f"Output: {outdir}")
    
    # Handle run_all flag
    if args.run_all:
        args.run_actin = True
        args.run_isg = True
        args.run_cell_shape = True
        args.run_nucleus = True
    
    config = {
        'run_actin': args.run_actin,
        'run_isg': args.run_isg,
        'run_cell_shape': args.run_cell_shape,
        'run_nucleus': args.run_nucleus,
        'visualization': args.visualization,
        'auto_save': args.auto_save
    }
    
    # Step 4: Save configuration
    config_data = {
        'config': config,
        'paths': paths,
        'seed': args.seed,
        'timestamp': time.strftime("%Y-%m-%d %H:%M:%S")
    }
    save_json(outdir, "config", config_data)
    
    # Step 5: Validate input files (optional - checked in each run function)
    check_exists(logger, 
                actin=paths["actin_path"] if args.run_actin or args.run_cell_shape else None,
                isg=paths["isg_path"] if args.run_isg else None,
                nucleus=paths["nucleus_path"] if args.run_nucleus else None)
    
    # Step 6: Define execution wrapper
    def timed_run(name: str, fn, **kwargs):
        t0 = time.time()
        try:
            logger.step(f"Starting {name}...")
            result = fn(**kwargs)
            dt = time.time() - t0
            if result is not None:
                logger.step(f"{name} completed in {dt:.2f}s")
            return result
        except FileNotFoundError as e:
            logger.step(f"{name} failed: file not found: {e}")
            return None
        except MemoryError:
            logger.step(f"{name} failed: out of memory")
            return None
        except Exception as e:
            logger.step(f"{name} failed: {type(e).__name__}: {e}")
            return None
    
    # Step 7: Conditional execution
    results = {}
    
    if args.run_actin:
        result = timed_run(
            "actin_skeletonization",
            run_actin_segmentation,
            paths=paths,
            outdir=outdir,
            visualization=args.visualization,
            logger=logger
        )
        if result:
            results['actin'] = result
    
    if args.run_isg:
        result = timed_run(
            "isg_segmentation",
            run_isg_segmentation,
            paths=paths,
            outdir=outdir,
            visualization=args.visualization,
            logger=logger
        )
        if result:
            results['isg'] = result
    
    if args.run_cell_shape:
        result = timed_run(
            "cell_shape_segmentation",
            run_cell_shape_segmentation,
            paths=paths,
            outdir=outdir,
            visualization=args.visualization,
            logger=logger
        )
        if result:
            results['cell_shape'] = result
    
    if args.run_nucleus:
        result = timed_run(
            "nucleus_segmentation",
            run_nucleus_segmentation,
            paths=paths,
            outdir=outdir,
            visualization=args.visualization,
            logger=logger
        )
        if result:
            results['nucleus'] = result
    
    # Save summary
    if results:
        save_json(outdir, "segmentation_results", results)
    
    # Step 8: Completion
    logger.step(f"Completed {len(results)} segmentation(s). Results: {outdir}")
    print(f"\nOutputs saved to: {outdir}")


# ═══════════════════════════════════════════════════════════════════
# Entry Point
# ═══════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    main()
