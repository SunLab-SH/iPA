#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
SIM Cell Segmentation - Inference for Cell Membrane and Nucleus

This script provides a production-ready interface for performing cell segmentation
on SIM imaging data using trained deep learning models. It segments cell membrane
and nucleus from SIM images.

Supported operations:
- Cell segmentation: Segment PM (plasma membrane) and NE (nuclear envelope)
- Batch processing: Process multiple datasets
- Result saving: Save segmentation masks

Usage:
    # Run segmentation on single dataset
    python demo_SIM_cell_segmentation.py --main_path /path/to/ipa --data_id 784_5 --run_segmentation
    
    # Run with custom model and output directory
    python demo_SIM_cell_segmentation.py --main_path /path/to/ipa --data_id 784_5 --run_segmentation --output_dir /custom/path --pool_processes 8
"""

# ═══════════════════════════════════════════════════════════════════
# Imports
# ═══════════════════════════════════════════════════════════════════
import os
import sys
import json
import time
import argparse
import warnings
import matplotlib
matplotlib.use('Agg')

# Filter MRC file warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, module="mrcfile")
warnings.filterwarnings("ignore", message=".*Unrecognised machine stamp.*")
warnings.filterwarnings("ignore", message=".*Map ID string not found.*")
warnings.filterwarnings("ignore", message=".*data block cannot be read.*")

from ipa.processing.segmentation import run_cell_segmentation
from ipa.data_loader import UniversalDataLoader, QuickLogger


# ═══════════════════════════════════════════════════════════════════
# Argument Parsing
# ═══════════════════════════════════════════════════════════════════
def parse_args():
    """Parse command line arguments for SIM cell segmentation"""
    p = argparse.ArgumentParser(
        description="SIM Cell Segmentation Inference",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    p.add_argument("--main_path", type=str, required=True,
                   help="Project main path (containing data/ directory)")
    
    # Optional arguments
    p.add_argument("--data_id", type=str, nargs='+', default=['784_5'],
                   help="One or more data IDs to process")
    p.add_argument("--output_dir", type=str, default=None,
                   help="Output directory (default: {main_path}/outputs/cell_segmentation)")
    p.add_argument("--seed", type=int, default=42,
                   help="Random seed for reproducibility")
    
    # Analysis selection flags
    p.add_argument("--run_segmentation", action="store_true",
                   help="Run cell segmentation")
    
    # Performance parameters
    p.add_argument("--pool_processes", type=int, default=6,
                   help="Number of parallel processes")
    
    # Domain-specific parameters
    p.add_argument("--image_dir", type=str, default=None,
                   help="Custom image directory (default: {main_path}/data/sxt)")
    
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


def build_paths(main_path: str, data_id: str, image_dir: str = None) -> dict:
    """Build file paths for segmentation"""
    if image_dir is None:
        image_dir = os.path.join(main_path, "data", "sxt")
    
    paths = {
        "image_path": os.path.join(image_dir, f"Stevens_pancreatic_INS_1E_{data_id}_pre_rec.mrc"),
        "data_id": data_id,
        "image_dir": image_dir
    }
    
    return paths


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
    
    if args.output_dir is None:
        args.output_dir = os.path.join(args.main_path, "outputs", "cell_segmentation")
    
    outdir = ensure_outdir(args.output_dir)
    
    # Step 2: Initialize logger
    logger = QuickLogger(name="sim_cell_segmentation", log_dir=outdir)
    logger.step("SIM Cell Segmentation")
    
    # Step 3: Build paths and config
    logger.step(f"Main path: {args.main_path}")
    logger.step(f"Output: {outdir}")
    logger.step(f"Processing {len(args.data_id)} dataset(s): {args.data_id}")
    
    config = {
        'pool_processes': args.pool_processes,
        'output_dir': outdir,
        'image_dir': args.image_dir if args.image_dir else os.path.join(args.main_path, "data", "sxt")
    }
    
    # Step 4: Save configuration
    config_data = {
        'config': config,
        'data_ids': args.data_id,
        'seed': args.seed,
        'timestamp': time.strftime("%Y-%m-%d %H:%M:%S")
    }
    save_json(outdir, "config", config_data)
    
    # Step 5: Load and validate input files
    image_data = {}
    valid_data_ids = []
    
    for data_id in args.data_id:
        paths = build_paths(args.main_path, data_id, args.image_dir)
        
        if not os.path.exists(paths["image_path"]):
            logger.step(f"Skip {data_id} - image not found: {paths['image_path']}")
            continue
        
        logger.step(f"Loading {data_id}")
        data = UniversalDataLoader.load_data(paths["image_path"])
        image_data[data_id] = data
        valid_data_ids.append(data_id)
        logger.step(f"  Shape: {data.shape}, dtype: {data.dtype}")
    
    if not image_data:
        logger.step("Error: No valid image data loaded")
        print(f"\nNo valid data found. Check paths in: {outdir}")
        return
    
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
            logger.step(f"{name} failed: out of memory. Try reducing pool_processes")
            return None
        except Exception as e:
            logger.step(f"{name} failed: {type(e).__name__}: {e}")
            return None
    
    # Step 7: Conditional execution
    if args.run_segmentation and len(image_data) > 0:
        result = timed_run(
            "cell_segmentation",
            run_cell_segmentation,
            save_dir=outdir,
            pool_processes=config['pool_processes'],
            dataid=valid_data_ids,
            image_data=image_data
        )
        
        if result is not None:
            save_json(outdir, "segmentation_results", {
                'processed_ids': valid_data_ids,
                'status': 'completed'
            })
    
    # Step 8: Completion
    logger.step(f"All operations completed. Results: {outdir}")
    print(f"\nOutputs saved to: {outdir}")


# ═══════════════════════════════════════════════════════════════════
# Entry Point
# ═══════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    main()
