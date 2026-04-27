"""
CryoET Multi-Organelle Interaction Analysis Demo

This script provides a production-ready interface for running various cryoET 
organelle interaction analyses with proper logging, error handling, and 
reproducibility features.

Supports 14 different analysis types including:
- Actin-related: actin-to-actin, actin-to-vesicle, actin-to-microtube, actin-to-mito, actin-to-ER
- Microtube-related: microtube-to-actin, microtube-to-vesicle, microtube-to-mito, microtube-to-ER
- Vesicle-related: vesicle-to-mito, vesicle-to-ER
- Organelle interactions: mito-to-ER
- Enhanced analyses: actin angle-distance pair enhanced analysis

Usage:
    # Run specific analyses
    python demo_cryoET_analysis.py --main_path /path/to/data --data_id 20220326_2.8_2 --run_a2a --run_a2v --visualization
    
    # Run ER-related analyses
    python demo_cryoET_analysis.py --main_path /path/to/data --data_id 20220326_2.8_2 --run_v2er --run_mito2er --run_a2er --run_mt2er
    
    # Run all analyses
    python demo_cryoET_analysis.py --main_path /path/to/data --data_id 20220326_2.8_2 --run_a2a --run_a2v --run_a2v_da --run_a2v_angle_dist_enhanced --run_a2mt --run_a2mito --run_mt2a --run_mt2v --run_mt2mito --run_mito2er --run_mt2er --run_a2er --run_v2mito --run_v2er
"""

import os
import sys
import json
import time
import argparse
import matplotlib
matplotlib.use('Agg')

from ipa.data_loader import QuickLogger
from ipa.analysis import (
    actin_to_actin_analysis,
    actin_to_vesicle_analysis,
    actin_to_vesicle_dist_angle_pair_analysis,
    actin_angle_distance_pair_enhanced_analysis,
    actin_to_microtube_analysis,
    actin_to_mito_analysis,
    actin_to_endoreticulum_analysis,
    microtube_to_actin_analysis,
    microtube_to_vesicle_analysis,
    microtube_to_mito_analysis,
    microtube_to_endoreticulum_analysis,
    mito_to_endoreticulum_analysis,
    vesicle_to_mito_analysis,
    vesicle_to_endoreticulum_analysis,
)


def parse_args():
    """Parse command line arguments for cryoET analysis"""
    p = argparse.ArgumentParser(description="cryoET Multi-Organelle Interaction Analysis Demo")
    
    # Required arguments
    p.add_argument("--main_path", type=str, required=True, 
                   help="Project main path (containing data/cryoET/...)")
    p.add_argument("--data_id", type=str, required=True, 
                   help="Data ID (e.g., 20220326_2.8_2)")
    
    # Optional arguments
    p.add_argument("--output_dir", type=str, default=None, 
                   help="Output directory for results (default: {main_path}/outputs/demo_cryoET_analysis)")
    p.add_argument("--seed", type=int, default=42, 
                   help="Random seed for reproducibility")
    
    # Analysis selection flags
    p.add_argument("--run_a2a", action="store_true", 
                   help="Run actin-to-actin analysis")
    p.add_argument("--run_a2v", action="store_true", 
                   help="Run actin-to-vesicle analysis")
    p.add_argument("--run_a2v_da", action="store_true", 
                   help="Run actin-to-vesicle distance-angle pair analysis")
    p.add_argument("--run_a2v_angle_dist_enhanced", action="store_true", 
                   help="Run enhanced actin angle-distance pair analysis")
    p.add_argument("--run_a2mt", action="store_true", 
                   help="Run actin-to-microtube analysis")
    p.add_argument("--run_a2mito", action="store_true", 
                   help="Run actin-to-mitochondria analysis")
    p.add_argument("--run_mt2a", action="store_true", 
                   help="Run microtube-to-actin analysis")
    p.add_argument("--run_mt2v", action="store_true", 
                   help="Run microtube-to-vesicle analysis")
    p.add_argument("--run_mt2mito", action="store_true", 
                   help="Run microtube-to-mitochondria analysis")
    p.add_argument("--run_mito2er", action="store_true", 
                   help="Run mitochondria-to-endoplasmic reticulum analysis")
    p.add_argument("--run_mt2er", action="store_true", 
                   help="Run microtube-to-endoplasmic reticulum analysis")
    p.add_argument("--run_a2er", action="store_true", 
                   help="Run actin-to-endoplasmic reticulum analysis")
    p.add_argument("--run_v2mito", action="store_true", 
                   help="Run vesicle-to-mitochondria analysis")
    p.add_argument("--run_v2er", action="store_true", 
                   help="Run vesicle-to-endoplasmic reticulum analysis")
    
    # Visualization and performance parameters
    p.add_argument("--visualization", action="store_true", 
                   help="Enable visualization")
    p.add_argument("--batch_size", type=int, default=512, 
                   help="Batch size for processing")
    p.add_argument("--zoom_degree", type=float, default=0.1, 
                   help="Zoom degree for downsampling")
    
    # Voxel size and shift bias
    p.add_argument("--voxel_size", type=float, nargs=3, 
                   default=[13.412, 13.412, 13.412], 
                   help="Voxel size [x, y, z] in nanometers")
    p.add_argument("--shift_bias", type=float, nargs=3, 
                   default=[3046.30, 3214.95, -6171.17], 
                   help="Shift bias [x, y, z]")
    
    return p.parse_args()


def set_seed(seed: int):
    """Set random seeds for reproducibility"""
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


def ensure_outdir(base: str):
    """Create timestamped output directory"""
    os.makedirs(base, exist_ok=True)
    run_id = time.strftime("%Y%m%d-%H%M%S")
    outdir = os.path.join(base, run_id)
    os.makedirs(outdir, exist_ok=True)
    return outdir


def build_paths(main_path: str, data_id: str):
    """Build file paths for all data files"""
    datapath = os.path.join(main_path, "data", "cryoET", data_id)
    
    paths = {
        "mask_file": os.path.join(datapath, f"{data_id}_vesicle.mrc"),
        "json_file": os.path.join(datapath, f"{data_id}_actin_filled_points.json"),
        "angle_file": os.path.join(datapath, f"{data_id}_actin_to_pm_angles.json"),
        "microtube_file": os.path.join(datapath, f"{data_id}_microtube_filled_points.json"),
        "mito_file": os.path.join(datapath, f"{data_id}_mito.mrc"),
        "er_file": os.path.join(datapath, f"{data_id}_endoreticulum.mrc"),
        "datapath": datapath,
        "data_id": data_id,
    }
    
    return paths


def check_exists(logger: QuickLogger, **kwargs):
    """Check if required files exist and log warnings for missing files"""
    for k, v in kwargs.items():
        if v and not os.path.exists(v):
            logger.error(f"Missing file for {k}: {v}")


def save_json(outdir: str, name: str, obj: dict):
    """Save analysis results to JSON file"""
    fp = os.path.join(outdir, f"{name}.json")
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)
    return fp


def main():
    """Main entry point for cryoET analysis demo"""
    args = parse_args()
    set_seed(args.seed)
    
    # Set default output directory under main_path if not specified
    if args.output_dir is None:
        args.output_dir = os.path.join(args.main_path, "outputs", "demo_cryoET_analysis")
    
    outdir = ensure_outdir(args.output_dir)
    
    # Initialize logger
    logger = QuickLogger(name="cryoET_demo", log_dir=outdir)
    
    # Build file paths
    paths = build_paths(args.main_path, args.data_id)
    
    logger.step(f"Main path: {args.main_path}")
    logger.step(f"Data ID: {args.data_id}")
    logger.step(f"Output directory: {outdir}")
    
    # Build configuration
    config = {
        "voxel_size_xyz": args.voxel_size,
        "shift_bias": args.shift_bias,
        "visualization": args.visualization,
        "save_intermediate": True,
        "output_dir": outdir,
        "zoom_degree": args.zoom_degree,
        "batch_size": args.batch_size,
    }
    
    # Save configuration and paths
    config_data = {
        "config": config,
        "paths": paths,
        "seed": args.seed,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
    }
    save_json(outdir, "config", config_data)
    logger.step(f"Configuration saved to {outdir}/config.json")
    
    # Check file existence
    check_exists(
        logger,
        mask_file=paths["mask_file"],
        json_file=paths["json_file"],
        angle_file=paths["angle_file"],
        microtube_file=paths["microtube_file"],
        mito_file=paths["mito_file"],
        er_file=paths["er_file"],
    )
    
    def timed_run(name, fn, **kwargs):
        """Run analysis function with timing and error handling"""
        t0 = time.time()
        try:
            logger.step(f"Starting {name}...")
            res = fn(**kwargs)
            dt = time.time() - t0
            fp = save_json(outdir, name, res)
            logger.step(f"{name}: completed in {dt:.2f}s -> {fp}")
        except FileNotFoundError as e:
            logger.error(f"{name}: file not found: {e}")
        except MemoryError:
            logger.error(f"{name}: out of memory. Try reducing resolution, cropping data, or decreasing batch_size.")
        except Exception as e:
            logger.error(f"{name}: {type(e).__name__}: {e}")
    
    # Run selected analyses
    if args.run_a2a and os.path.exists(paths["json_file"]) and os.path.exists(paths["mask_file"]):
        timed_run("actin_to_actin_analysis", actin_to_actin_analysis,
                  data_id=paths["data_id"], mask_file=paths["mask_file"],
                  actin_file=paths["json_file"], config=config)
    
    if args.run_a2v and os.path.exists(paths["json_file"]) and os.path.exists(paths["mask_file"]):
        timed_run("actin_to_vesicle_analysis", actin_to_vesicle_analysis,
                  data_id=paths["data_id"], mask_file=paths["mask_file"],
                  json_file=paths["json_file"], config=config)
    
    if args.run_a2v_da and all(os.path.exists(p) for p in [paths["json_file"], paths["mask_file"], paths["angle_file"]]):
        timed_run("actin_to_vesicle_dist_angle_pair_analysis", actin_to_vesicle_dist_angle_pair_analysis,
                  data_id=paths["data_id"], mask_file=paths["mask_file"],
                  json_file=paths["json_file"], angle_file=paths["angle_file"], config=config)
    
    if args.run_a2v_angle_dist_enhanced and all(os.path.exists(p) for p in [paths["json_file"], paths["mask_file"]]):
        timed_run("actin_angle_distance_pair_enhanced_analysis", actin_angle_distance_pair_enhanced_analysis,
                  data_id=paths["data_id"], actin_file=paths["json_file"],
                  vesicle_file=paths["mask_file"], config=config)
    
    if args.run_a2mt and all(os.path.exists(p) for p in [paths["json_file"], paths["microtube_file"]]):
        timed_run("actin_to_microtube_analysis", actin_to_microtube_analysis,
                  data_id=paths["data_id"], actin_file=paths["json_file"],
                  microtube_file=paths["microtube_file"], config=config)
    
    if args.run_a2mito and all(os.path.exists(p) for p in [paths["json_file"], paths["mito_file"]]):
        timed_run("actin_to_mito_analysis", actin_to_mito_analysis,
                  data_id=paths["data_id"], actin_file=paths["json_file"],
                  mito_file=paths["mito_file"], config=config)
    
    if args.run_mt2a and all(os.path.exists(p) for p in [paths["json_file"], paths["microtube_file"]]):
        timed_run("microtube_to_actin_analysis", microtube_to_actin_analysis,
                  data_id=paths["data_id"], actin_file=paths["json_file"],
                  microtube_file=paths["microtube_file"], config=config)
    
    if args.run_mt2v and all(os.path.exists(p) for p in [paths["microtube_file"], paths["mask_file"]]):
        timed_run("microtube_to_vesicle_analysis", microtube_to_vesicle_analysis,
                  data_id=paths["data_id"], microtube_file=paths["microtube_file"],
                  vesicle_file=paths["mask_file"], config=config)
    
    if args.run_mt2mito and all(os.path.exists(p) for p in [paths["microtube_file"], paths["mito_file"]]):
        timed_run("microtube_to_mito_analysis", microtube_to_mito_analysis,
                  data_id=paths["data_id"], microtube_file=paths["microtube_file"],
                  mito_file=paths["mito_file"], config=config)
    
    if args.run_mito2er and all(os.path.exists(p) for p in [paths["mito_file"], paths["er_file"]]):
        timed_run("mito_to_endoreticulum_analysis", mito_to_endoreticulum_analysis,
                  data_id=paths["data_id"], mito_file=paths["mito_file"],
                  er_file=paths["er_file"], config=config)
    
    if args.run_mt2er and all(os.path.exists(p) for p in [paths["microtube_file"], paths["er_file"]]):
        timed_run("microtube_to_endoreticulum_analysis", microtube_to_endoreticulum_analysis,
                  data_id=paths["data_id"], microtube_file=paths["microtube_file"],
                  er_file=paths["er_file"], config=config)
    
    if args.run_a2er and all(os.path.exists(p) for p in [paths["json_file"], paths["er_file"]]):
        timed_run("actin_to_endoreticulum_analysis", actin_to_endoreticulum_analysis,
                  data_id=paths["data_id"], actin_file=paths["json_file"],
                  er_file=paths["er_file"], config=config)
    
    if args.run_v2mito and all(os.path.exists(p) for p in [paths["mask_file"], paths["mito_file"]]):
        timed_run("vesicle_to_mito_analysis", vesicle_to_mito_analysis,
                  data_id=paths["data_id"], vesicle_file=paths["mask_file"],
                  mito_file=paths["mito_file"], config=config)
    
    if args.run_v2er and all(os.path.exists(p) for p in [paths["mask_file"], paths["er_file"]]):
        timed_run("vesicle_to_endoreticulum_analysis", vesicle_to_endoreticulum_analysis,
                  data_id=paths["data_id"], vesicle_file=paths["mask_file"],
                  er_file=paths["er_file"], config=config)
    
    logger.step(f"All analyses completed. Results saved to: {outdir}")
    print(f"\nOutputs saved to: {outdir}")


if __name__ == "__main__":
    main()
