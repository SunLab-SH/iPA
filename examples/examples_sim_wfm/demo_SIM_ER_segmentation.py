#!/usr/bin/env python3
"""
SIM ER Segmentation - Endoplasmic Reticulum Segmentation with ERNet

This script provides a production-ready interface for performing ER segmentation
on SIM imaging data using the ERNet deep learning model. It provides inference
functionality for segmenting ER structures from fluorescence microscopy images.

Supported operations:
- ER segmentation: Segment endoplasmic reticulum structures
- Custom image processing: Support for various image formats
- Model inference: Use pretrained ERNet model

Usage:
    # Run ER segmentation with default paths
    python demo_SIM_ER_segmentation.py --main_path /path/to/ipa --run_segmentation
    
    # Run with custom image and model
    python demo_SIM_ER_segmentation.py --main_path /path/to/ipa --run_segmentation --image_path /custom/image.png --weights_path /custom/model.pth
    
    # Run with CPU (no GPU)
    python demo_SIM_ER_segmentation.py --main_path /path/to/ipa --run_segmentation --use_cpu
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
import matplotlib
matplotlib.use('Agg')

from types import SimpleNamespace
from PIL import Image
import torch
import torch.nn as nn

from ipa.processing.segmentation import ernet_GetModel as GetModel
from ipa.processing.segmentation import toTensor, toPIL
from ipa.data_loader import QuickLogger


# ═══════════════════════════════════════════════════════════════════
# Argument Parsing
# ═══════════════════════════════════════════════════════════════════
def parse_args():
    """Parse command line arguments for SIM ER segmentation"""
    p = argparse.ArgumentParser(
        description="SIM ER Segmentation with ERNet",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    p.add_argument("--main_path", type=str, required=True,
                   help="Project main path (containing ipa/ directory)")
    
    # Optional arguments
    p.add_argument("--output_dir", type=str, default=None,
                   help="Output directory (default: {main_path}/outputs/er_segmentation)")
    p.add_argument("--seed", type=int, default=42,
                   help="Random seed for reproducibility")
    
    # Analysis selection flags
    p.add_argument("--run_segmentation", action="store_true",
                   help="Run ER segmentation")
    
    # Model parameters
    p.add_argument("--weights_path", type=str, default=None,
                   help="Path to model weights (default: ipa/processing/segmentation/.../pretrained_model/final.pth)")
    p.add_argument("--image_path", type=str, default=None,
                   help="Path to input image (default: ERNet testdata)")
    p.add_argument("--image_size", type=int, default=1000,
                   help="Image size for processing")
    p.add_argument("--use_cpu", action="store_true",
                   help="Use CPU instead of GPU")
    
    # Domain-specific parameters
    p.add_argument("--nch_in", type=int, default=1,
                   help="Number of input channels")
    p.add_argument("--nch_out", type=int, default=2,
                   help="Number of output channels")
    p.add_argument("--n_resgroups", type=int, default=5,
                   help="Number of residual groups")
    p.add_argument("--n_resblocks", type=int, default=10,
                   help="Number of residual blocks")
    
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


def create_pred_config(model_path, image_size, nch_in, nch_out, n_resgroups, n_resblocks, use_cpu):
    """Create configuration object for ERNet prediction"""
    opt = SimpleNamespace()
    opt.weights = model_path
    opt.imageSize = image_size
    opt.model = 'rcan'
    opt.nch_in = nch_in
    opt.nch_out = nch_out
    opt.n_resgroups = n_resgroups
    opt.n_resblocks = n_resblocks
    opt.cpu = use_cpu
    opt.undomulti = False
    opt.n_feats = 64
    opt.reduction = 16
    opt.narch = 0
    opt.multigpu = False
    return opt


def remove_dataparallel_wrapper(state_dict):
    """Remove DataParallel wrapper from state dict"""
    from collections import OrderedDict
    new_state_dict = OrderedDict()
    for k, v in state_dict.items():
        name = k[7:] if k.startswith('module.') else k
        new_state_dict[name] = v
    return new_state_dict


def load_model(config, logger):
    """Load the trained ERNet model"""
    if not os.path.exists(config.weights):
        raise FileNotFoundError(f"Model weights not found: {config.weights}")
    
    net = GetModel(config)
    checkpoint = torch.load(config.weights, map_location='cpu' if config.cpu else None)
    
    if config.undomulti:
        checkpoint['state_dict'] = remove_dataparallel_wrapper(checkpoint['state_dict'])
    
    net.load_state_dict(checkpoint['state_dict'])
    net.eval()
    
    if config.cpu:
        net = net.cpu()
    else:
        net = net.cuda()
    
    logger.step("Model loaded successfully")
    return net


def predict_single_image(model, img, config, save_path, logger):
    """Predict segmentation for a single image"""
    # Handle color images
    if len(img.shape) > 2:
        img = img[:, :, 0]
    
    h, w = img.shape
    imageSize = config.imageSize
    
    # Auto-adjust image size if needed
    if imageSize == 0:
        imageSize = 250
        while imageSize + 250 < h and imageSize + 250 < w:
            imageSize += 250
    
    # Create image patches
    images = [
        img[:imageSize, :imageSize],
        img[h-imageSize:, :imageSize],
        img[:imageSize, w-imageSize:],
        img[h-imageSize:, w-imageSize:]
    ]
    
    # Process patches
    proc_images = []
    for idx, sub_img in enumerate(images):
        pil_sub_img = Image.fromarray((sub_img * 255).astype('uint8'))
        sub_tensor = toTensor(pil_sub_img).unsqueeze(0)
        
        with torch.no_grad():
            sr = model(sub_tensor.cuda() if not config.cpu else sub_tensor)
            sr = sr.cpu()
            m = nn.LogSoftmax(dim=1)
            sr = m(sr).argmax(dim=1, keepdim=True)
            sr_squeezed = sr.squeeze(0)
            pil_sr_img = toPIL(sr_squeezed.float() / (config.nch_out - 1))
            proc_images.append(pil_sr_img)
    
    # Stitch patches
    img1, img2, img3, img4 = proc_images
    woffset = (2 * imageSize - w) // 2
    hoffset = (2 * imageSize - h) // 2
    
    img1 = np.array(img1)[:imageSize-hoffset, :imageSize-woffset]
    img3 = np.array(img3)[:imageSize-hoffset, woffset:]
    top = np.concatenate((img1, img3), axis=1)
    
    img2 = np.array(img2)[hoffset:, :imageSize-woffset]
    img4 = np.array(img4)[hoffset:, woffset:]
    bot = np.concatenate((img2, img4), axis=1)
    
    result = np.concatenate((top, bot), axis=0)
    
    # Clean borders
    result[:10, :] = 0
    result[-10:, :] = 0
    result[:, :10] = 0
    result[:, -10:] = 0
    
    # Save result
    Image.fromarray(result.astype('uint8')).save(save_path)
    logger.step(f"Result saved: {save_path}")
    
    return result


# ═══════════════════════════════════════════════════════════════════
# Main Function
# ═══════════════════════════════════════════════════════════════════

def main():
    """Main entry point"""
    
    # Step 1: Initialize and parse arguments
    args = parse_args()
    set_seed(args.seed)
    
    if args.output_dir is None:
        args.output_dir = os.path.join(args.main_path, "outputs", "er_segmentation")
    
    outdir = ensure_outdir(args.output_dir)
    
    # Step 2: Initialize logger
    logger = QuickLogger(name="sim_er_segmentation", log_dir=outdir)
    logger.step("SIM ER Segmentation")
    
    # Step 3: Build paths and config
    ernet_dir = os.path.join(args.main_path, 'ipa', 'processing', 'segmentation', 
                             'segmentation_sim_wfm', 'ERNet')
    
    if args.weights_path is None:
        args.weights_path = os.path.join(ernet_dir, 'pretrained_model', 'final.pth')
    
    if args.image_path is None:
        args.image_path = os.path.join(ernet_dir, 'testdata', 'input', 
                                       '190624_COS7_3_LAMP1-mCherry_ER-YFP_optokin_2s_40x_1pt5_01_2_MMStack_Pos0-10000.png')
    
    logger.step(f"Main path: {args.main_path}")
    logger.step(f"Output: {outdir}")
    logger.step(f"Model: {args.weights_path}")
    logger.step(f"Image: {args.image_path}")
    logger.step(f"Device: {'CPU' if args.use_cpu else 'GPU'}")
    
    config = create_pred_config(
        model_path=args.weights_path,
        image_size=args.image_size,
        nch_in=args.nch_in,
        nch_out=args.nch_out,
        n_resgroups=args.n_resgroups,
        n_resblocks=args.n_resblocks,
        use_cpu=args.use_cpu
    )
    
    # Step 4: Save configuration
    config_data = {
        'config': {
            'weights_path': args.weights_path,
            'image_path': args.image_path,
            'image_size': args.image_size,
            'use_cpu': args.use_cpu,
            'nch_in': args.nch_in,
            'nch_out': args.nch_out,
        },
        'seed': args.seed,
        'timestamp': time.strftime("%Y-%m-%d %H:%M:%S")
    }
    save_json(outdir, "config", config_data)
    
    # Step 5: Validate input files
    check_exists(logger, weights=args.weights_path, image=args.image_path)
    
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
            logger.step(f"{name} failed: out of memory")
            return None
        except Exception as e:
            logger.step(f"{name} failed: {type(e).__name__}: {e}")
            return None
    
    # Step 7: Conditional execution
    if args.run_segmentation:
        if not os.path.exists(args.weights_path):
            logger.step(f"Error: Model weights not found: {args.weights_path}")
        elif not os.path.exists(args.image_path):
            logger.step(f"Error: Image not found: {args.image_path}")
        else:
            # Load model
            model = timed_run("model_loading", load_model, config=config, logger=logger)
            
            if model is not None:
                # Load image
                img = np.array(Image.open(args.image_path)) / 255.0
                logger.step(f"Image loaded: shape={img.shape}")
                
                # Run prediction
                output_path = os.path.join(outdir, 'er_segmentation_result.png')
                result = timed_run(
                    "er_segmentation",
                    predict_single_image,
                    model=model,
                    img=img,
                    config=config,
                    save_path=output_path,
                    logger=logger
                )
                
                if result is not None:
                    save_json(outdir, "segmentation_results", {
                        'output_path': output_path,
                        'result_shape': result.shape,
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
