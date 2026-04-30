import os
import sys
import numpy as np
import mrcfile
import torch
import argparse

# Add iPA module path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
sys.path.insert(0, PROJECT_ROOT)

from ipa.processing.denoising.noise2void.n2v_wrapper import N2V
from ipa.data_loader import UniversalDataLoader
from ipa.data_loader import QuickLogger

def main():
    """
    Main function for SXT N2V denoising training
    """
    # Initialize logger
    log_dir = os.path.join(PROJECT_ROOT, 'logs')
    logger = QuickLogger("sxt_denoising_train", log_dir=log_dir)
    
    logger.step("Starting SXT N2V Denoising Training Demo")
    
    # 1. Setup Data Paths
    # Default SXT example file
    default_input = os.path.join(PROJECT_ROOT, 'data', 'sxt', 'Stevens_pancreatic_INS_1E_784_5_pre_rec.mrc')
    
    parser = argparse.ArgumentParser(description="Train N2V model (Reproduction)")
    parser.add_argument('--input', type=str, default=default_input, help='Path to input data')
    parser.add_argument('--output_dir', type=str, default=None, help='Directory to save model')
    parser.add_argument('--epochs', type=int, default=5, help='Number of epochs')
    parser.add_argument('--batch_size', type=int, default=2, help='Batch size')
    parser.add_argument('--lr', type=float, default=0.0004, help='Learning rate')
    parser.add_argument('--n_filters', type=int, default=32, help='Number of filters')
    
    args, unknown = parser.parse_known_args()
    
    input_file = args.input
    if not os.path.exists(input_file):
        logger.step(f"Error: Input file not found: {input_file}")
        return

    logger.file_in(input_file)
    logger.step(f"Input file: {input_file}")

    # 2. Load Data
    logger.step("Reading input image...")
    img = UniversalDataLoader.load_data(input_file)
    
    logger.step(f"Image shape: {img.shape}")
    logger.step(f"Image data type: {img.dtype}")
    
    # Ensure float32
    if img.dtype != np.float32:
        img = img.astype(np.float32)

    # 3. Setup Model Output
    if args.output_dir:
        checkpoints_dir = args.output_dir
    else:
        # Default to n2v_3D directory
        checkpoints_dir = os.path.join(mainpath, 'ipa', 'processing', 'denoising', 'models', 'n2v_3D')
    
    logger.step(f"Model will be saved to: {checkpoints_dir}")
    
    # 4. Train
    # Reproducing n2v_3D parameters
    n2v = N2V(
        model_name='sxt_repro_run', 
        checkpoint_dir=checkpoints_dir, 
        n_channels=1,
        n_filters=args.n_filters 
    )
    
    # Split Data
    num_slices = img.shape[0]
    split_idx = int(num_slices * 0.9)
    train_data = img[:split_idx]
    val_data = img[split_idx:]
    
    logger.step(f"Training on {train_data.shape[0]} slices, validating on {val_data.shape[0]} slices.")
    logger.step("Starting training...")
    
    history = n2v.train(
        train_data, 
        val_data, 
        batch_size=args.batch_size,   
        epochs=args.epochs,       
        lr=args.lr       
    )
    
    logger.step("Training finished.")

if __name__ == "__main__":
    main()
