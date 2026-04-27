
import os
import sys
import argparse

# Add project root to path to import ipa
current_dir = os.path.dirname(os.path.abspath(__file__))
ipa_root = os.path.dirname(os.path.dirname(current_dir))
sys.path.insert(0, ipa_root)

from ipa.processing.denoising.n2n_package.src.noise2noise_wrapper import train_noise2noise
from ipa.data_loader import QuickLogger
from parsers import get_args

def main():
    """
    Main function for SXT N2N denoising training
    """
    # Get main path and initialize logger
    # Use parse_known_args in get_args or handle conflict manually. 
    # Here we can just manually set mainpath if we want to avoid conflict, or use get_args() first.
    # The issue is get_args() calls parser.parse_args() which fails on unknown args.
    
    # Simple fix: Re-implement necessary path logic locally or use parse_known_args in a wrapper
    # For this demo, let's just get the path directly without using the shared parser that is strict.
    current_dir = os.path.dirname(os.path.abspath(__file__))
    mainpath = os.path.dirname(os.path.dirname(current_dir))
    
    log_dir = f'{mainpath}/logs'
    logger = QuickLogger("sxt_n2n_denoising_train", log_dir=log_dir)
    
    logger.step("Starting SXT N2N Denoising Training Demo")
    
    # 1. Setup Data Paths
    # Since N2N typically requires a directory of images, we point to the directory containing SXT data
    # For this demo, we'll assume the sxt_images folder is the training source
    default_data_dir = os.path.join(mainpath, 'data', 'sxt_images')
    
    parser = argparse.ArgumentParser(description="Train N2N model")
    parser.add_argument('--train_dir', type=str, default=default_data_dir, help='Path to training data directory')
    parser.add_argument('--valid_dir', type=str, default=default_data_dir, help='Path to validation data directory')
    parser.add_argument('--output_dir', type=str, default=None, help='Directory to save model')
    parser.add_argument('--epochs', type=int, default=50, help='Number of epochs')
    parser.add_argument('--batch_size', type=int, default=4, help='Batch size')
    parser.add_argument('--lr', type=float, default=0.001, help='Learning rate')
    
    args, unknown = parser.parse_known_args()
    
    # 2. Setup Model Output
    if args.output_dir:
        checkpoints_dir = args.output_dir
    else:
        checkpoints_dir = os.path.join(mainpath, 'ipa', 'processing', 'denoising', 'models', 'n2n')
    
    logger.step(f"Training data dir: {args.train_dir}")
    print(f"DEBUG: mainpath={mainpath}")
    print(f"DEBUG: train_dir={args.train_dir}")
    print(f"DEBUG: valid_dir={args.valid_dir}")
    print(f"DEBUG: checkpoints_dir={checkpoints_dir}")
    logger.step(f"Model will be saved to: {checkpoints_dir}")
    
    # 3. Train
    # Using Gaussian noise augmentation as default strategy for single noisy inputs if needed,
    # or assuming the folder contains pairs if formatted correctly.
    # The wrapper handles basic loading.
    
    logger.step("Starting training...")
    
    train_noise2noise(
        train_dir=args.train_dir,
        valid_dir=args.valid_dir,
        ckpt_save_path=checkpoints_dir,
        nb_epochs=args.epochs,
        batch_size=args.batch_size,
        learning_rate=args.lr,
        loss='l2',
        cuda=True,
        noise_type='gaussian', # Strategy for self-supervised/simulated noise
        noise_param=25,
        crop_size=64,
        report_interval=1
    )
    
    logger.step("Training finished.")

if __name__ == "__main__":
    main()
