
import os
import argparse
import numpy as np
import mrcfile
import tifffile
import torch
from .n2v_wrapper import N2V

def load_data(path):
    """Load MRC, TIF or NPY data."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")
        
    if path.endswith('.mrc'):
        with mrcfile.open(path, permissive=True) as mrc:
            data = mrc.data
    elif path.endswith(('.tif', '.tiff')):
        data = tifffile.imread(path)
    elif path.endswith('.npy'):
        data = np.load(path)
    else:
        raise ValueError("Unsupported format. Use .mrc, .tif/tiff or .npy")
        
    return data

def main():
    parser = argparse.ArgumentParser(description="Train N2V model")
    parser.add_argument('--input', type=str, required=True, help='Path to input data (mrc, tif, npy)')
    parser.add_argument('--name', type=str, default='n2v_model', help='Model name')
    parser.add_argument('--output_dir', type=str, default='./checkpoints', help='Directory to save model')
    parser.add_argument('--epochs', type=int, default=20, help='Number of epochs')
    parser.add_argument('--batch_size', type=int, default=4, help='Batch size')
    parser.add_argument('--lr', type=float, default=1e-3, help='Learning rate')
    parser.add_argument('--n_filters', type=int, default=64, help='Number of UNet filters (default: 64)')
    parser.add_argument('--n_channels', type=int, default=1, help='Number of channels (default: 1)')
    parser.add_argument('--val_split', type=float, default=0.1, help='Validation split ratio (default: 0.1)')

    args = parser.parse_args()
    
    # Load Data
    print(f"Loading data from {args.input}...")
    data = load_data(args.input)
    print(f"Data shape: {data.shape}")
    
    # Ensure float32
    if data.dtype != np.float32:
        if data.max() > 255:
            # Assume 16-bit or float
            pass
        elif data.max() > 1:
            data = data.astype(np.float32) / 255.0
        else:
            data = data.astype(np.float32)

    # Split Train/Val
    num_samples = data.shape[0]
    split_idx = int(num_samples * (1 - args.val_split))
    train_data = data[:split_idx]
    val_data = data[split_idx:]
    
    print(f"Training samples: {train_data.shape[0]}, Validation samples: {val_data.shape[0]}")
    
    # Initialize Model
    # Use absolute path for checkpoint_dir if provided relative path starts with ./
    if args.output_dir.startswith('./'):
         args.output_dir = os.path.abspath(args.output_dir)
         
    n2v = N2V(model_name=args.name, checkpoint_dir=args.output_dir, n_channels=args.n_channels, n_filters=args.n_filters)
    
    # Train
    print("Starting training...")
    n2v.train(train_data, val_data, batch_size=args.batch_size, epochs=args.epochs, lr=args.lr)
    print("Training finished.")

if __name__ == '__main__':
    main()
