
import os
import sys
import argparse
import numpy as np
import mrcfile
import matplotlib.pyplot as plt

# Add project root to path to import ipa
current_dir = os.path.dirname(os.path.abspath(__file__))
ipa_root = os.path.dirname(os.path.dirname(current_dir))
sys.path.insert(0, ipa_root)

from ipa.processing.denoising.n2n_package import predict_noise2noise
from ipa.data_loader import QuickLogger
from parsers import get_args

def main():
    """
    Main function for SXT N2N denoising prediction
    """
    # Get main path and initialize logger
    mainpath = get_args().main_path
    log_dir = f'{mainpath}/logs'
    logger = QuickLogger("sxt_n2n_denoising_pred", log_dir=log_dir)
    
    logger.step("Starting SXT N2N Denoising Prediction Demo")
    
    # 1. Setup Data Paths
    default_data_dir = os.path.join(mainpath, 'data', 'sxt_images')
    
    parser = argparse.ArgumentParser(description="Predict N2N model")
    parser.add_argument('--input_dir', type=str, default=default_data_dir, help='Path to input data directory')
    parser.add_argument('--model_path', type=str, help='Path to trained model .pt file')
    
    args, unknown = parser.parse_known_args()
    
    # 2. Setup Model Path
    if args.model_path:
        model_path = args.model_path
    else:
        # Try to find a default trained model
        model_path = os.path.join(mainpath, 'ipa', 'processing', 'denoising', 'models', 'n2n_sxt_trained.pt')
    
    if not os.path.exists(model_path):
        logger.step(f"Error: Model file not found at {model_path}. Please run training script first.")
        return
        
    logger.step(f"Loading model from: {model_path}")
    logger.step(f"Processing images in: {args.input_dir}")
    
    # 3. Predict
    # predict_noise2noise saves results directly to the data dir with suffix
    predict_noise2noise(
        data_dir=args.input_dir,
        model_path=model_path,
        cuda=True,
        show_output=-1 
    )
    
    logger.step("Prediction finished. Results saved in input directory.")

if __name__ == "__main__":
    main()
