
import os
import sys
import numpy as np
import mrcfile
import matplotlib.pyplot as plt
import argparse

# Add iPA module path
current_dir = os.path.dirname(os.path.abspath(__file__))
ipa_root = os.path.dirname(os.path.dirname(current_dir))
sys.path.insert(0, ipa_root)

from ipa.processing.denoising.noise2void.n2v_wrapper import N2V
from ipa.data_loader import UniversalDataLoader
from ipa.data_loader import QuickLogger
from parsers import get_args

def visualize_results(original, denoised, slice_idx=None):
    """
    Visualize original and denoised images
    """
    if slice_idx is None:
        slice_idx = original.shape[0] // 2
    
    plt.figure(figsize=(16, 8))
    
    # Display original image
    plt.subplot(1, 2, 1)
    plt.imshow(original[slice_idx], 
               cmap='magma',
               vmin=np.percentile(original, 0.1),
               vmax=np.percentile(original, 99.9))
    plt.title('Original (Noisy)')
    plt.colorbar()
    
    # Display denoised image
    plt.subplot(1, 2, 2)
    plt.imshow(denoised[slice_idx],
               cmap='magma', 
               vmin=np.percentile(denoised, 0.1),
               vmax=np.percentile(denoised, 99.9))
    plt.title('Denoised (N2V)')
    plt.colorbar()
    
    plt.tight_layout()
    plt.show()

def main():
    """
    Main function for SXT N2V denoising prediction
    """
    # Get main path and initialize logger
    mainpath = get_args().main_path
    log_dir = f'{mainpath}/logs'
    logger = QuickLogger("sxt_denoising_pred", log_dir=log_dir)
    
    logger.step("Starting SXT N2V Denoising Prediction Demo")
    
    # 1. Setup Data Paths
    # Default SXT example file
    default_input = os.path.join(mainpath, 'data', 'sxt_images', 'Stevens_pancreatic_INS_1E_784_5_pre_rec.mrc')
    
    # Parse input from command line or use default
    parser = argparse.ArgumentParser(description="Predict N2V model")
    parser.add_argument('--input', type=str, default=default_input, help='Path to input data')
    parser.add_argument('--model_path', type=str, help='Path to trained model .pth file')
    parser.add_argument('--output_dir', type=str, default=None, help='Output directory')
    parser.add_argument('--n_filters', type=int, default=32, help='Number of filters (must match training)')
    
    # We need to parse known args to avoid conflict with get_args() if it parses sys.argv
    # But get_args() uses argparse too. Let's merge logic or just use what we need.
    # To keep it simple and compatible with existing structure:
    args, unknown = parser.parse_known_args()
    
    input_file = args.input
    if not os.path.exists(input_file):
        logger.step(f"Error: Input file not found: {input_file}")
        return

    # 2. Setup Output
    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = os.path.dirname(input_file)
    os.makedirs(output_dir, exist_ok=True)

    logger.file_in(input_file)
    logger.step(f"Input file: {input_file}")
    logger.step(f"Output directory: {output_dir}")

    # 3. Load Data
    logger.step("Reading input image...")
    img = UniversalDataLoader.load_data(input_file)
    
    logger.step(f"Image shape: {img.shape}")
    logger.step(f"Image data type: {img.dtype}")
    
    # Ensure float32
    if img.dtype != np.float32:
        img = img.astype(np.float32)

    # 4. Setup Model
    # Default model path if not provided
    if args.model_path:
        model_path = args.model_path
    else:
        # Try to find the reproduced model
        model_path = os.path.join(mainpath, 'ipa', 'processing', 'denoising', 'models', 'n2v_3D', 'sxt_repro_run_final.pth')
    
    if not os.path.exists(model_path):
        logger.step(f"Error: Model file not found at {model_path}. Please run training script first.")
        return
        
    logger.step(f"Loading model from: {model_path}")
    
    # Initialize N2V Wrapper
    # checkpoint_dir is not used for prediction, but required by init
    n2v = N2V(model_name='inference', checkpoint_dir='./', n_channels=1, n_filters=args.n_filters)
    n2v.load_model(model_path)
    
    # 5. Predict
    logger.step("Performing denoising...")
    denoised_img = n2v.predict(img, batch_size=2)
    
    logger.step(f"Denoised image shape: {denoised_img.shape}")
    
    # 6. Save Results
    input_basename = os.path.basename(input_file)
    output_filename = input_basename.replace('.mrc', '_denoised.mrc')
    if output_filename == input_basename: # handle case where .mrc not in name
        output_filename = input_basename + '_denoised.mrc'
        
    output_filepath = os.path.join(output_dir, output_filename)
    
    logger.step(f"Saving results to: {output_filepath}")
    with mrcfile.new(output_filepath, overwrite=True) as mrc:
        mrc.set_data(denoised_img)
    
    logger.file_out(output_filepath)
    logger.step("Save successful!")
    
    # 7. Visualize
    logger.step("Generating visualization results...")
    visualize_results(img, denoised_img)

if __name__ == "__main__":
    main()
