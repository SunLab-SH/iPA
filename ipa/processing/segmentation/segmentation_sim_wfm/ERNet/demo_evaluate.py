#!/usr/bin/env python3
"""
Demo script for ERNet evaluation with external data and parameters
"""

import os
import glob
from types import SimpleNamespace
from evaluate import EvaluateModel

def create_opt_config(
    input_path, 
    weights_path, 
    output_path, 
    image_size=1000,
    model='rcan',
    nch_in=1,
    nch_out=2,
    n_resgroups=5,
    n_resblocks=10,
    use_cpu=False,
    undomulti=False
):
    """
    Create configuration object for ERNet evaluation
    
    Args:
        input_path: Path to input images directory or single image file
        weights_path: Path to model weights file
        output_path: Path to output directory
        image_size: Size for image processing (default: 1000)
        model: Model type (default: 'rcan')
        nch_in: Number of input channels (default: 1)
        nch_out: Number of output channels (default: 2)
        n_resgroups: Number of residual groups (default: 5)
        n_resblocks: Number of residual blocks (default: 10)
        use_cpu: Use CPU instead of GPU (default: False)
        undomulti: Remove DataParallel wrapper (default: False)
    
    Returns:
        SimpleNamespace: Configuration object
    """
    
    # Create configuration object
    opt = SimpleNamespace()
    
    # Core parameters
    opt.root = input_path
    opt.weights = weights_path
    opt.out = output_path
    opt.imageSize = image_size
    opt.model = model
    opt.nch_in = nch_in
    opt.nch_out = nch_out
    opt.n_resgroups = n_resgroups
    opt.n_resblocks = n_resblocks
    opt.cpu = use_cpu
    opt.undomulti = undomulti
    
    # Additional required parameters with defaults
    opt.n_feats = 64
    opt.reduction = 16
    opt.narch = 0
    opt.multigpu = False
    opt.dataset = 'imagedataset'
    opt.batchSize = 16
    opt.batchSize_test = 1
    opt.lr = 0.0001
    opt.nepoch = 10
    opt.ntrain = 0
    opt.ntest = 10
    opt.test = False
    opt.log = False
    opt.saveinterval = 10
    opt.testinterval = 1
    opt.plotinterval = 1
    opt.scheduler = ''
    opt.workers = 6
    
    return opt


def run_evaluation_demo():
    """
    Demo function showing how to use ERNet evaluation with external parameters
    """
    
    # Get script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Define paths - you can modify these for your data
    input_images_path = os.path.join(script_dir, 'testdata', 'input')
    weights_file_path = os.path.join(script_dir, 'pretrained_model', 'final.pth')
    output_results_path = os.path.join(script_dir, 'demo_output')
    
    # Check if input data exists
    if not os.path.exists(input_images_path):
        print(f"Warning: Input path '{input_images_path}' does not exist")
        print("Please make sure you have test images in the input directory")
        return
    
    if not os.path.exists(weights_file_path):
        print(f"Warning: Weights file '{weights_file_path}' does not exist")
        print("Please make sure you have the pretrained model file")
        return
    
    # Create output directory if it doesn't exist
    os.makedirs(output_results_path, exist_ok=True)
    
    # Create configuration with your parameters
    config = create_opt_config(
        input_path=input_images_path,
        weights_path=weights_file_path,
        output_path=output_results_path,
        image_size=1000,  # You can adjust this
        model='rcan',
        nch_in=1,         # Grayscale input
        nch_out=2,        # Binary segmentation output
        n_resgroups=5,
        n_resblocks=10,
        use_cpu=False,    # Set to True if you want to use CPU
        undomulti=False   # Set to True if model was saved with DataParallel
    )
    
    print("=" * 60)
    print("ERNet Evaluation Demo")
    print("=" * 60)
    print(f"Input images: {config.root}")
    print(f"Model weights: {config.weights}")
    print(f"Output directory: {config.out}")
    print(f"Image size: {config.imageSize}")
    print(f"Model: {config.model}")
    print(f"Input channels: {config.nch_in}")
    print(f"Output channels: {config.nch_out}")
    print(f"Using CPU: {config.cpu}")
    print("=" * 60)
    
    # Check available images
    image_extensions = ['*.jpg', '*.png', '*.tif']
    all_images = []
    for ext in image_extensions:
        all_images.extend(glob.glob(os.path.join(input_images_path, ext)))
        all_images.extend(glob.glob(os.path.join(input_images_path, '**', ext), recursive=True))
    
    print(f"Found {len(all_images)} images to process")
    if all_images:
        print("Sample images:")
        for img in all_images[:3]:  # Show first 3 images
            print(f"  - {os.path.basename(img)}")
    print("=" * 60)
    
    # Run evaluation
    try:
        EvaluateModel(config)
        print("\nEvaluation completed successfully!")
        print(f"Results saved to: {output_results_path}")
    except Exception as e:
        print(f"Error during evaluation: {str(e)}")
        raise


def evaluate_custom_data(input_path, weights_path, output_path, **kwargs):
    """
    Convenience function to evaluate custom data
    
    Args:
        input_path: Path to input images
        weights_path: Path to model weights
        output_path: Path to save results
        **kwargs: Additional parameters (image_size, model, etc.)
    
    Example usage:
        evaluate_custom_data(
            input_path='/path/to/images',
            weights_path='/path/to/model.pth',
            output_path='/path/to/results',
            image_size=512,
            use_cpu=True
        )
    """
    
    # Create configuration
    config = create_opt_config(
        input_path=input_path,
        weights_path=weights_path,
        output_path=output_path,
        **kwargs
    )
    
    # Run evaluation
    EvaluateModel(config)
    print(f"Evaluation completed. Results saved to: {output_path}")


if __name__ == '__main__':
    # Run the demo
    run_evaluation_demo()
    
    # Example of how to use with custom data:
    # evaluate_custom_data(
    #     input_path='path/to/your/images',
    #     weights_path='path/to/your/model.pth',
    #     output_path='path/to/your/results',
    #     image_size=512,
    #     use_cpu=True
    # )