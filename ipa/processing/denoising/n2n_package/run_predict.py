
import argparse
import os
import sys

# Remove sys.path hacking and use relative imports
try:
    from .src.noise2noise_wrapper import predict_noise2noise
except ImportError:
    # If run as a script, we might need to add the project root to path
    # assuming we are running from the file location
    import sys
    import os
    
    # Try to resolve relative import by making this directory a package context?
    # Better: Assume user runs this via -m or has ipa installed.
    # Fallback for direct script execution without package context (discouraged but common)
    # We will try to find the package root
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(current_dir))))
    if parent_dir not in sys.path:
        sys.path.append(parent_dir)
    
    from ipa.processing.denoising.n2n_package.src.noise2noise_wrapper import predict_noise2noise


def parse_args():
    """Command-line argument parser"""
    parser = argparse.ArgumentParser(description='Noise2Noise prediction')
    
    # Data parameters
    parser.add_argument('-d', '--data-dir', help='test data path', required=True)
    parser.add_argument('--load-ckpt', help='load model checkpoint path', required=True)
    parser.add_argument('--show-output', help='show output results', default=0, type=int)
    parser.add_argument('--cuda', help='use CUDA', action='store_true')
    
    return parser.parse_args()

def main():
    params = parse_args()
    param_dict = vars(params)
    
    predict_noise2noise(
        data_dir=params.data_dir,
        model_path=params.load_ckpt,
        **param_dict
    )

if __name__ == '__main__':
    main()
