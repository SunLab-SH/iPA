
import argparse
import os
import sys

# Remove sys.path hacking and use relative imports
try:
    from .src.noise2noise_wrapper import train_noise2noise
except ImportError:
    import sys
    import os
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(current_dir))))
    if parent_dir not in sys.path:
        sys.path.append(parent_dir)
    from ipa.processing.denoising.n2n_package.src.noise2noise_wrapper import train_noise2noise


def parse_args():
    """Command-line argument parser"""
    parser = argparse.ArgumentParser(description='Noise2Noise training')
    
    # Data parameters
    parser.add_argument('-t', '--train-dir', help='training set path', required=True)
    parser.add_argument('-v', '--valid-dir', help='validation set path', required=True)
    parser.add_argument('--ckpt-save-path', help='model save path', default='./ckpts')
    parser.add_argument('--ckpt-overwrite', help='overwrite saved model', action='store_true')
    parser.add_argument('--report-interval', help='report interval', default=100, type=int)
    parser.add_argument('-ts', '--train-size', help='training set size', type=int)
    parser.add_argument('-vs', '--valid-size', help='validation set size', type=int)
    
    # Training hyperparameters
    parser.add_argument('-lr', '--learning-rate', help='learning rate', default=0.001, type=float)
    parser.add_argument('-a', '--adam', help='adam parameters', nargs='+', default=[0.9, 0.99, 1e-8], type=list)
    parser.add_argument('-b', '--batch-size', help='batch size', default=4, type=int)
    parser.add_argument('-e', '--nb-epochs', help='number of epochs', default=50, type=int)
    parser.add_argument('-l', '--loss', help='loss function', choices=['l1', 'l2', 'hdr'], default='l1', type=str)
    parser.add_argument('--cuda', help='use CUDA', action='store_true')
    parser.add_argument('--plot-stats', help='plot statistics', action='store_true')
    
    # Noise parameters
    parser.add_argument('-n', '--noise-type', help='noise type',
                       choices=['gaussian', 'poisson', 'text', 'mc'], default='gaussian', type=str)
    parser.add_argument('-p', '--noise-param', help='noise parameter (e.g. std for gaussian)', default=25, type=float)
    parser.add_argument('-s', '--seed', help='random seed', type=int)
    parser.add_argument('-c', '--crop-size', help='crop size', default=128, type=int)
    parser.add_argument('--clean-targets', help='use clean targets for training', action='store_true')
    
    return parser.parse_args()

def main():
    params = parse_args()
    param_dict = vars(params)
    
    train_noise2noise(
        train_dir=params.train_dir,
        valid_dir=params.valid_dir,
        **param_dict
    )

if __name__ == '__main__':
    main()
