import argparse
import os

def get_args():
    parser = argparse.ArgumentParser(description='example Analysis')
    
    current_dir = os.path.dirname(os.path.abspath(__file__))
    default_main_path = os.path.dirname(os.path.dirname(current_dir))

    parser.add_argument('--main_path',
                       type=str,
                       default=default_main_path,
                       help='Main path to iPA directory')
    
    return parser.parse_args()
