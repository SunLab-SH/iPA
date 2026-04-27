import glob
import tifffile
import pandas as pd
import os
import sys
import numpy as np


ipa_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ipa_dir not in sys.path:
    sys.path.append(ipa_dir)

from ipa.analysis import analyze_vesicle_clusters
from parsers import get_args
from ipa.data_loader import UniversalDataLoader
from ipa.data_loader import QuickLogger
from ipa.analysis import reconstruct_4d_mask




RESCALE_TABLE = {
    "16.7-30_P21": [2.1, 1.7],
}

def load_sample_data(file_path, rescale_table, coord_base_path, info_base_path, logger):
    name = os.path.basename(file_path)[:-4]
    logger.step(f"Loading sample: {name}")
    logger.file_in(file_path)
    
    img = UniversalDataLoader.load_data(file_path)
    img_ch2 = img[:, :, 1, :, :]
    logger.step(f"Image loaded with shape: {img.shape}")
    logger.step(f"Channel 2 extracted with shape: {img_ch2.shape}")

    coord_file_path = os.path.join(coord_base_path, f"{name}_Position.csv")
    logger.file_in(coord_file_path)
    logger.step(f"Loading coordinate data from: {coord_file_path}")
    data = pd.read_csv(coord_file_path, skiprows=3)
    logger.step(f"Coordinate data loaded with {len(data)} entries")

    res_file = os.path.join(info_base_path, f"{name}.txt")
    logger.file_in(res_file)
    logger.step(f"Loading resolution info from: {res_file}")
    with open(res_file, 'r') as f:
        for line in f:
            if 'Width:' in line:
                width = int(line.split('(')[1].split(')')[0])
            elif 'Height:' in line:
                height = int(line.split('(')[1].split(')')[0])
            elif 'Depth:' in line:
                depth = int(line.split('(')[1].split(')')[0])
    res = (width, height, depth)
    logger.step(f"Resolution extracted: {res}")
    
    rescalex, rescaley = rescale_table[name]
    logger.step(f"Rescale factors - X: {rescalex}, Y: {rescaley}")


    logger.step("Reconstructing image data...")
    img_reconstructed, ves_coords = reconstruct_4d_mask(data, res, rescalex, rescaley)
    img_reconstructed = img_reconstructed.transpose((0, 3, 2, 1))
    logger.step(f"Vesicle coordinates extracted: {len(ves_coords)} vesicles")

    return name, img_ch2, img_reconstructed, ves_coords

def main():
    args = get_args()
    
    # Initialize logger
    log_dir = f'{args.main_path}/logs'
    logger = QuickLogger("wfm_isg_cluster_analysis", log_dir=log_dir)
    logger.step("Starting WFM ISG Cluster Analysis Demo")
    
    data_path = os.path.join(args.main_path, 'data', 'wfm_images')
    file_path = os.path.join(data_path, '16.7-30_P21.tif')
    
    logger.step(f"Main path: {args.main_path}")
    logger.step(f"Data path: {data_path}")
    logger.step("Processing sample: 16.7-30_P21")

    name, img_ch2, img_reconstructed, ves_coords = load_sample_data(file_path, RESCALE_TABLE, data_path, data_path, logger)

    logger.step("Performing vesicle cluster analysis...")
    results = analyze_vesicle_clusters(name, img_ch2, img_reconstructed, ves_coords)
    logger.step("Cluster analysis completed successfully")

    # Create results directory if it doesn't exist
    results_dir = f'./example/results'
    os.makedirs(results_dir, exist_ok=True)

    save_path = f'{results_dir}/{name}_find_clusters.npy'
    np.save(save_path, results['final_clusters'])
    
    logger.file_out(save_path)
    logger.step(f"Results saved for {name} at {save_path}")
    logger.step(f"Number of clusters found: {len(results['final_clusters'])}")
    logger.step("ISG cluster analysis demo completed successfully")


if __name__ == "__main__":
    main()

