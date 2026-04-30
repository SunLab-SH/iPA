"""
SXT Spatial Arrangement Analysis Demo

This script demonstrates spatial arrangement analysis using SXT data.
It calculates Radial Distribution Functions (RDF) based on NE-PM partitions.
"""

import os
import sys
import numpy as np
import tifffile
import pandas as pd

# Add project root to path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
sys.path.insert(0, PROJECT_ROOT)

from ipa.analysis import calculate_rdf_from_xvg
from ipa.processing.partitioning import Partitioning
from ipa.data_loader import QuickLogger

def main():
    # Initialize logger
    output_dir = os.path.join(PROJECT_ROOT, 'results', 'sxt_spatial_arrangement_demo')
    os.makedirs(output_dir, exist_ok=True)
    
    logger = QuickLogger("sxt_spatial_arrangement", log_dir=os.path.join(PROJECT_ROOT, 'logs'))
    logger.step("Starting SXT Spatial Arrangement Analysis Demo")

    # --- Configuration ---
    data_id = '784_5'
    data_root = os.path.join(PROJECT_ROOT, 'data', 'sxt')
    
    isg_mask_path = os.path.join(data_root, f'{data_id}_isg_label.tiff')
    pm_mask_path = os.path.join(data_root, f'{data_id}_wholecell_label.tiff')
    ne_mask_path = os.path.join(data_root, f'{data_id}_NC_label.tiff')
    
    if not all(os.path.exists(p) for p in [pm_mask_path, ne_mask_path, isg_mask_path]):
        logger.step("Error: Required mask files not found.")
        return

    pm_mask = tifffile.imread(pm_mask_path)
    ne_mask = tifffile.imread(ne_mask_path)
    isg_mask = tifffile.imread(isg_mask_path)

    # Create partitions
    logger.step("Creating NE-PM radial partitions...")
    partitioner = Partitioning(root_dir=output_dir, n_slices=8)
    center, ne_edge, pm_edge = partitioner.extract_ne_pm_edges(pm_mask, ne_mask)
    partition_mask = partitioner.create_nepm_radial_partitions(
        ne_edge, pm_edge, shape=pm_mask.shape, n_slices=8, 
        pm_mask=pm_mask, ne_mask=ne_mask
    )
    
    # Prepare inputs for RDF calculation
    partition_coords = {}
    radial_positions = {}
    
    for i in range(1, 9):
        coords = np.array(np.where(partition_mask == i)).T
        partition_coords[i] = coords
        norm_pos = (i - 0.5) / 8.0
        for coord in coords:
            radial_positions[tuple(coord)] = norm_pos
        
    # Calculate RDF
    logger.step("Calculating RDF...")
    rdf_results = calculate_rdf_from_xvg(isg_mask, partition_coords, radial_positions, bins=8)
    
    if rdf_results:
        peak_pos = rdf_results['radii'][np.argmax(rdf_results['rdf'])]
        logger.step(f"   Peak density at radial position: {peak_pos:.2f}")
        
        # Save RDF results
        rdf_df = pd.DataFrame({
            'Radial_Position': rdf_results['radii'],
            'g(r)': rdf_results['rdf']
        })
        rdf_df.to_csv(os.path.join(output_dir, f'{data_id}_rdf_results.csv'), index=False)
        logger.file_out(os.path.join(output_dir, f'{data_id}_rdf_results.csv'))

    logger.step("Analysis Completed!")

if __name__ == '__main__':
    main()
