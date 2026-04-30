"""
SIM Radial Distribution Function (RDF) Analysis Demo

This script calculates radial distribution functions (RDF) for organelles within 
cellular partitions using SIM data. It analyzes the spatial distribution of 
organelles (e.g., ISG) relative to cellular compartments (nucleus to plasma membrane).
"""

import os
import sys
import numpy as np

# Add project root to path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
sys.path.insert(0, PROJECT_ROOT)

from ipa.analysis import (
    load_partition_coords_from_xvg, 
    calculate_rdf_from_xvg, 
    save_radial_rdf_results
)
from ipa.processing.partitioning import plot_radial_rdf
from ipa.data_loader import QuickLogger, UniversalDataLoader

def main():
    # Initialize logger
    output_dir = os.path.join(PROJECT_ROOT, 'results', 'sim_rdf_analysis_demo')
    os.makedirs(output_dir, exist_ok=True)

    logger = QuickLogger("sim_rdf_analysis", log_dir=os.path.join(PROJECT_ROOT, 'logs'))
    logger.step("SIM Radial Distribution Function (RDF) Analysis Demo")
    
    dataid = '20220909_30-2-1-SIM'
    logger.step(f"Processing dataset: {dataid}")
    
    # Paths
    results_dir = os.path.join(PROJECT_ROOT, 'data', 'SIM', 'results')
    isg_file_path = os.path.join(PROJECT_ROOT, 'data', 'sim', f'{dataid}_ISG_seg.tif')
    partition_coords_file = os.path.join(results_dir, f'{dataid}_partition_coords.xvg')

    if not os.path.exists(partition_coords_file):
        logger.step(f"Error: Partition coordinates file not found: {partition_coords_file}")
        return
    
    logger.file_in(partition_coords_file)
    logger.step("Loading partition coordinates and radial positions from XVG...")
    partition_coords, radial_positions = load_partition_coords_from_xvg(partition_coords_file)
    
    logger.step("Loading ISG organelle data...")
    if not os.path.exists(isg_file_path):
        logger.step(f"Error: ISG file not found: {isg_file_path}")
        return
    
    logger.file_in(isg_file_path)
    isg_mask = UniversalDataLoader.load_data(isg_file_path)
    
    # Flip the ISG data to match coordinate system
    isg_mask = np.flip(isg_mask, axis=0)
    
    logger.step(f"Loaded ISG mask with shape: {isg_mask.shape}")
    
    # Calculate RDF directly from XVG data
    logger.step("Calculating RDF...")
    rdf_results = calculate_rdf_from_xvg(isg_mask, partition_coords, radial_positions, bins=8)
    
    if rdf_results is None:
        logger.step("Error: Failed to calculate RDF")
        return
    
    logger.step("RDF RESULTS SUMMARY")
    for i, (pos, rdf_val) in enumerate(zip(rdf_results['radii'], rdf_results['rdf'])):
        logger.step(f"Bin {i+1}: Position {pos:.3f} -> RDF {rdf_val:.4f}")
    
    # Plot results
    logger.step("Generating plots...")
    rdf_plot_path = os.path.join(output_dir, f'{dataid}_radial_rdf_analysis.png')
    plot_radial_rdf(rdf_results, save_path=rdf_plot_path)
    logger.file_out(rdf_plot_path)
    
    # Save results
    logger.step("Saving results...")
    save_radial_rdf_results(rdf_results, output_dir, dataid)
    
    logger.step("Radial RDF analysis completed successfully!")

if __name__ == '__main__':
    main()
