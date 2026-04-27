"""
Radial Distribution Function (RDF) Analysis for Cellular Partitions

This script calculates radial distribution functions for organelles within 
cellular partitions, using normalized radial position distribution method.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import tifffile
from scipy.spatial.distance import cdist
import pickle
import sys
from parsers import get_args
from scipy.spatial import cKDTree
from ipa.analysis import load_partition_coords_from_xvg, calculate_rdf_from_xvg, save_radial_rdf_results
from ipa.processing.partitioning import plot_radial_rdf
from ipa.data_loader import QuickLogger


def main():
    mainpath = get_args().main_path
    
    # Initialize logger
    log_dir = f'{mainpath}/logs'
    logger = QuickLogger("sxt_rdf_analysis", log_dir=log_dir)
    logger.step("Starting SXT RDF Analysis Demo")
    
    dataid = '784_5'
    logger.step(f"Processing dataset: {dataid}")
    
    # Paths
    results_dir = f'{mainpath}/data/sxt_images/results/'
    organelle_path = f'{mainpath}/data/sxt_images/{dataid}_isg_label.tiff'

    # Load partition coordinates from XVG file
    partition_coords_file = f'{results_dir}{dataid}_partition_coords.xvg'

    if not os.path.exists(partition_coords_file):
        logger.step(f"Error: Partition coordinates file not found: {partition_coords_file}")
        logger.step("Please run the partitioning analysis first")
        return
    
    logger.file_in(partition_coords_file)
    logger.step("Loading partition coordinates and radial positions from XVG...")
    partition_coords, radial_positions = load_partition_coords_from_xvg(partition_coords_file)
    
    logger.step("Loading organelle data...")
    # Check if organelle file exists
    if not os.path.exists(organelle_path):
        logger.step(f"Error: Organelle file not found: {organelle_path}")
        logger.step("Available files in directory:")
        sxt_dir = f'{mainpath}/data/sxt_images/'
        if os.path.exists(sxt_dir):
            for file in os.listdir(sxt_dir):
                if dataid in file:
                    logger.step(f"  {file}")
        return
    
    # Use tifffile to read TIFF format
    logger.file_in(organelle_path)
    organelle_mask = tifffile.imread(organelle_path)
    logger.step(f"Loaded organelle mask with shape: {organelle_mask.shape}")
    logger.step(f"Organelle mask dtype: {organelle_mask.dtype}")
    logger.step(f"Organelle mask value range: {organelle_mask.min()} to {organelle_mask.max()}")
    
    logger.step(f"Loaded partition coordinates for partitions: {list(partition_coords.keys())}")
    logger.step(f"Total radial positions: {len(radial_positions)}")
    
    # Calculate RDF directly from XVG data
    logger.step("Calculating RDF directly from XVG data...")
    rdf_results = calculate_rdf_from_xvg(organelle_mask, partition_coords, radial_positions, bins=8)
    
    if rdf_results is None:
        logger.step("Error: Failed to calculate RDF")
        return
    
    # Print RDF values before plotting
    logger.step("RDF RESULTS SUMMARY")
    logger.step(f"Dataset: {rdf_results['dataset']}")

    logger.step("RDF Values by Radial Position:")
    for i, (pos, rdf_val) in enumerate(zip(rdf_results['radii'], rdf_results['rdf'])):
        logger.step(f"Bin {i+1}: Position {pos:.3f} -> RDF {rdf_val:.4f}")
    
    # Plot results
    logger.step("Generating plots...")
    rdf_plot_path = f'{results_dir}{dataid}_radial_rdf_analysis.png'
    plot_radial_rdf(rdf_results, save_path=rdf_plot_path)
    logger.file_out(rdf_plot_path)
    
    # Save results
    logger.step("Saving results...")
    save_radial_rdf_results(rdf_results, results_dir, dataid)
    
    logger.step("Radial RDF analysis completed successfully!")
    logger.step(f"Final RDF info: {rdf_results}")

if __name__ == '__main__':
    main()
