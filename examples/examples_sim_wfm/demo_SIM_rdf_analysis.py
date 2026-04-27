#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
SIM Radial Distribution Function (RDF) Analysis Demo

This script calculates radial distribution functions (RDF) for organelles within 
cellular partitions using SIM data. It analyzes the spatial distribution of 
organelles (e.g., ISG) relative to cellular compartments (nucleus to plasma membrane).

Usage:
    python demo_SIM_rdf_analysis.py

Requirements:
    - Input: 
        * Partition coordinate file (.xvg) from partitioning analysis
        * Organelle segmentation mask (.tif format)
    - Output: RDF plots and statistical results

The script performs:
1. Loads partition coordinates and radial positions
2. Loads organelle mask data
3. Calculates RDF across radial positions
4. Generates plots and saves results
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from parsers import get_args
from ipa.analysis import (
    load_partition_coords_from_xvg, 
    calculate_rdf_from_xvg, 
    save_radial_rdf_results
)
from ipa.processing.partitioning import plot_radial_rdf
from ipa.data_loader import QuickLogger, UniversalDataLoader


def main():
    """
    Main function to perform RDF analysis on SIM organelle data.
    
    This function:
    1. Loads partition coordinates from XVG file
    2. Loads organelle segmentation data
    3. Calculates radial distribution function
    4. Generates visualizations and saves results
    """
    mainpath = get_args().main_path
    
    # Initialize logger
    log_dir = f'{mainpath}/logs'
    logger = QuickLogger("sim_rdf_analysis", log_dir=log_dir)
    logger.step("=" * 60)
    logger.step("SIM Radial Distribution Function (RDF) Analysis Demo")
    logger.step("=" * 60)
    
    dataid = '20220909_30-2-1-SIM'
    logger.step(f"Processing dataset: {dataid}")
    
    # Paths
    results_dir = f'{mainpath}/data/SIM/results/'
    isg_file_path = f'{mainpath}/data/sim_images/{dataid}_ISG_seg.tif'

    # Load partition coordinates from XVG file
    partition_coords_file = f'{results_dir}{dataid}_partition_coords.xvg'

    if not os.path.exists(partition_coords_file):
        logger.step(f"Error: Partition coordinates file not found: {partition_coords_file}")
        logger.step("Please run the SIM partitioning analysis first")
        return
    
    logger.file_in(partition_coords_file)
    logger.step("Loading partition coordinates and radial positions from XVG...")
    partition_coords, radial_positions = load_partition_coords_from_xvg(partition_coords_file)
    
    logger.step("Loading ISG organelle data...")
    # Check if organelle file exists
    if not os.path.exists(isg_file_path):
        logger.step(f"Error: ISG file not found: {isg_file_path}")
        logger.step("Available files in directory:")
        sim_dir = f'{mainpath}/data/sim_images/'
        if os.path.exists(sim_dir):
            for file in os.listdir(sim_dir):
                if dataid in file:
                    logger.step(f"  {file}")
        return
    
    # Use UniversalDataLoader to read the ISG data
    logger.file_in(isg_file_path)
    isg_mask = UniversalDataLoader.load_data(isg_file_path)
    
    # Flip the ISG data to match coordinate system (as done in SIM visualization)
    isg_mask = np.flip(isg_mask, axis=0)
    
    logger.step(f"Loaded ISG mask with shape: {isg_mask.shape}")
    logger.step(f"ISG mask dtype: {isg_mask.dtype}")
    logger.step(f"ISG mask value range: {isg_mask.min()} to {isg_mask.max()}")
    
    logger.step(f"Loaded partition coordinates for partitions: {list(partition_coords.keys())}")
    logger.step(f"Total radial positions: {len(radial_positions)}")
    
    # Calculate RDF directly from XVG data
    logger.step("Calculating RDF directly from XVG data...")
    rdf_results = calculate_rdf_from_xvg(isg_mask, partition_coords, radial_positions, bins=8)
    
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
    
    logger.step("=" * 60)
    logger.step("Radial RDF analysis completed successfully!")
    logger.step(f"Dataset: {dataid}")
    logger.step(f"Number of radial bins: {len(rdf_results['radii'])}")
    logger.step("=" * 60)


if __name__ == '__main__':
    main()
