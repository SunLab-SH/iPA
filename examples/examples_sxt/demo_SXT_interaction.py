"""
SXT Interaction Analysis Demo

This script demonstrates interaction analysis between organelles (e.g., ISG and Mito)
using SXT data. It calculates contact RDF and contact probability.
"""

import os
import sys
import numpy as np
import tifffile
import pandas as pd
from scipy.ndimage import distance_transform_edt, label

# Add project root to path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
sys.path.insert(0, PROJECT_ROOT)

from ipa.analysis import calculate_contact_rdf, calculate_contact_probability
from ipa.processing.partitioning import Partitioning
from ipa.data_loader import QuickLogger

def main():
    output_dir = os.path.join(PROJECT_ROOT, 'results', 'sxt_interaction_demo')
    os.makedirs(output_dir, exist_ok=True)

    logger = QuickLogger("sxt_interaction", log_dir=os.path.join(PROJECT_ROOT, 'logs'))
    logger.step("Starting SXT Interaction Analysis Demo")

    # --- Configuration ---
    data_id = '784_5'
    data_root = os.path.join(PROJECT_ROOT, 'data', 'sxt')
    
    isg_path = os.path.join(data_root, f'{data_id}_isg_label.tiff')
    mito_path = os.path.join(data_root, f'{data_id}_mito_label.tiff')
    pm_path = os.path.join(data_root, f'{data_id}_wholecell_label.tiff')
    ne_path = os.path.join(data_root, f'{data_id}_NC_label.tiff')

    if not all(os.path.exists(p) for p in [isg_path, mito_path, pm_path, ne_path]):
        logger.step("Error: Required mask files not found.")
        return

    # Load Masks
    isg_mask = tifffile.imread(isg_path).astype(bool)
    mito_mask = tifffile.imread(mito_path).astype(bool)
    pm_mask = tifffile.imread(pm_path)
    ne_mask = tifffile.imread(ne_path)

    logger.step(f"   Loaded full-size data shape: {isg_mask.shape}")
    
    # Use center slice for faster demo (optional) - DISABLED for full test
    USE_CENTER_SLICE = False
    if USE_CENTER_SLICE:
        mid_z = isg_mask.shape[0] // 2
        isg_mask = isg_mask[mid_z-10:mid_z+10]
        mito_mask = mito_mask[mid_z-10:mid_z+10]
        pm_mask = pm_mask[mid_z-10:mid_z+10]
        ne_mask = ne_mask[mid_z-10:mid_z+10]
        logger.step(f"   Using center 20 slices for faster demo: {isg_mask.shape}")

    logger.step("1. Identifying contact points (ISG-Mito distance < 3 voxels)...")
    mito_dist = distance_transform_edt(~mito_mask)
    contact_mask = (isg_mask & (mito_dist <= 3))
    
    contact_coords = np.array(np.where(contact_mask)).T
    logger.step(f"   Found {len(contact_coords)} contact voxels.")

    logger.step("2. Creating spatial partitions...")
    partitioner = Partitioning(root_dir=output_dir, n_slices=8)
    center, ne_edge, pm_edge = partitioner.extract_ne_pm_edges(pm_mask, ne_mask)
    partition_mask = partitioner.create_nepm_radial_partitions(
        ne_edge, pm_edge, shape=pm_mask.shape, n_slices=8, 
        pm_mask=pm_mask, ne_mask=ne_mask
    )

    logger.step("3. Calculating Contact RDF...")
    rdf_results = calculate_contact_rdf(contact_coords, partition_mask, n_shells=8)
    logger.step(f"   Contact RDF values: {rdf_results['rdf']}")

    logger.step("4. Calculating Contact Probability...")
    isg_labeled = label(isg_mask)[0]
    isg_counts_per_shell = np.zeros(8)
    for i in range(1, 9):
        shell_mask = (partition_mask == i)
        unique_isgs = np.unique(isg_labeled[shell_mask])
        isg_counts_per_shell[i-1] = len(unique_isgs[unique_isgs != 0])
        
    contact_probs = calculate_contact_probability(rdf_results['contact_counts'], isg_counts_per_shell)
    logger.step(f"   Contact Probabilities: {contact_probs}")

    # Save Results
    res_df = pd.DataFrame({
        'Shell_ID': range(1, 9),
        'Contact_Count': rdf_results['contact_counts'],
        'Contact_RDF': rdf_results['rdf'],
        'Contact_Probability': contact_probs
    })
    output_csv = os.path.join(output_dir, f'{data_id}_interaction_analysis.csv')
    res_df.to_csv(output_csv, index=False)
    logger.file_out(output_csv)
    logger.step("Demo Completed!")

if __name__ == '__main__':
    main()
