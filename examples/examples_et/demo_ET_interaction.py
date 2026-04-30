"""
Cryo-ET Interaction Analysis Demo

This script demonstrates interaction analysis between organelles in Cryo-ET data.
It supports multiple analysis types including actin-vesicle, mito-ER, etc.

Note: Interaction analysis can be computationally intensive and may take several minutes
to complete depending on data size. For quick testing, consider running only one analysis type.
"""

import os
import sys
import json
import time
import numpy as np

# Add project root to path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
sys.path.insert(0, PROJECT_ROOT)

from ipa.data_loader import QuickLogger
from ipa.analysis import (
    actin_to_actin_analysis,
    actin_to_vesicle_analysis,
    actin_to_microtube_analysis,
    actin_to_mito_analysis,
    microtube_to_mito_analysis,
    mito_to_endoreticulum_analysis,
    vesicle_to_mito_analysis,
)

def main():
    # --- Configuration ---
    DATA_ID = '20220326_2.8_2'  # Dataset with complete organelle data (actin, mito, ER, vesicle)
    OUTPUT_DIR = os.path.join(PROJECT_ROOT, 'results', 'et_interaction_demo')
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Cryo-ET voxel size: 1.34 nm (from Manuscript)
    VOXEL_SIZE = [1.34, 1.34, 1.34]
    SHIFT_BIAS = [0, 0, 0]  # No shift needed for this dataset

    logger = QuickLogger("et_interaction", log_dir=os.path.join(PROJECT_ROOT, 'logs'))
    logger.step("Starting Cryo-ET Interaction Analysis Demo")
    logger.step("Note: This analysis may take several minutes to complete...")

    # Build paths
    data_root = os.path.join(PROJECT_ROOT, 'data', 'cryoET', DATA_ID)
    paths = {
        "vesicle_mask": os.path.join(data_root, f"{DATA_ID}_vesicle.mrc"),
        "actin_coords": os.path.join(data_root, f"{DATA_ID}_actin_filled_points.json"),
        "microtube_coords": os.path.join(data_root, f"{DATA_ID}_microtube_filled_points.json"),
        "mito_mask": os.path.join(data_root, f"{DATA_ID}_mito.mrc"),
        "er_mask": os.path.join(data_root, f"{DATA_ID}_endoreticulum.mrc"),
    }

    config = {
        "voxel_size_xyz": VOXEL_SIZE,
        "shift_bias": SHIFT_BIAS,
        "visualization": False,
        "output_dir": OUTPUT_DIR,
        "zoom_degree": 0.1,
        "batch_size": 512,
    }

    def run_analysis(name, func, **kwargs):
        t0 = time.time()
        try:
            logger.step(f"Running {name}...")
            res = func(**kwargs)
            dt = time.time() - t0
            output_path = os.path.join(OUTPUT_DIR, f"{name}.json")
            with open(output_path, "w") as f:
                json.dump(res, f, indent=2)
            logger.step(f"   Completed in {dt:.2f}s -> {output_path}")
        except Exception as e:
            logger.error(f"   Failed: {e}")

    # 1. Actin-to-Vesicle Interaction
    if all(os.path.exists(p) for p in [paths["actin_coords"], paths["vesicle_mask"]]):
        run_analysis("actin_to_vesicle", actin_to_vesicle_analysis,
                     data_id=DATA_ID, mask_file=paths["vesicle_mask"],
                     json_file=paths["actin_coords"], config=config)

    # 2. Mito-to-ER Interaction
    if all(os.path.exists(p) for p in [paths["mito_mask"], paths["er_mask"]]):
        run_analysis("mito_to_er", mito_to_endoreticulum_analysis,
                     data_id=DATA_ID, mito_file=paths["mito_mask"],
                     er_file=paths["er_mask"], config=config)

    # 3. Actin-to-Mito Interaction
    if all(os.path.exists(p) for p in [paths["actin_coords"], paths["mito_mask"]]):
        run_analysis("actin_to_mito", actin_to_mito_analysis,
                     data_id=DATA_ID, actin_file=paths["actin_coords"],
                     mito_file=paths["mito_mask"], config=config)

    logger.step(f"Analysis completed. Results saved to: {OUTPUT_DIR}")

if __name__ == '__main__':
    main()
