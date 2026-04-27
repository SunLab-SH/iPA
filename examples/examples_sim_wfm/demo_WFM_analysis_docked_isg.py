"""
WFM Docked ISG Analysis Demo

This script demonstrates how to use the iPA library to analyze docked insulin 
secretory granules (ISG) from WFM images.
"""
import numpy as np
import pandas as pd
import os
from parsers import get_args
from ipa.analysis import analyze_docked_granules, extract_isg_features
from ipa.data_loader import UniversalDataLoader
from ipa.data_loader import QuickLogger
def main():
    """Main function to run WFM docked ISG analysis pipeline."""
    # Initialize logger
    mainpath = get_args().main_path
    log_dir = f'{mainpath}/logs'
    logger = QuickLogger("wfm_docked_isg_analysis", log_dir=log_dir)
    logger.step("Starting WFM Docked ISG Analysis Demo")
    
    # Setup data paths
    datapath = f'{mainpath}/data/wfm_images'
    dataid = '20220909_30-2-2-SIM'
    
    # File paths
    isg_tif_path = f"{datapath}/{dataid}_ISG_seg.tif"
    pm_mrc_path = f"{datapath}/{dataid}_PM.mrc"
    csv_path = f"{datapath}/{dataid}_ISG_seg_result3D_generated.csv"
    
    logger.step(f"Working with dataset: {dataid}")
    logger.file_in(isg_tif_path)
    logger.file_in(pm_mrc_path)
    

    # Load image data
    logger.step("Loading image data")
    img_isg = UniversalDataLoader.load_data(isg_tif_path)
    mem_mask = UniversalDataLoader.load_data(pm_mrc_path)
    
    logger.step(f"ISG data loaded: {img_isg.shape}")
    logger.step(f"Membrane mask loaded: {mem_mask.shape}")
    
    # Perform 3D segmentation
    logger.step("Performing 3D segmentation")
    df_results, segmented_labels = extract_isg_features(
        img_isg, 
        min_distance=5, 
        min_size=20, 
        save_csv=True, 
        output_path=csv_path
    )
    
    logger.step(f"Segmented {len(df_results)-1} granule regions")
    logger.file_out(csv_path)
    
    # Read CSV data for analysis
    logger.step("Reading CSV data for analysis")
    isg_data = np.genfromtxt(csv_path, delimiter=',')
    
    # Analyze docked granules
    logger.step("Analyzing docked granules")
    result = analyze_docked_granules(
        segmented_labels, 
        mem_mask, 
        isg_data, 
        show_visualization=True
    )
    
    # Save results
    outputpath = f'{mainpath}/data/results'
    os.makedirs(outputpath, exist_ok=True)
    
    df = pd.DataFrame([result])
    output_file = f"{outputpath}/cell_docking_statistics_{dataid}.csv"
    df.to_csv(output_file, index=False)
    
    logger.file_out(output_file)
    logger.step(f"Docking statistics: {result}")
        

if __name__ == "__main__":
    main()
