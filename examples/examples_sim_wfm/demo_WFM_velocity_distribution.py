import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')  # Set non-GUI backend before importing pyplot
import matplotlib.pyplot as plt
import os
from ipa.analysis.dynamic import create_velocity_violin_plot
from ipa.data_loader import QuickLogger
from parsers import get_args

mainpath = get_args().main_path

# Set global style parameters
plt.rcParams.update({
    'font.size': 8,
    'font.sans-serif': 'Arial',
    'axes.unicode_minus': False,
    'axes.spines.top': False,
    'axes.spines.right': False
})

def load_velocity_data(valid_sequences, data_folder=f'{mainpath}/data/wfm_images', logger=None):
    """
    Load velocity data from .npy files for multiple sequences
    
    Parameters:
    -----------
    valid_sequences : list
        List of sequence names to process
    data_folder : str, optional
        Path to folder containing velocity data files
    logger : QuickLogger, optional
        Logger instance for recording operations
        
    Returns:
    --------
    list : List of tuples (speed, category) for all processed sequences
    """
    speed_data = []
    
    for seq in valid_sequences:
        # Determine category based on sequence name prefix
        category = 'LG' if seq.startswith('2.8') else 'HG'
        
        try:
            # Load velocity data file
            file_path = os.path.join(data_folder, f"{seq}_Velocity_3dvelocity.npy")
            if logger:
                logger.file_in(file_path)
            data = np.load(file_path)
            
            # Process velocity data if it has sufficient dimensions
            if data.shape[0] >= 3:
                # Calculate velocity vectors between consecutive time points
                velocity_vectors = np.diff(data[:3, :], axis=1)
                # Calculate 3D speed (magnitude of velocity vectors)
                speeds = np.linalg.norm(velocity_vectors, axis=0)
                # Add speeds with category labels to data collection
                speed_data.extend([(s, category) for s in speeds])
                
            if logger:
                logger.step(f"Successfully processed {seq}: {len(speeds)} velocity points")
            print(f"Successfully processed {seq}: {len(speeds)} velocity points")
            
        except Exception as e:
            error_msg = f"Error processing {seq}: {str(e)}"
            if logger:
                logger.step(error_msg)
            print(error_msg)
    
    return speed_data

def main():
    """
    Main function to orchestrate velocity data analysis and visualization
    """
    # Initialize logger
    log_dir = f'{mainpath}/logs'
    logger = QuickLogger("wfm_velocity_distribution_analysis", log_dir=log_dir)
    logger.step("Starting WFM Velocity Distribution Analysis Demo")
    
    # Define valid sequences for analysis
    valid_sequences = [
        '16.7-30_P21', 
    ]
    
    logger.step("Starting velocity data analysis...")
    logger.step(f"Processing {len(valid_sequences)} sequences")
    print("Starting velocity data analysis...")
    print(f"Processing {len(valid_sequences)} sequences")
    
    # Load velocity data from all sequences
    logger.step("Loading velocity data from sequences")
    speed_data = load_velocity_data(valid_sequences, logger=logger)
    
    if not speed_data:
        error_msg = "No velocity data loaded. Please check file paths and data availability."
        logger.step(error_msg)
        print(error_msg)
        return
    
    logger.step(f"Total velocity data points collected: {len(speed_data)}")
    print(f"Total velocity data points collected: {len(speed_data)}")
    
    # Create and display violin plot
    logger.step("Creating velocity violin plot")
    fig = create_velocity_violin_plot(
        speed_data=speed_data,
        max_speed=30,
        y_limit=0.35,
        figure_size=(1.5, 1.2)
    )
    plt.show()
    
    # Save output
    outputpath = mainpath + '/results/data'
    if not os.path.exists(outputpath):
        os.makedirs(outputpath)
    
    output_file = '{}/velocity_violin_plot.png'.format(outputpath)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    
    logger.file_out(output_file)
    logger.step("Velocity violin plot saved as 'velocity_violin_plot.png'")
    
    # Save the plot instead of showing it to avoid Qt issues
    print("Velocity violin plot saved as 'velocity_violin_plot.png'")

if __name__ == "__main__":
    main()


