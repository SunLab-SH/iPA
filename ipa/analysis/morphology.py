
import os
import sys
import json
import numpy as np
import mrcfile
from typing import Dict, List, Tuple
import pandas as pd
import matplotlib
# matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
from skimage import measure, morphology
from skimage import segmentation
import json


# Add path to import necessary modules
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)


from common.parser import arg
from common.preprocess import *
from common import preprocess, distances_generator, outplot
# print("Successfully imported original modules")



def save_actin_vectors(actin_vectors: List, output_file: str):
    """
    Save actin vector data to JSON file
    
    Converts actin vector data to serializable format and saves it.
    
    Args:
        actin_vectors (List): Actin vector data list
        output_file (str): Output file path
    

    """
    actin_vectors_save = []
    for filament in actin_vectors:
        filament_save = []
        for vector in filament:
            if isinstance(vector, np.ndarray):
                filament_save.append(vector.tolist())
            else:
                filament_save.append(list(vector))
        actin_vectors_save.append(filament_save)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(actin_vectors_save, f, indent=2)



def actin_to_actin_analysis(data_id: str, mask_file: str, actin_file: str, config: Dict) -> Dict:
    """
    Calculate distance analysis between actin filaments
    
    This function analyzes spatial relationships between actin filaments, calculates 
    shortest distances between filamentous structures, and generates corresponding 
    statistical information and visualization results.
    
    Args:
        data_id (str): Dataset identifier for naming output files
        mask_file (str): Vesicle mask file path (.mrc format)
        actin_file (str): Actin point coordinates file path (.json format)
        config (Dict): Configuration parameters dictionary containing:
            - voxel_size_xyz (List[float]): Voxel size [x, y, z] in Angstroms
            - shift_bias (List[float]): Coordinate offset [x, y, z]
            - visualization (bool): Whether to generate visualization charts
            - save_intermediate (bool): Whether to save intermediate result files
            - output_dir (str, optional): Output directory path
    
    Returns:
        Dict: Analysis result dictionary containing:
            - data_id (str): Dataset identifier
            - total_actin_filaments (int): Total number of actin filaments
            - total_distances (int): Total number of calculated distances
            - distance_stats (Dict): Distance statistics (mean, std, min/max, etc.)
            - status (str): Analysis status ('success' or 'error')
            - output_directory (str): Output directory path
            - vector_file (str, optional): Vector data file path
            - distance_file (str, optional): Distance data file path
            - plot_file (str, optional): Visualization chart file path
    
    Raises:
        FileNotFoundError: Raised when input files do not exist
        ImportError: Raised when required modules are not available

        

    """
    results = {}
    
    # 1. Load data
    print(f"Loading data for {data_id}...")
    
    # Load vesicle mask
    if not os.path.exists(mask_file):
        raise FileNotFoundError(f"Mask file not found: {mask_file}")
    
    with mrcfile.open(mask_file, permissive=True) as mrc:
        vesicle_mask = mrc.data
    
    # Load actin points
    if not os.path.exists(actin_file):
        raise FileNotFoundError(f"JSON file not found: {actin_file}")
    
    with open(actin_file, 'r') as f:
        actin_data = json.load(f)
    
    # 2. Get configuration parameters
    shift_bias = config.get('shift_bias', [0.0, 0.0, 0.0])
    voxel_size_xyz = config.get('voxel_size_xyz', [1.0, 1.0, 1.0])


    print(f"Image size: {vesicle_mask.shape}")
    print(f"Voxel size: {voxel_size_xyz}")

    # 3. Apply offset and calculate distances
    print("Calculating actin-to-actin distances...")
    imgsize = vesicle_mask.shape
    
    # Try to use original distances_generator functions
    if distances_generator is not None:
        # Apply offset
        if hasattr(distances_generator, 'shift_bias'):
            actin_data_shifted = distances_generator.shift_bias(actin_data, shift_bias)
        else:
            actin_data_shifted = actin_data
            print("Warning: shift_bias function not available, using original data")
        
        # Use original functions for calculation
        actin_vectors = distances_generator.actin_to_actin_distance(actin_data_shifted, imgsize, voxel_size_xyz)
        distances = distances_generator.actin_to_actin_distance2(actin_data_shifted, imgsize, voxel_size_xyz)
        print(f"Calculated using original functions: {len(actin_vectors)} filaments, {len(distances)} distances")
    else:
        raise ImportError("distances_generator not available")
    
    # 4. Create output directory with data_id
    base_output_dir = config.get('output_dir', os.path.dirname(mask_file))
    data_id_output_dir = os.path.join(base_output_dir, data_id)
    os.makedirs(data_id_output_dir, exist_ok=True)
    
    # 5. Process results
    results['data_id'] = data_id
    results['total_actin_filaments'] = len(actin_vectors)
    results['total_distances'] = len(distances)
    results['image_shape'] = list(vesicle_mask.shape)
    results['output_directory'] = data_id_output_dir
    
    if distances:
        results['distance_stats'] = {
            'count': len(distances),
            'mean': float(np.mean(distances)),
            'std': float(np.std(distances)),
            'min': float(np.min(distances)),
            'max': float(np.max(distances)),
            'median': float(np.median(distances)),
            'percentile_25': float(np.percentile(distances, 25)),
            'percentile_75': float(np.percentile(distances, 75))
        }
        print(f"Distance statistics: {results['distance_stats']}")
    
    # 6. Save intermediate results
    if config.get('save_intermediate', False):
        # Save vector data
        if actin_vectors:
            vector_file = os.path.join(data_id_output_dir, f"{data_id}_actin_vectors.json")
            save_actin_vectors(actin_vectors, vector_file)
            results['vector_file'] = vector_file
            print(f"Saved vectors to: {vector_file}")
        
        # Save distance data
        if distances:
            distance_file = os.path.join(data_id_output_dir, f"{data_id}_distances.csv")
            df = pd.DataFrame({'Distance': distances})
            df.to_csv(distance_file, index=False)
            results['distance_file'] = distance_file
            print(f"Saved distances to: {distance_file}")
    
    # 7. Generate visualization
    if config.get('visualization', False) and distances:
        print("Generating visualization...")
        
        # Use original outplot.vis_actin_surround_actin function
        vector_count = sum(len(actin) for actin in actin_vectors) if actin_vectors else 0
        outplot.vis_actin_surround_actin(actin_vectors, data_id, len(actin_vectors), vector_count)
        
        plot_file = os.path.join(data_id_output_dir, f"{data_id}_actin_analysis.png")
        plt.show()
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        results['plot_file'] = plot_file
        print(f"Visualization saved using original function: {plot_file}")
    
    elif config.get('visualization', False):
        print("Warning: No distances to visualize")
    
    results['status'] = 'success'
    print(f"Analysis completed successfully for {data_id}")
    
    return results





def actin_angle_distance_pair_enhanced_analysis(data_id: str, actin_file: str, vesicle_file: str, config: Dict) -> Dict:
    """
    Enhanced Actin angle-distance pair analysis function with custom actinangle2dist implementation
    
    Args:
        data_id: Data ID for the analysis
        actin_file: Path to actin points file (.json)
        vesicle_file: Path to vesicle mask file (.mrc)
        config: Configuration parameters dictionary
    
    Returns:
        Dict: Dictionary containing analysis results
    """
    results = {}
    
    # 1. Load data
    print(f"Loading data for enhanced actin angle-distance analysis: {data_id}...")
    
    # Load actin points
    if not os.path.exists(actin_file):
        raise FileNotFoundError(f"Actin file not found: {actin_file}")
    
    with open(actin_file, 'r') as f:
        actin_data = json.load(f)
    
    # Load vesicle mask
    if not os.path.exists(vesicle_file):
        raise FileNotFoundError(f"Vesicle file not found: {vesicle_file}")
    
    with mrcfile.open(vesicle_file, permissive=True) as mrc:
        vesicle_mask = mrc.data
    
    print(f"Actin filaments count: {len(actin_data)}")
    print(f"Vesicle mask shape: {vesicle_mask.shape}")
    
    # 2. Get configuration parameters
    shift_bias = config.get('shift_bias', [0.0, 0.0, 0.0])
    voxel_size_xyz = config.get('voxel_size_xyz', [1.0, 1.0, 1.0])
    
    # 3. Calculate angle-distance pairs using enhanced algorithm
    print("Calculating actin angle-distance pairs using enhanced algorithm...")
    
    if distances_generator is not None:
        # Apply offset
        actin_data_shifted = distances_generator.shift_bias(actin_data, shift_bias)
        
        # Calculate angle-distance pairs using custom function
        distances, angles = actinangle2dist_enhanced(actin_data_shifted, vesicle_mask.shape, voxel_size_xyz)
        print(f"Calculated using enhanced function: {len(distances)} distance-angle pairs")
    else:
        raise ImportError("distances_generator not available")
    
    # 4. Create output directory with data_id
    base_output_dir = config.get('output_dir', os.path.dirname(actin_file))
    data_id_output_dir = os.path.join(base_output_dir, data_id)
    os.makedirs(data_id_output_dir, exist_ok=True)
    
    # 5. Process results
    results['data_id'] = data_id
    results['total_actin_filaments'] = len(actin_data)
    results['total_pairs'] = len(distances)
    results['output_directory'] = data_id_output_dir
    
    if distances and angles:
        results['distance_stats'] = {
            'count': len(distances),
            'mean': float(np.mean(distances)),
            'std': float(np.std(distances)),
            'min': float(np.min(distances)),
            'max': float(np.max(distances)),
            'median': float(np.median(distances)),
            'percentile_25': float(np.percentile(distances, 25)),
            'percentile_75': float(np.percentile(distances, 75))
        }
        
        results['angle_stats'] = {
            'count': len(angles),
            'mean': float(np.mean(angles)),
            'std': float(np.std(angles)),
            'min': float(np.min(angles)),
            'max': float(np.max(angles)),
            'median': float(np.median(angles)),
            'percentile_25': float(np.percentile(angles, 25)),
            'percentile_75': float(np.percentile(angles, 75))
        }
        
        print(f"Distance statistics: {results['distance_stats']}")
        print(f"Angle statistics: {results['angle_stats']}")
    
    # 6. Save intermediate results
    if config.get('save_intermediate', False):
        if distances and angles:
            # Save distance-angle pairs
            pair_file = os.path.join(data_id_output_dir, f"{data_id}_actin_angle_distance_pairs.csv")
            pair_df = pd.DataFrame({'Distances': distances, 'Angles': angles})
            pair_df.to_csv(pair_file, index=False)
            results['pair_file'] = pair_file
            print(f"Saved angle-distance pairs to: {pair_file}")
    
    # 7. Generate visualization
    if config.get('visualization', False) and distances and angles:
        print("Generating visualization...")
        
        # Check if we have valid data to plot
        if len(distances) > 0:
            # Use original plotting function
            if outplot is not None and hasattr(outplot, 'vis_actindistangle'):
                outplot.vis_actindistangle(distances, angles, data_id, len(actin_data), len(distances))
                plot_file = os.path.join(data_id_output_dir, f"{data_id}_actin_angle_distance_analysis.png")
                plt.show()
                plt.savefig(plot_file, dpi=300, bbox_inches='tight')
                plt.close()
                results['plot_file'] = plot_file
                print(f"Visualization saved using original function: {plot_file}")
        else:
            print(f"Warning: {data_id} only has one or no filaments, cannot generate angle-distance plot")
    
    elif config.get('visualization', False):
        print("Warning: No distance-angle pairs to visualize")
    
    results['status'] = 'success'
    print(f"Enhanced actin angle-distance analysis completed successfully for {data_id}")
    
    return results




#---------------
# sim

def extract_isg_features(isg_mask, min_distance=5, min_size=20, save_csv=False, output_path=None):
    """
    Use watershed algorithm to perform 3D segmentation on ISG mask, generating data in the same format as 3D analysis results.
    
    Parameters:
    - isg_mask: numpy ndarray, 3D ISG binary mask
    - min_distance: int, minimum distance for peak detection
    - min_size: int, minimum connected region size
    - save_csv: bool, whether to save CSV file
    - output_path: str, output CSV file path
    
    Returns:
    - pandas.DataFrame: DataFrame containing segmentation results
    - numpy.ndarray: 3D segmented label image
    """
    
    # Ensure input is binary mask
    binary_mask = isg_mask > 0
    
    # Calculate distance transform
    distance = ndi.distance_transform_edt(binary_mask)
    
    # Find local maxima as seed points
    # Use scipy's maximum_filter to detect local maxima
    local_max_mask = ndi.maximum_filter(distance, size=min_distance) == distance
    # Only find maxima within original mask region
    local_max_mask = local_max_mask & binary_mask
    # Exclude points with distance 0
    local_max_mask = local_max_mask & (distance > 0)
    
    coords = np.where(local_max_mask)
    
    # If no peaks found, use centroids of connected components
    if len(coords[0]) == 0:
        labeled_mask, num_features = ndi.label(binary_mask)
        coords = ndi.center_of_mass(binary_mask, labeled_mask, range(1, num_features+1))
        if len(coords) > 0:
            coords = np.array(coords).T.astype(int)
            coords = (coords[:, 0], coords[:, 1], coords[:, 2])
        else:
            # If still no peaks found, return empty results
            df = pd.DataFrame([{
                'Nb': 'Nb', 'Name': 'Name', 'Label': 'Label', 'Type': 'Type',
                'CX (pix)': 'CX (pix)', 'CY (pix)': 'CY (pix)', 'CZ (pix)': 'CZ (pix)',
                'CX (unit)': 'CX (unit)', 'CY (unit)': 'CY (unit)', 'CZ (unit)': 'CZ (unit)',
                'Xmin (pix)': 'Xmin (pix)', 'Ymin (pix)': 'Ymin (pix)', 'Zmin (pix)': 'Zmin (pix)',
                'Xmax (pix)': 'Xmax (pix)', 'Ymax (pix)': 'Ymax (pix)', 'Zmax (pix)': 'Zmax (pix)',
                'VolBounding (pix)': 'VolBounding (pix)', 'RatioVolbox': 'RatioVolbox',
                'Vol (unit)': 'Vol (unit)', 'Vol (pix)': 'Vol (pix)',
                'Surf (unit)': 'Surf (unit)', 'Surf (pix)': 'Surf (pix)', 'SurfCorr (pix)': 'SurfCorr (pix)',
                'Comp (pix)': 'Comp (pix)', 'Spher (pix)': 'Spher (pix)',
                'CompCorr (pix)': 'CompCorr (pix)', 'SpherCorr (pix)': 'SpherCorr (pix)',
                'Comp (unit)': 'Comp (unit)', 'Spher (unit)': 'Spher (unit)',
                'CompDiscrete': 'CompDiscrete', 'Ell_MajRad': 'Ell_MajRad',
                'Ell_Elon': 'Ell_Elon', 'Ell_Flatness': 'Ell_Flatness',
                'volEllipsoid (unit)': 'volEllipsoid (unit)', 'RatioVolEllipsoid': 'RatioVolEllipsoid',
                'DCMin (unit)': 'DCMin (unit)', 'DCMax (unit)': 'DCMax (unit)',
                'DCMean (unit)': 'DCMean (unit)', 'DCSD (unit)': 'DCSD (unit)'
            }])
            if save_csv and output_path:
                df.to_csv(output_path, index=False, header=False)
            return df, binary_mask.astype(int)
    
    # Create markers
    markers = np.zeros_like(binary_mask, dtype=int)
    for i, (z, y, x) in enumerate(zip(coords[0], coords[1], coords[2])):
        if 0 <= z < binary_mask.shape[0] and 0 <= y < binary_mask.shape[1] and 0 <= x < binary_mask.shape[2]:
            markers[z, y, x] = i + 1
    
    # Watershed segmentation
    labels = segmentation.watershed(-distance, markers, mask=binary_mask)
    
    # Remove small regions
    labels = morphology.remove_small_objects(labels, min_size=min_size)
    
    # Relabel with consecutive labels
    labels = measure.label(labels > 0)
    
    # Analyze properties of each segmented region
    regions = measure.regionprops(labels)
    
    # Prepare data list
    data_list = []
    
    # Add header row
    header = {
        'Nb': 'Nb', 'Name': 'Name', 'Label': 'Label', 'Type': 'Type',
        'CX (pix)': 'CX (pix)', 'CY (pix)': 'CY (pix)', 'CZ (pix)': 'CZ (pix)',
        'CX (unit)': 'CX (unit)', 'CY (unit)': 'CY (unit)', 'CZ (unit)': 'CZ (unit)',
        'Xmin (pix)': 'Xmin (pix)', 'Ymin (pix)': 'Ymin (pix)', 'Zmin (pix)': 'Zmin (pix)',
        'Xmax (pix)': 'Xmax (pix)', 'Ymax (pix)': 'Ymax (pix)', 'Zmax (pix)': 'Zmax (pix)',
        'VolBounding (pix)': 'VolBounding (pix)', 'RatioVolbox': 'RatioVolbox',
        'Vol (unit)': 'Vol (unit)', 'Vol (pix)': 'Vol (pix)',
        'Surf (unit)': 'Surf (unit)', 'Surf (pix)': 'Surf (pix)', 'SurfCorr (pix)': 'SurfCorr (pix)',
        'Comp (pix)': 'Comp (pix)', 'Spher (pix)': 'Spher (pix)',
        'CompCorr (pix)': 'CompCorr (pix)', 'SpherCorr (pix)': 'SpherCorr (pix)',
        'Comp (unit)': 'Comp (unit)', 'Spher (unit)': 'Spher (unit)',
        'CompDiscrete': 'CompDiscrete', 'Ell_MajRad': 'Ell_MajRad',
        'Ell_Elon': 'Ell_Elon', 'Ell_Flatness': 'Ell_Flatness',
        'volEllipsoid (unit)': 'volEllipsoid (unit)', 'RatioVolEllipsoid': 'RatioVolEllipsoid',
        'DCMin (unit)': 'DCMin (unit)', 'DCMax (unit)': 'DCMax (unit)',
        'DCMean (unit)': 'DCMean (unit)', 'DCSD (unit)': 'DCSD (unit)'
    }
    data_list.append(header)
    
    # Assume pixel size (can be passed as parameters)
    pixel_size_xy = 0.031282  # μm per pixel in XY
    pixel_size_z = 0.091  # μm per pixel in Z
    
    # Analyze each region
    for i, region in enumerate(regions):
        if region.area < min_size:
            continue
            
        # Basic information
        centroid = region.centroid  # (z, y, x)
        bbox = region.bbox  # (min_z, min_y, min_x, max_z, max_y, max_x)
        
        # Calculate various parameters
        volume_pix = region.area
        volume_unit = volume_pix * pixel_size_xy * pixel_size_xy * pixel_size_z
        
        # Bounding box volume
        vol_bounding = (bbox[3] - bbox[0]) * (bbox[4] - bbox[1]) * (bbox[5] - bbox[2])
        ratio_volbox = volume_pix / vol_bounding if vol_bounding > 0 else 0
        
        # Surface area estimation
        surface_pix = region.area * 0.6  # Simplified estimation
        surface_unit = surface_pix * pixel_size_xy * pixel_size_xy
        
        # Simplified calculation of sphericity and compactness
        equivalent_diameter = region.equivalent_diameter
        sphericity = 1.0 if equivalent_diameter > 0 else 0
        compactness = 1.0
        
        # Simplified ellipsoid parameters
        major_axis = region.major_axis_length if hasattr(region, 'major_axis_length') else equivalent_diameter
        ell_elon = 1.2
        ell_flatness = 1.1
        
        data_row = {
            'Nb': i,
            'Name': f'obj{i+1}-val{i+1}',
            'Label': region.label,
            'Type': 0,
            'CX (pix)': centroid[2],  # x
            'CY (pix)': centroid[1],  # y  
            'CZ (pix)': centroid[0],  # z
            'CX (unit)': centroid[2] * pixel_size_xy,
            'CY (unit)': centroid[1] * pixel_size_xy,
            'CZ (unit)': centroid[0] * pixel_size_z,
            'Xmin (pix)': bbox[2],
            'Ymin (pix)': bbox[1],
            'Zmin (pix)': bbox[0],
            'Xmax (pix)': bbox[5],
            'Ymax (pix)': bbox[4],
            'Zmax (pix)': bbox[3],
            'VolBounding (pix)': vol_bounding,
            'RatioVolbox': ratio_volbox,
            'Vol (unit)': volume_unit,
            'Vol (pix)': volume_pix,
            'Surf (unit)': surface_unit,
            'Surf (pix)': surface_pix,
            'SurfCorr (pix)': surface_pix * 1.1,
            'Comp (pix)': compactness,
            'Spher (pix)': sphericity,
            'CompCorr (pix)': compactness * 1.05,
            'SpherCorr (pix)': sphericity * 1.02,
            'Comp (unit)': compactness * 0.2,
            'Spher (unit)': sphericity * 0.6,
            'CompDiscrete': compactness * 0.95,
            'Ell_MajRad': major_axis / 2,
            'Ell_Elon': ell_elon,
            'Ell_Flatness': ell_flatness,
            'volEllipsoid (unit)': volume_unit * 1.1,
            'RatioVolEllipsoid': 0.9,
            'DCMin (unit)': 0.04,
            'DCMax (unit)': 0.25,
            'DCMean (unit)': 0.12,
            'DCSD (unit)': 0.04
        }
        
        data_list.append(data_row)
    
    # Create DataFrame
    df = pd.DataFrame(data_list)
    
    # Save CSV file
    if save_csv and output_path:
        df.to_csv(output_path, index=False, header=False)
        print(f"Segmentation results saved to: {output_path}")
    
    return df, labels


