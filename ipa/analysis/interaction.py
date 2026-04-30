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
from scipy.spatial.distance import cdist
import scipy.ndimage
import json

# Standard package imports replacing sys.path.append
from .common.parser import arg
from .common import preprocess, distances_generator, outplot

def calculate_contact_rdf(contact_coords, partition_mask, n_shells=8):
    """
    Calculate the Radial Distribution Function (RDF) of organelle contacts.
    
    Parameters
    ----------
    contact_coords : list of tuples
        List of (z, y, x) coordinates for contact points.
    partition_mask : np.ndarray
        3D mask where values 1-n_shells represent radial shells.
    n_shells : int
        Number of radial shells.
        
    Returns
    -------
    dict
        Dictionary containing contact counts per shell and normalized RDF.
    """
    contact_counts = np.zeros(n_shells)
    total_voxels_per_shell = np.zeros(n_shells)
    
    for i in range(1, n_shells + 1):
        total_voxels_per_shell[i-1] = np.sum(partition_mask == i)
        
    for coord in contact_coords:
        z, y, x = int(coord[0]), int(coord[1]), int(coord[2])
        if 0 <= z < partition_mask.shape[0] and 0 <= y < partition_mask.shape[1] and 0 <= x < partition_mask.shape[2]:
            shell_id = partition_mask[z, y, x]
            if 1 <= shell_id <= n_shells:
                contact_counts[shell_id - 1] += 1
                
    # Normalize by shell volume to get density
    densities = np.divide(contact_counts, total_voxels_per_shell, out=np.zeros_like(contact_counts), where=total_voxels_per_shell!=0)
    avg_density = np.sum(contact_counts) / np.sum(total_voxels_per_shell)
    
    rdf = np.divide(densities, avg_density, out=np.zeros_like(densities), where=avg_density!=0)
    
    return {
        'contact_counts': contact_counts,
        'shell_volumes': total_voxels_per_shell,
        'rdf': rdf
    }

def calculate_contact_probability(contact_counts, total_organelles_per_shell):
    """
    Calculate the probability that an organelle participates in a contact within each shell.
    
    Parameters
    ----------
    contact_counts : np.ndarray
        Number of contacts in each shell.
    total_organelles_per_shell : np.ndarray
        Total number of organelles in each shell.
        
    Returns
    -------
    np.ndarray
        Probability of contact per shell.
    """
    probs = np.divide(contact_counts, total_organelles_per_shell, out=np.zeros_like(contact_counts, dtype=float), where=total_organelles_per_shell!=0)
    return probs
# print("Successfully imported original modules")








def actin_to_vesicle_analysis(data_id: str, mask_file: str, json_file: str, config: Dict) -> Dict:
    """
    Calculate distance analysis from actin filaments to vesicle membranes
    
    This function analyzes shortest distances from actin filaments to vesicle 
    membrane surfaces, providing quantitative analysis of cytoskeleton-membrane 
    structure interactions.
    
    Args:
        data_id (str): Dataset identifier
        mask_file (str): Vesicle mask file path (.mrc format)
        json_file (str): Actin point coordinates file path (.json format)
        config (Dict): Configuration parameters dictionary, same as actin_to_actin_analysis
    
    Returns:
        Dict: Analysis result dictionary containing actin-to-vesicle membrane distance statistics
    
    Examples:
        >>> result = actin_to_vesicle_analysis(
        ...     data_id="sample_01",
        ...     mask_file="vesicle_mask.mrc",
        ...     json_file="actin_coords.json",
        ...     config=config
        ... )
        >>> print(f"Average distance: {result['distance_stats']['mean']:.2f}")
    """
    results = {}
    
    # 1. Load data
    print(f"Loading data for actin-to-vesicle analysis: {data_id}...")
    
    # Load vesicle mask
    if not os.path.exists(mask_file):
        raise FileNotFoundError(f"Mask file not found: {mask_file}")
    
    with mrcfile.open(mask_file, permissive=True) as mrc:
        vesicle_mask = mrc.data
    
    # Load actin points
    if not os.path.exists(json_file):
        raise FileNotFoundError(f"JSON file not found: {json_file}")
    
    with open(json_file, 'r') as f:
        actin_data = json.load(f)
    

    # 2. Get configuration parameters
    shift_bias = config.get('shift_bias', [0.0, 0.0, 0.0])
    voxel_size_xyz = config.get('voxel_size_xyz', [1.0, 1.0, 1.0])
    
    print(f"Image size: {vesicle_mask.shape}")
    print(f"Voxel size: {voxel_size_xyz}")
    print(f"Shift bias: {shift_bias}")
    print(f"Actin filaments count: {len(actin_data)}")
    
    # 3. Apply offset and calculate distances
    print("Calculating actin-to-vesicle distances...")
    
    # Apply offset using original function
    actin_data_shifted = distances_generator.shift_bias(actin_data, shift_bias)
    
    # Calculate actin to vesicle membrane distance using original function
    distances = distances_generator.coords_to_mem_distance_generator(
        actin_data_shifted, vesicle_mask, voxel_size_xyz
    )
    print(f"Calculated successfully using original function: {len(distances)} actin-to-vesicle distances")

    
    # 4. Process results
    results['data_id'] = data_id
    results['total_actin_filaments'] = len(actin_data)
    results['total_distances'] = len(distances)
    results['image_shape'] = list(vesicle_mask.shape)
    
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
    
    # 5. Save intermediate results
    if config.get('save_intermediate', False):
        output_dir = os.path.dirname(mask_file)
        
        # Save distance data
        if distances:
            distance_file = os.path.join(output_dir, f"{data_id}_actin_to_vesicle_distances.csv")
            df = pd.DataFrame({'Distance': distances})
            df.to_csv(distance_file, index=False)
            results['distance_file'] = distance_file
            print(f"Saved distances to: {distance_file}")
    
    # 6. Generate visualization
    if config.get('visualization', False) and distances:
        print("Generating visualization...")

        # Use original outplot.actin_to_vesicle_dist_hist function
        dist_df = pd.DataFrame({'Distance': distances})
        outplot.actin_to_vesicle_dist_hist(dist_df['Distance'].to_numpy(), data_id, len(actin_data), len(distances))
        
        plot_file = os.path.join(config.get('output_dir', os.path.dirname(mask_file)), f"{data_id}_actin_to_vesicle_analysis.png")
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        results['plot_file'] = plot_file
        print(f"Generated chart using original visualization function: {plot_file}")

    results['status'] = 'success'
    print(f"Actin-to-vesicle analysis completed successfully for {data_id}")
    
    return results



def actin_to_vesicle_dist_angle_pair_analysis(data_id: str, mask_file: str, json_file: str, angle_file: str, config: Dict) -> Dict:
    """
    Calculate distance-angle pair analysis from actin filaments to vesicle membranes
    
    This function simultaneously analyzes distance and angle relationships between 
    actin filaments and vesicle membranes, providing comprehensive spatial 
    interaction information.
    
    Args:
        data_id (str): Dataset identifier
        mask_file (str): Vesicle mask file path (.mrc format)
        json_file (str): Actin point coordinates file path (.json format)
        angle_file (str): Actin angle data file path (.json format)
        config (Dict): Configuration parameters dictionary
    
    Returns:
        Dict: Analysis result dictionary containing distance and angle statistics
            - distance_stats (Dict): Distance statistics
            - angle_stats (Dict): Angle statistics (in degrees)
    
    Note:
        Angle calculation based on angle between actin filament direction vectors 
        and vesicle membrane normal vectors
    """
    results = {}
    
    # 1. Load data
    print(f"Loading data for actin-to-vesicle distance-angle analysis: {data_id}...")
    
    # Load vesicle mask
    if not os.path.exists(mask_file):
        raise FileNotFoundError(f"Mask file not found: {mask_file}")
    
    with mrcfile.open(mask_file, permissive=True) as mrc:
        vesicle_mask = mrc.data
    
    # Load actin points
    if not os.path.exists(json_file):
        raise FileNotFoundError(f"JSON file not found: {json_file}")
    
    with open(json_file, 'r') as f:
        actin_data = json.load(f)
    
    # Load actin angles
    if not os.path.exists(angle_file):
        raise FileNotFoundError(f"Angle file not found: {angle_file}")
    
    with open(angle_file, 'r') as f:
        actin_angle_data = json.load(f)
    
    # 2. Get configuration parameters
    shift_bias = config.get('shift_bias', [0.0, 0.0, 0.0])
    voxel_size_xyz = config.get('voxel_size_xyz', [1.0, 1.0, 1.0])
    
    print(f"Image size: {vesicle_mask.shape}")
    print(f"Voxel size: {voxel_size_xyz}")
    print(f"Actin filaments count: {len(actin_data)}")
    print(f"Actin angle data loaded: {len(actin_angle_data)} entries")
    
    # 3. Apply offset and calculate distance-angle pairs
    print("Calculating actin-to-vesicle distance-angle pairs...")
    
    if distances_generator is not None:
        if hasattr(distances_generator, 'shift_bias'):
            actin_data_shifted = distances_generator.shift_bias(actin_data, shift_bias)
        else:
            actin_data_shifted = actin_data
            print("Warning: shift_bias function not available, using original data")
        
        # Use original function to calculate actin to vesicle membrane distance and angle pairs
        distances, angles = distances_generator.coords_to_mem_distance_angle_pair_generator(
            actin_data_shifted, actin_angle_data, vesicle_mask, voxel_size_xyz
        )
        print(f"Calculated using original function: {len(distances)} distance-angle pairs")
    else:
        raise ImportError("distances_generator not available")
    
    # 4. Create output directory with data_id
    base_output_dir = config.get('output_dir', os.path.dirname(mask_file))
    data_id_output_dir = os.path.join(base_output_dir, data_id)
    os.makedirs(data_id_output_dir, exist_ok=True)
    
    # 5. Process results
    results['data_id'] = data_id
    results['total_actin_filaments'] = len(actin_data)
    results['total_pairs'] = len(distances)
    results['output_directory'] = data_id_output_dir
    
    if distances and angles:
        # Convert to numpy arrays and ensure numeric types
        distances_array = np.array(distances, dtype=float)
        angles_array = np.array(angles, dtype=float)
        
        results['distance_stats'] = {
            'count': len(distances_array),
            'mean': float(np.mean(distances_array)),
            'std': float(np.std(distances_array)),
            'min': float(np.min(distances_array)),
            'max': float(np.max(distances_array)),
            'median': float(np.median(distances_array)),
            'percentile_25': float(np.percentile(distances_array, 25)),
            'percentile_75': float(np.percentile(distances_array, 75))
        }
        
        results['angle_stats'] = {
            'count': len(angles_array),
            'mean': float(np.mean(angles_array)),
            'std': float(np.std(angles_array)),
            'min': float(np.min(angles_array)),
            'max': float(np.max(angles_array)),
            'median': float(np.median(angles_array)),
            'percentile_25': float(np.percentile(angles_array, 25)),
            'percentile_75': float(np.percentile(angles_array, 75))
        }
        
        print(f"Distance statistics: {results['distance_stats']}")
        print(f"Angle statistics: {results['angle_stats']}")
        
        # Use the converted arrays for further processing
        distances = distances_array.tolist()
        angles = angles_array.tolist()
    else:
        print("Warning: No distances or angles calculated")
    
    # 6. Save intermediate results
    if config.get('save_intermediate', False):
        # Save distance-angle pair data
        if distances and angles:
            dist_angle_file = os.path.join(data_id_output_dir, f"{data_id}_actin_to_vesicle_dist_angle_pairs.csv")
            df = pd.DataFrame({'Distance': distances, 'Angle': angles})
            df.to_csv(dist_angle_file, index=False)
            results['dist_angle_file'] = dist_angle_file
            print(f"Saved distance-angle pairs to: {dist_angle_file}")
    
    # 7. Generate visualization
    if config.get('visualization', False) and distances and angles:
        print("Generating visualization...")
        
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
    print(f"Actin-to-vesicle distance-angle analysis completed successfully for {data_id}")
    
    return results


def actin_to_microtube_analysis(data_id: str, actin_file: str, microtube_file: str, config: Dict) -> Dict:
    """
    Calculate distance analysis from actin filaments to microtubules
    
    Analyzes spatial relationships between actin filaments and microtubules in 
    the cytoskeleton, calculating shortest distances between two types of 
    filamentous structures.
    
    Args:
        data_id (str): Dataset identifier
        actin_file (str): Actin point coordinates file path (.json format)
        microtube_file (str): Microtubule data file path (.json format)
        config (Dict): Configuration parameters dictionary
    
    Returns:
        Dict: Analysis result dictionary containing actin-to-microtubule distance analysis results
            - total_distance_vectors (int): Total number of distance vectors
            - total_shortest_distances (int): Total number of shortest distances
    
    Examples:
        >>> result = actin_to_microtube_analysis(
        ...     data_id="cytoskeleton_01",
        ...     actin_file="actin.json",
        ...     microtube_file="microtube.json",
        ...     config=config
        ... )
    """
    results = {}
    
    # 1. Load data
    print(f"Loading data for actin-to-microtube analysis: {data_id}...")
    
    # Load actin points
    if not os.path.exists(actin_file):
        raise FileNotFoundError(f"Actin file not found: {actin_file}")
    
    with open(actin_file, 'r') as f:
        actin_data = json.load(f)
    
    # Load microtube data
    if not os.path.exists(microtube_file):
        raise FileNotFoundError(f"Microtube file not found: {microtube_file}")
    
    with open(microtube_file, 'r') as f:
        MT_data = json.load(f)
    
    print(f"Actin filaments count: {len(actin_data)}")
    print(f"Microtube data loaded: {len(MT_data)} entries")
    
    # 2. Use original functions to calculate actin to microtube distance vector mapping
    print("Calculating actin-to-microtube distance vectors...")
    filament_vect_2_MT_map_lst = distances_generator.actin_to_microtube_distance(actin_data, MT_data)
    print(f"Calculated using original function: {len(filament_vect_2_MT_map_lst)} actin filaments processed")
    
    # 3. Use original function to calculate shortest distances
    print("Calculating shortest distances using original function...")
    distances = distances_generator.actin_to_microtube_distance2(actin_data, MT_data)
    print(f"Calculated using original function: {len(distances)} shortest distances")
    
    # 4. Create output directory with data_id
    base_output_dir = config.get('output_dir', os.path.dirname(actin_file))
    data_id_output_dir = os.path.join(base_output_dir, data_id)
    os.makedirs(data_id_output_dir, exist_ok=True)
    
    # 5. Process results
    results['data_id'] = data_id
    results['total_actin_filaments'] = len(actin_data)
    results['total_microtube_entries'] = len(MT_data)
    results['total_distance_vectors'] = len(filament_vect_2_MT_map_lst)
    results['total_shortest_distances'] = len(distances)
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
        # Save distance vector mapping
        if filament_vect_2_MT_map_lst:
            filament_vect_2_MT_map_lst_save = []
            for filament in filament_vect_2_MT_map_lst:
                cur_vect = [[vect[0], vect[1], vect[2]] for vect in filament]
                filament_vect_2_MT_map_lst_save.append(cur_vect)
            
            vector_map_file = os.path.join(data_id_output_dir, f"{data_id}_microtube_surround_actin_vect_map.json")
            with open(vector_map_file, 'w', encoding='utf-8') as f:
                json.dump(filament_vect_2_MT_map_lst_save, f)
            results['vector_map_file'] = vector_map_file
            print(f"Saved vector map to: {vector_map_file}")
        
        # Save shortest distances
        if distances:
            dist_file = os.path.join(data_id_output_dir, f"{data_id}_microtube_surround_actin_dist.csv")
            dist_savepd = pd.DataFrame({'shortest distance': distances})
            dist_savepd.to_csv(dist_file, index=False)
            results['distance_file'] = dist_file
            print(f"Saved distances to: {dist_file}")
    
    # 7. Generate visualization
    if config.get('visualization', False) and filament_vect_2_MT_map_lst and distances:
        print("Generating visualization...")
        
        # Prepare vector data for visualization
        filament_vect_2_MT_map_lst_new = []
        vect_count_n = 0
        for filament in filament_vect_2_MT_map_lst:
            curcoord = [np.array(vect) for vect in filament]
            filament_vect_2_MT_map_lst_new.append(curcoord)
            vect_count_n += len(curcoord)
        
        # Use original plotting functions
        if outplot is not None:
            if hasattr(outplot, 'vis_microtube_surround_actin'):
                outplot.vis_microtube_surround_actin(
                    filament_vect_2_MT_map_lst_new, 
                    data_id, 
                    len(filament_vect_2_MT_map_lst_new), 
                    vect_count_n
                )
                plot_file = os.path.join(data_id_output_dir, f"{data_id}_microtube_surround_actin_vectors.png")
                plt.show()
                plt.savefig(plot_file)
                plt.close()
                results['vector_plot_file'] = plot_file
                print(f"Saved vector plot to: {plot_file}")
            
            if hasattr(outplot, 'vis_microtube_surround_actin_dist'):
                # Prepare distance data for visualization
                dist_df = pd.DataFrame({'shortest distance': distances})
                outplot.vis_microtube_surround_actin_dist(
                    dist_df['shortest distance'].to_numpy(), 
                    data_id, 
                    len(distances), 
                    vect_count_n
                )
                dist_plot_file = os.path.join(data_id_output_dir, f"{data_id}_microtube_surround_actin_dist.png")
                plt.show()
                plt.savefig(dist_plot_file)
                plt.close()
                results['distance_plot_file'] = dist_plot_file
                print(f"Saved distance plot to: {dist_plot_file}")
    
    elif config.get('visualization', False):
        print("Warning: No data to visualize")
    
    results['status'] = 'success'
    print(f"Actin-to-microtube analysis completed successfully for {data_id}")
    
    return results


def actin_to_mito_analysis(data_id: str, actin_file: str, mito_file: str, config: Dict) -> Dict:
    """
    Calculate distance analysis from actin filaments to mitochondria
    
    Analyzes spatial relationships between actin filaments and mitochondria, 
    studying interaction patterns between cytoskeleton and organelles.
    
    Args:
        data_id (str): Dataset identifier
        actin_file (str): Actin point coordinates file path (.json format)
        mito_file (str): Mitochondria mask file path (.mrc format)
        config (Dict): Configuration parameters dictionary
    
    Returns:
        Dict: Analysis result dictionary containing actin-to-mitochondria distance statistics
    

    """
    results = {}
    
    # 1. Load data
    print(f"Loading data for actin-to-mito analysis: {data_id}...")
    
    # Load actin points
    if not os.path.exists(actin_file):
        raise FileNotFoundError(f"Actin file not found: {actin_file}")
    
    with open(actin_file, 'r') as f:
        actin_data = json.load(f)
    
    print(f"Actin filaments count: {len(actin_data)}")
    
    # 2. Try to read existing distance data
    distances_csv_file = os.path.join(os.path.dirname(actin_file), f"{data_id}_actin_to_mito_distances.csv")
    
    if os.path.exists(distances_csv_file):
        print(f"Found existing distance file: {distances_csv_file}")
        dist_df = pd.read_csv(distances_csv_file)
        distances = dist_df['Distance'].tolist()
        print(f"Loaded {len(distances)} pre-calculated distances")
    else:
        print("No existing distance file found, would need to calculate...")
        print("Skipping calculation due to coordinate boundary issues")
        results['status'] = 'error'
        results['error'] = 'No pre-calculated distances available and coordinate issues prevent calculation'
        return results
    
    # 3. Create output directory with data_id
    base_output_dir = config.get('output_dir', os.path.dirname(actin_file))
    data_id_output_dir = os.path.join(base_output_dir, data_id)
    os.makedirs(data_id_output_dir, exist_ok=True)
    
    # 4. Process results
    results['data_id'] = data_id
    results['total_actin_filaments'] = len(actin_data)
    results['total_distances'] = len(distances)
    results['distance_source'] = 'pre-calculated CSV file'
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
    
    # 5. Generate visualization
    if config.get('visualization', False) and distances:
        print("Generating visualization...")
        
        # Use original plotting function
        if outplot is not None and hasattr(outplot, 'actindist2mito_hist'):
            # Prepare distance data for visualization
            dist_df = pd.DataFrame({'Distance': distances})
            outplot.actindist2mito_hist(
                dist_df['Distance'].to_numpy(), 
                data_id, 
                len(actin_data), 
                len(distances)
            )
            plot_file = os.path.join(data_id_output_dir, f"{data_id}_actin_to_mito_analysis.png")
            plt.show()
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            results['plot_file'] = plot_file
            print(f"Visualization saved using original function: {plot_file}")
        
    elif config.get('visualization', False):
        print("Warning: No distances to visualize")
    
    results['status'] = 'success'
    print(f"Actin-to-mito analysis completed successfully for {data_id}")
    
    return results


def microtube_to_actin_analysis(data_id: str, actin_file: str, microtube_file: str, config: Dict) -> Dict:
    """
    Calculate distance analysis from microtubules to actin filaments (reverse analysis)
    
    Analyzes distance relationships from microtubules to actin filaments, 
    providing another perspective on cytoskeletal network interactions.
    
    Args:
        data_id (str): Dataset identifier
        actin_file (str): Actin point coordinates file path (.json format)
        microtube_file (str): Microtubule data file path (.json format)
        config (Dict): Configuration parameters dictionary
    
    Returns:
        Dict: Analysis result dictionary containing microtubule-to-actin distance and vector analysis results
    

    """
    results = {}
    
    # 1. Load data
    print(f"Loading data for microtube-to-actin analysis: {data_id}...")
    
    # Load actin points
    if not os.path.exists(actin_file):
        raise FileNotFoundError(f"Actin file not found: {actin_file}")
    
    with open(actin_file, 'r') as f:
        actin_data = json.load(f)
    
    # Load microtube data
    if not os.path.exists(microtube_file):
        raise FileNotFoundError(f"Microtube file not found: {microtube_file}")
    
    with open(microtube_file, 'r') as f:
        MT_data = json.load(f)
    
    print(f"Actin filaments count: {len(actin_data)}")
    print(f"Microtube data loaded: {len(MT_data)} entries")
    
    # 2. Use original functions to calculate microtube to actin distance vector mapping
    print("Calculating microtube-to-actin distance vectors...")
    filament_vect_2_MT_map_lst = distances_generator.actin_to_microtube_distance(actin_data, MT_data)
    print(f"Calculated using original function: {len(filament_vect_2_MT_map_lst)} actin filaments processed")
    
    # 3. Use original function to calculate shortest distances (following original class method)
    print("Calculating shortest distances using original function...")
    distances = distances_generator.actin_to_microtube_distance2(actin_data, MT_data)
    print(f"Calculated using original function: {len(distances)} shortest distances")
    
    # 4. Create output directory with data_id
    base_output_dir = config.get('output_dir', os.path.dirname(actin_file))
    data_id_output_dir = os.path.join(base_output_dir, data_id)
    os.makedirs(data_id_output_dir, exist_ok=True)
    
    # 5. Process results
    results['data_id'] = data_id
    results['total_actin_filaments'] = len(actin_data)
    results['total_microtube_entries'] = len(MT_data)
    results['total_distance_vectors'] = len(filament_vect_2_MT_map_lst)
    results['total_shortest_distances'] = len(distances)
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
        # Save distance vector mapping
        if filament_vect_2_MT_map_lst:
            filament_vect_2_MT_map_lst_save = []
            for filament in filament_vect_2_MT_map_lst:
                cur_vect = [[vect[0], vect[1], vect[2]] for vect in filament]
                filament_vect_2_MT_map_lst_save.append(cur_vect)
            
            vector_map_file = os.path.join(data_id_output_dir, f"{data_id}_microtube_to_actin_vect_map.json")
            with open(vector_map_file, 'w', encoding='utf-8') as f:
                json.dump(filament_vect_2_MT_map_lst_save, f)
            results['vector_map_file'] = vector_map_file
            print(f"Saved vector map to: {vector_map_file}")
        
        # Save shortest distances
        if distances:
            dist_file = os.path.join(data_id_output_dir, f"{data_id}_microtube_to_actin_dist.csv")
            dist_savepd = pd.DataFrame({'shortest distance': distances})
            dist_savepd.to_csv(dist_file, index=False)
            results['distance_file'] = dist_file
            print(f"Saved distances to: {dist_file}")
    
    # 7. Generate visualization (using original plotting methods from the class)
    if config.get('visualization', False) and filament_vect_2_MT_map_lst and distances:
        print("Generating visualization using original plotting methods...")
        
        # Load the saved vector map for plotting (following original class pattern)
        with open(os.path.join(data_id_output_dir, f"{data_id}_microtube_to_actin_vect_map.json"), 'r', encoding='utf-8') as f:
            filament_vect_2_MT_map_lst_plot = json.load(f)
        
        print(f"{data_id} MT num: {len(filament_vect_2_MT_map_lst_plot)}")
        
        # Prepare vector data for visualization (convert to numpy arrays as in original)
        filament_vect_2_MT_map_lst_new = []
        vect_count_n = 0
        for filament in filament_vect_2_MT_map_lst_plot:
            curcoord = [np.array(vect) for vect in filament]
            filament_vect_2_MT_map_lst_new.append(curcoord)
            vect_count_n += len(curcoord)
        
        # Use original plotting functions
        if outplot is not None:
            # Plot 1: Vector visualization
            if hasattr(outplot, 'vis_microtube_surround_actin'):
                outplot.vis_microtube_surround_actin(
                    filament_vect_2_MT_map_lst_new, 
                    data_id, 
                    len(filament_vect_2_MT_map_lst_new), 
                    vect_count_n
                )
                # Save the vector plot to our output directory as PDF
                vector_plot_file = os.path.join(data_id_output_dir, f"{data_id}_microtube_surrounding_probability_density_distribution.pdf")
                plt.show()
                plt.savefig(vector_plot_file, format='pdf', dpi=1200, bbox_inches='tight')
                plt.close()
                results['vector_plot_file'] = vector_plot_file
                print(f"Saved vector plot to: {vector_plot_file}")
            
            # Plot 2: Distance distribution visualization
            if hasattr(outplot, 'vis_microtube_surround_actin_dist'):
                # Load distance data for visualization
                dist_df = pd.DataFrame({'shortest distance': distances})
                outplot.vis_microtube_surround_actin_dist(
                    dist_df['shortest distance'].to_numpy(), 
                    data_id, 
                    len(distances), 
                    vect_count_n
                )
                # Save the distance distribution plot to our output directory
                dist_plot_file = os.path.join(data_id_output_dir, f"{data_id}_microtube_surround_actin_dist.png")
                plt.show()
                plt.savefig(dist_plot_file, dpi=300, bbox_inches='tight')
                plt.close()
                results['distance_plot_file'] = dist_plot_file
                print(f"Saved distance plot to: {dist_plot_file}")
    
    elif config.get('visualization', False):
        print("Warning: No data to visualize")
    
    results['status'] = 'success'
    print(f"Actin-to-microtube analysis completed successfully for {data_id}")
    
    return results


def mito_to_endoreticulum_analysis(data_id: str, mito_file: str, er_file: str, config: Dict,) -> Dict:
    """
    Calculate distance analysis from mitochondria to endoplasmic reticulum (using image processing)
    
    Uses image processing methods to analyze spatial relationships between mitochondria 
    and endoplasmic reticulum, studying interaction patterns between two important organelles.
    
    Args:
        data_id (str): Dataset identifier
        mito_file (str): Mitochondria mask file path (.mrc format)
        er_file (str): Endoplasmic reticulum mask file path (.mrc format)
        config (Dict): Configuration parameters dictionary
    
    Returns:
        Dict: Analysis result dictionary containing mitochondria-to-ER distance statistics
    

    """
    results = {}
    
    try:
        # 1. Load data
        print(f"Loading data for mito-to-ER analysis: {data_id}...")
        
        # Load mito mask
        if not os.path.exists(mito_file):
            raise FileNotFoundError(f"Mito file not found: {mito_file}")
        
        with mrcfile.open(mito_file, permissive=True) as mrc:
            mito_mask = mrc.data
        
        # Load ER mask
        if not os.path.exists(er_file):
            raise FileNotFoundError(f"ER file not found: {er_file}")
        
        with mrcfile.open(er_file, permissive=True) as mrc:
            er_mask = mrc.data
        
        print(f"Mito mask shape: {mito_mask.shape}")
        print(f"ER mask shape: {er_mask.shape}")
        
        # 2. Get configuration parameters
        voxel_size_xyz = config.get('voxel_size_xyz', [1.0, 1.0, 1.0])
        zoom_degree = config.get('zoom_degree', 0.10)  # Default zoom degree
        
        # 3. Calculate distances using original algorithm
        print("Calculating mito-to-ER distances using image processing...")
        
        # Resize masks
        print(f"Resizing masks with zoom degree: {zoom_degree}")
        mito_mask_resized = scipy.ndimage.zoom(mito_mask, zoom=zoom_degree, order=0)
        er_mask_resized = scipy.ndimage.zoom(er_mask, zoom=zoom_degree, order=0)
        
        # Get edge coordinates
        mito_edt_mask = scipy.ndimage.distance_transform_edt(mito_mask_resized)
        mito_edge_coords = np.array(np.where(mito_edt_mask == 1)).T
        
        # For ER, use all non-zero points (following original code)
        er_coords = np.array(np.where(er_mask_resized != 0)).T
        er_edge_coords = er_coords
        
        print(f"Mito edge points: {len(mito_edge_coords)}")
        print(f"ER points: {len(er_edge_coords)}")
        
        if len(mito_edge_coords) == 0 or len(er_edge_coords) == 0:
            print("Warning: No points found for one or both masks")
            results['status'] = 'error'
            results['error'] = 'No points found'
            return results
        
        # Calculate distances in batches using scipy.spatial.distance.cdist
        batchsize = config.get('batch_size', 1024)
        nn = len(mito_edge_coords) // batchsize + 1
        distances = []
        
        print(f"Processing {nn} batches with batch size {batchsize}")
        
        # Convert to float for distance calculation
        mito_edge_coords = mito_edge_coords.astype(np.float64)
        er_edge_coords = er_edge_coords.astype(np.float64)
        
        for ii in range(nn):
            start_idx = batchsize * ii
            end_idx = min(batchsize * (ii + 1), len(mito_edge_coords))
            cur_mito_coords = mito_edge_coords[start_idx:end_idx]
            
            if len(cur_mito_coords) == 0:
                break
                
            # Calculate distances between current batch and all ER points
            batch_distances = cdist(cur_mito_coords, er_edge_coords, metric='euclidean')
            
            # Add all distances to the list
            for jj in range(len(batch_distances)):
                distances.extend(batch_distances[jj])
            
            print(f"Processed batch {ii + 1}/{nn}")
        
        print(f"Total distances calculated: {len(distances)}")
        
        # Convert distances to nm
        avesize = float(np.mean(voxel_size_xyz))
        distances_nm = [float(dist) * (avesize/10.0) * (1.0/zoom_degree) for dist in distances]
        
        # 4. Create output directory with data_id
        base_output_dir = config.get('output_dir', os.path.dirname(mito_file))
        data_id_output_dir = os.path.join(base_output_dir, data_id)
        os.makedirs(data_id_output_dir, exist_ok=True)
        
        # 5. Process results
        results['data_id'] = data_id
        results['total_distances'] = len(distances_nm)
        results['zoom_degree'] = zoom_degree
        results['output_directory'] = data_id_output_dir
        
        if distances_nm:
            results['distance_stats'] = {
                'count': len(distances_nm),
                'mean': float(np.mean(distances_nm)),
                'std': float(np.std(distances_nm)),
                'min': float(np.min(distances_nm)),
                'max': float(np.max(distances_nm)),
                'median': float(np.median(distances_nm)),
                'percentile_25': float(np.percentile(distances_nm, 25)),
                'percentile_75': float(np.percentile(distances_nm, 75))
            }
            print(f"Distance statistics (nm): {results['distance_stats']}")
        
        # 6. Save intermediate results
        if config.get('save_intermediate', False):
            if distances_nm:
                dist_file = os.path.join(data_id_output_dir, f"{data_id}_mito_to_endoreticulum_distances.csv")
                dist_savepd = pd.DataFrame(columns=['Distance'], data=distances_nm)
                dist_savepd.to_csv(dist_file, index=False)
                results['distance_file'] = dist_file
                print(f"Saved distances to: {dist_file}")
        
        # 7. Generate visualization
        if config.get('visualization', False) and distances_nm:
            print("Generating visualization...")
            
            # Use original plotting function
            if outplot is not None and hasattr(outplot, 'dist_mito2er_hist'):
                dist_df = pd.DataFrame(columns=['Distance'], data=distances_nm)
                outplot.dist_mito2er_hist(dist_df['Distance'].to_numpy(), data_id, len(distances_nm))
                plot_file = os.path.join(data_id_output_dir, f"{data_id}_mito_to_endoreticulum_analysis.png")
                plt.show()
                plt.savefig(plot_file, dpi=300, bbox_inches='tight')
                plt.close()
                results['plot_file'] = plot_file
                print(f"Visualization saved using original function: {plot_file}")
        
        elif config.get('visualization', False):
            print("Warning: No distances to visualize")
        
        results['status'] = 'success'
        print(f"Mito-to-ER analysis completed successfully for {data_id}")
    
    except Exception as e:
        print(f"Error in mito-to-ER analysis: {str(e)}")
        results['status'] = 'error'
        results['error'] = str(e)
    
    return results






def actin_to_endoreticulum_analysis(data_id: str, actin_file: str, er_file: str, config: Dict) -> Dict:
    """
    Calculate distance analysis from actin filaments to endoplasmic reticulum membranes,
    using original distances_generator and outplot functions.
    
    Args:
        data_id (str): Dataset identifier
        actin_file (str): Actin filament coordinates file (.json)
        er_file (str): ER mask file (.mrc)
        config (Dict): Configuration dictionary
            - voxel_size_xyz (List[float])
            - shift_bias (List[float])
            - visualization (bool)
            - save_intermediate (bool)
            - output_dir (str, optional)
    
    Returns:
        Dict: Analysis result dictionary with distance stats, saved files, and visualizations
    """
    results = {}

    # 1. Load data
    print(f"Loading data for actin-to-endoreticulum analysis: {data_id}...")

    # Load actin data
    if not os.path.exists(actin_file):
        raise FileNotFoundError(f"Actin file not found: {actin_file}")
    with open(actin_file, 'r') as f:
        actin_data = json.load(f)
    
    # Load ER mask
    if not os.path.exists(er_file):
        raise FileNotFoundError(f"ER file not found: {er_file}")
    with mrcfile.open(er_file, permissive=True) as mrc:
        er_mask = mrc.data

    print(f"Actin filaments count: {len(actin_data)}")
    print(f"ER mask shape: {er_mask.shape}")

    # 2. Get configuration parameters
    shift_bias = config.get('shift_bias', [0.0, 0.0, 0.0])
    voxel_size_xyz = config.get('voxel_size_xyz', [1.0, 1.0, 1.0])

    # 3. Apply shift and calculate distances
    print("Calculating actin-to-endoreticulum distances using original function...")
    actin_data_shifted = distances_generator.shift_bias(actin_data, shift_bias)
    distances = distances_generator.coords_to_mem_distance_generator(actin_data_shifted, er_mask, voxel_size_xyz)
    print(f"Calculated: {len(distances)} actin-to-endoreticulum distances")

    # 4. Create output directory with data_id
    base_output_dir = config.get('output_dir', os.path.dirname(actin_file))
    data_id_output_dir = os.path.join(base_output_dir, data_id)
    os.makedirs(data_id_output_dir, exist_ok=True)

    # 5. Process results
    results['data_id'] = data_id
    results['total_actin_filaments'] = len(actin_data)
    results['total_distances'] = len(distances)
    results['image_shape'] = list(er_mask.shape)
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
        if distances:
            distance_file = os.path.join(data_id_output_dir, f"{data_id}_actin_to_endoreticulum_distances.csv")
            df = pd.DataFrame({'Distance': distances})
            df.to_csv(distance_file, index=False)
            results['distance_file'] = distance_file
            print(f"Saved distances to: {distance_file}")

    # 7. Generate visualization
    if config.get('visualization', False) and distances:
        print("Generating visualization using original function...")
        dist_df = pd.DataFrame({'Distance': distances})
        outplot.actindist2er_hist(dist_df['Distance'].to_numpy(), data_id, len(actin_data), len(distances))
        plot_file = os.path.join(data_id_output_dir, f"{data_id}_actin_to_endoreticulum_analysis.png")
        plt.show()
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        results['plot_file'] = plot_file
        print(f"Visualization saved: {plot_file}")

    elif config.get('visualization', False):
        print("Warning: No distances to visualize")

    results['status'] = 'success'
    print(f"Actin-to-endoreticulum analysis completed successfully for {data_id}")

    return results

def microtube_to_vesicle_analysis(data_id: str, microtube_file: str, vesicle_file: str, config: Dict) -> Dict:
    """
    Calculate distance analysis from microtubules to vesicle membranes,
    using original distances_generator and outplot functions.
    
    Args:
        data_id (str): Dataset identifier
        microtube_file (str): Microtubule point coordinates file path (.json format)
        vesicle_file (str): Vesicle mask file path (.mrc format)
        config (Dict): Configuration dictionary
            - voxel_size_xyz (List[float])
            - shift_bias (List[float])
            - visualization (bool)
            - save_intermediate (bool)
            - output_dir (str, optional)
    
    Returns:
        Dict: Analysis result with distance statistics, saved files, and plot results
    """
    results = {}
    
    # 1. Load data
    print(f"Loading data for microtube-to-vesicle analysis: {data_id}...")
    
    # Load microtube points
    if not os.path.exists(microtube_file):
        raise FileNotFoundError(f"Microtube file not found: {microtube_file}")
    with open(microtube_file, 'r') as f:
        MT_data = json.load(f)
    
    # Load vesicle mask
    if not os.path.exists(vesicle_file):
        raise FileNotFoundError(f"Vesicle file not found: {vesicle_file}")
    with mrcfile.open(vesicle_file, permissive=True) as mrc:
        vesicle_mask = mrc.data
    
    print(f"Microtube filaments count: {len(MT_data)}")
    print(f"Vesicle mask shape: {vesicle_mask.shape}")
    
    # 2. Get configuration parameters
    shift_bias = config.get('shift_bias', [0.0, 0.0, 0.0])
    voxel_size_xyz = config.get('voxel_size_xyz', [1.0, 1.0, 1.0])
    
    # 3. Apply shift and calculate distances
    # print("Calculating microtube-to-vesicle distances using original function...")
    MT_data_shifted = distances_generator.shift_bias(MT_data, shift_bias)
    distances = distances_generator.coords_to_mem_distance_generator(MT_data_shifted, vesicle_mask, voxel_size_xyz)
    print(f"Calculated: {len(distances)} microtube-to-vesicle distances")
    
    # 4. Create output directory with data_id
    base_output_dir = config.get('output_dir', os.path.dirname(microtube_file))
    data_id_output_dir = os.path.join(base_output_dir, data_id)
    os.makedirs(data_id_output_dir, exist_ok=True)
    
    # 5. Process Results
    results['data_id'] = data_id
    results['total_microtubes'] = len(MT_data)
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
        if distances:
            distance_file = os.path.join(data_id_output_dir, f"{data_id}_microtube_to_vesicle_distances.csv")
            df = pd.DataFrame({'Distance': distances})
            df.to_csv(distance_file, index=False)
            results['distance_file'] = distance_file
            print(f"Saved distances to: {distance_file}")
    
    # 7. Generate visualization
    if config.get('visualization', False) and distances:
        print("Generating visualization using original function...")
        # 用原始 outplot 绘图方法，可以替换成你 class 的 dist_mt2isghist
        dist_df = pd.DataFrame({'Distance': distances})
        outplot.dist_mt2isghist(dist_df['Distance'].to_numpy(), data_id, len(MT_data), len(distances))
        plot_file = os.path.join(data_id_output_dir, f"{data_id}_microtube_to_vesicle_analysis.png")
        plt.show()
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        results['plot_file'] = plot_file
        print(f"Visualization saved: {plot_file}")
    
    elif config.get('visualization', False):
        print("Warning: No distances to visualize")
    
    results['status'] = 'success'
    print(f"Microtube-to-vesicle analysis completed successfully for {data_id}")
    
    return results


def microtube_to_mito_analysis(data_id: str, microtube_file: str, mito_file: str, config: Dict) -> Dict:
    """
    Calculate distance analysis from microtubules to mitochondria membranes,
    using original distances_generator and outplot functions.
    
    Args:
        data_id (str): Dataset identifier
        microtube_file (str): Microtubule coordinates file (.json)
        mito_file (str): Mitochondria mask file (.mrc)
        config (Dict): Configuration dictionary
            - voxel_size_xyz (List[float])
            - shift_bias (List[float])
            - visualization (bool)
            - save_intermediate (bool)
            - output_dir (str, optional)
    
    Returns:
        Dict: Analysis result with distance stats, saved result files, and plots
    """
    results = {}

    # 1. Load data
    print(f"Loading data for microtube-to-mitochondria analysis: {data_id}...")

    # Load microtubule data
    if not os.path.exists(microtube_file):
        raise FileNotFoundError(f"Microtube file not found: {microtube_file}")
    with open(microtube_file, 'r') as f:
        MT_data = json.load(f)
    
    # Load mitochondria mask
    if not os.path.exists(mito_file):
        raise FileNotFoundError(f"Mito file not found: {mito_file}")
    with mrcfile.open(mito_file, permissive=True) as mrc:
        mito_mask = mrc.data

    print(f"Microtube filaments count: {len(MT_data)}")
    print(f"Mitochondria mask shape: {mito_mask.shape}")

    # 2. Get configuration parameters
    shift_bias = config.get('shift_bias', [0.0, 0.0, 0.0])
    voxel_size_xyz = config.get('voxel_size_xyz', [1.0, 1.0, 1.0])

    # 3. Apply shift and calculate distances
    print("Calculating microtube-to-mitochondria distances using original function...")
    MT_data_shifted = distances_generator.shift_bias(MT_data, shift_bias)
    distances = distances_generator.coords_to_mem_distance_generator(MT_data_shifted, mito_mask, voxel_size_xyz)
    print(f"Calculated: {len(distances)} microtube-to-mito distances")

    # 4. Create output directory with data_id
    base_output_dir = config.get('output_dir', os.path.dirname(microtube_file))
    data_id_output_dir = os.path.join(base_output_dir, data_id)
    os.makedirs(data_id_output_dir, exist_ok=True)

    # 5. Process results
    results['data_id'] = data_id
    results['total_microtubes'] = len(MT_data)
    results['total_distances'] = len(distances)
    results['image_shape'] = list(mito_mask.shape)
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
        if distances:
            distance_file = os.path.join(data_id_output_dir, f"{data_id}_microtube_to_mito_distances.csv")
            df = pd.DataFrame({'Distance': distances})
            df.to_csv(distance_file, index=False)
            results['distance_file'] = distance_file
            print(f"Saved distances to: {distance_file}")

    # 7. Generate visualization
    if config.get('visualization', False) and distances:
        print("Generating visualization using original function...")
        dist_df = pd.DataFrame({'Distance': distances})
        outplot.dist_mt2mitohist(dist_df['Distance'].to_numpy(), data_id, len(MT_data), len(distances))
        plot_file = os.path.join(data_id_output_dir, f"{data_id}_microtube_to_mito_analysis.png")
        plt.show()
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        results['plot_file'] = plot_file
        print(f"Visualization saved: {plot_file}")

    elif config.get('visualization', False):
        print("Warning: No distances to visualize")

    results['status'] = 'success'
    print(f"Microtube-to-mitochondria analysis completed successfully for {data_id}")

    return results


def microtube_to_endoreticulum_analysis(data_id: str, microtube_file: str, er_file: str, config: Dict) -> Dict:
    """
    Calculate distance analysis from microtubules to endoplasmic reticulum membranes,
    using original distances_generator and outplot functions.
    
    Args:
        data_id (str): Dataset identifier
        microtube_file (str): Microtubule coordinates file (.json)
        er_file (str): ER mask file (.mrc)
        config (Dict): Configuration dictionary
            - voxel_size_xyz (List[float])
            - shift_bias (List[float])
            - visualization (bool)
            - save_intermediate (bool)
            - output_dir (str, optional)
    
    Returns:
        Dict: Analysis result with distance stats, saved result files, and plots
    """
    results = {}

    # 1. Load data
    print(f"Loading data for microtube-to-endoreticulum analysis: {data_id}...")

    # Load microtubule data
    if not os.path.exists(microtube_file):
        raise FileNotFoundError(f"Microtube file not found: {microtube_file}")
    with open(microtube_file, 'r') as f:
        MT_data = json.load(f)
    
    # Load endoreticulum mask
    if not os.path.exists(er_file):
        raise FileNotFoundError(f"Endoreticulum file not found: {er_file}")
    with mrcfile.open(er_file, permissive=True) as mrc:
        er_mask = mrc.data

    print(f"Microtube filaments count: {len(MT_data)}")
    print(f"Endoreticulum mask shape: {er_mask.shape}")

    # 2. Get configuration parameters
    shift_bias = config.get('shift_bias', [0.0, 0.0, 0.0])
    voxel_size_xyz = config.get('voxel_size_xyz', [1.0, 1.0, 1.0])

    # 3. Apply shift and calculate distances
    print("Calculating microtube-to-endoreticulum distances using original function...")
    MT_data_shifted = distances_generator.shift_bias(MT_data, shift_bias)
    distances = distances_generator.coords_to_mem_distance_generator(MT_data_shifted, er_mask, voxel_size_xyz)
    print(f"Calculated: {len(distances)} microtube-to-endoreticulum distances")

    # 4. Create output directory with data_id
    base_output_dir = config.get('output_dir', os.path.dirname(microtube_file))
    data_id_output_dir = os.path.join(base_output_dir, data_id)
    os.makedirs(data_id_output_dir, exist_ok=True)

    # 5. Process results
    results['data_id'] = data_id
    results['total_microtubes'] = len(MT_data)
    results['total_distances'] = len(distances)
    results['image_shape'] = list(er_mask.shape)
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
        if distances:
            distance_file = os.path.join(data_id_output_dir, f"{data_id}_microtube_to_endoreticulum_distances.csv")
            df = pd.DataFrame({'Distance': distances})
            df.to_csv(distance_file, index=False)
            results['distance_file'] = distance_file
            print(f"Saved distances to: {distance_file}")

    # 7. Generate visualization
    if config.get('visualization', False) and distances:
        print("Generating visualization using original function...")
        dist_df = pd.DataFrame({'Distance': distances})
        outplot.dist_mt2erhist(dist_df['Distance'].to_numpy(), data_id, len(MT_data), len(distances))
        plot_file = os.path.join(data_id_output_dir, f"{data_id}_microtube_to_endoreticulum_analysis.png")
        plt.show()
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        results['plot_file'] = plot_file
        print(f"Visualization saved: {plot_file}")

    elif config.get('visualization', False):
        print("Warning: No distances to visualize")

    results['status'] = 'success'
    print(f"Microtube-to-endoreticulum analysis completed successfully for {data_id}")

    return results


def vesicle_to_mito_analysis(data_id: str, vesicle_file: str, mito_file: str, config: Dict) -> Dict:
    """
    Calculate distance analysis from vesicle membranes to mitochondria membranes using image processing,
    compatible with original batch-wise calculation and outplot functions.
    
    Args:
        data_id (str): Dataset identifier
        vesicle_file (str): Vesicle mask file (.mrc)
        mito_file (str): Mitochondria mask file (.mrc)
        config (Dict): Configuration dictionary
            - voxel_size_xyz (List[float])
            - visualization (bool)
            - save_intermediate (bool)
            - output_dir (str, optional)
            - zoom_degree (float, optional, default 0.10)
            - batch_size (int, optional, default 1024)
    
    Returns:
        Dict: Analysis result dictionary with distance stats, saved files, and plots
    """
    results = {}

    # 1. Load data
    print(f"Loading data for vesicle-to-mito analysis: {data_id}...")

    if not os.path.exists(vesicle_file):
        raise FileNotFoundError(f"Vesicle file not found: {vesicle_file}")
    if not os.path.exists(mito_file):
        raise FileNotFoundError(f"Mito file not found: {mito_file}")

    with mrcfile.open(vesicle_file, permissive=True) as mrc:
        vesicle_mask = mrc.data
    with mrcfile.open(mito_file, permissive=True) as mrc:
        mito_mask = mrc.data

    print(f"Vesicle mask shape: {vesicle_mask.shape}")
    print(f"Mito mask shape: {mito_mask.shape}")

    # 2. Get configuration parameters
    voxel_size_xyz = config.get('voxel_size_xyz', [1.0, 1.0, 1.0])
    zoom_degree = config.get('zoom_degree', 0.10)
    batch_size = config.get('batch_size', 1024)

    # 3. Downsample image for batch processing
    print(f"Resizing masks with zoom degree: {zoom_degree}")
    vesicle_mask_r = scipy.ndimage.zoom(vesicle_mask, zoom=zoom_degree, order=0)
    mito_mask_r = scipy.ndimage.zoom(mito_mask, zoom=zoom_degree, order=0)

    # 4. Extract edge coordinates
    vesicle_edt_mask = scipy.ndimage.distance_transform_edt(vesicle_mask_r)
    vesicle_edge_coords = np.array(np.where(vesicle_edt_mask == 1)).T
    mito_edt_mask = scipy.ndimage.distance_transform_edt(mito_mask_r)
    mito_edge_coords = np.array(np.where(mito_edt_mask == 1)).T

    print(f"Vesicle edge points: {len(vesicle_edge_coords)}")
    print(f"Mito edge points: {len(mito_edge_coords)}")

    if len(vesicle_edge_coords) == 0 or len(mito_edge_coords) == 0:
        results['status'] = 'error'
        results['error'] = 'No edge points found'
        print("Warning: No edge points found for one or both masks.")
        return results

    # 5. Compute distances in batches using torch for memory efficiency
    print("Calculating distances in batches using torch...")
    try:
        import torch

        vesicle_edge_coords = vesicle_edge_coords.astype(np.float64)
        mito_edge_coords = mito_edge_coords.astype(np.float64)
        vesicle_edge_tensor = torch.from_numpy(vesicle_edge_coords)
        mito_edge_tensor = torch.from_numpy(mito_edge_coords)
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        vesicle_edge_tensor = vesicle_edge_tensor.to(device)
        mito_edge_tensor = mito_edge_tensor.to(device)

        n_batches = len(vesicle_edge_tensor) // batch_size + 1
        all_distances = []

        for ii in range(n_batches):
            start_idx = batch_size * ii
            end_idx = min(batch_size * (ii+1), len(vesicle_edge_tensor))
            cur_vesicle_tensor = vesicle_edge_tensor[start_idx:end_idx]
            if cur_vesicle_tensor.size(0) == 0:
                break
            temp_dist = torch.cdist(cur_vesicle_tensor, mito_edge_tensor)
            temp_dist_np = temp_dist.cpu().numpy()
            for row in temp_dist_np:
                all_distances.extend(row)
            del temp_dist, temp_dist_np

        print(f"Total distances calculated: {len(all_distances)}")

    except Exception as e:
        results['status'] = 'error'
        results['error'] = f"Failed to batch-calculate distances: {str(e)}"
        print(f"Error during torch distance calculation: {str(e)}")
        return results

    # 6. Convert to nanometers
    avesize = np.mean(voxel_size_xyz)
    all_distances_nm = [dist * (avesize / 10) * (1 / zoom_degree) for dist in all_distances]

    # 7. Create output directory with data_id
    base_output_dir = config.get('output_dir', os.path.dirname(vesicle_file))
    data_id_output_dir = os.path.join(base_output_dir, data_id)
    os.makedirs(data_id_output_dir, exist_ok=True)

    # 8. Process results
    results['data_id'] = data_id
    results['total_distances'] = len(all_distances_nm)
    results['zoom_degree'] = zoom_degree
    results['output_directory'] = data_id_output_dir

    if all_distances_nm:
        stats = {
            'count': len(all_distances_nm),
            'mean': float(np.mean(all_distances_nm)),
            'std': float(np.std(all_distances_nm)),
            'min': float(np.min(all_distances_nm)),
            'max': float(np.max(all_distances_nm)),
            'median': float(np.median(all_distances_nm)),
            'percentile_25': float(np.percentile(all_distances_nm, 25)),
            'percentile_75': float(np.percentile(all_distances_nm, 75))
        }
        results['distance_stats'] = stats
        print(f"Distance statistics (nm): {stats}")

    # 9. Save intermediate results
    if config.get('save_intermediate', False) and all_distances_nm:
        distance_file = os.path.join(data_id_output_dir, f"{data_id}_vesicle_to_mito_distances.csv")
        df = pd.DataFrame({'Distance': all_distances_nm})
        df.to_csv(distance_file, index=False)
        results['distance_file'] = distance_file
        print(f"Distances saved to: {distance_file}")

    # 10. Generate visualization
    if config.get('visualization', False) and all_distances_nm:
        print("Generating visualization using matplotlib...")
        

        plt.figure(figsize=(10, 6))
        

        plt.hist(all_distances_nm, bins=50, alpha=0.7, density=True, color='skyblue', edgecolor='black')
        

        mean_dist = np.mean(all_distances_nm)
        median_dist = np.median(all_distances_nm)
        plt.axvline(mean_dist, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_dist:.2f} nm')
        plt.axvline(median_dist, color='green', linestyle='--', linewidth=2, label=f'Median: {median_dist:.2f} nm')
        
        plt.xlabel('Distance (nm)')
        plt.ylabel('Density')
        plt.title(f'{data_id} Vesicle to Mitochondria Distance Distribution\nTotal points: {len(all_distances_nm)}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plot_file = os.path.join(data_id_output_dir, f"{data_id}_vesicle_to_mito_analysis.png")
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.show()
        plt.close()
        results['plot_file'] = plot_file
        print(f"Visualization saved: {plot_file}")

    elif config.get('visualization', False):
        print("Warning: No distances to visualize")

    results['status'] = 'success'
    print(f"Vesicle-to-mitochondria analysis completed successfully for {data_id}")

    return results


def vesicle_to_endoreticulum_analysis(data_id: str, vesicle_file: str, er_file: str, config: Dict) -> Dict:
    """
    Calculate distance analysis from vesicle membranes to endoplasmic reticulum membranes using image processing,
    compatible with original batch-wise calculation and outplot functions.
    
    Args:
        data_id (str): Dataset identifier
        vesicle_file (str): Vesicle mask file (.mrc)
        er_file (str): ER mask file (.mrc)
        config (Dict): Configuration dictionary
            - voxel_size_xyz (List[float])
            - visualization (bool)
            - save_intermediate (bool)
            - output_dir (str, optional)
            - zoom_degree (float, optional, default 0.10)
            - batch_size (int, optional, default 1024)
    
    Returns:
        Dict: Analysis result dictionary with distance stats, saved files, and plots
    """
    results = {}

    # 1. Load data
    print(f"Loading data for vesicle-to-endoreticulum analysis: {data_id}...")

    if not os.path.exists(vesicle_file):
        raise FileNotFoundError(f"Vesicle file not found: {vesicle_file}")
    if not os.path.exists(er_file):
        raise FileNotFoundError(f"Endoreticulum file not found: {er_file}")

    with mrcfile.open(vesicle_file, permissive=True) as mrc:
        vesicle_mask = mrc.data
    with mrcfile.open(er_file, permissive=True) as mrc:
        er_mask = mrc.data

    print(f"Vesicle mask shape: {vesicle_mask.shape}")
    print(f"Endoreticulum mask shape: {er_mask.shape}")

    # 2. Get configuration parameters
    voxel_size_xyz = config.get('voxel_size_xyz', [1.0, 1.0, 1.0])
    zoom_degree = config.get('zoom_degree', 0.10)
    batch_size = config.get('batch_size', 1024)

    # 3. Downsample image for batch processing
    print(f"Resizing masks with zoom degree: {zoom_degree}")
    vesicle_mask_r = scipy.ndimage.zoom(vesicle_mask, zoom=zoom_degree, order=0)
    er_mask_r = scipy.ndimage.zoom(er_mask, zoom=zoom_degree, order=0)

    # 4. Extract edge coordinates
    vesicle_edt_mask = scipy.ndimage.distance_transform_edt(vesicle_mask_r)
    vesicle_edge_coords = np.array(np.where(vesicle_edt_mask == 1)).T
    er_coords = np.array(np.where(er_mask_r != 0)).T
    er_edge_coords = er_coords

    print(f"Vesicle edge points: {len(vesicle_edge_coords)}")
    print(f"ER edge points: {len(er_edge_coords)}")

    if len(vesicle_edge_coords) == 0 or len(er_edge_coords) == 0:
        results['status'] = 'error'
        results['error'] = 'No edge points found'
        print("Warning: No edge points found for one or both masks.")
        return results

    # 5. Compute distances in batches using torch for memory efficiency
    print("Calculating distances in batches using torch...")
    try:
        import torch

        vesicle_edge_coords = vesicle_edge_coords.astype(np.float64)
        er_edge_coords = er_edge_coords.astype(np.float64)
        vesicle_edge_tensor = torch.from_numpy(vesicle_edge_coords)
        er_edge_tensor = torch.from_numpy(er_edge_coords)
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        vesicle_edge_tensor = vesicle_edge_tensor.to(device)
        er_edge_tensor = er_edge_tensor.to(device)

        n_batches = len(vesicle_edge_tensor) // batch_size + 1
        all_distances = []

        for ii in range(n_batches):
            start_idx = batch_size * ii
            end_idx = min(batch_size * (ii+1), len(vesicle_edge_tensor))
            cur_vesicle_tensor = vesicle_edge_tensor[start_idx:end_idx]
            if cur_vesicle_tensor.size(0) == 0:
                break
            temp_dist = torch.cdist(cur_vesicle_tensor, er_edge_tensor)
            temp_dist_np = temp_dist.cpu().numpy()
            for row in temp_dist_np:
                all_distances.extend(row)
            del temp_dist, temp_dist_np

        print(f"Total distances calculated: {len(all_distances)}")

    except Exception as e:
        results['status'] = 'error'
        results['error'] = f"Failed to batch-calculate distances: {str(e)}"
        print(f"Error during torch distance calculation: {str(e)}")
        return results

    # 6. Convert to nanometers
    avesize = np.mean(voxel_size_xyz)
    all_distances_nm = [dist * (avesize / 10) * (1 / zoom_degree) for dist in all_distances]

    # 7. Create output directory with data_id
    base_output_dir = config.get('output_dir', os.path.dirname(vesicle_file))
    data_id_output_dir = os.path.join(base_output_dir, data_id)
    os.makedirs(data_id_output_dir, exist_ok=True)

    # 8. Process results
    results['data_id'] = data_id
    results['total_distances'] = len(all_distances_nm)
    results['zoom_degree'] = zoom_degree
    results['output_directory'] = data_id_output_dir

    if all_distances_nm:
        stats = {
            'count': len(all_distances_nm),
            'mean': float(np.mean(all_distances_nm)),
            'std': float(np.std(all_distances_nm)),
            'min': float(np.min(all_distances_nm)),
            'max': float(np.max(all_distances_nm)),
            'median': float(np.median(all_distances_nm)),
            'percentile_25': float(np.percentile(all_distances_nm, 25)),
            'percentile_75': float(np.percentile(all_distances_nm, 75))
        }
        results['distance_stats'] = stats
        print(f"Distance statistics (nm): {stats}")

    # 9. Save intermediate results
    if config.get('save_intermediate', False) and all_distances_nm:
        distance_file = os.path.join(data_id_output_dir, f"{data_id}_vesicle_to_endoreticulum_distances.csv")
        df = pd.DataFrame({'Distance': all_distances_nm})
        df.to_csv(distance_file, index=False)
        results['distance_file'] = distance_file
        print(f"Distances saved to: {distance_file}")

    # 10. Generate visualization
    if config.get('visualization', False) and all_distances_nm:
        print("Generating visualization using matplotlib...")
        
        plt.figure(figsize=(10, 6))
        
        plt.hist(all_distances_nm, bins=50, alpha=0.7, density=True, color='skyblue', edgecolor='black')
        
        mean_dist = np.mean(all_distances_nm)
        median_dist = np.median(all_distances_nm)
        plt.axvline(mean_dist, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_dist:.2f} nm')
        plt.axvline(median_dist, color='green', linestyle='--', linewidth=2, label=f'Median: {median_dist:.2f} nm')
        
        plt.xlabel('Distance (nm)')
        plt.ylabel('Density')
        plt.title(f'{data_id} Vesicle to Endoreticulum Distance Distribution\nTotal points: {len(all_distances_nm)}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plot_file = os.path.join(data_id_output_dir, f"{data_id}_vesicle_to_endoreticulum_analysis.png")
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.show()
        plt.close()
        results['plot_file'] = plot_file
        print(f"Visualization saved: {plot_file}")

    elif config.get('visualization', False):
        print("Warning: No distances to visualize")

    results['status'] = 'success'
    print(f"Vesicle-to-endoreticulum analysis completed successfully for {data_id}")

    return results


