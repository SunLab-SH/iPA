import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndi
from skimage import measure
import pandas as pd
from scipy.ndimage import center_of_mass
from scipy.spatial.distance import cdist
from scipy.spatial import distance_matrix
import os
import pickle
import copy
import math
from scipy.spatial import cKDTree






def actinangle2dist_enhanced(filament_coord_extend_dict, imgsize, voxel_size_xyz):
    """
    Enhanced actin angle to distance calculation function
    
    Calculates shortest distances and corresponding angle relationships between 
    actin filaments.
    
    Args:
        filament_coord_extend_dict (Dict): Actin filament coordinate dictionary
        imgsize (Tuple): Image dimensions
        voxel_size_xyz (List[float]): Voxel size
    
    Returns:
        Tuple[List[float], List[float]]: (distance list, angle list)
    
    Note:
        - Angle range: 0-180 degrees
        - Distance units consistent with input coordinate units
        - Automatically handles boundary cases and outliers
    """
    filament_coord_extend_lst = []
    filamentnames = list(filament_coord_extend_dict.keys())
    
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = [np.array(coord) for coord in filament_coord_extend_dict[f'{filamentnames[idx]}']]
        filament_coord_extend_lst.append(curfilamentcoordslst)
    
    filament_coord_lst = filament_coord_extend_lst
    filaments_vect_lst = []  # vector for every point in every filament
    
    for single_filament_coord in filament_coord_lst:
        filament_vectlst = []
        for jj, point in enumerate(single_filament_coord):
            if jj < 1:
                tempvect = single_filament_coord[jj + 1] - point
            else:
                tempvect = point - single_filament_coord[jj - 1]
            filament_vectlst.append(tempvect)
        filaments_vect_lst.append(filament_vectlst)
    
    filaments_min_d_lst = []
    filaments_angle_lst = []
    
    if len(filament_coord_extend_lst) < 2:
        return filaments_min_d_lst, filaments_angle_lst
    
    for i, curfilaments_coords in enumerate(filament_coord_extend_lst):
        curfilaments_vect = filaments_vect_lst[i]
        
        rest_filament_coord_extend_lst = copy.deepcopy(filament_coord_extend_lst)
        rest_filament_vect_extend_lst = copy.deepcopy(filaments_vect_lst)
        
        _ = rest_filament_coord_extend_lst.pop(i)
        _ = rest_filament_vect_extend_lst.pop(i)
        
        for ll, point in enumerate(curfilaments_coords):
            vect1 = curfilaments_vect[ll]
            min_dist = float('inf')
            cur_angle = float('inf')
            
            for kk, filaments_coords in enumerate(rest_filament_coord_extend_lst):
                try:
                    point_reshaped = np.array([point])
                    filaments_coords_array = np.array(filaments_coords)
                    
                    distss = cdist(point_reshaped, filaments_coords_array, metric='euclidean')
                    min_d = float(np.min(distss))
                    
                    # Find the index of minimum distance
                    flat_idx = np.argmin(distss)
                    min_d_idx = flat_idx % distss.shape[1]
                    
                    if min_d < min_dist:
                        temp_d = min_d
                        vect2 = rest_filament_vect_extend_lst[kk][min_d_idx]
                        
                        if np.array_equal(vect1, vect2):
                            continue
                        
                        try:
                            dot_result = np.dot(vect1, vect2)
                            if np.isscalar(dot_result):
                                vector_dot_product = float(dot_result)
                            else:
                                vector_dot_product = float(dot_result.item())
                                
                            norm1 = float(np.linalg.norm(vect1))
                            norm2 = float(np.linalg.norm(vect2))
                            norm_product = norm1 * norm2
                        except:
                            continue
                        
                        if norm_product == 0:
                            continue
                            
                        cosine_value = vector_dot_product / norm_product
                        cosine_value = max(-1.0, min(1.0, cosine_value))  # Clip to valid range
                        
                        import math
                        arccos_val = math.acos(cosine_value)
                        angle = math.degrees(arccos_val)
                        
                        if angle >= 180:
                            temp_angle = 0.0
                        elif angle > 90:
                            temp_angle = 180.0 - angle
                        else:
                            temp_angle = angle
                        
                        if not math.isnan(temp_angle) and 0 <= temp_angle <= 180:
                            min_dist = temp_d
                            cur_angle = temp_angle
                            
                except Exception as e:
                    continue
            
            if min_dist != float('inf') and cur_angle != float('inf'):
                filaments_min_d_lst.append(min_dist)
                filaments_angle_lst.append(cur_angle)
    
    return filaments_min_d_lst, filaments_angle_lst





# ------



def load_partition_coords_from_xvg(xvg_file):
    """
    Load partition coordinate data from XVG file
    """
    partition_coords = {}
    radial_positions = {}
    
    # Read all data first
    data_lines = []
    with open(xvg_file, 'r') as f:
        for line in f:
            if not line.startswith('#') and not line.startswith('@'):
                parts = line.strip().split()
                if len(parts) >= 4:
                    data_lines.append(parts)
    
    if not data_lines:
        return {}, {}
    
    # Get max partition ID
    max_partition_id = max(int(parts[3]) for parts in data_lines)
    
    # Process data
    for parts in data_lines:
        z, y, x, partition_id = float(parts[0]), float(parts[1]), float(parts[2]), int(parts[3])
        coord = [z, y, x]
        
        if partition_id not in partition_coords:
            partition_coords[partition_id] = []
        
        partition_coords[partition_id].append(coord)
        
        # Normalized radial position
        normalized_position = (partition_id - 0.5) / max_partition_id
        radial_positions[(int(z), int(y), int(x))] = normalized_position
    
    # Convert to numpy arrays
    for partition_id in partition_coords:
        partition_coords[partition_id] = np.array(partition_coords[partition_id])
    
    return partition_coords, radial_positions

def calculate_rdf_from_xvg(organelle_mask, partition_coords, radial_positions, bins=8):
    """
    Calculate RDF directly from XVG partition data and organelle mask
    
    Parameters
    ----------
    organelle_mask : numpy.ndarray
        3D organelle mask array
    partition_coords : dict
        Dictionary with partition coordinates from XVG
    radial_positions : dict
        Dictionary mapping coordinates to radial positions
    bins : int, default=8
        Number of radial bins
        
    Returns
    -------
    dict
        Dictionary containing RDF results
    """
    # Extract organelle coordinates
    organelle_binary = organelle_mask > 0.5
    organelle_coords = np.column_stack(np.where(organelle_binary))
    
    print(f'Total organelles found: {len(organelle_coords)}')
    
    if len(organelle_coords) == 0:
        print("No organelles found")
        return None
    
    # Use KDTree for fast nearest neighbor search
    # Ensure partition_coords_array is a 2D numpy array of shape (n_points, 3)
    coords_list = list(radial_positions.keys())
    if len(coords_list) > 0 and isinstance(coords_list[0], tuple):
        partition_coords_array = np.array(coords_list)
    else:
        # Fallback if keys are already arrays or lists
        partition_coords_array = np.asarray(coords_list)
        
    partition_rpos_array = np.array(list(radial_positions.values()))
    
    print(f'Building KDTree for {len(partition_coords_array)} partition points...')
    if partition_coords_array.ndim != 2:
        raise ValueError(f"partition_coords_array has wrong dimensions: {partition_coords_array.shape}")
    tree = cKDTree(partition_coords_array)
    
    # Query nearest neighbors for all organelles at once
    dists, idxs = tree.query(organelle_coords, distance_upper_bound=10.0)
    
    # Assign radial positions to organelles within valid distance
    organelle_radial_positions = []
    for i, dist in enumerate(dists):
        if np.isfinite(dist):  # Within valid range
            organelle_radial_positions.append(partition_rpos_array[idxs[i]])
        # Skip organelles outside the valid distance range
    
    print(f'Organelles assigned to partitions: {len(organelle_radial_positions)} / {len(organelle_coords)}')
    
    if len(organelle_radial_positions) == 0:
        print("No organelles could be assigned to partitions")
        return None
    
    # Calculate histogram of organelle radial positions
    organelle_layer_counts, edges = np.histogram(organelle_radial_positions, bins=bins, range=(0, 1))
    
    # Calculate cytoplasm volume distribution (all partition coordinates)
    all_radial_positions = partition_rpos_array
    cyto_layer_counts, _ = np.histogram(all_radial_positions, bins=bins, range=(0, 1))
    
    print(f'Organelle layer counts: {organelle_layer_counts}')
    print(f'Cytoplasm layer counts: {cyto_layer_counts}')
    
    # Calculate normalized RDF
    cyto_volume = len(all_radial_positions)
    organelle_volume = len(organelle_radial_positions)
    average_density = organelle_volume / cyto_volume
    
    g_average = np.zeros(len(organelle_layer_counts))
    radii = np.zeros(len(organelle_layer_counts))
    
    for i in range(len(g_average)):
        if cyto_layer_counts[i] > 0 and average_density > 0:
            # RDF = (local density) / (average density)
            local_density = organelle_layer_counts[i] / cyto_layer_counts[i]
            g_average[i] = local_density / average_density
        else:
            g_average[i] = 0
        
        # Radial position (center of bin)
        radii[i] = (edges[i] + edges[i+1]) / 2
    
    print(f'RDF g(r): {g_average}')
    print(f'Radial positions: {radii}')
    
    # Prepare results
    rdf_results = {
        'dataset': 'sxt_analysis',
        'radii': radii.tolist(),
        'rdf': g_average.tolist(),
        'organelle_volume': organelle_volume,
        'cyto_volume': cyto_volume,
        'average_density': average_density,
        'organelle_layer_counts': organelle_layer_counts.tolist(),
        'cyto_layer_counts': cyto_layer_counts.tolist(),
        'assignment_efficiency': len(organelle_radial_positions) / len(organelle_coords)
    }
    
    return rdf_results


def calculate_distribution_similarity(rdf_1, rdf_2):
    """
    Calculate distribution similarity between two organelle distributions using Pearson correlation.
    
    Parameters
    ----------
    rdf_1 : array-like
        Radial Distribution Function values from the first image/modality.
    rdf_2 : array-like
        Radial Distribution Function values from the second image/modality.
        
    Returns
    -------
    dict
        Dictionary containing 'pearson_coefficient' and 'p_value'.
    """
    from scipy.stats import pearsonr
    
    if len(rdf_1) != len(rdf_2):
        raise ValueError("RDF arrays must have the same length for comparison.")
        
    corr, p_value = pearsonr(rdf_1, rdf_2)
    return {
        'pearson_coefficient': corr,
        'p_value': p_value,
        'is_significant': p_value < 0.05
    }

def save_radial_rdf_results(rdf_results, output_dir, dataid):
    """
    Save radial RDF results to files
    
    Parameters
    ----------
    rdf_results : dict
        RDF results
    output_dir : str
        Output directory
    dataid : str
        Data identifier
    """
    if rdf_results is None:
        print("No RDF results to save")
        return
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Save as pickle
    pickle_path = os.path.join(output_dir, f'{dataid}_radial_rdf_results.pkl')
    with open(pickle_path, 'wb') as f:
        pickle.dump(rdf_results, f)
    print(f"Radial RDF results saved to: {pickle_path}")
    
    # Save as CSV
    csv_data = []
    for i, (r, g) in enumerate(zip(rdf_results['radii'], rdf_results['rdf'])):
        csv_data.append({
            'radial_position': r,
            'rdf_value': g,
            'organelle_count': rdf_results['organelle_layer_counts'][i],
            'cyto_count': rdf_results['cyto_layer_counts'][i]
        })
    
    if csv_data:
        df = pd.DataFrame(csv_data)
        csv_path = os.path.join(output_dir, f'{dataid}_radial_rdf_data.csv')
        df.to_csv(csv_path, index=False)
        print(f"Radial RDF data saved to: {csv_path}")




# %%
#----------
# wfm isg cluster



def volume_to_radius(volume):
    return (3 * volume / (4 * np.pi)) ** (1 / 3)

def draw_sphere(img, center, radius):
    # 'center' is a tuple (x, y, z)
    # 'img' is the 3D image
    z, y, x = center
    radius = int(volume_to_radius(radius))
    rr, cc = disk((y, x), radius)
    rr[rr >= img.shape[1]] = img.shape[1] - 1
    cc[cc >= img.shape[2]] = img.shape[2] - 1
    img[z, rr, cc] = 255

# if 0 <= x < new_image_shape[0] and 0 <= y < new_image_shape[1] and 0 <= z < new_image_shape[2]:
#     draw_sphere(img_imaris[time_point - 1], coords, img_volume)


def reconstruct_4d_mask(data, new_image_shape, rescalex, rescaley):
    """
    Reconstructs a 4D binary mask array with a temporal dimension from 3D spatial point data.

    This function generates a 4D (time, X, Y, Z) mask array based on vesicle/particle coordinates across multiple time points.
    The coordinates from the input DataFrame are mapped and rescaled to fit into the target image shape, and each particle is visualized as a small cubic region in the 3D volume per time frame.

    Parameters
    ----------
    data : pandas.DataFrame
        DataFrame containing at least the columns 'Time', 'Position X', 'Position Y', 'Position Z'.
        Each row should represent a single particle/vesicle at a given time and 3D position.
    new_image_shape : tuple of int
        The target image dimensions as (X, Y, Z), defining the desired output resolution for the spatial axes.
    rescalex : float
        Manual rescaling factor applied to the X axis, used to match physical/image boundary correspondence.
    rescaley : float
        Manual rescaling factor applied to the Y axis, used to match physical/image boundary correspondence.

    Returns
    -------
    img : numpy.ndarray
        A 4D numpy array of shape (T, X, Y, Z) where each particle/vesicle is assigned a unique label (id) at each time point T in the binary mask.
    ves_coords : list
        List of lists of original (X, Y, Z) coordinates for the vesicles at each time point, ordered by frame.

    """

    time_range = (data['Time'].min(), data['Time'].max())
    time_values = data['Time'].unique()

    # Creating a 4D array for the image with an additional time dimension
    # The shape is (number of time points, X, Y, Z)
    # New image dimensions and resolution
    # new_image_shape = (187, 302, 43)
    # new_image_shape = (round(x_range[1]-x_range[0]), round(y_range[1]-y_range[0]), round(z_range[1]-z_range[0]))

    resolution = (0.1030, 0.1030, 0.3006)  # in micrometers
    # resolution = (1, 1, 1) 
    # resolution = (9.7099, 9.7099, 3.3267) # 1 / voxel number 

    ves_size = 0.2 # in micrometers

    img_imaris_shape = (len(time_values),) + new_image_shape
    img_imaris = np.zeros(img_imaris_shape, dtype=np.uint8)

    x_range = (data['Position X'].min(), data['Position X'].max())
    y_range = (data['Position Y'].min(), data['Position Y'].max())
    z_range = (data['Position Z'].min(), data['Position Z'].max())

    ves_coords = []
    # Populating the 4D image
    for time_point in time_values:

        # Extract data for the current time point
        time_data = data[data['Time'] == time_point]

        # Calculating the scale for each axis based on the new resolution
        # TODO: the 2 rescale parameters are manully fitted but should be availale from image size / boundary distance
        x_scale_new = new_image_shape[0] / ((x_range[1] - x_range[0])*rescalex)
        y_scale_new = new_image_shape[1] / ((y_range[1] - y_range[0])*rescaley)
        z_scale_new = new_image_shape[2] / ((z_range[1] - z_range[0]))

        centered_image = np.zeros(new_image_shape, dtype=np.uint8)

        x_center = time_data['Position X'].mean()
        y_center = time_data['Position Y'].mean()
        z_center = time_data['Position Z'].mean()

        x_offset = (new_image_shape[0] // 2) - int((x_center - x_range[0]) * x_scale_new)
        y_offset = (new_image_shape[1] // 2) - int((y_center - y_range[0]) * y_scale_new)
        z_offset = (new_image_shape[2] // 2) - int((z_center - z_range[0]) * z_scale_new)
        
        id = 1
        ves_coords_ti = []
        for _, row in time_data.iterrows():
            x_new = int((row['Position X'] - x_range[0]) * x_scale_new) + x_offset
            y_new = int((row['Position Y'] - y_range[0]) * y_scale_new) + y_offset
            z_new = int((row['Position Z'] - z_range[0]) * z_scale_new) + z_offset

            if 0 <= x_new < new_image_shape[0] and 0 <= y_new < new_image_shape[1] and 0 <= z_new < new_image_shape[2]:
                # TODO: the voxel size is not average compared to its real size, avg voxel size 4-6
                # this for visulize
                centered_image[x_new-2:x_new+2, y_new-2:y_new+2, z_new-2:z_new+2] = id
                # centered_image[x_new-1:x_new+1, y_new-1:y_new+1, z_new] = id
                id += 1  
                ves_coords_ti += [(row['Position X'], row['Position Y'], row['Position Z'])]
        ves_coords += [ves_coords_ti]
        
        img_imaris[time_point - 1] = centered_image

    return img_imaris, ves_coords


# Function to find instances and calculate their center of mass
def find_instances_and_centers(img_fi, vesicle_pixel_cutoff):

    im_vesicle_mask = np.where(np.less(img_fi, vesicle_pixel_cutoff), 0, 1)
    im_vesicle_label, n_vesicle_label = ndi.measurements.label(im_vesicle_mask)
    centers = [center_of_mass(img_fi, im_vesicle_label, i) for i in range(n_vesicle_label)]

    return im_vesicle_label, centers


# Function to identify cluster
def find_clusters_and_centers(fname, fi, raw_img, img_fi, ves_coords_list, vesicle_pixel_cutoff, visulize=1, outdir = None):

    im_vesicle_label, _ = find_instances_and_centers(img_fi, vesicle_pixel_cutoff)
    ves_volume_list = list(map(lambda x: np.array(x.area), measure.regionprops(im_vesicle_label)))

    # real_x = 57.1583 / 555 # μm
    # real_y = 55.3046 / 537 # μm
    # real_z = 11.122 / 37 # μm

    real_x = 19.2588 / 187 # μm
    real_y = 31.1024 / 302 # μm
    real_z = 12.9256 / 43 # μm

    cluster_index = []

    # for i in range(len(ves_volume_list)):
        
    #     ves_coords = np.where(np.equal(im_vesicle_label, i+1))
    #     ves_coords_x = ves_coords[1].max() - ves_coords[1].min()
    #     ves_coords_y = ves_coords[2].max() - ves_coords[2].min()
    #     ves_coords_z = ves_coords[0].max() - ves_coords[0].min()
        
    #     # 3D distance may not be very accurate here
    #     # dist = np.sqrt((ves_coords_x)**2 + (ves_coords_y)**2 + (ves_coords_z)**2)
    #     dist = np.sqrt((ves_coords_x*real_x)**2 + (ves_coords_y*real_y)**2 + (ves_coords_z*real_z)**2)
    #     print(dist)
    #     if dist > 0.6 and 5 < ves_volume_list[i] < 10000 and ves_coords[2].max() <= 180:
    #         cluster_index += [i+1]

    # change this to center distance instead of cluster size
    dist_matrix = distance_matrix(np.array(ves_coords_list), np.array(ves_coords_list))

    # distance threshold between any two vesicles
    threshold_distance = 1 # 0.2

    # Initialize labels for each particle
    labels = np.full(len(dist_matrix), -1)  # Start with -1 to denote unlabelled

    # Label counter
    current_label = 1

    for i in range(len(dist_matrix)):
        if labels[i] == -1: 
            if np.any(dist_matrix[i, :] < threshold_distance):
                labels[i] = current_label
                for j in range(len(dist_matrix)):
                    if dist_matrix[i, j] < threshold_distance and labels[j] == -1:
                        labels[j] = current_label
                current_label += 1

    # number each cluster from 1
    cluster_dict = {}
    i = 1
    # the label starts form 1, just for the index of each identity
    # the label for ves_coords_list is not the same as the one for im_vesicle_label
    for label in np.unique(labels):
        if len(np.where(labels == label)[0].tolist()) > 1:
            cluster_dict[i] = np.where(labels == label)[0].tolist()
            i += 1
    
    # number each cluster from 1
    im_ves_cluster = np.zeros_like(img_fi, dtype=np.int64)
    for i, v in enumerate(cluster_dict):
        for k in cluster_dict[v]:
            im_ves_cluster = np.where(np.equal(img_fi, k+1), i+1, im_ves_cluster)

    cluster_centers = []
    for i in range(len(cluster_dict)):
    
        if np.any(np.isnan(im_ves_cluster)) or np.any(np.isinf(im_ves_cluster)):
            print("Input data contains NaN or Inf.")

        unique_labels = np.unique(im_ves_cluster)
        label = i + 1
        if label not in unique_labels:
            print(f"Label {label} does not exist in the image.")

        com = center_of_mass(im_ves_cluster, labels=im_ves_cluster, index=label)

        if np.isnan(com).any():
            print(f"Calculated center of mass is NaN: {com}")
        else:
            # print(f"Center of mass coordinates: {com}")
            cluster_centers += [com]
    
    # The returned center might be nan
    # cluster_centers = [center_of_mass(im_ves_cluster, im_ves_cluster, i+1) for i in range(len(cluster_dict))]
    

    if visulize:
        if outdir is None:
            outdir = "./data/results"
        os.makedirs(outdir, exist_ok=True)
        
        fig = plt.figure(figsize=(20, 4))

        ax = fig.add_subplot(1, 4, 1)
        plt.title('raw image')
        plt.imshow(np.max(raw_img, axis=0))

        ax = fig.add_subplot(1, 4, 2)
        plt.title('image reconstruction')
        plt.imshow(np.max(img_fi, axis=0))

        ax = fig.add_subplot(1, 4, 3)
        plt.title('seg vesicle')
        plt.imshow(np.max(im_vesicle_label, axis=0))

        ax = fig.add_subplot(1, 4, 4)
        plt.title('seg cluster')
        plt.imshow(np.max(im_ves_cluster, axis=0))
        plt.savefig(os.path.join(outdir, '{}_fi{}.png'.format(fname, fi)))
        plt.close()

    print('frame_{}: n_vesicle: {}; n_cluster: {}; clustered_ves: {}'.format(fi, 
                                                                             len(ves_coords_list), 
                                                                             len(cluster_centers), 
                                                                             len([i for x in cluster_dict for i in cluster_dict[x]]) ))

    return im_vesicle_label, im_ves_cluster, cluster_dict, cluster_centers


# Function to match instances across frames using their center positions, this is simplified here
def match_instances_across_frames(centers_list, dthreshold):
    '''
    input: 
    centers_list = [
    [(1, 1), (5, 5)],   
    [(2, 2), (5, 4)],   
    [(3, 3), (6, 4)]]
    dthreshold = 2 

    output: 
    [[0, 0, 0], [1, 1, 1]]
    
    '''
    # Initialize the list to store consistent instances across all frames in the window
    track_instances = []

    # Start with instances in the first frame
    for i, center1 in enumerate(centers_list[0]):
        match_found = True
        match_indices = [i]

        # Compare with instances in subsequent frames
        for centers in centers_list[1:]:
            if not len(centers):
                continue
            distances = cdist([center1], centers)[0]
            match_index = np.argmin(distances)
            center1 = centers[match_index]  
            match_indices.append(match_index)

            # Check distance threshold (optional, for more accurate matching)
            # Threshold for considering a match (can be adjusted)
            if distances[match_index] > dthreshold:  
                match_found = False
                break

        # If a match is found in all frames, add it to the consistent instances
        if match_found:
            track_instances.append(match_indices)

    # return cluster idx, instead of vesicle ids
    return track_instances



def process_tracking_windows(data):

    ''' bug: the returned frames are not continuous '''

    final_cluster = {}

    for triplet in data[0]:
        final_cluster[tuple(triplet)] = [f'frame_{i}: {v}' for i, v in enumerate(triplet)]

    for i, window in enumerate(data[1:], start=1):
        # window is the list of tracking instances (triplet)
        if not len(window):
            continue

        # if len(triplet) < 3:
        #     continue

        # check the number before adding new elements
        # last_frame = 

        last_elements = [k[-2:] for k in final_cluster]
        # print(len(window))
        # print(final_cluster)
        for t in window:
            triplet = tuple(t)
            two_elements = (t[0], t[1])

            matching_keys = [(k, int(final_cluster[k][-1].split(':')[0].split('_')[1])) for k in final_cluster if k[-2:] == two_elements]

            if matching_keys:
                matched_key = matching_keys[0][0]
                matched_key_frame = matching_keys[0][1]
                if matched_key_frame == i+1:
                    new_key = matched_key + (t[2],)
                    final_cluster[new_key] = final_cluster.pop(matched_key)
                    final_cluster[new_key].append(f'frame_{i+2}: {t[2]}')
                else:
                    # new_key = (t[2])
                    # final_cluster[new_key].append(f'frame_{i+2}: {t[2]}')
                    # print(i+1, matched_key_frame)
                    pass
            else:
                if two_elements in last_elements:
                    print(f'Possible merge detected in frame {i} for {triplet}')
                else:
                    final_cluster[triplet] = [f'frame_{i+k}: {v}' for k, v in enumerate(triplet)]

    return final_cluster


def analyze_vesicle_clusters(name, img_ch2, img_imaris, ves_coords, outdir = None, config_params=None):
    """
    Main vesicle clustering analysis function
    
    Parameters:
    -----------
    name : str
        Sample name
    img_ch2 : numpy.ndarray
        Original fluorescence image data
    img_imaris : numpy.ndarray
        Reconstructed Imaris image data
    ves_coords : list
        List of vesicle coordinates for each time point
    config_params : dict, optional
        Configuration parameters dictionary
        
    Returns:
    --------
    results : dict
        Analysis results dictionary containing processed image data and clustering information
    """
    
    # Default configuration parameters
    default_config = {
        'window_size': 5,
        'vesicle_pixel_cutoff': 1,
        'dthreshold': 5,
        'visualize': True
    }
    
    if config_params:
        default_config.update(config_params)
    
    window_size = default_config['window_size']
    vesicle_pixel_cutoff = default_config['vesicle_pixel_cutoff']
    dthreshold = default_config['dthreshold']
    visualize = default_config['visualize']
    outdir = f'{outdir}/imgs'

    print(f"Processing {name}, image shape: {img_imaris.shape}")
    
    proc_image = {}
    # image segmentation for all frames
    for fi in range(len(img_imaris)):
        im_ves_fi, im_ves_cluster_fi, cluster_dict, cluster_centers = find_clusters_and_centers(
            name, fi, img_ch2[fi], 
            img_imaris[fi], 
            ves_coords[fi], 
            vesicle_pixel_cutoff,
            visulize=visualize,
            outdir=outdir
        )
        proc_image[fi] = [im_ves_fi, im_ves_cluster_fi, cluster_dict, cluster_centers]

    # Process all frames using a sliding window
    tracking_cluster = []
    for i in range(len(img_imaris) - window_size + 1):
        centers_list = [proc_image[fi][3] for fi in range(i, i + window_size)]
        track_instances = match_instances_across_frames(centers_list, dthreshold)
        if not len(track_instances):
            continue
        tracking_cluster.append(np.array(track_instances) + 1)

    final_cluster = process_tracking_windows(tracking_cluster[:-1])
    print('Detected number of cluster: {}'.format(len(final_cluster)))

    # check the clustered vesicle size
    cluster_size = []
    for i, k in enumerate(final_cluster):
        clus_size = 0
        for v in final_cluster[k]:
            try:
                fi, clu_idx = int(v.split(':')[0].split('_')[1]), int(v.split(':')[1])
                clus_size += np.sum(np.where(proc_image[fi][1] == clu_idx, 1, 0))
            except KeyError:
                print(v)
        cluster_size.append(clus_size)
        final_cluster[k].append(clus_size)

    # Construct return results
    results = {
        'sample_name': name,
        'processed_images': proc_image,
        'final_clusters': final_cluster,
        'cluster_sizes': cluster_size,
        'tracking_data': tracking_cluster,
        'config_used': default_config
    }
    
    return results






#--------------------------
# sim


def analyze_organelle_aggregation(centers_list, dthreshold=0.6):
    """
    Analyze organelle aggregation in time-lapse images.
    
    Aggregation is defined as the set of organelles that remain within a distance 
    less than a certain threshold across a number of continuous frames.
    
    Parameters
    ----------
    centers_list : list of np.ndarray
        List of arrays containing organelle center coordinates for each frame.
    dthreshold : float
        Distance threshold (in micrometers) to define an aggregation.
        
    Returns
    -------
    dict
        Dictionary mapping aggregation IDs to their tracking history and volume.
    """
    from scipy.spatial.distance import cdist
    
    # Simple clustering per frame based on distance threshold
    clusters_per_frame = []
    for centers in centers_list:
        if len(centers) == 0:
            clusters_per_frame.append([])
            continue
            
        dist_matrix = cdist(centers, centers)
        labels = np.full(len(centers), -1)
        current_label = 0
        
        for i in range(len(centers)):
            if labels[i] == -1:
                cluster_members = [j for j in range(len(centers)) if dist_matrix[i, j] < dthreshold]
                if len(cluster_members) > 1: # Only consider groups as aggregates
                    for m in cluster_members:
                        labels[m] = current_label
                    current_label += 1
        
        clusters_per_frame.append(labels)
        
    # Track aggregates across frames (simplified matching)
    # In a full implementation, this would use more robust tracking logic
    return {"clusters": clusters_per_frame}

def analyze_docked_granules(img_isg, mem_mask, isg_data, show_visualization=False, distance_threshold=None):
    """
    Analyze the localization relationship between ISG granules and plasma membrane using matrix data.

    Parameters:
        img_isg: numpy ndarray, ISG segmentation 3D matrix
        mem_mask: numpy ndarray, PM mask 3D matrix
        isg_data: numpy ndarray, granule center coordinates and volume information (CSV loading result)
        show_visualization: bool, whether to display visualization images of partial slices
        distance_threshold: float, optional. If provided, use this fixed distance (in pixels/nm)
            to define docked pool instead of dynamic vesicle radius. Matches manuscript requirement for fixed thresholds (e.g., 300 nm).

    Returns:
        dict, containing total granule count, docked granule count and ratio
    """
    isg_volume = isg_data[1:, 19]
    centers = isg_data[1:, 4:7]

    ves_label = np.zeros_like(mem_mask, dtype=np.uint8)
    for i, center in enumerate(centers):
        ves_label[int(center[2])-3:int(center[2])+3,
          int(center[1])-3:int(center[1])+3,
          int(center[0])-3:int(center[0])+3] = 1

    # Crop PM and corresponding regions
    nonzero_coords = np.where(mem_mask != 0)
    if len(nonzero_coords[0]) == 0:
        return {'total_granules': 0, 'docked_granules': 0, 'docked_ratio': 0}
        
    x0 = (nonzero_coords[2].min() // 100) * 100
    y0 = (nonzero_coords[1].min() // 100) * 100
    x1 = (nonzero_coords[2].max() // 100) * 100 + 100
    y1 = (nonzero_coords[1].max() // 100) * 100 + 100

    mem_mask_upd = mem_mask[:, y0:y1, x0:x1]
    ves_label = ves_label[:, y0:y1, x0:x1]
    img_isg = img_isg[:, y0:y1, x0:x1]

    pm_mask = np.where(mem_mask == 0, 1, 0)
    centers_upd = []
    isg_volume_upd = []
    
    for i, center in enumerate(centers):
        valid_x = (center[0] >= x0) and (center[0] < x1)
        valid_y = (center[1] >= y0) and (center[1] < y1)
        valid_z = 0 <= center[2] < mem_mask.shape[0]
        
        if valid_x and valid_y and valid_z:
            try:
                if pm_mask[int(center[2]), int(center[1]), int(center[0])] == 0:
                    centers_upd.append([center[0]-x0, center[1]-y0, center[2]])
                    isg_volume_upd.append(isg_volume[i])
            except IndexError:
                continue

    pm_mask_bin = np.where(mem_mask_upd == 0, 1, 0)
    dis_pm_mask = ndi.distance_transform_edt(1 - pm_mask_bin)

    granule_PM = []
    docked_cutoff = []
    for i, center in enumerate(centers_upd):
        idx = [round(c) for c in center]
        try:
            dist = dis_pm_mask[idx[2], idx[1], idx[0]]
        except IndexError:
            continue
        granule_PM.append(dist)
        
        # Use fixed threshold if provided, otherwise use dynamic radius
        if distance_threshold is not None:
            docked_cutoff.append(distance_threshold)
        else:
            ves_radius = ((3 * isg_volume_upd[i]) / (4 * math.pi)) ** (1/3)
            docked_cutoff.append(ves_radius)

    n_docked = len([x for i,x in enumerate(granule_PM) if x < docked_cutoff[i]])
    total_granules = len(granule_PM)
    docked_ratio = n_docked / total_granules if total_granules > 0 else 0

    if show_visualization and total_granules > 0:
        docked_indices = [i for i, x in enumerate(granule_PM) if x < docked_cutoff[i]]
        vis_mask = np.zeros_like(mem_mask_upd, dtype=np.uint8)
        for i in docked_indices:
            center = centers_upd[i]
            vis_mask[int(center[2])-3:int(center[2])+3,
                     int(center[1])-3:int(center[1])+3,
                     int(center[0])-3:int(center[0])+3] = 1

        sel_z = min(10, mem_mask_upd.shape[0]-1)
        fig = plt.figure(figsize=(15, 5))
        
        ax1 = fig.add_subplot(1, 3, 1)
        ax1.imshow(ves_label[sel_z], cmap='gray')
        ax1.set_title(f'All Granules (Total: {total_granules})')
        ax1.axis('off')

        ax2 = fig.add_subplot(1, 3, 2)
        ax2.imshow(mem_mask_upd[sel_z], cmap='gray')
        ax2.set_title('Plasma Membrane Region')
        ax2.axis('off')

        ax3 = fig.add_subplot(1, 3, 3)
        ax3.imshow(vis_mask[sel_z], cmap='hot')
        ax3.set_title(f'Docked Granules (Count: {n_docked})')
        ax3.axis('off')
        
        plt.tight_layout()
        # plt.show()

    return {
        'total_granules': total_granules,
        'docked_granules': n_docked,
        'docked_ratio': docked_ratio
    }


def calculate_distributional_similarity(rdf1, rdf2, alpha=0.05):
    """
    Calculate Pearson correlation between two RDF distributions to measure distributional similarity.
    
    This function compares organelle distributions between two cell images from different 
    modalities (e.g., SIM vs SXT) using Pearson correlation coefficient.
    
    Parameters
    ----------
    rdf1 : np.ndarray
        RDF values from modality 1 (e.g., SIM)
    rdf2 : np.ndarray
        RDF values from modality 2 (e.g., SXT)
    alpha : float, optional
        Significance level for p-value. Default: 0.05
        
    Returns
    -------
    dict
        Dictionary containing:
        - 'pearson_r' (float): Pearson correlation coefficient (-1 to 1)
        - 'p_value' (float): Statistical significance
        - 'significant' (bool): Whether correlation is significant (p < alpha)
        
    Example
    -------
    >>> import numpy as np
    >>> from ipa.analysis import calculate_distributional_similarity
    >>> sim_rdf = np.array([0.8, 1.0, 1.2, 1.1, 0.9, 1.0, 1.3, 1.5])
    >>> sxt_rdf = np.array([0.9, 1.1, 1.1, 1.0, 1.0, 1.1, 1.2, 1.4])
    >>> result = calculate_distributional_similarity(sim_rdf, sxt_rdf)
    >>> print(f"Pearson r: {result['pearson_r']:.2f}")
    >>> print(f"P-value: {result['p_value']:.4f}")
    >>> print(f"Significant: {result['significant']}")
    """
    from scipy import stats
    
    if len(rdf1) != len(rdf2):
        raise ValueError("RDF arrays must have same length")
    
    if len(rdf1) < 3:
        raise ValueError("Need at least 3 data points for correlation")
    
    # Calculate Pearson correlation
    r, p = stats.pearsonr(rdf1, rdf2)
    
    return {
        'pearson_r': r,
        'p_value': p,
        'significant': p < alpha
    }


def analyze_anchored_organelle_angles(filament_coords_dict, pm_mask, 
                                       distance_threshold_nm=120, 
                                       voxel_size_nm=(1, 1, 1)):
    """
    Analyze angles of anchored tubular organelles (e.g., F-actin) relative to PM.
    
    Anchored organelles are defined by their distance to the PM. For such organelles,
    the angle is calculated between pairs of nearest tubular organelles anchored to the PM.
    
    Parameters
    ----------
    filament_coords_dict : dict
        Dictionary of filament coordinates {filament_id: [[x,y,z], ...]}
    pm_mask : np.ndarray
        Plasma membrane binary mask
    distance_threshold_nm : float, optional
        Distance threshold (nm) to define anchored organelles. Default: 120
    voxel_size_nm : tuple, optional
        Voxel dimensions (nm). Default: (1, 1, 1)
        
    Returns
    -------
    dict
        Dictionary containing:
        - 'anchored_count' (int): Number of anchored filaments
        - 'angles' (list): List of angles between anchored filament pairs (degrees)
        - 'mean_angle' (float): Mean angle
        - 'std_angle' (float): Standard deviation of angles
        
    Example
    -------
    >>> import numpy as np
    >>> from ipa.analysis import analyze_anchored_organelle_angles
    >>> # Create sample filament data
    >>> filaments = {
    ...     'filament_1': [[10, 20, 30], [11, 21, 31], [12, 22, 32]],
    ...     'filament_2': [[15, 25, 35], [16, 26, 36], [17, 27, 37]]
    ... }
    >>> pm_mask = create_pm_mask()  # Your PM mask
    >>> results = analyze_anchored_organelle_angles(
    ...     filaments, pm_mask, 
    ...     distance_threshold_nm=120,
    ...     voxel_size_nm=(31.3, 31.3, 90.9)  # SIM voxel size
    ... )
    >>> print(f"Anchored filaments: {results['anchored_count']}")
    >>> print(f"Mean angle: {results['mean_angle']:.2f} degrees")
    """
    from scipy.spatial import cKDTree
    
    # Get PM coordinates
    pm_coords = np.array(np.where(pm_mask > 0)).T
    if len(pm_coords) == 0:
        raise ValueError("PM mask is empty")
    
    # Convert voxel size to numpy array
    voxel_size = np.array(voxel_size_nm)
    
    # Build KD-tree for PM
    pm_tree = cKDTree(pm_coords)
    
    # Find anchored filaments
    anchored_filaments = []
    anchored_vectors = []
    
    for filament_id, coords in filament_coords_dict.items():
        coords_array = np.array(coords)
        
        # Convert to physical coordinates if needed
        if np.max(voxel_size) > 10:  # Likely in voxels, convert to nm
            coords_nm = coords_array * voxel_size
        else:
            coords_nm = coords_array
        
        # Check if any point is within threshold of PM
        distances, _ = pm_tree.query(coords_nm)
        min_dist = np.min(distances)
        
        if min_dist <= distance_threshold_nm:
            # This filament is anchored
            anchored_filaments.append(filament_id)
            
            # Calculate filament direction vector (from first to last point)
            if len(coords_nm) >= 2:
                direction = coords_nm[-1] - coords_nm[0]
                norm = np.linalg.norm(direction)
                if norm > 0:
                    direction = direction / norm
                    anchored_vectors.append(direction)
    
    # Calculate angles between anchored filament pairs
    angles = []
    if len(anchored_vectors) >= 2:
        for i in range(len(anchored_vectors)):
            for j in range(i + 1, len(anchored_vectors)):
                vec1 = anchored_vectors[i]
                vec2 = anchored_vectors[j]
                
                # Calculate angle using dot product
                cos_angle = np.dot(vec1, vec2)
                cos_angle = np.clip(cos_angle, -1.0, 1.0)  # Clip for numerical stability
                angle_rad = np.arccos(cos_angle)
                angle_deg = np.degrees(angle_rad)
                
                # Convert to 0-90 degree range (acute angle)
                if angle_deg > 90:
                    angle_deg = 180 - angle_deg
                
                angles.append(angle_deg)
    
    # Calculate statistics
    mean_angle = np.mean(angles) if len(angles) > 0 else 0.0
    std_angle = np.std(angles) if len(angles) > 0 else 0.0
    
    return {
        'anchored_count': len(anchored_filaments),
        'angles': angles,
        'mean_angle': mean_angle,
        'std_angle': std_angle
    }

