from skimage.morphology import skeletonize_3d, skeletonize
import numpy as np

from skimage.morphology import skeletonize_3d
from skimage import filters, morphology
from scipy import ndimage
import json

def connect_skeleton_segments(skeleton, max_distance=3):
    """
    Connect nearby skeleton segments by finding close endpoints and bridging them
    
    Args:
        skeleton: 3D binary array of skeleton
        max_distance: maximum distance to connect segments
        
    Returns:
        connected_skeleton: skeleton with connected segments
    """
    from scipy.spatial import cKDTree
    
    # Find skeleton endpoints (pixels with only one neighbor)
    skeleton_coords = np.where(skeleton)
    skeleton_points = np.array(list(zip(*skeleton_coords)))
    
    if len(skeleton_points) == 0:
        return skeleton
    
    print(f"Processing all skeleton points: {len(skeleton_points)}")
    
    # Find endpoints more efficiently using vectorized operations
    endpoints = []
    
    # Pre-compute neighbor offsets for 26-connectivity
    neighbor_offsets = []
    for dz in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            for dx in [-1, 0, 1]:
                if dz == 0 and dy == 0 and dx == 0:
                    continue
                neighbor_offsets.append((dz, dy, dx))
    neighbor_offsets = np.array(neighbor_offsets)
    
    # Vectorized endpoint detection with progress reporting
    for i, point in enumerate(skeleton_points):
        if i % 10000 == 0:  # Progress indicator for large datasets
            print(f"Processed {i}/{len(skeleton_points)} points for endpoint detection")
            
        z, y, x = point
        
        # Check all 26 neighbors at once
        neighbor_coords = point + neighbor_offsets
        
        # Filter out coordinates that are outside the volume
        valid_mask = (
            (neighbor_coords[:, 0] >= 0) & (neighbor_coords[:, 0] < skeleton.shape[0]) &
            (neighbor_coords[:, 1] >= 0) & (neighbor_coords[:, 1] < skeleton.shape[1]) &
            (neighbor_coords[:, 2] >= 0) & (neighbor_coords[:, 2] < skeleton.shape[2])
        )
        
        valid_neighbors = neighbor_coords[valid_mask]
        
        # Count neighbors that are skeleton points
        neighbors = 0
        for nz, ny, nx in valid_neighbors:
            if skeleton[nz, ny, nx]:
                neighbors += 1
        
        if neighbors == 1:  # This is an endpoint
            endpoints.append(point)
    
    if len(endpoints) < 2:
        print("Not enough endpoints found for connection")
        return skeleton
    
    print(f"Found {len(endpoints)} endpoints from {len(skeleton_points)} skeleton points")
    
    # Use KDTree for efficient nearest neighbor search - built once for all endpoints
    endpoints = np.array(endpoints)
    tree = cKDTree(endpoints)
    
    # Find all pairs within distance threshold using KDTree
    print("Finding all endpoint pairs within connection distance...")
    pairs_to_connect = []
    
    # Use sparse_distance_matrix or query_pairs for efficient pair finding
    try:
        # Method 1: Use query_pairs (more memory efficient)
        from scipy.spatial.distance import pdist
        from scipy.sparse import coo_matrix
        
        # Find all pairs within max_distance
        pair_indices = tree.query_pairs(max_distance, output_type='ndarray')
        print(f"Found {len(pair_indices)} potential endpoint pairs to connect")
        
        pairs_to_connect = pair_indices
        
    except Exception as e:
        print(f"Using fallback method due to: {e}")
        # Fallback: process in chunks but ensure global connectivity
        chunk_size = 2000
        processed_pairs = set()
        
        for i in range(0, len(endpoints), chunk_size):
            chunk_end = min(i + chunk_size, len(endpoints))
            print(f"Processing chunk {i//chunk_size + 1}/{(len(endpoints) + chunk_size - 1)//chunk_size}")
            
            for j in range(i, chunk_end):
                # Query all endpoints within distance (not just current chunk)
                neighbor_indices = tree.query_ball_point(endpoints[j], max_distance)
                
                for k in neighbor_indices:
                    if j < k and (j, k) not in processed_pairs:  # Avoid duplicates and self-pairs
                        pairs_to_connect.append([j, k])
                        processed_pairs.add((j, k))
    
    # Connect the endpoint pairs
    print(f"Connecting {len(pairs_to_connect)} endpoint pairs...")
    connected_skeleton = skeleton.copy()
    connection_count = 0
    
    for pair_idx, (i, j) in enumerate(pairs_to_connect):
        if pair_idx % 1000 == 0:
            print(f"Connected {pair_idx}/{len(pairs_to_connect)} pairs")
        
        # Get the actual endpoint coordinates
        p1, p2 = endpoints[i], endpoints[j]
        distance = np.linalg.norm(p2 - p1)
        
        if distance <= max_distance and distance > 0:
            # Simple line drawing in 3D using Bresenham-like algorithm
            steps = max(1, int(np.ceil(distance * 2)))  # More steps for smoother lines
            for step in range(steps + 1):
                t = step / steps
                point = p1 + t * (p2 - p1)
                z, y, x = np.round(point).astype(int)
                if (0 <= z < skeleton.shape[0] and 
                    0 <= y < skeleton.shape[1] and 
                    0 <= x < skeleton.shape[2]):
                    connected_skeleton[z, y, x] = True
            connection_count += 1
    
    print(f"Successfully connected {connection_count} endpoint pairs")
    return connected_skeleton


def skeletonization_et_segmentation(image_matrix, 
                                gaussian_sigma=1.0,
                                threshold_multiplier=1.2,
                                erosion_radius=1,
                                min_object_size=200,
                                dilation_radius=1,
                                min_skeleton_size=30,
                                max_connect_distance=5,
                                final_min_size=20):
    """
    Extract 3D skeleton using global thresholding and 3D skeletonization
    
    Args:
        image_matrix: 3D numpy array, input SIM image data
        gaussian_sigma: float, sigma for Gaussian smoothing (default: 1.0)
        threshold_multiplier: float, multiplier for threshold calculation (default: 1.2)
        erosion_radius: int, radius for morphological erosion (default: 1)
        min_object_size: int, minimum size for object removal (default: 200)
        dilation_radius: int, radius for morphological dilation (default: 1)
        min_skeleton_size: int, minimum size for skeleton cleaning (default: 30)
        max_connect_distance: int, maximum distance for skeleton connection (default: 5)
        final_min_size: int, final minimum size for skeleton fragments (default: 20)
    
    Returns:
        skeleton_final: 3D numpy array, final skeletonized result
    """
    print(f"Input image stats: min={np.min(image_matrix)}, max={np.max(image_matrix)}, mean={np.mean(image_matrix):.2f}")

    # Gaussian smoothing for noise reduction
    smoothed = filters.gaussian(image_matrix, sigma=gaussian_sigma)
    
    # Use higher threshold to keep only bright internal regions, avoiding edges
    threshold = np.mean(smoothed) + threshold_multiplier * np.std(smoothed)
    binary_image = smoothed > threshold
    print(f"Global threshold: {threshold:.2f}")

    print(f"Foreground pixel ratio after thresholding: {np.sum(binary_image) / binary_image.size * 100:.2f}%")

    # Morphological erosion to further remove shape edges, keep only central regions
    binary_image = morphology.binary_erosion(binary_image, morphology.ball(erosion_radius))
    
    # Remove small objects
    binary_image = morphology.remove_small_objects(binary_image, min_size=min_object_size)
    
    # Slight dilation to recover some structures
    binary_image = morphology.binary_dilation(binary_image, morphology.ball(dilation_radius))

    # Extract 3D skeleton
    skeleton = skeletonize_3d(binary_image)

    # Minimal cleaning, only remove tiny skeleton fragments
    skeleton_cleaned = morphology.remove_small_objects(skeleton, min_size=min_skeleton_size)
    
    # Connect segmented skeleton parts with improved memory management
    print("Connecting skeleton segments...")
    skeleton_connected = connect_skeleton_segments(skeleton_cleaned, max_distance=max_connect_distance)
    
    # Clean again after connection to remove possible small fragments
    skeleton_final = morphology.remove_small_objects(skeleton_connected, min_size=final_min_size)

    print(f"Skeleton voxel count after cleaning and connecting: {np.sum(skeleton_final)}")

    return skeleton_final


def save_filament_branches_json(skeleton_result, data_id, output_dir, logger, interpolate=True, points_per_branch=100):
    """
    Group skeleton mask into individual filaments and save as JSON, 
    each branch as high-precision interpolated floating-point coordinate trajectory
    """
    structure = ndimage.generate_binary_structure(3, 2)  # 26-connectivity
    labeled, num_features = ndimage.label(skeleton_result, structure=structure)
    logger.step(f"Identified {num_features} actin branches in total")

    filaments = {}
    for label in range(1, num_features + 1):
        coords = np.argwhere(labeled == label)  # shape: [N,3]
        if len(coords) < 2:
            continue
        # ---------- Interpolation (fit trajectory, uniform sampling) ----------
        # Sort by spatial distance, can also use farthest point method
        # Here we use first point to last point for trajectory ordering
        coords = coords[np.argsort(coords[:,0])]  # Initial sort by z, adjust according to actual needs
        if interpolate and len(coords) > 2:
            # Calculate total distance length
            lengths = np.linalg.norm(coords[1:] - coords[:-1], axis=1)
            cumulative = np.concatenate([[0], np.cumsum(lengths)])
            total_length = cumulative[-1]
            interp_points = np.linspace(0, total_length, points_per_branch)
            interp_xyz = []
            for p in interp_points:
                idx = np.searchsorted(cumulative, p)
                if idx == 0:
                    interp_xyz.append(coords[0])
                elif idx >= len(coords):
                    interp_xyz.append(coords[-1])
                else:
                    prev, next = coords[idx-1], coords[idx]
                    ratio = (p - cumulative[idx-1]) / (cumulative[idx] - cumulative[idx-1] + 1e-6)
                    pt = prev + ratio * (next - prev)
                    interp_xyz.append(pt.tolist())
            filaments[f"filament_{label-1}"] = interp_xyz
        else:
            filaments[f"filament_{label-1}"] = coords.tolist()

    # Save as JSON
    json_file = output_dir / f"{data_id}_actin_filament_branches.json"
    with open(json_file, 'w') as f:
        json.dump(filaments, f, indent=2)
    logger.file_out(str(json_file))
    logger.step(f"Branch interpolated trajectory JSON saved to {json_file}")


# save_filament_branches_json(skeleton_result, data_id, output_dir, logger, interpolate=True, points_per_branch=100)