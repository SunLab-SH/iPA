"""
Cellular Partitioning Module

This module implements radial cytoplasmic partitioning analysis, a novel approach for 
spatial zoning of cellular structures from nucleus to plasma membrane.
Through creating concentric radial partitions, it enables quantitative analysis of 
organelle distribution patterns across different cytoplasmic regions.
"""

import os
import numpy as np
import glob
from scipy import ndimage
import matplotlib.pyplot as plt
from tqdm import tqdm

# Add performance optimization dependencies
try:
    from numba import njit, prange
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    print("Warning: Numba not available. Performance optimizations will be limited.")

try:
    from joblib import Parallel, delayed
    JOBLIB_AVAILABLE = True
except ImportError:
    JOBLIB_AVAILABLE = False
    print("Warning: Joblib not available. Parallel processing will be disabled.")

# Numba-accelerated core computation functions
if NUMBA_AVAILABLE:
    @njit(parallel=True)
    def _vectorized_partition_assignment(coords, radial_dist, n_partitions, shape):
        """
        Numba-accelerated vectorized partition assignment with corrected logic
        """
        partition_ids = np.zeros(len(coords), dtype=np.int32)
        
        for i in prange(len(coords)):
            z, y, x = coords[i]
            if 0 <= z < shape[0] and 0 <= y < shape[1] and 0 <= x < shape[2]:
                base_radial_pos = radial_dist[z, y, x]
                
                # Corrected partition assignment logic
                # radial_dist ranges from 0 (at nucleus) to 1 (at PM)
                # Partition 1 = innermost (near nucleus), Partition n = outermost (near PM)
                if base_radial_pos <= 0.0:
                    partition_id = 1  # Very close to nucleus -> innermost partition
                elif base_radial_pos >= 1.0:
                    partition_id = n_partitions  # Very close to PM -> outermost partition
                else:
                    # Linear mapping from (0,1) to (1, n_partitions)
                    partition_id = int(base_radial_pos * n_partitions) + 1
                    partition_id = min(max(1, partition_id), n_partitions)
                
                partition_ids[i] = partition_id
            else:
                partition_ids[i] = 0
                
        return partition_ids

    @njit
    def _fast_distance_calculation(p1, p2):
        """
        Fast distance calculation
        """
        return np.sqrt(np.sum((p1 - p2) ** 2))

    @njit(parallel=True)
    def _vectorized_distance_filter(pairs, min_dist, max_dist):
        """
        Vectorized distance filtering
        """
        n_pairs = len(pairs)
        valid_mask = np.zeros(n_pairs, dtype=np.bool_)
        
        for i in prange(n_pairs):
            ne_pt = pairs[i, :3]
            pm_pt = pairs[i, 3:6]
            dist = _fast_distance_calculation(ne_pt, pm_pt)
            valid_mask[i] = min_dist <= dist <= max_dist
            
        return valid_mask
else:
    # Fallback functions for non-Numba version
    def _vectorized_partition_assignment(coords, radial_dist, n_partitions, shape):
        partition_ids = np.zeros(len(coords), dtype=np.int32)
        valid_mask = (
            (coords[:, 0] >= 0) & (coords[:, 0] < shape[0]) &
            (coords[:, 1] >= 0) & (coords[:, 1] < shape[1]) &  
            (coords[:, 2] >= 0) & (coords[:, 2] < shape[2])
        )
        
        valid_coords = coords[valid_mask]
        if len(valid_coords) > 0:
            radial_values = radial_dist[valid_coords[:, 0], valid_coords[:, 1], valid_coords[:, 2]]
            
            # Corrected partition assignment logic
            # radial_dist ranges from 0 (at nucleus) to 1 (at PM)
            # Partition 1 = innermost (near nucleus), Partition n = outermost (near PM)
            partition_values = np.zeros_like(radial_values, dtype=int)
            
            # Handle boundary cases
            inner_mask = radial_values <= 0.0
            outer_mask = radial_values >= 1.0
            middle_mask = (radial_values > 0.0) & (radial_values < 1.0)
            
            partition_values[inner_mask] = 1  # Innermost partition
            partition_values[outer_mask] = n_partitions  # Outermost partition
            
            # Linear mapping for middle values
            if np.any(middle_mask):
                middle_values = radial_values[middle_mask]
                partition_values[middle_mask] = np.clip(
                    (middle_values * n_partitions).astype(int) + 1,
                    1, n_partitions
                )
            
            partition_ids[valid_mask] = partition_values
            
        return partition_ids

    def _vectorized_distance_filter(pairs, min_dist, max_dist):
        distances = np.linalg.norm(pairs[:, :3] - pairs[:, 3:6], axis=1)
        return (distances >= min_dist) & (distances <= max_dist)

class Partitioning:
    """
    This class provides a comprehensive toolkit for cellular partitioning analysis,
    including boundary extraction, partition generation, feature calculation, and 
    visualization. Primarily used for analyzing organelle spatial distribution 
    between nuclear envelope and plasma membrane.
    
    Parameters
    
    root_dir : str
        Root directory path for result outputs
    n_slices : int, default=8
        Number of radial partitions, dividing nucleus-membrane distance into n_slices regions
    num_cores : int, optional
        Number of CPU cores for parallel processing, defaults to 1
        
    Attributes
    
    root_dir : str
        Output directory path
    n_slices : int
        Number of partitions
    num_cores : int
        Number of CPU cores
    """
    
    def __init__(self, root_dir, n_slices=8, num_cores=None):
        """
        Initialize Partitioning instance.
        :param root_dir: Directory to save outputs.
        :param n_slices: Number of shell partitions.
        :param num_cores: Number of CPU cores for parallel processing.
        """
        self.root_dir = root_dir
        self.n_slices = n_slices
        self.num_cores = num_cores
        os.makedirs(root_dir, exist_ok=True)
        os.makedirs(os.path.join(root_dir, 'results'), exist_ok=True)

    def load_mask_data(self, ne_path, pm_path):
        """
        Load nuclear envelope (NE) and plasma membrane (PM) mask data from files
        
        Supports multiple file formats including numpy arrays (.npy), TIFF images 
        (.tif/.tiff), MRC format (.mrc), and other common biological imaging formats.
        
        Parameters
        ----------
        ne_path : str
            Path to nuclear envelope mask file
        pm_path : str
            Path to plasma membrane mask file
            
        Returns
        -------
        tuple of numpy.ndarray
            (ne_mask, pm_mask) - Binary mask arrays for nuclear envelope and plasma membrane
            
        Raises
        ------
        FileNotFoundError
            When specified file paths do not exist
        ValueError
            When file format is unsupported or data format is incorrect
            
        Examples
        --------
        >>> ne_mask, pm_mask = partitioner.load_mask_data("nuclear.npy", "membrane.npy")
        >>> print(f"NE mask shape: {ne_mask.shape}, PM mask shape: {pm_mask.shape}")
        """
        # Implement loading logic here (e.g., using tifffile or mrcfile)
        # Placeholder:
        ne_mask = np.load(ne_path)
        pm_mask = np.load(pm_path)
        return ne_mask, pm_mask

    def extract_ne_pm_edges(self, pm_mask, ne_mask):
        """
        Extract nuclear envelope and plasma membrane boundary points and cellular center from 3D mask data
        
        Uses morphological operations to extract proper boundaries instead of EDT method
        which may produce incorrect boundary representations.
        
        Parameters
        ----------
        pm_mask : numpy.ndarray
            3D binary mask of plasma membrane, shape (Z, Y, X)
        ne_mask : numpy.ndarray
            3D binary mask of nuclear envelope, shape (Z, Y, X)
            
        Returns
        -------
        center : numpy.ndarray
            Cell center coordinates, shape (3,), order (Z, Y, X)
        ne_edge : numpy.ndarray
            Nuclear envelope boundary point coordinates array, shape (N, 3)
        pm_edge : numpy.ndarray
            Plasma membrane boundary point coordinates array, shape (M, 3)
        """
        from scipy.ndimage import binary_erosion, binary_fill_holes
        
        # Ensure binary masks
        pm_binary = pm_mask > 0
        ne_binary = ne_mask > 0
        
        # Fill holes to get solid regions
        pm_filled = binary_fill_holes(pm_binary)
        ne_filled = binary_fill_holes(ne_binary)
        
        # Extract boundaries using morphological operations
        # NE boundary: original mask minus eroded mask
        ne_eroded = binary_erosion(ne_filled, iterations=1)
        ne_boundary = ne_filled & ~ne_eroded
        
        # PM boundary: filled mask minus eroded mask (inner boundary)
        pm_eroded = binary_erosion(pm_filled, iterations=2)
        pm_boundary = pm_filled & ~pm_eroded
        
        # Get boundary coordinates
        ne_coords = np.column_stack(np.where(ne_boundary))
        pm_coords = np.column_stack(np.where(pm_boundary))
        
        # Calculate center as centroid of nuclear region
        ne_center_coords = np.column_stack(np.where(ne_filled))
        center = np.mean(ne_center_coords, axis=0)
        
        print(f"Extracted NE boundary points: {len(ne_coords)}")
        print(f"Extracted PM boundary points: {len(pm_coords)}")
        print(f"Cell center: {center}")
        
        return center, ne_coords, pm_coords

    @staticmethod
    def min_angle(ne_points, pm_point, center):
        """
        Find nuclear envelope point with minimal angle to plasma membrane point relative to center
        
        Calculates cosine angles between vectors to find the directionally closest 
        nuclear envelope point. Used for establishing optimal NE-PM point pair matching.
        
        Parameters
        ----------
        ne_points : numpy.ndarray
            Nuclear envelope point coordinates array, shape (N, 3)
        pm_point : numpy.ndarray
            Single plasma membrane point coordinates, shape (3,)
        center : numpy.ndarray
            Cell center coordinates, shape (3,)
            
        Returns
        -------
        numpy.ndarray
            Nuclear envelope point coordinates with minimal angle, shape (3,)
            
        Notes
        -----
        Uses cosine similarity to calculate vector angles:
        cos(θ) = (v1 · v2) / (|v1| × |v2|)
        Selects the point with maximum cos(θ), i.e., minimal angle.
        
        Examples
        --------
        >>> closest_ne = Partitioning.min_angle(ne_points, pm_point, center)
        >>> print(f"Closest NE point: {closest_ne}")
        """
        vectors_ne = ne_points - center
        vector_pm = pm_point - center
        
        # Vectorized calculation of all angles
        norms_ne = np.linalg.norm(vectors_ne, axis=1)
        norm_pm = np.linalg.norm(vector_pm)
        
        # Avoid division by zero
        valid_mask = norms_ne > 1e-8
        cos_angles = np.full(len(vectors_ne), -1.0)  # Initialize with minimum value
        
        if np.any(valid_mask):
            cos_angles[valid_mask] = np.dot(vectors_ne[valid_mask], vector_pm) / (norms_ne[valid_mask] * norm_pm + 1e-8)
        
        max_idx = np.argmax(cos_angles)
        return ne_points[max_idx]

    def find_ne_pm_pairs(self, center, ne_edge, pm_edge, 
                         step=8, save_txt=False, dataid=None):
        """
        Pair nuclear envelope and plasma membrane boundary points to establish radial connections
        
        Reduces computational complexity through controlled sampling strategies while maintaining 
        analysis accuracy. Algorithm prioritizes keeping NE point count within reasonable range, 
        then finds optimal matching PM points for each NE point.
        
        Parameters
        ----------
        center : numpy.ndarray
            Cell center coordinates
        ne_edge : numpy.ndarray
            Nuclear envelope boundary points array
        pm_edge : numpy.ndarray
            Plasma membrane boundary points array
        step : int, default=8
            Sampling step size for controlling point density (automatically adjusted based on target count)
        save_txt : bool, default=False
            Whether to save pairing results to text file
        dataid : str, optional
            Data identifier for file naming
            
        Returns
        -------
        numpy.ndarray
            Paired points array, shape (K, 6), each row contains [ne_z, ne_y, ne_x, pm_z, pm_y, pm_x]
            
        Notes
        -----
        Optimized version using KDTree for accelerated NE-PM pairing algorithm with distance filtering
        """
        from scipy.spatial import cKDTree
        
        print(f"Original NE edge points: {len(ne_edge)}")
        print(f"Original PM edge points: {len(pm_edge)}")

        # Intelligent sampling control
        target_ne_points = min(50000, len(ne_edge))
        if len(ne_edge) > target_ne_points:
            ne_step = len(ne_edge) // target_ne_points
            ne_sampled = ne_edge[::ne_step]
        else:
            ne_sampled = ne_edge
        
        target_pm_points = min(len(ne_sampled) * 8, len(pm_edge))
        if len(pm_edge) > target_pm_points:
            pm_step = len(pm_edge) // target_pm_points
            pm_sampled = pm_edge[::pm_step]
        else:
            pm_sampled = pm_edge
        
        print(f"Sampled NE edge points: {len(ne_sampled)}")
        print(f"Sampled PM edge points: {len(pm_sampled)}")

        # Optimized batch vector computation
        pm_vectors = pm_sampled - center
        pm_norms = np.linalg.norm(pm_vectors, axis=1)
        pm_unit_vectors = pm_vectors / (pm_norms[:, np.newaxis] + 1e-8)
        
        ne_vectors = ne_sampled - center
        ne_norms = np.linalg.norm(ne_vectors, axis=1)
        ne_unit_vectors = ne_vectors / (ne_norms[:, np.newaxis] + 1e-8)
        
        # Batch KDTree queries (key optimization)
        print("Building KDTree and performing batch queries...")
        pm_tree = cKDTree(pm_unit_vectors)
        
        # Batch query all NE points at once
        distances, indices = pm_tree.query(ne_unit_vectors, k=1)
        
        # Build initial pairing array
        pm_matched = pm_sampled[indices]
        pairs_array = np.column_stack([ne_sampled, pm_matched])
        
        # Vectorized distance calculation and filtering
        print("Applying vectorized distance filtering...")
        min_distance_threshold = 15.0
        max_distance_threshold = 200.0
        
        if NUMBA_AVAILABLE:
            valid_mask = _vectorized_distance_filter(
                pairs_array, min_distance_threshold, max_distance_threshold
            )
        else:
            euclidean_dists = np.linalg.norm(pairs_array[:, :3] - pairs_array[:, 3:6], axis=1)
            valid_mask = (euclidean_dists >= min_distance_threshold) & (euclidean_dists <= max_distance_threshold)
        
        pairs = pairs_array[valid_mask]
        
        if len(pairs) > 0:
            # Calculate distance statistics
            if NUMBA_AVAILABLE:
                pair_distances = np.array([_fast_distance_calculation(pair[:3], pair[3:6]) for pair in pairs])
            else:
                pair_distances = np.linalg.norm(pairs[:, :3] - pairs[:, 3:6], axis=1)
            
            print(f"Distance statistics: min={np.min(pair_distances):.2f}, "
                  f"max={np.max(pair_distances):.2f}, mean={np.mean(pair_distances):.2f}")
            
            # Further filtering based on distribution
            distance_threshold_low = np.percentile(pair_distances, 5)
            distance_threshold_high = np.percentile(pair_distances, 95)
            
            final_valid_mask = (pair_distances >= distance_threshold_low) & (pair_distances <= distance_threshold_high)
            pairs = pairs[final_valid_mask]
            
            print(f"After distance filtering: {len(pairs)} pairs retained")

        # Save results
        if save_txt and dataid is not None:
            os.makedirs(f"{self.root_dir}/results", exist_ok=True)
            out_path = f"{self.root_dir}/results/{dataid}_pairs.txt"
            np.savetxt(out_path, pairs, fmt='%.2f')
            print(f"[INFO] Saved NE-PM pairs to {out_path}")

        print(f"[INFO] NE-PM pairs generated: {len(pairs)}")
        return pairs

    @staticmethod
    def points_along_vector(p1, p2, n_slices):
        """
        Generate equally spaced division points between two points
        
        Creates n_slices+1 equally spaced points along the vector from p1 to p2, 
        including start and end points. These points are used to create boundary 
        surfaces for radial partitions.
        
        Parameters
        ----------
        p1 : numpy.ndarray
            Starting point coordinates, shape (3,)
        p2 : numpy.ndarray
            Ending point coordinates, shape (3,)
        n_slices : int
            Number of division segments
            
        Returns
        -------
        numpy.ndarray
            Equally spaced points array, shape (n_slices+1, 3)
            
        Examples
        --------
        >>> points = Partitioning.points_along_vector([0,0,0], [10,10,10], 4)
        >>> print(f"Generated {len(points)} division points")
        """
        return np.array([p1 + (p2 - p1) * i / n_slices for i in range(n_slices + 1)])

    def shell_partition(self, pairs, cell_mask=None):
        """
        Generate shell partition points from paired NE-PM points with adaptive partitioning
        
        Creates equally spaced partition points along the connection line for each 
        NE-PM point pair, with adaptive subdivision for better handling of irregular cell shapes.
        
        Parameters
        ----------
        pairs : numpy.ndarray
            Paired points array, shape (K, 6)
        cell_mask : numpy.ndarray, optional
            Binary mask indicating the active cell regions, shape (H, W, D)

        Returns
        -------
        list of numpy.ndarray
            Shell points list, each element is a coordinate array of all points in one shell
            
        Notes
        -----
        Optimized version with adaptive partitioning for irregular cell shapes
        """
        print(f"Creating {self.n_slices} adaptive shells from {len(pairs)} pairs...")
        
        all_points = []
        pair_distances = []
        
        # Calculate distances for adaptive partitioning
        for pair in pairs:
            ne_pt = pair[:3]
            pm_pt = pair[3:]
            distance = np.linalg.norm(ne_pt - pm_pt)
            pair_distances.append(distance)
        
        pair_distances = np.array(pair_distances)
        median_distance = np.median(pair_distances)
        
        print(f"Pair distance statistics: min={np.min(pair_distances):.2f}, "
              f"median={median_distance:.2f}, max={np.max(pair_distances):.2f}")
        
        for i, pair in enumerate(pairs):
            ne_pt = pair[:3]
            pm_pt = pair[3:]
            current_distance = pair_distances[i]
            
            # Adaptive subdivision based on distance
            if current_distance < median_distance * 0.5:
                # For very short distances (indentations), use fewer subdivisions
                adaptive_n_slices = max(3, self.n_slices // 2)
            elif current_distance > median_distance * 2.0:
                # For very long distances, use more subdivisions
                adaptive_n_slices = self.n_slices + 2
            else:
                # Normal distances use standard subdivisions
                adaptive_n_slices = self.n_slices
            
            # Generate adaptive subdivision points
            slice_pts = self.points_along_vector(ne_pt, pm_pt, adaptive_n_slices)
            all_points.append((slice_pts, adaptive_n_slices))

        # Organize points into shells with weighted contribution
        shell_points = [[] for _ in range(self.n_slices + 1)]
        
        for pts_info in all_points:
            slice_pts, adaptive_n_slices = pts_info
            
            # Map adaptive points to standard shell indices
            for i, pt in enumerate(slice_pts):
                # Calculate relative position (0 to 1)
                relative_pos = i / adaptive_n_slices
                
                # Map to standard shell index
                shell_idx = min(int(relative_pos * self.n_slices), self.n_slices)
                shell_points[shell_idx].append(pt)

        cleaned_shells = []
        for i, shell in enumerate(shell_points):
            if len(shell) == 0:
                cleaned_shells.append(np.array([]))
                continue
                
            shell_arr = np.unique(np.array(shell), axis=0)
            
            if cell_mask is not None:
                # Use vectorized operations to filter points
                shell_int = shell_arr.astype(int)
                valid_mask = (
                    (shell_int[:, 0] >= 0) & (shell_int[:, 0] < cell_mask.shape[0]) &
                    (shell_int[:, 1] >= 0) & (shell_int[:, 1] < cell_mask.shape[1]) &
                    (shell_int[:, 2] >= 0) & (shell_int[:, 2] < cell_mask.shape[2])
                )
                
                if np.any(valid_mask):
                    valid_shell_int = shell_int[valid_mask]
                    mask_values = cell_mask[valid_shell_int[:, 0], valid_shell_int[:, 1], valid_shell_int[:, 2]]
                    final_mask = mask_values > 0
                    shell_arr = shell_arr[valid_mask][final_mask]
                else:
                    shell_arr = np.array([])
            
            cleaned_shells.append(shell_arr)
            print(f"Shell {i}: {len(shell_arr)} points")

        print(f"Generated {len(cleaned_shells)} adaptive shells total")
        return cleaned_shells

    def _convert_shells_to_partition_mask_pure_pairs(self, shells, shape, pm_mask=None, ne_mask=None):
        """
        Convert shell point cloud data to continuous 3D partition masks using pure NE-PM pair information
        
        This method only uses the shell points generated from NE-PM pairs without any EDT guidance.
        It assigns partition IDs based on proximity to shell boundaries and interpolation between shells.
        
        Parameters
        ----------
        shells : list of numpy.ndarray
            List of shell point arrays, each containing 3D coordinates
        shape : tuple
            Shape of the output partition mask (Z, Y, X)
        pm_mask : numpy.ndarray, optional
            Original plasma membrane mask for defining cell boundaries
        ne_mask : numpy.ndarray, optional
            Original nuclear envelope mask for defining nuclear boundaries
            
        Returns
        -------
        numpy.ndarray
            3D partition mask where each voxel value indicates partition ID (0=background)
        """
        from scipy.ndimage import binary_fill_holes
        from scipy.spatial import cKDTree
        import numpy as np
        
        print(f"Converting {len(shells)} shells to partition mask...")
        
        # Create cell regions from masks
        if pm_mask is not None and ne_mask is not None:
            pm_binary = pm_mask > 0
            ne_binary = ne_mask > 0
            
            cell_interior = binary_fill_holes(pm_binary)
            nuclear_interior = binary_fill_holes(ne_binary)
            cytoplasm_region = cell_interior & ~nuclear_interior
            
            print(f"Cytoplasm region voxels: {np.sum(cytoplasm_region)}")
        else:
            print("Warning: No masks provided, using full volume")
            cytoplasm_region = np.ones(shape, dtype=bool)
        
        # Filter out empty shells and build KDTrees for non-empty shells
        valid_shells = []
        shell_trees = []
        shell_indices = []
        
        for i, shell in enumerate(shells):
            if len(shell) > 0:
                valid_shells.append(shell)
                shell_trees.append(cKDTree(shell))
                shell_indices.append(i + 1)  # Partition IDs start from 1
                print(f"Shell {i}: {len(shell)} points -> Partition {i + 1}")
        
        if len(valid_shells) < 2:
            print("Error: Need at least 2 valid shells for partitioning")
            return np.zeros(shape, dtype=int)
        
        print(f"Built KDTrees for {len(valid_shells)} valid shells")
        
        # Get cytoplasm coordinates
        cytoplasm_coords = np.column_stack(np.where(cytoplasm_region))
        print(f"Processing {len(cytoplasm_coords)} cytoplasm voxels...")
        
        partition_mask = np.zeros(shape, dtype=int)
        
        # Process in batches for memory efficiency
        batch_size = 50000
        n_batches = (len(cytoplasm_coords) + batch_size - 1) // batch_size
        
        for batch_idx in range(n_batches):
            start_idx = batch_idx * batch_size
            end_idx = min((batch_idx + 1) * batch_size, len(cytoplasm_coords))
            batch_coords = cytoplasm_coords[start_idx:end_idx]
            
            # Find distances to all shells for the batch
            batch_distances = []
            for tree in shell_trees:
                distances, _ = tree.query(batch_coords, k=1)
                batch_distances.append(distances)
            
            batch_distances = np.array(batch_distances).T  # Shape: (batch_size, n_shells)
            
            # Assign partitions based on corrected logic
            batch_partition_ids = np.zeros(len(batch_coords), dtype=int)
            
            for i, coord_distances in enumerate(batch_distances):
                # Find the closest shell
                closest_shell_idx = np.argmin(coord_distances)
                closest_shell_id = shell_indices[closest_shell_idx]
                
                # Check if closest shell is innermost (shell 1) or outermost (shell n)
                if closest_shell_id == 1:
                    # Closest to innermost shell -> directly assign to partition 1
                    batch_partition_ids[i] = 1
                elif closest_shell_id >= len(shell_indices):
                    # Closest to outermost shell -> directly assign to outermost partition
                    batch_partition_ids[i] = max(shell_indices)
                else:
                    # Point is closest to a middle shell -> need interpolation between shells
                    sorted_indices = np.argsort(coord_distances)
                    closest_shell_idx = sorted_indices[0]
                    second_closest_shell_idx = sorted_indices[1] if len(sorted_indices) > 1 else closest_shell_idx
                    
                    closest_dist = coord_distances[closest_shell_idx]
                    second_closest_dist = coord_distances[second_closest_shell_idx]
                    
                    if len(sorted_indices) == 1:
                        # Only one shell available, assign to its partition
                        batch_partition_ids[i] = shell_indices[closest_shell_idx]
                    else:
                        # Multiple shells: determine which partition region this voxel belongs to
                        shell_id_1 = shell_indices[closest_shell_idx]
                        shell_id_2 = shell_indices[second_closest_shell_idx]
                        
                        # Determine partition based on relative position between shells
                        total_dist = closest_dist + second_closest_dist
                        if total_dist > 0:
                            # Relative position: closer to which shell
                            relative_pos = second_closest_dist / total_dist
                            
                            # If significantly closer to first shell, assign to its partition
                            if relative_pos > 0.6:  # >60% means much closer to first shell
                                batch_partition_ids[i] = shell_id_1
                            # If significantly closer to second shell, assign to its partition  
                            elif relative_pos < 0.4:  # <40% means much closer to second shell
                                batch_partition_ids[i] = shell_id_2
                            else:
                                # In between -> assign to the partition between these shells
                                # Use the smaller shell ID (closer to nucleus)
                                batch_partition_ids[i] = min(shell_id_1, shell_id_2)
                        else:
                            # Fallback: use closest shell's partition
                            batch_partition_ids[i] = shell_indices[closest_shell_idx]
            
            # Assign batch results to partition mask
            for i, coord in enumerate(batch_coords):
                z, y, x = coord
                partition_mask[z, y, x] = batch_partition_ids[i]
            
            if (batch_idx + 1) % 20 == 0:
                print(f"Processed batch {batch_idx + 1}/{n_batches}")
        
        # Post-processing: Fill gaps and smooth boundaries
        print("Post-processing partition mask...")
        
        # Fill small gaps within partitions
        for partition_id in shell_indices:
            partition_region = partition_mask == partition_id
            if np.sum(partition_region) > 0:
                # Small morphological closing to fill gaps
                from scipy.ndimage import binary_closing
                filled_region = binary_closing(partition_region, iterations=1)
                partition_mask[filled_region & cytoplasm_region] = partition_id
        
        # Verify partition results
        unique_partitions = np.unique(partition_mask)
        print("Final partition distribution:")
        for partition_id in unique_partitions:
            if partition_id != 0:
                count = np.sum(partition_mask == partition_id)
                print(f"Partition {partition_id}: {count} voxels")
        
        print(f"Successfully created {len(unique_partitions)-1} pure pair-based partitions")
        return partition_mask

    def create_nepm_radial_partitions_with_edt(self, ne_edge, pm_edge, shape, n_slices=None, pm_mask=None, ne_mask=None):
        """
        Optimized version of radial partition generation
        
        Parameters
        ----------
        ne_edge : numpy.ndarray
            Nuclear envelope boundary coordinates array, shape (N, 3)
        pm_edge : numpy.ndarray
            Plasma membrane boundary coordinates array, shape (M, 3)
        shape : tuple
            Shape of original image volume (Z, Y, X)
        n_slices : int, optional
            Number of radial partitions, defaults to self.n_slices
        pm_mask : numpy.ndarray, optional
            Original plasma membrane mask for improved accuracy
        ne_mask : numpy.ndarray, optional
            Original nuclear envelope mask for improved accuracy
            
        Returns
        -------
        numpy.ndarray
            3D partition mask array where each voxel value indicates partition ID (0=background)
            
        Notes
        -----
        Uses local EDT-guided pair matching for better handling of irregular cell shapes
        """
        if n_slices is None:
            n_slices = self.n_slices
            
        print(f"Creating {n_slices} radial partitions using optimized algorithms...")
        
        # Calculate cell center
        center = np.mean(ne_edge, axis=0)
        
        # Generate optimized NE-PM pairs
        pairs = self.find_ne_pm_pairs(center, ne_edge, pm_edge, step=8, save_txt=False)
        print(f"Generated {len(pairs)} NE-PM pairs")
        
        # Set correct n_slices
        original_n_slices = self.n_slices
        self.n_slices = n_slices - 1
        
        # Use optimized shell generation
        if len(pairs) > 1000:
            shells = self.shell_partition_optimized(pairs, cell_mask=pm_mask)
        else:
            shells = self.shell_partition(pairs, cell_mask=pm_mask)
        
        print(f"Generated {len(shells)} shells for {n_slices} partitions")
        
        # Restore original n_slices
        self.n_slices = original_n_slices
        
        # Use optimized partition mask conversion
        partition_mask = self._convert_shells_to_partition_mask_optimized(shells, shape, pm_mask, ne_mask)
        
        return partition_mask

    def _convert_shells_to_partition_mask_optimized(self, shells, shape, pm_mask=None, ne_mask=None):
        """
        Optimized version of partition mask conversion with corrected shell assignment logic
        """
        from scipy.ndimage import binary_fill_holes, distance_transform_edt, binary_erosion
        from scipy.spatial import cKDTree
        
        print(f"Converting {len(shells)} shells to partition mask using corrected assignment logic...")
        
        # Create cytoplasm region mask
        if pm_mask is not None and ne_mask is not None:
            pm_binary = pm_mask > 0
            ne_binary = ne_mask > 0
            
            cell_interior = binary_fill_holes(pm_binary)
            nuclear_interior = binary_fill_holes(ne_binary)
            
            if cell_interior is not None and nuclear_interior is not None:
                cytoplasm_region = np.asarray(cell_interior, dtype=bool) & (~np.asarray(nuclear_interior, dtype=bool))
                
                pm_eroded = binary_erosion(pm_binary, iterations=2)
                pm_boundary = np.asarray(pm_binary, dtype=bool) & (~np.asarray(pm_eroded, dtype=bool))
                
                ne_eroded = binary_erosion(ne_binary, iterations=1)
                ne_boundary = np.asarray(ne_binary, dtype=bool) & (~np.asarray(ne_eroded, dtype=bool))
                
                print(f"Cytoplasm region voxels: {np.sum(cytoplasm_region)}")
            else:
                print("Error: Failed to create interior regions")
                return np.zeros(shape, dtype=int)
        else:
            print("Error: Original masks are required for optimized method")
            return np.zeros(shape, dtype=int)
        
        # Get cytoplasm coordinates
        cytoplasm_coords = np.column_stack(np.where(cytoplasm_region))
        print(f"Processing {len(cytoplasm_coords)} cytoplasm voxels...")
        
        n_partitions = len([s for s in shells if len(s) > 0])
        partition_mask = np.zeros(shape, dtype=int)
        
        if n_partitions < 2:
            print("Error: Need at least 2 valid shells for partitioning")
            return partition_mask
        
        # Build KDTrees for boundary detection
        ne_boundary_coords = np.column_stack(np.where(ne_boundary))
        pm_boundary_coords = np.column_stack(np.where(pm_boundary))
        
        if len(ne_boundary_coords) > 0 and len(pm_boundary_coords) > 0:
            ne_tree = cKDTree(ne_boundary_coords)
            pm_tree = cKDTree(pm_boundary_coords)
            print("Built boundary KDTrees for direct assignment")
        else:
            print("Warning: Cannot build boundary trees, falling back to shell-based method")
            ne_tree = pm_tree = None
        
        # Process in batches for memory efficiency
        batch_size = 50000
        n_batches = (len(cytoplasm_coords) + batch_size - 1) // batch_size
        
        # Define distance thresholds for direct assignment
        direct_assignment_threshold = 3.0  # pixels
        
        for batch_idx in range(n_batches):
            start_idx = batch_idx * batch_size
            end_idx = min((batch_idx + 1) * batch_size, len(cytoplasm_coords))
            batch_coords = cytoplasm_coords[start_idx:end_idx]
            
            batch_partition_ids = np.zeros(len(batch_coords), dtype=int)
            direct_assignment_mask = np.zeros(len(batch_coords), dtype=bool)
            
            # Step 1: Direct assignment for points very close to boundaries
            if ne_tree is not None and pm_tree is not None:
                # Distance to nuclear boundary
                ne_distances, _ = ne_tree.query(batch_coords, k=1)
                # Distance to plasma membrane boundary  
                pm_distances, _ = pm_tree.query(batch_coords, k=1)
                
                # Direct assignment to innermost shell (shell 1) for points very close to nucleus
                inner_mask = (ne_distances <= direct_assignment_threshold) & (ne_distances < pm_distances)
                batch_partition_ids[inner_mask] = 1
                direct_assignment_mask[inner_mask] = True
                
                # Direct assignment to outermost shell for points very close to PM
                outer_mask = (pm_distances <= direct_assignment_threshold) & (pm_distances < ne_distances)
                batch_partition_ids[outer_mask] = n_partitions
                direct_assignment_mask[outer_mask] = True
                
                print(f"Batch {batch_idx}: Direct assignment - {np.sum(inner_mask)} to inner, {np.sum(outer_mask)} to outer")
            
            # Step 2: Shell-based assignment for remaining points
            remaining_mask = ~direct_assignment_mask
            remaining_coords = batch_coords[remaining_mask]
            
            if len(remaining_coords) > 0:
                # Find distances to all shells for remaining points
                shell_distances = []
                valid_shell_indices = []
                
                for i, shell in enumerate(shells):
                    if len(shell) > 0:
                        shell_tree = cKDTree(shell)
                        distances, _ = shell_tree.query(remaining_coords, k=1)
                        shell_distances.append(distances)
                        valid_shell_indices.append(i + 1)  # Partition IDs start from 1
                
                if len(shell_distances) > 0:
                    shell_distances = np.array(shell_distances).T  # Shape: (n_remaining, n_shells)
                    
                    # For each remaining point, find the closest shell
                    closest_shell_indices = np.argmin(shell_distances, axis=1)
                    
                    # Assign to closest shell's partition
                    for i, shell_idx in enumerate(closest_shell_indices):
                        batch_partition_ids[remaining_mask][i] = valid_shell_indices[shell_idx]
            
            # Assign batch results to partition mask
            for i, coord in enumerate(batch_coords):
                if batch_partition_ids[i] > 0:
                    z, y, x = coord
                    partition_mask[z, y, x] = batch_partition_ids[i]
            
            if (batch_idx + 1) % 20 == 0:
                print(f"Processed batch {batch_idx + 1}/{n_batches}")
        
        # Post-processing: Ensure proper shell ordering and fill gaps
        print("Post-processing partition mask...")
        
        # Verify and correct partition ordering
        unique_partitions = np.unique(partition_mask)
        valid_partitions = unique_partitions[unique_partitions > 0]
        
        print("Partition distribution after corrected assignment:")
        for partition_id in valid_partitions:
            count = np.sum(partition_mask == partition_id)
            print(f"Partition {partition_id}: {count} voxels")
        
        # Fill small gaps within partitions
        from scipy.ndimage import binary_closing
        for partition_id in valid_partitions:
            partition_region = partition_mask == partition_id
            if np.sum(partition_region) > 0:
                filled_region = binary_closing(partition_region, iterations=1)
                combined_region = np.asarray(filled_region, dtype=bool) & cytoplasm_region
                partition_mask[combined_region] = partition_id
        
        print(f"Successfully created {len(valid_partitions)} corrected partitions")
        return partition_mask

    def _parallel_shell_processing(self, pairs_chunk, n_slices):
        """
        Parallel processing worker function for shell generation
        """
        shell_points = [[] for _ in range(n_slices + 1)]
        
        for pair in pairs_chunk:
            ne_pt = pair[:3]
            pm_pt = pair[3:]
            
            # Generate subdivision points
            slice_pts = self.points_along_vector(ne_pt, pm_pt, n_slices)
            
            for i, pt in enumerate(slice_pts):
                shell_points[i].append(pt)
        
        return shell_points

    def shell_partition_optimized(self, pairs, cell_mask=None):
        """
        Optimized version of shell partition generation with parallel processing support
        """
        print(f"Creating {self.n_slices + 1} shells from {len(pairs)} pairs ...")
        
        # Initialize shell_points for all code paths
        shell_points = [[] for _ in range(self.n_slices + 1)]
        
        # Use parallel processing if available and data is large
        try:
            if len(pairs) > 5000:
                from joblib import Parallel, delayed
                print("Using parallel processing for shell generation...")
                
                # Chunk processing
                n_cores = min(4, os.cpu_count() or 1)  # Limit cores to avoid over-parallelization
                chunk_size = len(pairs) // n_cores + 1
                pairs_chunks = [pairs[i:i+chunk_size] for i in range(0, len(pairs), chunk_size)]
                
                # Parallel processing
                results = Parallel(n_jobs=n_cores)(
                    delayed(self._parallel_shell_processing)(chunk, self.n_slices) 
                    for chunk in pairs_chunks
                )
                
                # Merge results
                if results is not None:
                    for result in results:
                        if result is not None:
                            for i, shell_pts in enumerate(result):
                                shell_points[i].extend(shell_pts)
            else:
                raise ImportError("Data size too small for parallel processing")
                
        except (ImportError, ModuleNotFoundError):
            print("Using sequential processing for shell generation...")
            # Original sequential processing logic
            for pair in pairs:
                ne_pt = pair[:3]
                pm_pt = pair[3:]
                slice_pts = self.points_along_vector(ne_pt, pm_pt, self.n_slices)
                
                for i, pt in enumerate(slice_pts):
                    shell_points[i].append(pt)

        # Cleanup and deduplication
        cleaned_shells = []
        for i, shell in enumerate(shell_points):
            if len(shell) == 0:
                cleaned_shells.append(np.array([]))
                continue
                
            shell_arr = np.unique(np.array(shell), axis=0)
            
            if cell_mask is not None:
                # Vectorized point validation
                shell_int = shell_arr.astype(int)
                valid_mask = (
                    (shell_int[:, 0] >= 0) & (shell_int[:, 0] < cell_mask.shape[0]) &
                    (shell_int[:, 1] >= 0) & (shell_int[:, 1] < cell_mask.shape[1]) &
                    (shell_int[:, 2] >= 0) & (shell_int[:, 2] < cell_mask.shape[2])
                )
                
                if np.any(valid_mask):
                    valid_shell_int = shell_int[valid_mask]
                    mask_values = cell_mask[valid_shell_int[:, 0], valid_shell_int[:, 1], valid_shell_int[:, 2]]
                    final_mask = mask_values > 0
                    shell_arr = shell_arr[valid_mask][final_mask]
                else:
                    shell_arr = np.array([])
            
            cleaned_shells.append(shell_arr)
            print(f"Shell {i}: {len(shell_arr)} points")

        print(f"Generated {len(cleaned_shells)} optimized shells total")
        return cleaned_shells

    def extract_partition_coordinates(self, partition_mask, sampling_density=0.1):
        """
        Extract coordinate points from continuous partition mask for saving in XVG format
        
        Parameters
        ----------
        partition_mask : numpy.ndarray
            3D partition mask array
        sampling_density : float, default=0.1
            Sampling density, 0.1 means sampling 10% of points
            
        Returns
        -------
        list of numpy.ndarray
            List of coordinate points for each partition
        """
        unique_partitions = np.unique(partition_mask)
        shell_coords = []
        
        for partition_id in unique_partitions:
            if partition_id == 0:  # Skip background
                continue
                
            # Get all coordinates of current partition
            coords = np.column_stack(np.where(partition_mask == partition_id))
            
            if len(coords) > 0:
                # Sample according to sampling density
                n_samples = max(1, int(len(coords) * sampling_density))
                if n_samples < len(coords):
                    indices = np.random.choice(len(coords), n_samples, replace=False)
                    sampled_coords = coords[indices]
                else:
                    sampled_coords = coords
                
                shell_coords.append(sampled_coords)
                print(f"Partition {partition_id}: {len(coords)} total points, {len(sampled_coords)} sampled")
            else:
                shell_coords.append(np.array([]))
        
        return shell_coords

    def save_partition_coords_to_xvg(self, partition_coords, dataid, output_dir):
        """
        Save partition coordinates to XVG format file
        
        Parameters
        ----------
        partition_coords : list of numpy.ndarray
            List of partition coordinates
        dataid : str
            Data identifier
        output_dir : str
            Output directory
            
        Returns
        -------
        str
            Output directory path
        """
        xvg_dir = output_dir
        os.makedirs(xvg_dir, exist_ok=True)
        
        # Save combined XVG file
        combined_xvg_path = os.path.join(xvg_dir, f'{dataid}_partition_coords.xvg')
        with open(combined_xvg_path, 'w') as f:
            # XVG header
            f.write("# XVG file generated by iPA - Partition Coordinates\n")
            f.write(f"# All partition coordinates from data ID {dataid}\n")
            f.write("# Contains 3D coordinates (z, y, x) and partition index\n")
            f.write("@    title \"Partition Coordinates\"\n")
            f.write("@    xaxis  label \"Z\"\n")
            f.write("@    yaxis  label \"Y\"\n")
            f.write("@TYPE xy\n")
            
            # Define legends for each partition
            for i, coords in enumerate(partition_coords):
                if len(coords) > 0:
                    f.write(f"@    s{i} legend \"Partition {i+1}\"\n")
            
            # Write data points
            for i, coords in enumerate(partition_coords):
                for j, point in enumerate(coords):
                    if len(point) == 2:  # 2D points
                        f.write(f"{point[0]:.3f} {point[1]:.3f} 0.000 {i+1}\n")
                    elif len(point) == 3:  # 3D points
                        f.write(f"{point[0]:.3f} {point[1]:.3f} {point[2]:.3f} {i+1}\n")
                    else:
                        raise ValueError(f"Unexpected point dimension: {len(point)}")

        print(f"[INFO] Saved partition coordinates as XVG to {combined_xvg_path}")
        
        return xvg_dir

    def create_nepm_radial_partitions(self, ne_edge, pm_edge, shape, n_slices=None, pm_mask=None, ne_mask=None):
        """
        Create radial partitions using pure NE-PM pair information without EDT guidance
        
        This method generates radial partitions based solely on the shell points created from
        NE-PM pairs, without using Euclidean Distance Transform for edge refinement.
        It's more straightforward and computationally efficient than the EDT-guided version.
        
        Parameters
        ----------
        ne_edge : numpy.ndarray
            Nuclear envelope boundary coordinates array, shape (N, 3)
        pm_edge : numpy.ndarray
            Plasma membrane boundary coordinates array, shape (M, 3)
        shape : tuple
            Shape of original image volume (Z, Y, X)
        n_slices : int, optional
            Number of radial partitions, defaults to self.n_slices
        pm_mask : numpy.ndarray, optional
            Original plasma membrane mask for defining cell boundaries
        ne_mask : numpy.ndarray, optional
            Original nuclear envelope mask for defining nuclear boundaries
            
        Returns
        -------
        numpy.ndarray
            3D partition mask array where each voxel value indicates partition ID (0=background)
            
        Notes
        -----
        This pure pair-based method is recommended when:
        - You want faster computation without EDT overhead
        - The cell shape is relatively regular
        - You prefer simpler, more interpretable partitioning logic
        
        Examples
        --------
        >>> partition_mask = partitioner.create_nepm_radial_partitions(
        ...     ne_coords, pm_coords, (100, 200, 200), n_slices=8
        ... )
        >>> print(f"Created partitions with shape: {partition_mask.shape}")
        """
        if n_slices is None:
            n_slices = self.n_slices
            
        print(f"Creating {n_slices} radial partitions using pure NE-PM pair method...")
        
        # Calculate cell center
        center = np.mean(ne_edge, axis=0)
        
        # Generate NE-PM pairs
        pairs = self.find_ne_pm_pairs(center, ne_edge, pm_edge, step=8, save_txt=False)
        print(f"Generated {len(pairs)} NE-PM pairs")
        
        # Store original n_slices and set for shell generation
        original_n_slices = self.n_slices
        self.n_slices = n_slices - 1
        
        # Generate shells from pairs
        if len(pairs) > 1000:
            shells = self.shell_partition_optimized(pairs, cell_mask=pm_mask)
        else:
            shells = self.shell_partition(pairs, cell_mask=pm_mask)
        
        print(f"Generated {len(shells)} shells for {n_slices} partitions")
        
        # Restore original n_slices
        self.n_slices = original_n_slices
        
        # Convert shells to partition mask using pure pair method
        partition_mask = self._convert_shells_to_partition_mask_pure_pairs(shells, shape, pm_mask, ne_mask)
        
        return partition_mask

    # Keep backward compatibility alias
    def create_shell_based_partitions(self, ne_edge, pm_edge, shape, n_slices=None, pm_mask=None, ne_mask=None):
        """
        Deprecated: Use create_nepm_radial_partitions_with_edt instead.
        
        This method is kept for backward compatibility.
        """
        print("Warning: create_shell_based_partitions is deprecated. Use create_nepm_radial_partitions_with_edt instead.")
        return self.create_nepm_radial_partitions_with_edt(ne_edge, pm_edge, shape, n_slices, pm_mask, ne_mask)
