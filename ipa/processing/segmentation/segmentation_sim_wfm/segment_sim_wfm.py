
from skimage import filters, morphology
from scipy import ndimage
from skimage.morphology import skeletonize_3d
from scipy import ndimage
import numpy as np
import tifffile
import matplotlib.pyplot as plt
from pathlib import Path




def remove_edges(img, edge=10):
    """
    Remove edge voxels from 3D image to avoid edge artifacts
    """
    img[:edge, :, :] = False
    img[-edge:, :, :] = False
    img[:, :edge, :] = False
    img[:, -edge:, :] = False
    img[:, :, :edge] = False
    img[:, :, -edge:] = False
    return img


def connect_skeleton_segments(skeleton, max_distance=3):
    """
    Connect nearby skeleton segments by finding close endpoints and bridging them
    
    Args:
        skeleton: 3D binary array of skeleton
        max_distance: maximum distance to connect segments
        
    Returns:
        connected_skeleton: skeleton with connected segments
    """
    from scipy.spatial.distance import cdist
    
    # Find skeleton endpoints (pixels with only one neighbor)
    skeleton_coords = np.where(skeleton)
    skeleton_points = np.array(list(zip(*skeleton_coords)))
    
    if len(skeleton_points) == 0:
        return skeleton
    
    endpoints = []
    for point in skeleton_points:
        z, y, x = point
        # Check 26-neighborhood
        neighbors = 0
        for dz in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                for dx in [-1, 0, 1]:
                    if dz == 0 and dy == 0 and dx == 0:
                        continue
                    nz, ny, nx = z + dz, y + dy, x + dx
                    if (0 <= nz < skeleton.shape[0] and 
                        0 <= ny < skeleton.shape[1] and 
                        0 <= nx < skeleton.shape[2] and
                        skeleton[nz, ny, nx]):
                        neighbors += 1
        
        if neighbors == 1:  # This is an endpoint
            endpoints.append(point)
    
    if len(endpoints) < 2:
        return skeleton
    
    # Calculate distances between all endpoints
    endpoints = np.array(endpoints)
    distances = cdist(endpoints, endpoints)
    
    # Connect nearby endpoints
    connected_skeleton = skeleton.copy()
    for i in range(len(endpoints)):
        for j in range(i + 1, len(endpoints)):
            if distances[i, j] <= max_distance:
                # Draw line between endpoints
                p1, p2 = endpoints[i], endpoints[j]
                
                # Simple line drawing in 3D
                steps = int(np.ceil(distances[i, j]))
                for step in range(steps + 1):
                    t = step / max(steps, 1)
                    point = p1 + t * (p2 - p1)
                    z, y, x = np.round(point).astype(int)
                    if (0 <= z < skeleton.shape[0] and 
                        0 <= y < skeleton.shape[1] and 
                        0 <= x < skeleton.shape[2]):
                        connected_skeleton[z, y, x] = True
    
    return connected_skeleton

def skeletonize_organelle(image_matrix):
    """
    Extract 3D skeleton using global thresholding and 3D skeletonization
    
    Args:
        image_matrix: 3D numpy array, input SIM image data
        max_distance: float, maximum distance for connecting skeleton segments

    Returns:
        skeleton_final: 3D numpy array, final skeletonized result
    """
    print(f"Input image stats: min={np.min(image_matrix)}, max={np.max(image_matrix)}, mean={np.mean(image_matrix):.2f}")

    # Gaussian smoothing for noise reduction
    smoothed = filters.gaussian(image_matrix, sigma=1.0)
    
    # Use higher threshold to keep only bright internal regions, avoiding edges
    threshold = np.mean(smoothed) + 1.2 * np.std(smoothed)
    binary_image = smoothed > threshold
    print(f"Global threshold: {threshold:.2f}")

    print(f"Foreground pixel ratio after thresholding: {np.sum(binary_image) / binary_image.size * 100:.2f}%")

    # Morphological erosion to further remove shape edges, keep only central regions
    binary_image = morphology.binary_erosion(binary_image, morphology.ball(1))
    
    # Remove small objects
    binary_image = morphology.remove_small_objects(binary_image, min_size=200)
    
    # Slight dilation to recover some structures
    binary_image = morphology.binary_dilation(binary_image, morphology.ball(1))

    # Extract 3D skeleton
    skeleton = skeletonize_3d(binary_image)

    # Minimal cleaning, only remove tiny skeleton fragments
    skeleton_cleaned = morphology.remove_small_objects(skeleton, min_size=30)
    
    # Connect segmented skeleton parts
    skeleton_connected = connect_skeleton_segments(skeleton_cleaned, max_distance=5)
    
    # Clean again after connection to remove possible small fragments
    skeleton_final = morphology.remove_small_objects(skeleton_connected, min_size=20)

    print(f"Skeleton voxel count after cleaning and connecting: {np.sum(skeleton_final)}")

    return skeleton_final



def segment_sphere_like_organelle(image_3d, threshold=None, min_size=100):
    """
    Simple segmentation of spherical structures in 3D images
    
    Args:
        image_3d: 3D numpy array, input image
        threshold: float or None, threshold value (if None, use Otsu automatic threshold)
        min_size: int, filter out connected components smaller than this voxel count
    
    Returns:
        labeled_spheres: 3D numpy integer array, connected component labeling result
        num_spheres: number of segmented spheres
    """
    # Gaussian filtering for noise reduction
    smoothed = filters.gaussian(image_3d, sigma=1)
    
    # Threshold segmentation
    if threshold is None:
        threshold = filters.threshold_otsu(smoothed)
    binary = smoothed > threshold
    
    # Morphological opening for noise removal, closing to connect small objects
    binary = morphology.binary_opening(binary, morphology.ball(1))
    binary = morphology.binary_closing(binary, morphology.ball(2))
    
    # Connected component labeling
    labeled_spheres, num_spheres = ndimage.label(binary)
    
    # Filter small fragments
    sizes = ndimage.sum(binary, labeled_spheres, range(num_spheres + 1))
    mask_sizes = sizes >= min_size
    mask_sizes[0] = 0  # Do not keep background
    cleaned = mask_sizes[labeled_spheres]
    
    # Re-label filtered connected components
    labeled_spheres, num_spheres = ndimage.label(cleaned)
    
    return labeled_spheres, num_spheres



def segment_cell_shape(cell_channel_3d, threshold=None, min_size=1000):
    """
    Basic cell shape segmentation
    
    Args:
        cell_channel_3d: 3D numpy array, input cell channel image
        threshold: float or None, threshold value
        min_size: int, minimum size for filtering
        
    Returns:
        labeled_cell: labeled cell mask
    """
    smoothed = filters.gaussian(cell_channel_3d, sigma=2)
    if threshold is None:
        threshold = filters.threshold_otsu(smoothed)
    binary = smoothed > threshold

    binary = morphology.binary_opening(binary, morphology.ball(3))
    binary = morphology.binary_closing(binary, morphology.ball(5))

    binary = ndimage.binary_fill_holes(binary)

    labeled_cell, num_cells = ndimage.label(binary)
    sizes = ndimage.sum(binary, labeled_cell, range(num_cells + 1))
    max_label = sizes.argmax()
    cleaned = labeled_cell == max_label
    labeled_cell, num_cells = ndimage.label(cleaned)

    return labeled_cell


def segment_nucleus(nucleus_channel_3d, threshold=None):
    """
    Basic nucleus segmentation

    Args:
        nucleus_channel_3d: 3D numpy array, input nucleus channel image
        threshold: float or None, threshold value (if None, use Otsu automatic threshold)

    Returns:
        labeled_nuclei: 3D numpy integer array, connected component labeling result
    """

    print(f"Input nucleus image stats: min={np.min(nucleus_channel_3d)}, max={np.max(nucleus_channel_3d)}, mean={np.mean(nucleus_channel_3d):.2f}")

    smoothed = filters.gaussian(nucleus_channel_3d, sigma=0.8)

    if threshold is None:
        threshold = filters.threshold_otsu(smoothed)
    print(f"Nucleus threshold: {threshold:.2f}")

    binary = smoothed > threshold
    print(f"Foreground ratio after thresholding: {np.sum(binary) / binary.size * 100:.2f}%")

    # 使用形态学开闭操作保持连贯性
    binary = morphology.binary_opening(binary, morphology.ball(1))
    binary = morphology.binary_closing(binary, morphology.ball(2))

    binary = ndimage.binary_fill_holes(binary)

    # 不做分水岭分割，直接把所有连通区域合并成一个====>只用一个label
    labeled_nuclei = np.zeros_like(binary, dtype=np.uint16)
    labeled_nuclei[binary] = 1
    num_nuclei = 1 if np.any(binary) else 0

    print(f"Number of nuclei found (merged as one): {num_nuclei}")

    return labeled_nuclei




# def main():
    # Data file path, change to actual file path
    # data_path = Path("d:/Gitspace/ipa_full/iPA/data/sim_images/20220909_30-2-1-SIM_raw_Actin.tif")

    # if not data_path.exists():
    #     print(f"Error: File not found {data_path}")
    #     return

    # print("Loading image data...")
    # image_data = tifffile.imread(str(data_path))
    # print(f"Image shape: {image_data.shape}")
    # print(f"Data type: {image_data.dtype}")

    # print("Performing 3D skeletonization...")
    # skeleton_result = skeletonize_image(image_data)

    # print("Generating visualization...")
    # fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # mid_slice = image_data.shape[0] // 2
    # mid_slice = 7 

    # # Original image three views
    # axes[0, 0].imshow(image_data[mid_slice], cmap='gray')
    # axes[0, 0].set_title(f'Original Image - Z slice {mid_slice}')
    # axes[0, 0].axis('off')

    # axes[0, 1].imshow(image_data[:, image_data.shape[1]//2], cmap='gray')
    # axes[0, 1].set_title('Original Image - Y slice')
    # axes[0, 1].axis('off')

    # axes[0, 2].imshow(image_data[:, :, image_data.shape[2]//2], cmap='gray')
    # axes[0, 2].set_title('Original Image - X slice')
    # axes[0, 2].axis('off')

    # # Skeleton result three views
    # axes[1, 0].imshow(skeleton_result[mid_slice], cmap='gray')
    # axes[1, 0].set_title(f'Skeleton Result - Z slice {mid_slice}')
    # axes[1, 0].axis('off')

    # axes[1, 1].imshow(skeleton_result[:, skeleton_result.shape[1]//2], cmap='gray')
    # axes[1, 1].set_title('Skeleton Result - Y slice')
    # axes[1, 1].axis('off')

    # axes[1, 2].imshow(skeleton_result[:, :, skeleton_result.shape[2]//2], cmap='gray')
    # axes[1, 2].set_title('Skeleton Result - X slice')
    # axes[1, 2].axis('off')

    # plt.tight_layout()
    # plt.show()

    # # Ask whether to save results
    # save_option = input("Save skeleton results? (y/n): ").lower().strip()
    # if save_option in ['y', 'yes']:
    #     output_dir = Path("d:/Gitspace/ipa_full/iPA/data/processed")
    #     output_dir.mkdir(exist_ok=True)
    #     output_path = output_dir / "skeleton_result.tif"
    #     tifffile.imwrite(str(output_path), skeleton_result.astype(np.uint8) * 255)
    #     print(f"Skeleton result saved to: {output_path}")

    #     vis_path = output_dir / "skeleton_visualization.png"
    #     fig.savefig(str(vis_path), dpi=300, bbox_inches='tight')
    #     print(f"Visualization saved to: {vis_path}")


# if __name__ == "__main__":
#     main()
