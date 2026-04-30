"""
Advanced Blob Separation for SXT ISG Instance Segmentation.
Adapted from post_processing_tool_for_instance_segmentation_on_SXT/vesicle_seg/blob_separation_v6_batch.py
"""
import numpy as np
import copy
from scipy import ndimage as ndi
from scipy.spatial.distance import cdist
from skimage.feature import peak_local_max
from scipy.ndimage import gaussian_filter


def discrete_sphere(max_r, min_r):
    """Generate discrete sphere radius and volume lookup table."""
    cube = np.ones([int(max_r * 2 + 1 + 2), int(max_r * 2 + 1 + 2), int(max_r * 2 + 1 + 2)])
    coords_ball = np.where(cube == 1)
    x_ball = int(np.mean(coords_ball[0]))
    y_ball = int(np.mean(coords_ball[1]))
    z_ball = int(np.mean(coords_ball[2]))
    cent_coord = np.array([[x_ball, y_ball, z_ball]])
    coords_allball = np.array(coords_ball).T
    dist_lst = np.array(cdist(cent_coord, coords_allball, metric='euclidean')[0])
    
    r_lst = np.array(sorted(list(set(dist_lst)), reverse=False))
    r_lst = r_lst[np.where(r_lst > min_r)]
    r_volume_lst = [len(np.where(dist_lst <= r)[0]) for r in r_lst]
    
    return np.concatenate((np.array([r_lst]), np.array([r_volume_lst])), axis=0)


def get_rad(centcoord, mask, r_seq, theta):
    """Calculate the optimal radius for a blob center."""
    r_lst = r_seq[0]
    r_volume = r_seq[1]
    
    all_coords = np.array(np.where(mask != 0)).T
    if len(all_coords) == 0:
        return 0
        
    dist = cdist(np.array([centcoord]), all_coords)[0]
    mask_volumes = [len(np.where(dist <= d)[0]) for d in r_lst]
    ratiolst = [mask_volumes[i] / r_volume[i] for i in range(len(r_lst))]
    
    if np.max(ratiolst) < theta:
        return 0
    else:
        idx = np.where(np.array(ratiolst) >= theta)[0].shape[0]
        # Ensure index is within bounds
        if idx >= len(r_lst):
            idx = len(r_lst) - 1
        r = r_lst[idx]
        return r


def blob_fit(mask, mrc, min_r=1.5, theta=0.8, check_=False):
    """
    Separate touching blobs in a 3D mask using intensity information from mrc.
    
    Args:
        mask: Binary or labeled 3D mask
        mrc: Raw intensity data (same shape as mask)
        min_r: Minimum radius for a blob
        theta: Threshold for sphere fitting ratio
        check_: If True, enable debug plotting (not implemented here)
        
    Returns:
        Labeled instance mask
    """
    mask = mask.copy()
    mask[np.where(mask != 0)] = 1
    
    sigmalst = list(np.arange(1, 10, 1))
    blobs_lst = []
    
    for sigma in sigmalst:
        filtered_mrc = gaussian_filter(mrc, sigma)
        temp_blobs_lst = peak_local_max(filtered_mrc, footprint=np.ones((3,) * (mask.ndim)), exclude_border=True)
        blobs_lst.extend(temp_blobs_lst)
    
    if not blobs_lst:
        if np.sum(mask) > 15:
            labeled, _ = ndi.label(mask)
            return labeled
        else:
            return np.zeros_like(mask)
    
    coord_intensity_lst = []
    coord_edt_lst = []
    edtmap = ndi.distance_transform_edt(ndi.binary_dilation(mask))
    
    for coord in blobs_lst:
        # Ensure coordinates are within bounds
        if (0 <= coord[0] < mrc.shape[0] and 
            0 <= coord[1] < mrc.shape[1] and 
            0 <= coord[2] < mrc.shape[2]):
            coord_intensity_lst.append(mrc[coord[0]][coord[1]][coord[2]])
            coord_edt_lst.append(edtmap[coord[0]][coord[1]][coord[2]])
        else:
            coord_intensity_lst.append(0)
            coord_edt_lst.append(0)

    # Sort coords by intensity (priority to brighter spots)
    indices = np.argsort(coord_intensity_lst)[::-1]
    coord_intensity_lst = [coord_intensity_lst[i] for i in indices]
    coord_edt_lst = [coord_edt_lst[i] for i in indices]
    blobs_lst = [blobs_lst[i] for i in indices]
    
    max_r = np.max(coord_edt_lst) // 1 + 2 + 1 if coord_edt_lst else 5
    if max_r < min_r:
        max_r = min_r + 1
        
    r_seq = discrete_sphere(max_r, min_r)
    
    final_cent_coord_lst = []
    final_cent_rad_lst = []
    iterated_mask = copy.deepcopy(mask)
    
    for idx, coord in enumerate(blobs_lst):
        if coord_edt_lst[idx] == 0:
            continue
            
        if len(final_cent_coord_lst) != 0:
            dist = cdist(np.array([coord]), np.array(final_cent_coord_lst))[0]
            diff = dist - np.array(final_cent_rad_lst)
            if np.any(diff < 0):
                continue
                
        r = get_rad(coord, mask, r_seq, theta)
        
        if r > 0:
            final_cent_coord_lst.append(coord)
            final_cent_rad_lst.append(r)
            
            temp_edt = np.ones_like(iterated_mask)
            temp_edt[coord[0]][coord[1]][coord[2]] = 0
            temp_edt = ndi.distance_transform_edt(temp_edt)
            iterated_mask[np.where(temp_edt < r)] = 0
            
    if not final_cent_coord_lst:
        if np.sum(mask) > 15:
            labeled, _ = ndi.label(mask)
            return labeled
        else:
            return np.zeros_like(mask)
    
    # Assign labels based on nearest center
    instance_mask = np.zeros_like(mask, dtype=np.uint16)
    coordlst = np.array(np.where(mask > 0)).T
    
    if len(coordlst) == 0:
        return instance_mask
        
    dists = cdist(coordlst, np.array(final_cent_coord_lst))
    # Normalize by radius to prefer closer centers relative to their size
    radii_array = np.array(final_cent_rad_lst)
    if np.any(radii_array > 0):
        ratio_d = dists / radii_array
    else:
        ratio_d = dists
        
    labels = np.argmin(ratio_d, axis=1) + 1
    
    for idx, coord in enumerate(coordlst):
        instance_mask[coord[0], coord[1], coord[2]] = labels[idx]
        
    return instance_mask
