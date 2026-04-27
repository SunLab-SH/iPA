import numpy as np


def anchors_cython(height, width, stride, base_anchors):
    """
    Generate anchors for RPN.
    Pure Python implementation as fallback for Cython version.
    
    Args:
        height: feature map height
        width: feature map width  
        stride: stride of the feature map
        base_anchors: base anchor boxes
        
    Returns:
        anchors: generated anchors
    """
    return generate_anchors_py(height, width, stride, base_anchors)


def generate_anchors_py(height, width, stride, base_anchors):
    """
    Generate anchors for RPN using pure Python.
    
    Args:
        height: feature map height
        width: feature map width
        stride: stride of the feature map
        base_anchors: base anchor boxes [N, 4]
        
    Returns:
        anchors: generated anchors [height*width*N, 4]
    """
    base_anchors = np.array(base_anchors, dtype=np.float32)
    num_anchors = base_anchors.shape[0]
    
    # Generate shifts for all positions
    shift_x = np.arange(0, width) * stride
    shift_y = np.arange(0, height) * stride
    shift_x, shift_y = np.meshgrid(shift_x, shift_y)
    
    shifts = np.vstack((shift_x.ravel(), shift_y.ravel(),
                       shift_x.ravel(), shift_y.ravel())).transpose()
    
    # Add A anchors (1, A, 4) to
    # cell K shifts (K, 1, 4) to get
    # shift anchors (K, A, 4)
    # reshape to (K*A, 4) shifted anchors
    A = num_anchors
    K = shifts.shape[0]
    anchors = (base_anchors.reshape((1, A, 4)) +
               shifts.reshape((1, K, 4)).transpose((1, 0, 2)))
    anchors = anchors.reshape((K * A, 4))
    
    return anchors


def generate_base_anchors(base_size=16, ratios=[0.5, 1, 2], scales=[8, 16, 32]):
    """
    Generate base anchor boxes.
    
    Args:
        base_size: base size of anchors
        ratios: aspect ratios
        scales: scales
        
    Returns:
        anchors: base anchor boxes [len(ratios)*len(scales), 4]
    """
    base_anchor = np.array([1, 1, base_size, base_size]) - 1
    ratio_anchors = _ratio_enum(base_anchor, ratios)
    anchors = np.vstack([_scale_enum(ratio_anchors[i, :], scales)
                        for i in range(ratio_anchors.shape[0])])
    return anchors


def _whctrs(anchor):
    """
    Return width, height, x center, and y center for an anchor (window).
    """
    w = anchor[2] - anchor[0] + 1
    h = anchor[3] - anchor[1] + 1
    x_ctr = anchor[0] + 0.5 * (w - 1)
    y_ctr = anchor[1] + 0.5 * (h - 1)
    return w, h, x_ctr, y_ctr


def _mkanchors(ws, hs, x_ctr, y_ctr):
    """
    Given a vector of widths (ws) and heights (hs) around a center
    (x_ctr, y_ctr), output a set of anchors (windows).
    """
    ws = ws[:, np.newaxis]
    hs = hs[:, np.newaxis]
    anchors = np.hstack((x_ctr - 0.5 * (ws - 1),
                        y_ctr - 0.5 * (hs - 1),
                        x_ctr + 0.5 * (ws - 1),
                        y_ctr + 0.5 * (hs - 1)))
    return anchors


def _ratio_enum(anchor, ratios):
    """
    Enumerate a set of anchors for each aspect ratio wrt an anchor.
    """
    w, h, x_ctr, y_ctr = _whctrs(anchor)
    size = w * h
    size_ratios = size / ratios
    ws = np.round(np.sqrt(size_ratios))
    hs = np.round(ws * ratios)
    anchors = _mkanchors(ws, hs, x_ctr, y_ctr)
    return anchors


def _scale_enum(anchor, scales):
    """
    Enumerate a set of anchors for each scale wrt an anchor.
    """
    w, h, x_ctr, y_ctr = _whctrs(anchor)
    ws = w * np.array(scales)
    hs = h * np.array(scales)
    anchors = _mkanchors(ws, hs, x_ctr, y_ctr)
    return anchors


# Vectorized anchor generation for better performance
def generate_anchors_vectorized(height, width, stride, base_anchors):
    """
    Vectorized version of anchor generation for better performance.
    """
    base_anchors = np.array(base_anchors, dtype=np.float32)
    
    # Create grid of shifts
    shifts_x = np.arange(0, width) * stride
    shifts_y = np.arange(0, height) * stride
    shift_x, shift_y = np.meshgrid(shifts_x, shifts_y)
    
    # Flatten and stack
    shifts = np.stack((shift_x.ravel(), shift_y.ravel(), 
                      shift_x.ravel(), shift_y.ravel()), axis=1)
    
    # Broadcasting to get all combinations
    num_base_anchors = base_anchors.shape[0]
    num_shifts = shifts.shape[0]
    
    # Reshape for broadcasting
    shifts = shifts.reshape((num_shifts, 1, 4))
    base_anchors = base_anchors.reshape((1, num_base_anchors, 4))
    
    # Generate all anchors
    all_anchors = shifts + base_anchors
    
    # Reshape to final form
    anchors = all_anchors.reshape((num_shifts * num_base_anchors, 4))
    
    return anchors