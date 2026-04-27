import numpy as np


def bbox_overlaps(boxes, query_boxes):
    """
    Calculate overlaps between boxes and query_boxes using pure Python.
    This is a fallback when Cython version is not available.
    
    :param boxes: n * 4 bounding boxes
    :param query_boxes: k * 4 bounding boxes
    :return: overlaps: n * k overlaps
    """
    return bbox_overlaps_py(boxes, query_boxes)


def bbox_overlaps_py(boxes, query_boxes):
    """
    Determine overlaps between boxes and query_boxes
    :param boxes: n * 4 bounding boxes
    :param query_boxes: k * 4 bounding boxes
    :return: overlaps: n * k overlaps
    """
    n_ = boxes.shape[0]
    k_ = query_boxes.shape[0]
    overlaps = np.zeros((n_, k_), dtype=np.float)
    
    for k in range(k_):
        query_box_area = (query_boxes[k, 2] - query_boxes[k, 0] + 1) * (query_boxes[k, 3] - query_boxes[k, 1] + 1)
        for n in range(n_):
            iw = min(boxes[n, 2], query_boxes[k, 2]) - max(boxes[n, 0], query_boxes[k, 0]) + 1
            if iw > 0:
                ih = min(boxes[n, 3], query_boxes[k, 3]) - max(boxes[n, 1], query_boxes[k, 1]) + 1
                if ih > 0:
                    box_area = (boxes[n, 2] - boxes[n, 0] + 1) * (boxes[n, 3] - boxes[n, 1] + 1)
                    all_area = float(box_area + query_box_area - iw * ih)
                    overlaps[n, k] = iw * ih / all_area
    return overlaps


def bbox_overlaps_vectorized(boxes, query_boxes):
    """
    Vectorized version of bbox overlaps calculation for better performance.
    """
    boxes = np.asarray(boxes, dtype=np.float32)
    query_boxes = np.asarray(query_boxes, dtype=np.float32)
    
    n = boxes.shape[0]
    k = query_boxes.shape[0]
    
    overlaps = np.zeros((n, k), dtype=np.float32)
    
    if n == 0 or k == 0:
        return overlaps
    
    # Compute areas
    box_areas = (boxes[:, 2] - boxes[:, 0] + 1) * (boxes[:, 3] - boxes[:, 1] + 1)
    query_areas = (query_boxes[:, 2] - query_boxes[:, 0] + 1) * (query_boxes[:, 3] - query_boxes[:, 1] + 1)
    
    for k_idx in range(k):
        query_area = query_areas[k_idx]
        
        # Intersection coordinates
        xx1 = np.maximum(boxes[:, 0], query_boxes[k_idx, 0])
        yy1 = np.maximum(boxes[:, 1], query_boxes[k_idx, 1])
        xx2 = np.minimum(boxes[:, 2], query_boxes[k_idx, 2])
        yy2 = np.minimum(boxes[:, 3], query_boxes[k_idx, 3])
        
        # Intersection area
        w = np.maximum(0., xx2 - xx1 + 1)
        h = np.maximum(0., yy2 - yy1 + 1)
        inter = w * h
        
        # Union area
        union = box_areas + query_area - inter
        
        # IoU
        overlaps[:, k_idx] = inter / union
    
    return overlaps