# filepath: d:\Gitspace\ipa_full\iPA\ipa\processing\segmentation\segmentation_sim_wfm\model_sim_isg_instance\simple_threshold_watershed.py
"""
简单的Threshold + Watershed ISG分割
1. 阈值分割筛选前景
2. Watershed分离粘连目标
3. 计算指标
"""

import cv2 as cv
import numpy as np
import os
import json
import argparse
from tqdm import tqdm
from skimage import measure
from scipy import ndimage


def threshold_watershed_segmentation(image):
    """阈值+Watershed分割"""
    # 转灰度图
    if len(image.shape) == 3:
        gray = cv.cvtColor(image, cv.COLOR_RGB2GRAY)
    else:
        gray = image.copy()
    
    # Step 1: 阈值分割
    _, binary = cv.threshold(gray, 0, 255, cv.THRESH_BINARY + cv.THRESH_OTSU)
    
    # 形态学开运算去噪
    kernel = cv.getStructuringElement(cv.MORPH_ELLIPSE, (3, 3))
    binary = cv.morphologyEx(binary, cv.MORPH_OPEN, kernel)
    
    # Step 2: Watershed分离
    if np.sum(binary) == 0:
        return np.zeros_like(gray, dtype=np.int32)
    
    # 距离变换
    dist = cv.distanceTransform(binary, cv.DIST_L2, 5)
    
    # 找局部最大值作为种子
    from scipy.ndimage import maximum_filter
    local_maxima = maximum_filter(dist, size=10) == dist
    local_maxima = local_maxima & (dist > 0.3 * dist.max())
    
    # 标记种子点
    markers, _ = ndimage.label(local_maxima)
    
    if np.max(markers) == 0:
        # 如果没有种子点，直接返回连通区域
        num_labels, labels = cv.connectedComponents(binary)
        return labels
    
    # 应用watershed
    img_3d = cv.cvtColor(gray, cv.COLOR_GRAY2BGR)
    watersheds = cv.watershed(img_3d, markers)
    
    # 只保留在原始二值图中的区域
    final_labels = np.zeros_like(watersheds)
    label_id = 1
    for label in np.unique(watersheds):
        if label <= 0:
            continue
        region_mask = (watersheds == label) & (binary > 0)
        if np.sum(region_mask) > 0:
            final_labels[region_mask] = label_id
            label_id += 1
    
    return final_labels


def filter_by_area(labeled_image, min_area=20, max_area=2000):
    """根据面积过滤"""
    filtered = np.zeros_like(labeled_image)
    regions = measure.regionprops(labeled_image)
    
    label_id = 1
    for region in regions:
        if min_area <= region.area <= max_area:
            mask = labeled_image == region.label
            filtered[mask] = label_id
            label_id += 1
    
    return filtered


def extract_instances(labeled_mask):
    """提取实例信息"""
    instances = []
    regions = measure.regionprops(labeled_mask)
    
    for region in regions:
        if region.label == 0:
            continue
        
        # 边界框
        y1, x1, y2, x2 = region.bbox
        bbox = [x1, y1, x2 - x1, y2 - y1]
        
        # 掩码
        mask = (labeled_mask == region.label).astype(np.uint8)
        
        # 简单置信度
        score = min(0.9, max(0.1, region.solidity))
        
        instances.append({
            'bbox': bbox,
            'mask': mask,
            'area': region.area,
            'score': score,
            'label': 1
        })
    
    return instances


class SimpleEvaluator:
    def __init__(self, data_root, coco_file):
        self.data_root = data_root
        self.coco_file = coco_file
        
        with open(coco_file, 'r', encoding='utf-8') as f:
            self.coco_data = json.load(f)
        
        self.images = {img['id']: img for img in self.coco_data['images']}
        self.annotations = {}
        for ann in self.coco_data['annotations']:
            img_id = ann['image_id']
            if img_id not in self.annotations:
                self.annotations[img_id] = []
            self.annotations[img_id].append(ann)
        
        self.image_ids = list(self.images.keys())
    
    def load_image(self, img_id):
        img_info = self.images[img_id]
        filename = img_info['file_name']
        
        images_base = os.path.join(self.data_root, 'images')
        for subdir in os.listdir(images_base):
            subdir_path = os.path.join(images_base, subdir)
            if os.path.isdir(subdir_path):
                img_path = os.path.join(subdir_path, filename)
                if os.path.exists(img_path):
                    image = cv.imread(img_path)
                    if image is not None:
                        return cv.cvtColor(image, cv.COLOR_BGR2RGB)
        raise FileNotFoundError(f"Cannot find image: {filename}")
    
    def get_gt_count(self, img_id):
        anns = self.annotations.get(img_id, [])
        return len(anns)
    
    def evaluate(self, data_ratio=0.01):
        print("=== Threshold + Watershed ISG Segmentation ===")
        print(f"Total images: {len(self.image_ids)}")
        
        # 随机选择数据
        import random
        random.seed(42)
        use_size = int(len(self.image_ids) * data_ratio)
        indices = random.sample(range(len(self.image_ids)), use_size)
        print(f"Using {use_size} images")
        
        total_pred = 0
        total_true = 0
        correct = 0
        
        for idx in tqdm(indices, desc="Processing"):
            img_id = self.image_ids[idx]
            
            try:
                # 加载图像
                image = self.load_image(img_id)
                
                # 真实ISG数量
                gt_count = self.get_gt_count(img_id)
                total_true += gt_count
                
                # 分割
                labeled = threshold_watershed_segmentation(image)
                
                # 面积过滤
                filtered = filter_by_area(labeled)
                
                # 提取实例
                instances = extract_instances(filtered)
                pred_count = len(instances)
                total_pred += pred_count
                
                # 简单匹配
                correct += min(pred_count, gt_count)
                
            except Exception as e:
                print(f"Error processing {img_id}: {e}")
                continue
        
        # 计算指标
        precision = correct / max(total_pred, 1)
        recall = correct / max(total_true, 1)
        f1 = 2 * precision * recall / max(precision + recall, 1e-6)
        
        print(f"\n=== Results ===")
        print(f"Total GT ISGs: {total_true}")
        print(f"Total Predicted: {total_pred}")
        print(f"Correct: {correct}")
        print(f"Precision: {precision:.4f}")
        print(f"Recall: {recall:.4f}")
        print(f"F1 Score: {f1:.4f}")
        
        return {
            'precision': precision,
            'recall': recall,
            'f1': f1,
            'total_true': total_true,
            'total_pred': total_pred,
            'correct': correct
        }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_root', default=r'D:\Gitspace\ipa_full\data\SIM\split_isg_img')
    parser.add_argument('--coco_file', default=r'D:\Gitspace\ipa_full\data\SIM\split_isg_img\split_sim_isg_instances.json')
    parser.add_argument('--data_ratio', default=0.8, type=float)
    
    args = parser.parse_args()
    
    evaluator = SimpleEvaluator(args.data_root, args.coco_file)
    results = evaluator.evaluate(args.data_ratio)


if __name__ == '__main__':
    main()