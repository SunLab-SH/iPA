# filepath: d:\Gitspace\ipa_full\iPA\ipa\processing\segmentation\segmentation_sim_wfm\model_sim_isg_instance\isg_watershed_threshold.py
"""
专门针对ISG的Watershed+Threshold分割方法
基于ISG的特征（小目标、圆形、边界清晰）进行优化
"""

import cv2 as cv
import numpy as np
import os
import json
import argparse
from tqdm import tqdm
from skimage import measure
from scipy import ndimage


def isg_optimized_segmentation(image, debug=False):
    """
    针对ISG优化的分割算法
    ISG特征：小目标、类圆形、边界清晰、颜色与背景有差异
    """
    # 转为灰度图
    if len(image.shape) == 3:
        gray = cv.cvtColor(image, cv.COLOR_RGB2GRAY)
    else:
        gray = image.copy()
    
    # 保存原始图像用于调试
    original_gray = gray.copy()
    
    # Step 1: 预处理 - 轻微去噪但保留边缘
    # 使用双边滤波去噪但保留边缘
    gray = cv.bilateralFilter(gray, 9, 75, 75)
    
    # Step 2: 多种阈值方法组合
    binary_masks = []
    
    # 2.1 Otsu阈值
    _, otsu_binary = cv.threshold(gray, 0, 255, cv.THRESH_BINARY + cv.THRESH_OTSU)
    binary_masks.append(otsu_binary)
    
    # 2.2 自适应阈值 - 小窗口适应局部变化
    adaptive_binary = cv.adaptiveThreshold(gray, 255, cv.ADAPTIVE_THRESH_GAUSSIAN_C,
                                         cv.THRESH_BINARY, 7, 2)
    binary_masks.append(adaptive_binary)
    
    # 2.3 基于直方图的阈值
    hist = cv.calcHist([gray], [0], None, [256], [0, 256])
    # 找到直方图的两个峰值之间的谷底作为阈值
    hist_smooth = ndimage.gaussian_filter1d(hist.flatten(), sigma=2)
    peaks = []
    for i in range(1, len(hist_smooth)-1):
        if hist_smooth[i] > hist_smooth[i-1] and hist_smooth[i] > hist_smooth[i+1]:
            peaks.append((i, hist_smooth[i]))
    
    if len(peaks) >= 2:
        # 找到两个最大的峰值
        peaks.sort(key=lambda x: x[1], reverse=True)
        peak1, peak2 = peaks[0][0], peaks[1][0]
        threshold_val = (peak1 + peak2) // 2
        _, hist_binary = cv.threshold(gray, threshold_val, 255, cv.THRESH_BINARY)
        binary_masks.append(hist_binary)
    
    # Step 3: 投票机制组合多个阈值结果
    vote_sum = np.sum(binary_masks, axis=0)
    # 至少有一半方法认为是前景
    combined_binary = (vote_sum >= len(binary_masks) * 255 * 0.5).astype(np.uint8) * 255
    
    # Step 4: 形态学操作 - 针对ISG的小目标特点
    # 小的开运算去除噪声点
    kernel_small = cv.getStructuringElement(cv.MORPH_ELLIPSE, (2, 2))
    combined_binary = cv.morphologyEx(combined_binary, cv.MORPH_OPEN, kernel_small)
    
    # 小的闭运算填充小孔
    kernel_close = cv.getStructuringElement(cv.MORPH_ELLIPSE, (3, 3))
    combined_binary = cv.morphologyEx(combined_binary, cv.MORPH_CLOSE, kernel_close)
    
    # Step 5: Watershed分离粘连目标
    if np.sum(combined_binary) > 0:
        # 距离变换
        dist_transform = cv.distanceTransform(combined_binary, cv.DIST_L2, 5)
        
        # 找局部最大值作为种子 - 调整参数适应ISG大小
        # ISG通常比较小，所以用较小的min_distance
        from scipy.ndimage import maximum_filter
        local_maxima = maximum_filter(dist_transform, size=7) == dist_transform
        local_maxima = local_maxima & (dist_transform > 0.3 * dist_transform.max())
        
        # 标记种子点
        markers = ndimage.label(local_maxima)[0]
        
        if np.max(markers) > 0:
            # 应用watershed
            img_for_watershed = cv.cvtColor(gray, cv.COLOR_GRAY2BGR)
            watersheds = cv.watershed(img_for_watershed, markers)
            
            # 创建最终mask
            final_mask = np.zeros_like(combined_binary)
            for label in np.unique(watersheds):
                if label <= 0:  # 跳过背景和边界
                    continue
                region_mask = (watersheds == label)
                # 只保留在原始二值图中的区域
                region_mask = region_mask & (combined_binary > 0)
                if np.sum(region_mask) > 0:
                    final_mask[region_mask] = 255
        else:
            final_mask = combined_binary
    else:
        final_mask = combined_binary
    
    # Step 6: 基于ISG特征的后处理过滤
    num_labels, labeled = cv.connectedComponents(final_mask)
    filtered_mask = np.zeros_like(final_mask)
    
    for label in range(1, num_labels + 1):
        region_mask = (labeled == label)
        region_props = measure.regionprops((region_mask).astype(int))[0]
        
        area = region_props.area
        
        # ISG面积范围 - 根据实际ISG大小调整
        if area < 20 or area > 2000:
            continue
        
        # 长宽比 - ISG通常接近圆形
        bbox_height = region_props.bbox[2] - region_props.bbox[0]
        bbox_width = region_props.bbox[3] - region_props.bbox[1]
        aspect_ratio = max(bbox_height, bbox_width) / max(min(bbox_height, bbox_width), 1)
        if aspect_ratio > 2.5:  # 不要太长条
            continue
        
        # 凸度 - ISG边界相对规整
        if region_props.solidity < 0.4:  # 降低要求，因为ISG可能有凹陷
            continue
        
        # 圆形度
        if region_props.perimeter > 0:
            circularity = 4 * np.pi * area / (region_props.perimeter ** 2)
            if circularity < 0.3:  # 不要太不规则
                continue
        
        # 通过所有检查，保留该区域
        filtered_mask[region_mask] = 255
    
    if debug:
        return filtered_mask, {
            'original': original_gray,
            'preprocessed': gray,
            'otsu': binary_masks[0] if len(binary_masks) > 0 else None,
            'adaptive': binary_masks[1] if len(binary_masks) > 1 else None,
            'combined': combined_binary,
            'final': filtered_mask
        }
    
    return filtered_mask


def extract_isg_instances(mask):
    """从分割mask中提取ISG实例"""
    instances = []
    num_labels, labeled = cv.connectedComponents(mask)
    
    for label in range(1, num_labels + 1):
        region_mask = (labeled == label).astype(np.uint8)
        region_props = measure.regionprops(region_mask)[0]
        
        # 边界框
        y1, x1, y2, x2 = region_props.bbox
        bbox = [x1, y1, x2 - x1, y2 - y1]
        
        # 置信度基于形状特征
        area = region_props.area
        solidity = region_props.solidity
        
        # 圆形度
        if region_props.perimeter > 0:
            circularity = 4 * np.pi * area / (region_props.perimeter ** 2)
        else:
            circularity = 0
        
        # 长宽比
        bbox_height = region_props.bbox[2] - region_props.bbox[0]
        bbox_width = region_props.bbox[3] - region_props.bbox[1]
        aspect_ratio = max(bbox_height, bbox_width) / max(min(bbox_height, bbox_width), 1)
        aspect_score = max(0, 1.0 - abs(aspect_ratio - 1.0) / 2.0)
        
        # 综合评分
        score = (solidity * 0.3 + circularity * 0.4 + aspect_score * 0.3)
        score = max(0.1, min(0.95, score))
        
        instances.append({
            'bbox': bbox,
            'mask': region_mask,
            'area': area,
            'score': score,
            'label': 1
        })
    
    return instances


class ISGWatershedEvaluator:
    """ISG Watershed+Threshold评估器"""
    
    def __init__(self, data_root, coco_file):
        self.data_root = data_root
        self.coco_file = coco_file
        
        # 加载COCO数据
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
        """加载图像"""
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
                        return cv.cvtColor(image, cv.COLOR_BGR2RGB), img_info
        
        raise FileNotFoundError(f"Cannot find image: {filename}")
    
    def get_ground_truth_count(self, img_id):
        """获取真实ISG数量"""
        anns = self.annotations.get(img_id, [])
        return len(anns)
    
    def evaluate(self, data_ratio=0.01, save_debug_images=False):
        """评估方法"""
        print("=== ISG Watershed+Threshold Evaluation ===")
        print(f"Total images: {len(self.image_ids)}")
        
        # 随机选择数据子集
        import random
        random.seed(42)
        use_size = int(len(self.image_ids) * data_ratio)
        indices = random.sample(range(len(self.image_ids)), use_size)
        print(f"Using {use_size} images")
        
        predictions = []
        ground_truths = []
        
        debug_dir = "debug_images" if save_debug_images else None
        if debug_dir:
            os.makedirs(debug_dir, exist_ok=True)
        
        for i, idx in enumerate(tqdm(indices, desc="Processing")):
            img_id = self.image_ids[idx]
            
            try:
                # 加载图像
                image, img_info = self.load_image(img_id)
                
                # 获取真实ISG数量
                gt_count = self.get_ground_truth_count(img_id)
                
                # 执行分割
                if debug_dir and i < 5:  # 只保存前5张图的调试信息
                    mask, debug_info = isg_optimized_segmentation(image, debug=True)
                    # 保存调试图像
                    self.save_debug_images(debug_info, img_id, debug_dir)
                else:
                    mask = isg_optimized_segmentation(image)
                
                # 提取实例
                pred_instances = extract_isg_instances(mask)
                
                predictions.append(pred_instances)
                ground_truths.append(gt_count)
                
                # 实时显示进度
                if i < 10:  # 前10张图显示详细信息
                    print(f"Image {img_id}: GT={gt_count}, Pred={len(pred_instances)}")
                
            except Exception as e:
                print(f"Error processing {img_id}: {e}")
                predictions.append([])
                ground_truths.append(0)
        
        # 计算指标
        total_pred = sum(len(pred) for pred in predictions)
        total_true = sum(ground_truths)
        correct = sum(min(len(pred), gt) for pred, gt in zip(predictions, ground_truths))
        
        precision = correct / max(total_pred, 1)
        recall = correct / max(total_true, 1)
        f1 = 2 * precision * recall / max(precision + recall, 1e-6)
        detection_rate = total_pred / max(total_true, 1)
        
        print(f"\nResults:")
        print(f"  Total GT ISGs: {total_true}")
        print(f"  Total Predicted: {total_pred}")
        print(f"  Correct matches: {correct}")
        print(f"  Precision: {precision:.4f}")
        print(f"  Recall: {recall:.4f}")
        print(f"  F1: {f1:.4f}")
        print(f"  Detection Rate: {detection_rate:.4f}")
        
        return {
            'precision': precision,
            'recall': recall,
            'f1': f1,
            'detection_rate': detection_rate,
            'total_true': total_true,
            'total_pred': total_pred
        }
    
    def save_debug_images(self, debug_info, img_id, debug_dir):
        """保存调试图像"""
        import matplotlib.pyplot as plt
        
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        
        axes[0, 0].imshow(debug_info['original'], cmap='gray')
        axes[0, 0].set_title('Original')
        axes[0, 0].axis('off')
        
        axes[0, 1].imshow(debug_info['preprocessed'], cmap='gray')
        axes[0, 1].set_title('Preprocessed')
        axes[0, 1].axis('off')
        
        if debug_info['otsu'] is not None:
            axes[0, 2].imshow(debug_info['otsu'], cmap='gray')
            axes[0, 2].set_title('Otsu')
            axes[0, 2].axis('off')
        
        if debug_info['adaptive'] is not None:
            axes[1, 0].imshow(debug_info['adaptive'], cmap='gray')
            axes[1, 0].set_title('Adaptive')
            axes[1, 0].axis('off')
        
        axes[1, 1].imshow(debug_info['combined'], cmap='gray')
        axes[1, 1].set_title('Combined')
        axes[1, 1].axis('off')
        
        axes[1, 2].imshow(debug_info['final'], cmap='gray')
        axes[1, 2].set_title('Final Result')
        axes[1, 2].axis('off')
        
        plt.tight_layout()
        plt.savefig(os.path.join(debug_dir, f'debug_{img_id}.png'), dpi=150, bbox_inches='tight')
        plt.close()


def main():
    parser = argparse.ArgumentParser(description='ISG Watershed+Threshold Segmentation')
    parser.add_argument('--data_root', default=r'D:\Gitspace\ipa_full\data\SIM\split_isg_img')
    parser.add_argument('--coco_file', default=r'D:\Gitspace\ipa_full\data\SIM\split_isg_img\split_sim_isg_instances.json')
    parser.add_argument('--data_ratio', default=0.01, type=float)
    parser.add_argument('--debug_images', action='store_true', help='保存调试图像')
    
    args = parser.parse_args()
    
    evaluator = ISGWatershedEvaluator(args.data_root, args.coco_file)
    results = evaluator.evaluate(args.data_ratio, args.debug_images)
    
    print(f"\n=== Final Results ===")
    print(f"ISG Watershed+Threshold Method Performance:")
    for key, value in results.items():
        if isinstance(value, float):
            print(f"  {key}: {value:.4f}")
        else:
            print(f"  {key}: {value}")


if __name__ == '__main__':
    main()