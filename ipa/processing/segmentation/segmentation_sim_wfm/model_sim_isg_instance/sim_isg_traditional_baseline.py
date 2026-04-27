# filepath: d:\Gitspace\ipa_full\iPA\ipa\processing\segmentation\segmentation_sim_wfm\model_sim_isg_instance\sim_isg_traditional_baseline.py
"""
SIM ISG传统图像处理基线方法
使用Watershed、阈值分割、形态学操作等传统方法进行ISG分割
用于与Mask R-CNN进行对比评估
"""

import cv2 as cv
import numpy as np
import os
import json
import argparse
import time
import logging
from pathlib import Path
from PIL import Image
from sklearn.metrics import precision_recall_fscore_support
from scipy import ndimage
from skimage import morphology, segmentation, measure, filters
from skimage.segmentation import watershed
import matplotlib.pyplot as plt
from tqdm import tqdm
import sys

# 导入配置
try:
    from traditional_configs import get_config, list_configs
except ImportError:
    print("Warning: traditional_configs.py not found, using default config")
    def get_config(name='base'):
        return {
            'gaussian_blur_kernel': 3, 'median_filter_size': 3,
            'threshold_method': 'otsu', 'manual_threshold': 128,
            'adaptive_block_size': 11, 'adaptive_c': 2,
            'watershed_connectivity': 2, 'min_distance': 10, 'watershed_threshold': 0.3,
            'opening_kernel_size': 3, 'closing_kernel_size': 5,
            'erosion_iterations': 1, 'dilation_iterations': 2,
            'min_area': 50, 'max_area': 10000, 'min_solidity': 0.3,
            'aspect_ratio_range': (0.3, 3.0), 'circularity_threshold': 0.4, 'use_shape_filter': True,
        }
    def list_configs():
        return ['base']


# 简单的局部最大值检测实现
def find_local_maxima(image, min_distance=1, threshold_abs=None):
    """局部最大值检测"""
    if threshold_abs is None:
        threshold_abs = 0.1 * image.max()
    
    try:
        # 尝试使用scipy的maximum_filter
        from scipy.ndimage import maximum_filter
        local_maxima = maximum_filter(image, size=min_distance*2+1) == image
        local_maxima = local_maxima & (image > threshold_abs)
        coords = np.where(local_maxima)
        return coords
    except ImportError:
        # 手动实现
        peaks = []
        h, w = image.shape
        half_dist = min_distance // 2
        
        for y in range(half_dist, h - half_dist):
            for x in range(half_dist, w - half_dist):
                if image[y, x] > threshold_abs:
                    # 检查是否为局部最大值
                    local_region = image[y-half_dist:y+half_dist+1, 
                                       x-half_dist:x+half_dist+1]
                    if image[y, x] == local_region.max():
                        peaks.append((y, x))
        
        if peaks:
            peaks = list(zip(*peaks))
            return (np.array(peaks[0]), np.array(peaks[1]))
        else:
            return (np.array([]), np.array([]))


class TraditionalISGSegmentation:
    """传统图像处理方法的ISG分割器"""
    
    def __init__(self, method='watershed_combined', config=None):
        """
        初始化分割器
        Args:
            method: 分割方法 ('threshold', 'watershed', 'watershed_combined', 'adaptive_threshold')
            config: 方法配置参数 (dict) 或配置名称 (str)
        """
        self.method = method
        
        # 处理配置参数
        if isinstance(config, str):
            # 如果config是字符串，从配置文件加载
            self.config = get_config(config)
        elif isinstance(config, dict):
            # 如果config是字典，直接使用
            self.config = config
        else:
            # 使用默认配置
            self.config = self._get_default_config()
        
    def _get_default_config(self):
        """获取默认配置参数"""
        return {
            # 预处理参数
            'gaussian_blur_kernel': 3,
            'median_filter_size': 3,
            
            # 阈值分割参数
            'threshold_method': 'otsu',  # 'otsu', 'adaptive', 'manual'
            'manual_threshold': 128,
            'adaptive_block_size': 11,
            'adaptive_c': 2,
            
            # Watershed参数
            'watershed_connectivity': 2,
            'min_distance': 10,  # 局部最大值之间的最小距离
            'watershed_threshold': 0.3,  # 距离变换的阈值比例
            
            # 形态学操作参数
            'opening_kernel_size': 3,
            'closing_kernel_size': 5,
            'erosion_iterations': 1,
            'dilation_iterations': 2,
            
            # 后处理参数
            'min_area': 50,  # 最小连通区域面积
            'max_area': 10000,  # 最大连通区域面积
            'min_solidity': 0.3,  # 最小凸度
            'aspect_ratio_range': (0.3, 3.0),  # 长宽比范围
            
            # ISG特异性参数
            'circularity_threshold': 0.4,  # 圆形度阈值
            'use_shape_filter': True,  # 是否使用形状过滤
        }
    
    def preprocess_image(self, image):
        """图像预处理"""
        if len(image.shape) == 3:
            # 转为灰度图
            gray = cv.cvtColor(image, cv.COLOR_RGB2GRAY)
        else:
            gray = image.copy()
        
        # 高斯模糊去噪
        if self.config['gaussian_blur_kernel'] > 0:
            gray = cv.GaussianBlur(gray, 
                                 (self.config['gaussian_blur_kernel'], self.config['gaussian_blur_kernel']), 
                                 0)
        
        # 中值滤波
        if self.config['median_filter_size'] > 0:
            gray = cv.medianBlur(gray, self.config['median_filter_size'])
        
        return gray
    
    def threshold_segmentation(self, gray_image):
        """阈值分割方法"""
        if self.config['threshold_method'] == 'otsu':
            _, binary = cv.threshold(gray_image, 0, 255, cv.THRESH_BINARY + cv.THRESH_OTSU)
        elif self.config['threshold_method'] == 'adaptive':
            binary = cv.adaptiveThreshold(gray_image, 255, cv.ADAPTIVE_THRESH_GAUSSIAN_C,
                                        cv.THRESH_BINARY, 
                                        self.config['adaptive_block_size'],
                                        self.config['adaptive_c'])
        else:  # manual
            _, binary = cv.threshold(gray_image, self.config['manual_threshold'], 255, cv.THRESH_BINARY)
        
        return binary
    
    def watershed_segmentation(self, gray_image):
        """Watershed分割方法"""
        # 阈值分割得到前景
        binary = self.threshold_segmentation(gray_image)
        
        # 形态学操作
        kernel = cv.getStructuringElement(cv.MORPH_ELLIPSE, 
                                        (self.config['opening_kernel_size'], 
                                         self.config['opening_kernel_size']))
        opening = cv.morphologyEx(binary, cv.MORPH_OPEN, kernel)
        
        # 如果没有前景像素，直接返回
        if np.sum(opening) == 0:
            return np.zeros_like(gray_image, dtype=np.int32), binary
        
        # 距离变换
        dist_transform = cv.distanceTransform(opening, cv.DIST_L2, 5)
        
        # 寻找局部最大值作为种子点
        local_maxima_coords = find_local_maxima(dist_transform, 
                                       min_distance=self.config['min_distance'],
                                       threshold_abs=dist_transform.max() * self.config['watershed_threshold'])
        
        # 创建标记
        markers = np.zeros_like(dist_transform, dtype=np.int32)
        if len(local_maxima_coords[0]) > 0:
            markers[local_maxima_coords[0], local_maxima_coords[1]] = np.arange(1, len(local_maxima_coords[0]) + 1)
        
        # 如果没有找到种子点，使用形态学方法
        if np.max(markers) == 0:
            # 确定背景
            sure_bg = cv.dilate(opening, kernel, iterations=3)
            # 确定前景
            if dist_transform.max() > 0:
                dist_transform_norm = (dist_transform / dist_transform.max() * 255).astype(np.uint8)
                _, sure_fg = cv.threshold(dist_transform_norm, int(0.7 * dist_transform_norm.max()), 255, 0)
                sure_fg = np.uint8(sure_fg)
                
                # 找到连通区域作为标记
                num_labels, labels = cv.connectedComponents(np.array(sure_fg, dtype=np.uint8))
                markers = labels
        
        # 如果仍然没有标记，返回简单的连通区域
        if np.max(markers) == 0:
            num_labels, watersheds = cv.connectedComponents(opening)
            return watersheds, binary
        
        # 应用watershed
        # 需要3通道图像用于watershed
        if len(gray_image.shape) == 2:
            img_for_watershed = cv.cvtColor(gray_image, cv.COLOR_GRAY2BGR)
        else:
            img_for_watershed = gray_image
            
        watersheds = cv.watershed(img_for_watershed, markers)
        
        return watersheds, binary
    
    def combined_watershed_threshold(self, gray_image):
        """组合Watershed和阈值的方法"""
        # 先用阈值分割得到粗略区域
        binary = self.threshold_segmentation(gray_image)
        
        # 形态学操作清理
        kernel_open = cv.getStructuringElement(cv.MORPH_ELLIPSE, 
                                              (self.config['opening_kernel_size'], 
                                               self.config['opening_kernel_size']))
        binary_clean = cv.morphologyEx(binary, cv.MORPH_OPEN, kernel_open)
        
        kernel_close = cv.getStructuringElement(cv.MORPH_ELLIPSE,
                                               (self.config['closing_kernel_size'],
                                                self.config['closing_kernel_size']))
        binary_clean = cv.morphologyEx(binary_clean, cv.MORPH_CLOSE, kernel_close)
        
        # 在清理后的二值图像上应用watershed进一步分离
        if np.sum(binary_clean) > 0:  # 确保有前景像素
            try:
                watersheds, _ = self.watershed_segmentation(gray_image * (binary_clean > 0).astype(np.uint8))
                
                # 将watershed结果与原始阈值结果结合
                final_mask = np.zeros_like(binary_clean)
                for label in np.unique(watersheds):
                    if label <= 0:  # 跳过背景和边界
                        continue
                    region_mask = (watersheds == label)
                    if np.sum(region_mask) >= self.config['min_area']:
                        final_mask[region_mask] = 255
            except Exception as e:
                print(f"Watershed failed, using binary mask: {e}")
                final_mask = binary_clean
        else:
            final_mask = binary_clean
            
        return final_mask
    
    def adaptive_threshold_method(self, gray_image):
        """自适应阈值方法"""
        # 多尺度自适应阈值
        results = []
        block_sizes = [7, 11, 15, 21]
        
        for block_size in block_sizes:
            adaptive_binary = cv.adaptiveThreshold(gray_image, 255, 
                                                 cv.ADAPTIVE_THRESH_GAUSSIAN_C,
                                                 cv.THRESH_BINARY, 
                                                 block_size, 2)
            results.append(adaptive_binary)
        
        # 投票机制
        vote_sum = np.sum(results, axis=0)
        final_binary = (vote_sum >= len(block_sizes) * 255 * 0.5).astype(np.uint8) * 255
        
        return final_binary
    
    def filter_by_shape_properties(self, labeled_image):
        """根据形状属性过滤连通区域"""
        filtered_labels = np.zeros_like(labeled_image)
        regions = measure.regionprops(labeled_image)
        
        valid_label = 1
        for region in regions:
            area = region.area
            
            # 面积过滤
            if area < self.config['min_area'] or area > self.config['max_area']:
                continue
            
            # 形状特征过滤
            if self.config['use_shape_filter']:
                # 长宽比
                bbox_height = region.bbox[2] - region.bbox[0]
                bbox_width = region.bbox[3] - region.bbox[1]
                aspect_ratio = max(bbox_height, bbox_width) / max(min(bbox_height, bbox_width), 1)
                
                if not (self.config['aspect_ratio_range'][0] <= aspect_ratio <= self.config['aspect_ratio_range'][1]):
                    continue
                
                # 凸度
                if region.solidity < self.config['min_solidity']:
                    continue
                
                # 圆形度 (4π*面积/周长²)
                if region.perimeter > 0:
                    circularity = 4 * np.pi * area / (region.perimeter ** 2)
                    if circularity < self.config['circularity_threshold']:
                        continue
            
            # 保留满足条件的区域
            mask = labeled_image == region.label
            filtered_labels[mask] = valid_label
            valid_label += 1
        
        return filtered_labels
    
    def segment(self, image):
        """执行分割"""
        # 预处理
        gray = self.preprocess_image(image)
        
        # 根据方法选择分割算法
        try:
            if self.method == 'threshold':
                mask = self.threshold_segmentation(gray)
            elif self.method == 'watershed':
                watersheds, binary = self.watershed_segmentation(gray)
                mask = (watersheds > 0).astype(np.uint8) * 255
            elif self.method == 'watershed_combined':
                mask = self.combined_watershed_threshold(gray)
            elif self.method == 'adaptive_threshold':
                mask = self.adaptive_threshold_method(gray)
            else:
                raise ValueError(f"Unknown method: {self.method}")
            
            # 连通区域标记
            if mask is None or np.sum(mask) == 0:
                return np.zeros_like(gray, dtype=np.int32)
            
            num_labels, labeled_mask = cv.connectedComponents(mask.astype(np.uint8))
            
            # 形状属性过滤
            filtered_labeled = self.filter_by_shape_properties(labeled_mask)
            
            return filtered_labeled
            
        except Exception as e:
            print(f"Segmentation failed: {e}")
            return np.zeros_like(gray, dtype=np.int32)
    
    def extract_instances(self, labeled_mask):
        """从标记图像中提取实例信息"""
        instances = []
        regions = measure.regionprops(labeled_mask)
        
        for region in regions:
            if region.label == 0:  # 跳过背景
                continue
            
            # 边界框 [x, y, width, height]
            y1, x1, y2, x2 = region.bbox
            bbox = [x1, y1, x2 - x1, y2 - y1]
            
            # 掩码
            mask = (labeled_mask == region.label).astype(np.uint8)
            
            # 置信度（基于形状特征计算）
            score = self._calculate_confidence_score(region)
            
            instances.append({
                'bbox': bbox,
                'mask': mask,
                'area': region.area,
                'score': score,
                'label': 1  # ISG类别
            })
        
        return instances
    
    def _calculate_confidence_score(self, region):
        """基于形状特征计算置信度分数"""
        # 基于多个形状特征的综合评分
        
        # 凸度评分 (0-1)
        solidity_score = min(region.solidity / 0.8, 1.0)
        
        # 圆形度评分 (0-1)
        if region.perimeter > 0:
            circularity = 4 * np.pi * region.area / (region.perimeter ** 2)
            circularity_score = min(circularity / 0.8, 1.0)
        else:
            circularity_score = 0.0
        
        # 长宽比评分 (0-1)
        bbox_height = region.bbox[2] - region.bbox[0]
        bbox_width = region.bbox[3] - region.bbox[1]
        aspect_ratio = max(bbox_height, bbox_width) / max(min(bbox_height, bbox_width), 1)
        aspect_score = max(0, 1.0 - abs(aspect_ratio - 1.0) / 2.0)  # 接近1.0的比例更好
        
        # 面积评分 (基于预期的ISG面积范围)
        ideal_area = (self.config['min_area'] + min(self.config['max_area'], 1000)) / 2
        area_score = 1.0 - abs(region.area - ideal_area) / ideal_area
        area_score = max(0, min(area_score, 1.0))
        
        # 综合评分
        total_score = (solidity_score * 0.3 + 
                      circularity_score * 0.3 + 
                      aspect_score * 0.2 + 
                      area_score * 0.2)
        
        return max(0.1, min(0.95, total_score))  # 限制在合理范围内


class TraditionalISGEvaluator:
    """传统方法评估器 - 与Mask R-CNN保持一致的评估指标"""
    
    def __init__(self, data_root, coco_file, method='watershed_combined', config=None):
        self.data_root = data_root
        self.coco_file = coco_file
        self.segmenter = TraditionalISGSegmentation(method, config)
        
        # 加载COCO数据
        with open(coco_file, 'r', encoding='utf-8') as f:
            self.coco_data = json.load(f)
        
        # 构建数据映射
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
        
        # 查找图像文件（与Mask R-CNN脚本相同的逻辑）
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
    
    def get_ground_truth(self, img_id):
        """获取真实标注"""
        anns = self.annotations.get(img_id, [])
        gt_masks = []
        gt_bboxes = []
        
        for ann in anns:
            # 边界框
            bbox = ann['bbox']
            gt_bboxes.append(bbox)
            
            # 掩码
            if 'segmentation' in ann:
                if isinstance(ann['segmentation'], dict):  # RLE格式
                    import pycocotools.mask as mask_util
                    rle = ann['segmentation']
                    if isinstance(rle.get('counts'), list):
                        compressed_rle = mask_util.frPyObjects(rle, rle['size'][0], rle['size'][1])
                        mask = mask_util.decode(compressed_rle)
                    else:
                        mask = mask_util.decode(rle)
                    gt_masks.append(mask)
        
        return gt_masks, gt_bboxes
    
    def compute_metrics(self, predictions, ground_truths, threshold=0.5):
        """计算评估指标 - 与Mask R-CNN保持一致"""
        total_pred_objects = 0
        total_true_objects = 0
        correct_detections = 0
        all_scores = []
        
        for pred_instances, (gt_masks, gt_bboxes) in zip(predictions, ground_truths):
            # 统计真实ISG数量
            n_true = len(gt_masks)
            total_true_objects += n_true
            
            # 统计预测ISG数量
            if len(pred_instances) > 0:
                pred_scores = [inst['score'] for inst in pred_instances]
                high_conf_mask = np.array(pred_scores) >= threshold
                n_pred = np.sum(high_conf_mask)
                total_pred_objects += n_pred
                
                # 简化的正确检测计算：min(预测数量, 真实数量)
                correct_detections += min(n_pred, n_true)
                
                # 收集所有分数
                all_scores.extend(pred_scores)
        
        # 计算基本指标
        if total_pred_objects > 0:
            precision = correct_detections / total_pred_objects
        else:
            precision = 0.0
        
        if total_true_objects > 0:
            recall = correct_detections / total_true_objects
            detection_rate = total_pred_objects / total_true_objects
        else:
            recall = 0.0
            detection_rate = 0.0
        
        if precision + recall > 0:
            f1 = 2 * (precision * recall) / (precision + recall)
        else:
            f1 = 0.0
        
        accuracy = recall  # 简化为召回率
        
        # AUC和AP
        auc = ap = 0.0
        if len(all_scores) > 0:
            ap = np.mean(all_scores)  # 简化为平均置信度
        
        return {
            'accuracy': float(accuracy),
            'precision': float(precision),
            'recall': float(recall),
            'f1': float(f1),
            'auc': float(auc),
            'ap': float(ap),
            'detection_rate': float(detection_rate)
        }
    
    def evaluate_subset(self, image_indices, show_progress=True):
        """评估图像子集"""
        predictions = []
        ground_truths = []
        
        iterator = tqdm(image_indices, desc="Processing images") if show_progress else image_indices
        
        for idx in iterator:
            img_id = self.image_ids[idx]
            
            try:
                # 加载图像
                image, img_info = self.load_image(img_id)
                
                # 获取真实标注
                gt_masks, gt_bboxes = self.get_ground_truth(img_id)
                
                # 执行分割
                labeled_mask = self.segmenter.segment(image)
                
                # 提取实例
                pred_instances = self.segmenter.extract_instances(labeled_mask)
                
                predictions.append(pred_instances)
                ground_truths.append((gt_masks, gt_bboxes))
                
            except Exception as e:
                print(f"Error processing image {img_id}: {e}")
                # 添加空预测
                predictions.append([])
                ground_truths.append(([], []))
                continue
        
        return predictions, ground_truths
    
    def evaluate_full(self, data_ratio=0.01, train_val_split=0.8, val_test_split=0.1):
        """完整评估 - 与Mask R-CNN相同的数据分割"""
        print(f"\n=== Traditional ISG Segmentation Evaluation ===")
        print(f"Method: {self.segmenter.method}")
        print(f"Total images: {len(self.image_ids)}")
        
        # 按比例选择数据
        full_size = len(self.image_ids)
        use_size = int(full_size * data_ratio)
        
        if use_size < full_size:
            import random
            random.seed(42)  # 固定随机种子确保可复现
            indices = random.sample(range(full_size), use_size)
        else:
            indices = list(range(full_size))
        
        print(f"Using {len(indices)}/{full_size} images ({data_ratio:.1%})")
        
        # 数据分割
        dataset_size = len(indices)
        train_size = int(train_val_split * dataset_size)
        val_size = int(val_test_split * dataset_size)
        test_size = dataset_size - train_size - val_size
        
        train_indices = indices[:train_size]
        val_indices = indices[train_size:train_size + val_size]
        test_indices = indices[train_size + val_size:]
        
        print(f"Train: {len(train_indices)}, Val: {len(val_indices)}, Test: {len(test_indices)}")
        
        # 评估验证集
        print("\n--- Validation Set Evaluation ---")
        val_predictions, val_ground_truths = self.evaluate_subset(val_indices)
        val_metrics = self.compute_metrics(val_predictions, val_ground_truths)
        
        print("Validation Metrics:")
        for key, value in val_metrics.items():
            print(f"  Val_{key}: {value:.4f}")
        
        # 评估测试集
        print("\n--- Test Set Evaluation ---")
        test_predictions, test_ground_truths = self.evaluate_subset(test_indices)
        test_metrics = self.compute_metrics(test_predictions, test_ground_truths)
        
        print("Test Metrics:")
        for key, value in test_metrics.items():
            print(f"  Test_{key}: {value:.4f}")
        
        return {
            'validation': val_metrics,
            'test': test_metrics,
            'method': self.segmenter.method,
            'config': self.segmenter.config
        }


def setup_logger(output_dir, name='traditional_isg'):
    """设置日志记录器"""
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    
    # 创建文件处理器
    log_file = os.path.join(output_dir, f'{name}.log')
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    
    # 创建控制台处理器
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    
    # 创建格式器
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)
    
    # 添加处理器
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return logger


def main():
    parser = argparse.ArgumentParser(description='Traditional ISG Segmentation Baseline')
    parser.add_argument('--data_root', default=r'D:\Gitspace\ipa_full\data\SIM\split_isg_img',
                       help='SIM ISG数据根目录')
    parser.add_argument('--coco_file', default=r'D:\Gitspace\ipa_full\data\SIM\split_isg_img\split_sim_isg_instances.json',
                       help='COCO格式标注文件')
    parser.add_argument('--output_dir', default=r'D:\Gitspace\ipa_full\data\SIM\split_isg_img\traditional_results',
                       help='结果输出目录')
    parser.add_argument('--method', default='watershed_combined', 
                       choices=['threshold', 'watershed', 'watershed_combined', 'adaptive_threshold'],
                       help='分割方法')
    parser.add_argument('--config', default='combined', type=str,
                       help='配置名称或配置文件路径')
    parser.add_argument('--data_ratio', default=0.8, type=float,
                       help='使用数据的比例 (0.0-1.0)')
    parser.add_argument('--train_val_split', default=0.8, type=float,
                       help='训练数据比例')
    parser.add_argument('--val_test_split', default=0.1, type=float,
                       help='验证数据比例')
    parser.add_argument('--compare_methods', action='store_true',
                       help='比较所有传统方法')
    parser.add_argument('--compare_configs', action='store_true',
                       help='比较不同配置')
    parser.add_argument('--save_visualizations', action='store_true',
                       help='保存可视化结果')
    parser.add_argument('--list_configs', action='store_true',
                       help='列出所有可用配置')
    
    args = parser.parse_args()
    
    if args.list_configs:
        list_configs()
        return
    
    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 设置日志
    logger = setup_logger(args.output_dir, 'traditional_isg')
    logger.info(f"Traditional ISG评估参数: {args}")
    
    start_time = time.time()
    
    if args.compare_configs:
        # 比较不同配置（固定方法）
        configs = ['base', 'small_isg', 'large_isg', 'high_precision', 'high_recall', 'combined']
        results = {}
        
        for config_name in configs:
            print(f"\n{'='*50}")
            print(f"Evaluating config: {config_name} with method: {args.method}")
            print(f"{'='*50}")
            
            evaluator = TraditionalISGEvaluator(args.data_root, args.coco_file, args.method, config_name)
            config_results = evaluator.evaluate_full(args.data_ratio, args.train_val_split, args.val_test_split)
            results[f"{args.method}_{config_name}"] = config_results
            
            logger.info(f"Config {config_name} results: {config_results}")
        
        # 保存比较结果
        comparison_file = os.path.join(args.output_dir, f'{args.method}_config_comparison.json')
        with open(comparison_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        
        # 打印比较总结
        print(f"\n{'='*70}")
        print(f"{args.method.upper()} CONFIG COMPARISON SUMMARY")
        print(f"{'='*70}")
        print(f"{'Config':<20} {'Test_F1':<10} {'Test_Precision':<15} {'Test_Recall':<12} {'Test_AP':<10}")
        print("-" * 70)
        
        for config_key, result in results.items():
            config_name = config_key.replace(f"{args.method}_", "")
            test_metrics = result['test']
            print(f"{config_name:<20} {test_metrics['f1']:<10.4f} {test_metrics['precision']:<15.4f} "
                  f"{test_metrics['recall']:<12.4f} {test_metrics['ap']:<10.4f}")
    
    elif args.compare_methods:
        # 比较所有方法（使用指定配置）
        methods = ['threshold', 'watershed', 'watershed_combined', 'adaptive_threshold']
        results = {}
        
        for method in methods:
            print(f"\n{'='*50}")
            print(f"Evaluating method: {method} with config: {args.config}")
            print(f"{'='*50}")
            
            evaluator = TraditionalISGEvaluator(args.data_root, args.coco_file, method, args.config)
            method_results = evaluator.evaluate_full(args.data_ratio, args.train_val_split, args.val_test_split)
            results[method] = method_results
            
            logger.info(f"Method {method} results: {method_results}")
        
        # 保存比较结果
        comparison_file = os.path.join(args.output_dir, f'{args.config}_method_comparison.json')
        with open(comparison_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        
        # 打印比较总结
        print(f"\n{'='*60}")
        print(f"METHOD COMPARISON SUMMARY (Config: {args.config})")
        print(f"{'='*60}")
        print(f"{'Method':<20} {'Test_F1':<10} {'Test_Precision':<15} {'Test_Recall':<12} {'Test_AP':<10}")
        print("-" * 70)
        
        for method, result in results.items():
            test_metrics = result['test']
            print(f"{method:<20} {test_metrics['f1']:<10.4f} {test_metrics['precision']:<15.4f} "
                  f"{test_metrics['recall']:<12.4f} {test_metrics['ap']:<10.4f}")
        
    else:
        # 单个方法和配置评估
        print(f"\nEvaluating: {args.method} with config: {args.config}")
        evaluator = TraditionalISGEvaluator(args.data_root, args.coco_file, args.method, args.config)
        results = evaluator.evaluate_full(args.data_ratio, args.train_val_split, args.val_test_split)
        
        # 保存结果
        results_file = os.path.join(args.output_dir, f'{args.method}_{args.config}_results.json')
        with open(results_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        
        logger.info(f"Results saved to: {results_file}")
    
    total_time = time.time() - start_time
    print(f"\nTotal evaluation time: {total_time:.2f} seconds")
    logger.info(f"Evaluation completed in {total_time:.2f} seconds")


if __name__ == '__main__':
    main()