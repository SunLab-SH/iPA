import os
import cv2 as cv
import numpy as np
import json
import argparse as arg
from collections import defaultdict
import pycocotools.mask as mask_util
from pathlib import Path


def load_split_mapping(mapping_path):
    """加载切分映射信息"""
    with open(mapping_path, 'r', encoding='utf-8') as f:
        mapping = json.load(f)
    return mapping


def load_predictions(predictions_path):
    """加载预测结果（COCO格式）"""
    with open(predictions_path, 'r', encoding='utf-8') as f:
        predictions = json.load(f)
    return predictions


def group_predictions_by_original_image(predictions, split_mapping):
    """按原始图像分组预测结果"""
    image_groups = defaultdict(list)
    
    for pred in predictions:
        image_id = pred['image_id']
        if str(image_id) in split_mapping:
            mapping_info = split_mapping[str(image_id)]
            original_image = mapping_info['original_image']
            
            # 添加映射信息到预测结果
            pred['split_info'] = mapping_info
            image_groups[original_image].append(pred)
    
    return image_groups


def merge_predictions_for_image(predictions_group, original_image_shape):
    """合并单个原始图像的所有预测结果"""
    original_h, original_w = original_image_shape[:2]
    merged_mask = np.zeros((original_h, original_w), dtype=np.uint16)
    
    instance_id = 1
    merged_annotations = []
    
    for pred in predictions_group:
        split_info = pred['split_info']
        split_x, split_y = split_info['split_position']
        split_w, split_h = split_info['split_size']
        
        # 解码预测的mask
        if 'segmentation' in pred:
            if isinstance(pred['segmentation'], dict):
                # RLE格式
                rle = pred['segmentation']
                binary_mask = mask_util.decode(rle)
            else:
                # Polygon格式，需要转换为mask
                from PIL import Image, ImageDraw
                img = Image.new('L', (split_w, split_h), 0)
                for polygon in pred['segmentation']:
                    polygon_coords = [(polygon[i], polygon[i+1]) for i in range(0, len(polygon), 2)]
                    ImageDraw.Draw(img).polygon(polygon_coords, outline=1, fill=1)
                binary_mask = np.array(img)
        else:
            # 如果没有segmentation，使用bbox创建mask
            bbox = pred['bbox']
            binary_mask = np.zeros((split_h, split_w), dtype=np.uint8)
            x, y, w, h = bbox
            binary_mask[y:y+h, x:x+w] = 1
        
        # 将预测mask映射到原始图像坐标
        if binary_mask.shape[0] > 0 and binary_mask.shape[1] > 0:
            # 确保不超出边界
            end_y = min(split_y + binary_mask.shape[0], original_h)
            end_x = min(split_x + binary_mask.shape[1], original_w)
            actual_h = end_y - split_y
            actual_w = end_x - split_x
            
            if actual_h > 0 and actual_w > 0:
                # 调整mask尺寸以匹配实际区域
                resized_mask = binary_mask[:actual_h, :actual_w]
                
                # 将mask添加到合并的mask中
                region = merged_mask[split_y:end_y, split_x:end_x]
                # 只在原mask为0的地方添加新的实例
                new_instances = (resized_mask > 0) & (region == 0)
                merged_mask[split_y:end_y, split_x:end_x][new_instances] = instance_id
                
                if np.sum(new_instances) > 0:  # 只有当实际添加了像素时才记录
                    # 计算合并后的bbox
                    y_indices, x_indices = np.where(merged_mask == instance_id)
                    if len(y_indices) > 0:
                        x_min, x_max = np.min(x_indices), np.max(x_indices)
                        y_min, y_max = np.min(y_indices), np.max(y_indices)
                        width = x_max - x_min + 1
                        height = y_max - y_min + 1
                        
                        merged_annotation = {
                            'id': instance_id,
                            'bbox': [int(x_min), int(y_min), int(width), int(height)],
                            'area': int(np.sum(merged_mask == instance_id)),
                            'category_id': pred.get('category_id', 1),
                            'score': pred.get('score', 1.0),
                            'original_split_info': split_info
                        }
                        merged_annotations.append(merged_annotation)
                        instance_id += 1
    
    return merged_mask, merged_annotations


def save_merged_results(merged_results, output_dir):
    """保存合并后的结果"""
    os.makedirs(output_dir, exist_ok=True)
    
    # 保存合并后的mask图像
    masks_dir = os.path.join(output_dir, 'merged_masks')
    os.makedirs(masks_dir, exist_ok=True)
    
    # 保存COCO格式的合并标注
    merged_annotations = []
    merged_images = []
    
    for idx, (original_image_path, (merged_mask, annotations)) in enumerate(merged_results.items()):
        # 保存mask图像
        image_name = Path(original_image_path).stem
        mask_filename = f"{image_name}_merged_mask.png"
        mask_path = os.path.join(masks_dir, mask_filename)
        
        # 将mask转换为可视化格式（每个实例不同颜色）
        max_instances = np.max(merged_mask)
        if max_instances > 0:
            # 创建彩色mask用于可视化
            colored_mask = np.zeros((merged_mask.shape[0], merged_mask.shape[1], 3), dtype=np.uint8)
            for instance_id in range(1, max_instances + 1):
                color = [(instance_id * 50) % 255, (instance_id * 100) % 255, (instance_id * 150) % 255]
                colored_mask[merged_mask == instance_id] = color
            
            cv.imwrite(mask_path, colored_mask)
            
            # 也保存原始mask（用于进一步处理）
            raw_mask_path = os.path.join(masks_dir, f"{image_name}_raw_mask.png")
            cv.imwrite(raw_mask_path, merged_mask.astype(np.uint16))
        
        # 添加到COCO格式
        # 读取原始图像以获取尺寸
        original_img = cv.imread(original_image_path)
        if original_img is not None:
            h, w = original_img.shape[:2]
            merged_images.append({
                'id': idx + 1,
                'file_name': Path(original_image_path).name,
                'width': w,
                'height': h,
                'original_path': original_image_path
            })
            
            # 更新标注的image_id
            for ann in annotations:
                ann['image_id'] = idx + 1
            merged_annotations.extend(annotations)
    
    # 保存COCO格式文件
    coco_output = {
        'images': merged_images,
        'annotations': merged_annotations,
        'categories': [
            {
                'id': 1,
                'name': 'granule',
                'supercategory': 'object'
            }
        ],
        'type': 'instances'
    }
    
    coco_path = os.path.join(output_dir, 'merged_predictions.json')
    with open(coco_path, 'w', encoding='utf-8') as f:
        json.dump(coco_output, f, indent=2, ensure_ascii=False)
    
    print(f"Merged results saved to: {output_dir}")
    print(f"Merged masks saved to: {masks_dir}")
    print(f"COCO annotations saved to: {coco_path}")
    print(f"Total merged images: {len(merged_images)}")
    print(f"Total merged annotations: {len(merged_annotations)}")


def merge_split_predictions(mapping_path, predictions_path, output_dir, original_images_dir=None):
    """主函数：合并切分图像的预测结果"""
    print("Loading split mapping...")
    split_mapping = load_split_mapping(mapping_path)
    
    print("Loading predictions...")
    predictions = load_predictions(predictions_path)
    
    print("Grouping predictions by original image...")
    image_groups = group_predictions_by_original_image(predictions, split_mapping)
    
    print(f"Found predictions for {len(image_groups)} original images")
    
    merged_results = {}
    
    for original_image_path, predictions_group in image_groups.items():
        print(f"Merging predictions for: {original_image_path}")
        
        # 尝试读取原始图像以获取尺寸
        if os.path.exists(original_image_path):
            original_img = cv.imread(original_image_path)
            if original_img is not None:
                original_shape = original_img.shape
            else:
                print(f"Warning: Cannot read original image {original_image_path}")
                continue
        else:
            # 如果原始路径不存在，尝试从其他地方查找
            if original_images_dir and os.path.exists(original_images_dir):
                image_name = Path(original_image_path).name
                alt_path = os.path.join(original_images_dir, image_name)
                if os.path.exists(alt_path):
                    original_img = cv.imread(alt_path)
                    if original_img is not None:
                        original_shape = original_img.shape
                        original_image_path = alt_path  # 更新路径
                    else:
                        print(f"Warning: Cannot read image from alternative path {alt_path}")
                        continue
                else:
                    print(f"Warning: Image not found in alternative directory {alt_path}")
                    continue
            else:
                print(f"Warning: Original image not found {original_image_path}")
                continue
        
        merged_mask, merged_annotations = merge_predictions_for_image(predictions_group, original_shape)
        merged_results[original_image_path] = (merged_mask, merged_annotations)
        
        print(f"  Merged {len(predictions_group)} split predictions into {len(merged_annotations)} instances")
    
    print("Saving merged results...")
    save_merged_results(merged_results, output_dir)
    
    return merged_results


def init():
    """初始化参数"""
    parser = arg.ArgumentParser(description='Merge split image predictions back to original images')
    
    parser.add_argument('-m', '--mapping', dest='mapping_path', required=True,
                       help='Path to split_mapping.json file')
    parser.add_argument('-p', '--predictions', dest='predictions_path', required=True,
                       help='Path to predictions file (COCO format)')
    parser.add_argument('-o', '--output', dest='output_dir', required=True,
                       help='Output directory for merged results')
    parser.add_argument('--original-images-dir', dest='original_images_dir',
                       help='Alternative directory to search for original images if paths in mapping are invalid')
    
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    cfg = init()
    
    # 创建输出目录
    os.makedirs(cfg.output_dir, exist_ok=True)
    
    # 合并预测结果
    results = merge_split_predictions(
        cfg.mapping_path, 
        cfg.predictions_path, 
        cfg.output_dir,
        cfg.original_images_dir
    )
    
    print(f"\nMerging completed! Results saved to: {cfg.output_dir}")