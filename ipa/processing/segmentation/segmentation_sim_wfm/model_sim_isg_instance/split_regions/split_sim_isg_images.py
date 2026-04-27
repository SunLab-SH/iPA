import os
import cv2 as cv
import numpy as np
import argparse as arg
from math import ceil
import json
from pathlib import Path
import tifffile
import re
import pycocotools.mask as mask_util
from itertools import groupby
import sys


def error(msg='Unknown'):
    print('\nERROR: ' + msg)
    exit(1)


def _fix_path(p):
    """确保路径以斜杠结尾"""
    if p[-1] not in ['/', '\\']:
        return p + os.sep
    else:
        return p


def get_sim_data_ids():
    """获取SIM数据的所有可用ID"""
    # SIM文件命名规则
    specified_files = [
        "_0-2-1-SIM", "_0-2-2-SIM", "_0-3-1-SIM", "_0-3-1-1-SIM", "_0-3-3-SIM",
        "_0-3-5-SIM", "_0-4-1-1-SIM", "_0-4-1-SIM", "_0-4-2-SIM",
        "_5-2-1-SIM", "_5-2-2-SIM", "_5-2-3-1-SIM", "_5-2-3-SIM",
        "_5-3-2-SIM", "_5-3-3-SIM", "_5-3-4-1-SIM", "_5-3-4-SIM",
        "_30-2-1-SIM", "_30-2-2-SIM", "_30-2-3-SIM", "_30-2-4-SIM",
        "_30-3-2-SIM", "_30-3-3-SIM", "_30-4-1-SIM", "_30-4-3-SIM",
        "_30-4-4-SIM", "_30-4-5-SIM"
    ]
    
    # 转换为不带下划线的ID格式
    data_ids = [file_id.strip('_') for file_id in specified_files]
    return sorted(data_ids)


def load_sim_data(img_path, mask_path):
    """加载SIM图像和ISG mask数据"""
    try:
        # 加载图像和mask - SIM数据通常是3D TIFF
        img = tifffile.imread(img_path)
        mask = tifffile.imread(mask_path)
        
        print(f"    Original image shape: {img.shape}")
        print(f"    Original mask shape: {mask.shape}")
        
        # SIM数据格式分析和修正
        if len(img.shape) == 3:
            # 检查各维度大小，确定合理的(H, W, Z)格式
            dims = img.shape
            print(f"    Dimensions: {dims[0]} x {dims[1]} x {dims[2]}")
            
            # 找出最可能是高度、宽度、Z轴的维度
            # 通常Z轴是最小的维度，高度和宽度相对较大且比较接近
            sorted_dims = sorted([(i, dim) for i, dim in enumerate(dims)], key=lambda x: x[1])
            
            print(f"    Sorted dimensions: {sorted_dims}")
            
            # 最小的维度很可能是Z轴
            z_axis = sorted_dims[0][0]
            remaining_axes = [i for i in range(3) if i != z_axis]
            
            print(f"    Identified Z-axis: {z_axis} (size: {dims[z_axis]})")
            print(f"    H,W axes: {remaining_axes} (sizes: {[dims[i] for i in remaining_axes]})")
            
            # 重新排列为 (H, W, Z) 格式
            if z_axis == 0:  # (Z, H, W) -> (H, W, Z)
                img = np.transpose(img, (1, 2, 0))
                mask = np.transpose(mask, (1, 2, 0))
                print(f"    Converted from (Z,H,W) to (H,W,Z)")
            elif z_axis == 1:  # (H, Z, W) -> (H, W, Z)
                img = np.transpose(img, (0, 2, 1))
                mask = np.transpose(mask, (0, 2, 1))
                print(f"    Converted from (H,Z,W) to (H,W,Z)")
            # z_axis == 2 的情况已经是 (H, W, Z) 格式，不需要转换
        
        print(f"    Final image shape: {img.shape}")
        print(f"    Final mask shape: {mask.shape}")
        print(f"    Image data type: {img.dtype}, range: {img.min():.2f} - {img.max():.2f}")
        
        # 打印ISG mask统计信息
        isg_pixels = np.sum(mask > 0)
        total_pixels = mask.size
        print(f"    ISG coverage: {isg_pixels}/{total_pixels} pixels ({isg_pixels/total_pixels*100:.2f}%)")
        
        return img, mask
    except Exception as e:
        print(f"    Error loading data: {e}")
        return None, None


def analyze_isg_distribution(cfg):
    """分析ISG分布，优化split参数"""
    print("Analyzing ISG distribution for optimal split parameters...")
    
    data_ids = get_sim_data_ids()[:3]  # 分析前3个数据
    
    total_coverage = 0
    total_slices = 0
    distribution_stats = []
    
    for data_id in data_ids:
        pairs = get_sim_image_mask_pairs(data_id, 
                                       raw_img_dir=cfg.raw_img_dir, 
                                       mask_dir=cfg.mask_dir)
        
        if pairs:
            img_path, mask_path, _ = pairs[0]
            img, mask = load_sim_data(img_path, mask_path)
            
            if img is not None and mask is not None:
                # 分析每个Z切片的ISG分布
                for z in range(mask.shape[2]):
                    mask_2d = mask[:, :, z]
                    isg_pixels = np.sum(mask_2d > 0)
                    
                    if isg_pixels > 0:
                        total_slices += 1
                        h, w = mask_2d.shape
                        coverage = isg_pixels / (h * w)
                        total_coverage += coverage
                        
                        # 找到ISG区域的边界框
                        y_indices, x_indices = np.where(mask_2d > 0)
                        if len(y_indices) > 0:
                            x_min, x_max = np.min(x_indices), np.max(x_indices)
                            y_min, y_max = np.min(y_indices), np.max(y_indices)
                            
                            distribution_stats.append({
                                'data_id': data_id,
                                'z_slice': z,
                                'coverage': coverage,
                                'bbox': (x_min, y_min, x_max - x_min + 1, y_max - y_min + 1),
                                'isg_pixels': isg_pixels
                            })
    
    if total_slices > 0:
        avg_coverage = total_coverage / total_slices
        print(f"\nISG Distribution Analysis:")
        print(f"  Average ISG coverage per slice: {avg_coverage*100:.2f}%")
        print(f"  Total slices with ISG: {total_slices}")
        
        # 分析边界框统计
        if distribution_stats:
            bbox_widths = [stat['bbox'][2] for stat in distribution_stats]
            bbox_heights = [stat['bbox'][3] for stat in distribution_stats]
            
            avg_bbox_w = np.mean(bbox_widths)
            avg_bbox_h = np.mean(bbox_heights)
            max_bbox_w = np.max(bbox_widths)
            max_bbox_h = np.max(bbox_heights)
            
            print(f"  ISG region bbox - avg: {avg_bbox_w:.0f}x{avg_bbox_h:.0f}, max: {max_bbox_w}x{max_bbox_h}")
            
            # 基于分析结果提供建议
            print(f"\nSplit Parameter Recommendations:")
            
            # 建议split尺寸
            if avg_coverage < 0.1:  # 覆盖率低于10%
                suggested_size = 256
                suggested_grid = 8
                print(f"  - Low ISG density detected")
            elif avg_coverage < 0.3:  # 覆盖率10-30%
                suggested_size = 384
                suggested_grid = 6
                print(f"  - Medium ISG density detected")
            else:  # 高密度
                suggested_size = 512
                suggested_grid = 4
                print(f"  - High ISG density detected")
            
            print(f"  - Suggested split size: {suggested_size}x{suggested_size}")
            print(f"  - Suggested grid: {suggested_grid}x{suggested_grid}")
            print(f"  - Current split size: {cfg.width}x{cfg.height}")
            print(f"  - Current grid: {cfg.xgrid}x{cfg.ygrid}")
            
            # 检查当前配置效率
            current_total_splits = cfg.xgrid * cfg.ygrid
            estimated_useful_splits = current_total_splits * avg_coverage
            
            print(f"\nCurrent Configuration Efficiency:")
            print(f"  - Total splits per slice: {current_total_splits}")
            print(f"  - Estimated useful splits: {estimated_useful_splits:.1f}")
            print(f"  - Efficiency: {estimated_useful_splits/current_total_splits*100:.1f}%")
            
            if estimated_useful_splits/current_total_splits < 0.3:
                print(f"  - WARNING: Low efficiency! Consider reducing grid size or increasing split size")


def test_sim_image_loading(cfg):
    """测试SIM图像加载和维度检测"""
    print("Testing SIM image loading...")
    
    data_ids = get_sim_data_ids()[:2]  # 只测试前2个
    
    for data_id in data_ids:
        print(f"\nTesting {data_id}:")
        pairs = get_sim_image_mask_pairs(data_id, 
                                       raw_img_dir=cfg.raw_img_dir, 
                                       mask_dir=cfg.mask_dir)
        
        if pairs:
            img_path, mask_path, _ = pairs[0]
            print(f"  Image file: {os.path.basename(img_path)}")
            print(f"  Mask file: {os.path.basename(mask_path)}")
            
            # 直接读取原始数据，不做任何转换
            try:
                img_raw = tifffile.imread(img_path)
                mask_raw = tifffile.imread(mask_path)
                
                print(f"  Raw image shape: {img_raw.shape}")
                print(f"  Raw mask shape: {mask_raw.shape}")
                print(f"  Raw image dtype: {img_raw.dtype}")
                print(f"  Raw mask dtype: {mask_raw.dtype}")
                
                # 分析维度含义
                if len(img_raw.shape) == 3:
                    dims = img_raw.shape
                    print(f"  Analyzing dimensions: {dims}")
                    
                    # 典型的SIM图像应该是 1280x1280 左右的正方形，Z维度较小
                    # 所以我们期望看到两个相近的大数（宽高）和一个小数（Z）
                    
                    max_dim = max(dims)
                    min_dim = min(dims)
                    mid_dim = sum(dims) - max_dim - min_dim
                    
                    print(f"  Min: {min_dim}, Mid: {mid_dim}, Max: {max_dim}")
                    
                    # 如果最小维度明显小于其他两个，那它可能是Z
                    if min_dim < max_dim * 0.1:  # Z维度应该远小于空间维度
                        z_axis = dims.index(min_dim)
                        print(f"  Likely Z-axis position: {z_axis} (size: {min_dim})")
                        
                        if z_axis == 0:
                            print("  Format appears to be (Z, H, W)")
                        elif z_axis == 1:
                            print("  Format appears to be (H, Z, W)")
                        else:
                            print("  Format appears to be (H, W, Z)")
                    else:
                        print("  WARNING: Cannot clearly identify Z-axis!")
                
            except Exception as e:
                print(f"  Error loading: {e}")


def get_sim_image_mask_pairs(data_id, base_date="20220909", raw_img_dir=None, mask_dir=None):
    """获取指定data_id下的SIM图像和mask文件对"""
    
    if raw_img_dir is None:
        raw_img_dir = r'F:\Salilab\Projects\IPA-toolbox\raw_image\SIM_data\Singlecolor'
    if mask_dir is None:
        mask_dir = r'F:\Salilab\Projects\IPA-toolbox\raw_image\SIM_data\from_bing\SIM\ISG_spot_seg'
    
    pairs = []
    
    # 构建文件路径
    img_path = os.path.join(raw_img_dir, f"{base_date}_{data_id}_Actin.tif")
    mask_path = os.path.join(mask_dir, f"{base_date}_{data_id}_ISG_seg.tif")
    
    # 检查文件是否存在
    if os.path.exists(img_path) and os.path.exists(mask_path):
        pairs.append((img_path, mask_path, f"{base_date}_{data_id}"))
        print(f"Found pair for {data_id}: image and mask files exist")
    else:
        missing = []
        if not os.path.exists(img_path):
            missing.append(f"image: {os.path.basename(img_path)}")
        if not os.path.exists(mask_path):
            missing.append(f"mask: {os.path.basename(mask_path)}")
        print(f"Warning: Missing files for {data_id}: {', '.join(missing)}")
    
    return pairs


def split_sim_image_with_3d_mask(img_path, mask_path, data_id, cfg):
    """使用3D SIM数据切分图像"""
    
    # 加载SIM数据
    img, mask = load_sim_data(img_path, mask_path)
    
    if img is None or mask is None:
        print(f"Error: Cannot load data for {data_id}")
        return []
    
    # 确保尺寸匹配
    if img.shape != mask.shape:
        print(f"    Warning: Shape mismatch - Image: {img.shape}, Mask: {mask.shape}")
        # 尝试调整mask尺寸以匹配图像
        if len(img.shape) == 3 and len(mask.shape) == 3:
            if img.shape[:2] == mask.shape[:2] and img.shape[2] != mask.shape[2]:
                # Z维度不同，取较小的
                min_z = min(img.shape[2], mask.shape[2])
                img = img[:, :, :min_z]
                mask = mask[:, :, :min_z]
                print(f"    Adjusted to common Z dimension: {img.shape}")
        
        if img.shape != mask.shape:
            print(f"    Cannot resolve shape mismatch, skipping {data_id}")
            return []
    
    print(f"    Processing 3D volume: {img.shape}")
    
    # 归一化图像到0-1范围
    if img.max() > img.min():
        img_norm = (img - img.min()) / (img.max() - img.min())
    else:
        img_norm = np.zeros_like(img)
    
    # 转换为uint8并创建3通道图像
    img_uint8 = (img_norm * 255).astype(np.uint8)
    
    # 将ISG mask转换为二值mask
    mask_binary = (mask > 0).astype(np.uint8)
    
    split_results = []
    
    # 找到包含ISG数据的Z切片
    z_with_isg = np.any(mask > 0, axis=(0, 1))
    z_indices = np.where(z_with_isg)[0]
    
    if len(z_indices) == 0:
        print(f"    Warning: No ISG data found in any slice for {data_id}")
        return []
    
    print(f"    Found ISG data in {len(z_indices)} slices: Z {z_indices[0]} to {z_indices[-1]}")
    
    # 处理每个包含ISG数据的Z切片
    for z_idx in z_indices:
        img_2d = img_uint8[:, :, z_idx]
        mask_2d = mask_binary[:, :, z_idx]
        
        # 创建3通道图像用于模型兼容性
        img_3ch = np.stack([img_2d, img_2d, img_2d], axis=2)
        
        raw_h, raw_w = img_2d.shape
        print(f"    Processing slice Z={z_idx}, size: {raw_w}x{raw_h}")
        
        # 使用split逻辑计算切分尺寸
        w, h = fix_wh(cfg, raw_w, raw_h)
        print(f"    Calculated split size: {w}x{h}, Grid: {cfg.xgrid}x{cfg.ygrid}")
        
        # 计算切分点
        if cfg.xgrid != 1:
            x_step = (raw_w - w) / (cfg.xgrid - 1)
        else:
            x_step = 0
        x_points = [int(i * x_step) for i in range(cfg.xgrid)]
        
        if cfg.ygrid != 1:
            y_step = (raw_h - h) / (cfg.ygrid - 1)
        else:
            y_step = 0
        y_points = [int(i * y_step) for i in range(cfg.ygrid)]
        
        # 执行切分
        for i, x in enumerate(x_points):
            for j, y in enumerate(y_points):
                # 确保不超出图像边界
                x_end = min(x + w, raw_w)
                y_end = min(y + h, raw_h)
                
                # 切分图像和mask
                split_img = img_3ch[y:y_end, x:x_end]
                split_mask = mask_2d[y:y_end, x:x_end]
                
                # 检查mask是否包含有效区域
                mask_area = np.sum(split_mask > 0)
                
                if cfg.min_mask_area > 0:
                    if mask_area < cfg.min_mask_area:
                        print(f"    Skipped split {i},{j} Z={z_idx} due to small mask area: {mask_area}")
                        continue
                
                split_info = {
                    'original_image': img_path,
                    'original_mask': mask_path,
                    'split_position': (x, y),
                    'split_size': (x_end - x, y_end - y),
                    'grid_index': (i, j),
                    'slice_index': z_idx,
                    'data_id': data_id
                }
                
                split_results.append((split_img, split_mask, split_info))
                print(f"    Generated split {i},{j} Z={z_idx}, mask area: {mask_area}")
    
    print(f"    Total splits generated for {data_id}: {len(split_results)}")
    return split_results


def fix_wh(cfg, raw_w, raw_h):
    """计算切分尺寸"""
    if cfg.overlap == -1:
        # 无重叠模式
        w = int(max(cfg.width, ceil(raw_w / cfg.xgrid)))
        h = int(max(cfg.height, ceil(raw_h / cfg.ygrid)))
    else:
        # 重叠模式
        w = int((2 * raw_w) / (2 * cfg.xgrid + cfg.overlap - cfg.overlap * cfg.xgrid))
        h = int((2 * raw_h) / (2 * cfg.ygrid + cfg.overlap - cfg.overlap * cfg.ygrid))
    return (w, h)


def binary_mask_to_rle(binary_mask):
    """将二值mask转换为RLE格式"""
    rle = {'counts': [], 'size': list(binary_mask.shape)}
    counts = rle.get('counts')
    for i, (value, elements) in enumerate(groupby(binary_mask.ravel(order='F'))):
        if i == 0 and value == 1:
            counts.append(0)
        counts.append(len(list(elements)))
    return rle


def mask_to_coco_annotations(mask, image_id, category_id=1):
    """将mask转换为COCO格式的标注"""
    annotations = []
    
    # 找到所有唯一的标签（排除背景0）
    unique_labels = np.unique(mask)
    unique_labels = unique_labels[unique_labels > 0]
    
    annotation_id = 1
    for label in unique_labels:
        # 创建当前标签的二值mask
        binary_mask = (mask == label).astype(np.uint8)
        
        # 计算面积
        area = np.sum(binary_mask)
        if area == 0:
            continue
        
        # 计算边界框
        y_indices, x_indices = np.where(binary_mask)
        if len(y_indices) == 0:
            continue
            
        x_min, x_max = np.min(x_indices), np.max(x_indices)
        y_min, y_max = np.min(y_indices), np.max(y_indices)
        width = x_max - x_min + 1
        height = y_max - y_min + 1
        bbox = [int(x_min), int(y_min), int(width), int(height)]
        
        # 转换为RLE格式
        rle = binary_mask_to_rle(binary_mask)
        
        annotation = {
            'id': annotation_id,
            'image_id': image_id,
            'category_id': category_id,
            'segmentation': rle,
            'area': int(area),
            'bbox': bbox,
            'iscrowd': 0
        }
        
        annotations.append(annotation)
        annotation_id += 1
    
    return annotations


def save_split_results(split_results, output_dir, data_id, img_id_start=1):
    """保存切分结果，包括COCO格式标注"""
    img_output_dir = os.path.join(output_dir, 'images', data_id)
    mask_output_dir = os.path.join(output_dir, 'masks', data_id)
    
    os.makedirs(img_output_dir, exist_ok=True)
    os.makedirs(mask_output_dir, exist_ok=True)
    
    saved_files = []
    coco_images = []
    coco_annotations = []
    
    current_img_id = img_id_start
    
    for idx, (split_img, split_mask, split_info) in enumerate(split_results):
        # 生成文件名
        grid_i, grid_j = split_info['grid_index']
        slice_z = split_info['slice_index']
        split_name = f"{data_id}_z{slice_z:03d}_split_{grid_i}_{grid_j}.png"
        
        img_save_path = os.path.join(img_output_dir, split_name)
        mask_save_path = os.path.join(mask_output_dir, split_name)
        
        # 保存图像和mask
        cv.imwrite(img_save_path, split_img)
        cv.imwrite(mask_save_path, split_mask)
        
        # 创建COCO图像信息
        h, w = split_img.shape[:2]
        coco_image = {
            'id': current_img_id,
            'file_name': split_name,
            'width': w,
            'height': h,
            'image_name': f'{data_id}_{current_img_id}'
        }
        coco_images.append(coco_image)
        
        # 创建COCO标注信息
        annotations = mask_to_coco_annotations(split_mask, current_img_id)
        coco_annotations.extend(annotations)
        
        saved_files.append({
            'image_path': img_save_path,
            'mask_path': mask_save_path,
            'split_info': split_info,
            'coco_image_id': current_img_id
        })
        
        current_img_id += 1
    
    return saved_files, coco_images, coco_annotations, current_img_id


def process_sim_data_id(data_id, cfg):
    """处理单个SIM data_id"""
    print(f"Processing SIM data_id: {data_id}")
    
    # 获取图像和mask文件对
    pairs = get_sim_image_mask_pairs(data_id, 
                                   raw_img_dir=cfg.raw_img_dir, 
                                   mask_dir=cfg.mask_dir)
    
    if not pairs:
        print(f"No valid image-mask pairs found for {data_id}")
        return [], [], []
    
    print(f"Found {len(pairs)} image-mask pairs")
    
    all_results = []
    all_coco_images = []
    all_coco_annotations = []
    current_img_id = 1
    
    for img_path, mask_path, filename in pairs:
        print(f"  Processing: {filename}")
        
        # 使用3D SIM数据切分图像
        split_results = split_sim_image_with_3d_mask(img_path, mask_path, data_id, cfg)
        
        if split_results:
            saved_files, coco_images, coco_annotations, current_img_id = save_split_results(
                split_results, cfg.output_dir, data_id, current_img_id)
            all_results.extend(saved_files)
            all_coco_images.extend(coco_images)
            all_coco_annotations.extend(coco_annotations)
            print(f"    Generated {len(split_results)} splits")
        else:
            print(f"    No valid splits generated")
    
    return all_results, all_coco_images, all_coco_annotations


def process_all_sim_data(cfg):
    """处理所有SIM数据"""
    data_ids = get_sim_data_ids()
    
    if not data_ids:
        error("No SIM data IDs found!")
    
    print(f"Found {len(data_ids)} SIM data IDs: {data_ids}")
    
    # 如果是测试模式，只处理第一个数据目录
    if hasattr(cfg, 'test') and cfg.test:
        data_ids = data_ids[:1]
        print("Test mode: processing only first data ID")
    
    # 限制处理的数据数量
    if cfg.limit > 0:
        data_ids = data_ids[:cfg.limit]
        print(f"Limited to first {cfg.limit} data IDs")
    
    all_results = []
    all_coco_images = []
    all_coco_annotations = []
    summary = {}
    
    # 全局图像ID计数器，确保ID唯一
    global_img_id = 1
    
    for data_id in data_ids:
        results, coco_images, coco_annotations = process_sim_data_id(data_id, cfg)
        
        # 更新图像ID为全局唯一
        for i, img in enumerate(coco_images):
            old_id = img['id']
            img['id'] = global_img_id + i
            
            # 更新对应的标注
            for ann in coco_annotations:
                if ann['image_id'] == old_id:
                    ann['image_id'] = global_img_id + i
        
        global_img_id += len(coco_images)
        
        all_results.extend(results)
        all_coco_images.extend(coco_images)
        all_coco_annotations.extend(coco_annotations)
        summary[data_id] = len(results)
    
    # 创建完整的COCO格式数据集
    coco_dataset = {
        'images': all_coco_images,
        'annotations': all_coco_annotations,
        'categories': [
            {
                'id': 1,
                'name': 'isg',
                'supercategory': 'vesicle'
            }
        ],
        'type': 'instances'
    }
    
    # 保存COCO格式标注文件
    coco_path = os.path.join(cfg.output_dir, 'split_sim_isg_instances.json')
    with open(coco_path, 'w', encoding='utf-8') as f:
        json.dump(coco_dataset, f, indent=2, ensure_ascii=False)
    
    # 保存切分信息映射文件（用于反推）
    split_mapping = {}
    for result in all_results:
        split_info = result['split_info']
        image_id = result['coco_image_id']
        split_mapping[str(image_id)] = {
            'original_image': split_info['original_image'],
            'original_mask': split_info['original_mask'],
            'split_position': [int(x) for x in split_info['split_position']],  # 转换为int
            'split_size': [int(x) for x in split_info['split_size']],  # 转换为int
            'grid_index': [int(x) for x in split_info['grid_index']],  # 转换为int
            'slice_index': int(split_info.get('slice_index', 0)),  # 转换为int
            'data_id': str(split_info.get('data_id', ''))  # 确保是字符串
        }
    
    mapping_path = os.path.join(cfg.output_dir, 'split_sim_isg_mapping.json')
    with open(mapping_path, 'w', encoding='utf-8') as f:
        json.dump(split_mapping, f, indent=2, ensure_ascii=False)
    
    # 保存处理摘要
    summary_path = os.path.join(cfg.output_dir, 'split_sim_isg_summary.json')
    with open(summary_path, 'w', encoding='utf-8') as f:
        json.dump({
            'total_splits': len(all_results),
            'data_summary': summary,
            'config': {
                'xgrid': cfg.xgrid,
                'ygrid': cfg.ygrid,
                'width': cfg.width,
                'height': cfg.height,
                'overlap': cfg.overlap,
                'min_mask_area': cfg.min_mask_area,
                'raw_img_dir': cfg.raw_img_dir,
                'mask_dir': cfg.mask_dir
            }
        }, f, indent=2, ensure_ascii=False)
    
    print(f"\nSIM ISG processing completed!")
    print(f"Total splits generated: {len(all_results)}")
    print(f"COCO dataset saved to: {coco_path}")
    print(f"Split mapping saved to: {mapping_path}")
    print(f"Summary saved to: {summary_path}")
    
    return all_results


def show_config(cfg):
    """显示配置信息"""
    print("=" * 50)
    print("SIM ISG Image Split Configuration")
    print("=" * 50)
    print(f"Raw image directory: {cfg.raw_img_dir}")
    print(f"Mask directory: {cfg.mask_dir}")
    print(f"Output directory: {cfg.output_dir}")
    print(f"Split grid: {cfg.xgrid} x {cfg.ygrid}")
    print(f"Split size: {cfg.width} x {cfg.height}")
    if cfg.overlap != -1:
        print(f"Overlap: {cfg.overlap}")
    else:
        print("No overlap")
    print(f"Minimum mask area: {cfg.min_mask_area}")
    if cfg.limit > 0:
        print(f"Limit to first {cfg.limit} data IDs")
    print("=" * 50)


def init():
    """初始化参数"""
    parser = arg.ArgumentParser(description='Split SIM images and corresponding ISG masks')
    
    # 基本参数
    parser.add_argument('-o', '--output', dest='output_dir', 
                       default=r'D:\Gitspace\ipa_full\data\SIM\split_isg_img',
                       help='Output directory for split ISG images')
    parser.add_argument('-l', '--limit', dest='limit', default=0, type=int,
                       help='Limit number of data IDs to process (0 = all)')
    parser.add_argument('--analyze-distribution', dest='analyze_distribution', action='store_true',
                       help='Analyze ISG distribution and suggest optimal split parameters')
    parser.add_argument('--test-image-loading', dest='test_image_loading', action='store_true',
                       help='Test image loading and dimension detection only')
    parser.add_argument('--test', dest='test', action='store_true',
                       help='Test mode: only process first data ID')
    
    # SIM数据路径
    parser.add_argument('--raw-img-dir', dest='raw_img_dir',
                       default=r'F:\Salilab\Projects\IPA-toolbox\raw_image\SIM_data\Singlecolor',
                       help='Raw SIM images directory')
    parser.add_argument('--mask-dir', dest='mask_dir',
                       default=r'F:\Salilab\Projects\IPA-toolbox\raw_image\SIM_data\from_bing\SIM\ISG_spot_seg',
                       help='ISG masks directory')
    
    # 切分参数 - 固定尺寸用于模型训练
    parser.add_argument('--width', dest='width', default=256, type=int,
                       help='Split width (fixed size for model training)')
    parser.add_argument('--height', dest='height', default=256, type=int,
                       help='Split height (fixed size for model training)')
    parser.add_argument('--xgrid', dest='xgrid', default=10, type=int,
                       help='Number of splits in x direction')
    parser.add_argument('--ygrid', dest='ygrid', default=10, type=int,
                       help='Number of splits in y direction')
    parser.add_argument('--overlap', dest='overlap', default=-1, type=float,
                       help='Overlap ratio between splits (-1 for no overlap)')
    
    # 过滤参数
    parser.add_argument('--min-mask-area', dest='min_mask_area', default=50, type=int,
                       help='Minimum mask area to keep the split')
    
    args = parser.parse_args()
    args = parser.parse_args()
    return args


def analyze_split_results(output_dir):
    """分析已生成的切分结果"""
    print("Analyzing generated split results...")
    
    # 加载COCO数据
    coco_path = os.path.join(output_dir, 'split_sim_isg_instances.json')
    mapping_path = os.path.join(output_dir, 'split_sim_isg_mapping.json')
    summary_path = os.path.join(output_dir, 'split_sim_isg_summary.json')
    
    if not os.path.exists(coco_path):
        print(f"COCO file not found: {coco_path}")
        return
    
    with open(coco_path, 'r') as f:
        coco_data = json.load(f)
    
    with open(summary_path, 'r') as f:
        summary = json.load(f)
    
    # 基本统计
    total_images = len(coco_data['images'])
    total_annotations = len(coco_data['annotations'])
    
    print(f"\n=== Split Results Analysis ===")
    print(f"Total images: {total_images}")
    print(f"Total annotations: {total_annotations}")
    print(f"Average annotations per image: {total_annotations/total_images:.2f}")
    
    # 按data_id分析
    print(f"\nPer data_id breakdown:")
    for data_id, count in summary['data_summary'].items():
        print(f"  {data_id}: {count} splits")
    
    # 分析Z切片分布
    z_slice_count = {}
    for img in coco_data['images']:
        filename = img['file_name']
        # 提取Z切片信息: 0-2-1-SIM_z000_split_2_3.png
        if '_z' in filename:
            z_part = filename.split('_z')[1].split('_')[0]
            z_idx = int(z_part)
            z_slice_count[z_idx] = z_slice_count.get(z_idx, 0) + 1
    
    print(f"\nZ-slice distribution:")
    for z_idx in sorted(z_slice_count.keys()):
        print(f"  Z-slice {z_idx:2d}: {z_slice_count[z_idx]:3d} splits")
    
    # 分析标注质量
    annotation_areas = []
    bbox_sizes = []
    
    for ann in coco_data['annotations']:
        if 'area' in ann and 'bbox' in ann:
            annotation_areas.append(ann['area'])
            bbox = ann['bbox']
            bbox_sizes.append((bbox[2], bbox[3]))  # width, height
    
    if annotation_areas:
        print(f"\nAnnotation statistics:")
        print(f"  Total annotations: {len(annotation_areas)}")
        print(f"  Area - min: {min(annotation_areas)}, max: {max(annotation_areas)}, avg: {np.mean(annotation_areas):.1f}")
        
        bbox_widths = [size[0] for size in bbox_sizes]
        bbox_heights = [size[1] for size in bbox_sizes]
        print(f"  BBox width - min: {min(bbox_widths)}, max: {max(bbox_widths)}, avg: {np.mean(bbox_widths):.1f}")
        print(f"  BBox height - min: {min(bbox_heights)}, max: {max(bbox_heights)}, avg: {np.mean(bbox_heights):.1f}")
    
    # 检查数据完整性
    print(f"\nData integrity check:")
    
    # 检查图像文件是否存在
    missing_images = 0
    for img in coco_data['images'][:10]:  # 只检查前10个
        data_id = img['file_name'].split('_')[0] + '-' + img['file_name'].split('_')[1] + '-' + img['file_name'].split('_')[2]
        img_path = os.path.join(output_dir, 'images', data_id, img['file_name'])
        if not os.path.exists(img_path):
            missing_images += 1
    
    if missing_images == 0:
        print(f"  ✅ Image files check passed (sample of 10)")
    else:
        print(f"  ❌ Found {missing_images} missing image files")
    
    # 验证COCO格式
    images_with_annotations = set()
    for ann in coco_data['annotations']:
        images_with_annotations.add(ann['image_id'])
    
    images_without_annotations = []
    for img in coco_data['images']:
        if img['id'] not in images_with_annotations:
            images_without_annotations.append(img['id'])
    
    print(f"  Images with annotations: {len(images_with_annotations)}/{total_images}")
    if images_without_annotations:
        print(f"  ⚠️  Images without annotations: {len(images_without_annotations)}")
        print(f"     (This is normal if some splits contain no ISG)")
    
    # 建议
    print(f"\n=== Recommendations ===")
    efficiency = len(images_with_annotations) / total_images
    print(f"Data efficiency: {efficiency*100:.1f}% (splits with ISG / total splits)")
    
    if efficiency < 0.3:
        print("❌ Low efficiency! Consider:")
        print("   - Reducing grid size (e.g., 8x8 instead of 10x10)")
        print("   - Increasing split size (e.g., 384x384 instead of 256x256)")
        print("   - Increasing min_mask_area threshold")
    elif efficiency < 0.6:
        print("⚠️  Medium efficiency. Could be optimized:")
        print("   - Fine-tune grid size or min_mask_area")
    else:
        print("✅ Good efficiency! Data is well balanced.")
    
    avg_ann_per_img = total_annotations / len(images_with_annotations) if images_with_annotations else 0
    if avg_ann_per_img < 1:
        print("⚠️  Low average annotations per image. ISG objects might be very small or sparse.")
    elif avg_ann_per_img > 5:
        print("⚠️  High average annotations per image. Consider checking for over-segmentation.")
    
    return coco_data, summary


if __name__ == "__main__":
    cfg = init()
    
    # 如果是测试图像加载模式
    if cfg.test_image_loading:
        print("Testing SIM image loading and dimension detection...")
        test_sim_image_loading(cfg)
        exit(0)
    
    # 如果是分析ISG分布模式
    if cfg.analyze_distribution:
        print("Analyzing ISG distribution...")
        analyze_isg_distribution(cfg)
        exit(0)
    
    # 如果是分析已生成结果模式
    if len(sys.argv) > 1 and sys.argv[1] == '--analyze-results':
        analyze_split_results(cfg.output_dir)
        exit(0)
    
    show_config(cfg)
    
    # 检查输入目录是否存在
    if not os.path.exists(cfg.raw_img_dir):
        error(f"Raw image directory not found: {cfg.raw_img_dir}")
    
    if not os.path.exists(cfg.mask_dir):
        error(f"Mask directory not found: {cfg.mask_dir}")
    
    # 创建输出目录
    os.makedirs(cfg.output_dir, exist_ok=True)
    
    # 处理所有SIM数据
    results = process_all_sim_data(cfg)
    
    print(f"\nAll done! Check SIM ISG results in: {cfg.output_dir}")