# filepath: d:\Gitspace\ipa_full\iPA\ipa\processing\segmentation\segmentation_sxt\model_sxt_mito_instance\split_regions\split_sxt_mito_images.py
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


def error(msg='Unknown'):
    print('\nERROR: ' + msg)
    exit(1)


def _fix_path(p):
    """确保路径以斜杠结尾"""
    if p[-1] not in ['/', '\\']:
        return p + os.sep
    else:
        return p


def get_data_ids():
    """获取所有数据ID"""
    base_path = r"D:\Gitspace\ipa_full\iPA\data\sxt_images\mrcslice_output_denoise"
    data_ids = []
    
    if os.path.exists(base_path):
        for item in os.listdir(base_path):
            item_path = os.path.join(base_path, item)
            if os.path.isdir(item_path):
                data_ids.append(item)
    
    return sorted(data_ids)


def load_slice_mapping():
    """加载切片映射信息"""
    mapping_path = r"D:\Gitspace\ipa_full\iPA\data\sxt_images\label_slice\slice_mapping_3d_to_2d.json"
    
    if not os.path.exists(mapping_path):
        print(f"Warning: Slice mapping file not found: {mapping_path}")
        return {}
    
    with open(mapping_path, 'r', encoding='utf-8') as f:
        mapping = json.load(f)
    
    # 创建从dataset_name到映射信息的字典
    dataset_mapping = {}
    for key, value in mapping.items():
        dataset_name = value['dataset_name']
        dataset_mapping[dataset_name] = value
    
    return dataset_mapping


def load_dataset_info():
    """加载数据集信息，包括需要翻转的数据集列表"""
    info_path = r"D:\Gitspace\ipa_full\iPA\data\sxt_images\label_slice\dataset_info.json"
    
    if not os.path.exists(info_path):
        print(f"Warning: Dataset info file not found: {info_path}")
        return {}
    
    with open(info_path, 'r', encoding='utf-8') as f:
        dataset_info = json.load(f)
    
    return dataset_info


def load_3d_mito_label(data_id, dataset_info=None):
    """加载3D mito标签文件，支持数据集别名映射"""
    # 检查是否有别名映射
    actual_data_id = data_id
    if dataset_info:
        aliases = dataset_info.get('dataset_aliases', {})
        if data_id in aliases:
            actual_data_id = aliases[data_id]
            print(f"    Using alias mapping: {data_id} -> {actual_data_id}")
    
    # 修改为mito标签文件路径
    label_path = rf"D:\Gitspace\ipa_full\iPA\data\sxt_images\label_slice\{actual_data_id}\combine_{actual_data_id}_mito_instances.tif"
    
    if not os.path.exists(label_path):
        print(f"Warning: 3D mito label file not found: {label_path}")
        return None
    
    try:
        # 加载3D mito标签文件
        label_3d = tifffile.imread(label_path)
        print(f"    Loaded 3D mito label: {label_3d.shape}")
        return label_3d
    except Exception as e:
        print(f"    Error loading 3D mito label: {e}")
        return None


def extract_slice_number(filename):
    """从文件名中提取切片号"""
    # 匹配文件名中的数字，通常在最后
    match = re.search(r'_(\d+)\.png$', filename)
    if match:
        return int(match.group(1))
    return None


def get_mask_from_3d_mito_label(label_3d, slice_mapping, filename, data_id, dataset_info):
    """从3D mito标签中提取对应切片的2D mask"""
    if label_3d is None or slice_mapping is None:
        return None
    
    # 提取切片号
    slice_num = extract_slice_number(filename)
    if slice_num is None:
        print(f"    Warning: Could not extract slice number from {filename}")
        return None
    
    # 检查是否有别名映射，获取实际的映射数据集ID
    mapping_data_id = data_id
    if dataset_info:
        aliases = dataset_info.get('dataset_aliases', {})
        if data_id in aliases:
            mapping_data_id = aliases[data_id]
            print(f"    Using mapping for alias: {data_id} -> {mapping_data_id}")
    
    # 获取映射信息
    if mapping_data_id not in slice_mapping:
        print(f"    Warning: No mapping found for dataset {mapping_data_id}")
        print(f"    Available datasets in mapping: {list(slice_mapping.keys())}")
        return None
    
    mapping_info = slice_mapping[mapping_data_id]
    slice_files = mapping_info.get('slice_files', [])
    print(f"    Mapping dataset {mapping_data_id} has {len(slice_files)} slice files")
    
    # 找到对应的切片文件在列表中的索引
    # 先尝试直接匹配
    target_filename = filename
    slice_index = None
    
    for i, slice_file in enumerate(slice_files):
        if slice_file == target_filename:
            slice_index = i
            break
    
    # 如果直接匹配失败，尝试大小写转换匹配
    if slice_index is None:
        # 尝试将 Ins -> INS
        alt_filename = filename.replace('_Ins_', '_INS_')
        for i, slice_file in enumerate(slice_files):
            if slice_file == alt_filename:
                slice_index = i
                print(f"    Found match with case conversion: {filename} -> {alt_filename}")
                break
    
    # 如果还是找不到，尝试更灵活的匹配（基于切片号）
    if slice_index is None:
        print(f"    Trying slice number matching for slice {slice_num}")
        for i, slice_file in enumerate(slice_files):
            file_slice_num = extract_slice_number(slice_file)
            if file_slice_num == slice_num:
                slice_index = i
                print(f"    Found match by slice number: {filename} (slice {slice_num}) -> {slice_file}")
                break
    
    # 如果还找不到，并且有别名映射，尝试基于数据集ID的文件名模式匹配
    if slice_index is None and mapping_data_id != data_id:
        print(f"    Trying ID pattern matching...")
        
        # 方案1: 直接替换整个data_id
        pattern_filename = filename.replace(data_id, mapping_data_id)
        print(f"    Trying pattern 1: {pattern_filename}")
        for i, slice_file in enumerate(slice_files):
            if slice_file == pattern_filename:
                slice_index = i
                print(f"    Found match with ID substitution: {filename} -> {pattern_filename}")
                break
        
        # 方案2: 替换 _data_id_ 模式
        if slice_index is None:
            pattern_filename = filename.replace(f"_{data_id}_", f"_{mapping_data_id}_")
            print(f"    Trying pattern 2: {pattern_filename}")
            for i, slice_file in enumerate(slice_files):
                if slice_file == pattern_filename:
                    slice_index = i
                    print(f"    Found match with ID pattern substitution: {filename} -> {pattern_filename}")
                    break
        
        # 方案3: 基于切片号和数据集ID重构文件名
        if slice_index is None:
            # 尝试从映射文件中找到相同切片号的文件
            for i, slice_file in enumerate(slice_files):
                file_slice_num = extract_slice_number(slice_file)
                if file_slice_num == slice_num:
                    slice_index = i
                    print(f"    Found match by reconstructed pattern: {filename} -> {slice_file}")
                    break
    
    if slice_index is None:
        print(f"    Warning: Slice file {filename} not found in mapping (tried all matching strategies)")
        print(f"    Looking for slice number: {slice_num}")
        print(f"    Original dataset: {data_id}, Mapping dataset: {mapping_data_id}")
        print(f"    Sample slice files from mapping:")
        for i, slice_file in enumerate(slice_files[:3]):
            print(f"      {i}: {slice_file}")
        return None
    
    # 从3D标签中提取对应的2D切片
    if slice_index >= label_3d.shape[0]:
        print(f"    Warning: Slice index {slice_index} out of range for 3D mito label shape {label_3d.shape}")
        return None
    
    mask_2d = label_3d[slice_index]
    
    # 检查是否需要翻转 - 使用原始data_id检查翻转规则
    need_flip = dataset_info.get('dataneedflipped', [])
    if data_id in need_flip:
        mask_2d = np.flipud(mask_2d)  # 上下翻转
        print(f"    Applied vertical flip to mito mask for dataset {data_id}")
    
    print(f"    Extracted 2D mito mask from slice {slice_index}, shape: {mask_2d.shape}")
    
    return mask_2d


def get_image_mask_pairs(data_id):
    """获取指定data_id下的图像文件列表，使用3D mito标签生成mask"""
    image_base = rf"D:\Gitspace\ipa_full\iPA\data\sxt_images\mrcslice_output_denoise\{data_id}\slices data"
    
    pairs = []
    
    if not os.path.exists(image_base):
        print(f"Warning: Image path not found for data_id {data_id}: {image_base}")
        return pairs
    
    # 获取所有图像文件
    image_files = [f for f in os.listdir(image_base) if f.endswith('.png')]
    
    for img_file in image_files:
        img_path = os.path.join(image_base, img_file)
        pairs.append((img_path, None, img_file))  # mask_path暂时为None，稍后从3D mito标签生成
    
    return pairs


def split_image_with_3d_mito_mask(img_path, label_3d, slice_mapping, data_id, filename, cfg):
    """使用3D mito标签切分单个图像"""
    # 读取图像
    img = cv.imread(img_path)
    
    if img is None:
        print(f"Error: Cannot read {img_path}")
        return []
    
    # 加载数据集信息
    dataset_info = load_dataset_info()
    
    # 从3D mito标签生成2D mask
    mask = get_mask_from_3d_mito_label(label_3d, slice_mapping, filename, data_id, dataset_info)
    
    if mask is None:
        print(f"    Warning: Could not generate mito mask for {filename}")
        return []
    
    raw_h, raw_w = img.shape[:2]
    print(f"    Original image size: {raw_w}x{raw_h}")
    print(f"    Mito mask shape: {mask.shape}")
    
    # 确保mask和图像尺寸匹配
    if mask.shape != (raw_h, raw_w):
        print(f"    Warning: Mito mask shape {mask.shape} doesn't match image shape {(raw_h, raw_w)}")
        # 尝试调整mask尺寸
        mask = cv.resize(mask.astype(np.uint8), (raw_w, raw_h), interpolation=cv.INTER_NEAREST)
        print(f"    Resized mito mask to: {mask.shape}")
    
    # 使用原始split_img.py的逻辑计算切分尺寸
    w, h = fix_wh(cfg, raw_w, raw_h)
    print(f"    Calculated split size: {w}x{h}, Grid: {cfg.xgrid}x{cfg.ygrid}")
    
    # 计算切分点 - 使用原始算法
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
    
    print(f"    X points: {x_points}")
    print(f"    Y points: {y_points}")
    
    split_results = []
    
    for i, x in enumerate(x_points):
        for j, y in enumerate(y_points):
            # 确保不超出图像边界
            x_end = min(x + w, raw_w)
            y_end = min(y + h, raw_h)
            
            print(f"    Split {i},{j}: ({x},{y}) to ({x_end},{y_end}), size: {x_end-x}x{y_end-y}")
            
            # 切分图像和mask
            split_img = img[y:y_end, x:x_end]
            split_mask = mask[y:y_end, x:x_end]
            
            # 检查mask是否包含有效区域
            mask_area = np.sum(split_mask > 0)
            print(f"    Mito mask area: {mask_area}, threshold: {cfg.min_mask_area}")
            
            if cfg.min_mask_area > 0:
                if mask_area < cfg.min_mask_area:
                    print(f"    Skipped due to small mito mask area")
                    continue
            
            split_info = {
                'original_image': img_path,
                'original_mask': f"3D_mito_label_slice_{extract_slice_number(filename)}",
                'split_position': (x, y),
                'split_size': (x_end - x, y_end - y),
                'grid_index': (i, j),
                'slice_index': extract_slice_number(filename)
            }
            
            split_results.append((split_img, split_mask, split_info))
    
    return split_results


def fix_wh(cfg, raw_w, raw_h):
    """计算切分尺寸 - 参考原始split_img.py的逻辑"""
    if cfg.overlap == -1:
        # 无重叠模式：确保切分尺寸不小于指定尺寸，且能完整覆盖原图
        w = int(max(cfg.width, ceil(raw_w / cfg.xgrid)))
        h = int(max(cfg.height, ceil(raw_h / cfg.ygrid)))
    else:
        # 重叠模式：根据重叠比例计算切分尺寸
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


def save_split_results(split_results, output_dir, data_id, img_filename, img_id_start=1):
    """保存切分结果，包括COCO格式标注"""
    img_output_dir = os.path.join(output_dir, 'images', data_id)
    mask_output_dir = os.path.join(output_dir, 'masks', data_id)
    
    os.makedirs(img_output_dir, exist_ok=True)
    os.makedirs(mask_output_dir, exist_ok=True)
    
    base_name = os.path.splitext(img_filename)[0]
    saved_files = []
    coco_images = []
    coco_annotations = []
    
    current_img_id = img_id_start
    
    for idx, (split_img, split_mask, split_info) in enumerate(split_results):
        # 生成文件名
        split_name = f"{base_name}_split_{idx:03d}.png"
        
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


def process_data_id(data_id, cfg):
    """处理单个data_id"""
    print(f"Processing data_id: {data_id}")
    
    # 加载数据集信息（包含别名映射）
    dataset_info = load_dataset_info()
    
    # 加载切片映射
    slice_mapping = load_slice_mapping()
    
    # 加载3D mito标签 - 传递dataset_info来处理别名映射
    label_3d = load_3d_mito_label(data_id, dataset_info)
    
    if label_3d is None:
        print(f"Cannot load 3D mito label for {data_id}, skipping...")
        return [], [], []
    
    pairs = get_image_mask_pairs(data_id)
    if not pairs:
        print(f"No valid image-mask pairs found for {data_id}")
        return [], [], []
    
    print(f"Found {len(pairs)} image files")
    
    all_results = []
    all_coco_images = []
    all_coco_annotations = []
    current_img_id = 1
    
    # 如果是测试模式，只处理第一个图像
    if hasattr(cfg, 'test') and cfg.test:
        pairs = pairs[:1]
        print("Test mode: processing only first image")
    
    for img_path, mask_path, img_filename in pairs:
        print(f"  Processing: {img_filename}")
        
        # 使用3D mito标签切分图像
        split_results = split_image_with_3d_mito_mask(img_path, label_3d, slice_mapping, data_id, img_filename, cfg)
        
        if split_results:
            saved_files, coco_images, coco_annotations, current_img_id = save_split_results(
                split_results, cfg.output_dir, data_id, img_filename, current_img_id)
            all_results.extend(saved_files)
            all_coco_images.extend(coco_images)
            all_coco_annotations.extend(coco_annotations)
            print(f"    Generated {len(split_results)} splits")
        else:
            print(f"    No valid splits generated")
    
    return all_results, all_coco_images, all_coco_annotations


def process_all_data(cfg):
    """处理所有数据"""
    data_ids = get_data_ids()
    
    if not data_ids:
        error("No data directories found!")
    
    print(f"Found {len(data_ids)} data directories: {data_ids}")
    
    # 如果是测试模式，只处理第一个数据目录
    if hasattr(cfg, 'test') and cfg.test:
        data_ids = data_ids[:1]
        print("Test mode: processing only first data directory")
    
    # 限制处理的数据数量
    if cfg.limit > 0:
        data_ids = data_ids[:cfg.limit]
        print(f"Limited to first {cfg.limit} data directories")
    
    all_results = []
    all_coco_images = []
    all_coco_annotations = []
    summary = {}
    
    # 全局图像ID计数器，确保ID唯一
    global_img_id = 1
    
    for data_id in data_ids:
        results, coco_images, coco_annotations = process_data_id(data_id, cfg)
        
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
                'name': 'mitochondria',
                'supercategory': 'organelle'
            }
        ],
        'type': 'instances'
    }
    
    # 保存COCO格式标注文件
    coco_path = os.path.join(cfg.output_dir, 'split_mito_instances.json')
    with open(coco_path, 'w', encoding='utf-8') as f:
        json.dump(coco_dataset, f, indent=2, ensure_ascii=False)
    
    # 保存切分信息映射文件（用于反推）
    split_mapping = {}
    for result in all_results:
        split_info = result['split_info']
        image_id = result['coco_image_id']
        split_mapping[image_id] = {
            'original_image': split_info['original_image'],
            'split_position': split_info['split_position'],
            'split_size': split_info['split_size'],
            'grid_index': split_info['grid_index'],
            'slice_index': split_info.get('slice_index')
        }
    
    mapping_path = os.path.join(cfg.output_dir, 'split_mito_mapping.json')
    with open(mapping_path, 'w', encoding='utf-8') as f:
        json.dump(split_mapping, f, indent=2, ensure_ascii=False)
    
    # 保存处理摘要
    summary_path = os.path.join(cfg.output_dir, 'split_mito_summary.json')
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
                'min_mask_area': cfg.min_mask_area
            }
        }, f, indent=2, ensure_ascii=False)
    
    print(f"\nProcessing completed!")
    print(f"Total splits generated: {len(all_results)}")
    print(f"COCO dataset saved to: {coco_path}")
    print(f"Split mapping saved to: {mapping_path}")
    print(f"Summary saved to: {summary_path}")
    
    return all_results


def show_config(cfg):
    """显示配置信息"""
    print("=" * 50)
    print("SXT Mito Image Split Configuration")
    print("=" * 50)
    print(f"Output directory: {cfg.output_dir}")
    print(f"Split grid: {cfg.xgrid} x {cfg.ygrid}")
    print(f"Split size: {cfg.width} x {cfg.height}")
    if cfg.overlap != -1:
        print(f"Overlap: {cfg.overlap}")
    else:
        print("No overlap")
    print(f"Minimum mask area: {cfg.min_mask_area}")
    if cfg.limit > 0:
        print(f"Limit to first {cfg.limit} data directories")
    print("=" * 50)


def init():
    """初始化参数"""
    parser = arg.ArgumentParser(description='Split SXT images and corresponding mito masks')
    
    # 基本参数 - 使用新的输出目录避免与ISG数据冲突
    parser.add_argument('-o', '--output', dest='output_dir', 
                       default=r'D:\Gitspace\ipa_full\iPA\data\sxt_images\split_mito_img',
                       help='Output directory for split mito images')
    parser.add_argument('-l', '--limit', dest='limit', default=0, type=int,
                       help='Limit number of data directories to process (0 = all)')
    parser.add_argument('--test', dest='test', action='store_true',
                       help='Test mode: only process first image of first data directory')
    
    # 切分参数 - 使用原始split_img.py的默认值
    parser.add_argument('--width', dest='width', default=100, type=int,
                       help='Split width')
    parser.add_argument('--height', dest='height', default=150, type=int,
                       help='Split height')
    parser.add_argument('--xgrid', dest='xgrid', default=4, type=int,
                       help='Number of splits in x direction')
    parser.add_argument('--ygrid', dest='ygrid', default=5, type=int,
                       help='Number of splits in y direction')
    parser.add_argument('--overlap', dest='overlap', default=-1, type=float,
                       help='Overlap ratio between splits (-1 for no overlap)')
    
    # 过滤参数
    parser.add_argument('--min-mask-area', dest='min_mask_area', default=100, type=int,
                       help='Minimum mask area to keep the split')
    
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    cfg = init()
    show_config(cfg)
    
    # 创建输出目录
    os.makedirs(cfg.output_dir, exist_ok=True)
    
    # 处理所有数据
    results = process_all_data(cfg)
    
    print(f"\nAll done! Check mito results in: {cfg.output_dir}")