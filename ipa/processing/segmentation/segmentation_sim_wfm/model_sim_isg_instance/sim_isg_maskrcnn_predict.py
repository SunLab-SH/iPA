# filepath: d:\Gitspace\ipa_full\iPA\ipa\processing\segmentation\segmentation_sim_wfm\model_sim_isg_instance\sim_isg_maskrcnn_predict.py
import torch
import torchvision
from torchvision.models.detection.faster_rcnn import FastRCNNPredictor
from torchvision.models.detection.mask_rcnn import MaskRCNNPredictor
from torch.utils.data import Dataset, DataLoader
import torchvision.transforms as T
import os
import json
import cv2 as cv
import numpy as np
from PIL import Image
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pathlib import Path
import pycocotools.mask as mask_util
from sim_isg_maskrcnn_train import SIMISGDataset, get_transform, get_model_instance_segmentation, collate_fn


def predict_single_image(model, image_path, device, threshold=0.5):
    """预测单张SIM ISG图像"""
    model.eval()
    
    # 读取图像
    image = Image.open(image_path).convert('RGB')
    
    # 转换为tensor
    transform = get_transform(train=False)
    image_tensor = transform(image).unsqueeze(0).to(device)
    
    with torch.no_grad():
        predictions = model(image_tensor)
    
    # 处理预测结果
    pred = predictions[0]
    
    # 过滤低置信度的预测
    scores = pred['scores'].cpu().numpy()
    boxes = pred['boxes'].cpu().numpy()
    labels = pred['labels'].cpu().numpy()
    masks = pred['masks'].cpu().numpy()
    
    # 应用阈值
    keep = scores >= threshold
    scores = scores[keep]
    boxes = boxes[keep]
    labels = labels[keep]
    masks = masks[keep]
    
    return {
        'scores': scores,
        'boxes': boxes,
        'labels': labels,
        'masks': masks,
        'image': np.array(image)
    }


def visualize_isg_predictions(prediction, save_path=None, show=True):
    """可视化SIM ISG预测结果"""
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
    
    image = prediction['image']
    scores = prediction['scores']
    boxes = prediction['boxes']
    labels = prediction['labels']
    masks = prediction['masks']
    
    # 显示原图
    ax1.imshow(image)
    ax1.set_title('Original SIM Image')
    ax1.axis('off')
    
    # 显示预测结果 - 边界框
    ax2.imshow(image)
    for i, (box, score, label) in enumerate(zip(boxes, scores, labels)):
        # 绘制边界框
        x1, y1, x2, y2 = box
        color = 'red' if label == 1 else 'blue'  # ISG为红色
        rect = patches.Rectangle((x1, y1), x2-x1, y2-y1, 
                               linewidth=2, edgecolor=color, facecolor='none')
        ax2.add_patch(rect)
        
        # 添加标签和分数
        label_text = 'ISG' if label == 1 else f'Class {label}'
        ax2.text(x1, y1-5, f'{label_text}: {score:.3f}', 
                bbox=dict(boxstyle="round,pad=0.3", facecolor=color, alpha=0.7),
                fontsize=8, color='white')
    
    ax2.set_title(f'ISG Detection (Found {len(scores)} objects)')
    ax2.axis('off')
    
    # 显示掩码结果
    ax3.imshow(image)
    
    # 创建组合掩码显示
    combined_mask = np.zeros_like(image)
    for i, (mask, score, label) in enumerate(zip(masks, scores, labels)):
        if len(mask.shape) == 3:
            mask = mask[0]  # 移除channel维度
        
        # 为每个ISG实例使用不同颜色
        color = plt.cm.Set3(i / max(len(masks), 1))[:3]  # 获取RGB颜色
        mask_colored = np.zeros_like(image)
        mask_binary = mask > 0.5
        mask_colored[mask_binary] = [int(c * 255) for c in color]
        
        # 叠加到组合掩码
        alpha = 0.6
        combined_mask = np.where(mask_binary[..., np.newaxis], 
                               mask_colored * alpha + combined_mask * (1 - alpha), 
                               combined_mask)
    
    # 显示组合掩码
    if len(masks) > 0:
        ax3.imshow(combined_mask.astype(np.uint8), alpha=0.7)
    
    ax3.set_title(f'ISG Segmentation Masks')
    ax3.axis('off')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"ISG预测可视化结果已保存到: {save_path}")
    
    if show:
        plt.show()
    
    plt.close()


def save_isg_predictions_coco_format(predictions_list, image_infos, output_file):
    """将ISG预测结果保存为COCO格式"""
    results = []
    
    for predictions, img_info in zip(predictions_list, image_infos):
        img_id = img_info['id']
        
        for i, (score, box, label, mask) in enumerate(zip(
            predictions['scores'], predictions['boxes'], 
            predictions['labels'], predictions['masks'])):
            
            # 转换mask为RLE格式
            if len(mask.shape) == 3:
                mask = mask[0]
            binary_mask = (mask > 0.5).astype(np.uint8)
            rle = mask_util.encode(np.asfortranarray(binary_mask))
            rle['counts'] = rle['counts'].decode('utf-8')
            
            result = {
                'image_id': img_id,
                'category_id': int(label),
                'bbox': [float(box[0]), float(box[1]), 
                        float(box[2] - box[0]), float(box[3] - box[1])],
                'score': float(score),
                'segmentation': rle,
                'area': float(np.sum(binary_mask))
            }
            results.append(result)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    
    print(f"ISG预测结果已保存到: {output_file}")


def analyze_isg_predictions(predictions_list, image_infos):
    """分析ISG预测结果统计信息"""
    total_images = len(predictions_list)
    images_with_isg = 0
    total_isg_detected = 0
    confidence_scores = []
    isg_sizes = []
    
    for predictions, img_info in zip(predictions_list, image_infos):
        scores = predictions['scores']
        masks = predictions['masks']
        
        if len(scores) > 0:
            images_with_isg += 1
            total_isg_detected += len(scores)
            confidence_scores.extend(scores)
            
            # 计算ISG大小
            for mask in masks:
                if len(mask.shape) == 3:
                    mask = mask[0]
                isg_area = np.sum(mask > 0.5)
                isg_sizes.append(isg_area)
    
    print(f"\n=== ISG预测结果分析 ===")
    print(f"总图像数: {total_images}")
    print(f"检测到ISG的图像: {images_with_isg} ({images_with_isg/total_images*100:.1f}%)")
    print(f"总检测到的ISG数量: {total_isg_detected}")
    print(f"平均每张图像ISG数量: {total_isg_detected/total_images:.2f}")
    if images_with_isg > 0:
        print(f"有ISG图像的平均ISG数量: {total_isg_detected/images_with_isg:.2f}")
    
    if confidence_scores:
        print(f"置信度统计:")
        print(f"  平均置信度: {np.mean(confidence_scores):.3f}")
        print(f"  置信度范围: {np.min(confidence_scores):.3f} - {np.max(confidence_scores):.3f}")
    
    if isg_sizes:
        print(f"ISG大小统计 (像素):")
        print(f"  平均大小: {np.mean(isg_sizes):.1f}")
        print(f"  大小范围: {int(np.min(isg_sizes))} - {int(np.max(isg_sizes))}")
        print(f"  中位数大小: {np.median(isg_sizes):.1f}")


def batch_predict_sim_data(model, data_root, coco_file, device, output_dir, 
                          threshold=0.5, num_samples=50, data_ids=None):
    """批量预测SIM数据的所有data_id"""
    
    # 加载映射文件以了解原始数据结构
    mapping_file = os.path.join(data_root, 'split_sim_isg_mapping.json')
    if os.path.exists(mapping_file):
        with open(mapping_file, 'r', encoding='utf-8') as f:
            mapping_data = json.load(f)
        print(f"加载映射文件: {len(mapping_data)} 个split记录")
    else:
        print("未找到映射文件，使用默认预测方式")
        mapping_data = {}
    
    # 创建数据集
    dataset = SIMISGDataset(data_root, coco_file, phase='pred', 
                           transforms=get_transform(train=False))
    
    # 如果指定了data_ids，只预测这些数据
    if data_ids:
        print(f"指定预测data_ids: {data_ids}")
        # 过滤数据集，只保留指定data_id的图像
        filtered_indices = []
        for i, img_id in enumerate(dataset.image_ids):
            img_info = dataset.images[img_id]
            filename = img_info['file_name']
            # 从文件名中提取data_id
            for data_id in data_ids:
                if data_id in filename:
                    filtered_indices.append(i)
                    break
        
        if filtered_indices:
            dataset.image_ids = [dataset.image_ids[i] for i in filtered_indices]
            print(f"过滤后的数据集大小: {len(dataset.image_ids)}")
        else:
            print("警告: 未找到匹配指定data_id的图像")
    
    # 随机选择样本进行预测
    import random
    sample_indices = random.sample(range(len(dataset)), 
                                 min(num_samples, len(dataset)))
    
    predictions_list = []
    image_infos = []
    
    print(f'开始预测 SIM ISG {len(sample_indices)} 个样本...')
    
    for i, idx in enumerate(sample_indices):
        print(f'预测进度: {i+1}/{len(sample_indices)}')
        
        image, target = dataset[idx]
        img_id = target['image_id'].item()
        img_info = dataset.images[img_id]
        image_infos.append(img_info)
        
        # 构建图像路径
        image_path = dataset._get_image_path(img_info['file_name'])
        
        # 预测
        prediction = predict_single_image(model, image_path, device, threshold)
        predictions_list.append(prediction)
        
        print(f"图像: {img_info['file_name']}")
        print(f"检测到 {len(prediction['scores'])} 个ISG")
        if len(prediction['scores']) > 0:
            avg_confidence = np.mean(prediction['scores'])
            print(f"平均置信度: {avg_confidence:.3f}")
        
        # 可视化前几个样本
        if i < 5:  # 只可视化前5个样本
            vis_path = os.path.join(output_dir, f'isg_prediction_{i:03d}_{img_info["file_name"]}')
            visualize_isg_predictions(prediction, save_path=vis_path, show=False)
    
    return predictions_list, image_infos


def reconstruct_full_image_predictions(predictions_list, image_infos, mapping_data, output_dir):
    """将split预测结果重构为完整图像的预测"""
    print("\n=== 重构完整图像预测结果 ===")
    
    # 按原始图像分组
    original_images = {}
    
    for pred, img_info in zip(predictions_list, image_infos):
        img_id = str(img_info['id'])
        if img_id in mapping_data:
            mapping_info = mapping_data[img_id]
            original_image = mapping_info['original_image']
            data_id = mapping_info['data_id']
            slice_index = mapping_info['slice_index']
            split_position = mapping_info['split_position']
            
            key = f"{data_id}_{slice_index}"
            if key not in original_images:
                original_images[key] = {
                    'data_id': data_id,
                    'slice_index': slice_index,
                    'original_image': original_image,
                    'predictions': []
                }
            
            # 调整预测坐标到原始图像坐标系
            adjusted_pred = {
                'scores': pred['scores'].copy(),
                'labels': pred['labels'].copy(),
                'boxes': pred['boxes'].copy(),
                'masks': pred['masks'].copy(),
                'split_position': split_position,
                'split_info': mapping_info
            }
            
            # 调整边界框坐标
            if len(adjusted_pred['boxes']) > 0:
                adjusted_pred['boxes'][:, [0, 2]] += split_position[0]  # x坐标
                adjusted_pred['boxes'][:, [1, 3]] += split_position[1]  # y坐标
            
            original_images[key]['predictions'].append(adjusted_pred)
    
    print(f"找到 {len(original_images)} 个原始图像slice")
    
    # 保存重构结果
    reconstruction_file = os.path.join(output_dir, 'reconstructed_isg_predictions.json')
    with open(reconstruction_file, 'w', encoding='utf-8') as f:
        # 转换numpy数组为列表以便JSON序列化
        json_data = {}
        for key, value in original_images.items():
            json_data[key] = {
                'data_id': value['data_id'],
                'slice_index': value['slice_index'],
                'original_image': value['original_image'],
                'total_predictions': len(value['predictions']),
                'total_isg_detected': sum(len(p['scores']) for p in value['predictions'])
            }
        json.dump(json_data, f, indent=2, ensure_ascii=False)
    
    print(f"重构结果已保存到: {reconstruction_file}")
    return original_images


def main():
    parser = argparse.ArgumentParser(description='SIM ISG Mask R-CNN Prediction')
    parser.add_argument('--data_root', default=r'D:\Gitspace\ipa_full\data\SIM\split_isg_img',
                       help='SIM ISG数据根目录')
    parser.add_argument('--coco_file', default=r'D:\Gitspace\ipa_full\data\SIM\split_isg_img\split_sim_isg_instances.json',
                       help='COCO格式标注文件')
    parser.add_argument('--model_path', default=None,
                       help='训练好的ISG模型路径（如果为None则使用预训练模型）')
    parser.add_argument('--output_dir', default=r'D:\Gitspace\ipa_full\data\SIM\split_isg_img\predictions',
                       help='预测结果输出目录')
    parser.add_argument('--threshold', type=float, default=0.5,
                       help='ISG检测置信度阈值')
    parser.add_argument('--num_samples', type=int, default=50,
                       help='预测样本数量')
    parser.add_argument('--visualize', action='store_true',
                       help='是否可视化预测结果')
    parser.add_argument('--save_coco', action='store_true',
                       help='是否保存COCO格式的预测结果')
    parser.add_argument('--data_ids', nargs='+', default=None,
                       help='指定要预测的data_id列表，如：0-2-1-SIM 0-2-2-SIM')
    parser.add_argument('--reconstruct', action='store_true',
                       help='是否重构完整图像的预测结果')
    
    args = parser.parse_args()
    
    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 设备
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'使用设备: {device}')
    
    # 创建模型
    num_classes = 2  # 背景 + ISG
    model = get_model_instance_segmentation(num_classes)
    
    # 加载训练好的模型（如果有）
    if args.model_path and os.path.exists(args.model_path):
        print(f'加载训练好的ISG模型: {args.model_path}')
        model.load_state_dict(torch.load(args.model_path, map_location=device))
    else:
        print('使用预训练模型进行ISG预测（建议先训练模型）')
    
    model.to(device)
    model.eval()
    
    # 批量预测
    predictions_list, image_infos = batch_predict_sim_data(
        model, args.data_root, args.coco_file, device, args.output_dir,
        threshold=args.threshold, num_samples=args.num_samples, data_ids=args.data_ids)
    
    # 分析预测结果
    analyze_isg_predictions(predictions_list, image_infos)
    
    # 保存COCO格式预测结果
    if args.save_coco:
        coco_output = os.path.join(args.output_dir, 'isg_predictions_coco_format.json')
        save_isg_predictions_coco_format(predictions_list, image_infos, coco_output)
    
    # 重构完整图像预测结果
    if args.reconstruct:
        mapping_file = os.path.join(args.data_root, 'split_sim_isg_mapping.json')
        if os.path.exists(mapping_file):
            with open(mapping_file, 'r', encoding='utf-8') as f:
                mapping_data = json.load(f)
            reconstruct_full_image_predictions(predictions_list, image_infos, mapping_data, args.output_dir)
        else:
            print("未找到映射文件，跳过重构步骤")
    
    print(f'\nSIM ISG预测完成！结果保存在: {args.output_dir}')


if __name__ == '__main__':
    main()