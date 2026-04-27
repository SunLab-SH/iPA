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


class SXTDataset(Dataset):
    """SXT切分图像数据集"""
    
    def __init__(self, data_root, coco_file, phase='train', transforms=None):
        self.data_root = data_root
        self.phase = phase
        self.transforms = transforms
        
        # 加载COCO格式标注
        with open(coco_file, 'r', encoding='utf-8') as f:
            self.coco_data = json.load(f)
        
        # 创建图像ID到图像信息的映射
        self.images = {img['id']: img for img in self.coco_data['images']}
        
        # 创建图像ID到标注的映射
        self.annotations = {}
        for ann in self.coco_data['annotations']:
            img_id = ann['image_id']
            if img_id not in self.annotations:
                self.annotations[img_id] = []
            self.annotations[img_id].append(ann)
        
        # 获取所有图像ID列表
        self.image_ids = list(self.images.keys())
        
        print(f"加载 {phase} 数据集: {len(self.image_ids)} 张图像")
        print(f"COCO数据总数: {len(self.coco_data['images'])} 张图像")
        if len(self.image_ids) > 0:
            print(f"图像ID范围: {min(self.image_ids)} - {max(self.image_ids)}")
        
        # 检查图像ID的连续性
        all_ids = [img['id'] for img in self.coco_data['images']]
        max_id = max(all_ids)
        expected_count = max_id
        actual_count = len(self.image_ids)
        if actual_count < expected_count:
            print(f"警告: 图像ID不连续，最大ID={max_id}，但只有{actual_count}个唯一ID")
        
    def __len__(self):
        return len(self.image_ids)
    
    def __getitem__(self, idx):
        img_id = self.image_ids[idx]
        img_info = self.images[img_id]
        
        # 构建图像路径
        img_path = self._get_image_path(img_info['file_name'])
        
        # 读取图像
        image = Image.open(img_path).convert('RGB')
        
        # 获取标注
        anns = self.annotations.get(img_id, [])
        
        # 准备target
        target = {}
        if anns:  # 如果有标注（训练阶段）
            boxes = []
            labels = []
            masks = []
            areas = []
            
            for ann in anns:
                # 边界框
                bbox = ann['bbox']  # [x, y, width, height]
                boxes.append([bbox[0], bbox[1], bbox[0] + bbox[2], bbox[1] + bbox[3]])
                
                # 标签
                labels.append(ann['category_id'])
                
                # 面积
                areas.append(ann['area'])
                
                # 掩码
                if 'segmentation' in ann:
                    if isinstance(ann['segmentation'], dict):  # RLE格式
                        rle = ann['segmentation']
                        # 确保RLE格式正确
                        if isinstance(rle.get('counts'), list):
                            # 如果counts是list，需要转换为压缩格式
                            compressed_rle = mask_util.frPyObjects(rle, rle['size'][0], rle['size'][1])
                            mask = mask_util.decode(compressed_rle)
                        else:
                            # 直接解码
                            mask = mask_util.decode(rle)
                    else:  # polygon格式
                        # 将polygon转换为mask
                        mask = np.zeros((img_info['height'], img_info['width']), dtype=np.uint8)
                        # 这里简化处理，实际可能需要更复杂的polygon转mask逻辑
                    
                    # 确保mask尺寸正确
                    expected_h, expected_w = img_info['height'], img_info['width']
                    if mask.shape != (expected_h, expected_w):
                        # 只在第一次遇到时警告，避免大量重复输出
                        if not hasattr(self, '_size_warning_shown'):
                            print(f"警告: 检测到mask尺寸不匹配，将自动调整。mask: {mask.shape}, 图像: {(expected_h, expected_w)}")
                            self._size_warning_shown = True
                        mask = cv.resize(mask, (expected_w, expected_h), interpolation=cv.INTER_NEAREST)
                    
                    masks.append(mask)
            
            if boxes:
                target['boxes'] = torch.as_tensor(boxes, dtype=torch.float32)
                target['labels'] = torch.as_tensor(labels, dtype=torch.int64)
                target['masks'] = torch.as_tensor(np.array(masks), dtype=torch.uint8)
                target['area'] = torch.as_tensor(areas, dtype=torch.float32)
                target['iscrowd'] = torch.zeros((len(anns),), dtype=torch.int64)
            else:
                # 空标注
                target['boxes'] = torch.zeros((0, 4), dtype=torch.float32)
                target['labels'] = torch.zeros((0,), dtype=torch.int64)
                target['masks'] = torch.zeros((0, img_info['height'], img_info['width']), dtype=torch.uint8)
                target['area'] = torch.zeros((0,), dtype=torch.float32)
                target['iscrowd'] = torch.zeros((0,), dtype=torch.int64)
        
        target['image_id'] = torch.tensor([img_id])
        
        # 应用transforms
        if self.transforms:
            image = self.transforms(image)
        
        return image, target
    
    def _get_image_path(self, filename):
        """根据文件名构建完整路径"""
        # 首先尝试在所有子文件夹中查找文件
        images_base = os.path.join(self.data_root, 'images')
        
        # 遍历所有子文件夹查找文件
        if os.path.exists(images_base):
            for subdir in os.listdir(images_base):
                subdir_path = os.path.join(images_base, subdir)
                if os.path.isdir(subdir_path):
                    img_path = os.path.join(subdir_path, filename)
                    if os.path.exists(img_path):
                        return img_path
        
        # 如果直接查找失败，尝试更灵活的文件名匹配
        print(f"Warning: 直接查找失败，尝试模糊匹配: {filename}")
        
        # 提取关键信息进行模糊匹配
        if os.path.exists(images_base):
            for subdir in os.listdir(images_base):
                subdir_path = os.path.join(images_base, subdir)
                if os.path.isdir(subdir_path):
                    # 列出该目录下的所有文件
                    for existing_file in os.listdir(subdir_path):
                        if existing_file.endswith('.png'):
                            # 简单的相似度匹配：如果大部分字符相同
                            if self._files_similar(filename, existing_file):
                                matched_path = os.path.join(subdir_path, existing_file)
                                print(f"找到相似文件: {filename} -> {existing_file}")
                                return matched_path
        
        # 最后尝试直接在images目录下查找
        direct_path = os.path.join(self.data_root, 'images', filename)
        if os.path.exists(direct_path):
            return direct_path
            
        raise FileNotFoundError(f"Cannot find image: {filename}")
    
    def _files_similar(self, file1, file2, threshold=0.8):
        """检查两个文件名是否相似"""
        # 移除扩展名
        name1 = os.path.splitext(file1)[0]
        name2 = os.path.splitext(file2)[0]
        
        # 简单的相似度计算：共同字符数 / 较长字符串长度
        common_chars = sum(1 for a, b in zip(name1, name2) if a == b)
        max_len = max(len(name1), len(name2))
        similarity = common_chars / max_len if max_len > 0 else 0
        
        return similarity >= threshold


def get_transform(train=False):
    """获取数据预处理transforms"""
    transforms = []
    transforms.append(T.PILToTensor())
    transforms.append(T.ConvertImageDtype(torch.float))
    if train:
        # 训练时可以添加数据增强
        pass
    return T.Compose(transforms)


def get_model_instance_segmentation(num_classes):
    """创建Mask R-CNN模型"""
    # 加载预训练的Mask R-CNN模型
    model = torchvision.models.detection.maskrcnn_resnet50_fpn(pretrained=True)
    
    # 获取分类器的输入特征数
    in_features = model.roi_heads.box_predictor.cls_score.in_features
    # 替换预训练的头部
    model.roi_heads.box_predictor = FastRCNNPredictor(in_features, num_classes)
    
    # 获取mask分类器的输入特征数
    in_features_mask = model.roi_heads.mask_predictor.conv5_mask.in_channels
    hidden_layer = 256
    # 替换mask预测器
    model.roi_heads.mask_predictor = MaskRCNNPredictor(in_features_mask, hidden_layer, num_classes)
    
    return model


def collate_fn(batch):
    """自定义collate函数"""
    return tuple(zip(*batch))


def predict_single_image(model, image_path, device, threshold=0.5):
    """预测单张图像"""
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


def visualize_predictions(prediction, save_path=None, show=True):
    """可视化预测结果"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))
    
    image = prediction['image']
    scores = prediction['scores']
    boxes = prediction['boxes']
    labels = prediction['labels']
    masks = prediction['masks']
    
    # 显示原图
    ax1.imshow(image)
    ax1.set_title('Original Image')
    ax1.axis('off')
    
    # 显示预测结果
    ax2.imshow(image)
    
    # 绘制边界框和掩码
    for i, (box, score, label, mask) in enumerate(zip(boxes, scores, labels, masks)):
        if len(mask.shape) == 3:
            mask = mask[0]  # 移除channel维度
        
        # 绘制掩码
        mask_colored = np.zeros_like(image)
        color = np.random.rand(3)
        mask_colored[mask > 0.5] = color * 255
        ax2.imshow(mask_colored, alpha=0.5)
        
        # 绘制边界框
        x1, y1, x2, y2 = box
        rect = patches.Rectangle((x1, y1), x2-x1, y2-y1, 
                               linewidth=2, edgecolor=color, facecolor='none')
        ax2.add_patch(rect)
        
        # 添加标签和分数
        ax2.text(x1, y1-5, f'Class {label}: {score:.3f}', 
                bbox=dict(boxstyle="round,pad=0.3", facecolor=color, alpha=0.7),
                fontsize=8, color='white')
    
    ax2.set_title(f'Predictions (Found {len(scores)} objects)')
    ax2.axis('off')
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"可视化结果已保存到: {save_path}")
    
    if show:
        plt.show()
    
    plt.close()


def save_predictions_coco_format(predictions_list, image_infos, output_file):
    """将预测结果保存为COCO格式"""
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
                'segmentation': rle
            }
            results.append(result)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    
    print(f"预测结果已保存到: {output_file}")


def main():
    parser = argparse.ArgumentParser(description='SXT Mask R-CNN Prediction')
    parser.add_argument('--data_root', default=r'D:\Gitspace\ipa_full\iPA\data\sxt_images\split_img',
                       help='数据根目录')
    parser.add_argument('--coco_file', default=r'D:\Gitspace\ipa_full\iPA\data\sxt_images\split_img\split_instances.json',
                       help='COCO格式标注文件')
    parser.add_argument('--model_path', default=None,
                       help='训练好的模型路径（如果为None则使用预训练模型）')
    parser.add_argument('--output_dir', default=r'D:\Gitspace\ipa_full\iPA\data\sxt_images\predictions',
                       help='预测结果输出目录')
    parser.add_argument('--threshold', type=float, default=0.5,
                       help='置信度阈值')
    parser.add_argument('--num_samples', type=int, default=10,
                       help='预测样本数量')
    parser.add_argument('--visualize', action='store_true',
                       help='是否可视化预测结果')
    parser.add_argument('--save_coco', action='store_true',
                       help='是否保存COCO格式的预测结果')
    
    args = parser.parse_args()
    
    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 设备
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'使用设备: {device}')
    
    # 创建模型
    num_classes = 2  # 背景 + granule
    model = get_model_instance_segmentation(num_classes)
    
    # 加载训练好的模型（如果有）
    if args.model_path and os.path.exists(args.model_path):
        print(f'加载训练好的模型: {args.model_path}')
        model.load_state_dict(torch.load(args.model_path, map_location=device))
    else:
        print('使用预训练模型进行预测')
    
    model.to(device)
    model.eval()
    
    # 创建数据集
    dataset = SXTDataset(args.data_root, args.coco_file, phase='pred', 
                        transforms=get_transform(train=False))
    
    # 随机选择样本进行预测
    import random
    sample_indices = random.sample(range(len(dataset)), 
                                 min(args.num_samples, len(dataset)))
    
    predictions_list = []
    image_infos = []
    
    print(f'开始预测 {len(sample_indices)} 个样本...')
    
    for i, idx in enumerate(sample_indices):
        print(f'预测进度: {i+1}/{len(sample_indices)}')
        
        image, target = dataset[idx]
        img_id = target['image_id'].item()
        img_info = dataset.images[img_id]
        image_infos.append(img_info)
        
        # 构建图像路径
        image_path = dataset._get_image_path(img_info['file_name'])
        
        # 预测
        prediction = predict_single_image(model, image_path, device, args.threshold)
        predictions_list.append(prediction)
        
        print(f"图像: {img_info['file_name']}")
        print(f"检测到 {len(prediction['scores'])} 个对象")
        
        # 可视化
        if args.visualize:
            vis_path = os.path.join(args.output_dir, f'prediction_{i:03d}_{img_info["file_name"]}')
            visualize_predictions(prediction, save_path=vis_path, show=False)
    
    # 保存COCO格式预测结果
    if args.save_coco:
        coco_output = os.path.join(args.output_dir, 'predictions_coco_format.json')
        save_predictions_coco_format(predictions_list, image_infos, coco_output)
    
    print(f'\n预测完成！结果保存在: {args.output_dir}')


if __name__ == '__main__':
    main()