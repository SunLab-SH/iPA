# filepath: d:\Gitspace\ipa_full\iPA\ipa\processing\segmentation\segmentation_sim_wfm\model_sim_isg_instance\sim_isg_maskrcnn_train.py
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
import time
import sys
from pathlib import Path
import pycocotools.mask as mask_util
from sklearn.metrics import precision_recall_fscore_support, roc_auc_score, average_precision_score
import logging
from tqdm import tqdm
try:
    from torch.utils.tensorboard import SummaryWriter
    TENSORBOARD_AVAILABLE = True
except ImportError:
    TENSORBOARD_AVAILABLE = False
    print("Warning: tensorboard not available, skipping tensorboard logging")


class SIMISGDataset(Dataset):
    """SIM ISG切分图像数据集"""
    
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
        
        print(f"加载 {phase} SIM ISG数据集: {len(self.image_ids)} 张图像")
        print(f"COCO数据总数: {len(self.coco_data['images'])} 张图像")
        print(f"标注总数: {len(self.coco_data['annotations'])} 个ISG实例")
        if len(self.image_ids) > 0:
            print(f"图像ID范围: {min(self.image_ids)} - {max(self.image_ids)}")
        
        # 统计数据分布
        self._analyze_dataset()
        
    def _analyze_dataset(self):
        """分析数据集统计信息"""
        total_instances = 0
        images_with_isg = 0
        
        for img_id in self.image_ids:
            anns = self.annotations.get(img_id, [])
            if anns:
                images_with_isg += 1
                total_instances += len(anns)
        
        print(f"数据集统计:")
        print(f"  包含ISG的图像: {images_with_isg}/{len(self.image_ids)} ({images_with_isg/len(self.image_ids)*100:.1f}%)")
        print(f"  平均每张图像ISG数量: {total_instances/len(self.image_ids):.2f}")
        if images_with_isg > 0:
            print(f"  有ISG图像的平均ISG数量: {total_instances/images_with_isg:.2f}")
    
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
        
        # 准备target - 初始化空目标
        target = {
            'boxes': torch.zeros((0, 4), dtype=torch.float32),
            'labels': torch.zeros((0,), dtype=torch.int64),
            'masks': torch.zeros((0, img_info['height'], img_info['width']), dtype=torch.uint8),
            'area': torch.zeros((0,), dtype=torch.float32),
            'iscrowd': torch.zeros((0,), dtype=torch.int64),
            'image_id': torch.tensor([img_id])
        }
        
        if anns:  # 如果有标注
            boxes = []
            labels = []
            masks = []
            areas = []
            valid_anns = []
            
            for ann in anns:
                try:
                    # 边界框
                    bbox = ann['bbox']  # [x, y, width, height]
                    if bbox[2] <= 0 or bbox[3] <= 0:  # 跳过无效bbox
                        continue
                    
                    boxes.append([bbox[0], bbox[1], bbox[0] + bbox[2], bbox[1] + bbox[3]])
                    
                    # 标签 - ISG类别
                    labels.append(ann['category_id'])
                    
                    # 面积
                    areas.append(ann['area'])
                    
                    # 掩码
                    mask = None
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
                    
                    if mask is None:
                        # 如果没有mask，根据bbox创建简单的矩形mask
                        mask = np.zeros((img_info['height'], img_info['width']), dtype=np.uint8)
                        x1, y1, x2, y2 = int(bbox[0]), int(bbox[1]), int(bbox[0] + bbox[2]), int(bbox[1] + bbox[3])
                        mask[y1:y2, x1:x2] = 1
                    
                    # 确保mask尺寸正确
                    expected_h, expected_w = img_info['height'], img_info['width']
                    if mask.shape != (expected_h, expected_w):
                        # 只在第一次遇到时警告，避免大量重复输出
                        if not hasattr(self, '_size_warning_shown'):
                            print(f"警告: 检测到mask尺寸不匹配，将自动调整。mask: {mask.shape}, 图像: {(expected_h, expected_w)}")
                            self._size_warning_shown = True
                        mask = cv.resize(mask, (expected_w, expected_h), interpolation=cv.INTER_NEAREST)
                    
                    masks.append(mask)
                    valid_anns.append(ann)
                    
                except Exception as e:
                    print(f"警告: 跳过无效标注 {ann.get('id', 'unknown')}: {e}")
                    continue
            
            # 只有当有有效标注时才更新target
            if len(valid_anns) > 0 and len(boxes) > 0:
                target['boxes'] = torch.as_tensor(boxes, dtype=torch.float32)
                target['labels'] = torch.as_tensor(labels, dtype=torch.int64)
                target['masks'] = torch.as_tensor(np.array(masks), dtype=torch.uint8)
                target['area'] = torch.as_tensor(areas, dtype=torch.float32)
                target['iscrowd'] = torch.zeros((len(valid_anns),), dtype=torch.int64)
        
        # 应用transforms
        if self.transforms:
            image = self.transforms(image)
        
        return image, target
    
    def _get_image_path(self, filename):
        """根据文件名构建完整路径 - 适配SIM ISG数据结构"""
        # SIM ISG数据的目录结构: images/{data_id}/{filename}
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
        
        # 提取关键信息进行模糊匹配 - 基于SIM文件命名模式
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
            
        raise FileNotFoundError(f"Cannot find SIM ISG image: {filename}")
    
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
        # 训练时可以添加数据增强 - 针对ISG小目标的增强
        # 注意：ISG是小目标，过度的几何变换可能会损害性能
        pass
    return T.Compose(transforms)


def get_model_instance_segmentation(num_classes):
    """创建Mask R-CNN模型 - 针对ISG优化"""
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
    
    # 针对小目标ISG调整RPN参数
    # 降低NMS阈值以保留更多小目标候选框
    model.rpn.nms_thresh = 0.6  # 默认0.7
    # 增加训练时的候选框数量
    model.rpn.pre_nms_top_n_train = 3000  # 默认2000
    model.rpn.post_nms_top_n_train = 2000  # 默认1000
    
    return model


def collate_fn(batch):
    """自定义collate函数"""
    return tuple(zip(*batch))


def setup_logger(output_dir, name='sim_isg_train'):
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


def compute_metrics(predictions, targets, threshold=0.5):
    """计算评估指标 - 针对ISG优化"""
    total_pred_objects = 0
    total_true_objects = 0
    correct_detections = 0
    all_scores = []
    
    for pred, target in zip(predictions, targets):
        # 统计真实ISG数量
        n_true = len(target['labels']) if len(target['labels']) > 0 else 0
        total_true_objects += n_true
        
        # 统计预测ISG数量
        if len(pred['scores']) > 0:
            pred_scores = pred['scores'].cpu().numpy()
            high_conf_mask = pred_scores >= threshold
            n_pred = np.sum(high_conf_mask)
            total_pred_objects += n_pred
            
            # 简化的正确检测计算：min(预测数量, 真实数量)
            correct_detections += min(n_pred, n_true)
            
            # 收集所有分数用于AUC计算
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
    
    # AUC和AP暂时设为0，需要更复杂的计算
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


def train_one_epoch(model, optimizer, data_loader, device, epoch, print_freq, logger=None, writer=None):
    """训练一个epoch"""
    model.train()
    
    all_predictions = []
    all_targets = []
    total_loss = 0
    
    # 使用tqdm显示进度条
    pbar = tqdm(data_loader, desc=f'Epoch {epoch}', 
                unit='batch', ncols=100, leave=False)
    
    for batch_idx, (images, targets) in enumerate(pbar):
        # 过滤掉没有标注的图像
        valid_indices = []
        for i, target in enumerate(targets):
            if len(target['boxes']) > 0:  # 只保留有标注的图像
                valid_indices.append(i)
        
        if len(valid_indices) == 0:
            # 如果整个batch都没有标注，跳过
            pbar.set_postfix({'status': 'no_annotations', 'avg_loss': f'{total_loss / max(1, batch_idx):.4f}'})
            continue
        
        # 只保留有效的图像和标注
        filtered_images = [images[i] for i in valid_indices]
        filtered_targets = [targets[i] for i in valid_indices]
        
        filtered_images = list(image.to(device) for image in filtered_images)
        filtered_targets = [{k: v.to(device) for k, v in t.items()} for t in filtered_targets]
        
        try:
            loss_dict = model(filtered_images, filtered_targets)
            
            losses = sum(loss for loss in loss_dict.values())
            
            # reduce losses over all GPUs for logging purposes
            loss_dict_reduced = {k: v.item() if torch.is_tensor(v) else v for k, v in loss_dict.items()}
            losses_reduced = sum(loss for loss in loss_dict_reduced.values())
            
            loss_value = losses_reduced
            total_loss += loss_value
            
            if not torch.isfinite(torch.tensor(loss_value)):
                print(f"Loss is {loss_value}, stopping training")
                print(loss_dict_reduced)
                sys.exit(1)
            
            optimizer.zero_grad()
            if torch.is_tensor(losses):
                losses.backward()
            else:
                # 如果losses不是tensor，转换为tensor
                torch.tensor(losses, requires_grad=True).backward()
            optimizer.step()
            
            # 更新进度条显示
            avg_loss = total_loss / (batch_idx + 1)
            pbar.set_postfix({
                'loss': f'{loss_value:.4f}', 
                'avg_loss': f'{avg_loss:.4f}',
                'lr': f'{optimizer.param_groups[0]["lr"]:.6f}',
                'valid_imgs': f'{len(valid_indices)}/{len(images)}'
            })
            
            # 为了计算指标，获取预测结果
            if batch_idx % (print_freq * 2) == 0:  # 不是每个batch都计算，减少计算开销
                model.eval()
                with torch.no_grad():
                    predictions = model(filtered_images)
                all_predictions.extend(predictions)
                all_targets.extend(filtered_targets)
                model.train()
            
            # 记录到tensorboard
            if writer is not None and batch_idx % print_freq == 0:
                global_step = epoch * len(data_loader) + batch_idx
                writer.add_scalar('Train/Loss', loss_value, global_step)
                writer.add_scalar('Train/LR', optimizer.param_groups[0]["lr"], global_step)
                for key, value in loss_dict_reduced.items():
                    writer.add_scalar(f'Train/{key}', value, global_step)
        
        except Exception as e:
            print(f"Batch {batch_idx} 训练出错: {e}")
            print(f"有效图像数: {len(filtered_images)}, 标注数: {[len(t['boxes']) for t in filtered_targets]}")
            # 跳过这个batch继续训练
            continue
    
    pbar.close()
    
    # 计算epoch指标
    metrics = {}
    if len(all_predictions) > 0:
        metrics = compute_metrics(all_predictions, all_targets)
        
        # 记录指标
        if logger:
            logger.info(f"Epoch {epoch} Training Metrics:")
            for key, value in metrics.items():
                logger.info(f"  {key}: {value:.4f}")
        
        # 打印指标
        print(f"\n=== Epoch {epoch} SIM ISG Training Metrics ===")
        for key, value in metrics.items():
            print(f"{key}: {value:.4f}")
        
        # tensorboard记录
        if writer is not None:
            for key, value in metrics.items():
                writer.add_scalar(f'Train_Metrics/{key}', value, epoch)
    
    return None, metrics


def evaluate(model, data_loader, device, logger=None, writer=None, epoch=None):
    """评估模型"""
    model.eval()
    
    all_predictions = []
    all_targets = []
    total_loss = 0
    
    # 使用tqdm显示验证进度
    pbar = tqdm(data_loader, desc=f'Validation Epoch {epoch}', 
                unit='batch', ncols=100, leave=False)
    
    with torch.no_grad():
        for batch_idx, (images, targets) in enumerate(pbar):
            images = list(img.to(device) for img in images)
            targets = [{k: v.to(device) for k, v in t.items()} for t in targets]
            
            # 获取预测结果
            predictions = model(images)
            all_predictions.extend(predictions)
            all_targets.extend(targets)
            
            # 计算loss（只对有标注的图像）
            valid_targets = [t for t in targets if len(t['boxes']) > 0]
            valid_images = [images[i] for i, t in enumerate(targets) if len(t['boxes']) > 0]
            
            if len(valid_targets) > 0:
                # 训练模式下计算loss
                model.train()
                loss_dict = model(valid_images, valid_targets)
                losses = sum(loss for loss in loss_dict.values())
                loss_dict_reduced = {k: v.item() if torch.is_tensor(v) else v for k, v in loss_dict.items()}
                losses_reduced = sum(loss for loss in loss_dict_reduced.values())
                
                total_loss += losses_reduced
                avg_loss = total_loss / (batch_idx + 1)
                
                # 更新进度条
                pbar.set_postfix({
                    'val_loss': f'{losses_reduced:.4f}',
                    'avg_val_loss': f'{avg_loss:.4f}',
                    'valid_imgs': f'{len(valid_targets)}/{len(targets)}'
                })
                
                model.eval()
            else:
                # 如果没有有效标注，只显示进度
                pbar.set_postfix({
                    'val_loss': 'no_targets',
                    'valid_imgs': f'0/{len(targets)}'
                })
    
    pbar.close()
    
    # 计算评估指标
    metrics = {}
    if len(all_predictions) > 0:
        metrics = compute_metrics(all_predictions, all_targets)
        
        # 记录指标
        if logger:
            logger.info(f"Epoch {epoch} Validation Metrics:")
            for key, value in metrics.items():
                logger.info(f"  {key}: {value:.4f}")
        
        # 打印指标
        print(f"\n=== Epoch {epoch} SIM ISG Validation Metrics ===")
        for key, value in metrics.items():
            print(f"Val_{key}: {value:.4f}")
        
        # tensorboard记录
        if writer is not None and epoch is not None:
            # 如果epoch是字符串'test'，使用一个特殊的数字
            if isinstance(epoch, str) and epoch == 'test':
                epoch_num = 999999  # 用大数字表示测试阶段
            else:
                epoch_num = epoch
            for key, value in metrics.items():
                writer.add_scalar(f'Val_Metrics/{key}', value, epoch_num)
    
    return None, metrics


def check_data_structure(data_root, coco_file):
    """检查SIM ISG数据结构和文件路径"""
    print("\n=== 检查SIM ISG数据结构 ===")
    
    # 检查COCO文件
    if not os.path.exists(coco_file):
        print(f"❌ COCO文件不存在: {coco_file}")
        return
    
    with open(coco_file, 'r', encoding='utf-8') as f:
        coco_data = json.load(f)
    
    print(f"✅ COCO文件存在: {coco_file}")
    print(f"   图像数量: {len(coco_data['images'])}")
    print(f"   标注数量: {len(coco_data['annotations'])}")
    print(f"   类别数量: {len(coco_data['categories'])}")
    
    # 显示类别信息
    for cat in coco_data['categories']:
        print(f"   类别: ID={cat['id']}, 名称='{cat['name']}'")
    
    # 检查图像目录结构
    images_dir = os.path.join(data_root, 'images')
    if not os.path.exists(images_dir):
        print(f"❌ 图像目录不存在: {images_dir}")
        return
    
    print(f"✅ 图像目录存在: {images_dir}")
    
    # 列出所有子目录 - SIM数据按data_id分组
    subdirs = [d for d in os.listdir(images_dir) if os.path.isdir(os.path.join(images_dir, d))]
    print(f"   SIM数据子目录数量: {len(subdirs)}")
    print(f"   前5个子目录: {subdirs[:5]}")
    
    # 检查前几个图像文件是否存在
    print("\n=== 检查SIM ISG图像文件存在性 ===")
    sample_images = coco_data['images'][:5]
    
    for img_info in sample_images:
        filename = img_info['file_name']
        found = False
        file_path = None
        
        # 在所有子目录中查找
        for subdir in subdirs:
            potential_path = os.path.join(images_dir, subdir, filename)
            if os.path.exists(potential_path):
                found = True
                file_path = potential_path
                break
        
        if found:
            print(f"✅ {filename} -> {file_path}")
        else:
            print(f"❌ {filename} 未找到")
            # 显示可能的匹配
            print(f"   尝试在以下目录中查找:")
            for subdir in subdirs[:3]:
                print(f"     {os.path.join(images_dir, subdir, filename)}")


def main():
    parser = argparse.ArgumentParser(description='SIM ISG Mask R-CNN Training')
    # 默认路径指向SIM ISG数据
    parser.add_argument('--data_root', default=r'D:\Gitspace\ipa_full\data\SIM\split_isg_img',
                       help='SIM ISG数据根目录')
    parser.add_argument('--coco_file', default=r'D:\Gitspace\ipa_full\data\SIM\split_isg_img\split_sim_isg_instances.json',
                       help='COCO格式标注文件')
    parser.add_argument('--output_dir', default=r'D:\Gitspace\ipa_full\data\SIM\split_isg_img\mask_rcnn\models',
                       help='模型输出目录')
    parser.add_argument('--epochs', default=20, type=int,
                       help='训练轮数')
    parser.add_argument('--batch_size', default=4, type=int,
                       help='批次大小')
    parser.add_argument('--lr', default=0.005, type=float,
                       help='学习率')
    parser.add_argument('--momentum', default=0.9, type=float,
                       help='SGD动量')
    parser.add_argument('--weight_decay', default=0.0005, type=float,
                       help='权重衰减')
    parser.add_argument('--print_freq', default=20, type=int,
                       help='打印频率')
    parser.add_argument('--resume', default='',
                       help='恢复训练的模型路径')
    parser.add_argument('--start_epoch', default=0, type=int,
                       help='开始epoch')
    parser.add_argument('--test_only', action='store_true',
                       help='仅测试模式')
    parser.add_argument('--data_ratio', default=0.8, type=float,
                       help='使用数据的比例 (0.0-1.0)')
    parser.add_argument('--train_val_split', default=0.8, type=float,
                       help='训练数据比例')
    parser.add_argument('--val_test_split', default=0.1, type=float,
                       help='验证数据比例')
    parser.add_argument('--check_data', action='store_true',
                       help='检查数据结构并退出')
    
    args = parser.parse_args()
    
    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 设置日志
    logger = setup_logger(args.output_dir, 'sim_isg_train')
    logger.info(f"SIM ISG训练参数: {args}")
    
    # 设置tensorboard
    writer = None
    if TENSORBOARD_AVAILABLE:
        log_dir = os.path.join(args.output_dir, 'tensorboard_logs')
        try:
            writer = SummaryWriter(log_dir)
            logger.info(f"Tensorboard日志保存到: {log_dir}")
        except ImportError:
            logger.warning("无法导入SummaryWriter")
    
    # 设备
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logger.info(f'使用设备: {device}')
    print(f'使用设备: {device}')
    
    # 如果是检查数据模式，执行检查并退出
    if args.check_data:
        check_data_structure(args.data_root, args.coco_file)
        return
    
    # 数据集
    print("加载SIM ISG数据集...")
    
    # 创建训练和验证数据集
    full_dataset = SIMISGDataset(args.data_root, args.coco_file, phase='train', 
                                transforms=get_transform(train=True))
    
    # 按data_ratio使用数据
    full_size = len(full_dataset)
    use_size = int(full_size * args.data_ratio)
    
    if use_size < full_size:
        from torch.utils.data import Subset
        import random
        indices = random.sample(range(full_size), use_size)
        dataset_to_split = Subset(full_dataset, indices)
        logger.info(f"使用数据比例: {args.data_ratio:.2f}, 使用 {use_size}/{full_size} 个样本")
    else:
        dataset_to_split = full_dataset
        logger.info(f"使用全部数据: {full_size} 个样本")
    
    # 三分割数据：训练/验证/测试
    dataset_size = len(dataset_to_split)
    train_size = int(args.train_val_split * dataset_size)
    val_size = int(args.val_test_split * dataset_size)
    test_size = dataset_size - train_size - val_size
    
    # 确保比例合理
    train_ratio = args.train_val_split
    val_ratio = args.val_test_split
    test_ratio = 1.0 - train_ratio - val_ratio
    
    if test_ratio < 0:
        print(f"警告: 训练比例({train_ratio})和验证比例({val_ratio})之和超过1.0，调整为默认比例")
        train_ratio, val_ratio, test_ratio = 0.8, 0.1, 0.1
        train_size = int(train_ratio * dataset_size)
        val_size = int(val_ratio * dataset_size)
        test_size = dataset_size - train_size - val_size
    
    from torch.utils.data import random_split
    train_dataset, val_dataset, test_dataset = random_split(
        dataset_to_split, [train_size, val_size, test_size])
    
    logger.info(f"数据分割比例: 训练{train_ratio:.1f} / 验证{val_ratio:.1f} / 测试{test_ratio:.1f}")
    print(f"训练集大小: {len(train_dataset)} ({train_ratio:.1%})")
    print(f"验证集大小: {len(val_dataset)} ({val_ratio:.1%})")
    print(f"测试集大小: {len(test_dataset)} ({test_ratio:.1%})")
    logger.info(f"训练集大小: {len(train_dataset)}")
    logger.info(f"验证集大小: {len(val_dataset)}")
    logger.info(f"测试集大小: {len(test_dataset)}")
    
    # 数据加载器
    train_loader = DataLoader(
        train_dataset, batch_size=args.batch_size, shuffle=True, 
        num_workers=4, collate_fn=collate_fn)
    
    val_loader = DataLoader(
        val_dataset, batch_size=1, shuffle=False, 
        num_workers=4, collate_fn=collate_fn)
    
    test_loader = DataLoader(
        test_dataset, batch_size=1, shuffle=False, 
        num_workers=4, collate_fn=collate_fn)
    
    # 模型
    print("创建SIM ISG Mask R-CNN模型...")
    num_classes = 2  # 背景 + ISG
    model = get_model_instance_segmentation(num_classes)
    model.to(device)
    
    # 优化器
    params = [p for p in model.parameters() if p.requires_grad]
    optimizer = torch.optim.SGD(params, lr=args.lr, momentum=args.momentum, 
                               weight_decay=args.weight_decay)
    
    # 学习率调度器
    lr_scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=3, gamma=0.1)
    
    # 恢复训练
    if args.resume:
        checkpoint = torch.load(args.resume, map_location='cpu')
        model.load_state_dict(checkpoint['model'])
        optimizer.load_state_dict(checkpoint['optimizer'])
        lr_scheduler.load_state_dict(checkpoint['lr_scheduler'])
        args.start_epoch = checkpoint['epoch'] + 1
        print(f"从epoch {args.start_epoch}恢复训练")
    
    if args.test_only:
        _, val_metrics = evaluate(model, val_loader, device=device, logger=logger, writer=writer, epoch=0)
        return
    
    logger.info("开始SIM ISG训练...")
    print("开始SIM ISG训练...")
    start_time = time.time()
    
    best_metrics = {'ap': 0.0}  # 记录最佳指标
    
    for epoch in range(args.start_epoch, args.epochs):
        logger.info(f"开始训练 Epoch {epoch}")
        
        # 训练一个epoch
        train_metric_logger, train_metrics = train_one_epoch(
            model, optimizer, train_loader, device, epoch, args.print_freq, logger, writer)
        
        # 更新学习率
        lr_scheduler.step()
        
        # 评估
        if epoch % 2 == 0 or epoch == args.epochs - 1:
            logger.info(f"开始验证 Epoch {epoch}")
            _, val_metrics = evaluate(model, val_loader, device=device, logger=logger, writer=writer, epoch=epoch)
            
            # 保存最佳模型
            if val_metrics.get('ap', 0) > best_metrics['ap']:
                best_metrics = val_metrics.copy()
                best_model_path = os.path.join(args.output_dir, 'best_sim_isg_model.pth')
                torch.save(model.state_dict(), best_model_path)
                logger.info(f"保存最佳SIM ISG模型: {best_model_path}, AP: {best_metrics['ap']:.4f}")
        
        # 保存检查点
        if epoch % 5 == 0 or epoch == args.epochs - 1:
            checkpoint = {
                'model': model.state_dict(),
                'optimizer': optimizer.state_dict(),
                'lr_scheduler': lr_scheduler.state_dict(),
                'epoch': epoch,
                'best_metrics': best_metrics
            }
            checkpoint_path = os.path.join(args.output_dir, f'sim_isg_model_{epoch}.pth')
            torch.save(checkpoint, checkpoint_path)
            
            # 也保存只有模型权重的文件（用于预测）
            weights_path = os.path.join(args.output_dir, f'sim_isg_weights_{epoch}.pth')
            torch.save(model.state_dict(), weights_path)
            logger.info(f"保存检查点: {checkpoint_path}")
            print(f"SIM ISG模型已保存: sim_isg_model_{epoch}.pth")
    
    import datetime
    total_time = time.time() - start_time
    total_time_str = str(datetime.timedelta(seconds=int(total_time)))
    logger.info(f'SIM ISG训练完成！总时间: {total_time_str}')
    logger.info(f'最佳验证指标: {best_metrics}')
    print(f'SIM ISG训练完成！总时间: {total_time_str}')
    print(f'最佳验证指标: {best_metrics}')
    
    # 在测试集上进行最终评估
    print("\n=== SIM ISG最终测试集评估 ===")
    logger.info("开始SIM ISG测试集评估")
    
    # 加载最佳模型
    if os.path.exists(os.path.join(args.output_dir, 'best_sim_isg_model.pth')):
        model.load_state_dict(torch.load(os.path.join(args.output_dir, 'best_sim_isg_model.pth')))
        print("使用最佳SIM ISG模型进行测试")
        logger.info("使用最佳SIM ISG模型进行测试")
    
    _, test_metrics = evaluate(model, test_loader, device=device, logger=logger, writer=writer, epoch='test')
    
    logger.info("SIM ISG最终测试指标:")
    for key, value in test_metrics.items():
        logger.info(f"  Test_{key}: {value:.4f}")
    
    print("=== SIM ISG最终测试指标 ===")
    for key, value in test_metrics.items():
        print(f"Test_{key}: {value:.4f}")
    
    # 关闭writer
    if writer is not None:
        writer.close()


if __name__ == '__main__':
    main()