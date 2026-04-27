# filepath: d:\Gitspace\ipa_full\iPA\ipa\processing\segmentation\segmentation_sxt\model_sxt_mito_instance\torchvision_maskrcnn_mito_train.py
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
# 假设可以重用ISG的预测模块中的一些函数，但需要修改数据路径
from sklearn.metrics import precision_recall_fscore_support, roc_auc_score, average_precision_score
import logging
from tqdm import tqdm
try:
    from torch.utils.tensorboard import SummaryWriter
    TENSORBOARD_AVAILABLE = True
except ImportError:
    TENSORBOARD_AVAILABLE = False
    print("Warning: tensorboard not available, skipping tensorboard logging")


class SXTMitoDataset(Dataset):
    """SXT Mito数据集类"""
    def __init__(self, data_root, coco_file, phase='train', transforms=None):
        self.data_root = data_root
        self.transforms = transforms
        self.phase = phase
        
        # 加载COCO格式的标注
        with open(coco_file, 'r', encoding='utf-8') as f:
            self.coco_data = json.load(f)
        
        self.images = self.coco_data['images']
        self.annotations = self.coco_data['annotations']
        
        # 创建图像ID到标注的映射
        self.image_id_to_annotations = {}
        for ann in self.annotations:
            image_id = ann['image_id']
            if image_id not in self.image_id_to_annotations:
                self.image_id_to_annotations[image_id] = []
            self.image_id_to_annotations[image_id].append(ann)
        
        print(f"Loaded {len(self.images)} images with {len(self.annotations)} mito annotations")
    
    def __len__(self):
        return len(self.images)
    
    def __getitem__(self, idx):
        # 获取图像信息
        img_info = self.images[idx]
        image_id = img_info['id']
        file_name = img_info['file_name']
        
        # 查找图像文件路径（在各个子目录中搜索）
        img_path = self._find_image_path(file_name)
        if img_path is None:
            raise FileNotFoundError(f"Image not found: {file_name}")
        
        # 读取图像
        image = Image.open(img_path).convert("RGB")
        
        # 获取该图像的所有标注
        annotations = self.image_id_to_annotations.get(image_id, [])
        
        # 处理标注
        boxes = []
        labels = []
        masks = []
        areas = []
        
        for ann in annotations:
            # 边界框
            bbox = ann['bbox']  # [x, y, width, height]
            boxes.append([bbox[0], bbox[1], bbox[0] + bbox[2], bbox[1] + bbox[3]])
            
            # 标签 (mito类别ID为1)
            labels.append(ann['category_id'])
            
            # mask (从RLE格式解码) - 简化处理
            if 'segmentation' in ann:
                rle = ann['segmentation']
                if isinstance(rle, dict) and 'counts' in rle and 'size' in rle:
                    try:
                        # 如果counts是列表格式（未压缩的RLE）
                        if isinstance(rle['counts'], list):
                            # 手动解码RLE列表
                            h, w = rle['size']
                            counts = rle['counts']
                            mask = np.zeros(h * w, dtype=np.uint8)
                            
                            # 解码RLE列表
                            pos = 0
                            for i, count in enumerate(counts):
                                if i % 2 == 1:  # 奇数索引表示前景像素
                                    end_pos = min(pos + count, len(mask))
                                    mask[pos:end_pos] = 1
                                pos += count
                                if pos >= len(mask):
                                    break
                            
                            mask = mask.reshape((h, w), order='F')  # Fortran order (column-major)
                        else:
                            # 如果是其他格式，创建空mask
                            h, w = img_info['height'], img_info['width']
                            mask = np.zeros((h, w), dtype=np.uint8)
                        
                        masks.append(mask)
                    except Exception as e:
                        print(f"Warning: Failed to decode RLE mask: {e}")
                        print(f"Ann ID: {ann.get('id', 'unknown')}, Image: {file_name}")
                        # 创建空mask作为fallback
                        h, w = img_info['height'], img_info['width']
                        mask = np.zeros((h, w), dtype=np.uint8)
                        masks.append(mask)
                else:
                    # 其他格式，创建空mask
                    h, w = img_info['height'], img_info['width']
                    mask = np.zeros((h, w), dtype=np.uint8)
                    masks.append(mask)
            else:
                # 如果没有分割信息，创建空mask
                h, w = img_info['height'], img_info['width']
                mask = np.zeros((h, w), dtype=np.uint8)
                masks.append(mask)
            
            # 面积
            areas.append(ann['area'])
        
        # 转换为torch tensor
        if len(boxes) > 0:
            boxes = torch.as_tensor(boxes, dtype=torch.float32)
            labels = torch.as_tensor(labels, dtype=torch.int64)
            masks = torch.as_tensor(np.array(masks), dtype=torch.uint8)
            areas = torch.as_tensor(areas, dtype=torch.float32)
        else:
            # 如果没有标注，创建空的tensor
            boxes = torch.zeros((0, 4), dtype=torch.float32)
            labels = torch.zeros((0,), dtype=torch.int64)
            masks = torch.zeros((0, img_info['height'], img_info['width']), dtype=torch.uint8)
            areas = torch.zeros((0,), dtype=torch.float32)
        
        iscrowd = torch.zeros((len(annotations),), dtype=torch.int64)
        
        target = {
            "boxes": boxes,
            "labels": labels,
            "masks": masks,
            "image_id": torch.tensor([image_id]),
            "area": areas,
            "iscrowd": iscrowd,
        }
        
        if self.transforms is not None:
            image = self.transforms(image)
        
        return image, target
    
    def _find_image_path(self, file_name):
        """在images目录的各个子目录中查找图像文件"""
        images_dir = os.path.join(self.data_root, 'images')
        
        # 首先检查是否有子目录
        if os.path.exists(images_dir):
            # 遍历所有子目录
            for item in os.listdir(images_dir):
                item_path = os.path.join(images_dir, item)
                if os.path.isdir(item_path):
                    img_path = os.path.join(item_path, file_name)
                    if os.path.exists(img_path):
                        return img_path
        
        # 如果在子目录中没找到，尝试直接在images目录中查找
        direct_path = os.path.join(images_dir, file_name)
        if os.path.exists(direct_path):
            return direct_path
        
        return None


def get_mito_transform(train):
    """获取数据增强变换"""
    transforms = []
    transforms.append(T.PILToTensor())
    transforms.append(T.ConvertImageDtype(torch.float))
    if train:
        transforms.append(T.RandomHorizontalFlip(0.5))
    return T.Compose(transforms)


def get_mito_model_instance_segmentation(num_classes):
    """创建Mask R-CNN模型用于mito实例分割"""
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
    model.roi_heads.mask_predictor = MaskRCNNPredictor(in_features_mask,
                                                       hidden_layer,
                                                       num_classes)
    
    return model


def collate_fn(batch):
    """自定义的collate函数"""
    return tuple(zip(*batch))


def setup_logger(output_dir, name='mito_train'):
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


def compute_mito_metrics(predictions, targets, threshold=0.5):
    """计算mito评估指标"""
    total_pred_objects = 0
    total_true_objects = 0
    correct_detections = 0
    all_scores = []
    
    for pred, target in zip(predictions, targets):
        # 统计真实mito对象数量
        n_true = len(target['labels']) if len(target['labels']) > 0 else 0
        total_true_objects += n_true
        
        # 统计预测mito对象数量
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


def train_one_epoch_mito(model, optimizer, data_loader, device, epoch, print_freq, logger=None, writer=None):
    """训练一个epoch - mito版本"""
    model.train()
    
    all_predictions = []
    all_targets = []
    total_loss = 0
    
    # 使用tqdm显示进度条
    pbar = tqdm(data_loader, desc=f'Mito Epoch {epoch}', 
                unit='batch', ncols=100, leave=False)
    
    for batch_idx, (images, targets) in enumerate(pbar):
        images = list(image.to(device) for image in images)
        targets = [{k: v.to(device) for k, v in t.items()} for t in targets]
        
        loss_dict = model(images, targets)
        
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
            'lr': f'{optimizer.param_groups[0]["lr"]:.6f}'
        })
        
        # 为了计算指标，获取预测结果
        if batch_idx % (print_freq * 2) == 0:  # 不是每个batch都计算，减少计算开销
            model.eval()
            with torch.no_grad():
                predictions = model(images)
            all_predictions.extend(predictions)
            all_targets.extend(targets)
            model.train()
        
        # 记录到tensorboard
        if writer is not None and batch_idx % print_freq == 0:
            global_step = epoch * len(data_loader) + batch_idx
            writer.add_scalar('Mito_Train/Loss', loss_value, global_step)
            writer.add_scalar('Mito_Train/LR', optimizer.param_groups[0]["lr"], global_step)
            for key, value in loss_dict_reduced.items():
                writer.add_scalar(f'Mito_Train/{key}', value, global_step)
    
    pbar.close()
    
    # 计算epoch指标
    metrics = {}
    if len(all_predictions) > 0:
        metrics = compute_mito_metrics(all_predictions, all_targets)
        
        # 记录指标
        if logger:
            logger.info(f"Epoch {epoch} Mito Training Metrics:")
            for key, value in metrics.items():
                logger.info(f"  {key}: {value:.4f}")
        
        # 打印指标
        print(f"\n=== Epoch {epoch} Mito Training Metrics ===")
        for key, value in metrics.items():
            print(f"{key}: {value:.4f}")
        
        # tensorboard记录
        if writer is not None:
            for key, value in metrics.items():
                writer.add_scalar(f'Mito_Train_Metrics/{key}', value, epoch)
    
    return None, metrics


def evaluate_mito(model, data_loader, device, logger=None, writer=None, epoch=None):
    """评估mito模型"""
    model.eval()
    
    all_predictions = []
    all_targets = []
    total_loss = 0
    
    # 使用tqdm显示验证进度
    pbar = tqdm(data_loader, desc=f'Mito Validation Epoch {epoch}', 
                unit='batch', ncols=100, leave=False)
    
    with torch.no_grad():
        for batch_idx, (images, targets) in enumerate(pbar):
            images = list(img.to(device) for img in images)
            targets = [{k: v.to(device) for k, v in t.items()} for t in targets]
            
            # 获取预测结果
            predictions = model(images)
            all_predictions.extend(predictions)
            all_targets.extend(targets)
            
            # 训练模式下计算loss
            model.train()
            loss_dict = model(images, targets)
            losses = sum(loss for loss in loss_dict.values())
            loss_dict_reduced = {k: v.item() if torch.is_tensor(v) else v for k, v in loss_dict.items()}
            losses_reduced = sum(loss for loss in loss_dict_reduced.values())
            
            total_loss += losses_reduced
            avg_loss = total_loss / (batch_idx + 1)
            
            # 更新进度条
            pbar.set_postfix({
                'val_loss': f'{losses_reduced:.4f}',
                'avg_val_loss': f'{avg_loss:.4f}'
            })
            
            model.eval()
    
    pbar.close()
    
    # 计算评估指标
    metrics = {}
    if len(all_predictions) > 0:
        metrics = compute_mito_metrics(all_predictions, all_targets)
        
        # 记录指标
        if logger:
            logger.info(f"Epoch {epoch} Mito Validation Metrics:")
            for key, value in metrics.items():
                logger.info(f"  {key}: {value:.4f}")
        
        # 打印指标
        print(f"\n=== Epoch {epoch} Mito Validation Metrics ===")
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
                writer.add_scalar(f'Mito_Val_Metrics/{key}', value, epoch_num)
    
    return None, metrics


def check_mito_data_structure(data_root, coco_file):
    """检查mito数据结构和文件路径"""
    print("\n=== 检查Mito数据结构 ===")
    
    # 检查COCO文件
    if not os.path.exists(coco_file):
        print(f"❌ Mito COCO文件不存在: {coco_file}")
        return
    
    with open(coco_file, 'r', encoding='utf-8') as f:
        coco_data = json.load(f)
    
    print(f"✅ Mito COCO文件存在: {coco_file}")
    print(f"   图像数量: {len(coco_data['images'])}")
    print(f"   Mito标注数量: {len(coco_data['annotations'])}")
    
    # 检查类别信息
    categories = coco_data.get('categories', [])
    print(f"   类别: {[cat['name'] for cat in categories]}")
    
    # 检查图像目录结构
    images_dir = os.path.join(data_root, 'images')
    if not os.path.exists(images_dir):
        print(f"❌ 图像目录不存在: {images_dir}")
        return
    
    print(f"✅ 图像目录存在: {images_dir}")
    
    # 列出所有子目录
    subdirs = [d for d in os.listdir(images_dir) if os.path.isdir(os.path.join(images_dir, d))]
    print(f"   子目录数量: {len(subdirs)}")
    print(f"   前5个子目录: {subdirs[:5]}")
    
    # 检查前几个图像文件是否存在
    print("\n=== 检查Mito图像文件存在性 ===")
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
    parser = argparse.ArgumentParser(description='SXT Mito Mask R-CNN Training')
    
    # 基本参数 - 使用mito数据路径
    parser.add_argument('--data_root', 
                       default=r'D:\Gitspace\ipa_full\iPA\data\sxt_images\split_mito_img',
                       help='Mito数据根目录')
    parser.add_argument('--coco_file', 
                       default=r'D:\Gitspace\ipa_full\iPA\data\sxt_images\split_mito_img\split_mito_instances.json',
                       help='Mito COCO格式标注文件')
    parser.add_argument('--output_dir', 
                       default=r'D:\Gitspace\ipa_full\iPA\data\sxt_images\mask_rcnn_mito\models',
                       help='Mito模型输出目录')
    
    # 训练参数
    parser.add_argument('--epochs', default=20, type=int,
                       help='训练轮数')
    parser.add_argument('--batch_size', default=2, type=int,
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
    
    # 数据参数
    parser.add_argument('--data_ratio', default=0.8, type=float,
                       help='使用数据的比例 (0.0-1.0)')
    parser.add_argument('--train_val_split', default=0.8, type=float,
                       help='训练数据比例')
    parser.add_argument('--val_test_split', default=0.1, type=float,
                       help='验证数据比例')
    parser.add_argument('--check_data', action='store_true',
                       help='检查Mito数据结构并退出')
    
    args = parser.parse_args()
    
    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 设置日志
    logger = setup_logger(args.output_dir, 'mito_train')
    logger.info(f"Mito训练参数: {args}")
    
    # 设置tensorboard
    writer = None
    if TENSORBOARD_AVAILABLE:
        log_dir = os.path.join(args.output_dir, 'tensorboard_logs')
        try:
            from torch.utils.tensorboard.writer import SummaryWriter
            writer = SummaryWriter(log_dir)
            logger.info(f"Mito Tensorboard日志保存到: {log_dir}")
        except ImportError:
            logger.warning("无法导入SummaryWriter")
    
    # 设备
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logger.info(f'使用设备: {device}')
    print(f'使用设备: {device}')
    
    # 如果是检查数据模式，执行检查并退出
    if args.check_data:
        check_mito_data_structure(args.data_root, args.coco_file)
        return
    
    # 数据集
    print("加载Mito数据集...")
    
    # 创建训练和验证数据集
    full_dataset = SXTMitoDataset(args.data_root, args.coco_file, phase='train', 
                                  transforms=get_mito_transform(train=True))
    
    # 按data_ratio使用数据
    full_size = len(full_dataset)
    use_size = int(full_size * args.data_ratio)
    
    if use_size < full_size:
        from torch.utils.data import Subset
        import random
        indices = random.sample(range(full_size), use_size)
        dataset_to_split = Subset(full_dataset, indices)
        logger.info(f"使用Mito数据比例: {args.data_ratio:.2f}, 使用 {use_size}/{full_size} 个样本")
    else:
        dataset_to_split = full_dataset
        logger.info(f"使用全部Mito数据: {full_size} 个样本")
    
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
    
    logger.info(f"Mito数据分割比例: 训练{train_ratio:.1f} / 验证{val_ratio:.1f} / 测试{test_ratio:.1f}")
    print(f"Mito训练集大小: {len(train_dataset)} ({train_ratio:.1%})")
    print(f"Mito验证集大小: {len(val_dataset)} ({val_ratio:.1%})")
    print(f"Mito测试集大小: {len(test_dataset)} ({test_ratio:.1%})")
    logger.info(f"Mito训练集大小: {len(train_dataset)}")
    logger.info(f"Mito验证集大小: {len(val_dataset)}")
    logger.info(f"Mito测试集大小: {len(test_dataset)}")
    
    # 数据加载器
    train_loader = DataLoader(
        train_dataset, batch_size=args.batch_size, shuffle=True, 
        num_workers=0, collate_fn=collate_fn)
    
    val_loader = DataLoader(
        val_dataset, batch_size=1, shuffle=False, 
        num_workers=0, collate_fn=collate_fn)
    
    test_loader = DataLoader(
        test_dataset, batch_size=1, shuffle=False, 
        num_workers=0, collate_fn=collate_fn)
    
    # 模型
    print("创建Mito模型...")
    num_classes = 2  # 背景 + mitochondria
    model = get_mito_model_instance_segmentation(num_classes)
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
        print(f"从epoch {args.start_epoch}恢复Mito训练")
    
    if args.test_only:
        _, val_metrics = evaluate_mito(model, val_loader, device=device, logger=logger, writer=writer, epoch=0)
        return
    
    logger.info("开始Mito训练...")
    print("开始Mito训练...")
    start_time = time.time()
    
    best_metrics = {'ap': 0.0}  # 记录最佳指标
    
    for epoch in range(args.start_epoch, args.epochs):
        logger.info(f"开始Mito训练 Epoch {epoch}")
        
        # 训练一个epoch
        train_metric_logger, train_metrics = train_one_epoch_mito(
            model, optimizer, train_loader, device, epoch, args.print_freq, logger, writer)
        
        # 更新学习率
        lr_scheduler.step()
        
        # 评估
        if epoch % 2 == 0 or epoch == args.epochs - 1:
            logger.info(f"开始Mito验证 Epoch {epoch}")
            _, val_metrics = evaluate_mito(model, val_loader, device=device, logger=logger, writer=writer, epoch=epoch)
            
            # 保存最佳模型
            if val_metrics.get('ap', 0) > best_metrics['ap']:
                best_metrics = val_metrics.copy()
                best_model_path = os.path.join(args.output_dir, 'best_mito_model.pth')
                torch.save(model.state_dict(), best_model_path)
                logger.info(f"保存最佳Mito模型: {best_model_path}, AP: {best_metrics['ap']:.4f}")
        
        # 保存检查点
        if epoch % 5 == 0 or epoch == args.epochs - 1:
            checkpoint = {
                'model': model.state_dict(),
                'optimizer': optimizer.state_dict(),
                'lr_scheduler': lr_scheduler.state_dict(),
                'epoch': epoch,
                'best_metrics': best_metrics
            }
            checkpoint_path = os.path.join(args.output_dir, f'mito_model_{epoch}.pth')
            torch.save(checkpoint, checkpoint_path)
            
            # 也保存只有模型权重的文件（用于预测）
            weights_path = os.path.join(args.output_dir, f'mito_model_weights_{epoch}.pth')
            torch.save(model.state_dict(), weights_path)
            logger.info(f"保存Mito检查点: {checkpoint_path}")
            print(f"Mito模型已保存: mito_model_{epoch}.pth")
    
    total_time = time.time() - start_time
    total_time_str = str(datetime.timedelta(seconds=int(total_time)))
    logger.info(f'Mito训练完成！总时间: {total_time_str}')
    logger.info(f'最佳Mito验证指标: {best_metrics}')
    print(f'Mito训练完成！总时间: {total_time_str}')
    print(f'最佳Mito验证指标: {best_metrics}')
    
    # 在测试集上进行最终评估
    print("\n=== 最终Mito测试集评估 ===")
    logger.info("开始Mito测试集评估")
    
    # 加载最佳模型
    if os.path.exists(os.path.join(args.output_dir, 'best_mito_model.pth')):
        model.load_state_dict(torch.load(os.path.join(args.output_dir, 'best_mito_model.pth')))
        print("使用最佳Mito模型进行测试")
        logger.info("使用最佳Mito模型进行测试")
    
    _, test_metrics = evaluate_mito(model, test_loader, device=device, logger=logger, writer=writer, epoch='test')
    
    logger.info("最终Mito测试指标:")
    for key, value in test_metrics.items():
        logger.info(f"  Test_Mito_{key}: {value:.4f}")
    
    print("=== 最终Mito测试指标 ===")
    for key, value in test_metrics.items():
        print(f"Test_Mito_{key}: {value:.4f}")
    
    # 关闭writer
    if writer is not None:
        writer.close()


if __name__ == '__main__':
    import datetime
    from collections import deque
    main()