#!/usr/bin/env python3
"""
简化的UPSNet训练脚本
使用虚拟数据和模型进行训练测试
"""

import os
import sys
import logging
import time
import torch
import torch.nn as nn
import torch.utils.data
import numpy as np
import cv2

# 设置随机种子
np.random.seed(235)
torch.manual_seed(235)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(235)

# 设置OpenCV
cv2.ocl.setUseOpenCL(False)

class SimpleConfig:
    """简化的配置类"""
    def __init__(self):
        # 基础配置
        self.output_path = "./output"
        self.model_prefix = "upsnet_"
        self.gpus = "0"
        
        # 训练配置
        self.batch_size = 2
        self.lr = 0.001
        self.momentum = 0.9
        self.weight_decay = 0.0001
        self.max_iteration = 1000  # 简化训练轮数
        self.display_iter = 10     # 每10轮显示一次
        self.snapshot_step = 100   # 每100轮保存一次
        
        # 网络配置
        self.has_rpn = True
        self.has_rcnn = True
        self.has_mask_head = True

class SimpleDataset(torch.utils.data.Dataset):
    """简化的数据集类"""
    def __init__(self, length=100):
        self.length = length
        
    def __len__(self):
        return self.length
    
    def __getitem__(self, idx):
        # 生成虚拟数据
        image = torch.randn(3, 256, 256)  # 简化图像尺寸
        
        # 虚拟标签
        data = {'image': image}
        label = {'gt_boxes': torch.tensor([[50, 50, 150, 150, 1]])}
        
        return data, label, idx

class SimpleUPSNet(nn.Module):
    """简化的UPSNet模型"""
    def __init__(self):
        super(SimpleUPSNet, self).__init__()
        print("创建简化的UPSNet模型...")
        
        # 简单的特征提取器
        self.backbone = nn.Sequential(
            nn.Conv2d(3, 64, 3, 1, 1),
            nn.ReLU(),
            nn.Conv2d(64, 128, 3, 2, 1),
            nn.ReLU(),
            nn.AdaptiveAvgPool2d((1, 1)),
            nn.Flatten(),
            nn.Linear(128, 256)
        )
        
        # 简单的分类头
        self.classifier = nn.Linear(256, 2)
        
    def forward(self, data, label=None):
        """前向传播"""
        if isinstance(data, dict):
            x = data['image']
        else:
            x = data
            
        # 特征提取
        features = self.backbone(x)
        
        # 分类
        cls_output = self.classifier(features)
        
        # 返回虚拟损失（模拟UPSNet的输出格式）
        device = x.device
        output = {
            'rpn_cls_loss': torch.tensor(0.5, requires_grad=True, device=device),
            'rpn_bbox_loss': torch.tensor(0.3, requires_grad=True, device=device),
            'cls_loss': torch.tensor(0.4, requires_grad=True, device=device),
            'bbox_loss': torch.tensor(0.2, requires_grad=True, device=device),
            'mask_loss': torch.tensor(0.3, requires_grad=True, device=device),
            'rcnn_accuracy': torch.tensor(0.8, device=device),
        }
        
        return output

def create_simple_logger(output_path):
    """创建简单的日志记录器"""
    os.makedirs(output_path, exist_ok=True)
    
    logger = logging.getLogger('SimpleUPSNet')
    logger.setLevel(logging.INFO)
    
    # 控制台输出
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # 文件输出
    log_file = os.path.join(output_path, 'train.log')
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    return logger

def simple_train():
    """简化的训练函数"""
    
    # 配置
    config = SimpleConfig()
    
    # 创建输出目录和日志
    logger = create_simple_logger(config.output_path)
    logger.info("开始简化训练...")
    
    # 设备设置
    device = torch.device(f'cuda:{config.gpus}' if torch.cuda.is_available() else 'cpu')
    logger.info(f"使用设备: {device}")
    
    # 创建模型
    model = SimpleUPSNet().to(device)
    
    # 创建优化器
    optimizer = torch.optim.SGD(
        model.parameters(), 
        lr=config.lr, 
        momentum=config.momentum, 
        weight_decay=config.weight_decay
    )
    
    # 创建数据集和数据加载器
    train_dataset = SimpleDataset(length=200)
    train_loader = torch.utils.data.DataLoader(
        train_dataset, 
        batch_size=config.batch_size, 
        shuffle=True,
        num_workers=0  # 简化为0避免多进程问题
    )
    
    # 训练循环
    model.train()
    total_loss_sum = 0.0
    
    for iteration in range(config.max_iteration):
        for batch_idx, (data, label, _) in enumerate(train_loader):
            if iteration >= config.max_iteration:
                break
                
            # 数据移到设备
            for k, v in data.items():
                if torch.is_tensor(v):
                    data[k] = v.to(device)
            
            # 前向传播
            optimizer.zero_grad()
            output = model(data, label)
            
            # 计算总损失
            total_loss = torch.tensor(0.0, device=device, requires_grad=True)
            loss_components = []
            
            if config.has_rpn:
                total_loss = total_loss + output['rpn_cls_loss'] + output['rpn_bbox_loss']
                loss_components.extend(['rpn_cls_loss', 'rpn_bbox_loss'])
                
            if config.has_rcnn:
                total_loss = total_loss + output['cls_loss'] + output['bbox_loss']
                loss_components.extend(['cls_loss', 'bbox_loss'])
                
            if config.has_mask_head:
                total_loss = total_loss + output['mask_loss']
                loss_components.append('mask_loss')
            
            # 反向传播
            total_loss.backward()
            optimizer.step()
            
            total_loss_sum += total_loss.item()
            
            # 显示进度
            if iteration % config.display_iter == 0:
                avg_loss = total_loss_sum / max(iteration + 1, 1)
                logger.info(f"Iteration [{iteration:4d}/{config.max_iteration}] "
                          f"Loss: {total_loss.item():.4f} "
                          f"Avg Loss: {avg_loss:.4f}")
            
            # 保存模型
            if iteration % config.snapshot_step == 0 and iteration > 0:
                model_path = os.path.join(config.output_path, f'{config.model_prefix}{iteration}.pth')
                torch.save(model.state_dict(), model_path)
                logger.info(f"模型已保存: {model_path}")
            
            iteration += 1
            if iteration >= config.max_iteration:
                break
    
    # 最终保存
    final_model_path = os.path.join(config.output_path, f'{config.model_prefix}final.pth')
    torch.save(model.state_dict(), final_model_path)
    logger.info(f"训练完成! 最终模型保存: {final_model_path}")
    
    return model

if __name__ == '__main__':
    print("=" * 50)
    print("简化的UPSNet训练脚本")
    print("=" * 50)
    
    try:
        model = simple_train()
        print("训练成功完成!")
    except Exception as e:
        print(f"训练出错: {e}")
        import traceback
        traceback.print_exc()