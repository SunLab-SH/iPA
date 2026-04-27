# ---------------------------------------------------------------------------
# Unified Panoptic Segmentation Network
#
# Copyright (c) 2018-2019 Uber Technologies, Inc.
#
# Licensed under the Uber Non-Commercial License (the "License");
# you may no# ----if is_master:
#     logger, final_output_path = create_logger(config.output_path, args.cfg, config.dataset.image_set)
#     # 使用简单的指标记录器替代TensorboardX
#     writer = SimpleMetricLogger(log_dir=os.path.join(config.output_path, 'metrics',
#                                                      os.path.basename(args.cfg).split('.')[0],
#                                                      '_'.join(config.dataset.image_set.split('+')),
#                                                      config.model_prefix+time.strftime('%Y-%m-%d-%H-%M')))
# else:
#     final_output_path = os.path.join(config.output_path, os.path.basename(args.cfg).split('.')[0], '{}'.format('_'.join([iset for iset in config.dataset.image_set.split('+')])))
#     writer = None---------------------------------------------------------------
# Unified Panoptic Segmentation Network
#
# Copyright (c) 2018-2019 Uber Technologies, Inc.
#
# Licensed under the Uber Non-Commercial License (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at the root directory of this project. 
#
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Written by Yuwen Xiong
# ---------------------------------------------------------------------------

from __future__ import print_function, division
import os
import sys
import logging
import pprint
import time
import numpy as np
import pickle
import torch
import torch.nn as nn
import torch.utils.data
import torch.backends.cudnn as cudnn
import cv2
import torch.utils.data.distributed as distributed


sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from collections import deque
from .upsnet.config.config import config
from .upsnet.config.parse_args import parse_args
from .utils.metrics import AvgMetric

args = parse_args()

# Update config from args
config.gpus = args.gpus if hasattr(args, 'gpus') else "0"
config.train.batch_size = args.batch_size if hasattr(args, 'batch_size') else 2
config.train.lr = args.lr if hasattr(args, 'lr') else 0.001
config.train.max_iteration = args.max_iter if hasattr(args, 'max_iter') else 10000
config.output_path = args.output if hasattr(args, 'output') else "./output"

# Single GPU training only
config.train.use_horovod = False
is_master = True

def create_logger(output_path, cfg_name, image_set):
    """创建日志记录器"""
    os.makedirs(output_path, exist_ok=True)
    
    # 创建日志记录器
    logger = logging.getLogger('UPSNet')
    logger.setLevel(logging.INFO)
    
    # 创建文件处理器
    log_file = os.path.join(output_path, 'train.log')
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    
    # 创建控制台处理器
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    
    # 创建格式器
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)
    
    # 添加处理器到记录器
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    final_output_path = os.path.join(output_path, 
                                   os.path.basename(cfg_name).split('.')[0], 
                                   '_'.join(image_set.split('+')))
    os.makedirs(final_output_path, exist_ok=True)
    
    return logger, final_output_path

# 简单的指标记录器，替代TensorboardX
class SimpleMetricLogger:
    """简单的指标记录器，将指标写入文件"""
    def __init__(self, log_dir):
        self.log_dir = log_dir
        os.makedirs(log_dir, exist_ok=True)
        self.metrics_file = os.path.join(log_dir, 'metrics.log')
        
        # 创建指标文件并写入头部
        with open(self.metrics_file, 'w') as f:
            f.write('iteration,metric_name,value\n')
    
    def add_scalar(self, metric_name, value, iteration):
        """记录标量指标"""
        with open(self.metrics_file, 'a') as f:
            f.write(f'{iteration},{metric_name},{value:.6f}\n')
    
    def close(self):
        """关闭记录器"""
        pass

# 简单的速度计
class Speedometer:
    def __init__(self, batch_size, display_freq):
        self.batch_size = batch_size
        self.display_freq = display_freq
        self.last_time = time.time()
        
    def __call__(self, curr_iter, metrics):
        if curr_iter % self.display_freq == 0:
            current_time = time.time()
            speed = self.batch_size * self.display_freq / (current_time - self.last_time)
            self.last_time = current_time
            
            log_str = f'Iter [{curr_iter}], Speed: {speed:.2f} samples/sec'
            if metrics:
                for metric in metrics:
                    name, value = metric.get()
                    log_str += f', {name}: {value:.4f}'
            print(log_str)

# 简单的SGD优化器包装器
class SGD:
    def __init__(self, params_lr, lr, momentum, weight_decay):
        self.param_groups = []
        for param_group in params_lr:
            group = {
                'params': param_group['params'],
                'lr': param_group['lr'] * lr,
                'momentum': momentum,
                'weight_decay': weight_decay
            }
            self.param_groups.append(group)
        
        # 使用第一个组的学习率作为默认lr
        default_lr = self.param_groups[0]['lr'] if self.param_groups else lr
        self.optimizer = torch.optim.SGD(
            [{'params': group['params'], 'lr': group['lr']} for group in self.param_groups],
            lr=default_lr,  # 添加默认lr参数
            momentum=momentum,
            weight_decay=weight_decay
        )
    
    def zero_grad(self):
        self.optimizer.zero_grad()
    
    def step(self, lr_scale):
        # 调整学习率
        for i, group in enumerate(self.optimizer.param_groups):
            group['lr'] = self.param_groups[i]['lr'] * lr_scale
        self.optimizer.step()
    
    def state_dict(self):
        return self.optimizer.state_dict()
    
    def load_state_dict(self, state_dict):
        self.optimizer.load_state_dict(state_dict)

# 虚拟数据集类（需要替换为实际实现）
class CocoDataset:
    def __init__(self, image_sets, flip=False, result_path=None, phase='train'):
        self.image_sets = image_sets
        self.flip = flip
        self.result_path = result_path
        self.phase = phase
        print(f"警告: CocoDataset 是虚拟实现，需要根据实际数据替换")
        # 这里需要实际的数据加载逻辑
        
    def __len__(self):
        return 100  # 虚拟长度
    
    def __getitem__(self, idx):
        # 返回虚拟数据
        data = {'image': torch.randn(3, 512, 512)}
        label = {'gt_boxes': torch.tensor([[100, 100, 200, 200, 1]])}
        return data, label, idx
    
    def collate(self, batch):
        # 简单的批处理函数
        return batch[0]  # 暂时返回第一个样本

# 虚拟模型类（需要替换为实际实现）
class resnet50_upsnet(nn.Module):
    def __init__(self):
        super(resnet50_upsnet, self).__init__()
        print(f"警告: resnet50_upsnet 是虚拟实现，需要加载实际模型")
        
        # 简单的虚拟网络
        self.backbone = nn.Sequential(
            nn.Conv2d(3, 64, 7, 2, 3),
            nn.ReLU(),
            nn.AdaptiveAvgPool2d((1, 1)),
            nn.Flatten(),
            nn.Linear(64, 512)
        )
        
        # 添加resnet_backbone属性用于冻结
        self.resnet_backbone = nn.Sequential(
            nn.Conv2d(3, 64, 7, 2, 3),
            nn.BatchNorm2d(64),
            nn.ReLU(inplace=True),
            nn.MaxPool2d(kernel_size=3, stride=2, padding=1)
        )
        
    def forward(self, data, label=None):
        # 虚拟前向传播
        if isinstance(data, dict):
            x = data.get('image', torch.randn(1, 3, 512, 512))
        else:
            x = data[0]['image'] if isinstance(data, (list, tuple)) else data
            
        features = self.backbone(x)
        
        # 确保损失在正确的设备上
        device = x.device if torch.is_tensor(x) else torch.device('cuda')
        
        # 返回虚拟损失
        output = {
            'rpn_cls_loss': torch.tensor(0.5, requires_grad=True, device=device),
            'rpn_bbox_loss': torch.tensor(0.3, requires_grad=True, device=device),
            'rcnn_accuracy': torch.tensor(0.8, device=device),
            'cls_loss': torch.tensor(0.4, requires_grad=True, device=device),
            'bbox_loss': torch.tensor(0.2, requires_grad=True, device=device),
            'mask_loss': torch.tensor(0.3, requires_grad=True, device=device),
        }
        return output
    
    def freeze_backbone(self, freeze_at):
        """冻结骨干网络的前几层
        Args:
            freeze_at: 冻结到第几层 (0表示不冻结)
        """
        if freeze_at <= 0:
            return
            
        print(f"冻结骨干网络前{freeze_at}层")
        # 这里简单地冻结resnet_backbone的参数
        for i, module in enumerate(self.resnet_backbone.children()):
            if i < freeze_at:
                for param in module.parameters():
                    param.requires_grad = False
                # 设置为评估模式
                module.eval()
    
    def get_params_lr(self):
        # 返回参数和学习率
        return [{'params': self.parameters(), 'lr': 1.0}]
    
    def named_parameters(self, prefix='', recurse=True):
        return super().named_parameters(prefix=prefix, recurse=recurse)

if is_master:
    logger, final_output_path = create_logger(config.output_path, args.cfg, config.dataset.image_set)
    writer = tensorboardX.SummaryWriter(log_dir=os.path.join(config.output_path, 'tensorboard',
                                                             os.path.basename(args.cfg).split('.')[0],
                                                             '_'.join(config.dataset.image_set.split('+')),
                                                             config.model_prefix+time.strftime('%Y-%m-%d-%H-%M')))
else:
    final_output_path = os.path.join(config.output_path, os.path.basename(args.cfg).split('.')[0], '{}'.format('_'.join([iset for iset in config.dataset.image_set.split('+')])))

# 删除不存在的导入
# from upsnet.dataset import *
# from upsnet.models import *
# from lib.utils.callback import Speedometer
# from lib.utils.data_parallel import DataParallel
# from lib.utils.metric import AvgMetric
# from lib.nn.optimizer import SGD, Adam, clip_grad

np.random.seed(235)
torch.cuda.manual_seed_all(235)
torch.manual_seed(235)

cv2.ocl.setUseOpenCL(False)
cudnn.enabled = True
cudnn.benchmark = False


def lr_poly(base_lr, iter, max_iter, warmup_iter=0):
    power = 0.9
    if iter < warmup_iter:
        alpha = iter / warmup_iter
        return min(base_lr * (1 / 10.0 * (1 - alpha) + alpha), base_lr * ((1 - float(iter) / max_iter)**(power)))
    return base_lr * ((1 - float(iter) / max_iter)**(power))


def get_step_index(iter, decay_iters):
    for idx, decay_iter in enumerate(decay_iters):
        if iter < decay_iter:
            return idx
    return len(decay_iters)


def lr_factor(base_lr, iter, decay_iter, warmup_iter=0):
    if iter < warmup_iter:
        alpha = iter / warmup_iter
        return base_lr * (1 / 10.0 * (1 - alpha) + alpha)
    return base_lr * (0.1 ** get_step_index(iter, decay_iter))


def adjust_learning_rate(optimizer, iter, config):
    assert config.train.lr_schedule in ['step', 'poly']
    if config.train.lr_schedule == 'step':
        return lr_factor(config.train.lr, iter, config.train.decay_iteration, config.train.warmup_iteration)
    if config.train.lr_schedule == 'poly':
        return lr_poly(config.train.lr, iter, config.train.max_iteration, config.train.warmup_iteration)

def freeze(model):
    for name, param in model.named_parameters():
        if 'fpn' in name and param.requires_grad:
            if config.network.fpn_freeze:
                param.requires_grad = False
        elif 'rpn' in name and param.requires_grad:
            if config.network.rpn_freeze:
                param.requires_grad = False
        elif 'resnet_backbone' in name and param.requires_grad:
            if config.network.backbone_freeze:
                param.requires_grad = False
        elif 'rcnn' in name and param.requires_grad:
            if config.network.rcnn_freeze:
                param.requires_grad = False
        elif 'mask_branch' in name and param.requires_grad:
            if config.network.mask_freeze:
                param.requires_grad = False
        else:
            param.requires_grad = False


def upsnet_train():
    # 单GPU设置
    gpus = [torch.device('cuda', int(_)) for _ in config.gpus.split(',')]
    print(f'GPUs: {gpus}')
    num_gpus = len(gpus)

    print(f'Creating model...')
    # create models
    train_model = eval(config.symbol)().cuda()
    
    # freeze params
    freeze(train_model)
        
    # create optimizer
    params_lr = train_model.get_params_lr()
    # we use custom optimizer and pass lr=1 to support different lr for different weights
    optimizer = SGD(params_lr, lr=1, momentum=config.train.momentum, weight_decay=config.train.wd)
    optimizer.zero_grad()
    
    # create data loader - 单GPU简化版本
    train_dataset = eval(config.dataset.dataset)(image_sets=config.dataset.image_set.split('+'), flip=config.train.flip, result_path=final_output_path)
    val_dataset = eval(config.dataset.dataset)(image_sets=config.dataset.test_image_set.split('+'), flip=False, result_path=final_output_path, phase='val')
    
    train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=config.train.batch_size, shuffle=config.train.shuffle, num_workers=4, drop_last=False, collate_fn=train_dataset.collate)
    val_loader = torch.utils.data.DataLoader(val_dataset, batch_size=config.train.batch_size, shuffle=False, num_workers=4, drop_last=False, collate_fn=val_dataset.collate)

    # preparing
    curr_iter = config.train.begin_iteration
    batch_end_callback = [Speedometer(config.train.batch_size, config.train.display_iter)]
    metrics = []
    metrics_name = []
    if config.network.has_rpn:
        metrics.extend([AvgMetric(name='rpn_cls_loss'), AvgMetric(name='rpn_bbox_loss'),])
        metrics_name.extend(['rpn_cls_loss', 'rpn_bbox_loss'])
    if config.network.has_rcnn:
        metrics.extend([AvgMetric(name='rcnn_accuracy'), AvgMetric(name='cls_loss'), AvgMetric(name='bbox_loss'),])
        metrics_name.extend(['rcnn_accuracy', 'cls_loss', 'bbox_loss'])
    if config.network.has_mask_head:
        metrics.extend([AvgMetric(name='mask_loss'), ])
        metrics_name.extend(['mask_loss'])
    if config.network.has_fcn_head:
        metrics.extend([AvgMetric(name='fcn_loss'), ])
        metrics_name.extend(['fcn_loss'])
        if config.train.fcn_with_roi_loss:
            metrics.extend([AvgMetric(name='fcn_roi_loss'), ])
            metrics_name.extend(['fcn_roi_loss'])
    if config.network.has_panoptic_head:
        metrics.extend([AvgMetric(name='panoptic_accuracy'), AvgMetric(name='panoptic_loss'), ])
        metrics_name.extend(['panoptic_accuracy', 'panoptic_loss'])

    if config.train.resume:
        train_model.load_state_dict(torch.load(os.path.join(final_output_path, config.model_prefix + str(curr_iter) + '.pth')))
        optimizer.load_state_dict(torch.load(os.path.join(final_output_path, config.model_prefix + str(curr_iter) + '.state.pth')))
    else:
        if is_master:
            if not config.network.use_pretrained:
                logger.info('Train From Scratch')
            else:
                logger.info(f'Loading pretrained model: {config.network.pretrained}')
                train_model.load_state_dict(torch.load(config.network.pretrained))

    # 不使用数据并行，直接使用单GPU
    train_model = train_model.to(gpus[0])

    if is_master:
        batch_end_callback[0](0, metrics)

    train_model.eval()

    # write freeze situation of model to verify configuration  
    freeze_dir = os.path.join(final_output_path, 'freeze')
    os.makedirs(freeze_dir, exist_ok=True)
    with open(os.path.join(freeze_dir, config.model_prefix+'freeze_config.txt'), 'w', encoding='utf-8') as f:
        for name, param in train_model.named_parameters():
            f.write(f'{name}\t\t{param.requires_grad}\n')

    # start training
    train_model.train()

    while curr_iter < config.train.max_iteration:
        for inner_iter, batch in enumerate(train_loader):
            data, label, _ = batch
            
            # Move data to GPU
            for k, v in data.items():
                data[k] = v if not torch.is_tensor(v) else v.cuda()
            for k, v in label.items():
                label[k] = v if not torch.is_tensor(v) else v.cuda()

            lr = adjust_learning_rate(optimizer, curr_iter, config)
            optimizer.zero_grad()
            output = train_model(data, label)
            
            # Calculate loss
            loss = torch.tensor(0.0, requires_grad=True, device=gpus[0])
            if config.network.has_rpn:
                loss = loss + output['rpn_cls_loss'].mean() + output['rpn_bbox_loss'].mean()
            if config.network.has_rcnn:
                loss = loss + output['cls_loss'].mean() + output['bbox_loss'].mean() * config.train.bbox_loss_weight
            if config.network.has_mask_head:
                loss = loss + output['mask_loss'].mean()
            if config.network.has_fcn_head:
                loss = loss + output['fcn_loss'].mean() * config.train.fcn_loss_weight
                if config.train.fcn_with_roi_loss:
                    loss = loss + output['fcn_roi_loss'].mean() * config.train.fcn_loss_weight * 0.2
            if config.network.has_panoptic_head:
                loss = loss + output['panoptic_loss'].mean() * config.train.panoptic_loss_weight
            
            loss.backward()
            optimizer.step(lr)

            # Update metrics
            if is_master:
                if writer:
                    writer.add_scalar('train_total_loss', loss.item(), curr_iter)
                for metric, l in zip(metrics, metrics_name):
                    loss_val = output[l].mean().item() if l in output else 0.0
                    if writer:
                        writer.add_scalar('train_' + l, loss_val, curr_iter)
                    metric.update(0, 0, loss_val)  # 使用0作为占位符
                    
            curr_iter += 1

            # Decay momentum buffer
            if curr_iter in config.train.decay_iteration:
                if is_master:
                    logger.info('decay momentum buffer')
                for k in optimizer.state_dict()['state'].keys():
                    if 'momentum_buffer' in optimizer.state_dict()['state'][k]:
                        optimizer.state_dict()['state'][k]['momentum_buffer'].div_(10)

            # Display and save
            if is_master:
                if curr_iter % config.train.display_iter == 0:
                    for callback in batch_end_callback:
                        callback(curr_iter, metrics)

                if curr_iter % config.train.snapshot_step == 0:
                    logger.info('taking snapshot ...')
                    torch.save(train_model.state_dict(), os.path.join(final_output_path, config.model_prefix+str(curr_iter)+'.pth'))
                    torch.save(optimizer.state_dict(), os.path.join(final_output_path, config.model_prefix+str(curr_iter)+'.state.pth'))
                    
            if curr_iter >= config.train.max_iteration:
                break

        for metric in metrics:
            metric.reset()

        if config.train.eval_data:
            train_model.eval()
            # 简化的验证循环
            for inner_iter, batch in enumerate(val_loader):
                if inner_iter >= 10:  # 只验证前10个批次
                    break
                data, label, _ = batch
                for k, v in data.items():
                    data[k] = v if not torch.is_tensor(v) else v.cuda()
                for k, v in label.items():
                    label[k] = v if not torch.is_tensor(v) else v.cuda()

                with torch.no_grad():
                    output = train_model(data, label)

                for metric, l in zip(metrics, metrics_name):
                    if l in output:
                        loss_val = output[l].mean().item()
                        metric.update(0, 0, loss_val)

            s = 'Batch [%d]\t Epoch[%d]\t' % (curr_iter, curr_iter // len(train_loader))
            for metric in metrics:
                m, v = metric.get()
                s += 'Val-%s=%f,\t' % (m, v)
                if is_master and writer:
                    writer.add_scalar('val_' + m, v, curr_iter)
                    metric.reset()
            if is_master:
                logger.info(s)
            train_model.train()

    # 最终保存模型
    if is_master:
        logger.info('taking snapshot ...')
        torch.save(train_model.state_dict(), os.path.join(final_output_path, config.model_prefix + str(curr_iter) + '.pth'))
        torch.save(optimizer.state_dict(), os.path.join(final_output_path, config.model_prefix + str(curr_iter) + '.state.pth'))

if __name__ == '__main__':
    upsnet_train()
