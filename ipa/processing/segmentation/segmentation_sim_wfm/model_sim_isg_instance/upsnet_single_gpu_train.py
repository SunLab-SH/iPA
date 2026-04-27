# ---------------------------------------------------------------------------
# Unified Panoptic Segmentation Network - Single GPU Training
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
# Modified for single GPU training
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

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from collections import deque
from upsnet.config.config import config
from upsnet.config.parse_args import parse_args
import glob
from PIL import Image
import random

# 导入本地utils
sys.path.append(os.path.join(os.path.dirname(__file__), 'utils'))
from utils.logging import create_logger

# 解析命令行参数，如果没有提供则使用默认配置
try:
    args = parse_args()
except SystemExit:
    # 如果没有提供命令行参数，创建默认参数
    import argparse
    from argparse import Namespace
    args = Namespace()
    args.cfg = os.path.join(os.path.dirname(__file__), 'configs', 'sxt_train.yaml')
    args.eval_only = False
    args.weight_path = None
    print(f"使用默认配置文件: {args.cfg}")

# 更新配置
if hasattr(args, 'cfg') and args.cfg and os.path.exists(args.cfg):
    from upsnet.config.config import update_config
    update_config(args.cfg)
    print(f"已加载配置文件: {args.cfg}")
else:
    print("使用默认配置")

# Single GPU setup
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
logger, final_output_path = create_logger(config.output_path, getattr(args, 'cfg', 'default'), config.dataset.image_set)

# 数据路径配置
DATA_ROOT = r"D:\Gitspace\ipa_full\iPA\data\sxt_images\split_img"
IMAGES_DIR = os.path.join(DATA_ROOT, "images")
MASKS_DIR = os.path.join(DATA_ROOT, "masks")

def check_data_samples():
    """检查3张随机样本确保图像和掩码对应"""
    print("\n=== 检查数据样本 ===")
    
    # 获取所有子文件夹
    image_folders = [f for f in os.listdir(IMAGES_DIR) if os.path.isdir(os.path.join(IMAGES_DIR, f))]
    
    if not image_folders:
        print("错误：未找到图像文件夹")
        return
    
    # 随机选择3个样本
    for i in range(3):
        folder = random.choice(image_folders)
        img_folder_path = os.path.join(IMAGES_DIR, folder)
        mask_folder_path = os.path.join(MASKS_DIR, folder)
        
        if not os.path.exists(mask_folder_path):
            print(f"警告：掩码文件夹不存在: {mask_folder_path}")
            continue
            
        # 获取该文件夹下的图像文件
        img_files = [f for f in os.listdir(img_folder_path) if f.endswith(('.png', '.jpg', '.jpeg'))]
        
        if not img_files:
            continue
            
        img_file = random.choice(img_files)
        img_path = os.path.join(img_folder_path, img_file)
        mask_path = os.path.join(mask_folder_path, img_file)
        
        print(f"\n样本 {i+1}:")
        print(f"图像: {img_path}")
        print(f"掩码: {mask_path}")
        print(f"图像存在: {os.path.exists(img_path)}")
        print(f"掩码存在: {os.path.exists(mask_path)}")
        
        if os.path.exists(img_path) and os.path.exists(mask_path):
            try:
                img = Image.open(img_path)
                mask = Image.open(mask_path)
                print(f"图像尺寸: {img.size}")
                print(f"掩码尺寸: {mask.size}")
                print(f"图像模式: {img.mode}")
                print(f"掩码模式: {mask.mode}")
            except Exception as e:
                print(f"读取错误: {e}")

from upsnet.dataset import *
from upsnet.models import *
# 使用本地utils
from utils.metrics import AvgMetric
from utils.callbacks import Speedometer
from utils.optimizers import SGD, Adam, clip_grad

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
    # 在训练开始前检查数据样本
    check_data_samples()
    
    logger.info('training config:{}\n'.format(pprint.pformat(config)))
    print(f'Using device: {device}')
    print(f'模型保存路径: {final_output_path}')
    
    # 确保输出目录存在
    os.makedirs(final_output_path, exist_ok=True)
    print(f'输出目录已创建: {final_output_path}')

    print(f'Creating model...')
    # create model
    train_model = eval(config.symbol)().to(device)
    
    # freeze params
    freeze(train_model)
        
    # create optimizer
    params_lr = train_model.get_params_lr()
    # we use custom optimizer and pass lr=1 to support different lr for different weights
    optimizer = SGD(params_lr, lr=1, momentum=config.train.momentum, weight_decay=config.train.wd)
    optimizer.zero_grad()

    # create data loader
    train_dataset = eval(config.dataset.dataset)(image_sets=config.dataset.image_set.split('+'), flip=config.train.flip, result_path=final_output_path)
    val_dataset = eval(config.dataset.dataset)(image_sets=config.dataset.test_image_set.split('+'), flip=False, result_path=final_output_path, phase='val')
    
    train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=config.train.batch_size, shuffle=config.train.shuffle, 
                                             num_workers=4, drop_last=False, collate_fn=train_dataset.collate)
    val_loader = torch.utils.data.DataLoader(val_dataset, batch_size=config.train.batch_size, shuffle=False, 
                                            num_workers=4, drop_last=False, collate_fn=val_dataset.collate)

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
        train_model.load_state_dict(torch.load(os.path.join(final_output_path, config.model_prefix + str(curr_iter) + '.pth')), resume=True)
        optimizer.load_state_dict(torch.load(os.path.join(final_output_path, config.model_prefix + str(curr_iter) + '.state.pth')))
    else:
        if not config.network.use_pretrained:
            logger.info('Train From Scratch')
        else:
            logger.info(f'Loading pretrained model: {config.network.pretrained}')
            train_model.load_state_dict(torch.load(config.network.pretrained))

    batch_end_callback[0](0, 0)

    # write freeze situation of model to verify configuration
    freeze_dir = os.path.join(final_output_path, 'freeze')
    os.makedirs(freeze_dir, exist_ok=True)
    with open(os.path.join(freeze_dir, config.model_prefix+'freeze_config.txt'), 'w', encoding='utf-8') as f:
        for name, param in train_model.named_parameters():
            f.write(f'{name}\t\t{param.requires_grad}\n')

    # start training
    while curr_iter < config.train.max_iteration:
        train_model.train()
        
        if config.network.use_syncbn:
            if config.network.backbone_freeze_at > 0:
                train_model.freeze_backbone(config.network.backbone_freeze_at)
            if config.network.backbone_fix_bn:
                train_model.resnet_backbone.eval()

        for inner_iter, batch in enumerate(train_loader):
            data, label, _ = batch
            
            # move data to device
            for k, v in data.items():
                data[k] = v if not torch.is_tensor(v) else v.to(device, non_blocking=True)
            for k, v in label.items():
                label[k] = v if not torch.is_tensor(v) else v.to(device, non_blocking=True)

            lr = adjust_learning_rate(optimizer, curr_iter, config)
            optimizer.zero_grad()
            output = train_model(data, label)
            
            loss = 0
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
            optimizer.step_with_lr(lr)

            # record losses for metrics
            for i, (metric, l) in enumerate(zip(metrics, metrics_name)):
                if l in output and output[l] is not None:
                    loss_val = output[l].mean().item()
                    metric.update(None, None, loss_val)
            
            curr_iter += 1

            if curr_iter in config.train.decay_iteration:
                logger.info('decay momentum buffer')
                for k in optimizer.state_dict()['state'].keys():
                    if 'momentum_buffer' in optimizer.state_dict()['state'][k]:
                        optimizer.state_dict()['state'][k]['momentum_buffer'].div_(10)

            if curr_iter % config.train.display_iter == 0:
                for callback in batch_end_callback:
                    callback(curr_iter, metrics)
                
                # 输出训练指标
                print(f"\n=== 训练迭代 {curr_iter} ===")
                for metric in metrics:
                    m, v = metric.get()
                    print(f"Train_{m}: {v:.6f}")
                    logger.info(f'Train_{m}: {v:.6f}')

            if curr_iter % config.train.snapshot_step == 0:
                logger.info('taking snapshot ...')
                model_path = os.path.join(final_output_path, config.model_prefix+str(curr_iter)+'.pth')
                optimizer_path = os.path.join(final_output_path, config.model_prefix+str(curr_iter)+'.state.pth')
                torch.save(train_model.state_dict(), model_path)
                torch.save(optimizer.state_dict(), optimizer_path)
                print(f"模型已保存到: {model_path}")
                print(f"优化器状态已保存到: {optimizer_path}")

            if curr_iter >= config.train.max_iteration:
                break

        for metric in metrics:
            metric.reset()

        # validation
        if config.train.eval_data:
            train_model.eval()

            for inner_iter, batch in enumerate(val_loader):
                data, label, _ = batch
                
                # move data to device
                for k, v in data.items():
                    data[k] = v if not torch.is_tensor(v) else v.to(device, non_blocking=True)
                for k, v in label.items():
                    label[k] = v if not torch.is_tensor(v) else v.to(device, non_blocking=True)

                with torch.no_grad():
                    output = train_model(data, label)

                for metric, l in zip(metrics, metrics_name):
                    if l in output and output[l] is not None:
                        loss_val = output[l].mean().item()
                        metric.update(None, None, loss_val)

            s = 'Batch [%d]\t Epoch[%d]\t' % (curr_iter, curr_iter // len(train_loader))

            for metric in metrics:
                m, v = metric.get()
                s += 'Val-%s=%f,\t' % (m, v)
                logger.info(f'Val_{m}: {v:.6f}')
                print(f"Val_{m}: {v:.6f}")  # 添加打印输出
                metric.reset()
            logger.info(s)
            print(s)  # 添加打印输出

    # final snapshot
    logger.info('taking final snapshot ...')
    torch.save(train_model.state_dict(), os.path.join(final_output_path, config.model_prefix + str(curr_iter) + '.pth'))
    torch.save(optimizer.state_dict(), os.path.join(final_output_path, config.model_prefix + str(curr_iter) + '.state.pth'))


if __name__ == '__main__':
    upsnet_train()