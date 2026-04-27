# -*- coding: utf-8 -*-
"""
简化的配置文件用于演示
Simplified configuration for demo
"""

import argparse
import os

class DemoArgs:
    """演示用参数配置类"""
    def __init__(self):
        # 基本设置
        self.exp = 'FS_mem_nu_demo'
        self.seed = 2019
        self.gpu = '0'
        self.has_dropout = True
        self.is_conv_downsample = False
        self.all_label = True
        
        # 数据设置
        self.root_dir = r'f:\Salilab\Projects\auto-segmentation_Aneesh\autoseg\semantic'
        self.data_root_dir = r'f:\Salilab\Projects\auto-segmentation_Aneesh\Data\raw_img'
        self.num_workers = 2
        self.num_classes = 2  # membrane + nucleus
        self.batch_size = 1
        self.label_batch = 1
        
        # 学习设置
        self.step_size = 15
        self.num_epochs = 60
        self.lr = 1e-4
        self.momentum = 0.99
        self.weight_decay = 1e-4
        self.ema_decay = 0.99
        self.consistency_type = "mse"
        self.consistency = 0.1
        self.consistency_rampup = 40.0
        
        # 检查点和评估
        self.print_freq = 50
        self.epoch_val = 10
        self.contrast = 0
        self.test_idx = 'demo'
        
        # 评估模式
        self.mode = 'mem_nu'  # 'mem_nu' 或 'mito'
        self.post = True

# 创建全局args对象供导入使用
args = DemoArgs()
