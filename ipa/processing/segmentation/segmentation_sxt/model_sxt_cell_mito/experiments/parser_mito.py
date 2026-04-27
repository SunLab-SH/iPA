# -*- coding: utf-8 -*-
# @Time    : 2020-06-01 12:08
# @Author  : Xiangyi Zhang
# @File    : parser.py
# @Email   : zhangxy9@shanghaitech.edu.cn
import argparse
import logging


LOG = logging.getLogger('main')

def create_parser():
    """创建参数解析器"""
    parser = argparse.ArgumentParser(description='PyTorch LENF Training')
    parser.add_argument('--exp', type=str, default='FS_0106_001', help='the experiment name')
    parser.add_argument('--seed', type=int, default=2019, help='random seed')
    parser.add_argument('--gpu', type=str, default='0', help='GPU to use')
    parser.add_argument('--has-dropout', type=bool, default=True, help='dropout or not(Unused)')
    # conv downsample maybe help the samll object, but failed. So disable it now.
    parser.add_argument('--is-conv-downsample', type=bool, default=False, help='conv-downsample or not')
    parser.add_argument('--all-label', type=bool, default=True, help='True: 5 label, False: 3 label')

    # dataset setting
    parser.add_argument('--root-dir', type=str, default='F:\\salilab\\salilab\\projects\\auto-segmentation_Aneesh\\autoseg\\semantic', help='root path to dataset')
    parser.add_argument('--data-root-dir', type=str, default='F:\\salilab\\salilab\\projects\\auto-segmentation_Aneesh\\Data\\raw_img', help='root path to dataset')
    parser.add_argument('--num-workers', default=4, type=int, help='number of data loading workers (default: 4)')
    parser.add_argument('--num-classes', default=2, type=int, help='number of class(default 5)')
    parser.add_argument('--batch-size', default=4, type=int, help='mini-batch size (default: 4)')
    parser.add_argument('--label-batch', default=2, type=int, help='how many label image in a batch (default: 2)')

    # learning setting
    parser.add_argument('--step-size', type=int, default=15, help='every step size to decay the learning rate')
    parser.add_argument('--num-epochs', default=60, type=int, help='num of epoch')
    parser.add_argument('--lr', '--learning-rate', default=1e-4, type=float, help='max learning rate')
    parser.add_argument('--momentum', default=0.99, type=float, metavar='M', help='momentum')
    parser.add_argument('--weight-decay', default=1e-4, type=float, help='weight decay (default: 1e-4)')
    parser.add_argument('--ema-decay', type=float,  default=0.99, help='ema decay')
    parser.add_argument('--consistency_type', type=str,  default="mse", help='consistency_type')
    parser.add_argument('--consistency', type=float,  default=0.1, help='consistency')
    parser.add_argument('--consistency_rampup', type=float,  default=40.0, help='consistency_rampup')

    # checkpoint and evaluation
    parser.add_argument('--print-freq', default=50, type=int, help='print frequency (default: 10)')
    parser.add_argument('--epoch_val', default=10, type=int, help='epoch_val')
    parser.add_argument('--contrast', default=0, type=int, help='contrast')
    parser.add_argument('--test_idx', default='all', type=str, help='tested data id')

    ## Evaluation
    parser.add_argument('--mode', default='mito', type=str, help='eval or postprocess')
    parser.add_argument('--post', default=True, type=bool, help='do postprocess or not(including 3d fusion and coarse label refine)')
    
    return parser

def get_base_args():
    """获取基础参数配置"""
    return ['--gpu', '0', 
            '--num-workers', '2', 
            '--batch-size', '1', 
            '--num-classes', '2', 
            '--test_idx', 'iso', 
            '--mode', 'mem_nu']

def get_args():
    """获取默认参数"""
    parser = create_parser()
    base_args = get_base_args()
    return parser.parse_args(args=base_args)

def get_args_nu():
    """获取nu模式参数"""
    parser = create_parser()
    base_args = get_base_args()
    args_nu_list = base_args + ['--exp', "FS_mem_nu"]
    return parser.parse_args(args=args_nu_list)

def get_args_mito():
    """获取mito模式参数"""
    parser = create_parser()
    base_args = get_base_args()
    args_mito_list = base_args + ['--exp', "FS_mito"]
    return parser.parse_args(args=args_mito_list)

# 为了向后兼容，保留这些变量
args = None
args_nu = None
args_mito = None

def init_args():
    """初始化参数（延迟初始化）"""
    global args, args_nu, args_mito
    if args is None:
        args = get_args()
        args_nu = get_args_nu()
        args_mito = get_args_mito()
    return args, args_nu, args_mito

