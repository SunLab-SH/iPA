# -*- coding: utf-8 -*-
# @Time    : 2023-06-01 15:58
# @File    : parser.py
# @Author   : Angdi

import argparse
import logging


LOG = logging.getLogger('main')
parser = argparse.ArgumentParser(description='autoseg vesicle on SIM data')

# parser.add_argument('--datapath', type=str, default='I:\\Bing\\fluorescence\\3D\\20230320MT-ISG\\03_versicle', help='path to dataset')
# parser.add_argument('--resultpath', type=str, default='I:\\Bing\\fluorescence\\3D\\20230320MT-ISG\\03_versicle', help='path to save results')
parser.add_argument('--datapath', type=str, default='I:\\Bing\\fluorescence\\3D\\03_versicle', help='path to dataset')
parser.add_argument('--resultpath', type=str, default='I:\\Bing\\fluorescence\\3D\\03_versicle', help='path to save results')


args = parser.parse_args(args=[])
# args = parser.parse_args(args=[ '--exp', "FS_mem_nu", 
#                                 '--gpu', '0', 
#                                 '--num-workers', '2', 
#                                 '--batch-size', '1', 
#                                 '--num-classes', '2', 
#                                 '--test_idx', 'iso', 
#                                 '--mode', 'mem_nu',])
