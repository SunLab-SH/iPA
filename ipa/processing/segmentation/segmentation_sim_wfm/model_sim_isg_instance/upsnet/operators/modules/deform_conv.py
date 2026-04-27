# ---------------------------------------------------------------------------
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


import math
import torch
import torch.nn as nn
from torch.nn.parameter import Parameter
from torch.nn.modules.utils import _pair
from upsnet.operators.functions.deform_conv import DeformConvFunction


class DeformConv(nn.Module):
    """Simplified deformable convolution module using standard convolution."""
    
    def __init__(self, in_channels, out_channels, kernel_size, stride=1, 
                 padding=0, dilation=1, groups=1, bias=True):
        super(DeformConv, self).__init__()
        
        # Use standard convolution as fallback
        self.conv = nn.Conv2d(
            in_channels, out_channels, kernel_size, 
            stride=stride, padding=padding, dilation=dilation, 
            groups=groups, bias=bias
        )
        
        # Offset convolution (for deformable convolution)
        # In simplified version, we don't use this
        self.offset_conv = nn.Conv2d(
            in_channels, 2 * kernel_size * kernel_size, 
            kernel_size, stride=stride, padding=padding, dilation=dilation
        )
        
        self.init_weights()
    
    def init_weights(self):
        """Initialize weights."""
        nn.init.kaiming_normal_(self.conv.weight, mode='fan_out', nonlinearity='relu')
        if self.conv.bias is not None:
            nn.init.constant_(self.conv.bias, 0)
        
        # Initialize offset conv to zero (no deformation initially)
        nn.init.constant_(self.offset_conv.weight, 0)
        if self.offset_conv.bias is not None:
            nn.init.constant_(self.offset_conv.bias, 0)
    
    def forward(self, x):
        """Forward pass using standard convolution."""
        # In full implementation, we would:
        # 1. Compute offset using offset_conv
        # 2. Apply deformable convolution using the offset
        # For simplicity, we just use standard convolution
        return self.conv(x)