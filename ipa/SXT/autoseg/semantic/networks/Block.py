# -*- coding: utf-8 -*-
# @Time    : 2020-06-02 23:08
# @Author  : Xiangyi Zhang
# @File    : Block.py
# @Email   : zhangxy9@shanghaitech.edu.cn

import torch.nn as nn
import torch.nn.functional as F
import torch


class Convblock(nn.Module):
    """
    Convolution Block
    """

    def __init__(self, in_channel, out_channel, downsample=False):
        super(Convblock, self).__init__()
        self.downsample = downsample
        self.conv = nn.Sequential(
            nn.Conv2d(in_channel, out_channel, kernel_size=3, stride=1, padding=1, bias=True),
            nn.BatchNorm2d(out_channel),
            nn.ReLU(inplace=True),
            nn.Conv2d(out_channel, out_channel, kernel_size=3, stride=1, padding=1, bias=True),
            nn.BatchNorm2d(out_channel),
            nn.ReLU(inplace=True))

        self.maxpool = nn.MaxPool2d(kernel_size=2)

    def forward(self, x):
        if self.downsample:
            x = self.maxpool(x)

        x = self.conv(x)
        return x


class Upconv(nn.Module):
    """
    Up Convolution Block
    """
    def __init__(self, in_channel, out_channel):
        super(Upconv, self).__init__()
        self.up = nn.Sequential(
            nn.Upsample(scale_factor=2),
            nn.Conv2d(in_channel, out_channel, kernel_size=3, stride=1, padding=1, bias=True),
            nn.BatchNorm2d(out_channel),
            nn.ReLU(inplace=True)
        )

    def forward(self, x):
        x = self.up(x)
        return x


class UpconvWithCat(nn.Module):
    def __init__(self, in_channel, out_channel, n_concat=2):
        super(UpconvWithCat, self).__init__()
        self.conv = Convblock(in_channel + (n_concat - 2) * out_channel, out_channel)
        self.up = Upconv(in_channel, out_channel)

    def forward(self, high_feature, *low_feature):
        outputs = self.up(high_feature)
        for feature in low_feature:
            outputs = torch.cat([outputs, feature], 1)
        return self.conv(outputs)