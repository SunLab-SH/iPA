# -*- coding: utf-8 -*-
# @Time    : 2020-06-02 23:08
# @Author  : Xiangyi Zhang
# @File    : Unet.py
# @Email   : zhangxy9@shanghaitech.edu.cn

import torch.nn as nn
import torch.nn.functional as F
import torch
from networks.Block import Convblock, Upconv, UpconvWithCat


class Unet(nn.Module):
    """
    UNet - Basic Implementation
    Paper : https://arxiv.org/abs/1505.04597
    """

    def __init__(self, in_channel=1, n_class=2, is_dropout=False):
        super(Unet, self).__init__()

        n1 = 64
        filters = [n1, n1 * 2, n1 * 4, n1 * 8, n1 * 16]
        self.is_dropout = is_dropout
        self.Maxpool = nn.MaxPool2d(kernel_size=2, stride=2)

        self.conv1 = Convblock(in_channel, filters[0])
        self.conv2 = Convblock(filters[0], filters[1], downsample=True)
        self.conv3 = Convblock(filters[1], filters[2], downsample=True)
        self.conv4 = Convblock(filters[2], filters[3], downsample=True)
        self.conv5 = Convblock(filters[3], filters[4], downsample=True)

        self.up_concat5 = UpconvWithCat(filters[4], filters[3])

        self.up_concat4 = UpconvWithCat(filters[3], filters[2])

        self.up_concat3 = UpconvWithCat(filters[2], filters[1])

        self.up_concat2 = UpconvWithCat(filters[1], filters[0])

        self.Conv = nn.Conv2d(filters[0], n_class, kernel_size=1, stride=1, padding=0)
        self.dropout = nn.Dropout(0.5)

    def encoder(self, x):
        x1 = self.conv1(x)
        x2 = self.conv2(x1)
        x3 = self.conv3(x2)
        x4 = self.conv4(x3)
        x5 = self.conv5(x4)
        if self.is_dropout:
            x5 = self.dropout(x5)
        # print(x1.shape, x2.shape, x3.shape, x4.shape, x5.shape)
        return [x1, x2, x3, x4, x5]

    def decoder(self, feat):
        x1, x2, x3, x4, x5 = feat

        out5 = self.up_concat5(x5, x4)
        out4 = self.up_concat4(out5, x3)
        out3 = self.up_concat3(out4, x2)
        out2 = self.up_concat2(out3, x1)
        # print(out2.shape, out3.shape, out4.shape, out5.shape)

        return out2

    def forward(self, inputs):
        feat = self.encoder(inputs)
        out = self.decoder(feat)
        if self.is_dropout:
            out = self.dropout(out)
        out = self.Conv(out)

        return out


class UnetNested(nn.Module):

    def __init__(self, in_channels=1, n_classes=2, is_deconv=True, is_deep_supervised=True):
        super(UnetNested, self).__init__()
        self.is_deconv = is_deconv
        self.is_deep_supervised = is_deep_supervised

        filters = [64, 128, 256, 512, 1024]

        # downsampling
        self.conv00 = Convblock(in_channels, filters[0])
        self.conv10 = Convblock(filters[0], filters[1], is_conv_downsample=False)
        self.conv20 = Convblock(filters[1], filters[2], is_conv_downsample=False)
        self.conv30 = Convblock(filters[2], filters[3], is_conv_downsample=False)
        self.conv40 = Convblock(filters[3], filters[4], is_conv_downsample=False)

        # upsampling
        self.up_concat01 = UpconvWithCat(filters[1], filters[0])
        self.up_concat11 = UpconvWithCat(filters[2], filters[1])
        self.up_concat21 = UpconvWithCat(filters[3], filters[2])
        self.up_concat31 = UpconvWithCat(filters[4], filters[3])

        self.up_concat02 = UpconvWithCat(filters[1], filters[0], n_concat=3)
        self.up_concat12 = UpconvWithCat(filters[2], filters[1], n_concat=3)
        self.up_concat22 = UpconvWithCat(filters[3], filters[2], n_concat=3)

        self.up_concat03 = UpconvWithCat(filters[1], filters[0], n_concat=4)
        self.up_concat13 = UpconvWithCat(filters[2], filters[1], n_concat=4)

        self.up_concat04 = UpconvWithCat(filters[1], filters[0], n_concat=5)

        # final conv (without any concat)
        self.final_1 = nn.Conv2d(filters[0], n_classes, 1)
        self.final_2 = nn.Conv2d(filters[0], n_classes, 1)
        self.final_3 = nn.Conv2d(filters[0], n_classes, 1)
        self.final_4 = nn.Conv2d(filters[0], n_classes, 1)

    def forward(self, inputs):
        # column : 0
        X_00 = self.conv00(inputs)
        X_10 = self.conv10(X_00)
        X_20 = self.conv20(X_10)
        X_30 = self.conv30(X_20)
        X_40 = self.conv40(X_30)
        # column : 1
        X_01 = self.up_concat01(X_10, X_00)
        X_11 = self.up_concat11(X_20, X_10)
        X_21 = self.up_concat21(X_30, X_20)
        X_31 = self.up_concat31(X_40, X_30)
        # column : 2
        X_02 = self.up_concat02(X_11, X_00, X_01)
        X_12 = self.up_concat12(X_21, X_10, X_11)
        X_22 = self.up_concat22(X_31, X_20, X_21)
        # column : 3
        X_03 = self.up_concat03(X_12, X_00, X_01, X_02)
        X_13 = self.up_concat13(X_22, X_10, X_11, X_12)
        # column : 4
        X_04 = self.up_concat04(X_13, X_00, X_01, X_02, X_03)

        # final layer
        final_1 = self.final_1(X_01)
        final_2 = self.final_2(X_02)
        final_3 = self.final_3(X_03)
        final_4 = self.final_4(X_04)

        final = (final_1 + final_2 + final_3 + final_4) / 4

        if self.is_deep_supervised:
            return final
        else:
            return final_4

