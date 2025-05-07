# -*- coding: utf-8 -*-
# @Time    : 2020-06-02 23:08
# @Author  : Xiangyi Zhang
# @File    : net.py
# @Email   : zhangxy9@shanghaitech.edu.cn
from torchvision.models.resnet import BasicBlock, ResNet, Bottleneck
import torch.nn as nn
import torch
import torch.nn.functional as F


class Encoder(ResNet):
    def __init__(self, name="resnet101", is_dropout=False):
        self.is_dropout = is_dropout
        component = {"resnet150": [3, 8, 36, 3], "resnet101": [3, 4, 23, 3], "resnet50": [3, 4, 6, 3],
                     "resnet34": [3, 4, 6, 3], "resnet18": [2, 2, 2, 2]}
        if name == "resnet34":
            super(Encoder, self).__init__(BasicBlock, component[name], 1000)
        else:
            super(Encoder, self).__init__(Bottleneck, component[name], 1000)

        self.dropout = nn.Dropout(0.5)

    def forward(self, x):
        x = self.conv1(x)
        x = self.bn1(x)
        x1 = self.relu(x)
        x = self.maxpool(x1)
        x2 = self.layer1(x)
        x3 = self.layer2(x2)
        x4 = self.layer3(x3)
        x5 = self.layer4(x4)
        if self.is_dropout:
            x5 = self.dropout(x5)
        return x5, x4, x3, x2, x1

"""
class Bottleneck(nn.Module):
    expansion = 4

    def __init__(self, inplanes, planes, stride=1):
        super(Bottleneck, self).__init__()
        self.conv1 = nn.Conv2d(inplanes, planes, kernel_size=1, bias=True)
        self.bn1 = nn.BatchNorm2d(planes)
        self.conv2 = nn.Conv2d(planes, planes, kernel_size=3, stride=stride,
                               padding=1, bias=False)
        self.bn2 = nn.BatchNorm2d(planes)
        self.conv3 = nn.Conv2d(planes, planes * 4, kernel_size=1, bias=True)
        self.bn3 = nn.BatchNorm2d(planes * 4)
        self.relu = nn.ReLU(inplace=True)
        self.stride = stride
        self.shortconv = nn.Conv2d(inplanes, planes*4, kernel_size=1, bias=True)

    def forward(self, x):
        residual = x

        out = self.conv1(x)
        out = self.bn1(out)
        out = self.relu(out)

        out = self.conv2(out)
        out = self.bn2(out)
        out = self.relu(out)

        out = self.conv3(out)
        out = self.bn3(out)
        out += self.shortconv(residual)
        out = self.relu(out)
        return out
"""


class Decoder1(nn.Module):
    def __init__(self, is_init_decoder=True, is_dropout=False, n_class=5):
        super(Decoder1, self).__init__()
        self.is_init_decoder = is_init_decoder
        self.uplayer4 = Bottleneck(inplanes=2048 + 1024, planes=512, stride=1)
        self.uplayer3 = Bottleneck(inplanes=2048 + 512, planes=256, stride=1)
        self.uplayer2 = Bottleneck(inplanes=1024 + 256, planes=128, stride=1)
        self.uplayer1 = Bottleneck(inplanes=512 + 64, planes=64, stride=1)
        self.uplayer0 = Bottleneck(inplanes=256, planes=32, stride=1)
        self.outconv = nn.Conv2d(in_channels=128, out_channels=n_class, kernel_size=1, bias=True, stride=1)
        self.is_dropout = is_dropout
        self.dropout = nn.Dropout(0.5)
        if self.is_init_decoder:
            for m in self.modules():
                if isinstance(m, nn.Conv2d):
                    torch.nn.init.kaiming_normal_(m.weight, mode='fan_in')
                elif isinstance(m, nn.BatchNorm2d):
                    m.weight.data.fill_(1)
            print('init the decoder parameters')

    def forward(self, feat, residual):
        x5, x4, x3, x2, x1 = feat[4], feat[3], feat[2], feat[1], feat[0]
        x5_up = F.interpolate(x5, size=x4.size()[2:], mode='bilinear')
        x5_cat = torch.cat((x4, x5_up), dim=1)
        x5_out = self.uplayer4(x5_cat)

        x4_up = F.interpolate(x5_out, size=x3.size()[2:], mode='bilinear')
        x4_cat = torch.cat((x3, x4_up), dim=1)
        x4_out = self.uplayer3(x4_cat)

        x3_up = F.interpolate(x4_out, size=x2.size()[2:], mode='bilinear')
        x3_cat = torch.cat((x2, x3_up), dim=1)
        x3_out = self.uplayer2(x3_cat)

        x2_up = F.interpolate(x3_out, size=x1.size()[2:], mode='bilinear')

        x2_cat = torch.cat((x1, x2_up), dim=1)
        x2_out = self.uplayer1(x2_cat)

        x1_up = F.interpolate(x2_out, size=residual.size()[2:], mode='bilinear')
        x1_out = self.uplayer0(x1_up)

        if self.is_dropout:
            x1_out = self.dropout(x1_out)
        out = self.outconv(x1_out)
        return out


class Decoder2(nn.Module):

    def __init__(self, is_init_decoder=True):
        super(Decoder2, self).__init__()
        n_class = 5
        self.is_init_decoder = is_init_decoder
        self.shortconv4 = nn.Conv2d(in_channels=1024, out_channels=1024, kernel_size=3, padding=1, stride=1, bias=True)
        self.shortbn4 = nn.BatchNorm2d(num_features=1024)
        self.shortconv3 = nn.Conv2d(in_channels=512, out_channels=512, kernel_size=3, padding=1, stride=1, bias=True)
        self.shortbn3 = nn.BatchNorm2d(num_features=512)
        self.shortconv2 = nn.Conv2d(in_channels=256, out_channels=256, kernel_size=3, padding=1, stride=1, bias=True)
        self.shortbn2 = nn.BatchNorm2d(num_features=256)
        self.shortconv1 = nn.Conv2d(in_channels=64, out_channels=64, kernel_size=3, padding=1, stride=1, bias=True)
        self.shortbn1 = nn.BatchNorm2d(num_features=64)

        self.uplayer4 = Bottleneck(inplanes=2048 + 1024, planes=512, stride=1)
        self.uplayer3 = Bottleneck(inplanes=2048 + 512, planes=256, stride=1)
        self.uplayer2 = Bottleneck(inplanes=1024 + 256, planes=128, stride=1)
        self.uplayer1 = Bottleneck(inplanes=512 + 64, planes=64, stride=1)
        self.uplayer0 = Bottleneck(inplanes=256, planes=32, stride=1)

        self.outconv0 = nn.Conv2d(in_channels=128, out_channels=128, kernel_size=3, bias=True, stride=1, padding=1)
        self.outconvs0 = nn.Conv2d(in_channels=128, out_channels=n_class, kernel_size=1, bias=True, stride=1)
        self.outconv1 = nn.Conv2d(in_channels=256, out_channels=128, kernel_size=3, bias=True, stride=1, padding=1)
        self.outconvs1 = nn.Conv2d(in_channels=128, out_channels=n_class, kernel_size=1, bias=True, stride=1)
        self.outconv2 = nn.Conv2d(in_channels=512, out_channels=128, kernel_size=3, bias=True, stride=1, padding=1)
        self.outconvs2 = nn.Conv2d(in_channels=128, out_channels=n_class, kernel_size=1, bias=True, stride=1)
        self.outconv3 = nn.Conv2d(in_channels=1024, out_channels=128, kernel_size=3, bias=True, stride=1, padding=1)
        self.outconvs3 = nn.Conv2d(in_channels=128, out_channels=n_class, kernel_size=1, bias=True, stride=1)
        self.outconv4 = nn.Conv2d(in_channels=2048, out_channels=128, kernel_size=3, bias=True, stride=1, padding=1)
        self.outconvs4 = nn.Conv2d(in_channels=128, out_channels=n_class, kernel_size=1, bias=True, stride=1)

        # init the decoder parameters
        if self.is_init_decoder:
            for m in self.modules():
                if isinstance(m, nn.Conv2d):
                    torch.nn.init.kaiming_normal_(m.weight, mode='fan_in')
                if isinstance(m, nn.BatchNorm2d):
                    m.weight.data.fill_(1)
            print('init the decoder parameters')

    def forward(self, feat, residual):

        x5, x4, x3, x2, x1 = feat[4], feat[3], feat[2], feat[1], feat[0]

        middle_output = []

        x5_up = F.interpolate(x5, size=x4.size()[2:], mode='nearest')
        x4_skip = self.shortbn4(self.shortconv4(x4))
        x5_cat = torch.cat((x4_skip, x5_up), dim=1)
        x5_out = self.uplayer4(x5_cat)
        x5_out1 = self.outconvs4(F.interpolate(self.outconv4(x5_out), size=residual.size()[2:], mode='nearest'))
        middle_output.append(x5_out1)

        x4_up = F.interpolate(x5_out, size=x3.size()[2:], mode='nearest')
        x3_skip = self.shortbn3(self.shortconv3(x3))
        x4_cat = torch.cat((x3_skip, x4_up), dim=1)
        x4_out = self.uplayer3(x4_cat)
        x4_out1 = self.outconvs3(F.interpolate(self.outconv3(x4_out), size=residual.size()[2:], mode='nearest'))
        middle_output.append(x4_out1)

        x3_up = F.interpolate(x4_out, size=x2.size()[2:], mode='nearest')
        x2_skip = self.shortbn2(self.shortconv2(x2))
        x3_cat = torch.cat((x2_skip, x3_up), dim=1)
        x3_out = self.uplayer2(x3_cat)
        x3_out1 = self.outconvs2(F.interpolate(self.outconv2(x3_out), size=residual.size()[2:], mode='nearest'))
        middle_output.append(x3_out1)

        x2_up = F.interpolate(x3_out, size=x1.size()[2:], mode='nearest')
        x1_skip = self.shortbn1(self.shortconv1(x1))
        x2_cat = torch.cat((x1_skip, x2_up), dim=1)
        x2_out = self.uplayer1(x2_cat)
        x2_out1 = self.outconvs1(F.interpolate(self.outconv1(x2_out), size=residual.size()[2:], mode='nearest'))
        middle_output.append(x2_out1)

        x1_up = F.interpolate(x2_out, size=residual.size()[2:], mode='nearest')
        x1_out = self.uplayer0(x1_up)
        out = self.outconvs0(self.outconv0(x1_out))

        middle_output.append(out)

        return middle_output