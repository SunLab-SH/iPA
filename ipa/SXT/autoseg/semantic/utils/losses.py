#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2019/10/27 下午8:31
# @Author  : chuyu zhang
# @File    : losses.py
# @Software: PyCharm

import torch.nn.functional as F
import torch.nn as nn
from torch.autograd import Variable
import torch


class FocalLoss(nn.Module):
    def __init__(self, class_num=5, alpha=None, gamma=2, size_average=True):
        super(FocalLoss, self).__init__()
        if alpha is None:
            self.alpha = Variable(torch.ones(class_num, 1))
        else:
            if isinstance(alpha, Variable):
                self.alpha = alpha
            else:
                self.alpha = Variable(alpha)
        self.gamma = gamma
        self.class_num = class_num
        self.size_average = size_average

    def forward(self, inputs, targets):
        """
        calculate the focalloss
        :param inputs: N*C*H*W
        :param targets: N*H*W
        :return:
        """
        N = inputs.size(0)
        C = inputs.size(1)
        # For numerical stability using log_softmax
        P = F.log_softmax(inputs, dim=1)
        if P.dim() > 2:
            P = P.view(P.size(0), P.size(1), -1)
            P = P.permute(0, 2, 1).contiguous()
            P = P.view(-1, P.size(-1))

        ids = targets.view(-1, 1)

        class_mask = inputs.data.new(ids.size(0), C).fill_(0)
        class_mask = Variable(class_mask)
        class_mask = class_mask.scatter_(1, ids.data, 1.)

        if class_mask.device != P.device:
            class_mask = class_mask.to(P.device)

        if inputs.is_cuda and not self.alpha.is_cuda:
            self.alpha = self.alpha.cuda()
        alpha = self.alpha[ids.data.view(-1)]

        probs = (P * class_mask).sum(axis=1).view(-1, 1)

        # log_p = probs.log()

        batch_loss = -alpha * (torch.pow((1 - probs.exp()), self.gamma)) * probs

        if self.size_average:
            loss = batch_loss.mean()
        else:
            loss = batch_loss.sum()
        return loss


class Entropy(nn.Module):

    def __init__(self, size_average=True, alpha=0.01):
        """
        This loss function is designed for unlabel data.
        :param size_average: if true, mean loss
        :param alpha: balance the label loss and unlabel loss, set it 0.01
        """
        super(Entropy, self).__init__()
        self.size_average = size_average
        self.alpha = alpha

    def forward(self, pred):
        probs = F.log_softmax(pred, dim=1)
        if probs.dim() > 2:
            probs = probs.view(probs.size(0), probs.size(1), -1)
            probs = probs.permute(0, 2, 1).contiguous()
            probs = probs.view(-1, probs.size(-1))

        batch_loss = -self.alpha * probs * probs.exp()
        batch_loss = batch_loss.sum(axis=1)

        if self.size_average:
            loss = batch_loss.mean()
        else:
            loss = batch_loss.sum()
        return loss


def softmax_mse_loss(input_logits, target_logits):
    """Takes softmax on both sides and returns MSE loss

    Note:
    - Returns the sum over all examples. Divide by the batch size afterwards
      if you want the mean.
    - Sends gradients to inputs but not the targets.
    """
    assert input_logits.size() == target_logits.size()
    input_softmax = F.softmax(input_logits, dim=1)
    target_softmax = F.softmax(target_logits, dim=1)
    mse_loss = (input_softmax-target_softmax)**2
    return mse_loss


def softmax_kl_loss(input_logits, target_logits):
    """Takes softmax on both sides and returns KL divergence

    Note:
    - Returns the sum over all examples. Divide by the batch size afterwards
      if you want the mean.
    - Sends gradients to inputs but not the targets.
    """
    assert input_logits.size() == target_logits.size()
    input_log_softmax = F.log_softmax(input_logits, dim=1)
    target_softmax = F.softmax(target_logits, dim=1)

    # return F.kl_div(input_log_softmax, target_softmax)
    kl_div = F.kl_div(input_log_softmax, target_softmax, reduction='none')
    # mean_kl_div = torch.mean(0.2*kl_div[:,0,...]+0.8*kl_div[:,1,...])
    return kl_div


if __name__ == "__main__":
    pass
