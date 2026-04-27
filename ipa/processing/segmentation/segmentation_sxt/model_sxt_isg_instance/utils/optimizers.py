"""
自定义优化器
支持不同参数组的不同学习率
"""

import torch
import torch.nn as nn
from torch.optim import Optimizer
import math


class SGD(Optimizer):
    """Custom SGD optimizer with learning rate scaling support."""
    
    def __init__(self, params, lr=1e-3, momentum=0, dampening=0,
                 weight_decay=0, nesterov=False):
        if not 0.0 <= lr:
            raise ValueError("Invalid learning rate: {}".format(lr))
        if momentum < 0.0:
            raise ValueError("Invalid momentum value: {}".format(momentum))
        if weight_decay < 0.0:
            raise ValueError("Invalid weight_decay value: {}".format(weight_decay))

        defaults = dict(lr=lr, momentum=momentum, dampening=dampening,
                       weight_decay=weight_decay, nesterov=nesterov)
        if nesterov and (momentum <= 0 or dampening != 0):
            raise ValueError("Nesterov momentum requires a momentum and zero dampening")
        
        # Handle parameter groups with different learning rates
        if not isinstance(params, list):
            params = list(params)
        
        # Convert to parameter groups if needed
        param_groups = []
        for param in params:
            if isinstance(param, dict):
                param_groups.append(param)
            else:
                param_groups.append({'params': param})
        
        super(SGD, self).__init__(param_groups, defaults)

    def __setstate__(self, state):
        super(SGD, self).__setstate__(state)
        for group in self.param_groups:
            group.setdefault('nesterov', False)

    def step_with_lr(self, lr_scale):
        """Step with learning rate scaling."""
        for group in self.param_groups:
            # Apply learning rate scaling
            scaled_lr = group['lr'] * lr_scale
            self._step_group(group, scaled_lr)

    def step(self, closure=None):
        """Standard step function."""
        loss = None
        if closure is not None:
            loss = closure()

        for group in self.param_groups:
            self._step_group(group, group['lr'])

        return loss

    def _step_group(self, group, lr):
        """Apply SGD update to a parameter group."""
        weight_decay = group['weight_decay']
        momentum = group['momentum']
        dampening = group['dampening']
        nesterov = group['nesterov']

        for p in group['params']:
            if p.grad is None:
                continue
            d_p = p.grad.data
            if weight_decay != 0:
                d_p = d_p.add(p.data, alpha=weight_decay)
            if momentum != 0:
                param_state = self.state[p]
                if 'momentum_buffer' not in param_state:
                    buf = param_state['momentum_buffer'] = torch.zeros_like(p.data)
                    buf.mul_(momentum).add_(d_p)
                else:
                    buf = param_state['momentum_buffer']
                    buf.mul_(momentum).add_(d_p, alpha=1 - dampening)
                if nesterov:
                    d_p = d_p.add(buf, alpha=momentum)
                else:
                    d_p = buf

            p.data.add_(d_p, alpha=-lr)


class Adam(Optimizer):
    """Custom Adam optimizer with learning rate scaling support."""
    
    def __init__(self, params, lr=1e-3, betas=(0.9, 0.999), eps=1e-8,
                 weight_decay=0, amsgrad=False):
        if not 0.0 <= lr:
            raise ValueError("Invalid learning rate: {}".format(lr))
        if not 0.0 <= eps:
            raise ValueError("Invalid epsilon value: {}".format(eps))
        if not 0.0 <= betas[0] < 1.0:
            raise ValueError("Invalid beta parameter at index 0: {}".format(betas[0]))
        if not 0.0 <= betas[1] < 1.0:
            raise ValueError("Invalid beta parameter at index 1: {}".format(betas[1]))
        if not 0.0 <= weight_decay:
            raise ValueError("Invalid weight_decay value: {}".format(weight_decay))

        defaults = dict(lr=lr, betas=betas, eps=eps,
                       weight_decay=weight_decay, amsgrad=amsgrad)
        
        # Handle parameter groups
        if not isinstance(params, list):
            params = list(params)
        
        param_groups = []
        for param in params:
            if isinstance(param, dict):
                param_groups.append(param)
            else:
                param_groups.append({'params': param})
        
        super(Adam, self).__init__(param_groups, defaults)

    def __setstate__(self, state):
        super(Adam, self).__setstate__(state)
        for group in self.param_groups:
            group.setdefault('amsgrad', False)

    def step_with_lr(self, lr_scale):
        """Step with learning rate scaling."""
        for group in self.param_groups:
            scaled_lr = group['lr'] * lr_scale
            self._step_group(group, scaled_lr)

    def step(self, closure=None):
        """Standard step function."""
        loss = None
        if closure is not None:
            loss = closure()

        for group in self.param_groups:
            self._step_group(group, group['lr'])

        return loss

    def _step_group(self, group, lr):
        """Apply Adam update to a parameter group."""
        for p in group['params']:
            if p.grad is None:
                continue
            grad = p.grad.data
            if grad.dtype in {torch.float16, torch.bfloat16}:
                grad = grad.float()

            amsgrad = group['amsgrad']

            state = self.state[p]

            # State initialization
            if len(state) == 0:
                state['step'] = 0
                state['exp_avg'] = torch.zeros_like(p.data).float()
                state['exp_avg_sq'] = torch.zeros_like(p.data).float()
                if amsgrad:
                    state['max_exp_avg_sq'] = torch.zeros_like(p.data).float()

            exp_avg, exp_avg_sq = state['exp_avg'], state['exp_avg_sq']
            if amsgrad:
                max_exp_avg_sq = state['max_exp_avg_sq']
            beta1, beta2 = group['betas']

            state['step'] += 1
            bias_correction1 = 1 - beta1 ** state['step']
            bias_correction2 = 1 - beta2 ** state['step']

            if group['weight_decay'] != 0:
                grad = grad.add(p.data, alpha=group['weight_decay'])

            # Exponential moving average of gradient values
            exp_avg.mul_(beta1).add_(grad, alpha=1 - beta1)

            # Exponential moving average of squared gradient values
            exp_avg_sq.mul_(beta2).addcmul_(grad, grad, value=1 - beta2)
            if amsgrad:
                # Maintains the maximum of all 2nd moment running averages
                torch.maximum(max_exp_avg_sq, exp_avg_sq, out=max_exp_avg_sq)
                # Use the max. for normalizing running avg. of gradient
                denom = (max_exp_avg_sq.sqrt() / math.sqrt(bias_correction2)).add_(group['eps'])
            else:
                denom = (exp_avg_sq.sqrt() / math.sqrt(bias_correction2)).add_(group['eps'])

            step_size = lr / bias_correction1

            p.data.addcdiv_(exp_avg, denom, value=-step_size)


def clip_grad(parameters, max_norm, norm_type=2):
    """Clip gradients by norm."""
    if isinstance(parameters, torch.Tensor):
        parameters = [parameters]
    parameters = list(filter(lambda p: p.grad is not None, parameters))
    max_norm = float(max_norm)
    norm_type = float(norm_type)
    
    if norm_type == float('inf'):
        total_norm = max(p.grad.data.abs().max() for p in parameters)
    else:
        total_norm = 0
        for p in parameters:
            param_norm = p.grad.data.norm(norm_type)
            total_norm += param_norm.item() ** norm_type
        total_norm = total_norm ** (1. / norm_type)
    
    clip_coef = max_norm / (total_norm + 1e-6)
    if clip_coef < 1:
        for p in parameters:
            p.grad.data.mul_(clip_coef)
    
    return total_norm