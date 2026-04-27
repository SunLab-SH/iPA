"""
训练回调函数
包含速度计、进度监控等功能
"""

import time
import logging


class Speedometer:
    """训练速度监控器"""
    def __init__(self, batch_size, frequent=100):
        self.batch_size = batch_size
        self.frequent = frequent
        self.init = False
        self.tic = 0
        self.last_count = 0
    
    def __call__(self, param_iter, metrics=None):
        """记录训练速度和指标"""
        count = param_iter
        if self.last_count > count:
            self.init = False
        self.last_count = count
        
        if self.init:
            if count % self.frequent == 0:
                speed = self.frequent * self.batch_size / (time.time() - self.tic)
                if metrics is not None:
                    metrics_str = ', '.join([f'{m.get()[0]}={m.get()[1]:.4f}' for m in metrics])
                    print(f'Iter[{count}] Speed: {speed:.2f} samples/sec, {metrics_str}')
                else:
                    print(f'Iter[{count}] Speed: {speed:.2f} samples/sec')
                self.tic = time.time()
        else:
            self.init = True
            self.tic = time.time()


class LoggingCallback:
    """日志记录回调"""
    def __init__(self, logger, frequent=100):
        self.logger = logger
        self.frequent = frequent
    
    def __call__(self, param_iter, metrics=None):
        """记录指标到日志"""
        if param_iter % self.frequent == 0 and metrics is not None:
            for metric in metrics:
                name, value = metric.get()
                self.logger.info(f'Iter[{param_iter}] {name}: {value:.6f}')


class EarlyStoppingCallback:
    """基于验证损失的提前停止回调"""
    def __init__(self, patience=10, min_delta=0.001):
        self.patience = patience
        self.min_delta = min_delta
        self.best_loss = float('inf')
        self.wait = 0
        self.stopped_epoch = 0
    
    def __call__(self, epoch, val_loss):
        """检查是否应该停止训练"""
        if val_loss < self.best_loss - self.min_delta:
            self.best_loss = val_loss
            self.wait = 0
        else:
            self.wait += 1
            if self.wait >= self.patience:
                self.stopped_epoch = epoch
                return True
        return False


class LearningRateScheduler:
    """学习率调度器回调"""
    def __init__(self, optimizer, schedule='step', **kwargs):
        self.optimizer = optimizer
        self.schedule = schedule
        self.kwargs = kwargs
        self.step_count = 0
    
    def __call__(self, epoch=None, iter_count=None):
        """更新学习率"""
        if self.schedule == 'step':
            self._step_schedule(epoch)
        elif self.schedule == 'poly':
            self._poly_schedule(iter_count)
        elif self.schedule == 'cosine':
            self._cosine_schedule(epoch)
    
    def _step_schedule(self, epoch):
        """阶梯学习率调度"""
        step_size = self.kwargs.get('step_size', 30)
        gamma = self.kwargs.get('gamma', 0.1)
        
        if epoch % step_size == 0 and epoch > 0:
            for param_group in self.optimizer.param_groups:
                param_group['lr'] *= gamma
    
    def _poly_schedule(self, iter_count):
        """多项式学习率调度"""
        base_lr = self.kwargs.get('base_lr', 0.001)
        max_iter = self.kwargs.get('max_iter', 10000)
        power = self.kwargs.get('power', 0.9)
        
        if iter_count is not None:
            lr = base_lr * ((1 - float(iter_count) / max_iter) ** power)
            for param_group in self.optimizer.param_groups:
                param_group['lr'] = lr
    
    def _cosine_schedule(self, epoch):
        """余弦退火学习率调度"""
        import math
        base_lr = self.kwargs.get('base_lr', 0.001)
        max_epochs = self.kwargs.get('max_epochs', 100)
        
        if epoch is not None:
            lr = base_lr * 0.5 * (1 + math.cos(math.pi * epoch / max_epochs))
            for param_group in self.optimizer.param_groups:
                param_group['lr'] = lr