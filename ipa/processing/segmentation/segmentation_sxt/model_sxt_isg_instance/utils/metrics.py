"""
训练和评估指标工具类
包含准确率、AUC、召回率等常用指标
"""

import numpy as np
import torch
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc


class AvgMetric:
    """平均指标计算器"""
    def __init__(self, name):
        self.name = name
        self.reset()
    
    def reset(self):
        """重置指标"""
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0
    
    def update(self, batch_idx, epoch, val, n=1):
        """更新指标值"""
        if val is None:
            return
        self.val = val
        self.sum += val * n
        self.count += n
        self.avg = self.sum / self.count if self.count > 0 else 0
    
    def get(self):
        """获取当前指标名称和值"""
        return self.name, self.avg
    
    def __str__(self):
        return f'{self.name}: {self.avg:.6f}'


class AccuracyMetric(AvgMetric):
    """准确率指标"""
    def __init__(self, name='accuracy'):
        super(AccuracyMetric, self).__init__(name)
    
    def update_with_logits(self, logits, targets, n=1):
        """使用logits和真实值更新准确率"""
        if logits is None or targets is None:
            return
        
        if isinstance(logits, torch.Tensor):
            predictions = torch.argmax(logits, dim=1)
            if isinstance(targets, torch.Tensor):
                correct = (predictions == targets).float().mean().item()
            else:
                correct = (predictions.cpu().numpy() == targets).mean()
        else:
            predictions = np.argmax(logits, axis=1)
            correct = (predictions == targets).mean()
        
        self.update(None, None, correct, n)


class LossMetric(AvgMetric):
    """损失值指标"""
    def __init__(self, name='loss'):
        super(LossMetric, self).__init__(name)
    
    def update_with_tensor(self, loss_tensor, n=1):
        """使用张量值更新损失"""
        if loss_tensor is None:
            return
        
        if isinstance(loss_tensor, torch.Tensor):
            loss_val = loss_tensor.item() if loss_tensor.numel() == 1 else loss_tensor.mean().item()
        else:
            loss_val = float(loss_tensor)
        
        self.update(None, None, loss_val, n)


class AUCMetric:
    """AUC指标"""
    def __init__(self, name='auc'):
        self.name = name
        self.reset()
    
    def reset(self):
        self.predictions = []
        self.targets = []
    
    def update(self, predictions, targets, mask=None):
        """更新AUC
        Args:
            predictions: 预测概率 (N, C) 或 (N,)
            targets: 真实标签 (N,)
            mask: 可选的掩码 (N,)
        """
        if torch.is_tensor(predictions):
            predictions = predictions.detach().cpu().numpy()
        if torch.is_tensor(targets):
            targets = targets.detach().cpu().numpy()
        if mask is not None and torch.is_tensor(mask):
            mask = mask.detach().cpu().numpy()
            
        # 如果是多类别，取正类概率
        if len(predictions.shape) > 1 and predictions.shape[1] > 1:
            predictions = predictions[:, 1]  # 假设第1列是正类概率
            
        if mask is not None:
            valid_idx = mask.astype(bool)
            predictions = predictions[valid_idx]
            targets = targets[valid_idx]
            
        self.predictions.extend(predictions.flatten())
        self.targets.extend(targets.flatten())
    
    def get(self):
        """计算AUC"""
        if len(set(self.targets)) < 2:
            return self.name, 0.5  # 只有一个类别时返回0.5
        try:
            auc_score = roc_auc_score(self.targets, self.predictions)
            return self.name, auc_score
        except:
            return self.name, 0.5


class RecallMetric:
    """召回率指标"""
    def __init__(self, name='recall'):
        self.name = name
        self.reset()
    
    def reset(self):
        self.true_positives = 0
        self.false_negatives = 0
    
    def update(self, predictions, targets, mask=None):
        """更新召回率
        Args:
            predictions: 预测值 (N, C) 或 (N,)
            targets: 真实标签 (N,)
            mask: 可选的掩码 (N,)
        """
        if torch.is_tensor(predictions):
            predictions = predictions.detach().cpu().numpy()
        if torch.is_tensor(targets):
            targets = targets.detach().cpu().numpy()
        if mask is not None and torch.is_tensor(mask):
            mask = mask.detach().cpu().numpy()
            
        if len(predictions.shape) > 1:
            predictions = np.argmax(predictions, axis=-1)
            
        if mask is not None:
            valid_idx = mask.astype(bool)
            predictions = predictions[valid_idx]
            targets = targets[valid_idx]
            
        self.true_positives += np.sum((predictions == 1) & (targets == 1))
        self.false_negatives += np.sum((predictions == 0) & (targets == 1))
    
    def get(self):
        """获取召回率"""
        recall = self.true_positives / max(self.true_positives + self.false_negatives, 1)
        return self.name, recall


class PrecisionMetric:
    """精确率指标"""
    def __init__(self, name='precision'):
        self.name = name
        self.reset()
    
    def reset(self):
        self.true_positives = 0
        self.false_positives = 0
    
    def update(self, predictions, targets, mask=None):
        """更新精确率"""
        if torch.is_tensor(predictions):
            predictions = predictions.detach().cpu().numpy()
        if torch.is_tensor(targets):
            targets = targets.detach().cpu().numpy()
        if mask is not None and torch.is_tensor(mask):
            mask = mask.detach().cpu().numpy()
            
        if len(predictions.shape) > 1:
            predictions = np.argmax(predictions, axis=-1)
            
        if mask is not None:
            valid_idx = mask.astype(bool)
            predictions = predictions[valid_idx]
            targets = targets[valid_idx]
            
        self.true_positives += np.sum((predictions == 1) & (targets == 1))
        self.false_positives += np.sum((predictions == 1) & (targets == 0))
    
    def get(self):
        """获取精确率"""
        precision = self.true_positives / max(self.true_positives + self.false_positives, 1)
        return self.name, precision


class F1Metric:
    """F1分数指标"""
    def __init__(self, name='f1'):
        self.name = name
        self.precision_metric = PrecisionMetric('precision')
        self.recall_metric = RecallMetric('recall')
    
    def reset(self):
        self.precision_metric.reset()
        self.recall_metric.reset()
    
    def update(self, predictions, targets, mask=None):
        """更新F1分数"""
        self.precision_metric.update(predictions, targets, mask)
        self.recall_metric.update(predictions, targets, mask)
    
    def get(self):
        """获取F1分数"""
        _, precision = self.precision_metric.get()
        _, recall = self.recall_metric.get()
        
        f1 = 2 * (precision * recall) / max(precision + recall, 1e-8)
        return self.name, f1


class MetricTracker:
    """指标跟踪器，管理多个指标"""
    def __init__(self):
        self.metrics = {}
    
    def add_metric(self, metric):
        """添加指标"""
        self.metrics[metric.name] = metric
    
    def update_all(self, predictions, targets, mask=None):
        """更新所有指标"""
        for metric in self.metrics.values():
            metric.update(predictions, targets, mask)
    
    def get_all(self):
        """获取所有指标结果"""
        results = {}
        for name, metric in self.metrics.items():
            _, value = metric.get()
            results[name] = value
        return results
    
    def reset_all(self):
        """重置所有指标"""
        for metric in self.metrics.values():
            metric.reset()