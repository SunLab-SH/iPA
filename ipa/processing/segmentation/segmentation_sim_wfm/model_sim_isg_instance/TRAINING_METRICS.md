# 训练和评估指标示例

## 新增功能

### 1. 完整的评估指标

现在训练脚本会在每个epoch打印以下指标：

**训练指标**:
- accuracy: 准确率
- precision: 精确率  
- recall: 召回率
- f1: F1分数
- auc: ROC曲线下面积
- ap: 平均精确率 (Average Precision)

**验证指标**:
- 相同的指标，前缀为`Val_`

### 2. 数据使用控制

- **data_ratio**: 控制使用数据的比例 (0.0-1.0)
- **train_val_split**: 控制训练/验证数据分割比例，默认0.8 (80%训练，20%验证)

### 3. 日志记录

- **文件日志**: 自动保存到 `{output_dir}/train.log`
- **TensorBoard**: 保存到 `{output_dir}/tensorboard_logs/`
- 移除了WandB，简化依赖

### 4. 使用示例

#### 基本训练（使用全部数据，80%训练20%验证）
```bash
python torchvision_maskrcnn_train.py \
    --epochs 20 \
    --batch_size 4 \
    --print_freq 10
```

#### 使用部分数据训练（比如只用50%的数据）
```bash
python torchvision_maskrcnn_train.py \
    --epochs 20 \
    --data_ratio 0.5 \
    --train_val_split 0.8
```

#### 调整训练/验证分割比例（90%训练10%验证）
```bash
python torchvision_maskrcnn_train.py \
    --epochs 20 \
    --train_val_split 0.9 \
    --data_ratio 1.0
```

#### 查看TensorBoard日志
```bash
tensorboard --logdir D:\Gitspace\ipa_full\iPA\data\sxt_images\models\tensorboard_logs
```

### 5. 训练输出示例

```
使用数据比例: 0.50, 使用 6906/13812 个样本
数据分割比例: 0.80
训练集大小: 5524
验证集大小: 1382

=== Epoch 0 Training Metrics ===
accuracy: 0.8234
precision: 0.7891
recall: 0.8567
f1: 0.8215
auc: 0.8934
ap: 0.8456

=== Epoch 0 Validation Metrics ===
Val_accuracy: 0.8123
Val_precision: 0.7734
Val_recall: 0.8345
Val_f1: 0.8023
Val_auc: 0.8723
Val_ap: 0.8234
```

### 6. 参数说明

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--data_ratio` | 1.0 | 使用数据的比例，1.0=100%，0.5=50% |
| `--train_val_split` | 0.8 | 训练/验证分割比例，0.8=80%训练20%验证 |
| `--epochs` | 10 | 训练轮数 |
| `--batch_size` | 2 | 批次大小 |
| `--lr` | 0.005 | 学习率 |

### 7. 最佳模型保存

训练会自动保存性能最好的模型到 `best_model.pth`，基于AP指标。

### 8. 安装依赖

```bash
pip install scikit-learn tensorboard matplotlib
```

### 9. 恢复训练

```bash
python torchvision_maskrcnn_train.py \
    --resume models/model_10.pth \
    --start_epoch 11 \
    --data_ratio 0.8 \
    --train_val_split 0.85
```

所有指标都会自动记录和可视化，现在可以灵活控制数据使用量和分割比例！