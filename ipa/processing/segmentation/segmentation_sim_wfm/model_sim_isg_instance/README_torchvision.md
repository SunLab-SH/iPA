# Torchvision Mask R-CNN for SXT Images

基于torchvision的Mask R-CNN模型，用于SXT切分图像的实例分割任务。

## 文件说明

- `torchvision_maskrcnn_train.py` - 训练脚本
- `torchvision_maskrcnn_predict.py` - 预测脚本
- 数据来源：`split_sxt_images.py`生成的COCO格式数据

## 快速开始

### 1. 训练模型

```bash
# 基本训练
python torchvision_maskrcnn_train.py

# 自定义参数
python torchvision_maskrcnn_train.py \
    --data_root D:\Gitspace\ipa_full\iPA\data\sxt_images\split_img \
    --coco_file D:\Gitspace\ipa_full\iPA\data\sxt_images\split_img\split_instances.json \
    --output_dir D:\Gitspace\ipa_full\iPA\data\sxt_images\models \
    --epochs 20 \
    --batch_size 4 \
    --lr 0.005
```

### 2. 预测

```bash
# 使用预训练模型预测
python torchvision_maskrcnn_predict.py --visualize --save_coco

# 使用训练好的模型预测
python torchvision_maskrcnn_predict.py \
    --model_path D:\Gitspace\ipa_full\iPA\data\sxt_images\models\model_weights_19.pth \
    --visualize \
    --save_coco \
    --num_samples 20 \
    --threshold 0.7
```

## 训练参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--data_root` | split_img目录 | 数据根目录 |
| `--coco_file` | split_instances.json | COCO格式标注文件 |
| `--output_dir` | models目录 | 模型输出目录 |
| `--epochs` | 10 | 训练轮数 |
| `--batch_size` | 2 | 批次大小 |
| `--lr` | 0.005 | 学习率 |
| `--momentum` | 0.9 | SGD动量 |
| `--weight_decay` | 0.0005 | 权重衰减 |

## 预测参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--model_path` | None | 训练好的模型路径 |
| `--threshold` | 0.5 | 置信度阈值 |
| `--num_samples` | 10 | 预测样本数量 |
| `--visualize` | False | 是否可视化结果 |
| `--save_coco` | False | 是否保存COCO格式结果 |

## 数据格式

### 输入数据结构
```
split_img/
├── images/
│   ├── 766_10/
│   │   ├── slice_766_10_Ins_0_split_000.png
│   │   ├── slice_766_10_Ins_0_split_001.png
│   │   └── ...
│   ├── 766_11/
│   └── ...
├── masks/
│   ├── 766_10/
│   │   ├── slice_766_10_Ins_0_split_000.png
│   │   └── ...
│   └── ...
├── split_instances.json    # COCO格式标注
├── split_mapping.json      # 切分映射信息
└── split_summary.json      # 处理摘要
```

### COCO格式标注
```json
{
  "images": [
    {
      "id": 1,
      "file_name": "slice_766_10_Ins_0_split_000.png",
      "width": 100,
      "height": 150,
      "image_name": "766_10_1"
    }
  ],
  "annotations": [
    {
      "id": 1,
      "image_id": 1,
      "category_id": 1,
      "segmentation": {"counts": "...", "size": [150, 100]},
      "area": 1234,
      "bbox": [x, y, width, height],
      "iscrowd": 0
    }
  ],
  "categories": [
    {
      "id": 1,
      "name": "granule",
      "supercategory": "object"
    }
  ]
}
```

## 输出结果

### 训练输出
```
models/
├── model_0.pth           # 完整检查点
├── model_weights_0.pth   # 仅模型权重
├── model_5.pth
├── model_weights_5.pth
└── ...
```

### 预测输出
```
predictions/
├── prediction_000_slice_766_10_Ins_0_split_000.png  # 可视化结果
├── prediction_001_slice_766_11_Ins_5_split_003.png
├── ...
└── predictions_coco_format.json                     # COCO格式预测结果
```

## 特性

1. **基于torchvision**: 使用官方实现的Mask R-CNN
2. **COCO格式支持**: 完全兼容COCO数据格式
3. **预训练模型**: 基于MS COCO预训练的权重
4. **可视化**: 自动生成预测结果可视化
5. **灵活配置**: 支持多种训练和预测参数

## 模型架构

- **Backbone**: ResNet-50 + FPN
- **RPN**: Region Proposal Network
- **分类器**: Fast R-CNN头部
- **分割器**: Mask预测头部
- **类别数**: 2 (背景 + granule)

## 性能优化建议

1. **批次大小**: 根据GPU内存调整batch_size
2. **学习率**: 建议从0.005开始，可根据收敛情况调整
3. **数据增强**: 可在训练时添加更多数据增强
4. **阈值调整**: 根据验证结果调整预测阈值

## 与原始切分数据的兼容性

此脚本完全兼容`split_sxt_images.py`生成的数据格式，可以直接使用生成的COCO格式标注文件进行训练和预测。预测结果也可以通过`merge_predictions.py`脚本合并回原始图像。