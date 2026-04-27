# SXT图像切分和预测结果合并工具使用说明

## 概述
这套工具包含两个主要脚本：
1. `split_sxt_images.py` - 将SXT图像切分成小块，生成COCO格式数据集
2. `merge_predictions.py` - 将切分图像的预测结果合并回原图

## 1. 图像切分 (split_sxt_images.py)

### 功能特点
- 自动从3D标签文件生成2D mask
- 支持数据集别名映射和翻转
- 生成COCO格式标注文件，支持后续反推
- 可配置切分网格、尺寸和重叠参数
- 生成切分映射文件用于反推

### 基本用法
```bash
# 处理所有数据（使用默认参数）
python split_sxt_images.py

# 测试模式（只处理第一个数据目录的第一张图）
python split_sxt_images.py --test

# 限制处理数据数量
python split_sxt_images.py -l 3

# 自定义切分参数
python split_sxt_images.py --xgrid 6 --ygrid 8 --width 120 --height 180

# 设置最小mask面积阈值
python split_sxt_images.py --min-mask-area 200
```

### 参数说明
- `-o, --output`: 输出目录 (默认: D:\Gitspace\ipa_full\iPA\data\sxt_images\split_img)
- `-l, --limit`: 限制处理的数据目录数量 (0=全部)
- `--test`: 测试模式，只处理第一个数据目录的第一张图
- `--width`: 切分宽度 (默认: 100)
- `--height`: 切分高度 (默认: 150)
- `--xgrid`: X方向切分数量 (默认: 4)
- `--ygrid`: Y方向切分数量 (默认: 5)
- `--overlap`: 重叠比例 (-1表示无重叠)
- `--min-mask-area`: 最小mask面积阈值 (默认: 100)

### 输出文件
```
output_dir/
├── images/                    # 切分后的图像
│   ├── 766_10/
│   ├── 766_11/
│   └── ...
├── masks/                     # 对应的mask
│   ├── 766_10/
│   ├── 766_11/
│   └── ...
├── split_instances.json       # COCO格式标注文件
├── split_mapping.json         # 切分映射信息（用于反推）
└── split_summary.json         # 处理摘要
```

## 2. 预测结果合并 (merge_predictions.py)

### 功能特点
- 将切分图像的预测结果合并回原始图像
- 支持COCO格式的预测结果
- 处理重叠区域的冲突
- 生成可视化的合并mask和COCO格式标注

### 基本用法
```bash
# 合并预测结果
python merge_predictions.py -m split_mapping.json -p predictions.json -o merged_results

# 如果原始图像路径发生变化，指定备用搜索目录
python merge_predictions.py -m split_mapping.json -p predictions.json -o merged_results --original-images-dir /path/to/images
```

### 参数说明
- `-m, --mapping`: split_mapping.json文件路径 (必需)
- `-p, --predictions`: 预测结果文件路径 (COCO格式, 必需)
- `-o, --output`: 输出目录 (必需)
- `--original-images-dir`: 原始图像的备用搜索目录

### 输出文件
```
output_dir/
├── merged_masks/              # 合并后的mask
│   ├── image1_merged_mask.png # 彩色可视化mask
│   ├── image1_raw_mask.png    # 原始mask（用于进一步处理）
│   └── ...
└── merged_predictions.json    # 合并后的COCO格式标注
```

## 3. 完整工作流程

### 步骤1: 切分图像
```bash
# 切分SXT图像
python split_sxt_images.py --test  # 先测试
python split_sxt_images.py         # 正式处理
```

### 步骤2: 训练模型
使用生成的`split_instances.json`和`images/`目录训练你的实例分割模型

### 步骤3: 运行预测
在切分后的图像上运行模型预测，生成COCO格式的预测结果

### 步骤4: 合并结果
```bash
# 将预测结果合并回原图
python merge_predictions.py \
    -m split_img/split_mapping.json \
    -p your_predictions.json \
    -o merged_results
```

## 4. 注意事项

### 数据要求
- SXT图像位于: `D:\Gitspace\ipa_full\iPA\data\sxt_images\mrcslice_output_denoise\{data_id}\slices data\`
- 3D标签位于: `D:\Gitspace\ipa_full\iPA\data\sxt_images\label_slice\{data_id}\combine_{data_id}_label.tif`
- 切片映射文件: `D:\Gitspace\ipa_full\iPA\data\sxt_images\label_slice\slice_mapping_3d_to_2d.json`
- 数据集信息文件: `D:\Gitspace\ipa_full\iPA\data\sxt_images\label_slice\dataset_info.json`

### 性能建议
- 使用`--test`参数先验证流程
- 根据GPU内存调整切分尺寸
- 合理设置`min_mask_area`避免产生过多小目标

### 故障排除
1. 如果找不到切片映射，检查文件名模式是否匹配
2. 如果mask尺寸不匹配，会自动调整尺寸
3. 如果合并时找不到原始图像，使用`--original-images-dir`参数

## 5. 与原始split_img.py的关系

`split_sxt_images.py`参考了`split_img.py`的核心算法：
- 相同的切分尺寸计算逻辑(`fix_wh`函数)
- 相同的网格切分算法
- 相同的重叠处理方式

主要改进：
- 适配SXT数据结构
- 自动从3D标签生成2D mask
- 完整的COCO格式支持
- 支持预测结果反推