# 数据准备指南

## 需要准备的数据

### 1. 原始图像 (Raw Images)
**存放位置**: `d:\Gitspace\ipa_full\data\raw_images\train\`

**数据要求**:
- 格式: .png 或 .jpg
- 类型: 灰度图像
- 内容: 胰岛素囊泡的显微镜图像
- 命名: 建议使用标准命名如 `Stevens_pancreatic_INS_1E_766_8_204.png`

### 2. 标注图像 (Labels)  
**存放位置**: `d:\Gitspace\ipa_full\data\labels\train\`

**数据要求**:
- 格式: .png (灰度图像)
- 内容: 分割mask，非零像素表示感兴趣区域
- 命名: 必须与原始图像完全一致
- 用途: 用于确定crop的边界

## 目录结构

请按以下结构准备数据:

```
d:\Gitspace\ipa_full\
├── data\
│   ├── raw_images\
│   │   └── train\
│   │       ├── Stevens_pancreatic_INS_1E_766_8_204.png
│   │       ├── Stevens_pancreatic_INS_1E_784_5_204.png
│   │       └── ... (更多原始图像)
│   └── labels\
│       └── train\
│           ├── Stevens_pancreatic_INS_1E_766_8_204.png
│           ├── Stevens_pancreatic_INS_1E_784_5_204.png
│           └── ... (对应的标注图像)
└── iPA\
    └── examples\
        └── test_isg model\
            └── instance\
                └── data\
                    └── crop_images\
                        └── train\  (输出目录)
```

## 数据获取建议

### 如果有现有数据集:
1. 将原始图像复制到 `raw_images/train/`
2. 将标注mask复制到 `labels/train/`

### 如果需要创建标注:
1. 使用 LabelMe 或 CVAT 创建分割mask
2. 确保mask是灰度图像，感兴趣区域为非零值

### 数据验证:
```bash
# 检查数据完整性
python -c "
import os
raw_dir = r'd:\Gitspace\ipa_full\data\raw_images\train'
label_dir = r'd:\Gitspace\ipa_full\data\labels\train'
raw_files = set(os.listdir(raw_dir)) if os.path.exists(raw_dir) else set()
label_files = set(os.listdir(label_dir)) if os.path.exists(label_dir) else set()
print(f'Raw images: {len(raw_files)}')
print(f'Label images: {len(label_files)}')
print(f'Missing labels: {raw_files - label_files}')
print(f'Missing raw images: {label_files - raw_files}')
"
```

## 运行Crop步骤

准备好数据后，运行:
```bash
cd d:\Gitspace\ipa_full\iPA\examples\test_isg model\instance\crop\
python crop.py
```

这将:
1. 读取raw_images和labels
2. 基于label确定crop边界
3. 输出cropped图像到 `instance\data\crop_images\train\`
4. 生成crop信息到 `crop_info_train.json`