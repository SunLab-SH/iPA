# SIM/WFM Examples 使用指南

本目录包含 Structured Illumination Microscopy (SIM) 和 Widefield Microscopy (WFM) 数据的分析示例脚本。

## 📁 目录结构

```
examples/examples_sim_wfm/
├── USAGE.md                          # 本文件
├── demo_SIM_isg_segmentation.py      # ISG 分割示例
├── demo_SIM_mito_segmentation.py     # Mito 分割示例
├── demo_SIM_Mito_predict.py          # Mito 预测（合成数据）
├── demo_SIM_ER_predict.py            # ER 预测（真实数据）
├── demo_SIM_cell_segmentation.py     # Cell 分割示例
├── demo_SIM_rdf_analysis.py          # RDF 分析示例
├── demo_SIM_visualization.py         # 可视化示例
├── demo_SIM_partitioning.py          # SIM 空间分区
├── demo_SIM_WFM_interaction.py       # SIM/WFM 相互作用分析
├── demo_WFM_dynamics.py              # WFM 动力学分析
├── demo_WFM_partitioning.py          # WFM 空间分区
├── demo_SIM_cell_seg_train.py        # Cell 分割训练
└── demo_SIM_ER_train.py              # ER 分割训练
```

## 📊 脚本分类

### 1. 分割与预测脚本 (6个)

#### ✅ 使用真实数据

| 脚本 | 功能 | 数据来源 | 说明 |
|------|------|---------|------|
| `demo_SIM_isg_segmentation.py` | ISG 实例分割 | `data/sim_images/sim_images/20220909_30-2-1-SIM_raw_ISG.tif` (99 MB) | 阈值法 + 形态学后处理 |
| `demo_SIM_ER_predict.py` | ER 网络分割 | `data/other_modalities/fluorescence/ER/*.png` | ERNet 模型预测，COS7 细胞 ER-YFP 荧光数据 |
| `demo_SIM_cell_segmentation.py` | 细胞分割 | U-Net 模型 | 细胞边界检测 |
| `demo_SIM_mito_segmentation.py` | 线粒体分割 | 阈值法 | 线粒体形态学分析 |

#### ⚠️ 使用合成数据

| 脚本 | 功能 | 说明 |
|------|------|------|
| `demo_SIM_Mito_predict.py` | Mito 预测演示 | 无真实 SIM Mito 数据，使用合成数据合理 |

### 2. 空间分析脚本 (4个)

| 脚本 | 功能 | 模块 | 数据来源 |
|------|------|------|---------|
| `demo_SIM_rdf_analysis.py` | 径向分布函数 (RDF) | arrangement | ISG 从 NE 到 PM 的分布 |
| `demo_SIM_partitioning.py` | SIM 空间分区分析 | partitioning | `data/wfm_images/*_PM.mrc`, `*_N.mrc` (SIM masks) |
| `demo_WFM_partitioning.py` | WFM 空间分区 | partitioning | WFM 时间序列数据的分区 |
| `demo_SIM_WFM_interaction.py` | 细胞器相互作用 | interaction | 距离分析、接触概率 |

### 3. 动力学分析脚本 (1个)

| 脚本 | 功能 | 数据来源 | 结果 |
|------|------|---------|------|
| `demo_WFM_dynamics.py` | 3D 速度分析 | `data/wfm_images/16.7-30_P21_Velocity_3dvelocity.npy` | 17,445 个速度测量值 |

**注意**: WFM 第一个轴是时间维度，通过前后帧差分计算速度。

### 4. 可视化脚本 (1个)

| 脚本 | 功能 | 说明 |
|------|------|------|
| `demo_SIM_visualization.py` | SIM 数据可视化 | 多通道叠加、3D 渲染 |

### 5. 训练脚本 (2个)

| 脚本 | 功能 | 说明 |
|------|------|------|
| `demo_SIM_cell_seg_train.py` | Cell 分割模型训练 | U-Net 训练演示 |
| `demo_SIM_ER_train.py` | ER 分割模型训练 | ERNet 训练演示 |

## 🔧 运行示例

### 基本运行

```bash
# 激活环境
conda activate ipa_env

# 进入项目根目录
cd /media/cuixi/data7/liad/gitspace/iPA

# 运行任意脚本
python examples/examples_sim_wfm/demo_SIM_isg_segmentation.py
```

### 后台运行（长时间任务）

```bash
# 后台运行并记录日志
nohup python examples/examples_sim_wfm/demo_WFM_dynamics.py > /tmp/wfm_dynamics.log 2>&1 &

# 查看进度
tail -f /tmp/wfm_dynamics.log
```

## 📈 典型输出

### 分割结果
```
results/sim_isg_segmentation_demo/
├── isg_mask.npy              # ISG 分割掩码
├── input_image.npy           # 输入图像
└── segmentation_result.png   # 可视化结果
```

### 分析结果
```
results/sim_rdf_analysis_demo/
├── rdf_results.json          # RDF 计算结果
└── rdf_plot.png              # RDF 分布图
```

### 动力学结果
```
results/wfm_dynamics_demo/
├── velocity_distribution.png # 速度分布图
├── velocity_statistics.json  # 速度统计
└── trajectory_plots.png      # 轨迹图
```

## 💡 关键特性

### 1. 真实数据优先

所有预测和分析脚本都优先使用真实数据：
- **SIM ISG**: 99 MB TIF 文件，真实胰岛素分泌颗粒
- **WFM Dynamics**: 17,445 个真实速度测量值
- **ER Prediction**: COS7 细胞 ER-YFP 荧光图像

### 2. 模块化设计

每个脚本独立运行，依赖清晰的输入输出：
```python
# 标准路径配置
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
sys.path.insert(0, PROJECT_ROOT)
```

### 3. 统一日志系统

使用 `QuickLogger` 记录运行状态：
```python
from ipa.data_loader import QuickLogger
logger = QuickLogger("sim_isg_segmentation", log_dir=f'{PROJECT_ROOT}/logs')
logger.step("Starting analysis...")
logger.file_out(output_path)
```

### 4. 非交互式绘图

所有脚本使用 `matplotlib.use('Agg')` 后端，适合服务器环境：
```python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.close()
```

## 📝 数据使用说明

### SIM 数据位置

```
data/sim/
├── 20220909_30-2-1-SIM_raw_ISG.tif      (99 MB)   # ISG 原始图像
├── 20220909_30-2-1-SIM_ISG_seg.tif      (1.5 MB)  # ISG 分割结果
├── 20220909_30-2-1-SIM_raw_Actin.tif    (171 MB)  # Actin 原始图像
├── 20220909_30-2-1-SIM_PM.mrc                       # SIM PM mask
├── 20220909_30-2-1-SIM_N.mrc                        # SIM NE mask
├── 20220909_30-2-2-SIM_PM.mrc                       # SIM PM mask (partitioning)
├── 20220909_30-2-2-SIM_N.mrc                        # SIM NE mask (partitioning)
└── 20220909_30-2-2-SIM_ISG_seg.tif                  # SIM ISG seg (partitioning)
```

### WFM 数据位置

```
data/wfm/
├── 16.7-30_P21_Velocity_3dvelocity.npy  # 速度跟踪数据
├── 16.7-30_P21_Memb.mrc                 # WFM membrane mask
├── 16.7-30_P21_Nuc.mrc                  # WFM nucleus mask
├── 2.8-30_P3S3-1_Memb.mrc               # WFM membrane mask
└── 2.8-30_P3S3-1_Nuc.mrc                # WFM nucleus mask
```

### ER 荧光数据位置

```
data/other/fluorescence/ER/
├── *.png                                 # COS7 细胞 ER-YFP 图像
└── README.md                             # 数据说明
```

## ⚠️ 注意事项

### 1. 内存管理

部分脚本处理大型 3D 数据，确保有足够内存：
- SIM 原始数据：最高 328 MB
- WFM 时间序列：可能达到 GB 级别

### 2. GPU 加速

Partitioning 模块支持 GPU 加速：
```python
# 自动检测并使用 GPU（如果可用）
from ipa.processing.partitioning import Partitioning
partitioner = Partitioning(root_dir="results/", n_slices=8)
```

### 3. 训练脚本

训练脚本仅提供 API 演示，需要自备标注数据：
```python
# 训练脚本需要用户自己的数据集
# demo_SIM_cell_seg_train.py
# demo_SIM_ER_train.py
```

## 🔗 相关文档

- [FINAL_LIBRARY_CHECKLIST.md](../../documentations/FINAL_LIBRARY_CHECKLIST.md) - 完整库检查清单
- [ANALYSIS_MODULE.md](../../documentations/ANALYSIS_MODULE.md) - 分析模块文档
- [SEGMENTATION_MODULE.md](../../documentations/SEGMENTATION_MODULE.md) - 分割模块文档

## 📞 问题反馈

如遇到问题，请检查：
1. 数据文件是否存在且路径正确
2. 依赖包是否安装完整 (`requirements.txt`)
3. 日志文件中的错误信息 (`logs/` 目录)

---

**最后更新**: 2026-04-30  
**维护者**: iPA 开发团队
