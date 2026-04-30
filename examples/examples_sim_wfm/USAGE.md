# SIM/WFM Demo Scripts - Quick Reference

⚠️ **注意**: 详细文档请查看 [README.md](README.md)

本文件提供 SIM (Structured Illumination Microscopy) 和 WFM (Widefield Microscopy) 示例脚本的快速参考。

## 📋 脚本列表

### 分割与预测 (6个)
- `demo_SIM_isg_segmentation.py` - ISG 分割（真实数据）
- `demo_SIM_mito_segmentation.py` - Mito 分割（阈值法）
- `demo_SIM_Mito_predict.py` - Mito 预测（合成数据）
- `demo_SIM_ER_predict.py` - ER 预测（真实荧光数据）
- `demo_SIM_cell_segmentation.py` - Cell 分割（U-Net）
- `demo_SIM_visualization.py` - 可视化示例

### 空间分析 (4个)
- `demo_SIM_rdf_analysis.py` - RDF 径向分布函数
- `demo_SIM_partitioning.py` - SIM 空间分区
- `demo_WFM_partitioning.py` - WFM 空间分区
- `demo_SIM_WFM_interaction.py` - 细胞器相互作用

### 动力学分析 (1个)
- `demo_WFM_dynamics.py` - WFM 3D 速度分析（17,445 measurements）

### 训练脚本 (2个)
- `demo_SIM_cell_seg_train.py` - Cell 分割训练
- `demo_SIM_ER_train.py` - ER 分割训练

## 🚀 快速开始

```bash
# 激活环境
conda activate ipa_env

# 进入项目根目录
cd /media/cuixi/data7/liad/gitspace/iPA

# 运行任意脚本
python examples/examples_sim_wfm/demo_SIM_isg_segmentation.py
```

## 📊 数据来源

| 数据类型 | 位置 | 大小 |
|---------|------|------|
| SIM ISG | `data/sim_images/sim_images/` | 99 MB |
| WFM tracking | `data/wfm_images/` | NPY 文件 |
| ER 荧光 | `data/other_modalities/fluorescence/ER/` | PNG 文件 |

## 💡 关键特性

✅ **真实数据优先** - 所有预测和分析使用真实数据  
✅ **模块化设计** - 每个脚本独立运行  
✅ **统一日志** - QuickLogger 记录运行状态  
✅ **非交互式绘图** - 适合服务器环境  

## 📖 详细文档

完整的用法说明、参数配置和故障排除，请查看：
- [README.md](README.md) - 完整使用指南
- [FINAL_LIBRARY_CHECKLIST.md](../../documentations/FINAL_LIBRARY_CHECKLIST.md) - 库检查清单

---

**最后更新**: 2026-04-30





