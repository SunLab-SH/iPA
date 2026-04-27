# SIM Demo Scripts - Usage Guide

Quick reference for using the standardized SIM (Structured Illumination Microscopy) demo scripts.

## Scripts Overview

| Script | Purpose | Primary Use Case |
|--------|---------|------------------|
| `demo_SIM_cell_seg_train.py` | Training PM/NE segmentation models | Model training and data preparation |
| `demo_SIM_cell_segmentation.py` | Cell segmentation inference | Segment plasma membrane and nucleus |
| `demo_SIM_ER_segmentation.py` | ER segmentation with ERNet | Segment endoplasmic reticulum |
| `demo_SIM_mito_segmentation.py` | Mitochondria segmentation | Segment mitochondria structures |
| `demo_SIM_organelle_segmentation.py` | Multi-organelle segmentation | Segment actin, ISG, cell shape, nuclei |

---

## 1. Cell Segmentation Training

Train deep learning models for PM (plasma membrane) and NE (nuclear envelope) segmentation.

### Basic Usage

```bash
# Prepare data and train model
python demo_SIM_cell_seg_train.py --main_path /path/to/ipa --run_data_prep --run_training

# Training only (data already prepared)
python demo_SIM_cell_seg_train.py --main_path /path/to/ipa --run_training
```

### Advanced Options

```bash
# Custom training parameters
python demo_SIM_cell_seg_train.py \
    --main_path /path/to/ipa \
    --run_training \
    --epochs 30 \
    --batch_size 16 \
    --learning_rate 0.0001 \
    --visualization
```

### Key Parameters
- `--epochs`: Number of training epochs (default: 20)
- `--batch_size`: Batch size (default: 8)
- `--learning_rate`: Learning rate (default: 1e-4)
- `--data_usage_ratio`: Fraction of data to use (default: 0.8)

---

## 2. Cell Segmentation Inference

Perform cell segmentation on SIM images using trained models.

### Basic Usage

```bash
# Segment single dataset
python demo_SIM_cell_segmentation.py \
    --main_path /path/to/ipa \
    --data_id 784_5 \
    --run_segmentation
```

### Multiple Datasets

```bash
# Process multiple datasets
python demo_SIM_cell_segmentation.py \
    --main_path /path/to/ipa \
    --data_id 784_5 784_6 784_7 \
    --run_segmentation
```

### Key Parameters
- `--data_id`: One or more dataset IDs
- `--pool_processes`: Number of parallel processes (default: 6)
- `--output_dir`: Custom output directory

---

## 3. ER Segmentation

Segment endoplasmic reticulum structures using ERNet model.

### Basic Usage

```bash
# Run with default paths
python demo_SIM_ER_segmentation.py \
    --main_path /path/to/ipa \
    --run_segmentation
```

### Custom Image

```bash
# Process custom image
python demo_SIM_ER_segmentation.py \
    --main_path /path/to/ipa \
    --run_segmentation \
    --image_path /path/to/image.png \
    --weights_path /path/to/model.pth
```

### CPU Mode

```bash
# Run on CPU (no GPU)
python demo_SIM_ER_segmentation.py \
    --main_path /path/to/ipa \
    --run_segmentation \
    --use_cpu
```

### Key Parameters
- `--image_size`: Processing image size (default: 1000)
- `--use_cpu`: Use CPU instead of GPU
- `--weights_path`: Custom model weights path

---

## 4. Mitochondria Segmentation

Segment mitochondria using vessel enhancement filters.

### Basic Usage

```bash
# High-quality segmentation
python demo_SIM_mito_segmentation.py \
    --main_path /path/to/ipa \
    --input_image /path/to/image.tif \
    --run_segmentation
```

### Quality Tuning

```bash
# Adjust quality parameters
python demo_SIM_mito_segmentation.py \
    --main_path /path/to/ipa \
    --input_image /path/to/image.tif \
    --run_segmentation \
    --downsample_factor 2 \
    --z_range_ratio 1.0 \
    --min_object_size_3d 80
```

### Key Parameters
- `--downsample_factor`: Downsampling factor (default: 2)
- `--z_range_ratio`: Z range to process, 0-1 (default: 1.0)
- `--sigma1`, `--sigma2`: Multi-scale vessel detection (default: 0.8, 2.5)
- `--min_object_size_3d`: Minimum object size (default: 80)

---

## 5. Multi-Organelle Segmentation

Segment multiple organelle types from SIM images.

### Basic Usage

```bash
# Run all segmentation types
python demo_SIM_organelle_segmentation.py \
    --main_path /path/to/ipa \
    --run_all \
    --visualization
```

### Specific Organelles

```bash
# Segment specific organelles
python demo_SIM_organelle_segmentation.py \
    --main_path /path/to/ipa \
    --run_actin \
    --run_nucleus \
    --visualization
```

### Key Parameters
- `--run_actin`: Actin filament skeletonization
- `--run_isg`: ISG sphere segmentation
- `--run_cell_shape`: Cell shape segmentation
- `--run_nucleus`: Nucleus segmentation
- `--run_all`: Run all segmentation types

---

## Common Parameters

All scripts support these common parameters:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--main_path` | Project main path (required) | - |
| `--output_dir` | Custom output directory | `{main_path}/outputs/{script_name}` |
| `--seed` | Random seed for reproducibility | 42 |
| `--visualization` | Enable visualization output | False |

---

## Output Structure

All scripts create timestamped output directories:

```
{main_path}/outputs/{script_name}/
└── YYYYMMDD-HHMMSS/
    ├── config.json              # Configuration and parameters
    ├── results/                 # Segmentation results
    ├── checkpoints/             # Model checkpoints (training only)
    ├── visualizations/          # Visualization images
    └── {script_name}_*.log      # Log files
```

---

## Troubleshooting

### Out of Memory Error

**Problem**: Script fails with "out of memory" error

**Solutions**:
- Training: Reduce `--batch_size` (e.g., from 8 to 4)
- Mito segmentation: Increase `--downsample_factor` or reduce `--crop_size`
- Cell segmentation: Reduce `--pool_processes`

### File Not Found Error

**Problem**: Script cannot find input files

**Solutions**:
1. Check file paths in error message
2. Verify `--main_path` is correct
3. For custom paths, use absolute paths
4. Ensure input files have correct naming format

### Model Loading Error

**Problem**: Cannot load pretrained model

**Solutions**:
1. Verify model weights exist at specified path
2. For ER segmentation, check ERNet pretrained model location
3. For cell segmentation training, ensure data preparation completed
4. Use `--weights_path` to specify custom model location

---

## Quick Start Examples

### Complete Training Pipeline

```bash
# 1. Prepare data and train model
python demo_SIM_cell_seg_train.py \
    --main_path /path/to/ipa \
    --run_data_prep \
    --run_training \
    --epochs 20 \
    --visualization

# 2. Run inference
python demo_SIM_cell_segmentation.py \
    --main_path /path/to/ipa \
    --data_id 784_5 \
    --run_segmentation
```

### Batch Processing Multiple Organelles

```bash
# Process all organelle types with visualization
python demo_SIM_organelle_segmentation.py \
    --main_path /path/to/ipa \
    --run_all \
    --visualization \
    --output_dir /custom/output
```

---

## Tips

1. **First-time use**: Run data preparation once, then training only for subsequent runs
2. **GPU usage**: Scripts automatically use GPU if available; use `--use_cpu` to force CPU mode
3. **Large datasets**: Use `--data_usage_ratio` to train on subset of data for faster iteration
4. **Reproducibility**: Set consistent `--seed` value across runs
5. **Monitoring**: Check log files in output directory for detailed progress

---

For more details, see individual script docstrings or run with `--help`:

```bash
python demo_SIM_cell_seg_train.py --help
```





