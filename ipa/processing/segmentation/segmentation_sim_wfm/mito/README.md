# SIM Mitochondria Segmentation - CLI Demo

Simple command-line tool for SIM mitochondria segmentation.

## Usage

```bash
python cli_demo.py --input path/to/image.tif [options]
```

## Examples

### Basic usage (default parameters):
```bash
python cli_demo.py --input image.tif --output results
```

### Fast processing:
```bash
python cli_demo.py --input image.tif --downsample 8 --max_z_layers 15 --max_crop_size 512
```

### High quality processing:
```bash
python cli_demo.py --input image.tif --downsample 2 --max_z_layers 0 --max_crop_size 2048 --min_size_3d 100
```

### Custom enhancement parameters:
```bash
python cli_demo.py --input image.tif --sigma1 1.2 --sigma2 2.5 --laplacian_weight 0.3
```

### Quiet mode:
```bash
python cli_demo.py --input image.tif --quiet
```

## Parameters

### Required:
- `--input`, `-i`: Input image file path

### Optional:
- `--output`, `-o`: Output directory (default: 'segmentation_results')

### Processing control:
- `--downsample`: Downsampling factor, higher = faster (default: 4)
- `--max_z_layers`: Max Z layers to process, 0 = all (default: 30)
- `--max_crop_size`: Max crop size for XY (default: 1024)

### Enhancement:
- `--sigma1`: Fine structure Gaussian sigma (default: 1.0)
- `--sigma2`: Coarse structure Gaussian sigma (default: 2.0)
- `--laplacian_weight`: Laplacian enhancement weight (default: 0.5)

### Segmentation:
- `--min_size_3d`: Minimum object size for 3D (default: 50)
- `--min_size_2d`: Minimum object size for 2D (default: 20)
- `--morph_size`: Morphological operations size (default: 1)

### Output:
- `--quiet`, `-q`: Suppress progress output

## Output Files

The tool creates the following files in the output directory:
- `mitochondria_mask.tif`: **Main result** - Binary mask (0-1 values, original image size)
- `middle_slice_original.tif`: Original middle slice for visualization
- `middle_slice_mask.tif`: Mask middle slice for visualization

## Parameter Recommendations

### For speed (preview/testing):
```bash
--downsample 8 --max_z_layers 15 --max_crop_size 512
```

### For balance (recommended):
```bash
--downsample 4 --max_z_layers 30 --max_crop_size 1024
```

### For quality (final results):
```bash
--downsample 2 --max_z_layers 0 --max_crop_size 2048 --min_size_3d 100
```

### For fine structures:
```bash
--sigma1 0.7 --sigma2 1.3 --laplacian_weight 0.7 --min_size_3d 20
```