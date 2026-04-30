# Integrated Processing and Analysis Toolkit (iPA)

<div align="center">

```
​**​*​**​*​**​*​**​*​**​*​**​*​**​*​**​**​**​**​
*    ██╗ ██████╗   █████╗     *
*    ██║ ██╔══██╗ ██╔══██╗    *
*    ██║ ██████╔╝ ███████║    *
*    ██║ ██╔═══╝  ██╔══██║    *
*    ██║ ██║      ██║  ██║    *
*    ╚═╝ ╚═╝      ╚═╝  ╚═╝    *
​**​*​**​*​**​*​**​*​**​*​**​*​**​*​**​**​**​**​
```





**A comprehensive toolkit for multi-modal cellular imaging analysis**
[![Python](https://img.shields.io/badge/Python-3.9-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Documentation](https://img.shields.io/badge/docs-available-brightgreen.svg)](docs/)
[![Docker Pulls](https://img.shields.io/docker/pulls/a13707621/ipa.svg)](https://hub.docker.com/r/a13707621/ipa)

</div>

---

**Developed by:** Angdi Li ([liad@shanghaitech.edu.cn](mailto:liad@shanghaitech.edu.cn))  
**Institution:** ShanghaiTech University



## Introduction

iPA is an open-source Python platform for comprehensive analysis of cellular and subcellular imaging data. Our toolkit streamlines the workflow from raw image processing to advanced quantitative analysis, enabling researchers to extract biologically meaningful features from diverse microscopy modalities.



For detailed documentation, tutorials, and examples, please visit:  
[Documentation on ReadTheDocs](https://ipa.readthedocs.io/en/latest/)

---

## Quick Start

### 1. Installation

#### Option A: Docker Installation (Recommended for Quick Start)

We provide an official Docker image that includes all dependencies and supports both GPU and CPU.

**Prerequisites:**
*   [Docker](https://docs.docker.com/get-docker/)
*   [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html) (for GPU support)

**Quick Start:**
```bash
# Pull the latest image
docker pull a13707621/ipa:latest

# Run with GPU support and mount your local data folder
docker run --gpus all -it -v $(pwd)/data:/app/data a13707621/ipa:latest

# If you don't have a GPU, simply remove the --gpus all flag:
# docker run -it -v $(pwd)/data:/app/data a13707621/ipa:latest
```

Inside the container, you can run any example script:
```bash
python examples/examples_et/demo_cryoET_filament_predict.py
```

#### Option B: Local Installation

```bash
# Clone the repository
git clone https://github.com/SunLab-SH/iPA.git
cd iPA

# Initialize Git LFS (required for pre-trained models)
git lfs install
git lfs pull

# Install dependencies and iPA in development mode
pip install -r requirements.txt
pip install -e .
```

### 2. Download Example Data

Due to their large size, example data files are hosted on **Zenodo**. 

1.  Download `data.zip` from: [Zenodo Record 19903896](https://zenodo.org/records/19903896)
2.  Extract the archive into the `data/` directory:
    ```bash
    unzip data.zip -d data/
    ```

**Note**: Some data files are gzip/zlib compressed. The `UniversalDataLoader` will automatically decompress them at runtime.

For more details, see [DATA_SETUP_GUIDE.md](DATA_SETUP_GUIDE.md).

### 3. Run Examples

```bash
# Cryo-ET filament segmentation
python examples/examples_et/demo_cryoET_filament_predict.py

# SXT cell analysis
python examples/examples_sxt/demo_SXT_composite_mask.py

# SIM ISG segmentation
python examples/examples_sim_wfm/demo_SIM_isg_segmentation.py
```

All example scripts are located in the `examples/` directory, organized by imaging modality.

---

#### Core Dependencies

*   numpy, scipy, pandas, matplotlib
*   scikit-image, tifffile, Pillow, mrcfile, imageio
*   torch, torchvision
*   tqdm, h5py

*Optional dependencies (for training or advanced visualization) can be installed via:* `pip install -e ".[all]"`

---

#### Example Data & Models

*   **Data**: Available at [Zenodo Record 19903896](https://zenodo.org/records/19903896). Unzip into the `data/` folder.
*   **Models**: Pre-trained weights (`.pth`) are managed via **Git LFS**. Ensure you run `git lfs pull` after cloning.





### Overview


![iPA Analysis Pipeline](./workflow_images/figure_1_v39.jpg)

*Figure 1: iPA Analysis Pipeline Overview*




#### Install from local source

**Requirements:** Python 3.9 or 3.10

```bash
# Create conda environment (Python 3.9 recommended)
conda create -n ipa_env python=3.9
conda activate ipa_env

# Install Git LFS (if not already installed)
# On Ubuntu/Debian:
sudo apt install git-lfs
# On macOS:
brew install git-lfs
# On Windows: download from https://git-lfs.github.io/

# Initialize Git LFS
git lfs install

# Clone the repository
git clone https://github.com/SunLab-SH/iPA.git
cd iPA

# Pull large files (including .pth model files)
git lfs pull

# Install dependencies and the package
pip install -r requirements.txt
pip install -e .
```

#### Verify Installation

Test that iPA is correctly installed by running a simple demo:

```bash
# Run a basic Cryo-ET filament prediction demo
python examples/examples_et/demo_cryoET_filament_predict.py
```

If the script runs without import errors and produces output, your installation is successful.



<!-- ## Citation

If you use iPA in your research, please cite:

```bibtex
@software{iPA2024,
  title={Integrated Processing and Analysis Toolkit for Multi-Scale Cellular Imaging},
  author={Li, Angdi and others},
  year={2024},
  url={https://github.com/SunLab-SH/iPA}
}
``` -->