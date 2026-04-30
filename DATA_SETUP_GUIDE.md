# iPA Data and Model Setup Guide

This guide explains how to prepare the data files and pre-trained models required to run iPA.

---

## 📋 Overview

Due to their large size, example data files are hosted on **Zenodo**, while pre-trained model files (`.pth`) are managed via **Git LFS** within the repository.

---

## 🔧 Setup Steps

### 1. Download Example Data from Zenodo

The complete dataset (approx. 800 MB) is available as a single archive.

**Steps:**
1.  Download `data.zip` from: [Zenodo Record 19903896](https://zenodo.org/records/19903896)
2.  Create a `data` directory in the project root:
    ```bash
    mkdir -p data
    ```
3.  Extract the archive into the `data` directory:
    ```bash
    unzip data.zip -d data/
    ```

**Note on Compressed Files:**
Some large files (e.g., `.mrc`, `.tif`) inside the archive are stored in gzip or zlib compressed formats to save space. The iPA `UniversalDataLoader` automatically detects and transparently decompresses these files at runtime. **You do not need to manually decompress them.**

### 2. Retrieve Pre-trained Models via Git LFS

Model files are tracked using Git Large File Storage (LFS).

**Steps:**
1.  Ensure Git LFS is installed:
    ```bash
    # Ubuntu/Debian
    sudo apt-get install git-lfs
    # macOS
    brew install git-lfs
    ```
2.  Initialize Git LFS (if not already done):
    ```bash
    git lfs install
    ```
3.  Pull the model files:
    ```bash
    git lfs pull
    ```

The models will be automatically placed in `ipa/processing/segmentation/models/`.

---

## 📊 Dataset Content

| Modality | Folder | Description |
| :--- | :--- | :--- |
| **Cryo-ET** | `data/cryoET/` | Denoised subvolumes and coordinate files for filaments/vesicles. |
| **SXT** | `data/sxt/` | Stevens pancreatic cell dataset with labeled masks. |
| **SIM** | `data/sim/` | Super-resolution images for ISG and Actin analysis. |
| **WFM** | `data/wfm/` | Time-lapse tracking data and membrane/nucleus masks. |
| **Other** | `data/other/` | Supplementary fluorescence microscopy images. |

---

## ❓ Troubleshooting

### Q1: "File not found" errors when running demos?
**A**: 
- Verify that you have extracted `data.zip` into the `data/` folder at the project root.
- Check that the directory structure matches the one described above.

### Q2: "Model not loaded" errors?
**A**: 
- Ensure you have run `git lfs pull`.
- Check if the `.pth` files in `ipa/processing/segmentation/models/` are actual model files (several MBs) and not small text files (which would indicate LFS pointers).

### Q3: Slow download speeds from Zenodo?
**A**: 
- Zenodo uses a global CDN. If speeds are low, try downloading during off-peak hours or use a download manager that supports multi-threading.

---

## 🔗 Useful Links

- **iPA Documentation**: [https://ipa.readthedocs.io/](https://ipa.readthedocs.io/)
- **Example Data (Zenodo)**: [https://zenodo.org/records/19903896](https://zenodo.org/records/19903896)
- **Issue Tracker**: GitHub Issues

---

## 📞 Support

For technical support or questions, please contact:
- **Developer**: Angdi Li
- **Email**: liad@shanghaitech.edu.cn
- **Institution**: ShanghaiTech University

---

**Last Updated**: 2026-04-30
