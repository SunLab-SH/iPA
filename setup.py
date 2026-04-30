from setuptools import setup, find_packages

setup(
    name="ipa",
    version="0.1.0",
    description="Multi-Modal Image Processing Toolkit",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    author="liad",
    packages=find_packages(),
    python_requires=">=3.9,<3.11",
    install_requires=[
        # Core scientific computing
        "numpy>=1.21.0,<2.0.0",
        "scipy>=1.9.0",
        "pandas>=1.5.0",
        "matplotlib>=3.6.0",
        
        # Image processing
        "scikit-image>=0.20.0",
        "tifffile>=2022.4.8",
        "Pillow>=9.0.0",
        "mrcfile>=1.3.0",
        "imageio>=2.25.0",
        
        # Deep learning
        "torch>=1.13.0",
        "torchvision>=0.14.0",
        
        # Utilities
        "tqdm>=4.60.0",
        "h5py>=3.7.0",
    ],
    extras_require={
        "cv": [
            "opencv-python>=4.6.0",
        ],
        "viz": [
            "plotly>=5.0.0",
            "seaborn>=0.12.0",
        ],
        "ml": [
            "scikit-learn>=1.2.0",
        ],
        "fileformats": [
            "readlif>=0.6.0",
            "czifile>=2019.7.2", 
            "nd2reader>=3.2.0",
            "aicsimageio>=4.10.0",
        ],
        "training": [
            "tensorboard>=2.10.0",
            "wandb>=0.13.0",
            "pycocotools>=2.0.6",
        ],
        "perf": [
            # numba removed as per project cleanup
        ],
        "jupyter": [
            "jupyter>=1.0.0",
            "jupyterlab>=3.0.0",
            "ipywidgets>=7.6.0",
        ],
        "dev": [
            "pytest>=6.0.0",
            "black>=21.0.0",
            "flake8>=3.8.0",
            "sphinx>=4.0.0",
            "sphinx-rtd-theme>=1.0.0",
        ],
        "all": [
            "opencv-python>=4.6.0",
            "plotly>=5.0.0",
            "seaborn>=0.12.0",
            "scikit-learn>=1.2.0",
            "readlif>=0.6.0",
            "czifile>=2019.7.2",
            "nd2reader>=3.2.0",
            "aicsimageio>=4.10.0",
            "pycocotools>=2.0.6",
            "jupyter>=1.0.0",
            "jupyterlab>=3.0.0",
            "ipywidgets>=7.6.0",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
    ],
)




