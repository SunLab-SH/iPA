from setuptools import setup, find_packages

setup(
    name="ipa",
    version="0.1.0",
    description="Multi-Modal Image Processing Toolkit",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    author="liad",
    packages=find_packages(),
    python_requires=">=3.8,<3.12",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0",
        "matplotlib>=3.3.0",
        "scikit-image>=0.18.0",
        "tifffile>=2021.7.2",
        "pandas>=1.3.0",
        "tqdm>=4.60.0",
        "plotly>=5.0.0",
        "Pillow>=8.0.0",
        "mrcfile>=1.3.0",
        "h5py>=3.0.0",
    ],
    extras_require={
        "ml": [
            "torch>=1.9.0",
            "torchvision>=0.10.0",
            "n2v>=0.3.0",
        ],
        "cv": [
            "opencv-python>=4.5.0",
        ],
        "viz": [
            "seaborn>=0.11.0",
        ],
        "fileformats": [
            "readlif>=0.6.0",
            "czifile>=2019.7.2", 
            "nd2reader>=3.2.0",
            "aicsimageio>=4.0.0",
            "imageio>=2.9.0",
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
            "torch>=1.9.0",
            "torchvision>=0.10.0",
            "opencv-python>=4.5.0",
            "seaborn>=0.11.0",
            "jupyter>=1.0.0",
            "jupyterlab>=3.0.0",
            "ipywidgets>=7.6.0",
            "scikit-learn>=1.0.0",
            "n2v>=0.3.0",
            "readlif>=0.6.0",
            "czifile>=2019.7.2",
            "nd2reader>=3.2.0", 
            "aicsimageio>=4.0.0",
            "imageio>=2.9.0",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
)




