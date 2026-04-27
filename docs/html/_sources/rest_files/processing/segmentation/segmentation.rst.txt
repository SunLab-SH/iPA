Segmentation Module
===================

The segmentation module provides computational tools for automated cellular structure segmentation, leveraging advanced deep learning architectures optimized for different microscopy modalities.

Overview
--------

This module contains specialized segmentation tools for various microscopy techniques:

* **SXT Segmentation**: Soft X-ray Tomography cellular structure identification
* **ET Segmentation**: Cryo-Electron Tomography structural analysis
* **SIM/WFM Segmentation**: Structured Illumination and Widefield Microscopy processing

The segmentation module is designed to handle data from multiple imaging modalities with modality-specific optimization and preprocessing pipelines.

Submodules
----------

.. toctree::
   :maxdepth: 2

   sxt_segmentation
   et_segmentation
   sim_wfm_segmentation
