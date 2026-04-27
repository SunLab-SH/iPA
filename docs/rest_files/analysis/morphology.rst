Morphology Analysis
===================

The morphology module provides tools for analyzing cellular and subcellular structure morphological features.

Overview
--------

This module contains specialized morphological analysis functions for different imaging modalities:


Functions
---------

Actin Analysis
~~~~~~~~~~~~~~

.. autofunction:: ipa.analysis.actin_to_actin_analysis

   Analyzes spatial relationships within actin filament networks in Cryo-ET data.
   
   **Parameters:**
   
   * `data_id` (str): Dataset identifier
   * `mask_file` (str): Path to vesicle mask file (.mrc)
   * `actin_file` (str): Path to actin coordinates file (.json)
   * `config` (dict): Configuration parameters including voxel size and offsets
   
   **Returns:**
   
   * Dictionary containing analysis results with statistical measures

.. autofunction:: ipa.analysis.actin_angle_distance_pair_enhanced_analysis

   Enhanced analysis of actin angle-distance pair relationships with advanced statistical measures.
   
   **Parameters:**
   
   * `data_id` (str): Dataset identifier  
   * `actin_file` (str): Path to actin coordinates file
   * `vesicle_file` (str): Path to vesicle mask file
   * `config` (dict): Configuration parameters
   
   **Returns:**
   
   * Dictionary with enhanced analysis metrics

ISG Feature Extraction
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: ipa.analysis.extract_isg_features

   Extracts morphological features from insulin secretory granules in SIM/WFM data.
   
   **Parameters:**
   
   * `img_isg` (ndarray): ISG image data
   * `min_distance` (int): Minimum distance between peaks (default: 5)
   * `min_size` (int): Minimum size threshold (default: 20)
   * `save_csv` (bool): Whether to save results to CSV
   * `output_path` (str): Output file path
   
   **Returns:**
   
   * Tuple of (DataFrame with features, segmented labels)

Usage Examples
--------------

Actin Network Analysis
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from ipa.analysis import actin_to_actin_analysis
   
   # Configuration for Cryo-ET analysis
   config = {
       'voxel_size_xyz': [13.412, 13.412, 13.412],
       'shift_bias': [3046.30, 3214.95, -6171.17],
       'visualization': True,
       'save_intermediate': True
   }
   
   results = actin_to_actin_analysis(
       data_id="20220326_2.8_2",
       mask_file="vesicle_mask.mrc",
       actin_file="actin_points.json",
       config=config
   )

ISG Morphological Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from ipa.analysis import extract_isg_features
   from ipa.data_loader import UniversalDataLoader
   
   # Load ISG image data
   img_isg = UniversalDataLoader.load_data("isg_data.tif")
   
   # Extract features
   df_results, segmented_labels = extract_isg_features(
       img_isg,
       min_distance=5,
       min_size=20,
       save_csv=True,
       output_path="isg_features.csv"
   )
