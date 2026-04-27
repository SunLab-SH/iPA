Interaction Analysis
====================

The interaction module provides comprehensive tools for analyzing spatial interactions between different cellular components in Cryo-ET data.


Functions
---------

Actin Interaction Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: ipa.analysis.actin_to_vesicle_analysis

   Analyzes spatial interactions between actin filaments and vesicles.
   
   **Parameters:**
   
   * `data_id` (str): Dataset identifier
   * `mask_file` (str): Path to vesicle mask file
   * `json_file` (str): Path to actin coordinates file
   * `config` (dict): Configuration parameters
   
   **Returns:**
   
   * Dictionary containing interaction analysis results

.. autofunction:: ipa.analysis.actin_to_vesicle_dist_angle_pair_analysis

   Combined distance-angle analysis for actin-vesicle interactions.
   
   **Parameters:**
   
   * `data_id` (str): Dataset identifier
   * `mask_file` (str): Path to vesicle mask file
   * `json_file` (str): Path to actin coordinates file
   * `angle_file` (str): Path to angle data file
   * `config` (dict): Configuration parameters
   
   **Returns:**
   
   * Dictionary with distance-angle pair analysis results

.. autofunction:: ipa.analysis.actin_to_microtube_analysis

   Analyzes spatial relationships between actin filaments and microtubules.

.. autofunction:: ipa.analysis.actin_to_mito_analysis

   Quantifies interactions between actin filaments and mitochondria.

.. autofunction:: ipa.analysis.actin_to_endoreticulum_analysis

   Analyzes actin-endoplasmic reticulum spatial interactions.

Microtubule Interaction Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: ipa.analysis.microtube_to_actin_analysis

   Analyzes microtubule-actin filament interactions from microtubule perspective.

.. autofunction:: ipa.analysis.microtube_to_vesicle_analysis

   Quantifies spatial relationships between microtubules and vesicles.

.. autofunction:: ipa.analysis.microtube_to_mito_analysis

   Analyzes microtubule-mitochondria interactions.

.. autofunction:: ipa.analysis.microtube_to_endoreticulum_analysis

   Studies microtubule-ER spatial relationships.

Organelle-Organelle Interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: ipa.analysis.vesicle_to_mito_analysis

   Analyzes direct vesicle-mitochondria spatial interactions.

.. autofunction:: ipa.analysis.vesicle_to_endoreticulum_analysis

   Quantifies vesicle-ER spatial relationships.

.. autofunction:: ipa.analysis.mito_to_endoreticulum_analysis

   Analyzes mitochondria-endoplasmic reticulum interactions.

Usage Examples
--------------

Basic Actin-Vesicle Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from ipa.analysis import actin_to_vesicle_analysis
   
   # Configuration for Cryo-ET data
   config = {
       'voxel_size_xyz': [13.412, 13.412, 13.412],
       'shift_bias': [3046.30, 3214.95, -6171.17],
       'visualization': True,
       'save_intermediate': True,
       'output_dir': "./results"
   }
   
   # Analyze actin-vesicle interactions
   results = actin_to_vesicle_analysis(
       data_id="20220326_2.8_2",
       mask_file="vesicle_mask.mrc",
       json_file="actin_coordinates.json",
       config=config
   )

Advanced Distance-Angle Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from ipa.analysis import actin_to_vesicle_dist_angle_pair_analysis
   
   # Combined distance-angle analysis
   results = actin_to_vesicle_dist_angle_pair_analysis(
       data_id="20220326_2.8_2",
       mask_file="vesicle_mask.mrc", 
       json_file="actin_coordinates.json",
       angle_file="actin_angles.json",
       config=config
   )

Multi-Organelle Interaction Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from ipa.analysis import (actin_to_microtube_analysis,
                            actin_to_mito_analysis,
                            vesicle_to_mito_analysis,
                            mito_to_endoreticulum_analysis)
   
   # Analyze multiple interaction pairs
   interactions = {}
   
   # Actin interactions
   interactions['actin_microtube'] = actin_to_microtube_analysis(
       data_id, actin_file, microtube_file, config
   )
   
   interactions['actin_mito'] = actin_to_mito_analysis(
       data_id, actin_file, mito_file, config  
   )
   
   # Organelle interactions
   interactions['vesicle_mito'] = vesicle_to_mito_analysis(
       data_id, vesicle_file, mito_file, config
   )
   
   interactions['mito_er'] = mito_to_endoreticulum_analysis(
       data_id, mito_file, er_file, config
   )

Configuration Parameters
------------------------

All interaction analysis functions accept a configuration dictionary with the following parameters:

* `voxel_size_xyz` (list): Voxel dimensions in xyz [nm]
* `shift_bias` (list): Coordinate offset correction [nm]  
* `visualization` (bool): Enable result visualization
* `save_intermediate` (bool): Save intermediate processing results
* `output_dir` (str): Directory for output files
* `zoom_degree` (float): Visualization zoom factor
* `batch_size` (int): Processing batch size for memory management
