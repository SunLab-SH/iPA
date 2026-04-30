Example Scripts
===============

iPA provides a comprehensive set of example scripts to demonstrate its capabilities across different imaging modalities. All scripts are located in the ``examples/`` directory.

Cryo-Electron Tomography (Cryo-ET)
----------------------------------

*   **Filament Prediction**: ``examples/examples_et/demo_cryoET_filament_predict.py``
    *   Demonstrates filament skeletonization and prediction on denoised subvolumes.
*   **Morphology Analysis**: ``examples/examples_et/demo_ET_morphology.py``
    *   Calculates tubular length and filament statistics.
*   **Spatial Arrangement**: ``examples/examples_et/demo_ET_arrangement.py``
    *   Analyzes the angular distribution of anchored filaments.
*   **Interaction Analysis**: ``examples/examples_et/demo_ET_interaction.py``
    *   Studies interactions between actin, vesicles, mitochondria, and ER.

Soft X-ray Tomography (SXT)
---------------------------

*   **Composite Mask Generation**: ``examples/examples_sxt/demo_SXT_composite_mask.py``
    *   Combines cell, nucleus, mito, and ISG masks into a multi-class label map.
*   **Spatial Partitioning**: ``examples/examples_sxt/demo_SXT_partitioning.py``
    *   Divides cytoplasm into concentric shells based on NE and PM boundaries.
*   **Spatial Arrangement (RDF)**: ``examples/examples_sxt/demo_SXT_spatial_arrangement.py``
    *   Computes Radial Distribution Functions for organelles.
*   **Interaction Analysis**: ``examples/examples_sxt/demo_SXT_interaction.py``
    *   Analyzes ISG-Mitochondria contacts and distances.
*   **Morphology Analysis**: ``examples/examples_sxt/demo_SXT_morphology.py``
    *   Extracts volume distributions for individual ISGs and Mitochondria.
*   **Denoising (N2N/N2V)**: ``examples/examples_sxt/demo_SXT_n2n_predict.py``, ``demo_SXT_n2v_predict.py``
    *   Applies Noise2Noise and Noise2Void for tomogram enhancement.

Structured Illumination Microscopy (SIM) & Widefield (WFM)
----------------------------------------------------------

*   **ISG Segmentation**: ``examples/examples_sim_wfm/demo_SIM_isg_segmentation.py``
    *   Threshold-based segmentation for insulin granules.
*   **ER Prediction**: ``examples/examples_sim_wfm/demo_SIM_ER_predict.py``
    *   Uses the ERNet transformer model for ER segmentation.
*   **RDF Analysis**: ``examples/examples_sim_wfm/demo_SIM_rdf_analysis.py``
    *   Calculates organelle distribution relative to cell boundaries.
*   **Partitioning**: ``examples/examples_sim_wfm/demo_SIM_partitioning.py``, ``demo_WFM_partitioning.py``
    *   Spatial partitioning for SIM and WFM datasets.
*   **Dynamics Analysis**: ``examples/examples_sim_wfm/demo_WFM_dynamics.py``
    *   Analyzes 3D and radial velocities from tracking data.
*   **Interaction**: ``examples/examples_sim_wfm/demo_SIM_WFM_interaction.py``
    *   Cross-organelle distance analysis in fluorescence data.

Cross-Modal Comparison
----------------------

*   **SIM vs SXT**: ``examples/examples_cross_modal/demo_SIM_SXT_comparison.py``
    *   Compares organelle distribution patterns across different imaging techniques using Pearson correlation.

Running Examples
----------------

All examples can be run directly from the command line:

.. code-block:: bash

    # Navigate to examples directory
    cd examples/examples_sxt
    
    # Run a demo script
    python demo_SXT_partitioning.py

Each example includes detailed comments and generates output in the ``results/`` directory with corresponding log files in ``logs/``.
