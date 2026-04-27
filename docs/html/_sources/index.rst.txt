.. iPA documentation master file, created by
   sphinx-quickstart on Sun Jul 28 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Integrated Processing and Analysis Toolkit (iPA)




iPA is an open-source Python toolkit for comprehensive analysis of cellular and subcellular imaging data. 
Our platform streamlines the workflow from raw image processing to advanced quantitative analysis, enabling 
researchers to extract biologically meaningful features from diverse microscopy modalities including 
cryo-electron tomography ( **cryo-ET**), soft X-ray tomography ( **SXT**), structured illumination microscopy ( **SIM**), and widefield microscopy( **WFM**).







Multi-Modal Cellular Imaging Analysis
======================================

* The included cellular imaging analysis pipeline can be modified to:
    * Load various microscopy data formats (MRC, TIFF, etc.)
    * Process and segment organelles (mitochondria, vesicles, nucleus and whole-cells)
    * Extract morphological and spatial features
    * Perform partitioning on cytoplasmic regions
    * Generate visualization results



.. image:: rest_files/resources/figure_1_v25.jpg
   :alt: iPA Overview
   :align: center
   :width: 600px


.. Radial Cytoplasmic Partitioning
.. ===============================

.. * Novel spatial zoning from nucleus to membrane

.. .. image:: ./rest_files/resources/partitioning_example.png 
..     :width: 600
..     :alt: Radial partitioning example
    

.. Organelle Spatial Analysis
.. ===========================

.. * Distance and angular analysis between cellular structures

.. .. image:: ./rest_files/resources/spatial_analysis.png 
..     :width: 600
..     :alt: Spatial analysis example

.. * Multi-organelle interaction mapping

.. .. image:: ./rest_files/resources/interaction_mapping.png 
..     :width: 600
..     :alt: Interaction mapping example


.. CryoET Data Processing
.. ======================

.. * Functions for processing cryo-electron tomography data are included in the `et_analysising` module

.. .. image:: ./rest_files/resources/cryoet_processing.png
..     :width: 900
..     :alt: CryoET processing pipeline



 

Dependencies 
------------

    * numpy
    * mrcfile
    * tifffile
    * matplotlib
    * pandas
    * scipy
    * scikit-image
    * tqdm
    * plotly


    if you use segmentation modules, you will also need:
    * pytorch



.. Examples 
.. --------

.. This toolkit contains comprehensive examples for various cellular imaging analysis tasks:
..     * CryoET organelle interaction analysis
..     * Multi-modal spatial correlation studies  
..     * Radial cytoplasmic partitioning workflows
..     * Cross-platform validation pipelines
    
    
..     `Examples On GitHub <https://github.com/SunLab-SH/iPA/tree/main/examples>`_


Supported Imaging Modalities
-----------------------------

* **Cryo-ET**: Cryo-electron tomography analysis
* **SIM**: Structured illumination microscopy
* **SXT**: Soft X-ray tomography  
* **WFM**: Wide-field microscopy


.. Cite 
.. -----

.. If you use iPA in your research, please cite:

.. .. code-block:: bibtex

.. ..    @software{iPA2024,
..       title={Integrated Processing and Analysis Toolkit for Multi-Scale Cellular Imaging},
..       author={Li, Angdi and others},
..       year={2024},
..       url={https://github.com/SunLab-SH/iPA}
..     }


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   rest_files/data_loader/data_loader
   rest_files/processing/processing
   rest_files/analysis/analysis





..    :caption: API Reference:
..    :hidden:
  
..    rest_files/data_loader
..    rest_files/denoising
..    rest_files/segmentation
..    rest_files/partitioning
..    rest_files/processing
..    rest_files/et_analysising
..    rest_files/visualization

.. .. toctree::
..    :maxdepth: 1
..    :caption: Tutorials:
..    :hidden:

..    rest_files/examples

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
* :ref:`search`
