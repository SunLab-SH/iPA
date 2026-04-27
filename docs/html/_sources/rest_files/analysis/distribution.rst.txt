Distribution Analysis
=====================

The distribution module provides tools for analyzing spatial distributions and radial distribution functions (RDF) of cellular components.


Functions
---------

Radial Distribution Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: ipa.analysis.calculate_rdf_from_xvg

   Calculates radial distribution function from partition coordinate data.
   
   **Parameters:**
   
   * `organelle_mask` (ndarray): 3D mask of organelle locations
   * `partition_coords` (dict): Dictionary of partition coordinates
   * `radial_positions` (list): List of radial position values
   * `bins` (int): Number of radial bins (default: 8)
   
   **Returns:**
   
   * Dictionary containing RDF results with radii and RDF values

.. autofunction:: ipa.analysis.load_partition_coords_from_xvg

   Loads partition coordinates and radial positions from XVG format files.
   
   **Parameters:**
   
   * `xvg_file` (str): Path to XVG coordinate file
   
   **Returns:**
   
   * Tuple of (partition_coords dict, radial_positions list)

.. autofunction:: ipa.analysis.save_radial_rdf_results

   Saves RDF analysis results to file.
   
   **Parameters:**
   
   * `rdf_results` (dict): RDF analysis results
   * `output_dir` (str): Output directory path  
   * `data_id` (str): Dataset identifier

Spatial Arrangement Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: ipa.analysis.reconstruct_4d_mask

   Reconstructs 4D mask data from coordinate information.
   
   **Parameters:**
   
   * `data` (DataFrame): Coordinate data
   * `resolution` (tuple): Image resolution (width, height, depth)
   * `rescale_x` (float): X-axis rescaling factor
   * `rescale_y` (float): Y-axis rescaling factor
   
   **Returns:**
   
   * Tuple of (reconstructed image, vesicle coordinates)

.. autofunction:: ipa.analysis.analyze_vesicle_clusters

   Analyzes vesicle clustering patterns and spatial organization.
   
   **Parameters:**
   
   * `sample_name` (str): Sample identifier
   * `img_ch2` (ndarray): Channel 2 image data
   * `img_reconstructed` (ndarray): Reconstructed image data
   * `ves_coords` (list): Vesicle coordinate list
   
   **Returns:**
   
   * Dictionary containing cluster analysis results

.. autofunction:: ipa.analysis.analyze_docked_granules

   Analyzes docking relationships between granules and membrane structures.
   
   **Parameters:**
   
   * `segmented_labels` (ndarray): Segmented granule labels
   * `mem_mask` (ndarray): Membrane mask
   * `isg_data` (ndarray): ISG coordinate data
   * `show_visualization` (bool): Whether to display visualizations
   
   **Returns:**
   
   * Dictionary with docking statistics

Enhanced Distribution Methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: ipa.analysis.actinangle2dist_enhanced

   Enhanced analysis of actin filament angle-distance distributions.

Usage Examples
--------------

RDF Analysis for SXT Data
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from ipa.analysis import (load_partition_coords_from_xvg, 
                            calculate_rdf_from_xvg, 
                            save_radial_rdf_results)
   from ipa.processing.partitioning import plot_radial_rdf
   import tifffile
   
   # Load partition data
   partition_coords, radial_positions = load_partition_coords_from_xvg(
       "partition_coords.xvg"
   )
   
   # Load organelle mask
   organelle_mask = tifffile.imread("organelle_mask.tiff")
   
   # Calculate RDF
   rdf_results = calculate_rdf_from_xvg(
       organelle_mask, 
       partition_coords, 
       radial_positions, 
       bins=8
   )
   
   # Plot and save results
   plot_radial_rdf(rdf_results, save_path="rdf_plot.png")
   save_radial_rdf_results(rdf_results, "./results", "dataset_id")

WFM Cluster Analysis
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from ipa.analysis import analyze_vesicle_clusters, reconstruct_4d_mask
   import pandas as pd
   
   # Load coordinate data
   data = pd.read_csv("vesicle_positions.csv", skiprows=3)
   
   # Reconstruct 4D data
   img_reconstructed, ves_coords = reconstruct_4d_mask(
       data, 
       resolution=(512, 512, 100),
       rescale_x=2.1, 
       rescale_y=1.7
   )
   
   # Analyze clusters
   results = analyze_vesicle_clusters(
       "sample_name",
       img_ch2, 
       img_reconstructed, 
       ves_coords
   )

Docked Granule Analysis
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from ipa.analysis import analyze_docked_granules, extract_isg_features
   
   # Extract ISG features first
   df_results, segmented_labels = extract_isg_features(
       img_isg, 
       min_distance=5, 
       min_size=20
   )
   
   # Analyze docking
   docking_results = analyze_docked_granules(
       segmented_labels,
       membrane_mask,
       isg_coordinate_data,
       show_visualization=True
   )
