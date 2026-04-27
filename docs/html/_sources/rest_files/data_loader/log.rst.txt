Log 
=================================

The log module provides logging utilities for data analysis workflows with two main components: the ``@log_analysis`` decorator for automatic function logging and the ``QuickLogger`` class for manual step-by-step logging.

.. autofunction:: ipa.data_loader.log_analysis

.. autoclass:: ipa.data_loader.QuickLogger
   :no-members:

Overview
--------

The log module offers two complementary approaches for tracking analysis workflows:

1. **@log_analysis decorator** - Automatic logging of function execution, parameters, and outputs
2. **QuickLogger class** - Manual logging for detailed step-by-step tracking

All logs are saved to timestamped files with UTF-8 encoding, creating a complete audit trail of analysis activities.

QuickLogger - Manual Step Logging
----------------------------------

The ``QuickLogger`` class provides a simple interface for logging analysis steps, input/output files, and errors. It has four main methods:

- **step(message)** - Log general analysis steps and information
- **file_in(file_path)** - Log input file paths
- **file_out(file_path)** - Log output file paths  
- **error(message)** - Log error messages

Basic Usage
^^^^^^^^^^^

.. code-block:: python

    from ipa.data_loader import QuickLogger
    
    # Initialize logger for SXT partitioning analysis
    logger = QuickLogger("sxt_partitioning", log_dir="logs")
    
    # Log analysis steps
    logger.step("Starting SXT Partitioning Demo")
    logger.step(f"Main path: {mainpath}")
    logger.step(f"Processing dataset: {dataid}")
    
    # Log input files
    logger.file_in(mask_pm_path)
    logger.file_in(mask_ne_path)
    
    # Log processing steps with details
    logger.step("Loading mask data...")
    logger.step(f"PM mask shape: {mask_data_pm.shape}, dtype: {mask_data_pm.dtype}")
    logger.step("Initializing partitioning processor...")
    
    # Log output files
    logger.file_out(xvg_path)
    logger.file_out(comparison_path)

Error Handling
^^^^^^^^^^^^^^

.. code-block:: python

    from ipa.data_loader import QuickLogger
    
    logger = QuickLogger("data_analysis")
    
    try:
        logger.step("Loading data from file")
        data = load_data(input_file)
        logger.step(f"Successfully loaded {data.shape} data")
        
    except Exception as e:
        logger.error(f"Failed to load data: {str(e)}")
        raise

@log_analysis - Automatic Function Logging
-------------------------------------------

The ``@log_analysis`` decorator automatically logs function execution, parameters, and results.

Function Decoration
^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from ipa.data_loader import log_analysis
    
    @log_analysis(log_dir="analysis_logs")
    def process_tomography_data(input_path, output_dir, threshold=0.5):
        """Process tomography data with automatic logging"""
        # Function logic here
        processed_files = []
        # ... processing ...
        return processed_files
    
    # Usage - logging happens automatically
    results = process_tomography_data("data/tomo/", "results/", threshold=0.7)

Parameter Detection
^^^^^^^^^^^^^^^^^^^

The decorator automatically detects and logs:

- **File paths**: Parameters containing 'path' or 'file' in the name
- **Output directories**: Parameters named 'output_dir', 'save_dir', 'result_path'
- **Return values**: File paths returned as strings, Paths, or lists

.. code-block:: python

    @log_analysis(log_dir="batch_logs")
    def batch_process_images(image_dir, mask_file, output_dir, n_slices=8):
        """Batch process images with comprehensive logging"""
        # Automatically logs:
        # - image_dir (contains 'dir')
        # - mask_file (contains 'file') 
        # - output_dir (recognized output parameter)
        # - Any returned file paths
        
        result_files = []
        # ... processing logic ...
        return result_files

Example from SXT Demo
---------------------------------

This example shows how logging is used in the SXT partitioning demo:

.. code-block:: python

    from ipa.data_loader import QuickLogger
    import tifffile
    import os
    
    def main():
        # Initialize logger with descriptive name
        log_dir = f'{mainpath}/logs'
        logger = QuickLogger("sxt_partitioning", log_dir=log_dir)
        logger.step("Starting SXT Partitioning Demo")
        
        # Log configuration and paths
        logger.step(f"Main path: {mainpath}")
        dataid = '784_5'
        logger.step(f"Processing dataset: {dataid}")
        
        # Log directory creation
        results_dir = f'{root_dir}/results/'
        os.makedirs(results_dir, exist_ok=True)
        logger.step(f"Results directory: {results_dir}")
        
        # Log input files
        mask_pm_path = f'{mainpath}/data/sxt_images/{dataid}_wholecell_label.tiff'
        mask_ne_path = f'{mainpath}/data/sxt_images/{dataid}_NC_label.tiff'
        logger.file_in(mask_pm_path)
        logger.file_in(mask_ne_path)
        
        # Log data loading and properties
        logger.step("Loading mask data...")
        mask_data_pm = tifffile.imread(mask_pm_path).astype(int)
        mask_data_ne = tifffile.imread(mask_ne_path).astype(int)
        logger.step(f"PM mask shape: {mask_data_pm.shape}, dtype: {mask_data_pm.dtype}")
        logger.step(f"NE mask shape: {mask_data_ne.shape}, dtype: {mask_data_ne.dtype}")
        
        # Log processing steps
        logger.step("Initializing partitioning processor...")
        logger.step("Extracting boundaries for analysis...")
        logger.step("Creating continuous radial partitions using distance transform...")
        
        # Log results and outputs
        xvg_path = partitioner.save_partition_coords_to_xvg(partition_coords, dataid, results_dir)
        logger.file_out(xvg_path)
        
        comparison_path = f"{results_dir}{dataid}_partition_comparison.png"
        logger.file_out(comparison_path)
        
        logger.step(f"Processing completed for dataset: {dataid}")

Log File Output
---------------

The logging generates timestamped files with structured output:

.. code-block:: text

    2024-01-15 14:30:25,123 - INFO - Starting SXT Partitioning Demo
    2024-01-15 14:30:25,124 - INFO - Main path: /data/sxt_analysis
    2024-01-15 14:30:25,125 - INFO - Processing dataset: 784_5
    2024-01-15 14:30:25,126 - INFO - Results directory: /data/sxt_analysis/results/
    2024-01-15 14:30:25,127 - INFO - Input file: /data/sxt_images/784_5_wholecell_label.tiff
    2024-01-15 14:30:25,128 - INFO - Input file: /data/sxt_images/784_5_NC_label.tiff
    2024-01-15 14:30:25,129 - INFO - Loading mask data...
    2024-01-15 14:30:26,456 - INFO - PM mask shape: (512, 1024, 1024), dtype: int32
    2024-01-15 14:30:26,457 - INFO - NE mask shape: (512, 1024, 1024), dtype: int32
    2024-01-15 14:30:26,458 - INFO - Output file: /data/results/784_5_partition_coords.xvg
    2024-01-15 14:30:27,123 - INFO - Processing completed for dataset: 784_5

