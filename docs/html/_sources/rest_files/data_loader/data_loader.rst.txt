Data loader Module
===================

The data loader module provides comprehensive support for loading multiple microscopy file formats commonly used in cellular imaging analysis. It offers a unified interface for handling various commercial and open-source imaging formats with automatic format detection and optional dependency management.

.. autoclass:: ipa.data_loader.UniversalDataLoader
   :members:
   :undoc-members:
   :show-inheritance:

.. .. automodule:: ipa.data_loader.data_loading
..    :members:
..    :undoc-members:
..    :show-inheritance:

Overview
--------

The ``UniversalDataLoader`` class provides a single entry point for loading various microscopy file formats. The module  handles missing optional dependencies by providing clear warnings and fallback options.

Supported File Formats
-----------------------

Core Formats (Always Available)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* **.mrc** - Medical Research Council format for cryo-EM and tomography data
* **.tif/.tiff** - Tagged Image File Format with multi-channel support
* **.npz** - Compressed NumPy array format for processed data
* **.lif** - Leica Image Format (requires ``readlif``)
* **.czi** - Carl Zeiss Image format (requires ``czifile``)  
* **.nd2** - Nikon ND2 format (requires ``nd2reader``)
.. * **Additional formats** - Extended support via ``aicsimageio`` (.oib, .vsi, etc.)


.. Dependency Installation
.. -----------------------

.. Optional dependencies for commercial formats:

.. .. code-block:: bash

..     # Individual packages
..     pip install readlif         # For Leica .lif files
..     pip install czifile         # For Zeiss .czi files  
..     pip install nd2reader       # For Nikon .nd2 files

..     # pip install aicsimageio     # Supports many additional formats



Usage Examples
-----------

Basic Loading
^^^^^^^^^^^^^

.. code-block:: python

    from ipa.data_loader.data_loading import UniversalDataLoader
    
    # Simple data loading
    data = UniversalDataLoader.load_data('experiment.tif')
    print(f"Loaded data shape: {data.shape}")
    
    # Load with metadata
    data = UniversalDataLoader.load_data('experiment.tif', )
    print(f"Format: {metadata['format']}, Shape: {metadata['shape']}")

Multi-channel Data Loading
^^^^^^^^^^^^^^^^^^

.. code-block:: python

    # Load specific channel from multi-channel data
    nuclei_channel = UniversalDataLoader.load_data('cells.nd2', channel=0)
    protein_channel = UniversalDataLoader.load_data('cells.nd2', channel=1)
    
    # Normalize intensity values
    normalized_data = UniversalDataLoader.load_data('cells.nd2', 
                                                   channel=0, 
                                                   normalize=True)



Batch File Loading
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    # Process multiple files
    file_list = [
        'experiment_01.mrc',
        'experiment_02.tif', 
        'experiment_03.czi'
    ]
    
    results = UniversalDataLoader.batch_load(
        file_list, 
        channel=0, 
        normalize=True
    )
    
    # Process results
    for filename, data in results.items():
        if data is not None:
            print(f"Successfully loaded {filename}: {data.shape}")
        else:
            print(f"Failed to load {filename}")

Log 
----------

The data_loader module includes a logging system for tracking analysis workflows:

.. toctree::
   :maxdepth: 2

   log




