SIM/WFM Segmentation
====================

Structured Illumination Microscopy (SIM) and Widefield Microscopy (WFM) segmentation module provides comprehensive tools for super-resolution and conventional fluorescence microscopy data analysis.

Core Functions
--------------

**skeletonize_organelle**

Extracts skeleton structures from SIM microscopy images for morphological analysis of filamentous organelles like actin networks.

.. autofunction:: ipa.processing.segmentation.skeletonize_organelle

**Parameters:**

* ``image_data`` (numpy.ndarray): Input 3D image volume

**Returns:**

* ``numpy.ndarray``: Skeletonized 3D image showing structural connectivity

**segment_sphere_like_organelle**

Detects and segments spherical structures in SIM images using advanced detection algorithms optimized for immunofluorescence spherical granules (ISG).

.. autofunction:: ipa.processing.segmentation.segment_sphere_like_organelle

**Parameters:**

* ``image_data`` (numpy.ndarray): Input 3D image volume
* ``threshold`` (float, optional): Detection threshold (None for automatic Otsu)
* ``min_size`` (int, optional): Minimum sphere size (default: 100)

**Returns:**

* ``tuple``: (labeled_image, number_of_objects)

**segment_cell_shape**

Performs whole-cell segmentation to define cellular boundaries using membrane or cytoplasmic markers.

.. autofunction:: ipa.processing.segmentation.segment_cell_shape

**Parameters:**

* ``image_data`` (numpy.ndarray): Input cell channel image
* ``threshold`` (float, optional): Segmentation threshold (None for automatic)
* ``min_size`` (int, optional): Minimum cell size (default: 1000)

**Returns:**

* ``numpy.ndarray``: Labeled cell regions

**segment_nucleus**

Specifically segments nuclear regions in fluorescence microscopy images with optimized parameters for nuclear staining.

.. autofunction:: ipa.processing.segmentation.segment_nucleus

**Parameters:**

* ``image_data`` (numpy.ndarray): Input fluorescence image with nuclear staining

**Returns:**

* ``numpy.ndarray``: Labeled nuclear regions

Example Usage
-------------

**Actin Filament Skeletonization**

.. code-block:: python

    from ipa.processing.segmentation import skeletonize_organelle
    from ipa.data_loader import UniversalDataLoader
    
    # Load SIM actin data
    actin_data = UniversalDataLoader.load_data("SIM_raw_Actin.tif")
    
    # Extract actin filament skeleton
    skeleton = skeletonize_organelle(actin_data)
    print(f"Skeleton points: {np.sum(skeleton > 0)}")

**ISG Sphere Detection**

.. code-block:: python

    from ipa.processing.segmentation import segment_sphere_like_organelle
    
    # Load ISG immunofluorescence data
    isg_data = UniversalDataLoader.load_data("SIM_raw_ISG.tif")
    
    # Detect spherical granules
    labeled_spheres, num_spheres = segment_sphere_like_organelle(
        isg_data, 
        threshold=None,  # Automatic Otsu threshold
        min_size=100
    )
    print(f"Detected {num_spheres} ISG spheres")

**Nucleus and Cell Segmentation**

.. code-block:: python

    from ipa.processing.segmentation import segment_cell_shape, segment_nucleus
    
    # Cell shape segmentation using actin
    cell_mask = segment_cell_shape(
        actin_data,
        threshold=None,  # Automatic threshold
        min_size=1000
    )
    
    # Nuclear segmentation
    nucleus_data = UniversalDataLoader.load_data("SIM_raw_N.tif")
    labeled_nuclei = segment_nucleus(nucleus_data)
    num_nuclei = np.max(labeled_nuclei)
    print(f"Detected {num_nuclei} nuclei")

