SXT ISG Mask R-CNN Instance Segmentation
========================================

Overview
--------

The SXT ISG Mask R-CNN segmenter provides **instance-level segmentation** of insulin secretory granules (ISGs) in Soft X-ray Tomography data. Unlike the U-Net based approach which produces semantic segmentation, Mask R-CNN enables:

* **Instance Segmentation**: Distinguish individual ISG instances with unique IDs
* **Precise Boundaries**: Pixel-level mask accuracy for each granule
* **Morphological Analysis**: Extract per-instance features (volume, shape, position)
* **Spatial Relationships**: Analyze ISG-ISG and ISG-organelle interactions at instance level

Model Performance
-----------------

**Trained Model**: ``isg_mask_rcnn_final.pth``

**Performance Metrics** (on validation set):

* mAP @ [0.50:0.95]: **0.2276**
* mAP @ 0.50: **0.6137**
* Recall @ 100 detections: **0.367**

**Architecture**:

* Backbone: ResNet-50 with Feature Pyramid Network (FPN)
* Head: Mask R-CNN with box and mask prediction branches
* Input: 2D slices from 3D SXT volumes (256×256 patches)
* Output: Instance masks with confidence scores and bounding boxes
* Classes: 2 (background + ISG)

**Training Configuration**:

* Dataset: COCO-format annotations from SXT ISG training data
* Epochs: 20 (with early stopping)
* Batch size: 4
* Learning rate: 0.005 (with scheduler)
* Device: CUDA (GPU acceleration)
* Evaluation frequency: Every 2 epochs

Quick Start
-----------

Using Unified API (Recommended):

.. code-block:: python

    from ipa.processing.segmentation import create_segmenter
    import mrcfile
    import numpy as np
    
    # Load SXT volume
    with mrcfile.open('data/sxt_images/Stevens_pancreatic_INS_1E_784_5_pre_rec.mrc', permissive=True) as mrc:
        volume = mrc.data.astype(np.float32)
    
    # Create segmenter and load model
    segmenter = create_segmenter(modality='sxt', task='isg_maskrcnn')
    segmenter.load_model()  # Automatically loads isg_mask_rcnn_final.pth
    
    # Predict instance masks
    instance_mask = segmenter.predict(volume, threshold=0.5)
    
    # Analyze results
    unique_ids = np.unique(instance_mask)
    num_instances = len(unique_ids[unique_ids > 0])
    print(f"Detected {num_instances} ISG instances")
    
    # Extract per-instance statistics
    if num_instances > 0:
        volumes = []
        for uid in unique_ids:
            if uid > 0:
                volumes.append(np.sum(instance_mask == uid))
        
        import numpy as np
        volumes = np.array(volumes)
        print(f"Average volume: {np.mean(volumes):.2f} voxels")
        print(f"Median volume:  {np.median(volumes):.2f} voxels")

Direct Class Usage:

.. code-block:: python

    from ipa.processing.segmentation.unified import SXTISGMaskRCNNSegmenter
    
    # Initialize segmenter
    segmenter = SXTISGMaskRCNNSegmenter(device='cuda')
    
    # Load trained model
    segmenter.load_model()
    
    # Predict
    results = segmenter.predict(volume, threshold=0.5)
    
    # results is a labeled 3D mask where each ISG has a unique integer ID
    # Background = 0, ISG instances = 1, 2, 3, ...

API Reference
-------------

SXTISGMaskRCNNSegmenter
^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: ipa.processing.segmentation.unified.SXTISGMaskRCNNSegmenter
   :members: __init__, _load_model_impl, predict
   :undoc-members:
   :show-inheritance:

**predict() Parameters**:

* ``data`` (np.ndarray): Input 3D volume with shape (D, H, W)
* ``threshold`` (float): Confidence threshold for detection (default: 0.5)
* Returns: Labeled 3D mask (np.ndarray, dtype=uint16) where each ISG has unique ID

**Processing Pipeline**:

1. Volume is processed slice-by-slice (2D Mask R-CNN)
2. Each slice is normalized to [0, 1] and converted to 3-channel RGB
3. Mask R-CNN predicts boxes, masks, and scores for each slice
4. High-confidence masks (>threshold) are combined into binary slice
5. Connected component labeling assigns unique IDs within each slice
6. Global instance IDs ensure uniqueness across all slices

Demo Script
-----------

A complete demo script is available at:

``examples/examples_sxt/demo_SXT_isg_maskrcnn.py``

**Run the demo**:

.. code-block:: bash

    cd /media/cuixi/data7/liad/gitspace/iPA
    python3 examples/examples_sxt/demo_SXT_isg_maskrcnn.py

**What the demo does**:

1. Loads real SXT volume (cell 784_5)
2. Initializes Mask R-CNN segmenter
3. Runs instance segmentation on the full 3D volume
4. Calculates per-instance statistics (count, volumes)
5. Generates visualization with 3 panels:
   - Original slice
   - Binary mask (all ISGs)
   - Instance mask (colored by instance ID)
6. Saves results as PNG and TIFF files

**Expected output**:

.. code-block:: text

    ============================================================
    SXT ISG Mask R-CNN Instance Segmentation Demo
    ============================================================
    Test Cell: 784_5
    Image Path: data/sxt_images/sxt_isg_training/for_24_datasets/...
    
    Note: This uses Mask R-CNN (mAP=0.23)
    For better results, use demo_SXT_isg_instance.py (blob_fit)
    
    [1/4] Loading 3D volume...
      Volume shape: (473, 486, 473)
    
    [2/4] Loading Mask R-CNN model...
      Model loaded successfully!
    
    [3/4] Running Mask R-CNN prediction...
    
    [4/4] Analyzing instances...
      Total instances detected: XX
      Average volume: XXX.XX voxels
      Median volume:  XXX.XX voxels
      Max volume:     XXX voxels
      Min volume:     XX voxels
    
    Generating visualization...
    Visualization saved to: sxt_isg_maskrcnn_result_784_5.png
    Instance mask saved to: sxt_isg_maskrcnn_784_5.tiff
    
    ============================================================
    Demo completed!
    ============================================================

Training Details
----------------

**Training Script**: ``examples/examples_sxt/train_SXT_isg_mask_rcnn_model.py``

**Dataset Preparation**:

The model was trained on COCO-format annotations created from SXT ISG ground truth masks. The dataset includes:

* Training set: ~XX patches (256×256)
* Validation set: ~XX patches
* Annotation format: COCO JSON with polygon masks
* Preprocessing: Intensity normalization, patch extraction

**Training Process**:

.. code-block:: bash

    python3 examples/examples_sxt/train_SXT_isg_mask_rcnn_model.py

**Training Configuration**:

* Optimizer: SGD with momentum (0.9)
* Learning rate: 0.005 (step decay scheduler)
* Weight decay: 0.0005
* Gradient clipping: max_norm=0.5
* Mixed precision: Enabled (torch.cuda.amp)
* Checkpoint saving: Best model + final model + every 5 epochs

**Model Files**:

After training, the following models are saved:

* ``isg_mask_rcnn_best.pth`` - Best validation performance
* ``isg_mask_rcnn_final.pth`` - Final epoch model (used by default)
* ``isg_mask_rcnn_ep{N}.pth`` - Checkpoints every 5 epochs

Model Location: ``ipa/processing/segmentation/models/sxt/isg_mask_rcnn_final.pth``

Comparison with Other Methods
------------------------------

**Mask R-CNN vs U-Net + Blob Fit**:

+---------------------+------------------+------------------+
| Feature             | Mask R-CNN       | U-Net + Blob Fit |
+=====================+==================+==================+
| Instance Separation | ✅ Native        | ⚠️ Post-process  |
+---------------------+------------------+------------------+
| mAP                 | 0.23             | Higher           |
+---------------------+------------------+------------------+
| Inference Speed     | Fast (~seconds)  | Slower           |
+---------------------+------------------+------------------+
| Boundary Precision  | Good             | Excellent        |
+---------------------+------------------+------------------+
| Complexity          | Simple API       | Multi-step       |
+---------------------+------------------+------------------+
| Recommended Use     | Quick analysis   | Production       |
+---------------------+------------------+------------------+

**When to use Mask R-CNN**:

* Quick prototyping and exploration
* When you need native instance IDs without post-processing
* For comparison with other instance segmentation methods
* When inference speed is critical

**When to use U-Net + Blob Fit**:

* Production-quality analysis
* Maximum accuracy requirements
* Publication-ready results
* Complex morphological analysis

For production use, we recommend: ``examples/examples_sxt/demo_SXT_isg_instance.py`` (blob_fit method)

See Also
--------

* :doc:`sxt_segmentation` - Complete SXT segmentation documentation
* :doc:`../partitioning` - Radial partitioning for spatial analysis
* :doc:`../../analysis/morphology` - Morphological feature extraction
