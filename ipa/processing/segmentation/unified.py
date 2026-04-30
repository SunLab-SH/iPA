"""
Unified Segmentation Interface

Provides a simplified, unified API for all segmentation tasks while maintaining
backward compatibility with existing implementations.

Supported Modalities:
- SXT (Soft X-ray Tomography): cell, nucleus, mitochondria
- SIM (Structured Illumination Microscopy): ER, mitochondria, ISG
- WFM (Widefield Microscopy): cell, nucleus, general organelles
- Cryo-ET (Cryo-Electron Tomography): filaments, membranes
"""

import os
import sys
import numpy as np
from typing import Optional, Dict, Any

from .base import BaseSegmenter

# Define model root relative to the project root (ipa directory)
# This ensures paths are robust regardless of where the script is executed from
PACKAGE_ROOT = os.path.dirname(os.path.abspath(__file__))
MODEL_ROOT = os.path.join(PACKAGE_ROOT, 'models')

DEFAULT_SXT_CELL_MODEL = os.path.join(MODEL_ROOT, 'sxt', 'cell_nucleus_unet_best.pth')
DEFAULT_SXT_MITO_MODEL = os.path.join(MODEL_ROOT, 'sxt', 'mito_unet_best.pth')
DEFAULT_SXT_ISG_MODEL = os.path.join(MODEL_ROOT, 'sxt', 'isg_unet_best.pth')
DEFAULT_SXT_ISG_MASKRCNN_MODEL = os.path.join(MODEL_ROOT, 'sxt', 'isg_mask_rcnn_final.pth')
DEFAULT_SIM_ER_MODEL = os.path.join(MODEL_ROOT, 'sim', 'ernet_final.pth')


class SXTCellSegmenter(BaseSegmenter):
    """SXT Cell and Nucleus Segmentation using U-Net."""
    
    def __init__(self, device: Optional[str] = None):
        super().__init__(modality='sxt', task='cell', device=device)
        self.model = None
        
    def _load_model_impl(self, path: str = None, **kwargs):
        """Load SXT cell segmentation model."""
        if path is None:
            path = DEFAULT_SXT_CELL_MODEL
        if not os.path.exists(path):
            raise FileNotFoundError(f"SXT Cell model not found at {path}.")
            
        import torch
        # Use the Unet class that matches the nu_best.pth structure
        from .segmentation_sxt.model_sxt_cell_mito.networks.Unet import Unet
        
        # Initialize model: n_class=3 (bg, mem, nu)
        self.model = Unet(n_class=3, is_dropout=True)
        
        # Load state dict (handle DataParallel 'module.' prefix if present)
        checkpoint = torch.load(path, map_location=self.device, weights_only=False)
        new_state_dict = {}
        for k, v in checkpoint.items():
            name = k.replace('module.', '')  # remove 'module.' prefix
            new_state_dict[name] = v
        
        self.model.load_state_dict(new_state_dict)
        self.model.to(self.device)
        self.model.eval()
        self.model_path = path
        
    def predict(self, data: np.ndarray, **kwargs) -> Dict[str, np.ndarray]:
        """
        Predict cell and nucleus masks.
        
        Args:
            data: Input 3D volume (D, H, W)
            
        Returns:
            Dictionary with 'cell_mask' and 'nucleus_mask'
        """
        if not self.is_loaded or self.model is None:
            raise RuntimeError("Model not loaded!")
        
        import torch
        from skimage.transform import resize
        import cv2
        
        # Normalize input to [0, 255] uint8 as expected by the legacy model
        data_norm = ((data - data.min()) / (data.max() - data.min()) * 255).astype(np.uint8)
        
        D, H, W = data.shape
        pred_volume = np.zeros((D, H, W), dtype=np.uint8)
        
        self.model.eval()
        with torch.no_grad():
            for i in range(D):
                slice_img = data_norm[i]
                # Legacy model expects specific preprocessing: normalize + resize to (288, 480)
                inputslice = cv2.resize(slice_img, (480, 288), interpolation=cv2.INTER_NEAREST)
                
                # Convert to tensor with normalization used in training
                transform_input = (inputslice.astype(np.float32) / 255.0 - 0.456) / 0.224
                input_tensor = torch.from_numpy(transform_input).float().unsqueeze(0).unsqueeze(0).to(self.device)
                
                output = self.model(input_tensor)
                probs = torch.softmax(output, dim=1)
                pred = torch.argmax(probs, dim=1).squeeze().cpu().numpy()
                
                # Resize back to original shape
                pred = cv2.resize(pred.astype(np.uint8), (W, H), interpolation=cv2.INTER_NEAREST)
                    
                pred_volume[i] = pred
        
        # Extract masks: 1 = membrane/cell, 2 = nucleus
        cell_mask = (pred_volume == 1).astype(np.uint8)
        nucleus_mask = (pred_volume == 2).astype(np.uint8)
        
        return {'cell_mask': cell_mask, 'nucleus_mask': nucleus_mask}


class SXTMitoSegmenter(BaseSegmenter):
    """SXT Mitochondria Segmentation using Mask R-CNN.
    
    Implements instance segmentation for mitochondria using a pretrained
    Mask R-CNN model (ResNet50-FPN backbone).
    """
    
    def __init__(self, device: Optional[str] = None):
        super().__init__(modality='sxt', task='mito', device=device)
        self.model = None
        
    def _load_model_impl(self, path: str = None, **kwargs):
        """Load SXT mitochondria U-Net model."""
        if path is None:
            path = DEFAULT_SXT_MITO_MODEL
        if not os.path.exists(path):
            raise FileNotFoundError(f"SXT Mito model not found at {path}.")
            
        import torch
        from .segmentation_sxt.model_sxt_isg.unet.unet_model import UNet
        
        # Initialize model (assuming 3 channels, 2 classes: bg + mito)
        self.model = UNet(n_channels=3, n_classes=2, bilinear=False)
        self.model.load_state_dict(torch.load(path, map_location=self.device, weights_only=False))
        self.model.to(self.device)
        self.model.eval()
        self.model_path = path
        
    def predict(self, data: np.ndarray, **kwargs) -> np.ndarray:
        """
        Predict mitochondria mask using U-Net.
        """
        if not self.is_loaded or self.model is None:
            raise RuntimeError("Model not loaded!")
        
        import torch
        from skimage.transform import resize
        import cv2
        
        # Normalize input to [0, 255] uint8
        data_norm = ((data - data.min()) / (data.max() - data.min()) * 255).astype(np.uint8)
        
        D, H, W = data.shape
        pred_volume = np.zeros((D, H, W), dtype=np.uint8)
        
        self.model.eval()
        with torch.no_grad():
            for i in range(D):
                slice_img = data_norm[i]
                # Preprocess: resize to (480, 288) and normalize
                inputslice = cv2.resize(slice_img, (480, 288), interpolation=cv2.INTER_NEAREST)
                transform_input = (inputslice.astype(np.float32) / 255.0 - 0.456) / 0.224
                # Convert to 3-channel tensor as expected by the model
                input_tensor = torch.from_numpy(np.stack([transform_input]*3, axis=0)).float().unsqueeze(0).to(self.device)
                
                output = self.model(input_tensor)
                probs = torch.softmax(output, dim=1)
                pred = torch.argmax(probs, dim=1).squeeze().cpu().numpy()
                
                # Resize back to original shape
                pred = cv2.resize(pred.astype(np.uint8), (W, H), interpolation=cv2.INTER_NEAREST)
                pred_volume[i] = pred
        
        # Class 1 is mitochondria
        return (pred_volume == 1).astype(np.uint8)


class SXTISGSegmenter(BaseSegmenter):
    """SXT ISG (Insulin Granule) Segmentation using U-Net."""
    
    def __init__(self, device: Optional[str] = None):
        super().__init__(modality='sxt', task='isg', device=device)
        self.model = None
        
    def _load_model_impl(self, path: str = None, **kwargs):
        """Load SXT ISG U-Net model."""
        if path is None:
            path = DEFAULT_SXT_ISG_MODEL
        if not os.path.exists(path):
            raise FileNotFoundError(f"SXT ISG model not found at {path}.")
            
        import torch
        from .segmentation_sxt.model_sxt_isg.unet.unet_model import UNet
        from PIL import Image
        
        # Initialize U-Net model
        self.model = UNet(n_channels=1, n_classes=2, bilinear=False)
        self.model.to(self.device)
        
        # Load weights
        checkpoint = torch.load(path, map_location=self.device, weights_only=False)
        self.model.load_state_dict(checkpoint)
        self.model.eval()
        
        self.PIL_Image = Image  # Store for preprocessing
        self.model_path = path
        
    def predict(self, data: np.ndarray, **kwargs) -> np.ndarray:
        """
        Predict ISG (insulin granule) mask.
        
        Args:
            data: Input 3D volume (D, H, W)
            **kwargs: Additional parameters
            
        Returns:
            ISG segmentation mask (binary)
        """
        if not self.is_loaded or self.model is None:
            raise RuntimeError("Model not loaded! Call load_model() first.")
        
        import torch
        from skimage.transform import resize
        
        D, H, W = data.shape
        pred_volume = np.zeros((D, H, W), dtype=np.uint8)
        
        self.model.eval()
        with torch.no_grad():
            for z in range(D):
                slice_data = data[z]
                
                # 1. Per-slice normalization (0-1)
                v_min, v_max = slice_data.min(), slice_data.max()
                if v_max - v_min > 0:
                    slice_norm = (slice_data - v_min) / (v_max - v_min)
                else:
                    slice_norm = slice_data
                
                # 2. Padding to square
                h, w = slice_norm.shape
                size = max(h, w)
                img_padded = np.zeros((size, size), dtype=np.float32)
                start_h = (size - h) // 2
                start_w = (size - w) // 2
                img_padded[start_h:start_h+h, start_w:start_w+w] = slice_norm
                
                # 3. Resize to 400x400
                pil_img = self.PIL_Image.fromarray((img_padded * 255).astype(np.uint8))
                pil_img = pil_img.resize((400, 400), resample=self.PIL_Image.BILINEAR)
                input_tensor = torch.from_numpy(np.asarray(pil_img) / 255.0).float().unsqueeze(0).unsqueeze(0).to(self.device)
                
                # 4. Predict
                output = self.model(input_tensor)
                pred_slice = torch.argmax(output, dim=1).squeeze().cpu().numpy()
                
                # 5. Resize back to original size
                pred_resized = (resize(pred_slice, (H, W), order=0) > 0.5).astype(np.uint8)
                pred_volume[z] = pred_resized
        
        return pred_volume


class SXTISGMaskRCNNSegmenter(BaseSegmenter):
    """SXT ISG Instance Segmentation using Mask R-CNN.
    
    Implements true instance segmentation where each ISG gets a unique ID.
    Uses a ResNet50-FPN backbone trained on COCO-format patches.
    """
    
    def __init__(self, device: Optional[str] = None):
        super().__init__(modality='sxt', task='isg_maskrcnn', device=device)
        self.model = None
        
    def _load_model_impl(self, path: str = None, **kwargs):
        """Load Mask R-CNN model."""
        if path is None:
            path = DEFAULT_SXT_ISG_MASKRCNN_MODEL
        if not os.path.exists(path):
            raise FileNotFoundError(f"SXT ISG Mask R-CNN model not found at {path}.")
            
        import torch
        import torchvision
        from torchvision.models.detection import maskrcnn_resnet50_fpn
        
        # Initialize model (2 classes: background + ISG)
        self.model = maskrcnn_resnet50_fpn(weights=None, num_classes=2)
        self.model.load_state_dict(torch.load(path, map_location=self.device, weights_only=False))
        self.model.to(self.device)
        self.model.eval()
        self.model_path = path
        
    def predict(self, data: np.ndarray, threshold: float = 0.5, **kwargs) -> np.ndarray:
        """
        Predict ISG instance mask for a 3D volume by processing 2D slices.
        
        Args:
            data: Input 3D volume (D, H, W)
            threshold: Confidence threshold for detection
            
        Returns:
            Labeled 3D mask where each ISG has a unique integer ID
        """
        import torch
        from torchvision.transforms import functional as F
        from scipy import ndimage
        
        D, H, W = data.shape
        # We'll collect all detected instances across slices
        # For true 3D instance consistency, we'd need 3D Mask R-CNN, 
        # but here we use 2D slices and label them.
        final_3d_mask = np.zeros((D, H, W), dtype=np.uint16)
        global_inst_id = 1
        
        self.model.eval()
        with torch.no_grad():
            for z in range(D):
                slice_data = data[z]
                # Normalize to 0-1 and convert to tensor
                v_min, v_max = slice_data.min(), slice_data.max()
                if v_max - v_min > 0:
                    slice_norm = (slice_data - v_min) / (v_max - v_min)
                else:
                    slice_norm = slice_data
                    
                # Convert to RGB tensor (Mask R-CNN expects 3 channels)
                img_tensor = torch.from_numpy(np.stack([slice_norm]*3, axis=0)).float().to(self.device)
                
                # Predict
                predictions = self.model([img_tensor])
                
                # Process predictions for this slice
                boxes = predictions[0]['boxes'].cpu().numpy()
                masks = predictions[0]['masks'].cpu().numpy()[0, 0] # (H, W)
                scores = predictions[0]['scores'].cpu().numpy()
                
                # Filter by threshold
                keep = scores > threshold
                filtered_masks = masks * keep[:, None, None] if len(keep) > 0 else np.zeros_like(masks)
                
                # Label the filtered masks in this slice
                if np.sum(filtered_masks) > 0:
                    # Combine all high-confidence masks into one binary slice
                    binary_slice = (np.max(filtered_masks, axis=0) > 0.5).astype(np.uint8)
                    labeled_slice, n = ndimage.label(binary_slice)
                    
                    # Assign global IDs
                    for inst_id in range(1, n + 1):
                        final_3d_mask[z][labeled_slice == inst_id] = global_inst_id
                        global_inst_id += 1
                        
        return final_3d_mask


class SIMERSegmenter(BaseSegmenter):
    """SIM Endoplasmic Reticulum Segmentation using ERNet."""
    
    def __init__(self, device: Optional[str] = None):
        super().__init__(modality='sim', task='er', device=device)
        
    def _load_model_impl(self, path: str = None, image_size: int = 1000, **kwargs):
        """Load ERNet model."""
        if path is None:
            path = DEFAULT_SIM_ER_MODEL
        if not os.path.exists(path):
            raise FileNotFoundError(f"SIM ER model not found at {path}. Please provide a valid model path.")
            
        import torch
        from types import SimpleNamespace
        from .segmentation_sim_wfm.ERNet.models import GetModel
        from .segmentation_sim_wfm.ERNet.datahandler import toTensor
        
        # Create config for ERNet
        opt = SimpleNamespace()
        opt.model = 'rcan'
        opt.nch_in = 1
        opt.nch_out = 2
        opt.n_resgroups = 5
        opt.n_resblocks = 10
        opt.n_feats = 64
        opt.reduction = 16
        opt.narch = 0
        opt.multigpu = False
        opt.cpu = (self.device.type == 'cpu')
        opt.undomulti = False
        
        # Load model with config
        self.model = GetModel(opt)
        checkpoint = torch.load(path, map_location=self.device, weights_only=False)
        
        # Handle DataParallel wrapper
        from collections import OrderedDict
        new_state_dict = OrderedDict()
        for k, v in checkpoint['state_dict'].items():
            name = k[7:] if k.startswith('module.') else k
            new_state_dict[name] = v
        
        self.model.load_state_dict(new_state_dict)
        self.model.to(self.device)
        self.model.eval()
        
        self.toTensor = toTensor
        self.image_size = image_size
        
    def predict(self, data: np.ndarray, **kwargs) -> np.ndarray:
        """
        Predict ER mask.
        
        Args:
            data: Input image (H, W) or (H, W, C)
            **kwargs: Additional parameters
            
        Returns:
            ER segmentation mask
        """
        import torch
        
        if not self.is_loaded:
            raise RuntimeError("Model not loaded! Call load_model() first.")
        
        # Preprocess
        data_normalized = self._normalize_data(data)
        
        # Convert to tensor and move to device
        input_tensor = self.toTensor(data_normalized)
        input_tensor = input_tensor.unsqueeze(0).to(self.device)
        
        # Predict
        with torch.no_grad():
            output = self.model(input_tensor)
        
        # Postprocess
        mask = output.squeeze().cpu().numpy()
        mask = (mask > 0.5).astype(np.uint8)
        
        return mask


class SIMMitoSegmenter(BaseSegmenter):
    """SIM Mitochondria Segmentation using thresholding-based approach.
    
    Implements a thresholding-based pipeline inspired by Lefebvre et al. (2021):
    1. Background correction (Gaussian blur + subtraction)
    2. Global/adaptive thresholding (Otsu)
    3. Morphological filtering (opening/closing)
    4. Connected-component labeling with size filtering
    """
    
    def __init__(self, device: Optional[str] = None):
        super().__init__(modality='sim', task='mito', device=device)
        # Thresholding-based method, no model needed
        self.is_loaded = True
        
    def predict(self, data: np.ndarray, threshold: Optional[float] = None,
                min_size: int = 50, max_size: Optional[int] = None,
                sigma_bg: float = 20.0, **kwargs) -> np.ndarray:
        """
        Segment mitochondria using thresholding-based pipeline.
        
        Args:
            data: Input image (H, W) or (D, H, W)
            threshold: Segmentation threshold (auto-calculated via Otsu if None)
            min_size: Minimum object size to keep (default: 50 pixels)
            max_size: Maximum object size (None for no limit)
            sigma_bg: Gaussian sigma for background estimation (default: 20)
            **kwargs: Additional parameters
            
        Returns:
            Binary segmentation mask (same shape as input)
        """
        from skimage import filters, morphology, measure
        from scipy import ndimage
        
        # Ensure float type
        data_float = data.astype(np.float32)
        
        # Step 1: Background correction
        # Estimate background with large Gaussian blur and subtract
        background = filters.gaussian(data_float, sigma=sigma_bg)
        corrected = data_float - background
        
        # Normalize to [0, 1]
        corrected = corrected - corrected.min()
        if corrected.max() > 0:
            corrected = corrected / corrected.max()
        
        # Step 2: Thresholding
        if threshold is None:
            # Use Otsu's method for automatic threshold
            threshold = filters.threshold_otsu(corrected)
        
        binary = corrected > threshold
        
        # Step 3: Morphological operations
        # Choose structuring element based on dimensionality
        if data.ndim == 2:
            elem_open = morphology.disk(1)
            elem_close = morphology.disk(2)
        else:
            elem_open = morphology.ball(1)
            elem_close = morphology.ball(2)
        
        # Opening to remove small noise
        binary = morphology.binary_opening(binary, elem_open)
        # Closing to fill small holes and smooth boundaries
        binary = morphology.binary_closing(binary, elem_close)
        
        # Fill holes
        binary = ndimage.binary_fill_holes(binary)
        
        # Step 4: Connected component labeling and size filtering
        labeled = measure.label(binary)
        props = measure.regionprops(labeled)
        
        # Filter by size
        filtered_mask = np.zeros_like(binary, dtype=np.uint8)
        for prop in props:
            area = prop.area
            if area >= min_size:
                if max_size is None or area <= max_size:
                    filtered_mask[labeled == prop.label] = 1
        
        return filtered_mask


class SIMISGSegmenter(BaseSegmenter):
    """SIM ISG Instance Segmentation using intensity-based watershed.
    
    Implements a robust pipeline for detecting individual insulin granules:
    1. Background subtraction (Top-hat filtering)
    2. LoG blob detection for seed points
    3. Marker-controlled Watershed for instance separation
    """
    
    def __init__(self, device: Optional[str] = None):
        super().__init__(modality='sim', task='isg', device=device)
        self.is_loaded = True # No model needed
        
    def predict(self, data: np.ndarray, sigma: float = 2.0, 
                threshold_rel: float = 0.1, min_size: int = 10, **kwargs) -> np.ndarray:
        """
        Segment individual ISGs.
        
        Args:
            data: Input image (H, W) or (D, H, W)
            sigma: Sigma for LoG filter (controls spot size)
            threshold_rel: Relative threshold for blob detection
            min_size: Minimum object size in pixels
            
        Returns:
            Labeled mask where each ISG has a unique integer ID
        """
        from skimage import filters, morphology, measure, segmentation
        from scipy import ndimage
        
        data_float = data.astype(np.float32)
        if data_float.max() > 0:
            data_float = data_float / data_float.max()
            
        # 1. Enhance spots using LoG
        log_filtered = ndimage.gaussian_laplace(data_float, sigma=sigma)
        # Invert because LoG produces dark spots on bright background for blobs
        spots = -log_filtered 
        spots = spots - spots.min()
        if spots.max() > 0:
            spots = spots / spots.max()
            
        # 2. Find markers (seeds) using local maxima
        from skimage.feature import peak_local_max
        threshold = filters.threshold_otsu(spots) * threshold_rel
        coordinates = peak_local_max(spots, min_distance=int(sigma*2), threshold_abs=threshold)
        markers = np.zeros_like(data_float, dtype=bool)
        if len(coordinates) > 0:
            # Handle both 2D and 3D coordinates
            if data_float.ndim == 2:
                markers[coordinates[:, 0], coordinates[:, 1]] = True
            else:
                markers[coordinates[:, 0], coordinates[:, 1], coordinates[:, 2]] = True
        markers = ndimage.label(markers)[0]
        
        # 3. Watershed segmentation
        # Use the original intensity as elevation map
        elevation_map = 1 - data_float 
        labeled_mask = segmentation.watershed(elevation_map, markers, mask=data_float > 0.1)
        
        # 4. Filter by size
        props = measure.regionprops(labeled_mask)
        final_mask = np.zeros_like(labeled_mask, dtype=np.uint16)
        current_id = 1
        for prop in props:
            if prop.area >= min_size:
                final_mask[labeled_mask == prop.label] = current_id
                current_id += 1
                
        return final_mask


class WFMSegmenter(BaseSegmenter):
    """WFM General Organelle Segmentation."""
    
    def __init__(self, task: str = 'cell', device: Optional[str] = None):
        super().__init__(modality='wfm', task=task, device=device)
        
    def predict(self, data: np.ndarray, threshold: Optional[float] = None, 
                min_size: int = 100, **kwargs) -> np.ndarray:
        """
        Segment organelles in WFM images using traditional methods.
        """
        from skimage import filters, morphology, measure
        from scipy import ndimage
        
        # Handle multi-dimensional data (e.g., 5D time series)
        if data.ndim > 3:
            print(f"Warning: Input has {data.ndim} dimensions. Extracting t=0 slice...")
            if data.ndim == 5:
                z_mid = data.shape[1] // 2
                data = data[0, z_mid, 0]
            elif data.ndim == 4:
                data = data[0, 0]
        
        data_normalized = self._normalize_data(data)
        
        # Auto threshold if not provided
        if threshold is None:
            threshold = filters.threshold_otsu(data_normalized)
            
        binary = data_normalized > threshold
        
        # Morphological cleaning
        elem = morphology.ball(1) if data.ndim == 3 else morphology.disk(1)
        binary = morphology.binary_opening(binary, elem)
        binary = ndimage.binary_fill_holes(binary)
        
        # Size filtering
        labeled = measure.label(binary)
        props = measure.regionprops(labeled)
        mask = np.zeros_like(binary, dtype=np.uint8)
        for prop in props:
            if prop.area >= min_size:
                mask[labeled == prop.label] = 1
                
        return mask


class ETFilamentSegmenter(BaseSegmenter):
    """Cryo-ET Filament Segmentation and Skeletonization."""
    
    def __init__(self, device: Optional[str] = None):
        super().__init__(modality='et', task='filament', device=device)
        
    def predict(self, data: np.ndarray, skeletonize: bool = True,
                threshold_multiplier: float = 1.2, **kwargs) -> Dict[str, np.ndarray]:
        """
        Segment and skeletonize filaments in Cryo-ET data.
        
        Args:
            data: Input 3D volume (D, H, W)
            skeletonize: Whether to extract skeleton (default: True)
            threshold_multiplier: Multiplier for automatic threshold (default: 1.2)
            **kwargs: Additional parameters passed to skeletonization_et_segmentation
            
        Returns:
            Dictionary with 'mask' and optionally 'skeleton'
        """
        from .segmentation_et.segment_et import (
            skeletonization_et_segmentation,
            save_filament_branches_json
        )
        
        # Normalize
        data_normalized = self._normalize_data(data)
        
        # Segment using threshold multiplier
        mask = skeletonization_et_segmentation(
            data_normalized,
            threshold_multiplier=threshold_multiplier
        )
        
        result = {'mask': mask.astype(np.uint8)}
        
        # Extract skeleton if requested
        if skeletonize:
            from skimage.morphology import skeletonize_3d
            skeleton = skeletonize_3d(mask > 0)
            result['skeleton'] = skeleton.astype(np.uint8)
        
        return result


class ETMembraneSegmenter(BaseSegmenter):
    """Cryo-ET Membrane Segmentation."""
    
    def __init__(self, device: Optional[str] = None):
        super().__init__(modality='et', task='membrane', device=device)
        
    def predict(self, data: np.ndarray, threshold: Optional[float] = None,
                **kwargs) -> np.ndarray:
        """
        Segment membranes in Cryo-ET data.
        
        Args:
            data: Input 3D volume (D, H, W)
            threshold: Segmentation threshold
            **kwargs: Additional parameters
            
        Returns:
            Membrane segmentation mask
        """
        # Use simple thresholding for now
        # Can be extended with more sophisticated methods
        data_normalized = self._normalize_data(data)
        
        if threshold is None:
            # Auto threshold using Otsu's method
            from skimage.filters import threshold_otsu
            threshold = threshold_otsu(data_normalized)
        
        mask = (data_normalized > threshold).astype(np.uint8)
        
        # Clean up small objects
        from skimage.morphology import remove_small_objects
        mask = remove_small_objects(mask.astype(bool), min_size=100).astype(np.uint8)
        
        return mask


class SXTMitoInstanceSegmenter(BaseSegmenter):
    """SXT Mitochondria Instance Segmentation.
    
    Implements instance segmentation for mitochondria using the same blob separation
    algorithm as ISGs, adapted for tubular/network structures.
    """
    
    def __init__(self, device: Optional[str] = None):
        super().__init__(modality='sxt', task='mito_instance', device=device)
        self.semantic_segmenter = SXTMitoSegmenter(device=device)
        
    def _load_model_impl(self, path: str = None, **kwargs):
        """Load the underlying semantic model."""
        self.semantic_segmenter.load_model(path=path, **kwargs)
        self.is_loaded = True
        
    def predict(self, data: np.ndarray, min_voxels: int = 50, **kwargs) -> np.ndarray:
        """
        Predict mitochondria instance mask using blob separation.
        
        Args:
            data: Input 3D volume (D, H, W)
            min_voxels: Minimum number of voxels to keep an instance (default: 50)
            **kwargs: Additional parameters passed to semantic segmenter
            
        Returns:
            Labeled mask where each mitochondrion has a unique integer ID
        """
        from .segmentation_sxt.blob_separation import blob_fit
        from scipy import ndimage
        
        # 1. Get semantic mask
        semantic_mask = self.semantic_segmenter.predict(data, **kwargs)
        
        # Prepare raw data (normalize to 0-1 range)
        raw_data = data.astype(np.float32)
        v_min, v_max = raw_data.min(), raw_data.max()
        if v_max - v_min > 0:
            raw_data = (raw_data - v_min) / (v_max - v_min)
        else:
            raw_data = np.zeros_like(raw_data)
        
        # 2. Get initial connected components
        labeled, num_features = ndimage.label(semantic_mask)
        if num_features == 0:
            return np.zeros_like(semantic_mask, dtype=np.uint16)
        
        final_mask = np.zeros_like(semantic_mask, dtype=np.uint16)
        current_global_id = 1
        
        print(f"Starting mito instance separation for {num_features} components...")
        
        # 3. Process each component with blob separation
        for local_id in range(1, num_features + 1):
            coords = np.where(labeled == local_id)
            if len(coords[0]) < min_voxels:
                continue
            
            # Extract ROI with padding
            x_min, x_max = np.min(coords[0]), np.max(coords[0])
            y_min, y_max = np.min(coords[1]), np.max(coords[1])
            z_min, z_max = np.min(coords[2]), np.max(coords[2])
            
            pad = 5
            roi = (
                slice(max(0, x_min - pad), min(data.shape[0], x_max + 1 + pad)),
                slice(max(0, y_min - pad), min(data.shape[1], y_max + 1 + pad)),
                slice(max(0, z_min - pad), min(data.shape[2], z_max + 1 + pad))
            )
            
            roi_mask = semantic_mask[roi].astype(np.uint8)
            roi_raw = raw_data[roi]
            
            # Run blob separation (same algorithm as ISG)
            separated_roi = blob_fit(roi_mask, roi_raw, min_r=1.5, theta=0.8, check_=False)
            
            # Map results back to global mask
            non_zero = np.where(separated_roi > 0)
            if len(non_zero[0]) > 0:
                unique_local_ids = np.unique(separated_roi[non_zero])
                for lid in unique_local_ids:
                    if lid == 0: continue
                    local_coords = np.where(separated_roi == lid)
                    gx = local_coords[0] + roi[0].start
                    gy = local_coords[1] + roi[1].start
                    gz = local_coords[2] + roi[2].start
                    final_mask[gx, gy, gz] = current_global_id
                    current_global_id += 1
        
        print(f"Mito instance separation complete. Found {current_global_id - 1} individual mitochondria.")
        return final_mask


class SXTISGInstanceSegmenter(BaseSegmenter):
    """SXT ISG Instance Segmentation using Organelle Separation.
    
    Implements instance segmentation for ISGs by performing connected component
    analysis on the semantic mask, with optional size filtering to remove noise.
    Uses a simplified blob separation algorithm based on intensity and morphology.
    """
    
    def __init__(self, device: Optional[str] = None):
        super().__init__(modality='sxt', task='isg_instance', device=device)
        self.semantic_segmenter = SXTISGSegmenter(device=device)
        
    def _load_model_impl(self, path: str = None, **kwargs):
        """Load the underlying semantic model."""
        self.semantic_segmenter.load_model(path=path, **kwargs)
        self.is_loaded = True
        
    def predict(self, data: np.ndarray, min_voxels: int = 8, use_advanced_separation: bool = True, **kwargs) -> np.ndarray:
        """
        Predict ISG instance mask.
        
        Args:
            data: Input 3D volume (D, H, W)
            min_voxels: Minimum number of voxels to keep an instance (default: 8)
            use_advanced_separation: Whether to use the advanced blob separation algorithm
            **kwargs: Additional parameters passed to semantic segmenter
            
        Returns:
            Labeled mask where each ISG has a unique integer ID
        """
        from scipy import ndimage
        from skimage.feature import peak_local_max
        from skimage.segmentation import watershed
        
        # 1. Get semantic mask
        semantic_mask = self.semantic_segmenter.predict(data, **kwargs)
        
        if not use_advanced_separation:
            # Simple connected components approach
            labeled, num_features = ndimage.label(semantic_mask)
            if num_features == 0:
                return np.zeros_like(semantic_mask, dtype=np.uint16)
                
            component_sizes = np.bincount(labeled.ravel())
            small_components = component_sizes < min_voxels
            small_components[0] = False  # Keep background
            
            filtered_labeled = labeled.copy()
            for label_id in np.where(small_components)[0]:
                filtered_labeled[labeled == label_id] = 0
                
            # Re-label to ensure consecutive IDs starting from 1
            final_mask = np.zeros_like(filtered_labeled, dtype=np.uint16)
            unique_labels = np.unique(filtered_labeled)
            unique_labels = unique_labels[unique_labels > 0]
            
            for new_id, old_id in enumerate(unique_labels, start=1):
                final_mask[filtered_labeled == old_id] = new_id
            return final_mask
        else:
            # Advanced blob separation using the integrated professional algorithm
            from .segmentation_sxt.blob_separation import blob_fit
            
            # Prepare raw data (normalize to 0-1 range)
            raw_data = data.astype(np.float32)
            v_min, v_max = raw_data.min(), raw_data.max()
            if v_max - v_min > 0:
                raw_data = (raw_data - v_min) / (v_max - v_min)
            else:
                raw_data = np.zeros_like(raw_data)
                
            # Get initial connected components to process individually
            labeled, num_features = ndimage.label(semantic_mask)
            if num_features == 0:
                return np.zeros_like(semantic_mask, dtype=np.uint16)
            
            final_mask = np.zeros_like(semantic_mask, dtype=np.uint16)
            current_global_id = 1
            
            print(f"Starting advanced blob separation for {num_features} components...")
            
            # Process each component
            for local_id in range(1, num_features + 1):
                coords = np.where(labeled == local_id)
                if len(coords[0]) < min_voxels:
                    continue
                    
                # Extract ROI with padding
                x_min, x_max = np.min(coords[0]), np.max(coords[0])
                y_min, y_max = np.min(coords[1]), np.max(coords[1])
                z_min, z_max = np.min(coords[2]), np.max(coords[2])
                
                pad = 5
                roi = (
                    slice(max(0, x_min - pad), min(data.shape[0], x_max + 1 + pad)),
                    slice(max(0, y_min - pad), min(data.shape[1], y_max + 1 + pad)),
                    slice(max(0, z_min - pad), min(data.shape[2], z_max + 1 + pad))
                )
                
                roi_mask = semantic_mask[roi].astype(np.uint8)
                roi_raw = raw_data[roi]
                
                # Run the professional blob separation algorithm
                separated_roi = blob_fit(roi_mask, roi_raw, min_r=1.5, theta=0.8, check_=False)
                
                # Map results back to global mask
                non_zero = np.where(separated_roi > 0)
                if len(non_zero[0]) > 0:
                    unique_local_ids = np.unique(separated_roi[non_zero])
                    for lid in unique_local_ids:
                        if lid == 0: continue
                        local_coords = np.where(separated_roi == lid)
                        gx = local_coords[0] + roi[0].start
                        gy = local_coords[1] + roi[1].start
                        gz = local_coords[2] + roi[2].start
                        final_mask[gx, gy, gz] = current_global_id
                        current_global_id += 1
                    
            print(f"Advanced separation complete. Found {current_global_id - 1} instances.")
            return final_mask


# Factory function for easy access
def create_segmenter(modality: str, task: str, device: Optional[str] = None) -> BaseSegmenter:
    """
    Create a segmenter instance based on modality and task.
    
    Args:
        modality: 'sxt', 'sim', 'wfm', or 'et'
        task: Segmentation task (see table below)
        device: 'cuda' or 'cpu'
        
    Returns:
        Segmenter instance
        
    Supported Modality-Task Combinations:
    
    | Modality | Task         | Description                    |
    |----------|--------------|--------------------------------|
    | sxt      | cell         | Cell and nucleus segmentation  |
    | sxt      | mito         | Mitochondria segmentation      |
    | sxt      | isg          | ISG (insulin granule) seg      |
    | sim      | er           | ER segmentation (ERNet)        |
    | sim      | mito         | Mitochondria segmentation      |
    | wfm      | cell         | Cell shape segmentation        |
    | wfm      | nucleus      | Nucleus segmentation           |
    | wfm      | sphere       | Spherical organelles           |
    | wfm      | vesicle      | Vesicle segmentation           |
    | et       | filament     | Filament skeletonization       |
    | et       | membrane     | Membrane segmentation          |
        
    Example:
        >>> # SXT cell segmentation
        >>> segmenter = create_segmenter('sxt', 'cell')
        >>> segmenter.load_model('model.pth')
        >>> mask = segmenter.predict(data)
        >>> 
        >>> # WFM nucleus segmentation (no model needed)
        >>> segmenter = create_segmenter('wfm', 'nucleus')
        >>> mask = segmenter.predict(data)
    """
    key = f"{modality.lower()}_{task.lower()}"
    
    segmenters = {
        # SXT
        'sxt_cell': SXTCellSegmenter,
        'sxt_mito': SXTMitoSegmenter,
        'sxt_isg': SXTISGSegmenter,
        'sxt_isg_instance': SXTISGInstanceSegmenter,
        'sxt_mito_instance': SXTMitoInstanceSegmenter,
        'sxt_isg_maskrcnn': SXTISGMaskRCNNSegmenter,
        
        # SIM
        'sim_er': SIMERSegmenter,
        'sim_mito': SIMMitoSegmenter,
        'sim_isg': SIMISGSegmenter,
        
        # WFM
        'wfm_cell': lambda device=None: WFMSegmenter(task='cell', device=device),
        'wfm_nucleus': lambda device=None: WFMSegmenter(task='nucleus', device=device),
        'wfm_sphere': lambda device=None: WFMSegmenter(task='sphere', device=device),
        'wfm_vesicle': lambda device=None: WFMSegmenter(task='vesicle', device=device),
        
        # Cryo-ET
        'et_filament': ETFilamentSegmenter,
        'et_membrane': ETMembraneSegmenter,
    }
    
    if key not in segmenters:
        available = list(segmenters.keys())
        raise ValueError(
            f"Unsupported modality/task combination: '{modality}/{task}'.\n"
            f"Key: '{key}'\n"
            f"Available combinations: {available}\n\n"
            f"Example usage:\n"
            f"  create_segmenter('sxt', 'cell')\n"
            f"  create_segmenter('sim', 'er')\n"
            f"  create_segmenter('wfm', 'nucleus')\n"
            f"  create_segmenter('et', 'filament')"
        )
    
    # Handle lambda functions for parameterized segmenters
    creator = segmenters[key]
    if callable(creator) and not hasattr(creator, '__class__'):
        return creator(device)
    else:
        return creator(device=device)
