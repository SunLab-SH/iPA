"""
Universal Data Loading Module for Multiple Imaging Formats

This module provides comprehensive support for loading various microscopy file formats
commonly used in cellular imaging analysis, including both open-source and commercial formats.

The module offers:
- Unified interface for multiple file formats
- Automatic format detection based on file extensions
- Optional dependency management with graceful fallbacks
- Channel selection for multi-channel data
- Intensity normalization capabilities
- Comprehensive metadata extraction
- Batch processing support
- Robust error handling

Supported formats include .mrc, .tif/.tiff, .lif, .czi, .nd2, and others
via the aicsimageio library.
"""

import os
import json
import warnings
import gzip
import tempfile
from typing import Union, Tuple, Dict, Optional, Any
import numpy as np

# Core imaging libraries
import mrcfile
import tifffile

# Optional libraries for commercial formats
try:
    from readlif.reader import LifFile
    LIF_AVAILABLE = True
except ImportError:
    LIF_AVAILABLE = False
    warnings.warn("readlif not available. .lif files cannot be loaded.")

try:
    import czifile
    CZI_AVAILABLE = True
except ImportError:
    CZI_AVAILABLE = False
    warnings.warn("czifile not available. .czi files cannot be loaded.")

try:
    import nd2reader
    ND2_AVAILABLE = True
except ImportError:
    ND2_AVAILABLE = False
    warnings.warn("nd2reader not available. .nd2 files cannot be loaded.")

# try:
#     from aicsimageio import AICSImage
#     AICS_AVAILABLE = True
# except ImportError:
#     AICS_AVAILABLE = False
#     warnings.warn("aicsimageio not available. Advanced format support limited.")


class UniversalDataLoader:
    """
    Universal data loader for multiple microscopy file formats.
    
    Supports .mrc, .tif/.tiff, .lif, .czi, .nd2 and other formats
    with automatic format detection and optional dependency handling.
    """
    
    def __init__(self):
        self._last_metadata = {}
        self._index = self._load_index()

    @staticmethod
    def _load_index():
        """Load the dataset index for simplified data access."""
        try:
            # Try to find index in the project's data root
            script_dir = os.path.dirname(os.path.abspath(__file__))
            project_root = os.path.abspath(os.path.join(script_dir, '..', '..'))
            index_path = os.path.join(project_root, 'data', 'dataset_index.json')
            if os.path.exists(index_path):
                with open(index_path, 'r') as f:
                    return json.load(f)
        except Exception:
            pass
        return {}
    
    @staticmethod
    def load_data(filepath: str, 
                  channel: Optional[int] = None,
                  normalize: bool = False) -> np.ndarray:
        """
        Load data from various microscopy file formats.
        
        Args:
            filepath: Path to the data file (supports absolute paths, relative paths, 
                      or indexed keys like 'sxt/784_5/raw')
            channel: Channel number for multi-channel data (0-indexed)
            normalize: Whether to normalize values to [0,1] range
            
        Returns:
            Data array
        """
        loader = UniversalDataLoader()
        return loader._load_with_metadata(filepath, channel, normalize)
    
    def _load_with_metadata(self, filepath: str, channel: Optional[int], normalize: bool) -> np.ndarray:
        """Internal method that loads data and stores metadata."""
        # Check if filepath is an indexed key
        if '/' in filepath and not os.path.exists(filepath):
            parts = filepath.split('/')
            current = self._index
            try:
                for part in parts:
                    current = current[part]
                if isinstance(current, str):
                    # Resolve relative path from project root
                    script_dir = os.path.dirname(os.path.abspath(__file__))
                    project_root = os.path.abspath(os.path.join(script_dir, '..', '..'))
                    filepath = os.path.join(project_root, 'data', current)
            except (KeyError, TypeError):
                pass  # Fall back to treating as normal path

        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")
        
        # Get file extension
        _, ext = os.path.splitext(filepath.lower())
        
        # Load data based on file extension
        if ext == '.mrc':
            data = self._load_mrc(filepath)
        elif ext in ['.tif', '.tiff']:
            data = self._load_tiff(filepath, channel)
        elif ext == '.lif':
            data = self._load_lif(filepath, channel)
        elif ext == '.czi':
            data = self._load_czi(filepath, channel)
        elif ext == '.nd2':
            data = self._load_nd2(filepath, channel)
        elif ext == '.npz':
            data = self._load_npz(filepath)
        elif ext == '.json':
            data = self._load_json(filepath)
        else:
            # Try using aicsimageio as fallback
            if AICS_AVAILABLE:
                data = self._load_aics(filepath, channel)
            else:
                raise ValueError(f"Unsupported file format: {ext}")
        
        # Normalize if requested
        if normalize and data.dtype in [np.uint8, np.uint16, np.uint32, np.float32, np.float64]:
            if data.dtype in [np.uint8, np.uint16, np.uint32]:
                data = data.astype(np.float32)
            data = (data - data.min()) / (data.max() - data.min())
        
        return data
    
    def get_last_metadata(self) -> Dict:
        """Get metadata from the last loaded file."""
        return self._last_metadata.copy()
    
    @staticmethod
    def _is_gzip_file(filepath: str) -> bool:
        """Check if a file is gzip compressed by reading magic bytes."""
        try:
            with open(filepath, 'rb') as f:
                magic = f.read(2)
                return magic == b'\x1f\x8b'
        except Exception:
            return False
    
    @staticmethod
    def _decompress_gzip(filepath: str) -> str:
        """Decompress a gzip file to a temporary file and return the temp path."""
        temp_fd, temp_path = tempfile.mkstemp(suffix=os.path.splitext(filepath)[1])
        try:
            with gzip.open(filepath, 'rb') as f_in:
                with os.fdopen(temp_fd, 'wb') as f_out:
                    import shutil
                    shutil.copyfileobj(f_in, f_out)
            return temp_path
        except Exception:
            os.close(temp_fd)
            if os.path.exists(temp_path):
                os.remove(temp_path)
            raise
    
    def _load_mrc(self, filepath: str) -> np.ndarray:
        """Load MRC format file (cryo-EM/tomography data).
        
        Automatically detects and decompresses gzip-compressed MRC files.
        """
        temp_file = None
        try:
            # Check if file is gzip compressed
            if self._is_gzip_file(filepath):
                temp_file = self._decompress_gzip(filepath)
                load_path = temp_file
            else:
                load_path = filepath
            
            with mrcfile.open(load_path, permissive=True) as mrc:
                data = mrc.data
                if data is None:
                    raise ValueError("MRC file contains no data")
                
                self._last_metadata = {
                    'format': 'MRC',
                    'shape': data.shape,
                    'dtype': str(data.dtype),
                    'voxel_size': getattr(mrc, 'voxel_size', None),
                    'header': str(mrc.header) if hasattr(mrc, 'header') and mrc.header is not None else {},
                    'compressed': temp_file is not None
                }
            return data
        except Exception as e:
            raise IOError(f"Failed to load MRC file {filepath}: {str(e)}")
        finally:
            # Clean up temporary file
            if temp_file and os.path.exists(temp_file):
                try:
                    os.remove(temp_file)
                except Exception:
                    pass
    
    def _load_tiff(self, filepath: str, channel: Optional[int] = None) -> np.ndarray:
        """Load TIFF format file with optional channel selection."""
        try:
            data = tifffile.imread(filepath)
            if not isinstance(data, np.ndarray):
                data = np.array(data)
            
            original_shape = data.shape
            
            # Handle channel selection for multi-channel data
            if channel is not None and data.ndim > 2:
                if data.ndim == 3 and data.shape[0] < data.shape[1] and data.shape[0] < data.shape[2]:
                    # Assume first dimension is channel
                    if channel < data.shape[0]:
                        data = data[channel]
                    else:
                        raise ValueError(f"Channel {channel} not available. Data has {data.shape[0]} channels.")
                elif data.ndim == 4:
                    # Assume TCZYX or TCYX format
                    if channel < data.shape[1]:
                        data = data[:, channel, ...]
                    else:
                        raise ValueError(f"Channel {channel} not available.")
            
            self._last_metadata = {
                'format': 'TIFF',
                'shape': data.shape,
                'dtype': str(data.dtype),
                'original_shape': original_shape
            }
            return data
        except Exception as e:
            raise IOError(f"Failed to load TIFF file {filepath}: {str(e)}")
    
    def _load_lif(self, filepath: str, channel: Optional[int] = None) -> np.ndarray:
        """Load Leica LIF format file."""
        if not LIF_AVAILABLE:
            raise ImportError("readlif package required for .lif files. Install with: pip install readlif")
        
        try:
            from readlif.reader import LifFile
            lif = LifFile(filepath)
            
            # Get first image if multiple images exist
            img_list = list(lif.get_iter_image())
            if not img_list:
                raise ValueError("No images found in LIF file")
            
            img = img_list[0]  # Take first image
            data = np.array(img)
            
            # Handle channel selection
            if channel is not None and data.ndim > 2:
                if data.ndim == 3:  # CZY or CYX
                    if channel < data.shape[0]:
                        data = data[channel]
                    else:
                        raise ValueError(f"Channel {channel} not available. Data has {data.shape[0]} channels.")
            
            self._last_metadata = {
                'format': 'LIF',
                'shape': data.shape,
                'dtype': str(data.dtype),
                'n_images': len(img_list),
                'channels': getattr(img, 'channels', None),
                'info': getattr(img, 'info', {})
            }
            return data
        except Exception as e:
            raise IOError(f"Failed to load LIF file {filepath}: {str(e)}")
    
    def _load_czi(self, filepath: str, channel: Optional[int] = None) -> np.ndarray:
        """Load Zeiss CZI format file."""
        if not CZI_AVAILABLE:
            raise ImportError("czifile package required for .czi files. Install with: pip install czifile")
        
        try:
            import czifile
            with czifile.CziFile(filepath) as czi:
                data = czi.asarray()
                
                # CZI data often has complex dimensionality (STCZYX, etc.)
                # Squeeze singleton dimensions
                data = np.squeeze(data)
                
                # Handle channel selection
                if channel is not None and data.ndim > 2:
                    # This is simplified - CZI format can be complex
                    # May need more sophisticated dimension parsing
                    if data.ndim >= 3:
                        data = data[channel] if channel < data.shape[0] else data[0]
                
                self._last_metadata = {
                    'format': 'CZI',
                    'shape': data.shape,
                    'dtype': str(data.dtype),
                    'metadata': getattr(czi, 'metadata', None)
                }
            return data
        except Exception as e:
            raise IOError(f"Failed to load CZI file {filepath}: {str(e)}")
    
    def _load_nd2(self, filepath: str, channel: Optional[int] = None) -> np.ndarray:
        """Load Nikon ND2 format file."""
        if not ND2_AVAILABLE:
            raise ImportError("nd2reader package required for .nd2 files. Install with: pip install nd2reader")
        
        try:
            import nd2reader
            with nd2reader.ND2Reader(filepath) as nd2:
                if channel is not None:
                    nd2.default_coords['c'] = channel
                
                data = np.array(nd2)
                
                self._last_metadata = {
                    'format': 'ND2',
                    'shape': data.shape,
                    'dtype': str(data.dtype),
                    'metadata': nd2.metadata,
                    'sizes': getattr(nd2, 'sizes', {}),
                    'axes': getattr(nd2, 'axes', [])
                }
            return data
        except Exception as e:
            raise IOError(f"Failed to load ND2 file {filepath}: {str(e)}")
    
    def _load_npz(self, filepath: str) -> np.ndarray:
        """Load compressed NumPy array file."""
        try:
            npz_data = np.load(filepath)
            
            # If single array, return it directly
            if len(npz_data.files) == 1:
                data = npz_data[npz_data.files[0]]
            else:
                # For multiple arrays, return the first one
                first_key = npz_data.files[0]
                data = npz_data[first_key]
            
            self._last_metadata = {
                'format': 'NPZ',
                'files': npz_data.files,
                'shape': data.shape,
                'dtype': str(data.dtype)
            }
            return data
        except Exception as e:
            raise IOError(f"Failed to load NPZ file {filepath}: {str(e)}")
    
    def _load_json(self, filepath: str) -> Any:
        """Load JSON coordinate/metadata file."""
        try:
            import json
            with open(filepath, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            self._last_metadata = {
                'format': 'JSON',
                'type': type(data).__name__,
                'size': len(data) if hasattr(data, '__len__') else 'unknown'
            }
            return data
        except Exception as e:
            raise IOError(f"Failed to load JSON file {filepath}: {str(e)}")
    
    def _load_aics(self, filepath: str, channel: Optional[int] = None) -> np.ndarray:
        """Load file using aicsimageio as fallback loader."""
        if not AICS_AVAILABLE:
            raise ImportError("aicsimageio package required. Install with: pip install aicsimageio")
        
        try:
            from aicsimageio import AICSImage
            img = AICSImage(filepath)
            data = img.data
            
            # Handle channel selection
            if channel is not None and img.dims.C > 1:
                data = img.get_image_data("ZYX", C=channel)
            else:
                data = np.squeeze(data)
            
            self._last_metadata = {
                'format': f'AICS_{img.reader.__class__.__name__}',
                'shape': data.shape,
                'dtype': str(data.dtype),
                'dims': str(img.dims),
                'physical_pixel_sizes': img.physical_pixel_sizes._asdict(),
                'channel_names': img.channel_names,
                'metadata': img.metadata
            }
            return data
        except Exception as e:
            raise IOError(f"Failed to load file with aicsimageio {filepath}: {str(e)}")
    
    @staticmethod
    def get_supported_formats() -> Dict[str, str]:
        """
        Get supported file formats.
        
        Returns:
            Dictionary mapping file extensions to descriptions
        """
        formats = {
            '.mrc': 'MRC - Cryo-EM and tomography data',
            '.tif/.tiff': 'TIFF - Standard image format',
            '.npz': 'NPZ - Compressed numpy arrays',
            '.json': 'JSON - Coordinate and metadata'
        }
        
        if LIF_AVAILABLE:
            formats['.lif'] = 'LIF - Leica Image Format'
        if CZI_AVAILABLE:
            formats['.czi'] = 'CZI - Zeiss Image Format'
        if ND2_AVAILABLE:
            formats['.nd2'] = 'ND2 - Nikon Image Format'
        if AICS_AVAILABLE:
            formats['others'] = 'Additional formats via aicsimageio'
        
        return formats
    
    @staticmethod
    def batch_load(file_list: list, 
                   channel: Optional[int] = None,
                   normalize: bool = False) -> Dict[str, np.ndarray]:
        """
        Load multiple files in batch.
        
        Args:
            file_list: List of file paths to load
            channel: Channel selection for all files
            normalize: Whether to normalize all data
            
        Returns:
            Dictionary mapping filenames to data arrays
        """
        results = {}
        for filepath in file_list:
            try:
                filename = os.path.basename(filepath)
                data = UniversalDataLoader.load_data(filepath, channel=channel, normalize=normalize)
                results[filename] = data
                print(f"✓ Loaded {filename}: {data.shape}")
            except Exception as e:
                print(f"✗ Failed to load {filepath}: {str(e)}")
                results[os.path.basename(filepath)] = None
        
        return results


