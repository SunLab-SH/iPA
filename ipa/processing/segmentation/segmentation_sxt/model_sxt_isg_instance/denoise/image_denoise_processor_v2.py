import os
import sys
from pathlib import Path
import numpy as np
from PIL import Image
import shutil
from tqdm import tqdm
import logging
import cv2

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def setup_n2v():
    """
    Setup N2V environment
    """
    try:
        from n2v.models import N2V
        logger.info("N2V module imported successfully")
        return True
    except ImportError as e:
        logger.error(f"Cannot import N2V: {e}")
        logger.info("Please ensure noise2void package is installed: pip install noise2void")
        return False

class ImageDenoiser:
    def __init__(self, input_dir, output_dir, model_name='2denoise_190920_071249', basedir=None):
        """
        Initialize image denoiser
        
        Args:
            input_dir: Input directory path
            output_dir: Output directory path
            model_name: Model name
            basedir: Model base directory
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.model_name = model_name
        
        # Set basedir to the same directory as this script if not provided
        if basedir is None:
            self.basedir = Path(__file__).parent
        else:
            self.basedir = basedir
            
        self.model = None
        self.use_simple_denoising = False
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Supported image formats
        self.image_extensions = {'.png', '.jpg', '.jpeg', '.tif', '.tiff'}
        
    def load_model(self):
        """
        Load N2V model - following predictn2v.py approach
        """
        try:
            from n2v.models import N2V
            import json
            
            # Try to load model following predictn2v.py approach
            model_path = Path(self.basedir) / self.model_name
            
            # Check if model directory exists
            if not model_path.exists():
                logger.warning(f"Model directory does not exist: {model_path}")
                self.use_simple_denoising = True
                return True
            
            # Check if required files exist
            config_file = model_path / "config.json"
            weights_file = model_path / "weights_last.h5"
            
            if not config_file.exists() or not weights_file.exists():
                logger.warning(f"Model files missing in {model_path}")
                self.use_simple_denoising = True
                return True
            
            try:
                # Try to load the model with version compatibility handling
                self.model = N2V(config=None, name=self.model_name, basedir=str(self.basedir))
                logger.info(f"Successfully loaded model: {self.model_name} from {self.basedir}")
                return True
                
            except Exception as model_error:
                logger.warning(f"Model loading failed due to version incompatibility: {model_error}")
                
                # Try to load config and check version compatibility
                try:
                    with open(config_file, 'r') as f:
                        config_data = json.load(f)
                    logger.info(f"Model config loaded, but incompatible with current N2V version")
                    logger.info("Falling back to simple denoising method")
                except:
                    logger.warning("Cannot read model config file")
                
                self.use_simple_denoising = True
                return True
                
        except Exception as e:
            logger.error(f"Model setup failed: {e}")
            logger.info("Will use simple denoising method")
            self.use_simple_denoising = True
            return True
    
    def simple_denoise(self, image_array):
        """
        Simple denoising method - using multiple approaches for better results
        """
        try:
            from scipy import ndimage
            from scipy.signal import medfilt2d
            
            # Multi-step denoising approach
            # Step 1: Median filter to remove salt-and-pepper noise
            denoised = medfilt2d(image_array, kernel_size=3)
            
            # Step 2: Gaussian filter for smoothing
            sigma = 0.8
            denoised = ndimage.gaussian_filter(denoised, sigma=sigma)
            
            # Step 3: Edge-preserving bilateral-like filtering using scipy
            # Apply another round of gentle gaussian filtering
            denoised = ndimage.gaussian_filter(denoised, sigma=0.5)
            
            # logger.info("Applied enhanced simple denoising (median + gaussian)")
            return denoised
            
        except ImportError:
            # If scipy is not available, use OpenCV
            try:
                import cv2
                # Convert to proper data type
                img_uint8 = (image_array * 255).astype(np.uint8) if image_array.max() <= 1.0 else image_array.astype(np.uint8)
                
                # Apply median filter
                denoised = cv2.medianBlur(img_uint8, 3)
                # Apply Gaussian filter
                denoised = cv2.GaussianBlur(denoised, (5, 5), 0.8)
                # Apply bilateral filter for edge preservation
                denoised = cv2.bilateralFilter(denoised, 9, 75, 75)
                
                result = denoised.astype(np.float32) / 255.0 if image_array.max() <= 1.0 else denoised.astype(np.float32)
                logger.info("Applied enhanced OpenCV denoising (median + gaussian + bilateral)")
                return result
                
            except:
                logger.warning("Cannot use advanced denoising, returning original image")
                return image_array
    
    def preprocess_image(self, image_path):
        """
        Preprocess image - following predictn2v.py approach
        """
        try:
            # Use cv2 to read image (grayscale)
            img_array = cv2.imread(str(image_path), 0)
            
            if img_array is None:
                # If cv2 fails, try PIL
                img = Image.open(image_path)
                img_array = np.array(img)
                if len(img_array.shape) == 3:
                    img_array = np.mean(img_array, axis=2)
            
            # Convert to float32 and normalize to [0, 1]
            img_array = img_array.astype(np.float32)
            if img_array.max() > 1.0:
                img_array = img_array / 255.0
            
            return img_array
            
        except Exception as e:
            logger.error(f"Image preprocessing failed {image_path}: {e}")
            return None
    
    def denoise_image(self, image_array):
        """
        Denoise image - following predictn2v.py approach
        """
        try:
            if self.use_simple_denoising or self.model is None:
                return self.simple_denoise(image_array)
            
            # Use N2V model for prediction - following predictn2v.py
            pred = self.model.predict(image_array, axes='YX', n_tiles=(1, 1))
            
            return pred
            
        except Exception as e:
            logger.error(f"N2V denoising failed: {e}, using simple denoising method")
            return self.simple_denoise(image_array)
    
    def save_image(self, image_array, output_path):
        """
        Save processed image - following predictn2v.py approach
        """
        try:
            # Ensure output directory exists
            output_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Convert image array to appropriate range
            if image_array.max() <= 1.0:
                # If in [0,1] range, convert to [0,255]
                image_uint8 = (np.clip(image_array, 0.0, 1.0) * 255).astype(np.uint8)
            else:
                # If already in [0,255] range
                image_uint8 = np.clip(image_array, 0, 255).astype(np.uint8)
            
            # Save using cv2
            cv2.imwrite(str(output_path), image_uint8)
            
        except Exception as e:
            logger.error(f"Image saving failed {output_path}: {e}")
    
    def process_folder(self, folder_name):
        """
        Process a single data folder
        """
        input_folder = self.input_dir / folder_name
        output_folder = self.output_dir / folder_name
        
        if not input_folder.exists():
            logger.warning(f"Input folder does not exist: {input_folder}")
            return
        
        # Create output folder
        output_folder.mkdir(parents=True, exist_ok=True)
        
        # Copy label file
        label_file = input_folder / f"combine_{folder_name}_label.tif"
        if label_file.exists():
            shutil.copy2(label_file, output_folder / label_file.name)
            logger.info(f"Copied label file: {label_file.name}")
        
        # Process slices data folder
        slices_input_dir = input_folder
        slices_output_dir = output_folder / "slices data"
        
        if not slices_input_dir.exists():
            logger.warning(f"Slices data folder does not exist: {slices_input_dir}")
            return
        
        slices_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Get all image files
        image_files = [f for f in slices_input_dir.iterdir() 
                      if f.is_file() and f.suffix.lower() in self.image_extensions]
        
        logger.info(f"Processing folder {folder_name}: found {len(image_files)} image files")
        
        # Process each image file
        for img_file in tqdm(image_files, desc=f"Processing {folder_name}"):
            try:
                # Preprocess image
                img_array = self.preprocess_image(img_file)
                if img_array is None:
                    continue
                
                # Denoise
                denoised_img = self.denoise_image(img_array)
                if denoised_img is None:
                    continue
                
                # Save result
                output_path = slices_output_dir / img_file.name
                self.save_image(denoised_img, output_path)
                
            except Exception as e:
                logger.error(f"Image processing failed {img_file}: {e}")
                continue
        
        logger.info(f"Completed processing folder: {folder_name}")
    
    def process_all_folders(self):
        """
        Process all data folders
        """
        if not self.load_model():
            logger.error("Model loading failed, cannot continue processing")
            return
        
        # Get all data folders
        data_folders = [f for f in self.input_dir.iterdir() 
                       if f.is_dir() and '_' in f.name and not f.name.startswith('.')]
        
        logger.info(f"Found {len(data_folders)} data folders")
        
        for folder in sorted(data_folders):
            logger.info(f"Starting to process folder: {folder.name}")
            self.process_folder(folder.name)
        
        logger.info("All folders processing completed!")

def main():
    """
    Main function
    """
    # Input and output paths
    input_dir = r"D:\Gitspace\ipa_full\iPA\data\sxt_images\mrcslice_output"
    output_dir = r"D:\Gitspace\ipa_full\iPA\data\sxt_images\mrcslice_output_denoise"
    
    # Model parameters - using the trained model you provided
    model_name = '2denoise_190920_071249'  # Use your trained model
    basedir = None  # Will use the same directory as the script
    
    # Check if N2V is available
    if not setup_n2v():
        logger.warning("N2V setup failed, will use simple denoising method")
    
    # Check if input directory exists
    if not Path(input_dir).exists():
        logger.error(f"Input directory does not exist: {input_dir}")
        return
    
    logger.info("Starting image denoising processing...")
    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Model name: {model_name}")
    logger.info(f"Model directory: script directory")
    
    # Create denoiser and process
    denoiser = ImageDenoiser(input_dir, output_dir, model_name, basedir)
    denoiser.process_all_folders()

if __name__ == "__main__":
    main()