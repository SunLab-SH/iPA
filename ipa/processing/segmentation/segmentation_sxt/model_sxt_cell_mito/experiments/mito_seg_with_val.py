#%%
# -*- coding: utf-8 -*-
# @Time    : 2020-06-02 22:10
# @Author  : Xiangyi Zhang
# @File    : eval_mem_nu.py
# @Email   : zhangxy9@shanghaitech.edu.cn

import os, sys
from unittest import TestLoader
# filepath = os.path.dirname(os.path.abspath(__file__))
# os.chdir(filepath)

# Simplified logging configuration
import logging
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger(__name__)

filepath = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(filepath)

sys.path.append(filepath)
sys.path.append(parent_dir) 



from experiments.parser import args_mito as args

from dataloaders.dataset import *


from networks.Unet import Unet
from torch import nn
import logging

import cv2
LOG = logging.getLogger('main')

import tifffile, mrcfile
from utils.src import *

os.environ['CUDA_VISIBLE_DEVICES'] = args.gpu
from pylab import plt
from multiprocessing.dummy import Pool
import time


# print(args)  
# print(type(args.num_classes))
# image_size = (288, 480) # (288, 480)

#%%
'''
1: mem
2: mito
3: nu
4: gr
test_idx = ['784_5', '766_8', '842_17']
val_idx = ['783_5', '766_5', '842_12']
train_idx = ['766_2', '766_7', '766_10', '766_11', '769_5', '769_7', '783_6', '783_12', '784_4', '784_6', '784_7', '785_7', '822_4', '822_6', '822_7', '842_13', '931_11', '931_14']
'''
def normalize(image):
    '''
    for 2d image to normalize 
    '''
    image = (image - np.min(image) )/ np.max(image) * 255
    return image.astype(np.uint8)


class load_data():
    def __init__(self, dir, dataid) -> None:
        self.dir = dir
        self.dataid = dataid

    def get_filelist(self, dir, Filelist):
        newDir = dir
        if os.path.isfile(dir):
            Filelist.append(dir)
            # Filelist.append(os.path.basename(dir))
        elif os.path.isdir(dir):
            for s in os.listdir(dir):
                # If you need to ignore certain folders, use the following code
                #if s == "xxx":
                    #continue
                newDir=os.path.join(dir,s)
                self.get_filelist(newDir, Filelist)

        return Filelist


    def get_rawimgname(self):
        namelist = self.get_filelist(self.dir, [])
        # print(namelist)
        for name in namelist:
            if f'{self.dataid}' in name and 'pre_rec.tif' in name:
                return name
            elif f'{self.dataid}' in name and 'pre_rec.mrc' in name:
                return name
        else:
            print(' raw img not exist')
            os._exit(0)


    def get_labelname(self):
        # Remove label-related functionality, only used for prediction
        return None


    def get_labelfile(self):
        # Remove label-related functionality, only used for prediction
        return None

class cellmapping_test():
    def __init__(self, type='imageall_x', dataid=None, pool=None, image_data=None):
        self.type = type
        # Remove dependency on pre-sliced image directory
        # self.image_root_dir = osp.join('data/image_xyz/', self.type)
        self.checkpoint_dir =  'data/unet_pth/mito_best.pth' # 'results/{}/best.pth'.format(args.exp)
        self.device = torch.device("cuda")
        # Remove dependency on pre-sliced label directory
        # self.gt_path = osp.join('data/mask_xyz')
        # self.coarse_label_dir = 'data/coarse_label'
        self.pool = pool
        self.image_data = image_data  # Pre-loaded image data

        self.data_idxs = dataid if dataid else []
        self.test_idx = self.data_idxs

    def save_tiff(self, pred, imagefile):
        savepath = f'{args.root_dir}//results//masks//{self.type}_{args.mode}_predtiff'
        os.makedirs(savepath, exist_ok=True)
        
        tifffile.imwrite(f'{savepath}/{imagefile}.tiff', np.array(pred))
        print(f"save {savepath}/{imagefile}.tiff")

    # Remove all evaluation-related methods
    # def evaluate_iou(self, x, y, num_class):
    # def resize(self, pred, image_ori_path=None):
    # def get_path(self, image_root_dir):
    # def turnGtBigLabel(self, imagefile):
    # def evaluate(self, pred_list, gt):

    def turnPredMitoLabel(self, pred):
        pred[pred != 1] = 0
        return pred

    ### Inference ###
    def test(self):
        self.model = Unet(n_class=args.num_classes + 1, is_dropout=True)
        self.model = nn.DataParallel(self.model).cuda()
        self.model.load_state_dict(torch.load(self.checkpoint_dir))
        self.model = self.model.module.to(torch.device('cpu'))  # run on cpu
        print('Resume from {}'.format(self.checkpoint_dir))

        print(self.data_idxs)
        for dataid in self.data_idxs:
            
            # Use pre-loaded image data
            if self.image_data is not None and dataid in self.image_data:
                image_x = self.image_data[dataid]
                print(f'Using preloaded 3D image for {dataid}')
            else:
                # If no pre-loaded data, fallback to original loading method
                data_a = load_data(os.path.join(args.data_root_dir), dataid)
                namee = data_a.get_rawimgname()
                print(f'Loading 3D image: {namee}')

                nameelst = namee.split( '.' )
                if nameelst[-1] == 'mrc':
                    image_x = mrcfile.open(namee, permissive=True).data
                elif nameelst[-1] == 'tif':
                    image_x = tifffile.imread(namee)

            # Rotate 3D image based on view type
            if self.type == 'imageall_x':
                image = image_x
            elif self.type == 'imageall_y':
                image = np.rot90(image_x, k=3, axes=(2, 0))
            elif self.type == 'imageall_z':
                image = np.rot90(image_x, k=1, axes=(1, 0))

            # Get 3D image shape and create output mask
            imgshape = image.shape
            print(f'Image shape: {imgshape}')
            outputmask = np.zeros_like(image)
            
            self.model.eval()
            
            logging.info(f'\nPredicting image {dataid} with view {self.type}...')

            # Process 3D image slice by slice
            for ii in range(imgshape[0]):
                imageslice = image[ii]
                shape_ = imageslice.shape
                
                inputslice = normalize(imageslice)
                image_size = (288, 480)
                inputslice = cv2.resize(inputslice, image_size, cv2.INTER_NEAREST)
                
                transform = Compose([tf.ToTensor(),
                        tf.Normalize(mean=[0.456], std=[0.224]),
                        ])
                sample = {'image': inputslice}
                imagedata = transform(sample)
                imageslice_tensor = imagedata['image']
                predslice = self.step(imageslice_tensor, dataid)
                predslice = predslice.cpu().numpy()
                predslice = cv2.resize(predslice.astype(np.uint8), (shape_[1], shape_[0]), cv2.INTER_NEAREST)
                
                outputmask[ii] = predslice

                if ii % 50 == 0:
                    print(f'Processing slice {ii}/{imgshape[0]}, max: {np.max(predslice)}, min: {np.min(predslice)}')

            print(f'Final output max: {np.max(outputmask)}')
            self.save_tiff(outputmask, dataid)
            print(f'Prediction completed for {dataid} with view {self.type}')

    def step(self, image, imagefile=None):
        self.imagefile= imagefile
        image = torch.unsqueeze(image, 0)
        images = image
        images = images.to(torch.device('cpu')) 
        outputs = self.model(images)
        pred = torch.argmax(torch.softmax(outputs, dim=1), dim=1)
        if args.mode == 'mito':
            pred = self.turnPredMitoLabel(pred)
        pred = torch.squeeze(pred, 0)
        return pred


class Post_process(cellmapping_test):
    def __init__(self, dataid=None):
        super(Post_process, self).__init__(type='imageall_x', dataid=dataid)
        # Remove dependency on pre-sliced data directory
        # self.path = f'data/mask_xyz/'
        if dataid:
            self.data_idxs = dataid if isinstance(dataid, list) else [dataid]
        else:
            self.data_idxs = []

    def getIsolateLabel(self, x, label):
        x_copy = x.copy()
        x_copy[x_copy != label] = 255
        x_copy[x_copy == label] = 1
        x_copy[x_copy == 255] = 0
        return x_copy

    def fix012(self, xyzLabel, xyz, x, y):
        if len(x[xyz==6] == 2) !=  len(x[xyz==6]) :
            print('012 appear !')

    def save_tiff_coarse(self, pred, imagefile):
        savepath = f'{args.root_dir}//results//masks'
        os.makedirs(savepath, exist_ok=True)
        tifffile.imwrite(f'{savepath}/{imagefile}_pred.tiff', np.array(pred))
        print(f"save {savepath}/{imagefile}_pred.tiff")

    def post_fusexyz(self):
        for imagefile in self.data_idxs:
            self.imagefile = imagefile
            print(f"Processing fusion for {self.imagefile}")
            mask_base_dir = f'{args.root_dir}//results//masks'
            x = tifffile.imread(os.path.join(mask_base_dir, f'imageall_x_{args.mode}_predtiff', f'{self.imagefile}.tiff'))
            y = tifffile.imread(os.path.join(mask_base_dir, f'imageall_y_{args.mode}_predtiff', f'{self.imagefile}.tiff'))
            z = tifffile.imread(os.path.join(mask_base_dir, f'imageall_z_{args.mode}_predtiff', f'{self.imagefile}.tiff'))
            
            turned_y = np.rot90(y, k=1, axes=(2, 0))
            turned_z = np.rot90(z, k=3, axes=(1, 0))
            print('x y z turnedy turnedz', x.shape, y.shape,z.shape, turned_y.shape, turned_z.shape)
            assert x.shape == turned_y.shape and x.shape == turned_z.shape

            xyzLabel = np.zeros_like(x)
            for p in [0, 1, 2]:
                xyz = x + turned_y + turned_z
                xi = self.getIsolateLabel(x, p)
                yi = self.getIsolateLabel(turned_y, p)
                zi = self.getIsolateLabel(turned_z, p)
                xyzi = xi + yi + zi
                xyzLabel[xyzi >= 2] = p

            self.fix012(xyzLabel, xyz, x, turned_y)
            self.save_tiff_coarse(xyzLabel, self.imagefile)
            print(f'Fusion completed for {imagefile}')


def run_mito_segmentation(data_root_dir=None, model_path=None, save_dir=None, 
                    gpu='0', mode='mito', 
                    num_classes=4, 
                    pool_processes=8, view_types=['imageall_x', 'imageall_y', 'imageall_z'],
                    dataid=None, image_data=None):
    """
    Main function for running cellular structure segmentation with multi-view prediction and fusion.
    
    This function performs deep learning-based segmentation of cellular structures (mitochondria, nucleus, membrane, etc.)
    using a U-Net model. It processes 3D images from multiple viewing angles (X, Y, Z axes) and fuses
    the predictions to generate accurate segmentation masks.
    
    Workflow:
        1. Load pre-trained U-Net model
        2. Process input 3D images slice-by-slice from three orthogonal views
        3. Apply normalization and resizing for model input
        4. Generate predictions for each view
        5. Perform post-processing fusion to combine multi-view results
        6. Save final segmentation masks
    
    Args:
        data_root_dir (str, optional): Root directory containing input data files.
            Used when image_data is None. Should contain .mrc or .tif files.
        model_path (str): Path to the pre-trained PyTorch model (.pth file).
            Model should be compatible with the specified num_classes.
        save_dir (str): Output directory for saving segmentation results.
            Will create subdirectories for intermediate and final results.
        gpu (str, optional): GPU device identifier. Defaults to '0'.
            Set to CPU-compatible value if CUDA unavailable.
        mode (str, optional): Segmentation target mode. Defaults to 'mito'.
            Options: 'mito', 'nucleus', 'membrane', 'mem_nu' (combined).
        num_classes (int, optional): Number of segmentation classes excluding background.
            Defaults to 2. Background class is automatically added.
        pool_processes (int, optional): Number of parallel processes for multiprocessing.
            Defaults to 8. Adjust based on available CPU cores.
        view_types (list, optional): List of viewing angles to process.
            Defaults to ['imageall_x', 'imageall_y', 'imageall_z'].
            Each view corresponds to slicing along different axes.
        dataid (str or list): Identifier(s) for specific datasets to process.
            Required parameter. Can be single string or list of strings.
        image_data (dict, optional): Pre-loaded image data dictionary.
            Format: {dataid: numpy_array}. If provided, bypasses file loading.
    
    Returns:
        str: Path to the directory containing final segmentation results.
            Results include both individual view predictions and fused outputs.
    
    Raises:
        ValueError: If dataid parameter is missing or improperly formatted.
        FileNotFoundError: If specified model_path or input data files cannot be located.
        RuntimeError: If GPU operations fail or insufficient memory for processing.
        TypeError: If image_data contains incompatible data types or shapes.
    
        
    Notes:
        - Input images are automatically normalized to [0, 255] range
        - Model predictions are resized back to original image dimensions
        - Multi-view fusion uses majority voting for robust segmentation
        - Intermediate results are saved for debugging and analysis
        - GPU memory usage scales with image size and number of classes
    """
    import os
    import sys
    from multiprocessing.dummy import Pool
    
    # Simplified logging - only show important messages
    logging.getLogger().setLevel(logging.WARNING)
    


    model_path = f'{parent_dir}/models/mito_best.pth'

    if data_root_dir:
        args.data_root_dir = data_root_dir
    args.gpu = gpu
    args.mode = mode
    args.num_classes = num_classes
    args.root_dir = save_dir
    
    if dataid is not None:
        if isinstance(dataid, str):
            processed_dataid = [dataid]
        elif isinstance(dataid, list):
            processed_dataid = dataid
        else:
            raise ValueError("dataid must be string or list of strings")
    else:
        raise ValueError("dataid parameter is required for prediction")
    
    os.makedirs(save_dir, exist_ok=True)
    results_dir = os.path.join(save_dir, 'results', 'masks')
    os.makedirs(results_dir, exist_ok=True)
    
    class ConfigurableCellMapping(cellmapping_test):
        def __init__(self, type, model_path, data_root, save_root, data_idxs, image_data=None, pool=None):
            super().__init__(type, data_idxs, pool, image_data)
            self.checkpoint_dir = model_path
            self.data_root_dir = data_root
            self.save_root = save_root
            self.data_idxs = data_idxs
            self.test_idx = self.data_idxs

    class ConfigurablePostProcess(Post_process):
        def __init__(self, data_root, save_root, model_path, data_idxs):
            self.type = 'imageall_x'
            self.data_root_dir = data_root
            self.save_root = save_root
            self.checkpoint_dir = model_path
            self.data_idxs = data_idxs
            self.test_idx = self.data_idxs
        
        def save_tiff_coarse(self, pred, imagefile):
            savepath = os.path.join(self.save_root, 'results', 'masks')
            os.makedirs(savepath, exist_ok=True)
            tifffile.imwrite(f'{savepath}/{imagefile}_pred.tiff', np.array(pred))
            print(f"save {savepath}/{imagefile}_pred.tiff")
    
    pool = Pool(pool_processes)

    
    print(f"Starting prediction, mode: {mode}, number of classes: {num_classes}")
    if image_data:
        print(f"Using pre-loaded image data, data IDs: {list(image_data.keys())}")
    else:
        print(f"Data directory: {data_root_dir}")
    print(f"Model path: {model_path}")
    print(f"Save directory: {save_dir}")
    print(f"Processing data IDs: {processed_dataid}")
    
    # Process different views
    for view_type in view_types:
        print(f'############# {view_type} ###############')
        test_instance = ConfigurableCellMapping(
            type=view_type, 
            model_path=model_path,
            data_root=data_root_dir,
            save_root=save_dir,
            data_idxs=processed_dataid,
            image_data=image_data,
            pool=pool
        )
        test_instance.test()

    # Post-processing fusion
    print("Starting post-processing fusion...")
    post_process = ConfigurablePostProcess(
        data_root=data_root_dir,
        save_root=save_dir,
        model_path=model_path,
        data_idxs=processed_dataid
    )
    post_process.post_fusexyz()
    
    final_results_path = os.path.join(save_dir, 'results', 'masks')
    print(f"Prediction completed! Results saved in: {final_results_path}")
    
    return final_results_path


# if __name__ == "__main__":
#     # Example usage
#     import os
#     mainpath = f'D:/Gitspace/ipa_full/iPA'

#     # Method 1: Pre-load data
#     dataid_list = ['784_5']
#     image_data = {}
    
#     # Manually load image data
#     for dataid in dataid_list:
#         data_loader = load_data(f'{mainpath}/data', dataid)
#         image_path = data_loader.get_rawimgname()
#         print(f'Loading {image_path}')
        
#         if image_path.endswith('.mrc'):
#             image_data[dataid] = mrcfile.open(image_path, permissive=True).data
#         elif image_path.endswith('.tif'):
#             image_data[dataid] = tifffile.imread(image_path)
    
#     model_path = f'{mainpath}/data/unet_pth/mito_best.pth'
#     save_dir = f'{mainpath}/results/mito_seg_with_val'

#     # Use pre-loaded data for segmentation
#     run_cell_segmentation(
#         model_path=model_path,
#         save_dir=save_dir,
#         mode='mito',
#         pool_processes=6,
#         dataid=dataid_list,
#         image_data=image_data  # Pass pre-loaded image data
#     )



# if __name__ == "__main__":
#     # Keep original execution logic unchanged
#     pool = Pool(8)
#     logging.basicConfig(level=logging.INFO)
#     # for type in ['imageall_y', 'imageall_z']:
#     for type in ['imageall_x', 'imageall_y', 'imageall_z']:
#         print('#############', type, '###############')
#         Test = cellmapping_test(type=type)
#         Test.test()
#     post_process = Post_process()
#     post_process.post_fusexyz()


# %%

#%% test

# %%

#%% test
#         Test = cellmapping_test(type=type)
#         Test.test()
#     post_process = Post_process()
#     post_process.post_fusexyz()


# %%

#%% test

# %%

#%% test

