'''
Author: Yu Bing
Data: 2023-05-04

1. build masks of nucleus and cell
2. align all cell to the central site
3. cut edge to reduce storage

'''
# Actin-ISG: 'I:\\Fluorescence\\SIM\\20220909_Actin-ISG\\Singlecolor'
# MT-ISG: 'I:\\Fluorescence\\SIM\\20230320MT-ISG\\Singlecolor\\'
# MT-ISG(Ex4): 'I:\\Fluorescence\\SIM\\20230425_MT-ISG+EX4\\Singlecolor\\'

#%% -1
from imaris_ims_file_reader.ims import ims
import matplotlib.pyplot as plt

file_name = 'I:\\Fluorescence\\SIM\Actin-ISG\\0min_3_5_SIM_2.tif'
a = ims(file_name)

print(a.ResolutionLevelLock)
print(a.ResolutionLevels)
print(a.TimePoints)
print(a.Channels)
print(a.shape)
print(a.chunks)
print(a.dtype)
print(a.ndim)




#%%
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from skimage import data,filters,segmentation,measure,morphology,color
import tifffile
import cv2
import scipy.ndimage as ndimage
from scipy import ndimage
import scipy.signal as signal
from scipy.optimize import linear_sum_assignment
from scipy.ndimage import gaussian_filter
from skimage.measure import label, regionprops
from skimage.filters import threshold_otsu
from skimage.segmentation import clear_border, watershed
from skimage.measure import label, regionprops
from skimage.morphology import binary_erosion, binary_dilation, binary_opening

from PIL import Image
from PIL import ImageEnhance
import mrcfile

#%%
file_dir = 'I:\\Fluorescence\\SIM\\20230425_MT-ISG+EX4\\Singlecolor\\3\\'
file_list = glob.glob(file_dir + '*_MT.tif')

output_dir = 'I:\\Bing\\fluorescence\\3D\\mask\\'
corre_x, corre_y, corre_z = 1280, 1280, 64
# %%
count = 0

for i in range(len(file_list)):
    seq = file_list[i][54:-7]
    print(f'file name: {seq}')
    img_MT = tifffile.imread(file_list[i])#.astype(np.uint8)  # np.array(56, 2560, 2560)
    # img_ISG = tifffile.imread(file_dir + seq + '_ISG.tif')
    img_N = tifffile.imread(file_dir + seq + '_N.tif')
    mask_all = np.zeros((img_MT.shape[0],2,2560,2560))
    z_label = img_MT.shape[0]

    for z in range(img_MT.shape[0]):
        thresh_MT = filters.threshold_otsu(img_MT[z]) 
        mask_PM = morphology.closing(img_MT[z] > thresh_MT, morphology.square(3))
        
        mask = np.copy(mask_PM)
        # mask[np.where(mask_ISG[:,:] == True)] = True

        # Remove small objects (voxels) with less than 100 pixels
        mask_PM_cleaned = morphology.remove_small_objects(mask, min_size=100)
        # plt.imshow(mask_PM_cleaned,plt.cm.gray)
        # plt.show()
        # mask_PM_filled = ndimage.binary_fill_holes(mask_PM).astype(int)
        # mask_PM_01[np.where(mask_PM[:,:] == True)] = 1
        # props = regionprops(mask_PM_01)

        # image_filledhole = ndimage.binary_fill_holes(binaried).astype(int)
        # image_bk = np.copy(image_filledhole)
        # plt.imshow(image_filledhole, plt.cm.gray)
        # plt.title('image_filledhole')

        # Label the objects in the stack
        label_stack = label(mask_PM_cleaned)
        # Create an array to store the properties of the objects
        areas = []
        # Loop through each object and calculate its properties
        for obj in regionprops(label_stack):
            # Ignore small objects
            if obj.area < 100:
                continue
            # Add the properties to the array
            areas.append(obj.area)
        
        if areas == []:
            print('blank')
            
        elif np.array(areas).mean() >= 400 or z > (img_MT.shape[0]/2-1):
            print(np.array(areas).mean())
        # # Sort the objects by size
        # properties.sort(key=lambda x: x.area, reverse=True)
        # # Select the largest object as the cell
        # cell = properties[0]
        # # Create a mask of the cell
        # cell_mask = np.zeros_like(mask_PM_cleaned, dtype=bool)
        # cell_mask[label_stack == cell.label] = True
        # plt.imshow(cell_mask,plt.cm.gray)
        # plt.show()
        # break

            thresh_N = filters.threshold_otsu(img_N[z]) 
            mask_N = morphology.closing(img_N[z] > thresh_N, morphology.square(3))

            # Erode the mask to remove the noise signal from other layers
            eroded_mask = binary_erosion(mask_PM_cleaned, footprint=np.ones((5,5)))
            # Dilate the eroded mask to recover the lost information
            dilated_mask = binary_dilation(eroded_mask, footprint=np.ones((5,5)))
            # Perform an opening operation to smooth the mask
            opened_mask = binary_opening(dilated_mask, footprint=np.ones((2,2))) # selen=
            # Apply morphological closing to smooth the mask
            kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (100, 100))
            smoothed_mask = cv2.morphologyEx(opened_mask.astype('uint8'), cv2.MORPH_CLOSE, kernel)
            # Fill holes
            image_filledhole = ndimage.binary_fill_holes(smoothed_mask).astype(int)
            # plt.imshow(image_filledhole,plt.cm.gray)
            # plt.title(f'pmZ{z}')
            # plt.show()

            # Add NE
            props = regionprops(image_filledhole)
            bbox = props[0].bbox # (min_row, min_col, max_row, max_col)
            if z == 20:
                center = props[0].centroid
            img_bbox = np.zeros_like(img_MT[z], dtype=bool)
            img_bbox[bbox[0]:bbox[2], bbox[1]:bbox[3]] = True

            mask_N_one = img_bbox & mask_N
            opened_mask = binary_opening(mask_N_one, selem=np.ones((2,2)))
            kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (40, 40))
            smoothed_mask = cv2.morphologyEx(opened_mask.astype('uint8'), cv2.MORPH_CLOSE, kernel)
            NE_mask = ndimage.binary_fill_holes(smoothed_mask).astype(int)

            label_stack = label(NE_mask)
            properties = []
            for obj in regionprops(label_stack):
                # Ignore small objects
                if obj.area < 1000:
                    continue
                # Add the properties to the array
                properties.append(obj)

            if properties == []:
                if count == 0: z_label = z
                count += 1
                NE_mask_final = np.zeros_like(mask_PM_cleaned, dtype=bool)
            elif z > z_label:
                NE_mask_final = np.zeros_like(mask_PM_cleaned, dtype=bool)
            else:
                # Sort the objects by size
                properties.sort(key=lambda x: x.area, reverse=True)
                # Select the largest object as the cell
                cell = properties[0]
                # Create a mask of the cell
                cell_mask = np.zeros_like(mask_PM_cleaned, dtype=bool)
                cell_mask[label_stack == cell.label] = True

                kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (20, 20))
                NE_mask_final = cv2.morphologyEx(cell_mask.astype('uint8'), cv2.MORPH_OPEN, kernel)
                # plt.imshow(NE_mask_final,plt.cm.gray)
                # plt.title(f'Z{z}')
                # plt.show()

            kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (40, 40))
            PM_mask_0 = cv2.morphologyEx((image_filledhole|NE_mask_final).astype('uint8'), cv2.MORPH_CLOSE, kernel)
            PM_mask = cv2.morphologyEx((PM_mask_0).astype('uint8'), cv2.MORPH_OPEN, kernel)

            # Find the contours of the masked shape
            contours = measure.find_contours(PM_mask, 0.5)
            print(len(contours))
            mask_all[z,0,:,:] = PM_mask
            mask_all[z,1,:,:] = NE_mask_final

            # smoothed_contour = filters.gaussian(contours[0], sigma=2, preserve_range=True)

            # fig, ax = plt.subplots(figsize=(5,5))
            # ax.plot(smoothed_contour[:, 0], smoothed_contour[:, 1], linewidth=2)
            # # Reverse the y-axis
            # ax.invert_yaxis()
            # # Set the range of all axes to start from zero
            # ax.set_xlim([0, 2560])
            # ax.set_ylim([0, 2560])
            # plt.title(f'pm_Z{z}')
            # plt.show()    

        else:
            print(np.array(areas).mean())
            image_filledhole = np.zeros_like(img_MT[z])

    # move cell to center and padding to the same size
    correct_image = np.full((corre_z, 2, corre_x, corre_y), 0) # Z, C, X, Y
    shift_x, shift_y = int(center[0]- 2560/2), int(center[1] - 2560/2)
    correct_image[:img_MT.shape[0], :,:,:] \
        = mask_all[:, :, int(center[0]-corre_x/2):int(center[0]+corre_x/2), int(center[1]-corre_y/2):int(center[1]+corre_y/2)]
    
    for z in range(0, img_MT.shape[0], 5):
        plt.imshow(correct_image[z][0],plt.cm.gray)
        plt.title(f'{seq}_Z{z}')
        plt.savefig(f'{output_dir}plots\\Ex4_{seq}_PM_z{z}_mask.png', transparent = True, dpi = 300)

        plt.imshow(correct_image[z][1],plt.cm.gray)
        plt.title(f'{seq}_Z{z}')
        plt.savefig(f'{output_dir}plots\\Ex4_{seq}_NE_z{z}_mask.png', transparent = True, dpi = 300)
            
    np.save(f'{output_dir}Ex4_{seq}_PM_NE_mask.npy', correct_image)

    # break
# %%
