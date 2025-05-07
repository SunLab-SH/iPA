#%%
import numpy as np
import math
import matplotlib.pyplot as plt
import tifffile
from skimage import data,filters,segmentation,measure,morphology,color
from skimage.measure import regionprops
import mrcfile
#%%
t = np.linspace(0, 2*math.pi)

x = 2*np.sin(t)
y = 3*np.cos(t)

plt.plot(x, y, c = 'red')
plt.axvline(x = 0)
plt.axhline(y = 0)
plt.text(0.2, 1, r'$\frac{x^2}{2^2}+\frac{y^2}{3^2}=1$',fontsize = 15)

# %%
PM_mask = np.zeros((128*128, 2))
matrix_size = (1280, 1280)  
center = (640, 640)  
radius_x = 200  
radius_y = 400  


matrix = np.zeros(matrix_size, dtype=np.int32)


for i in range(matrix_size[0]):
    for j in range(matrix_size[1]):
        if ((i - center[0]) / radius_x) ** 2 + ((j - center[1]) / radius_y) ** 2 <= 1:
            matrix[i, j] = 1


plt.imshow(matrix, cmap='gray')
plt.axis('off')
plt.show()
# %%
file_dir = 'I:\\Fluorescence\\SIM\\20230425_MT-ISG+EX4\\Singlecolor\\'
mask_dir = 'I:\\Bing\\fluorescence\\3D\\mask\\'
output_dir = 'I:\\Bing\\fluorescence\\3D\\03_versicle\\'

file_list=[]
with open(f'../file_list.txt','r') as f:
	for line in f:
		file_list.append(line.strip('\n'))

for file_name in file_list[12:28]:
    print(file_name)

    isg_out = open(f'{output_dir}{file_name}_ne_index.xvg', 'w') # output
    PM_mask = mrcfile.read(f'{mask_dir}{file_name}_MT.labels.mrc')
    img_ISG = tifffile.imread(f'{file_dir}{file_name}_ISG.tif') #.astype(np.uint8)  # np.array(56, 2560, 2560)

    for z in range(10, img_ISG.shape[0]):
        thresh = filters.threshold_otsu(img_ISG[z]) 
        ISG_voxel = morphology.closing(img_ISG[z] > thresh, morphology.square(3))

        downsampled_matrix = ISG_voxel.reshape(2560 // 2, 2, 2560 // 2, 2).mean(axis=(1, 3))
        binary_matrix = np.where(downsampled_matrix >= 0.5, 1, 0)
        binary_matrix = binary_matrix & PM_mask[z, 640:1920, 640:1920]

        plt.imshow(matrix & binary_matrix, cmap='gray')
        plt.axis('off')
        plt.show()
        
        # ISG_ndx = np.array(np.where(ISG_voxel[:,:] == 1)).T

        # np.savetxt(isg_out, ISG_ndx, fmt='%d')
    break
        
# %%
