# coding = 'utf-8'
#%%
import numpy as np
import csv
import pandas as pd
from datetime import date,datetime
import sys, os, glob
import matplotlib.pyplot as plt
from scipy.spatial import distance
import math
import re

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
'''
Read the data file, calculate the radial velocity of insulin vesicles projected along the vector from this insulin vesicle to 
the nearest point of the membrane.

Checked the vesicle XYZ and mask XYZ, all coordinates should be in accordance with the XYZ in the cell mask.
'''


#%%
file_dir = 'I:\\Fluorescence\\SIM\\20230425_MT-ISG+EX4\\Singlecolor\\'    # Modified
file_dir2 = 'I:\\Bing\\fluorescence\\3D\\00_index\\'    # Modified
output_dir = 'I:\\Bing\\fluorescence\\3D\\06_stat\\'    # Modified

file_list=[]
with open(f'../file_list_MT_ex4.txt','r') as f:     # Modified
	for line in f:
		file_list.append(line.strip('\n'))

for file_name in file_list:
    print(file_name)

    isg_out = open(f'{output_dir}{file_name}_bin2_mt_index.xvg', 'w') # output
    img_ISG = mrcfile.open(f'{output_dir}{file_name}_volumn_bin2_MT.mrc', mode='r').data
    
    # thresh = filters.threshold_otsu(img_ISG[12])
    for z in range(img_ISG.shape[0]):
        thresh = filters.threshold_otsu(img_ISG[z]) 
        ISG_voxel = morphology.closing(img_ISG[z] > thresh, morphology.square(3))
        ISG_ndx = np.array(np.where(ISG_voxel[:,:] == 1)).T
        ISG_z = np.full((len(ISG_ndx),1),z)
        ISG = np.hstack((ISG_z,ISG_ndx))
        np.savetxt(isg_out, ISG, fmt='%d') # , fmt='%d'


#%%
########################
#   Plot
########################
file_list=[]
with open(f'../file_list_MT.txt','r') as f:    # Modified
	for line in f:
		file_list.append(line.strip('\n'))

for file_name in file_list:
    print(file_name)
    
    ne = np.loadtxt(f'{file_dir2}{file_name}_bin2_ne_index.xvg')
    xne = [item[0] for item in ne]
    yne = [item[1] for item in ne]
    zne = [item[2] for item in ne]
    pm = np.loadtxt(f'{file_dir2}{file_name}_bin2_pm_index.xvg')
    xpm = [item[0] for item in pm]
    ypm = [item[1] for item in pm]
    zpm = [item[2] for item in pm]
    isg = np.loadtxt(f'{output_dir}{file_name}_bin2_isg_index.xvg')
    xisg = [item[0] for item in isg]
    yisg = [item[1] for item in isg]
    zisg = [item[2] for item in isg]

    sep_tiff01 = np.ones((64,800,800,3))

    for i in range(0, len(ne)):
        sep_tiff01[int(xne[i]),int(yne[i]),int(zne[i]),:] = [0,0,1]

    for i in range(0, len(pm)):
        sep_tiff01[int(xpm[i]),int(ypm[i]),int(zpm[i]),:] = [0,1,0]

    for i in range(0, len(isg)):
        sep_tiff01[int(isg[i][0]),int(isg[i][1]),int(isg[i][2]),:] = [0.8,0.1,0]

    fig, (ax1) = plt.subplots(1,1, figsize=(7,7))
    #ax1.imshow(sep_tiff01[48, :, :])
    ax1.imshow(sep_tiff01[15, :, :])
    #ax1.imshow(sep_tiff02[232, :, :])
    #sys.exit()

    # Plot a basic wireframe.
    #ax.scatter(xs, ys, zs, c='k', marker='.')
    #ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)
    plt.title(file_name)
    plt.show()
    fig.savefig(f'{output_dir}plots\\{file_name}_isg_ne_pm_z15_8slice.png', transparent=True, dpi=60)





# %%
