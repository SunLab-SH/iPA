#!/usr/bin/env python
# coding: utf-8
#%%
import os
import sys
import numpy as np
import tifffile
import time, os, math, json, copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tifffile
import mrcfile
import scipy.ndimage as ndimage
from scipy.ndimage import gaussian_filter, gaussian_laplace
from skimage.feature import peak_local_max
from scipy import spatial
# from skimage.morphology import watershed
from scipy import ndimage as ndi
import random
import skimage
import time
import glob

'''
This tif file only has three compartments:

exterior - 0
cytoplasm - 2
nucleus - 3

'''

#%%
# root_dir = 'F:\\Bing\\Features\\Results\\'
file = 'I:\\Bing\\fluorescence\\Results\\SIM\\0min_3_5_SIM_PM_NE_mask.npy'
start_time = time.time() # Measure time elapsed in Python

## read npy
# file_list = glob.glob(f'{root_dir}masks\\*.npy')
# for file in file_list:
mask_t = np.load(file) # T, z, C(PM/NE), x, y
# out_name = file[31:-4]


pm = mask_t[:, 0, :, :]
ne = mask_t[:, 1, :, :]

ne_edt = ndimage.distance_transform_edt(ne) # calculate the distance from non-zero to the nearest zero.
ne_ndx = np.array(np.where(ne_edt[:,:,:] == 1.0))
center = [np.mean(ne_ndx[0]), np.mean(ne_ndx[1]), np.mean(ne_ndx[2])] # the coordinate of the NE center

ne_edge=ne_ndx.T
ne_out = open('I:\\Bing\\fluorescence\\Results\\SIM\\0min_3_5_SIM' + '_ne_index'+ '.xvg', 'w') # output
np.savetxt(ne_out, [center], fmt='%.2f')
np.savetxt(ne_out, ne_edge, fmt='%d')

pm_edt = ndimage.distance_transform_edt(pm) # calculate the distance from non-zero to the nearest zero.
pm_ndx = np.array(np.where(pm_edt[:,:,:] == 1.0))

pm_edge=pm_ndx.T 
pm_out=open('I:\\Bing\\fluorescence\\Results\\SIM\\0min_3_5_SIM' + '_pm_index'+ '.xvg', 'w') # output
np.savetxt(pm_out, pm_edge, fmt='%d')
    
print("--- %s seconds ---" % (time.time() - start_time))

# tiff_index = tiff[:,:,:]
# shape3d = tiff_index.shape

# print(shape3d)
# # extend y axis both from top and bottom for r, g, b channels
# tiff_extend = np.zeros((shape3d[0]+2,shape3d[1],shape3d[2]))
# tiff_extend[1:shape3d[0]+1,:,:] = tiff[:,:,:]

# print(tiff_extend.shape)
# # check number in this tiff, three channels for rgb color, _r, _g, _b are arbitary names
# tiff_index = tiff_extend[:,:,:]
# print(tiff[1,200,100:300])
# print(tiff[21,240,100:300])
# print(tiff[29,240,100:300])

# %%
