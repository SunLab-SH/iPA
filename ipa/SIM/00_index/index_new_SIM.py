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
import time
import glob
from skimage.measure import regionprops
from skimage import io, morphology


'''
This tif file only has three compartments:

exterior - 0
cytoplasm - 2
nucleus - 3

'''
#%%

########################
# for all SIM loop
########################
start_time = time.time() # Measure time elapsed in Python

file_dir = 'I:\\Bing\\fluorescence\\3D\\mask\\'
output_dir = 'I:\\Bing\\fluorescence\\3D\\00_index\\'
corre_x, corre_y, corre_z = 1280, 1280, 64

## read mrc
file_list=[]
with open('../file_list.txt','r') as f:
	for line in f:
		file_list.append(line.strip('\n'))

for file_name in file_list[11:15]:
	print(file_name)
	ne_out = open(f'{output_dir}{file_name}_ne_index.xvg', 'w') # output
	pm_out=open(f'{output_dir}{file_name}_pm_index.xvg', 'w') # output

	PM_mask = mrcfile.read(f'{file_dir}{file_name}_MT.labels.mrc') # Z, X, Y
	NE_mask = mrcfile.read(f'{file_dir}{file_name}_N.labels.mrc')
	props = regionprops(PM_mask[12])
	center = props[0].centroid
	correct_image = np.full((corre_z, 2, corre_x, corre_y), 0)
	correct_image[:PM_mask.shape[0], 0, :, :] \
		= PM_mask[:, int(center[0]-corre_x/2):int(center[0]+corre_x/2), int(center[1]-corre_y/2):int(center[1]+corre_y/2)]
	correct_image[:PM_mask.shape[0], 1, :, :] \
		= NE_mask[:, int(center[0]-corre_x/2):int(center[0]+corre_x/2), int(center[1]-corre_y/2):int(center[1]+corre_y/2)]

	ne_edt_3d = ndimage.distance_transform_edt(ne)
	ne_ndx_3d = np.array(np.where(ne_edt_3d[:,:,:] == 1.0))
	center_3d = [np.mean(ne_ndx_3d[0]), np.mean(ne_ndx_3d[1]), np.mean(ne_ndx_3d[2])]
	np.savetxt(ne_out, [center_3d], fmt='%.2f')

	pm = correct_image[:, 0, :, :]
	ne = correct_image[:, 1, :, :]

	for z in range(1,64,1):
		if (np.unique(ne[z-1]) == 0).all() and (np.unique(ne[z]) != 0).all():
			ne_ndx = np.array(np.where(ne[z][:,:] == 1)).T
			ne_edge_3d = np.full((ne_ndx.shape[0], 3), z)
			ne_edge_3d[:, 1:3] = ne_ndx
			np.savetxt(ne_out, ne_edge_3d, fmt='%d')
			continue
		elif (np.unique(ne[z-1]) != 0).all() and (np.unique(ne[z]) == 0).all():
			ne_ndx = np.array(np.where(ne[z-1][:,:] == 1)).T
			ne_edge_3d = np.full((ne_ndx.shape[0], 3), z)
			ne_edge_3d[:, 1:3] = ne_ndx
			np.savetxt(ne_out, ne_edge_3d, fmt='%d')
			break

		ne_edt = ndimage.distance_transform_edt(ne[z])
		ne_ndx = np.array(np.where(ne_edt[:,:] == 1.0))
		ne_edge=ne_ndx.T
		if ne_edge!=[]: 
			ne_edge_3d = np.full((ne_edge.shape[0], 3), z)
			ne_edge_3d[:, 1:3] = ne_edge
			np.savetxt(ne_out, ne_edge_3d, fmt='%d')


	for z in range(63):
		if (np.unique(ne[z-1]) == 0).all() and (np.unique(ne[z]) != 0).all():
			pm_ndx = np.array(np.where(pm[z][:,:] == 1)).T
			pm_edge_3d = np.full((pm_ndx.shape[0], 3), z)
			pm_edge_3d[:, 1:3] = pm_ndx
			np.savetxt(ne_out, pm_edge_3d, fmt='%d')
			continue
		elif (np.unique(ne[z-1]) != 0).all() and (np.unique(ne[z]) == 0).all():
			pm_ndx = np.array(np.where(ne[z-1][:,:] == 1)).T
			pm_edge_3d = np.full((pm_ndx.shape[0], 3), z)
			pm_edge_3d[:, 1:3] = pm_ndx
			np.savetxt(ne_out, pm_edge_3d, fmt='%d')
			break

		pm_edt = ndimage.distance_transform_edt(pm[z]) # calculate the distance from non-zero to the nearest zero.
		pm_ndx = np.array(np.where(pm_edt[:,:] == 1.0))
		pm_edge=pm_ndx.T 

		if pm_edge != []:
			pm_edge_3d = np.full((pm_edge.shape[0], 3), z)
			pm_edge_3d[:, 1:3] = pm_edge
			np.savetxt(pm_out, pm_edge_3d, fmt='%d')

	print("--- %s seconds ---" % (time.time() - start_time))
# %%
