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
import time
import glob
from skimage.measure import regionprops
from skimage import io, morphology

import torchvision.transforms as T
import torch


'''
This tif file only has three compartments:

exterior - 0
cytoplasm - 2
nucleus - 3

'''
file_dir = 'I:\\Bing\\fluorescence\\3D\\mask\\'

#%%
#########################
# generate file list
#########################
# files = glob.glob(f'{file_dir}*_MT.labels.mrc')
# for file in files:
# 	print(file[44:-14])

file_list=[]
with open('../file_list.txt','r') as f:
	for line in f:
		file_list.append(line.strip('\n'))

#%%

########################
# for all SIM loop
########################
start_time = time.time() # Measure time elapsed in Python


output_dir = 'I:\\Bing\\fluorescence\\3D\\00_index\\'

####################### save for he junjie ##########################
tif_dir = 'I:\\Fluorescence\\SIM\\20230425_MT-ISG+EX4\\Singlecolor\\'
isg_dir = 'I:\\Bing\\fluorescence\\3D\\03_versicle\\'
# img_dir = 'I:\\Bing\\fluorescence\\3D\\processed_img\\'
corre_x, corre_y, corre_z = 1600, 1600, 64

## read mrc

def resize(img):
		img = torch.tensor(img)
		transform = T.Resize((800, 800))
		img = transform(img)
		img = img.numpy()
		return img

for file_name in file_list:
	print(file_name)
	
	PM_mask = mrcfile.read(f'{file_dir}{file_name}_MT.labels.mrc') # Z, X, Y
	NE_mask = mrcfile.read(f'{file_dir}{file_name}_N.labels.mrc')

	####################### save for he junjie ##########################
	# img_MT = tifffile.imread(tif_dir + file_name + '_MT.tif')#.astype(np.uint8)  # np.array(56, 2560, 2560)
	img_ISG = tifffile.imread(tif_dir + file_name + '_ISG.tif')
	# img_N = tifffile.imread(tif_dir + file_name + '_N.tif')

	props = regionprops(PM_mask[12])
	center = props[0].centroid
	correct_image = np.full((corre_z, 2, corre_x, corre_y), 0)

	####################### save for he junjie ##########################
	# x_low_bnd, x_high_bnd = np.max((0, int(center[0]-corre_x/2))), np.min((corre_x, int(center[0]-corre_x/2)))
	# y_low_bnd, y_high_bnd = np.max((0, int(center[1]-corre_y/2))), np.min((corre_y, int(center[1]-corre_y/2)))
	if 0 > int(center[0]-corre_x/2):
		x_low_bnd = 0
		corr_x_low_bnd = - int(center[0]-corre_x/2) + 1
	else:
		x_low_bnd = int(center[0]-corre_x/2)
		corr_x_low_bnd = 0
	if 0 > int(center[1]-corre_y/2):
		y_low_bnd = 0
		corr_y_low_bnd = - int(center[1]-corre_y/2) + 1
	else:
		y_low_bnd = int(center[1]-corre_y/2)
		corr_y_low_bnd = 0

	if 2560 < int(center[0]+corre_x/2):
		x_high_bnd = 2560
		corr_x_high_bnd = 2560 + corre_x - int(center[0]+corre_x/2)
	else:
		x_high_bnd = int(center[0]+corre_x/2)
		corr_x_high_bnd = corre_x
	if 2560 < int(center[1]+corre_y/2):
		y_high_bnd = 2560
		corr_y_high_bnd = 2560 + corre_y - int(center[1]+corre_y/2)
	else:
		y_high_bnd = int(center[1]+corre_y/2)
		corr_y_high_bnd = corre_y


	correct_image[:PM_mask.shape[0], 0, corr_x_low_bnd:corr_x_high_bnd, corr_y_low_bnd:corr_y_high_bnd] \
		= PM_mask[:, x_low_bnd:x_high_bnd, y_low_bnd:y_high_bnd]
	correct_image[:PM_mask.shape[0], 1, corr_x_low_bnd:corr_x_high_bnd, corr_y_low_bnd:corr_y_high_bnd] \
		= NE_mask[:, x_low_bnd:x_high_bnd, y_low_bnd:y_high_bnd]

	############################
	# Save for he junjie
	############################
	# saved_img = np.full((4, corre_z, corre_x, corre_y), 0)
	# saved_img[0,:PM_mask.shape[0], corr_x_low_bnd:corr_x_high_bnd, corr_y_low_bnd:corr_y_high_bnd] \
	# 	= img_MT[:, x_low_bnd:x_high_bnd, y_low_bnd:y_high_bnd]
	# saved_img[1,:PM_mask.shape[0],corr_x_low_bnd:corr_x_high_bnd, corr_y_low_bnd:corr_y_high_bnd] \
	# 	= img_N[:, x_low_bnd:x_high_bnd, y_low_bnd:y_high_bnd]
	# saved_img[2,:PM_mask.shape[0], corr_x_low_bnd:corr_x_high_bnd, corr_y_low_bnd:corr_y_high_bnd] \
	# 	= img_ISG[:, x_low_bnd:x_high_bnd, y_low_bnd:y_high_bnd]
	# cyto_mask = PM_mask - NE_mask
	# saved_img[3,:PM_mask.shape[0], corr_x_low_bnd:corr_x_high_bnd, corr_y_low_bnd:corr_y_high_bnd] \
	# 	= cyto_mask[:, x_low_bnd:x_high_bnd, y_low_bnd:y_high_bnd]

	# np.save(f'{img_dir}{file_name}_MT_N_isg_mask_center.npy', saved_img)
	correct_isg = np.full((64, 1600, 1600), 0)
	correct_isg[:PM_mask.shape[0], corr_x_low_bnd:corr_x_high_bnd, corr_y_low_bnd:corr_y_high_bnd] \
		= img_ISG[:, x_low_bnd:x_high_bnd, y_low_bnd:y_high_bnd]
	isg_saved = []
	for i in range(64):
		isg_bin2 = resize(correct_isg[i].reshape(-1, 1600, 1600))
		isg_saved.append(isg_bin2)
	isg_saved = np.array(isg_saved, dtype=np.int16).reshape(64,800,800)
	# with open(f'{isg_dir}{file_name}_volumn_bin2_ISG{i}.xvg', 'w') as outfile:
	# 	np.savetxt(outfile, isg_saved, fmt='%.4f')

	pm = correct_image[:, 0, :, :]
	ne = correct_image[:, 1, :, :]
	



	##########################
	# downsampling bin2
	##########################
	
	
	pm_total = []
	ne_total = []
	for i in range(64):
		pm_bin2 = resize(pm[i].reshape(-1, 1600, 1600))
		pm_total.append(pm_bin2)
		ne_bin2 = resize(ne[i].reshape(-1, 1600, 1600))
		ne_total.append(ne_bin2)

	pm_total = np.array(pm_total, dtype=np.int16).reshape(64,800,800)
	ne_total = np.array(ne_total, dtype=np.int16).reshape(64,800,800)

	with mrcfile.new(f'{isg_dir}{file_name}_volumn_bin2_ISG.mrc', overwrite=True) as mrc:
		mrc.set_data(isg_saved*pm_total)

	ne_edt = ndimage.distance_transform_edt(ne_total) # calculate the distance from non-zero to the nearest zero.
	ne_ndx = np.array(np.where(ne_edt[:,:,:] == 1.0))
	center = [np.mean(ne_ndx[0]), np.mean(ne_ndx[1]), np.mean(ne_ndx[2])] # the coordinate of the NE center
	ne_edge=ne_ndx.T
	ne_out = open(f'{output_dir}{file_name}_bin2_ne_index.xvg', 'w') # output
	np.savetxt(ne_out, [center], fmt='%.2f')
	np.savetxt(ne_out, ne_edge, fmt='%d')
	pm_edt = ndimage.distance_transform_edt(pm_total) # calculate the distance from non-zero to the nearest zero.
	pm_ndx = np.array(np.where(pm_edt[:,:,:] == 1.0))
	pm_edge=pm_ndx.T 
	pm_out=open(f'{output_dir}{file_name}_bin2_pm_index.xvg', 'w') # output
	np.savetxt(pm_out, pm_edge, fmt='%d')
	

print("--- %s seconds ---" % (time.time() - start_time))
















#%% BACKUP
########################
# for all SIM loop
########################
start_time = time.time() # Measure time elapsed in Python

file_dir = 'I:\\Bing\\fluorescence\\3D\\mask\\'
output_dir = 'I:\\Bing\\fluorescence\\3D\\00_index\\'
tif_dir = 'I:\\Fluorescence\\SIM\\20230425_MT-ISG+EX4\\Singlecolor\\'
img_dir = 'I:\\Bing\\fluorescence\\3D\\processed_img\\'
corre_x, corre_y, corre_z = 1600, 1600, 64

## read mrc
file_list=[]
with open('../file_list.txt','r') as f:
	for line in f:
		file_list.append(line.strip('\n'))

for file_name in file_list:
	print(file_name)
	
	PM_mask = mrcfile.read(f'{file_dir}{file_name}_MT.labels.mrc') # Z, X, Y
	NE_mask = mrcfile.read(f'{file_dir}{file_name}_N.labels.mrc')
	img_MT = tifffile.imread(tif_dir + file_name + '_MT.tif')#.astype(np.uint8)  # np.array(56, 2560, 2560)
	img_ISG = tifffile.imread(tif_dir + file_name + '_ISG.tif')
	img_N = tifffile.imread(tif_dir + file_name + '_N.tif')

	props = regionprops(PM_mask[12])
	center = props[0].centroid
	correct_image = np.full((corre_z, 2, corre_x, corre_y), 0)
	correct_image[:PM_mask.shape[0], 0, :, :] \
		= PM_mask[:, int(center[0]-corre_x/2):int(center[0]+corre_x/2), int(center[1]-corre_y/2):int(center[1]+corre_y/2)]
	correct_image[:PM_mask.shape[0], 1, :, :] \
		= NE_mask[:, int(center[0]-corre_x/2):int(center[0]+corre_x/2), int(center[1]-corre_y/2):int(center[1]+corre_y/2)]

	saved_img = np.full((4, corre_z, corre_x, corre_y), 0)
	saved_img[0,:PM_mask.shape[0], :, :] \
		= img_MT[:, int(center[0]-corre_x/2):int(center[0]+corre_x/2), int(center[1]-corre_y/2):int(center[1]+corre_y/2)]
	saved_img[1,:PM_mask.shape[0], :, :] \
		= img_N[:, int(center[0]-corre_x/2):int(center[0]+corre_x/2), int(center[1]-corre_y/2):int(center[1]+corre_y/2)]
	saved_img[2,:PM_mask.shape[0], :, :] \
		= img_ISG[:, int(center[0]-corre_x/2):int(center[0]+corre_x/2), int(center[1]-corre_y/2):int(center[1]+corre_y/2)]
	cyto_mask = PM_mask - NE_mask
	saved_img[3,:PM_mask.shape[0], :, :] \
		= cyto_mask[:, int(center[0]-corre_x/2):int(center[0]+corre_x/2), int(center[1]-corre_y/2):int(center[1]+corre_y/2)]

	np.save(f'{img_dir}{file_name}_MT_N_isg_mask_center.npy', saved_img)


	pm = correct_image[:, 0, :, :]
	ne = correct_image[:, 1, :, :]
	
	ne_edt = ndimage.distance_transform_edt(ne) # calculate the distance from non-zero to the nearest zero.
	ne_ndx = np.array(np.where(ne_edt[:,:,:] == 1.0))
	center = [np.mean(ne_ndx[0]), np.mean(ne_ndx[1]), np.mean(ne_ndx[2])] # the coordinate of the NE center
	ne_edge=ne_ndx.T
	ne_out = open(f'{output_dir}{file_name}_ne_index.xvg', 'w') # output
	np.savetxt(ne_out, [center], fmt='%.2f')
	np.savetxt(ne_out, ne_edge, fmt='%d')
	pm_edt = ndimage.distance_transform_edt(pm) # calculate the distance from non-zero to the nearest zero.
	pm_ndx = np.array(np.where(pm_edt[:,:,:] == 1.0))
	pm_edge=pm_ndx.T 
	pm_out=open(f'{output_dir}{file_name}_pm_index.xvg', 'w') # output
	np.savetxt(pm_out, pm_edge, fmt='%d')
	

print("--- %s seconds ---" % (time.time() - start_time))
#%% 
########################
# for single SIM example
########################

# root_dir = 'F:\\Bing\\Features\\Results\\'
# 0+10-1-2_SIM_MT.labels.mrc
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
np.savetxt(pm_out, pm_edge)
	
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
