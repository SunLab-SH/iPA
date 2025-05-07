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
'''
Read the data file, calculate the radial velocity of insulin vesicles projected along the vector from this insulin vesicle to 
the nearest point of the membrane.

Checked the vesicle XYZ and mask XYZ, all coordinates should be in accordance with the XYZ in the cell mask.
'''


#%%
root_dir = 'I:\\Bing\\fluorescence\\Results\\'

# Voxel size: 0.1030x0.1030x0.3006
voxel_size_X = 0.3006
voxel_size_Y = 0.1030
voxel_size_Z = 0.1030

# max coords
#max_coords  = [85*0.125, 480*0.0805940, 480*0.0805940]
# define the function to find the the nearest index in list A for the point b.
def closest_node(node, nodes):
    closest_index = distance.cdist([node], nodes).argmin()
    return nodes[closest_index]

timeslice = "t8"

# seq = "HighGlu"
# dataset = pd.read_csv(f'C2-20210107_G4_N5_522_2.8mM_3D_3.9s_1min_002_D3D_ALX_position.csv')

conditions = ['2.8-30', '16.7-5', '16.7-30']

#%%
shift = np.array(pd.read_csv(root_dir + 'max.csv'))

shift_X = {}
shift_Y = {}
shift_Z = {}
for i in range(len(shift)):
    shift_X[shift[i,0]]= shift[i,1]
    shift_Y[shift[i,0]]= shift[i,2]
    shift_Z[shift[i,0]]= shift[i,3]

# 3D coordinates of ISG, um

for i in conditions:
    tag = root_dir + i + '\\'
    file_dir_list = glob.glob(tag + '*')
    for file_dir_dir in file_dir_list:
        ## 1. import data
        file_name = glob.glob(file_dir_dir + '\\' + '*_Position.csv')
        # print(file_name[0])
        position_file = pd.read_csv(file_name[0], header=2) 
        ret = re.match(r'.*[P,t](.*)_Statistics\\(.*).csv', str(file_name[0]))
        seq = ret.group(2)[:-9]
        print(seq)

        k1 = np.loadtxt(f'{root_dir}01_index\\{seq}_PM_NE_mask_pm_index_{timeslice}.xvg')
        pm_coord_voxel = k1
        k2 = np.loadtxt(f'{root_dir}01_index\\{seq}_PM_NE_mask_ne_index_{timeslice}.xvg')
        ne_center_coord_voxel = k2[0]

        ## 2. extract dt from the info file from Fiji
        with open(f'I:\\Bing\\fluorescence\\3+1D\\Info for {seq}.txt', "r") as f:
            for idx, line in enumerate(f.readlines()):
                line = line.strip('\n')  # 
                if 'Frame interval: ' in line: 
                    line_dt = re.match(r'Frame interval: (.*) sec', line)
                    dt = float(line_dt.group(1))
                    print(dt) # s

        #threshold = 0 # threshold in degrees to determine if a voxel is moving 'towards' the plasma membrane in an irregular cell shape
        threshold = 15/180*math.pi # threshold in degrees to determine if a voxel is moving 'towards' the plasma membrane in an irregular cell shape
        # for the subsequent vzsicle rdf calculation
        coord_voxel_all = []

        # for the subsequent vesicle velocity rdf calculation
        coord_voxel_vel = []
        velocity = []
        velocity_3d_all = []
        direction = []

        ID = list(set(position_file[position_file.columns[7]].tolist()))
        # idall = dataset[dataset.columns[7]].tolist()

        for id in ID: # calculate by each ISG 
            singleisgdata = position_file.loc[(position_file[position_file.columns[7]] == id)]
            singleisgdata = singleisgdata.sort_values([position_file.columns[6]],ascending=True)
            #print(len(singleisgdata))
            single_coord_voxel_all = []
            single_coord_all = []
            
            # Save the coordinates of the vesicles
            X = (singleisgdata[position_file.columns[2]] - shift_Z[seq]).tolist()
            Y = (singleisgdata[position_file.columns[1]] - shift_Y[seq]).tolist()
            Z = (singleisgdata[position_file.columns[0]] - shift_X[seq]).tolist()
            X_voxel = [x/voxel_size_X for x in X] 
            Y_voxel = [x/voxel_size_Y for x in Y] 
            Z_voxel = [x/voxel_size_Z for x in Z] 
            for t in range(len(singleisgdata)):
                coord_voxel_all.append(np.array([X_voxel[t],Y_voxel[t],Z_voxel[t]]))
                single_coord_voxel_all.append(np.array([X_voxel[t],Y_voxel[t],Z_voxel[t]]))
                single_coord_all.append(np.array([X[t],Y[t],Z[t]]))

            # Calculate the velocity
            if len(singleisgdata) > 1:
                timepoints = singleisgdata[position_file.columns[6]].tolist()
                for i in range(len(timepoints) - 1):
                    if timepoints[i+1] - timepoints[i] == 1:
                        velocity_3d = (single_coord_all[i+1]-single_coord_all[i])/dt # um/s
                        velocity_3d_length = np.linalg.norm(velocity_3d) # norm of the vector (length)
                        velocity_3d_normalized = velocity_3d/velocity_3d_length
                        #vector_to_pm = closest_node(single_coord_voxel_all[i],pm_coord_voxel) - single_coord_voxel_all[i]
                        vector_to_pm = (single_coord_voxel_all[i] - ne_center_coord_voxel)
                        vector_to_pm_length = np.linalg.norm(vector_to_pm)
                        vector_to_pm_normlized = vector_to_pm/vector_to_pm_length
                        dot_product = np.dot(velocity_3d_normalized, vector_to_pm_normlized)
                        angle = np.arccos(dot_product)
                        if angle <= threshold:
                            velocity_proj = velocity_3d_length
                        else:
                            velocity_proj = velocity_3d_length*np.cos(angle - threshold)

                        coord_voxel_vel.append(single_coord_voxel_all[i])
                        velocity.append(velocity_proj)
                        velocity_3d_all.append(velocity_3d)

        coord_radial_velocity = np.column_stack([coord_voxel_vel, velocity, velocity, velocity]) # ISG coordinates, radial velocity
        coord_3d_velocity = np.column_stack([coord_voxel_vel, velocity_3d_all]) # ISG coordinates, 3d velocity
        
        with open(f'{root_dir}04_versicle\\{seq}_PM_NE_mask_{timeslice}_coord_voxel_velocity.xvg', 'w') as outfile:
            np.savetxt(outfile, coord_radial_velocity, fmt='%.4f')

        with open(f'{root_dir}04_versicle\\{seq}_PM_NE_mask_{timeslice}_coord_voxel_all.xvg', 'w') as outfile:
            np.savetxt(outfile, coord_voxel_all, fmt='%.4f')

        with open(f'{root_dir}04_versicle\\{seq}_PM_NE_mask_{timeslice}_coord_voxel_3d_velocity.xvg', 'w') as outfile:
            np.savetxt(outfile, coord_3d_velocity, fmt='%.4f')


# %%
#######################
#     plot
#######################
shift = np.array(pd.read_csv(root_dir + 'max.csv'))

shift_X = {}
shift_Y = {}
shift_Z = {}
for i in range(len(shift)):
    shift_X[shift[i,0]]= shift[i,1]
    shift_Y[shift[i,0]]= shift[i,2]
    shift_Z[shift[i,0]]= shift[i,3]

for i in conditions:
    tag = root_dir + i + '\\'
    file_dir_list = glob.glob(tag + '*')
    for file_dir_dir in file_dir_list:
        ## 1. import data
        file_name = glob.glob(file_dir_dir + '\\' + '*_Position.csv')
        print(file_name[0])
        position_file = pd.read_csv(file_name[0], header=2) 
        ret = re.match(r'.*[P,t](.*)_Statistics\\(.*).csv', str(file_name[0]))
        seq = ret.group(2)[:-9]
        print(seq)

        k1 = np.loadtxt(f'{root_dir}01_index\\{seq}_PM_NE_mask_pm_index_{timeslice}.xvg')
        pm_coord_voxel = k1
        k2 = np.loadtxt(f'{root_dir}01_index\\{seq}_PM_NE_mask_ne_index_{timeslice}.xvg')
        ne_center_coord_voxel = k2[0]

        ## 2. extract dt from the info file from Fiji
        with open(f'I:\\Bing\\fluorescence\\3+1D\\Info for {seq}.txt', "r") as f:
            for idx, line in enumerate(f.readlines()):
                line = line.strip('\n')  # 
                if 'Frame interval: ' in line: 
                    line_dt = re.match(r'Frame interval: (.*) sec', line)
                    dt = float(line_dt.group(1))
                    print(dt) # s

        #threshold = 0 # threshold in degrees to determine if a voxel is moving 'towards' the plasma membrane in an irregular cell shape
        threshold = 15/180*math.pi # threshold in degrees to determine if a voxel is moving 'towards' the plasma membrane in an irregular cell shape
        # for the subsequent vzsicle rdf calculation
        coord_voxel_all = []
        single_coord_all = []

        # for the subsequent vesicle velocity rdf calculation
        coord_voxel_vel = []
        velocity = []
        velocity_3d_all = []
        direction = []

        ID = list(set(position_file[position_file.columns[7]].tolist()))
        # idall = dataset[dataset.columns[7]].tolist()

        for id in ID: # calculate by each ISG 
            singleisgdata = position_file.loc[(position_file[position_file.columns[7]] == id)]
            singleisgdata = singleisgdata.sort_values([position_file.columns[6]],ascending=True)
            #print(len(singleisgdata))
            single_coord_voxel_all = []
            
            
            # Save the coordinates of the vesicles
            X = (singleisgdata[position_file.columns[2]] - shift_Z[seq]).tolist()
            Y = (singleisgdata[position_file.columns[1]] - shift_Y[seq]).tolist()
            Z = (singleisgdata[position_file.columns[0]] - shift_X[seq]).tolist()
            X_voxel = [x/voxel_size_X for x in X] 
            Y_voxel = [x/voxel_size_Y for x in Y] 
            Z_voxel = [x/voxel_size_Z for x in Z] 
            
            for t in range(len(singleisgdata)):
                ISG = np.array([X_voxel[t],Y_voxel[t],Z_voxel[t]])
                # for t in range(len(singleisgdata)):
                coord_voxel_all.append(ISG)
                single_coord_all.append(np.array([X[t],Y[t],Z[t]]))

        # print('voxel coord:',np.array(coord_voxel_all, dtype=int)[:,0].max(),
        #       np.array(coord_voxel_all, dtype=int)[:,1].max(), 
        #       np.array(coord_voxel_all, dtype=int)[:,2].max())
        
        # mask_t = np.load(f'{root_dir}masks\\{seq}_PM_NE_mask.npy')
        # sep_tiff01 = np.ones((mask_t.shape[1], mask_t.shape[3], mask_t.shape[4],3))
        # print(sep_tiff01.shape)

        # print('coord: ',np.array(single_coord_all)[:,0].max(),
        #       np.array(single_coord_all)[:,1].max(),
        #       np.array(single_coord_all)[:,2].max())
        
        
        ne = np.loadtxt(f'{root_dir}01_index\\{seq}_PM_NE_mask_ne_index_{timeslice}.xvg')
        xne = [item[0] for item in ne]
        yne = [item[1] for item in ne]
        zne = [item[2] for item in ne]
        pm = np.loadtxt(f'{root_dir}01_index\\{seq}_PM_NE_mask_pm_index_{timeslice}.xvg')
        xpm = [item[0] for item in pm]
        ypm = [item[1] for item in pm]
        zpm = [item[2] for item in pm]

        mask_t = np.load(f'{root_dir}masks\\{seq}_PM_NE_mask.npy')
        sep_tiff01 = np.ones((mask_t.shape[1], mask_t.shape[3], mask_t.shape[4],3))

        for i in range(0, len(ne)):
            sep_tiff01[int(xne[i]),int(yne[i]),int(zne[i]),:] = [0,0,1]

        for i in range(0, len(pm)):
            sep_tiff01[int(xpm[i]),int(ypm[i]),int(zpm[i]),:] = [0,1,0]

        for i in range(0, len(coord_voxel_all)):
            sep_tiff01[int(coord_voxel_all[i][0]),int(coord_voxel_all[i][1]),int(coord_voxel_all[i][2]),:] = [0.8,0,0]

        fig, (ax1) = plt.subplots(1,1, figsize=(7,7))
        #ax1.imshow(sep_tiff01[48, :, :])
        ax1.imshow(sep_tiff01[15, :, :])
        #ax1.imshow(sep_tiff02[232, :, :])
        #sys.exit()

        # Plot a basic wireframe.
        #ax.scatter(xs, ys, zs, c='k', marker='.')
        #ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)
        plt.title(seq)
        plt.show()
        fig.savefig(root_dir + "04_versicle\\plot\\plot_"+ seq + "_isg_ne_pm_z15_t8_" + "16slice.png", transparent=True, dpi=300)


# %%
