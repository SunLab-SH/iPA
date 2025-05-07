# coding = utf-8


'''
To calculate filament distance to vesicle membrane.
To calculate filaments points around individual surrounding spaces.

By angdi

'''
#%%

import numpy as np
import scipy.ndimage
import scipy
import math 
from scipy.spatial.distance import cdist
import copy 
import pandas as pd

#%%
def shift_bias(filament_coords, shift):
    # shift should be minused 
    filament_coords_modified ={}    
    filamentnames = list(filament_coords.keys())
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords[f'{filamentnames[idx]}']
        filament_coords_modified[f'{filamentnames[idx]}'] =[]

        for i, coord in enumerate(curfilamentcoordslst):
            coord_new = [coord[0] - shift[0], coord[1] - shift[1], coord[2] - shift[2]]
            # coord_new = [coord_new[0]/voxelsize[0], coord_new[1]/voxelsize[1], coord_new[2]/voxelsize[2]]
            filament_coords_modified[f'{filamentnames[idx]}'].append(coord_new)

    return filament_coords_modified



def coords_to_mem_distance_generator(filament_coords, mask, voxelsize):
    '''
    input: 
        filament coord: dict
        mask: array 
    '''

    filament_coords_modified ={}    
    filamentnames = list(filament_coords.keys())
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords[f'{filamentnames[idx]}']
        filament_coords_modified[f'{filamentnames[idx]}'] =[]
        for i, coord_new in enumerate(curfilamentcoordslst):
            coord_new = [coord_new[0]/voxelsize[0], coord_new[1]/voxelsize[1], coord_new[2]/voxelsize[2]]
            filament_coords_modified[f'{filamentnames[idx]}'].append(coord_new)


    coords_lst = []
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords_modified[f'{filamentnames[idx]}']
        coords_lst.extend(curfilamentcoordslst)

    coords = np.array(coords_lst).reshape(-1,3)
    # print('line56', coords)
    voxel_coords = coords

    coords_x = [ voxel_coords[i][0] for i in range(len(voxel_coords))]
    coords_y = [ voxel_coords[i][1] for i in range(len(voxel_coords))]
    coords_z = [ voxel_coords[i][2] for i in range(len(voxel_coords))]

    # print(sorted(coords_x))
    
    # print(np.max(coords_x), np.max(coords_y), np.max(coords_z))
    # print(np.min(coords_x), np.min(coords_y), np.min(coords_z))
    
    # print('coords range', [[np.round(np.min(coords_x)), np.round(np.max(coords_x))],[np.round(np.min(coords_y)), np.round(np.max(coords_y))],[np.round(np.min(coords_z)), np.round(np.max(coords_z))]] )
    # print('mask size', mask.shape)
    # print('coords in total', len(voxel_coords))
    # if True:
    #     print('line78')
    #     print('x',np.max(coords_x),np.min(coords_x), mask.shape[0])
    #     print('y',np.max(coords_y),np.min(coords_y), mask.shape[1])
    #     print('z',np.max(coords_z),np.min(coords_z), mask.shape[2])


    if not (np.max(coords_x) <= mask.shape[0] and np.min(coords_x)>= 0) or not (np.max(coords_y) <= mask.shape[1] and np.min(coords_y)>= 0) or not (np.max(coords_z) <= mask.shape[2] and np.min(coords_z)>= 0):
        print('x',np.max(coords_x),np.min(coords_x), mask.shape[0])
        print('y',np.max(coords_y),np.min(coords_y), mask.shape[1])
        print('z',np.max(coords_z),np.min(coords_z), mask.shape[2])

    if len(coords_x) > 0:
        assert np.max(coords_x) <= mask.shape[0] and np.min(coords_x)>= 0
        assert np.max(coords_y) <= mask.shape[1] and np.min(coords_y)>= 0
        assert np.max(coords_z) <= mask.shape[2] and np.min(coords_z)>= 0
    else:
        print('no filament')

    img01 = mask
    reverse_img01 = img01 * (-1) + 1
    img01_rev_edt = scipy.ndimage.distance_transform_edt(reverse_img01)

    avesize = np.average(voxelsize)
    dist = [ img01_rev_edt[int(voxel_coords[i][0])][int(voxel_coords[i][1])][int(voxel_coords[i][2])]  for i in range(len(voxel_coords)) ]
    dist = [dist[i] * avesize/10 for i in range(len(dist))] # nm
    # print('distance num', len(dist))

    return dist





def actin_to_vesicle_distance_generator(filament_coords, mask, voxelsize):
    '''
    input: 
        filament coord: dict
        mask: array 
    '''

    filament_coords_modified ={}    
    filamentnames = list(filament_coords.keys())
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords[f'{filamentnames[idx]}']
        filament_coords_modified[f'{filamentnames[idx]}'] =[]
        for i, coord_new in enumerate(curfilamentcoordslst):
            coord_new = [coord_new[0]/voxelsize[0], coord_new[1]/voxelsize[1], coord_new[2]/voxelsize[2]]
            filament_coords_modified[f'{filamentnames[idx]}'].append(coord_new)


    coords_lst = []
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords_modified[f'{filamentnames[idx]}']
        coords_lst.extend(curfilamentcoordslst)

    coords = np.array(coords_lst).reshape(-1,3)
    # print('line56', coords)
    voxel_coords = coords

    coords_x = [ voxel_coords[i][0] for i in range(len(voxel_coords))]
    coords_y = [ voxel_coords[i][1] for i in range(len(voxel_coords))]
    coords_z = [ voxel_coords[i][2] for i in range(len(voxel_coords))]

    # print(sorted(coords_x))
    
    # print(np.max(coords_x), np.max(coords_y), np.max(coords_z))
    # print(np.min(coords_x), np.min(coords_y), np.min(coords_z))
    
    # print('coords range', [[np.round(np.min(coords_x)), np.round(np.max(coords_x))],[np.round(np.min(coords_y)), np.round(np.max(coords_y))],[np.round(np.min(coords_z)), np.round(np.max(coords_z))]] )
    # print('mask size', mask.shape)
    # print('coords in total', len(voxel_coords))


    assert np.max(coords_x) <= mask.shape[0] and np.min(coords_x)>= 0
    assert np.max(coords_y) <= mask.shape[1] and np.min(coords_y)>= 0
    assert np.max(coords_z) <= mask.shape[2] and np.min(coords_z)>= 0


    img01 = mask
    reverse_img01 = img01 * (-1) + 1
    img01_rev_edt = scipy.ndimage.distance_transform_edt(reverse_img01)

    avesize = np.average(voxelsize)
    dist = [ img01_rev_edt[int(voxel_coords[i][0])][int(voxel_coords[i][1])][int(voxel_coords[i][2])]  for i in range(len(voxel_coords)) ]
    dist = [dist[i] * avesize/10 for i in range(len(dist))] # nm
    # print('distance num', len(dist))

    return dist








def actin_to_actin_distance(filament_coord_extend_dict, imgsize, voxel_size_xyz):
    'input: filaments distance dict'
    edge = 500
    x_bound = [0 + edge, imgsize[0]*voxel_size_xyz[0] - edge]
    y_bound = [0 + edge, imgsize[1]*voxel_size_xyz[1] - edge]
    z_bound = [0 + edge, imgsize[2]*voxel_size_xyz[2] - edge]

    filament_coord_extend_lst = []
    filamentnames = list(filament_coord_extend_dict.keys())

    for idx in range(len(filamentnames)):
        # curfilamentcoordslst = filament_coord_extend_dict[f'{filamentnames[idx]}']
        curfilamentcoordslst = [ np.array(coord) for coord in filament_coord_extend_dict[f'{filamentnames[idx]}'] ]
        filament_coord_extend_lst.append(curfilamentcoordslst)


    filament_coord_lst = filament_coord_extend_lst
    filaments_vect_lst = [] 

    for single_filament_coord in filament_coord_lst:
        # coords1, coord2 = single_filament_coord[0], single_filament_coord[-1]
        vect = np.array(single_filament_coord[-1]) - np.array(single_filament_coord[0])
        unit_vect = vect / np.linalg.norm(vect)
        # unit_vect = vect
        filaments_vect_lst.append(unit_vect)

    filament_img_vect = np.array([0,0,0])
    for vect1 in filaments_vect_lst:
        filament_img_vect = filament_img_vect + vect1



    filament_vect_map_lst = []

    for i, curfilaments_coords in enumerate(filament_coord_extend_lst):
        rest_filament_coord_extend_lst = copy.deepcopy(filament_coord_extend_lst)
        _ = rest_filament_coord_extend_lst.pop(i)
        cur_filament_surr_coords = []
        vect_single_filament_map_lst = []
        vect_single_filament_angle_map_lst = []
        check_normvect = []
        curfilaments_coords_filtered= []
        for coords in curfilaments_coords:
            # print(coords)
            # print(x_bound,y_bound,z_bound)
            if (coords[0] >= x_bound[0] and  coords[0] <= x_bound[1]) or \
                (coords[1] >= y_bound[0] and  coords[1] <= y_bound[1]) or \
                (coords[2] >= z_bound[0] and  coords[2] <= z_bound[1]):

                curfilaments_coords_filtered.append(coords)
        # print('filtered num', len(curfilaments_coords), len(curfilaments_coords_filtered), len(curfilaments_coords) - len(curfilaments_coords_filtered))
        curfilaments_coords = curfilaments_coords_filtered
        if len(curfilaments_coords) > 0:

            for j, rest_single_filament_coords in enumerate(rest_filament_coord_extend_lst):
            # print(len(rest_filament_coord_extend_lst), len(filament_coord_extend_lst))  # n-1, n
            # print(rest_filament_coord_extend_lst[:5], '\n', filament_coord_extend_lst[:5])  
                dists = cdist(curfilaments_coords, rest_single_filament_coords, metric='euclidean') # return array, size [a,b]
                # print(dists)
                filtered_coords = np.where(dists <= edge ) #500 = 50nm  #100nm
                surrfilament_idx_check_lst = []
                for ii in range(len(filtered_coords[0])): 
                    if filtered_coords[1][ii] not in surrfilament_idx_check_lst:
                        surrfilament_idx_check_lst.append(filtered_coords[1][ii])
                        cur_filament_surr_coords.append(rest_single_filament_coords[filtered_coords[1][ii]])

            # get all surround coords and set vector
            for k, surround_coord in enumerate(cur_filament_surr_coords):
                distss = cdist([surround_coord], curfilaments_coords, metric='euclidean')
                min_d = np.min(distss)
                min_d_idx = np.where(distss == min_d)[1][0]

                cur_fila_coord = curfilaments_coords[min_d_idx] 
                cur_fila_correspond_coord = surround_coord

                if min_d_idx+1 >= len(curfilaments_coords):
                    cur_direction_vect = curfilaments_coords[min_d_idx] - curfilaments_coords[min_d_idx-1]

                elif min_d_idx - 1 >= 0:
                    cur_direction_vect = curfilaments_coords[min_d_idx] - curfilaments_coords[min_d_idx-1]
                else:
                    cur_direction_vect = curfilaments_coords[min_d_idx+1] - curfilaments_coords[min_d_idx]
                

                if np.linalg.norm(cur_direction_vect) > 0:  ## check if coord is not nan
                    # depend on single filament direction to )filaments overall direction
                    vect_direction = np.dot(cur_direction_vect,filament_img_vect) 
                    temp_vect_in_single_vect = cur_fila_correspond_coord - cur_fila_coord
                        
                    if vect_direction >= 0:
                        normed_mornvect, temp_vect_in_single_vect = rotate_normvect2z(cur_direction_vect, temp_vect_in_single_vect)
                    else:
                        normed_mornvect, temp_vect_in_single_vect = rotate_normvect2nz(cur_direction_vect, temp_vect_in_single_vect)
                    
                    vect_single_filament_map_lst.append(temp_vect_in_single_vect) # rotate among cur filament coord


        filament_vect_map_lst.append(vect_single_filament_map_lst) # vect to surrounding in a flat

    # print(filament_coord_extend_dict)
    # print(filament_vect_map_lst)
    # print('filament num', len(filament_vect_map_lst))


    return filament_vect_map_lst



def actin_to_actin_distance2(filament_coord_extend_dict, imgsize, voxel_size_xyz):
    'input: filaments distance dict'


    edge = 500
    x_bound = [0 + edge, imgsize[0]*voxel_size_xyz[0] - edge]
    y_bound = [0 + edge, imgsize[1]*voxel_size_xyz[1] - edge]
    z_bound = [0 + edge, imgsize[2]*voxel_size_xyz[2] - edge]

    filament_coord_extend_lst = []
    filamentnames = list(filament_coord_extend_dict.keys())

    for idx in range(len(filamentnames)):
        # curfilamentcoordslst = filament_coord_extend_dict[f'{filamentnames[idx]}']
        curfilamentcoordslst = [ np.array(coord) for coord in filament_coord_extend_dict[f'{filamentnames[idx]}'] ]
        filament_coord_extend_lst.append(curfilamentcoordslst)


    filament_coord_lst = filament_coord_extend_lst
    filaments_vect_lst = [] 

    for single_filament_coord in filament_coord_lst:
        # coords1, coord2 = single_filament_coord[0], single_filament_coord[-1]
        vect = np.array(single_filament_coord[-1]) - np.array(single_filament_coord[0])
        unit_vect = vect / np.linalg.norm(vect)
        # unit_vect = vect
        filaments_vect_lst.append(unit_vect)

    filament_img_vect = np.array([0,0,0])
    for vect1 in filaments_vect_lst:
        filament_img_vect = filament_img_vect + vect1



    filament_vect_map_lst = []
    dist = []

    for i, curfilaments_coords in enumerate(filament_coord_extend_lst):
        rest_filament_coord_extend_lst = copy.deepcopy(filament_coord_extend_lst)
        _ = rest_filament_coord_extend_lst.pop(i)
        cur_filament_surr_coords = []
        vect_single_filament_map_lst = []
        vect_single_filament_angle_map_lst = []
        check_normvect = []
        curfilaments_coords_filtered= []
        for coords in curfilaments_coords:
            # print(coords)
            # print(x_bound,y_bound,z_bound)
            if (coords[0] >= x_bound[0] and  coords[0] <= x_bound[1]) or \
                (coords[1] >= y_bound[0] and  coords[1] <= y_bound[1]) or \
                (coords[2] >= z_bound[0] and  coords[2] <= z_bound[1]):

                curfilaments_coords_filtered.append(coords)
        # print('filtered num', len(curfilaments_coords), len(curfilaments_coords_filtered), len(curfilaments_coords) - len(curfilaments_coords_filtered))
        curfilaments_coords = curfilaments_coords_filtered
        if len(curfilaments_coords) > 0:

            for j, rest_single_filament_coords in enumerate(rest_filament_coord_extend_lst):
            # print(len(rest_filament_coord_extend_lst), len(filament_coord_extend_lst))  # n-1, n
            # print(rest_filament_coord_extend_lst[:5], '\n', filament_coord_extend_lst[:5])  
                dists = cdist(curfilaments_coords, rest_single_filament_coords, metric='euclidean') # return array, size [a,b]
                # print(dists)
                filtered_coords = np.where(dists <= edge ) #500 = 50nm  #100nm
                surrfilament_idx_check_lst = []
                for ii in range(len(filtered_coords[0])): 
                    if filtered_coords[1][ii] not in surrfilament_idx_check_lst:
                        surrfilament_idx_check_lst.append(filtered_coords[1][ii])
                        cur_filament_surr_coords.append(rest_single_filament_coords[filtered_coords[1][ii]])

            # get all surround coords and set vector
            for k, surround_coord in enumerate(cur_filament_surr_coords):
                distss = cdist([surround_coord], curfilaments_coords, metric='euclidean')
                min_d = np.min(distss)
                min_d_idx = np.where(distss == min_d)[1][0]

                cur_fila_coord = curfilaments_coords[min_d_idx] 
                cur_fila_correspond_coord = surround_coord

                if min_d_idx+1 >= len(curfilaments_coords):
                    cur_direction_vect = curfilaments_coords[min_d_idx] - curfilaments_coords[min_d_idx-1]

                elif min_d_idx - 1 >= 0:
                    cur_direction_vect = curfilaments_coords[min_d_idx] - curfilaments_coords[min_d_idx-1]
                else:
                    cur_direction_vect = curfilaments_coords[min_d_idx+1] - curfilaments_coords[min_d_idx]
                

                if np.linalg.norm(cur_direction_vect) > 0:  ## check if coord is not nan
                    # depend on single filament direction to )filaments overall direction
                    vect_direction = np.dot(cur_direction_vect,filament_img_vect) 
                    temp_vect_in_single_vect = cur_fila_correspond_coord - cur_fila_coord
                        
                    if vect_direction >= 0:
                        normed_mornvect, temp_vect_in_single_vect = rotate_normvect2z(cur_direction_vect, temp_vect_in_single_vect)
                    else:
                        normed_mornvect, temp_vect_in_single_vect = rotate_normvect2nz(cur_direction_vect, temp_vect_in_single_vect)
                    
                    vect_single_filament_map_lst.append(temp_vect_in_single_vect) # rotate among cur filament coord

                    dist.append(np.linalg.norm(temp_vect_in_single_vect))




    return dist



def actin_to_microtube_distance(filament_coord_extend_dict, MT_coord_extend_dict):
    'input: filaments distance dict'
    
    filament_coord_extend_lst = []
    filamentnames = list(filament_coord_extend_dict.keys())

    for idx in range(len(filamentnames)):
        # curfilamentcoordslst = filament_coord_extend_dict[f'{filamentnames[idx]}']
        curfilamentcoordslst = [ np.array(coord) for coord in filament_coord_extend_dict[f'{filamentnames[idx]}'] ]
        filament_coord_extend_lst.append(curfilamentcoordslst)




    MT_coord_extend_lst = [] # [mt1 coords, mt2 coords, ..., mt3 coords ]
    MTnames = list(MT_coord_extend_dict.keys())

    for idx in range(len(MTnames)):
        # curfilamentcoordslst = filament_coord_extend_dict[f'{filamentnames[idx]}']
        curMTcoordslst = [ np.array(coord) for coord in MT_coord_extend_dict[f'{MTnames[idx]}'] ]
        MT_coord_extend_lst.append(curMTcoordslst) 

    



    MT_coord_lst = MT_coord_extend_lst
    MT_vect_lst = [] 

    for single_MT_coord in MT_coord_lst:
        # coords1, coord2 = single_filament_coord[0], single_filament_coord[-1]
        vect = np.array(single_MT_coord[-1]) - np.array(single_MT_coord[0])
        unit_vect = vect / np.linalg.norm(vect)
        # unit_vect = vect
        MT_vect_lst.append(unit_vect)

    MT_img_vect = np.array([0,0,0])
    for vect1 in MT_vect_lst:
        MT_img_vect = MT_img_vect + vect1





    # calculate dist map
    filament_vect_2_MT_map_lst = []

    for i, curMT_coords in enumerate(MT_coord_extend_lst):
        # rest_filament_coord_extend_lst = copy.deepcopy(filament_coord_extend_lst)
        # _ = rest_filament_coord_extend_lst.pop(i)
        cur_filament_surr_coords = []
        vect_single_filament_map_lst = []
        vect_single_filament_angle_map_lst = []
        check_normvect = []


        for j, single_filament_coords in enumerate(filament_coord_extend_lst):
        # print(len(rest_filament_coord_extend_lst), len(filament_coord_extend_lst))  # n-1, n
        # print(rest_filament_coord_extend_lst[:5], '\n', filament_coord_extend_lst[:5])  
            dists = cdist(curMT_coords, single_filament_coords, metric='euclidean') # return array, size [a,b]
            # print(dists)
            filtered_coords = np.where(dists <= 1000 ) #100nm
            surrfilament_idx_check_lst = []
            for ii in range(len(filtered_coords[0])): 
                if filtered_coords[1][ii] not in surrfilament_idx_check_lst:
                    surrfilament_idx_check_lst.append(filtered_coords[1][ii])
                    cur_filament_surr_coords.append(single_filament_coords[filtered_coords[1][ii]])

        # get all surround coords and set vector
        for k, surround_coord in enumerate(cur_filament_surr_coords):
            distss = cdist([surround_coord], curMT_coords, metric='euclidean')
            min_d = np.min(distss)
            min_d_idx = np.where(distss == min_d)[1][0]

            cur_fila_coord = curMT_coords[min_d_idx] 
            cur_fila_correspond_coord = surround_coord

            if min_d_idx > 1:
                cur_direction_vect = curMT_coords[min_d_idx] - curMT_coords[min_d_idx-1]
            else:
                cur_direction_vect = curMT_coords[min_d_idx+1] - curMT_coords[min_d_idx]
            
            if np.linalg.norm(cur_direction_vect) > 0:  ## check if coord is not nan
                # depend on single filament direction to )filaments overall direction
                vect_direction = np.dot(cur_direction_vect, MT_img_vect) 
                temp_vect_in_single_vect = cur_fila_correspond_coord - cur_fila_coord
                

                if vect_direction >= 0:
                    normed_mornvect, temp_vect_in_single_vect = rotate_normvect2z(cur_direction_vect, temp_vect_in_single_vect)
                else:
                    normed_mornvect, temp_vect_in_single_vect = rotate_normvect2nz(cur_direction_vect, temp_vect_in_single_vect)
                check_normvect.append(normed_mornvect)

                vect_single_filament_map_lst.append(temp_vect_in_single_vect) # rotate among cur filament coord


        filament_vect_2_MT_map_lst.append(vect_single_filament_map_lst) # vect to surrounding in a flat

    # print(filament_coord_extend_dict)
    # print(filament_vect_map_lst)
    # print('filament num', len(filament_vect_map_lst))


    return filament_vect_2_MT_map_lst


def actin_to_microtube_distance2(filament_coord_extend_dict, MT_coord_extend_dict):
    'input: filaments distance dict'
    
    filament_coord_extend_lst = []
    filamentnames = list(filament_coord_extend_dict.keys())

    for idx in range(len(filamentnames)):
        # curfilamentcoordslst = filament_coord_extend_dict[f'{filamentnames[idx]}']
        curfilamentcoordslst = [ np.array(coord) for coord in filament_coord_extend_dict[f'{filamentnames[idx]}'] ]
        filament_coord_extend_lst.append(curfilamentcoordslst)




    MT_coord_extend_lst = [] # [mt1 coords, mt2 coords, ..., mt3 coords ]
    MTnames = list(MT_coord_extend_dict.keys())

    for idx in range(len(MTnames)):
        # curfilamentcoordslst = filament_coord_extend_dict[f'{filamentnames[idx]}']
        curMTcoordslst = [ np.array(coord) for coord in MT_coord_extend_dict[f'{MTnames[idx]}'] ]
        MT_coord_extend_lst.append(curMTcoordslst) 

    



    MT_coord_lst = MT_coord_extend_lst
    MT_vect_lst = [] 

    for single_MT_coord in MT_coord_lst:
        # coords1, coord2 = single_filament_coord[0], single_filament_coord[-1]
        vect = np.array(single_MT_coord[-1]) - np.array(single_MT_coord[0])
        unit_vect = vect / np.linalg.norm(vect)
        # unit_vect = vect
        MT_vect_lst.append(unit_vect)

    MT_img_vect = np.array([0,0,0])
    for vect1 in MT_vect_lst:
        MT_img_vect = MT_img_vect + vect1





    # calculate dist map
    filament_vect_2_MT_map_lst = []
    dist = []

    for i, curMT_coords in enumerate(MT_coord_extend_lst):
        # rest_filament_coord_extend_lst = copy.deepcopy(filament_coord_extend_lst)
        # _ = rest_filament_coord_extend_lst.pop(i)
        cur_filament_surr_coords = []
        vect_single_filament_map_lst = []
        vect_single_filament_angle_map_lst = []
        check_normvect = []


        for j, single_filament_coords in enumerate(filament_coord_extend_lst):
        # print(len(rest_filament_coord_extend_lst), len(filament_coord_extend_lst))  # n-1, n
        # print(rest_filament_coord_extend_lst[:5], '\n', filament_coord_extend_lst[:5])  
            dists = cdist(curMT_coords, single_filament_coords, metric='euclidean') # return array, size [a,b]
            # print(dists)
            filtered_coords = np.where(dists <= 1000 ) #100nm
            surrfilament_idx_check_lst = []
            for ii in range(len(filtered_coords[0])): 
                if filtered_coords[1][ii] not in surrfilament_idx_check_lst:
                    surrfilament_idx_check_lst.append(filtered_coords[1][ii])
                    cur_filament_surr_coords.append(single_filament_coords[filtered_coords[1][ii]])

        # get all surround coords and set vector
        for k, surround_coord in enumerate(cur_filament_surr_coords):
            distss = cdist([surround_coord], curMT_coords, metric='euclidean')
            min_d = np.min(distss)
            min_d_idx = np.where(distss == min_d)[1][0]

            cur_fila_coord = curMT_coords[min_d_idx] 
            cur_fila_correspond_coord = surround_coord

            if min_d_idx > 1:
                cur_direction_vect = curMT_coords[min_d_idx] - curMT_coords[min_d_idx-1]
            else:
                cur_direction_vect = curMT_coords[min_d_idx+1] - curMT_coords[min_d_idx]
            
            if np.linalg.norm(cur_direction_vect) > 0:  ## check if coord is not nan
                # depend on single filament direction to )filaments overall direction
                vect_direction = np.dot(cur_direction_vect, MT_img_vect) 
                temp_vect_in_single_vect = cur_fila_correspond_coord - cur_fila_coord
                

                if vect_direction >= 0:
                    normed_mornvect, temp_vect_in_single_vect = rotate_normvect2z(cur_direction_vect, temp_vect_in_single_vect)
                else:
                    normed_mornvect, temp_vect_in_single_vect = rotate_normvect2nz(cur_direction_vect, temp_vect_in_single_vect)
                check_normvect.append(normed_mornvect)

                vect_single_filament_map_lst.append(temp_vect_in_single_vect) # rotate among cur filament coord
                dist.append(np.linalg.norm(temp_vect_in_single_vect))

        # filament_vect_2_MT_map_lst.append(vect_single_filament_map_lst) # vect to surrounding in a flat

    # print(filament_coord_extend_dict)
    # print(filament_vect_map_lst)
    # print('filament num', len(filament_vect_map_lst))


    return dist





def FA_pairdist(FAcoords, voxel_size_xyz):
    '''return dist '''
    avesize = np.average(voxel_size_xyz)
    dists = []
    for ii, cur_FAcoords in enumerate(FAcoords):
        rest_Facoords = FAcoords[ii:]
        if len(rest_Facoords) > 1:  # exclude last num
            rest_Facoords = rest_Facoords[1:]
            cur_dist = cdist([cur_FAcoords], rest_Facoords, metric='euclidean')
            cur_dist = list(cur_dist.reshape(1,-1))[0]
            cur_dist = [ i * avesize/10 for i in cur_dist ] #nm
            dists.extend(cur_dist)

    return dists



def FA2ca_dist(FAcoords, cacoords, voxel_size_xyz):
    '''return dist '''

    avesize = np.average(voxel_size_xyz)
    Dists = cdist(FAcoords, cacoords, metric='euclidean')

    Dists = list(Dists.reshape(1,-1))[0]
    Dists = [ i * avesize/10 for i in Dists ]  # nm
    # print(Dists)
    return Dists







#%%




def MTdistance2mem(mt_coords, mask, voxelsize):
    '''
    input: 
        mt coord: dict
        mask: array 
    '''
    filament_coords = mt_coords
    filament_coords_modified ={}    
    filamentnames = list(filament_coords.keys())
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords[f'{filamentnames[idx]}']
        filament_coords_modified[f'{filamentnames[idx]}'] =[]
        for i, coord_new in enumerate(curfilamentcoordslst):
            coord_new = [coord_new[0]/voxelsize[0], coord_new[1]/voxelsize[1], coord_new[2]/voxelsize[2]]
            filament_coords_modified[f'{filamentnames[idx]}'].append(coord_new)




    coords_lst = []
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords_modified[f'{filamentnames[idx]}']
        coords_lst.extend(curfilamentcoordslst)

    coords = np.array(coords_lst).reshape(-1,3)
    # print('line56', coords)
    voxel_coords = coords

    coords_x = [ voxel_coords[i][0] for i in range(len(voxel_coords))]
    coords_y = [ voxel_coords[i][1] for i in range(len(voxel_coords))]
    coords_z = [ voxel_coords[i][2] for i in range(len(voxel_coords))]

    # print(sorted(coords_x))
    
    # print(np.max(coords_x), np.max(coords_y), np.max(coords_z))
    # print(np.min(coords_x), np.min(coords_y), np.min(coords_z))
    
    print('mt coords range', [[np.round(np.min(coords_x)), np.round(np.max(coords_x))],[np.round(np.min(coords_y)), np.round(np.max(coords_y))],[np.round(np.min(coords_z)), np.round(np.max(coords_z))]] )
    print('mask size', mask.shape)
    print('coords in total', len(voxel_coords))


    assert np.max(coords_x) <= mask.shape[0] and np.min(coords_x)>= 0
    assert np.max(coords_y) <= mask.shape[1] and np.min(coords_y)>= 0
    assert np.max(coords_z) <= mask.shape[2] and np.min(coords_z)>= 0


    img01 = mask
    reverse_img01 = img01 * (-1) + 1
    img01_rev_edt = scipy.ndimage.distance_transform_edt(reverse_img01)

    avesize = np.average(voxelsize)
    dist = [ img01_rev_edt[int(voxel_coords[i][0])][int(voxel_coords[i][1])][int(voxel_coords[i][2])]  for i in range(len(voxel_coords)) ]
    dist = [dist[i] * avesize/10 for i in range(len(dist))] # nm
    print('distance num', len(dist))

    return dist


def MTdistance2mitomem(mt_coords, mask, voxelsize):
    '''
    input: 
        mt coord: dict
        mask: array 
    '''
    filament_coords = mt_coords
    filament_coords_modified ={}    
    filamentnames = list(filament_coords.keys())
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords[f'{filamentnames[idx]}']
        filament_coords_modified[f'{filamentnames[idx]}'] =[]
        for i, coord_new in enumerate(curfilamentcoordslst):
            coord_new = [coord_new[0]/voxelsize[0], coord_new[1]/voxelsize[1], coord_new[2]/voxelsize[2]]
            filament_coords_modified[f'{filamentnames[idx]}'].append(coord_new)




    coords_lst = []
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords_modified[f'{filamentnames[idx]}']
        coords_lst.extend(curfilamentcoordslst)

    coords = np.array(coords_lst).reshape(-1,3)
    # print('line56', coords)
    voxel_coords = coords

    coords_x = [ voxel_coords[i][0] for i in range(len(voxel_coords))]
    coords_y = [ voxel_coords[i][1] for i in range(len(voxel_coords))]
    coords_z = [ voxel_coords[i][2] for i in range(len(voxel_coords))]

    # print(sorted(coords_x))
    
    # print(np.max(coords_x), np.max(coords_y), np.max(coords_z))
    # print(np.min(coords_x), np.min(coords_y), np.min(coords_z))
    
    print('mt coords range', [[np.round(np.min(coords_x)), np.round(np.max(coords_x))],[np.round(np.min(coords_y)), np.round(np.max(coords_y))],[np.round(np.min(coords_z)), np.round(np.max(coords_z))]] )
    print('mask size', mask.shape)
    print('coords in total', len(voxel_coords))


    assert np.max(coords_x) <= mask.shape[0] and np.min(coords_x)>= 0
    assert np.max(coords_y) <= mask.shape[1] and np.min(coords_y)>= 0
    assert np.max(coords_z) <= mask.shape[2] and np.min(coords_z)>= 0


    img01 = mask
    reverse_img01 = img01 * (-1) + 1
    img01_rev_edt = scipy.ndimage.distance_transform_edt(reverse_img01)

    avesize = np.average(voxelsize)
    dist = [ img01_rev_edt[int(voxel_coords[i][0])][int(voxel_coords[i][1])][int(voxel_coords[i][2])]  for i in range(len(voxel_coords)) ]
    dist = [dist[i] * avesize/10 for i in range(len(dist))] # nm
    print('distance num', len(dist))

    return dist





#%%
# for distance to filaments

def rotate_ingrid(vect, 
            angle, 
            ):  # Euler-Rodrigues formula
    '''
    rotate in 3d 
    angle = [theta1, theta2, theta3]  shunshizhen  = [angle in yz, angle in xz, angle in xy]
     + shunshizhen
     - nishizhen
    '''

    sin1, cos1 = round(math.sin(angle[0]),3), round(math.cos(angle[0]),3)
    sin2, cos2 = round(math.sin(angle[1]),3), round(math.cos(angle[1]),3)
    sin3, cos3 = round(math.sin(angle[2]),3), round(math.cos(angle[2]),3)

    
    
    Rx = np.array([[ 1,     0,     0,   ],
                   [ 0,     cos1,  sin1,],
                   [ 0,    -sin1,  cos1,]])
    Ry = np.array([[ cos2,  0,    -sin2,],
                   [ 0,     1,     0,   ],
                   [ sin2,  0,     cos2,]])
    Rz = np.array([[ cos3,  sin3,  0,   ],
                   [-sin3,  cos3,  0,   ],
                   [ 0,     0,     1,   ]])    

    
    return np.dot(Rx, np.dot(Ry, np.dot(Rz, vect)))


def rotate_normvect2z(coord0, vect2):
    '''
    rotate normal vector to z axis, then the vect in the same way 
    '''
    # print(coord0,vect2)
    threshold_d = 0.05
    angle1 = [-np.arctan((coord0[1])/(coord0[2]+1e-6)), 0, 0]
    aa = rotate_ingrid(coord0, angle1)
    angle2 = [0, np.arctan(aa[0] / (aa[2]+1e-6)), 0]
    aa = rotate_ingrid(aa, angle2)
    aa = aa/ np.linalg.norm(aa)
    # print('line347',aa)
    if not abs(aa[0])<threshold_d or not abs(aa[1])<threshold_d :
        print('line347-1',aa)
    # print('line347-0-1',aa)
    assert abs(aa[0])<threshold_d
    assert abs(aa[1])<threshold_d
    bb = rotate_ingrid(vect2, angle1)
    bb = rotate_ingrid(bb, angle2)
    return aa, bb

def rotate_normvect2nz(coord0, vect2):
    '''
    rotate normal vector to z axis, then the vect in the same way 
    '''
    # print(coord0,vect2)
    threshold_d = 0.05
    angle1 = [-np.arctan((coord0[1])/(coord0[2]+1e-6)), 0, 0]
    aa = rotate_ingrid(coord0, angle1)
    angle2 = [0, np.arctan(aa[0] / (aa[2]+1e-6))+np.pi, 0]
    aa = rotate_ingrid(aa, angle2)
    aa = aa/ np.linalg.norm(aa)
    # print('line367',aa)
    if not abs(aa[0])<threshold_d or not abs(aa[1])<threshold_d :
        print('line367-2',aa)
    assert abs(aa[0])<threshold_d
    assert abs(aa[1])<threshold_d
    bb = rotate_ingrid(vect2, angle1)
    bb = rotate_ingrid(bb, angle2)

    return aa, bb




def actinangle2dist(filament_coord_extend_dict, imgsize, voxel_size_xyz):
    'input: filaments distance dict'
    edge = 0.
    x_bound = [0 + edge, imgsize[0]*voxel_size_xyz[0] - edge]
    y_bound = [0 + edge, imgsize[1]*voxel_size_xyz[1] - edge]
    z_bound = [0 + edge, imgsize[2]*voxel_size_xyz[2] - edge]

    filament_coord_extend_lst = []
    filamentnames = list(filament_coord_extend_dict.keys())

    for idx in range(len(filamentnames)):
        # curfilamentcoordslst = filament_coord_extend_dict[f'{filamentnames[idx]}']
        curfilamentcoordslst = [ np.array(coord) for coord in filament_coord_extend_dict[f'{filamentnames[idx]}'] ]
        filament_coord_extend_lst.append(curfilamentcoordslst)

    filament_coord_lst = filament_coord_extend_lst
    filaments_vect_lst = []  # vect for every point in every filament

    for single_filament_coord in filament_coord_lst:
        filament_vectlst = []
        for jj, point in enumerate(single_filament_coord):
            if jj < 1:
                tempvect = single_filament_coord[jj+1] - point
            else:
                tempvect = point- single_filament_coord[jj-1]
            # unit_vect = vect / np.linalg.norm(vect)
            filament_vectlst.append(tempvect)  
        filaments_vect_lst.append(filament_vectlst)  # not unit vect 


    filaments_min_d_lst = []
    filaments_angle_lst = []

    for i, curfilaments_coords in enumerate(filament_coord_extend_lst):
        curfilaments_vect = filaments_vect_lst[i]

        rest_filament_coord_extend_lst = copy.deepcopy(filament_coord_extend_lst)
        rest_filament_vect_extend_lst = copy.deepcopy(filaments_vect_lst)
        _ = rest_filament_coord_extend_lst.pop(i)
        _ = rest_filament_vect_extend_lst.pop(i)

        # curfilaments_coords = [point1, point2, ...]
        filament_dist_lst = []
        filament_angle_lst = []
        for ll, point in enumerate(curfilaments_coords):
            vect1 = curfilaments_vect[ll]
            min_dist = np.inf
            cur_angle = np.inf
            for kk, filaments_coords in enumerate(rest_filament_coord_extend_lst):
                distss = cdist([point], filaments_coords, metric='euclidean')
                min_d = np.min(distss)
                min_d_idx = np.where(distss == min_d)[1][0]

                if min_d < min_dist:
                    min_dist = min_d
                    vect2 = rest_filament_vect_extend_lst[kk][min_d_idx]
                    vector_dot_product = np.dot(vect1, vect2)
                    arccos = np.arccos(vector_dot_product / (np.linalg.norm(vect1) * np.linalg.norm(vect2)))
                    angle = np.degrees(arccos)
                    if angle >= 180:
                        cur_angle = 0
                    elif angle > 90:
                        cur_angle = 180 - angle
                    else:
                        cur_angle = angle

            filament_dist_lst.append(min_dist)
            filament_angle_lst.append(cur_angle)

        filaments_min_d_lst.extend(filament_dist_lst)
        filaments_angle_lst.extend(filament_angle_lst)


    return  filaments_min_d_lst, filaments_angle_lst


   
def actinangle2dist(filament_coord_extend_dict, imgsize, voxel_size_xyz):
    'input: filaments distance dict'
    edge = 0.
    x_bound = [0 + edge, imgsize[0]*voxel_size_xyz[0] - edge]
    y_bound = [0 + edge, imgsize[1]*voxel_size_xyz[1] - edge]
    z_bound = [0 + edge, imgsize[2]*voxel_size_xyz[2] - edge]

    filament_coord_extend_lst = []
    filamentnames = list(filament_coord_extend_dict.keys())

    for idx in range(len(filamentnames)):
        # curfilamentcoordslst = filament_coord_extend_dict[f'{filamentnames[idx]}']
        curfilamentcoordslst = [ np.array(coord) for coord in filament_coord_extend_dict[f'{filamentnames[idx]}'] ]
        filament_coord_extend_lst.append(curfilamentcoordslst)

    filament_coord_lst = filament_coord_extend_lst
    filaments_vect_lst = []  # vect for every point in every filament

    for single_filament_coord in filament_coord_lst:
        filament_vectlst = []
        for jj, point in enumerate(single_filament_coord):
            if jj < 1:
                tempvect = single_filament_coord[jj+1] - point
            else:
                tempvect = point- single_filament_coord[jj-1]
            # unit_vect = vect / np.linalg.norm(vect)
            filament_vectlst.append(tempvect)  
        filaments_vect_lst.append(filament_vectlst)  # not unit vect 


    filaments_min_d_lst = []
    filaments_angle_lst = []
    if len(filament_coord_extend_lst) < 2:
        return  filaments_min_d_lst, filaments_angle_lst
    else:
        for i, curfilaments_coords in enumerate(filament_coord_extend_lst):
            curfilaments_vect = filaments_vect_lst[i]

            rest_filament_coord_extend_lst = copy.deepcopy(filament_coord_extend_lst)
            rest_filament_vect_extend_lst = copy.deepcopy(filaments_vect_lst)
            # rest_filament_coord_extend_lst = [ filament_coord_extend_lst[ii] for ii in range(len(filament_coord_extend_lst)) if ii != i ]
            # rest_filament_vect_extend_lst = [ filaments_vect_lst[ii] for ii in range(len(filaments_vect_lst)) if ii != i ]

            
            _ = rest_filament_coord_extend_lst.pop(i)
            _ = rest_filament_vect_extend_lst.pop(i)

            # curfilaments_coords = [point1, point2, ...]
            filament_dist_lst = []
            filament_angle_lst = []
            for ll, point in enumerate(curfilaments_coords):
                vect1 = curfilaments_vect[ll]
                min_dist = np.inf
                cur_angle = np.inf
                for kk, filaments_coords in enumerate(rest_filament_coord_extend_lst):
                    distss = cdist([point], filaments_coords, metric='euclidean')
                    min_d = np.min(distss)
                    min_d_idx = np.where(distss == min_d)[1][0]

                    if min_d < min_dist:
                        min_dist = min_d
                        vect2 = rest_filament_vect_extend_lst[kk][min_d_idx]
                        vector_dot_product = np.dot(vect1, vect2)
                        arccos = np.arccos(vector_dot_product / (np.linalg.norm(vect1) * np.linalg.norm(vect2)))
                        angle = np.degrees(arccos)
                        if angle >= 180:
                            cur_angle = 0
                        elif angle > 90:
                            cur_angle = 180 - angle
                        else:
                            cur_angle = angle
                
                filament_dist_lst.append(min_dist)
                filament_angle_lst.append(cur_angle)

            filaments_min_d_lst.extend(filament_dist_lst)
            filaments_angle_lst.extend(filament_angle_lst)


    return  filaments_min_d_lst, filaments_angle_lst



def actin2Fa_angle(filament_coords, fa_data, voxelsize, dist_threshold, voxel_size_xyz):
    
    '''
    input: 
        filament coord: dict
        facoord: dataframe, on image 
        mask: array 
        dist_threshold:  nm
    '''

    angleslst = []
    avevoxelsize = np.average(voxelsize)

    fa_coord = [[fa_data['X'][i], fa_data['Y'][i], fa_data['Z'][i] ] for i in range(fa_data.shape[0]) ]
    fa_coord_idx = [i for i in range(len(fa_coord))]

    filament_coords_modified ={}    
    filamentnames = list(filament_coords.keys())
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords[f'{filamentnames[idx]}']
        filament_coords_modified[f'{filamentnames[idx]}'] =[]
        for i, coord_new in enumerate(curfilamentcoordslst):
            coord_new = [coord_new[0]/voxelsize[0], coord_new[1]/voxelsize[1], coord_new[2]/voxelsize[2]]
            filament_coords_modified[f'{filamentnames[idx]}'].append(coord_new)  ## coord in axis


    filament_coords_byfilament = list(filament_coords_modified.values())

    filament_coord_extend_lst = []
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords_modified[f'{filamentnames[idx]}']
        filament_coord_extend_lst.extend(curfilamentcoordslst)

    filament_extend_coords = np.array(filament_coord_extend_lst).reshape(-1,3)

    # print('line56', coords)
    voxel_coords = filament_extend_coords

    coords_x = [ voxel_coords[i][0] for i in range(len(voxel_coords))]
    coords_y = [ voxel_coords[i][1] for i in range(len(voxel_coords))]
    coords_z = [ voxel_coords[i][2] for i in range(len(voxel_coords))]



    fa_actin_match_dict = dict()

    for i in range(len(fa_coord)):
        fa_actin_match_dict[f'fa_{i}'] = list()


    for i, curfilaments_coords in enumerate(filament_coords_byfilament):
        'assign actin to focal adhesion'
        cur_filamentname = filamentnames[i]
        curfilaments_coords = [curfilaments_coords[0], curfilaments_coords[-1]]  # obtain two ends


        dists_matrix = cdist(curfilaments_coords, fa_coord, metric='euclidean')
        # print('matrix shape',dists_matrix.shape)
        min_d = np.min(dists_matrix)
        
        # print(np.where(dists_matrix == min_d))
        min_d_idx = np.where(dists_matrix == min_d)[1][0]  ## get row for fa and num

        if min_d <= (dist_threshold/np.average(voxel_size_xyz)*10):
        # if True:
            cur_fa_coord = fa_coord[min_d_idx]
            cur_fa_coord_idx = fa_coord_idx[min_d_idx]
            fa_actin_match_dict[f'fa_{cur_fa_coord_idx}'].append(i)
        

    fa_actin_angle_lst =[]
    fa_actin_dist_lst = []
    fa_actin_actincoords = []
    fa_actin_actinidx = []
    fa_actin_faidx = []
    for i ,cur_fa_coord in enumerate(fa_coord):
        '''calculate angle for each filament'''
        curfa_filamentlsts = fa_actin_match_dict[f'fa_{i}']
        if len(curfa_filamentlsts) > 0:  
            for fila_idx in curfa_filamentlsts:
                cur_filamentcoords = filament_coords_byfilament[fila_idx]
                curfilaments_coords = cur_filamentcoords
                
                dists_matrix = cdist([curfilaments_coords[0], curfilaments_coords[-1]], [cur_fa_coord], metric='euclidean')
                min_d = np.min(dists_matrix)
                min_d_idx = np.where(dists_matrix == min_d)[0][0]


                if min_d_idx == 0:
                    filament_end_point = curfilaments_coords[0]
                    filament_end_vect = np.array(curfilaments_coords[1]) - np.array(curfilaments_coords[0])
                elif min_d_idx == 1:
                    filament_end_point = curfilaments_coords[-1]
                    filament_end_vect = np.array(curfilaments_coords[-2]) - np.array(curfilaments_coords[-1])


                fila_to_fa_vect = filament_end_point - np.array(cur_fa_coord)

                vector_dot_product = np.dot(fila_to_fa_vect, filament_end_vect)
                arccos = np.arccos(vector_dot_product / (np.linalg.norm(fila_to_fa_vect) * np.linalg.norm(filament_end_vect)))
                angle = np.degrees(arccos)
                if angle >= 180:  # include nan
                    cur_angle = 0
                # elif angle > 90:
                #     cur_angle = 180 - angle
                else:
                    cur_angle = angle

                fa_actin_angle_lst.append(cur_angle)
                fa_actin_dist_lst.append(min_d * avevoxelsize/10) # nm
                fa_actin_actincoords.append([curfilaments_coords[0], curfilaments_coords[-1]][min_d_idx])
                fa_actin_actinidx.append(fila_idx)
                fa_actin_faidx.append(i)
            # currespond_angle_for_filament[i] = angle



    return fa_actin_angle_lst, fa_actin_dist_lst




def actin2Fa_angle_check(filament_coords, fa_data, voxelsize, dist_threshold):
    
    '''
    input: 
        filament coord: dict
        facoord: dataframe, on image 
        mask: array 
    '''
    angleslst = []

    fa_coord = [[fa_data['X'][i], fa_data['Y'][i], fa_data['Z'][i] ] for i in range(fa_data.shape[0]) ]
    fa_coord_idx = [i for i in range(len(fa_coord))]

    filament_coords_modified ={}    
    filamentnames = list(filament_coords.keys())
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords[f'{filamentnames[idx]}']
        filament_coords_modified[f'{filamentnames[idx]}'] =[]
        for i, coord_new in enumerate(curfilamentcoordslst):
            coord_new = [coord_new[0]/voxelsize[0], coord_new[1]/voxelsize[1], coord_new[2]/voxelsize[2]]
            filament_coords_modified[f'{filamentnames[idx]}'].append(coord_new)  ## coord in axis


    filament_coords_byfilament = list(filament_coords_modified.values())

    filament_coord_extend_lst = []
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords_modified[f'{filamentnames[idx]}']
        filament_coord_extend_lst.extend(curfilamentcoordslst)

    filament_extend_coords = np.array(filament_coord_extend_lst).reshape(-1,3)

    # print('line56', coords)
    voxel_coords = filament_extend_coords

    coords_x = [ voxel_coords[i][0] for i in range(len(voxel_coords))]
    coords_y = [ voxel_coords[i][1] for i in range(len(voxel_coords))]
    coords_z = [ voxel_coords[i][2] for i in range(len(voxel_coords))]



    fa_actin_match_dict = dict()

    for i in range(len(fa_coord)):
        fa_actin_match_dict[f'fa_{i}'] = list()


    for i, curfilaments_coords in enumerate(filament_coords_byfilament):
        'assign actin to focal adhesion'
        cur_filamentname = filamentnames[i]
        curfilaments_coords = curfilaments_coords

        dists_matrix = cdist(curfilaments_coords, fa_coord, metric='euclidean')

        min_d = np.min(dists_matrix)
        # print(np.where(dists_matrix == min_d))
        min_d_idx = np.where(dists_matrix == min_d)[1][0]  ## get row for fa and num

        # if min_d <= (dist_threshold/np.mean(voxel_size_xyz)):
        if True:
            cur_fa_coord = fa_coord[min_d_idx]
            cur_fa_coord_idx = fa_coord_idx[min_d_idx]
            fa_actin_match_dict[f'fa_{cur_fa_coord_idx}'].append(i)
        

    fa_actin_angle_lst =[]
    fa_actin_dist_lst = []
    fa_actin_actincoords = []
    fa_actin_actinidx = []
    fa_actin_faidx = []
    for i ,cur_fa_coord in enumerate(fa_coord):
        '''calculate angle for each filament'''
        curfa_filamentlsts = fa_actin_match_dict[f'fa_{i}']

        if len(curfa_filamentlsts) > 0:
            
            for fila_idx in curfa_filamentlsts:
                cur_filamentcoords = filament_coords_byfilament[fila_idx]
                
                dists_matrix = cdist([cur_filamentcoords[0], cur_filamentcoords[-1]], [cur_fa_coord], metric='euclidean')
                min_d = np.min(dists_matrix)


                min_d_idx = np.where(dists_matrix == min_d)[0][0]


                if min_d_idx == 0:
                    filament_end_point = curfilaments_coords[0]
                    filament_end_vect = np.array(curfilaments_coords[1]) - np.array(curfilaments_coords[0])
                elif min_d_idx == 1:
                    filament_end_point = curfilaments_coords[-1]
                    filament_end_vect = np.array(curfilaments_coords[-2]) - np.array(curfilaments_coords[-1])


                fila_to_fa_vect = filament_end_point - np.array(cur_fa_coord)

                vector_dot_product = np.dot(fila_to_fa_vect, filament_end_vect)
                arccos = np.arccos(vector_dot_product / (np.linalg.norm(fila_to_fa_vect) * np.linalg.norm(filament_end_vect)))
                angle = np.degrees(arccos)
                if angle >= 180:  # include nan
                    cur_angle = 0
                # elif angle > 90:
                #     cur_angle = 180 - angle
                else:
                    cur_angle = angle

                fa_actin_angle_lst.append(cur_angle)
                fa_actin_dist_lst.append(min_d)
                fa_actin_actincoords.append([curfilaments_coords[0], curfilaments_coords[-1]][min_d_idx])
                fa_actin_actinidx.append(fila_idx)
                fa_actin_faidx.append(i)
            # currespond_angle_for_filament[i] = angle



    return fa_actin_angle_lst, fa_actin_dist_lst





def Fa_aggragated_actinactin_angle(filament_coords, filament_angles, fa_data, voxelsize, dist_threshold):
    
    '''
    input: 
        filament coord: dict
        filament_angle_coords: dict 
        facoord: dataframe, coords on image 
        mask: array 
    '''

    angleslst = []

    fa_coord = [[fa_data['X'][i], fa_data['Y'][i], fa_data['Z'][i] ] for i in range(fa_data.shape[0]) ]
    fa_coord_idx = [i for i in range(len(fa_coord))]

    filament_coords_modified ={}    
    filamentnames = list(filament_coords.keys())
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords[f'{filamentnames[idx]}']
        filament_coords_modified[f'{filamentnames[idx]}'] =[]
        for i, coord_new in enumerate(curfilamentcoordslst):
            coord_new = [coord_new[0]/voxelsize[0], coord_new[1]/voxelsize[1], coord_new[2]/voxelsize[2]]
            filament_coords_modified[f'{filamentnames[idx]}'].append(coord_new)  ## coord in axis


    filament_coords_byfilament = list(filament_coords_modified.values())
    filament_angle_byfilament = list(filament_angles.values())

    filament_coord_extend_lst = []
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords_modified[f'{filamentnames[idx]}']
        filament_coord_extend_lst.extend(curfilamentcoordslst)




    fa_actin_match_dict = dict()

    for i in range(len(fa_coord)):
        fa_actin_match_dict[f'fa_{i}'] = list()


    for i, curfilaments_coords in enumerate(filament_coords_byfilament):
        'assign actin to focal adhesion'
        cur_filamentname = filamentnames[i]
        curfilaments_coords = curfilaments_coords

        dists_matrix = cdist(curfilaments_coords, fa_coord, metric='euclidean')

        min_d = np.min(dists_matrix)
        # print(np.where(dists_matrix == min_d))
        min_d_idx = np.where(dists_matrix == min_d)[1][0]  ## get row for fa and num

        if min_d <= (dist_threshold/np.average(voxelsize)*10):
        # if True:
            cur_fa_coord = fa_coord[min_d_idx]
            cur_fa_coord_idx = fa_coord_idx[min_d_idx]
            fa_actin_match_dict[f'fa_{cur_fa_coord_idx}'].append(i)
        



    actin_actin_angle_lst =[]
    actin_mem_angle_lst = []
    actin_actin_actinidx = []

    for i ,cur_fa_coord in enumerate(fa_coord):
        '''calculate angle for each filament'''
        curfa_filamentlsts = fa_actin_match_dict[f'fa_{i}']  # {'fa_1': filament1, filament2}

        if len(curfa_filamentlsts) == 1:
            fila_idx = curfa_filamentlsts[0]
            cur_filamentcoords = filament_coords_byfilament[fila_idx]
            cur_actin_mem_angle = filament_angle_byfilament[fila_idx]
            
            # actin_actin_angle_lst.append(cur_actin_mem_angle)
            # actin_mem_angle_lst.append(cur_actin_mem_angle)
            # actin_actin_actinidx.append(fila_idx)

        if len(curfa_filamentlsts) > 1:
            

            for fila_idx in range(len(curfa_filamentlsts)-1):

                cur_filamentcoords1 = filament_coords_byfilament[fila_idx]
                cur_actin_mem_angle = filament_angle_byfilament[fila_idx]
                vect1 = np.array(cur_filamentcoords1[0]) - np.array(cur_filamentcoords1[-1])


                for fila2_idx in range(fila_idx, len(curfa_filamentlsts)):
                    # print(fila_idx, fila2_idx, len(curfa_filamentlsts))
                    cur_filamentcoords2 = filament_coords_byfilament[fila2_idx]
                    vect2 = np.array(cur_filamentcoords2[0]) - np.array(cur_filamentcoords2[-1])

                    vector_dot_product = np.dot(vect1, vect2)
                    arccos = np.arccos(vector_dot_product / (np.linalg.norm(vect1) * np.linalg.norm(vect2)))
                    angle = np.degrees(arccos)
                    if angle >= 180:  # include nan
                        # print(111, angle)
                        cur_angle = 0
                    elif angle > 90:
                        cur_angle = 180 - angle
                    else:
                        cur_angle = angle
                    
                    actin_actin_angle_lst.append(cur_angle)
                    actin_mem_angle_lst.append(cur_actin_mem_angle)
                    actin_actin_actinidx.append(fila_idx)


                    


    return actin_actin_angle_lst, actin_mem_angle_lst





# def filament angle and dist 

#%%

def filadistance2mem_vis(filament_coords, mask, voxelsize):
    '''
    input: 
        filament coord: dict
        mask: array 
    '''

    filament_coords_modified ={}    
    filamentnames = list(filament_coords.keys())
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords[f'{filamentnames[idx]}']
        filament_coords_modified[f'{filamentnames[idx]}'] =[]
        for i, coord_new in enumerate(curfilamentcoordslst):
            coord_new = [coord_new[0]/voxelsize[0], coord_new[1]/voxelsize[1], coord_new[2]/voxelsize[2]]
            filament_coords_modified[f'{filamentnames[idx]}'].append(coord_new)




    coords_lst = []
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords_modified[f'{filamentnames[idx]}']
        coords_lst.extend(curfilamentcoordslst)

    coords = np.array(coords_lst).reshape(-1,3)
    # print('line56', coords)
    voxel_coords = coords

    coords_x = [ voxel_coords[i][0] for i in range(len(voxel_coords))]
    coords_y = [ voxel_coords[i][1] for i in range(len(voxel_coords))]
    coords_z = [ voxel_coords[i][2] for i in range(len(voxel_coords))]

    # print(sorted(coords_x))
    
    # print(np.max(coords_x), np.max(coords_y), np.max(coords_z))
    # print(np.min(coords_x), np.min(coords_y), np.min(coords_z))
    
    print('coords range', [[np.round(np.min(coords_x)), np.round(np.max(coords_x))],[np.round(np.min(coords_y)), np.round(np.max(coords_y))],[np.round(np.min(coords_z)), np.round(np.max(coords_z))]] )
    print('mask size', mask.shape)
    print('coords in total', len(voxel_coords))


    assert np.max(coords_x) <= mask.shape[0] and np.min(coords_x)>= 0
    assert np.max(coords_y) <= mask.shape[1] and np.min(coords_y)>= 0
    assert np.max(coords_z) <= mask.shape[2] and np.min(coords_z)>= 0


    img01 = mask
    reverse_img01 = img01 * (-1) + 1
    img01_rev_edt = scipy.ndimage.distance_transform_edt(reverse_img01)

    filtered_map = np.zeros_like(img01_rev_edt)
    for i in range(len(voxel_coords)):
        filtered_map[int(voxel_coords[i][0])][int(voxel_coords[i][1])][int(voxel_coords[i][2])] = 1  
    print('should be 1,', np.max(filtered_map))

    filtered_map = scipy.ndimage.morphology.binary_dilation(filtered_map, iterations= 3)
    img01_rev_edt2 = img01_rev_edt + 50
    
    fila_dist_map = filtered_map * img01_rev_edt2
    # fila_dist_map = fila_dist_map + img01 * (-1)


    # avesize = np.average(voxelsize)
    # dist = [ img01_rev_edt[int(voxel_coords[i][0])][int(voxel_coords[i][1])][int(voxel_coords[i][2])]  for i in range(len(voxel_coords)) ]
    # dist = [dist[i] * avesize/10 for i in range(len(dist))] # nm
    # print('distance num', len(dist))

    return fila_dist_map


