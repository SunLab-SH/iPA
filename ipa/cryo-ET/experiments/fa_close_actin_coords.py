# %%
import os, sys
import re
from tabnanny import check
from importlib_metadata import method_cache
import seaborn as sns

filedir = os.path.dirname(os.path.abspath(__file__))
os.chdir(filedir)
sys.path.append(os.path.pardir)
print(os.getcwd())

import numpy as np
import threading
from common.parser import arg
from common.preprocess import *
import glob 
from multiprocessing.dummy import Pool as ThreadPool
from common import preprocess, distances_generator, outplot
from scipy.spatial.distance import cdist
import torch 
from process_organelle_data import process_organelles
# import cv2

import seaborn as sns
print(torch.__version__)
print(torch.cuda.is_available())

device = torch.device('cuda')

#%%


def getsingledatafilename():
    filename =  'focal_aggregated_actin_coord.csv' 
    plotname =  'focal_aggregated_actin_coord.png'

    return filename, plotname



def actin2Fa_coord(filament_coords, fa_data, voxelsize, dist_threshold, voxel_size_xyz):
    
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
        

    fa_coordlst = list()
    actin_coordlst = list()
    min_distlst =list()

    for i ,cur_fa_coord in enumerate(fa_coord):
        '''calculate distance for each filament'''
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
                elif min_d_idx == 1:
                    filament_end_point = curfilaments_coords[-1]

                fa_coordlst.append(cur_fa_coord)
                actin_coordlst.append(filament_end_point)
                min_distlst.append(min_d)


    #sort
    for ii in range(len(min_distlst)):
        for jj in range(ii, len(min_distlst)):
            if min_distlst[ii] > min_distlst[jj]:
                min_distlst[ii], min_distlst[jj] = min_distlst[jj], min_distlst[ii]
                fa_coordlst[ii], fa_coordlst[jj] = fa_coordlst[jj], fa_coordlst[ii]
                actin_coordlst[ii], actin_coordlst[jj] = actin_coordlst[jj], actin_coordlst[ii]



    return fa_coordlst, actin_coordlst, min_distlst





class focal_aggregated_actin_coord(process_organelles):
    '''

    '''
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'focal_aggregated_actin_coord'
        self.dist_threshold = 100 # nm
        self.focal_aggregated_actin_coord_filename,  self.focal_aggregated_actin_coord_plotname = getsingledatafilename()



    def step(self):
        organelletypelst = [self.focal_adhesion_type, self.actin_type]
        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} {self.name} data not exist.')
        else:
            if not self.check_data_updated(organelletypelst):
                # print(11)
                self.generate_processed_data()

            # self.focal_aggregated_actin_coord_check()
            self.focal_aggregated_actin_coord_plot(overlap = self.plotdata_overlap)
        print(f'Done with {self.name}.')

    def focal_aggregated_actin_coord_plot(self, overlap):
        # print(1111)
        self.focal_aggregated_actin_coord_map_filepath = f'{self.datanewpath}/{self.dataid}_{self.focal_aggregated_actin_coord_filename}'
        self.curdatapath = self.focal_aggregated_actin_coord_map_filepath
        if not os.path.exists(self.curdatapath) or overlap == True:
            self.focal_aggregated_actin_coord_map()


        # focal_aggregated_actin_map_filename

        with open(f'{self.actinfilepath}', 'r') as f:
            actin_data = json.load(f)
        fa_data = pd.read_csv(self.focaladhesionfilepath)


        coords_pairs = pd.read_csv(self.curdatapath) ## angles

        numlst =  coords_pairs['facoord']
        print('num', coords_pairs.shape[0])
        print(coords_pairs)

        # outplot.fa_actin_angle_dist_map_plot(angle_dist_pairs, self.dataid,)
        # plotname = f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.focal_aggregated_actin_coord_plotname}'
        # plt.savefig(plotname)
        # plt.show()
        # plt.close()
     

    def focal_aggregated_actin_coord_map(self):
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)
        
        with open(self.actinfilepath, 'r') as f:
            actin_data = json.load(f)
        fa_data = pd.read_csv(self.focaladhesionfilepath)

        actin_data = distances_generator.shift_bias(actin_data, shift_xml_xyz)

        facoord, actincoord, min_distlst = actin2Fa_coord(actin_data, fa_data, voxel_size_xyz, self.dist_threshold, voxel_size_xyz)  # in nm
        # anglelst, distlst = self.actin2Fa_angle(actin_data, fa_data, voxel_size_xyz, self.dist_threshold)

        coord_savepd = pd.DataFrame({'facoord':facoord, 'actincoord':actincoord, 'dist':min_distlst})  

        # filament2filament_dist_lst_savepd.loc['distance'] = filament2filament_dist_lst
        coord_savepd.to_csv(f'{self.datanewpath}/{self.dataid}_{self.focal_aggregated_actin_coord_filename}', index=False)




#%%


def single_main():
    dataidlst = preprocess.read_dataid()
    print(dataidlst)
    # dataidlst = dataidlst[:1]

    for dataid in dataidlst:
        print(dataid)

        faan = focal_aggregated_actin_coord(dataid)
        faan.plotdata_overlap = True
        faan.step()


if __name__ == '__main__':
    single_main()






#%%