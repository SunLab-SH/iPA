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
    filename =  'fa_surround_actin_num.csv' 
    plotname =  'fa_surround_actin_num.png'

    return filename, plotname



def actin2Fa_num(filament_coords, fa_data, voxelsize, dist_threshold, voxel_size_xyz):
    
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
        

    
    fa_actin_num = list()

    for i ,cur_fa_coord in enumerate(fa_coord):
        '''calculate angle for each filament'''
        curfa_filamentlsts = fa_actin_match_dict[f'fa_{i}']
        cur_actin_num = len(curfa_filamentlsts)
        fa_actin_num.append(cur_actin_num)


    return fa_actin_num





class focal_aggregated_actin_num(process_organelles):
    '''
    Actins aggregated within one focal adhesion are regareded as one cluster.
    calculate the angle distribution of the actin end to focal adhesion angle.
    distributions.
    '''
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'focal_aggregated_actin_num'
        self.dist_threshold = 100 # nm
        self.focal_aggregated_actin_num_filename,  self.focal_aggregated_actin_num_plotname = getsingledatafilename()



    def step(self):
        organelletypelst = [self.focal_adhesion_type, self.actin_type]
        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} {self.name} data not exist.')
        else:
            if not self.check_data_updated(organelletypelst):
                # print(11)
                self.generate_processed_data()

            # self.focal_aggregated_actin_num_check()
            self.focal_aggregated_actin_num_plot(overlap = self.plotdata_overlap)
        print(f'Done with {self.name}.')

    def focal_aggregated_actin_num_plot(self, overlap):
        # print(1111)
        self.focal_aggregated_actin_num_map_filepath = f'{self.datanewpath}/{self.dataid}_{self.focal_aggregated_actin_num_filename}'
        self.curdatapath = self.focal_aggregated_actin_num_map_filepath
        if not os.path.exists(self.curdatapath) or overlap == True:
            self.focal_aggregated_actin_num_map()


        # focal_aggregated_actin_map_filename

        with open(f'{self.actinfilepath}', 'r') as f:
            actin_data = json.load(f)
        fa_data = pd.read_csv(self.focaladhesionfilepath)


        num_dist_pairs = pd.read_csv(self.curdatapath) ## angles

        numlst =  num_dist_pairs['Num']
        print(num_dist_pairs)

        # outplot.fa_actin_angle_dist_map_plot(angle_dist_pairs, self.dataid,)
        # plotname = f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.focal_aggregated_actin_num_plotname}'
        # plt.savefig(plotname)
        # plt.show()
        # plt.close()
     

    def focal_aggregated_actin_num_map(self):
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)
        
        with open(self.actinfilepath, 'r') as f:
            actin_data = json.load(f)
        fa_data = pd.read_csv(self.focaladhesionfilepath)

        actin_data = distances_generator.shift_bias(actin_data, shift_xml_xyz)

        numlst = actin2Fa_num(actin_data, fa_data, voxel_size_xyz, self.dist_threshold, voxel_size_xyz)  # in nm
        # anglelst, distlst = self.actin2Fa_angle(actin_data, fa_data, voxel_size_xyz, self.dist_threshold)

        angle_savepd = pd.DataFrame({'Num':numlst})  

        # filament2filament_dist_lst_savepd.loc['distance'] = filament2filament_dist_lst
        angle_savepd.to_csv(f'{self.datanewpath}/{self.dataid}_{self.focal_aggregated_actin_num_filename}', index=False)




#%%


def single_main():
    dataidlst = preprocess.read_dataid()
    print(dataidlst)

    for dataid in dataidlst:
        print(dataid)

        faan = focal_aggregated_actin_num(dataid)
        faan.plotdata_overlap = True
        faan.step()


if __name__ == '__main__':
    single_main()




#%%