#%%
from concurrent.futures import process
from math import nan
import os, sys
import re
from tabnanny import check
from importlib_metadata import method_cache
import seaborn as sns
# from .process_organelle_data import focal_to_ca

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
import copy
# import cv2

print(torch.__version__)
print(torch.cuda.is_available())

device = torch.device('cuda')

#%%

class process_organelles():
    def __init__(self, dataid) -> None:
        self.dataid = dataid

        self.vesicle_type = 0
        self.actin_type = 1
        self.microtube_type = 2
        self.mito_type = 3
        self.endo_reticulum_type = 4
        self.focal_adhesion_type = 5
        self.ca_type = 6

        self.vesiclefilename = 'vesicle.mrc'
        self.actinfilename ='actin_filled_points.json'
        self.microtubefilename = 'microtube_filled_points.json'
        self.mitofilename ='mito.mrc'
        self.endoreticulumfilename = 'endoreticulum.mrc'
        self.focaladhesionfilename = 'focaladhesion_points.csv'
        self.focaladhesion_imagefilename =  'focaladhesion.mrc'
        self.cafilename = 'calcium_points.csv'
        self.ca_imagefilename = 'calcium.mrc'
        self.focal_to_ca_distance_filename = 'focal_to_ca_distance.csv'
        
        self.actinanglefilename = 'actin_to_pm_angles.json'
        self.plotdata_overlap = False

        self.get_datanewpath()
        self.load_plotdataname()

    def load_datatype(self, file_idx):
        'load data type to see what type of organelle is available for calculation'
        df = pd.read_excel(os.path.join(arg.root_dir, arg.parameter_dir, 'Organelle_Sum.xlsx'))
        # print(df)
        organelles = df.values
        # print(organelles)
        for line in organelles:
            if file_idx in line:
                organelle_label = line[1:]
        assert len(organelle_label) > 0
        ## ISG	Actin	MT	Mito  ER
        return organelle_label
    
    def check_type(self, organelle1, organelle2):
        'check if organelle type is able to make further calculation'
        organelle_type = self.load_datatype(self.dataid)
        print(organelle_type)
        if organelle_type[organelle1] == 1 and organelle_type[organelle2] == 1:
            return True
        else: return False

    def check_organelle_type(self,organelletype):
        # ISG	Actin	MT	Mito ER fa ca
        return self.check_type(organelletype[0], organelletype[1])

    def get_filelist(self, dir, Filelist):
        newDir = dir
        if os.path.isfile(dir):
            Filelist.append(dir)
            # # 若只是要返回文件文，使用这个
            # Filelist.append(os.path.basename(dir))
        elif os.path.isdir(dir):
            for s in os.listdir(dir):
                # 如果需要忽略某些文件夹，使用以下代码
                #if s == "xxx":
                    #continue
                newDir=os.path.join(dir,s)
                self.get_filelist(newDir, Filelist)
        return Filelist


    def check_data_updated(self, organelletypelist):
        '''
        check if data in updated files
        '''

        self.organelletypelist = organelletypelist
        # organellenamelst = [self.vesiclefilename,self.actinfilename, self.microtubefilename, self.mitofilename, self.endoreticulumfilename]
        filelist = self.get_filelist(os.path.join(arg.root_dir, arg.processed_data_dir),[])
        filelist = [ file for file in filelist if self.dataid in file.lower()] 
        # print(filelist)

        self.vesiclefilepath, self.actinfilepath, self.microtubefilepath, self.mitofilepath, self.endoreticulumfilepath = str(),str(),str(),str(),str()
        self.focaladhesionfilepath, self.actinanglefilepath = str(), str()
        self.cafilepath = str()

        if 0 in self.organelletypelist: # vesicle 
            for name in filelist:
                if self.vesiclefilename in name:            self.vesiclefilepath = name
        if 1 in self.organelletypelist: # actin 
            for name in filelist:
                if self.actinfilename in name:              self.actinfilepath = name
                elif self.actinanglefilename in name:       self.actinanglefilepath = name
        if 2 in self.organelletypelist: # microtube 
            for name in filelist:
                if self.microtubefilename in name:          self.microtubefilepath = name
        if 3 in self.organelletypelist: # mito 
            for name in filelist:
                if self.mitofilename in name:               self.mitofilepath = name
        if 4 in self.organelletypelist: # endo_reticulum 
            for name in filelist:
                if self.endoreticulumfilename in name:      self.endoreticulumfilepath = name
        if 5 in self.organelletypelist: #  focal_adhesion
            for name in filelist:
                if self.focaladhesionfilename in name:      self.focaladhesionfilepath = name
        if 6 in self.organelletypelist: # calcium
            for name in filelist:
                if self.cafilename in name:                 self.cafilepath = name


        # check data
        if 0 in self.organelletypelist: # vesicle 
            if not os.path.exists(self.vesiclefilepath):        return False
        if 1 in self.organelletypelist: # actin 
            if not os.path.exists(self.actinfilepath):          return False
            if not os.path.exists(self.actinanglefilepath):     return False
        if 2 in self.organelletypelist: # microtube 
            if not os.path.exists(self.microtubefilepath):      return False
        if 3 in self.organelletypelist: # mito 
            if not os.path.exists(self.microtubefilepath):      return False
        if 4 in self.organelletypelist: # endo_reticulum 
            if not os.path.exists(self.endoreticulumfilepath):  return False
        if 5 in self.organelletypelist: # endo_reticulum 
            if not os.path.exists(self.focaladhesionfilepath):  return False
        if 6 in self.organelletypelist: # calcium
            if not os.path.exists(self.cafilepath):             return False

        print('organelletypelist', organelletypelist)
        print('filepath',self.vesiclefilepath, self.actinfilepath, self.microtubefilepath, self.mitofilepath, self.endoreticulumfilepath, self.focaladhesionfilepath)
        print(self.actinanglefilepath)
        print(self.cafilepath)
        return True


    def check_data_available(self, organelleidxlist):
        'check whether raw data exists'

        rawdata_dir = os.path.join(arg.root_dir, arg.data_root_dir)
        filelist = self.get_filelist(rawdata_dir,[])
        filelist = [ file for file in filelist if self.dataid in file.lower()] 
        datapath = '//'.join(filelist[0].split('\\')[:-1])
        filenamelist = preprocess.load_files(datapath)
        
        self.get_datanewpath()

        # print(self.condition)

        # print(145,filenamelist)

        if self.vesicle_type in organelleidxlist:
            if not os.path.exists(filenamelist[0]): sys.exit('Vesicle data lacked.')
            else: self.vesicle_rawdata_name = filenamelist[0]
        if self.actin_type in organelleidxlist:
            if not os.path.exists(filenamelist[2]): sys.exit('Actin data lacked.')
            else: self.actin_rawdata_name = filenamelist[2]
            if not os.path.exists(filenamelist[0]): sys.exit('Vesicle data lacked.')
            else: self.vesicle_rawdata_name = filenamelist[0]
        if self.microtube_type in organelleidxlist:
            if not os.path.exists(filenamelist[4]): sys.exit('Microtube data lacked.')
            else: self.microtube_rawdata_name = filenamelist[4]
            if not os.path.exists(filenamelist[0]): sys.exit('Vesicle data lacked.')
            else: self.vesicle_rawdata_name = filenamelist[0]
        if self.mito_type in organelleidxlist:
            if not os.path.exists(filenamelist[6]): sys.exit('Mtio data lacked.')
            else: self.mito_rawdata_name = filenamelist[6]
        if self.endo_reticulum_type in organelleidxlist:
            if not os.path.exists(filenamelist[7]): sys.exit('Endo reticulum data lacked.')
            else: self.endo_reticulum_rawdata_name = filenamelist[7]

        if self.focal_adhesion_type in organelleidxlist:
            self.check_fa_data_available()

        if self.ca_type in organelleidxlist:
            self.check_calcium_available()

    def check_fa_data_available(self):
        filelist = preprocess.load_focaladhesion_data(arg.fa_raw_data_dir, self.dataid)
        self.fa_raw_datapath = ''
        self.false_fa_raw_datapath = ''
        self.fa_raw_img_datapath = ''

        for name in filelist:
            if 'pl_list_all.xml' in name:
                self.fa_raw_datapath = name

            if 'false_FA.csv' in name:
                self.false_fa_raw_datapath = name
            elif 'FA_false.csv' in name:
                self.false_fa_raw_datapath = name

            if '_map.mrc' in name:
                self.fa_raw_img_datapath = name

        if not os.path.exists(self.fa_raw_datapath): sys.exit('focal adhesion data lacked.')
        if not os.path.exists(self.false_fa_raw_datapath):
            print('false_fa_raw_datapath not exists.')
        if not os.path.exists(self.fa_raw_img_datapath):
            sys.exit('fa_raw_img_datapath data lacked.')       

    def check_calcium_available(self):
        filelist = preprocess.load_ca_data(arg.ca_raw_data_dir, self.dataid)
        self.ca_raw_datapath = ''
        self.false_ca_raw_datapath = ''
        self.ca_raw_img_datapath = ''

        for name in filelist:
            if 'pl_list_all.xml' in name:
                self.ca_raw_datapath = name

            if 'false_Ca.csv' in name:
                self.false_ca_raw_datapath = name
            elif 'Ca_false.csv' in name:
                self.false_ca_raw_datapath = name
            if '_map.mrc' in name and 'ca' in name.lower():
                self.ca_raw_img_datapath = name

        print(self.ca_raw_datapath, self.false_ca_raw_datapath)
        if not os.path.exists(self.ca_raw_datapath): sys.exit('calcium data lacked.')
        if not os.path.exists(self.false_ca_raw_datapath):
            print('false_ca_raw_datapath not exists.')
        if not os.path.exists(self.ca_raw_img_datapath):
            sys.exit('ca_raw_img_datapath data lacked.')    


    def get_datanewpath(self):
        rawdata_dir = os.path.join(arg.root_dir, arg.data_root_dir)
        filelist = self.get_filelist(rawdata_dir,[])
        filelist = [ file for file in filelist if self.dataid in file.lower()] 
        print(filelist)
        datapath = '//'.join(filelist[0].split('\\')[:-1])
        filenamelist = preprocess.load_files(datapath)
        self.condition = filenamelist[0].split('//')[-2]
        self.location = filenamelist[0].split('//')[-3]
        self.datanewpath = f'{arg.root_dir}/{arg.processed_data_dir}/{self.location}/{self.condition}/{self.dataid}'


    def generate_processed_data(self):
        self.check_data_available(self.organelletypelist)

        self.datanewpath = f'{arg.root_dir}/{arg.processed_data_dir}/{self.location}/{self.condition}/{self.dataid}'
        if not os.path.exists(self.datanewpath): os.makedirs(self.datanewpath)
        
        if self.actin_type in self.organelletypelist:           self.preprocess_actin_data()
        if self.vesicle_type in self.organelletypelist:         self.preprocess_vesicle_data()
        if self.microtube_type in self.organelletypelist:       self.preprocess_microtube_data()
        if self.mito_type in self.organelletypelist:            self.preprocess_mito_data()
        if self.endo_reticulum_type in self.organelletypelist:  self.preprocess_endoreticulum_data()
        if self.focal_adhesion_type in self.organelletypelist:  self.preprocess_focaladhesion_data()
        if self.ca_type in self.organelletypelist:              self.preprocess_calcium_data()
        
        
        self.check_data_updated(self.organelletypelist)

    def preprocess_vesicle_data(self):
        rawdata_dir = os.path.join(arg.root_dir, arg.data_root_dir)
        vesicleoutputname = f'{self.datanewpath}/{self.dataid}_{self.vesiclefilename}'
        vesicle_data = mrcfile.open(self.vesicle_rawdata_name, permissive=True).data
        with mrcfile.new(vesicleoutputname, overwrite=True) as mrc:
            mrc.set_data(vesicle_data )

    def preprocess_actin_data(self):
        self.preprocess_vesicle_data()
        rawdata_dir = os.path.join(arg.root_dir, arg.data_root_dir, self.condition, self.dataid)
        # print(rawdata_dir)
        # filelist = self.get_filelist(rawdata_dir,[])
        # filelist = [ file for file in filelist if self.dataid in file.lower()] 
        actindict = preprocess.fill_actin_point_gap(self.actin_rawdata_name, )
        actinoutputpath = f'{self.datanewpath}/{self.dataid}_{self.actinfilename}'
        

        with open(actinoutputpath, 'w') as f:
            json.dump(actindict, f)
        print('saved', actinoutputpath)

        self.preprocess_actin_angle_data()

    def preprocess_actin_angle_data(self):

        # angles are coming from actin rawdata file

        actinangledict = preprocess.get_actin_to_pm_angle(self.actin_rawdata_name, )

        actinangleoutputpath = f'{self.datanewpath}/{self.dataid}_{self.actinanglefilename}'

        with open(actinangleoutputpath, 'w') as f:
            json.dump(actinangledict, f)
        print('saved', actinangleoutputpath)

    def preprocess_microtube_data(self):
        self.preprocess_vesicle_data()
        microtubedict = preprocess.fill_microtube_point_gap(self.microtube_rawdata_name,)
        microtubeoutputpath = f'{self.datanewpath}/{self.dataid}_{self.microtubefilename}'
        with open(microtubeoutputpath, 'w') as f:
            json.dump(microtubedict, f)
        print('saved', microtubeoutputpath)

    def preprocess_mito_data(self):
        mitooutputname = f'{self.datanewpath}/{self.dataid}_{self.mitofilename}'
        mito_data = mrcfile.open(self.mito_rawdata_name, permissive=True).data
        with mrcfile.new(mitooutputname, overwrite=True) as mrc:
            mrc.set_data(mito_data )

    def preprocess_endoreticulum_data(self):
        endoreticulumoutputname = f'{self.datanewpath}/{self.dataid}_{self.endoreticulumfilename}'
        endoreticulum_data = mrcfile.open(self.endo_reticulum_rawdata_name, permissive=True).data
        with mrcfile.new(endoreticulumoutputname, overwrite=True) as mrc:
            mrc.set_data(endoreticulum_data )        

    def preprocess_focaladhesion_data(self):

        facoords, fa_img = preprocess.filter_fa_point(self.fa_raw_datapath, self.false_fa_raw_datapath, self.fa_raw_img_datapath)

        confirmed_facoord = np.array(facoords)
        coordpd = pd.DataFrame({'X': confirmed_facoord[:,0], 'Y':confirmed_facoord[:,1] ,'Z':confirmed_facoord[:,2]})
        Fa_data = coordpd

        focaladhesion_outputpath = f'{self.datanewpath}/{self.dataid}_{self.focaladhesionfilename}'
        Fa_data.to_csv(focaladhesion_outputpath, index=False)

        focaladhesion_img_outputpath = f'{self.datanewpath}/{self.dataid}_{self.focaladhesion_imagefilename}'
        with mrcfile.new(focaladhesion_img_outputpath, overwrite = True) as mrc:
            mrc.set_data(fa_img) 

    def preprocess_calcium_data(self):
        cacoords, ca_img = preprocess.filter_fa_point(self.ca_raw_datapath, self.false_ca_raw_datapath, self.ca_raw_img_datapath)
        confirmed_cacoord = np.array(cacoords)
        coordpd = pd.DataFrame({'X': confirmed_cacoord[:,0], 'Y':confirmed_cacoord[:,1] ,'Z':confirmed_cacoord[:,2]})
        Ca_data = coordpd

        calcium_outputpath = f'{self.datanewpath}/{self.dataid}_{self.cafilename}'
        Ca_data.to_csv(calcium_outputpath, index=False)

        Ca_img_outputpath = f'{self.datanewpath}/{self.dataid}_{self.ca_imagefilename}'
        # with mrcfile.new(Ca_img_outputpath, overwrite = True) as mrc:
        #     mrc.set_data(ca_img) 


        print(calcium_outputpath, Ca_img_outputpath)


    def load_plotdataname(self):

        self.singledataoutplotpath = f'{arg.root_dir}/{arg.output_dir}/singledata/'
        self.conditiondataoutplotpath = f'{arg.root_dir}/{arg.output_dir}/conditions/'
        self.singledatavis_outplotpath = f'{arg.root_dir}/{arg.output_dir}/singlevisdata/'

        self.actin_surround_actin_vect_map_filename = f'actin_surround_actin_vect_map.csv'
        self.actin_surround_actin_plotname = f'actin_surrounding_actin_probability_density_ditribution.png'

        self.actin_surround_actin_vect_dist_map_filename = f'actin_surround_actin_vect_dist_map.csv'
        self.actin_surround_actin_dist_plotname = f'actin_surrounding_actin_dist_probability_density_ditribution.png'


        self.actin_to_vesicle_filename = f'actin_to_vesicle_distances.csv' 
        self.actin_to_vesicle_plotname = f'actin_distance_to_vesicle_distance_distribution.png' 

        self.actin_vesicle_vis_filename = f'space_vis_actin_no_vesicle.mrc'
        self.microtube_vesicle_vis_filename = f'space_vis_microtube_no_vesicle.mrc'

        self.microtube_surround_actin_vect_map_filename = f'microtube_surround_actin_vect_map.csv'
        self.microtube_surround_actin_plotname = f'microtube_surrounding_actin_probability_density_ditribution.png'

        self.microtube_surround_actin_vect_dist_map_filename = f'microtube_surround_actin_vect_dist_map.csv'
        self.microtube_surround_actin_dist_plotname = f'microtube_surrounding_actin_dist_probability_density_ditribution.png'


        self.actin_to_mito_filename = f'actin_to_mito_distances.csv'
        self.actin_to_mito_plotname = f'actin_to_mito_distance_distribution.png'

        self.actin_to_endoreticulum_filename = f'actin_to_endoreticulum_distances.csv'
        self.actin_to_endoreticulum_plotname = f'actin_to_endoreticulum_distance_distribution.png'

        self.microtube_to_vesicle_filename = f'microtube_to_vesicle_distances.csv'
        self.microtube_to_vesicle_plotname = f'microtube_to_vesicle_distance_distribution.png'

        self.microtube_to_mito_filename = f'microtube_to_mito_distances.csv'
        self.microtube_to_mito_plotname = f'microtube_to_mito_distance_distribution.png'

        self.microtube_to_endoreticulum_filename = f'microtube_to_endoreticulum_distances.csv'
        self.microtube_to_endoreticulum_plotname = f'microtube_to_endoreticulum_distance_distribution.png'

        self.vesicle_to_mito_filename = f'vesicle_to_mito_distances.csv'
        self.vesicle_to_mito_plotname = f'vesicle_to_mito_distance_distribution.png'

        self.vesicle_to_endoreticulum_filename = f'vesicle_to_endoreticulum_distances.csv'
        self.vesicle_to_endoreticulum_plotname = f'vesicle_to_endoreticulum_distances_distribution.png'

        self.mito_to_endoreticulum_filename = f'mito_to_endoreticulum_distances.csv'
        self.mito_to_endoreticulum_plotname = f'mito_to_endoreticulum_distances_distribution.png'

        self.actin_angle_distance_pair_filename = f'actin_angle_distance_pair.csv'
        self.actin_angle_distance_pair_plotname = f'actin_angle_distance_pair_distribution.png'

        self.focal_to_focal_distance_filename = f'focal_to_focal_distance.csv'
        self.focal_to_focal_distance_plotname = f'focal_to_focal_distance_distribution.png'

        self.focal_aggregated_actin_angle_map_filename = f'focal_aggregated_actin_angle_map.csv'
        self.focal_aggregated_actin_angle_map_plotname = f'focal_aggregated_actin_angle_map_distribution.png'

        self.focal_aggregated_actin_angles_pair_filename = f'focal_aggregated_actin_angles_pair.csv'
        self.focal_aggregated_actin_angles_pair_plotname = f'focal_aggregated_actin_angles_pair_distribution.png'
        
        self.focal_to_ca_distance_filename = f'focal_to_ca_distance.csv'
        self.focal_to_ca_distance_plotname = f'focal_to_ca_distance_distribution.png'




class check_data(process_organelles):

    def __init__(self, dataid) -> None:
        super().__init__(dataid)
    
    def check_organelle(self):
        # ISG	Actin	MT	Mito ER
        return self.check_type(self.actin_type, self.vesicle_type)


    def step(self):
        # ISG	Actin	MT	Mito ER 
        print('cur dataid', self.dataid)

        if self.check_type(self.vesicle_type, self.actin_type):
            organelletypelst = [self.actin_type, self.vesicle_type]
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()

            self.match_vesicle_actin()


    def match_vesicle_actin(self):
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)  
        with open(f'{self.actinfilepath}', 'r') as f:
            actin_data = json.load(f)

        # print(111,len(actin_data))
        if len(actin_data) == 0:
            print('no filament')
        else:
            vesicle_mem_mask = mrcfile.open(self.vesiclefilepath, permissive=True).data
            # print(self.dataid, 'vesicle_mem_mask shape:', vesicle_mem_mask.shape)

            if np.sum(vesicle_mem_mask) == 0:
                print('no vesicle')
            else:
                actin_data_new = distances_generator.shift_bias(actin_data, shift_xml_xyz)

                # dist = distances_generator.coords_to_mem_distance_generator(actin_data_new, vesicle_mem_mask, voxel_size_xyz)

                filament_coords, mask, voxelsize = actin_data_new, vesicle_mem_mask, voxel_size_xyz

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

                voxel_coords = coords

                coords_x = [ voxel_coords[i][0] for i in range(len(voxel_coords))]
                coords_y = [ voxel_coords[i][1] for i in range(len(voxel_coords))]
                coords_z = [ voxel_coords[i][2] for i in range(len(voxel_coords))]



                if not (np.max(coords_x) <= mask.shape[0] and np.min(coords_x)>= 0) or not (np.max(coords_y) <= mask.shape[1] and np.min(coords_y)>= 0) or not (np.max(coords_z) <= mask.shape[2] and np.min(coords_z)>= 0):
                    print('x',np.max(coords_x),np.min(coords_x), mask.shape[0])
                    print('y',np.max(coords_y),np.min(coords_y), mask.shape[1])
                    print('z',np.max(coords_z),np.min(coords_z), mask.shape[2])

                assert np.max(coords_x) <= mask.shape[0] and np.min(coords_x)>= 0
                assert np.max(coords_y) <= mask.shape[1] and np.min(coords_y)>= 0
                assert np.max(coords_z) <= mask.shape[2] and np.min(coords_z)>= 0



class actin_to_actin(process_organelles):
    def __init__(self,dataid):
        super(actin_to_actin, self).__init__(dataid)


    def check_actintype(self):
        # ISG	Actin	MT	Mito ER

        return self.check_type(self.actin_type, self.actin_type)

    def step(self):

        if not self.check_actintype():
            print(f'{self.dataid} actin_to_actin data not exist.')
        else:
            organelletypelst = [self.actin_type]
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()
            self.actin_surround_actin_map_plot(overlap = self.plotdata_overlap)



    def actin_surround_actin_map_plot(self, overlap):
        self.get_datanewpath()
        self.actin_surround_actin_vect_map_filepath = f'{self.datanewpath}/{self.dataid}_{self.actin_surround_actin_vect_map_filename}'
        self.actin_surround_actin_vect_dist_map_filepath = f'{self.datanewpath}/{self.dataid}_{self.actin_surround_actin_vect_dist_map_filename}'
        if not os.path.exists(self.actin_surround_actin_vect_map_filepath) or not os.path.exists(self.actin_surround_actin_vect_dist_map_filepath) or overlap == True:
            self.actin_surround_actin_map()


        with open(f'{self.actin_surround_actin_vect_map_filepath}', 'r', encoding='utf-8') as f:
            filament_vect_map_lst = json.load(f)
        print(self.dataid, 'filament num,', len(filament_vect_map_lst))  ## list[x,y,z]

        filament_vect_map_lst_new = []
        vect_count_n = 0
        for filament in filament_vect_map_lst:
            curcoord = [ np.array(vect) for vect in filament ]
            filament_vect_map_lst_new.append(curcoord)
            vect_count_n += len(curcoord)
        filament_vect_map_lst = filament_vect_map_lst_new



        outplot.vis_actin_surround_actin(filament_vect_map_lst, self.dataid, len(filament_vect_map_lst), vect_count_n) 
        # plt.clim()
         
        plt.savefig(f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.actin_surround_actin_plotname}')
        # plt.savefig(f'{arg.root_dir}/{arg.output_dir}/singledata/{location}_{condition}_{dataid}_actin_surrounding_probability_density_ditribution.pdf')
        # plt.show()
        plt.close()
        print(f'saved {self.dataid} {self.actin_surround_actin_plotname}')


        # ----
        filament_vect_dist_map_lst = pd.read_csv(self.actin_surround_actin_vect_dist_map_filepath)

        print(self.dataid, 'filament num,', len(filament_vect_map_lst))  ## list[x,y,z]

        outplot.vis_actin_surround_actin_dist(filament_vect_dist_map_lst,self.dataid, len(filament_vect_map_lst), vect_count_n)
        plt.savefig(f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.actin_surround_actin_dist_plotname}')
        plt.show()
        plt.close()



    def actin_surround_actin_map(self):
        '''
        1. shift bias
        2. generate map
        '''
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)

        actinoutputpath = f'{self.datanewpath}/{self.dataid}_{self.actinfilename}'
        with open(f'{actinoutputpath}', 'r') as f:
            actin_data = json.load(f)

        actin_data_new = distances_generator.shift_bias(actin_data, shift_xml_xyz)

        vesicle_mem_mask = mrcfile.open(f'{self.datanewpath}/{self.dataid}_{self.vesiclefilename}', permissive=True).data

        imgsize = vesicle_mem_mask.shape
        print(self.dataid, imgsize)

        actin_vect_map_lst = distances_generator.actin_to_actin_distance(actin_data_new, imgsize, voxel_size_xyz)
        actin_vect_map_filename = self.actin_surround_actin_vect_map_filepath
        actin_vect_map_lst_save = []


        for actin in actin_vect_map_lst:
            cur_vect = [ [vect[0], vect[1], vect[2]]  for vect in actin ]
            actin_vect_map_lst_save.append(cur_vect)


        with open(self.actin_surround_actin_vect_map_filepath, 'w', encoding='utf-8') as f:
                json.dump(actin_vect_map_lst_save, f)

        self.actin_surround_actin_dist_map()


    def actin_surround_actin_dist_map(self):
        
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)

        actinoutputpath = f'{self.datanewpath}/{self.dataid}_{self.actinfilename}'
        with open(f'{actinoutputpath}', 'r') as f:
            actin_data = json.load(f)

        actin_data_new = distances_generator.shift_bias(actin_data, shift_xml_xyz)

        vesicle_mem_mask = mrcfile.open(f'{self.datanewpath}/{self.dataid}_{self.vesiclefilename}', permissive=True).data

        imgsize = vesicle_mem_mask.shape
        print(self.dataid, imgsize)

        dist = distances_generator.actin_to_actin_distance2(actin_data_new, imgsize, voxel_size_xyz)

        
        dist_savepd = pd.DataFrame(columns=['Distance'], data= dist)
        dist_savepd.to_csv(f'{self.datanewpath}/{self.dataid}_{self.actin_surround_actin_vect_dist_map_filename}', index=False)




class actin_to_vesicle(process_organelles):
    def __init__(self, dataid):
        super().__init__(dataid)
    
    def check_actin_vesicle_type(self):
        # ISG	Actin	MT	Mito ER
        return self.check_type(self.actin_type, self.vesicle_type)

    def step(self):

        if not self.check_actin_vesicle_type():
            print(f'{self.dataid} actin_to_vesicle data not exist.')

        else:
            organelletypelst = [self.actin_type, self.vesicle_type]
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()
            self.actin_to_vesicle_plot(overlap = self.plotdata_overlap)



    def actin_to_vesicle_plot(self, overlap):
        self.get_datanewpath()
        self.actin_to_vesicle_filepath = f'{self.datanewpath}/{self.dataid}_{self.actin_to_vesicle_filename}'
        if not os.path.exists(self.actin_to_vesicle_filepath) or overlap == True:
            self.actin_to_vesicle_distance()


        with open(f'{self.actinfilepath}', 'r') as f:
            actin_data = json.load(f)

        dist = pd.read_csv(self.actin_to_vesicle_filepath)
        outplot.actin_to_vesicle_dist_hist(dist, self.dataid, len(actin_data), len(dist))
        plt.savefig(f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.actin_to_vesicle_plotname}')
        # plt.show()
        plt.close()
        


    def actin_to_vesicle_distance(self):
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)
        
        with open(f'{self.actinfilepath}', 'r') as f:
            actin_data = json.load(f)

        # print(self.dataid, 'actin num,', len(actin_data))

        vesicle_mem_mask = mrcfile.open(self.vesiclefilepath, permissive=True).data
        # print(self.dataid, 'vesicle_mem_mask shape:', vesicle_mem_mask.shape)


        actin_data_new = distances_generator.shift_bias(actin_data, shift_xml_xyz)

        dist = distances_generator.coords_to_mem_distance_generator(actin_data_new, vesicle_mem_mask, voxel_size_xyz)

        dist_savepd = pd.DataFrame(columns=['Distance'], data= dist)
        # filament2filament_dist_lst_savepd.loc['distance'] = filament2filament_dist_lst
        dist_savepd.to_csv(f'{self.datanewpath}/{self.dataid}_{self.actin_to_vesicle_filename}', index=False)


class actin_to_vesicle_dist_angle_pair(process_organelles):
    def __init__(self, dataid):
        super().__init__(dataid)
    
    def check_actin_vesicle_type(self):
        # ISG	Actin	MT	Mito ER
        return self.check_type(self.actin_type, self.vesicle_type)



    def step(self):

        if not self.check_actin_vesicle_type():
            print(f'{self.dataid} actin_to_vesicle data not exist.')

        else:
            organelletypelst = [self.actin_type, self.vesicle_type]
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()
            self.actin_to_vesicle_dist_angle_pair_plot(overlap = self.plotdata_overlap)

        


    def actin_to_vesicle_dist_angle_pair_plot(self, overlap):
        self.get_datanewpath()
        self.actin_to_vesicle_filepath = f'{self.datanewpath}/{self.dataid}_{self.actin_to_vesicle_filename}'
        if not os.path.exists(self.actin_to_vesicle_filepath) or overlap == True:
            self.actin_to_vesicle_distance()


        with open(f'{self.actinfilepath}', 'r') as f:
            actin_data = json.load(f)


        dist = pd.read_csv(self.actin_to_vesicle_filepath)
        outplot.actin_to_vesicle_dist_angle_pair_plot(dist, actin_angle_data,  self.dataid, len(actin_data), len(dist))
        
        # plt.savefig(f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.actin_to_vesicle_plotname}')
        # # plt.show()
        # plt.close()
        


    def actin_to_vesicle_distance(self):
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)
        
        with open(f'{self.actinfilepath}', 'r') as f:
            actin_data = json.load(f)

        # print(self.dataid, 'actin num,', len(actin_data))

        vesicle_mem_mask = mrcfile.open(self.vesiclefilepath, permissive=True).data
        # print(self.dataid, 'vesicle_mem_mask shape:', vesicle_mem_mask.shape)


        actin_data_new = distances_generator.shift_bias(actin_data, shift_xml_xyz)

        dist = distances_generator.coords_to_mem_distance_generator(actin_data_new, vesicle_mem_mask, voxel_size_xyz)

        dist_savepd = pd.DataFrame(columns=['Distance'], data= dist)
        # filament2filament_dist_lst_savepd.loc['distance'] = filament2filament_dist_lst
        dist_savepd.to_csv(f'{self.datanewpath}/{self.dataid}_{self.actin_to_vesicle_filename}', index=False)

        # ----

        with open(self.actinanglefilepath, 'r') as f:
            actin_angle_data = json.load(f)    

    
        anglelst1, anglelst2 = distances_generator.Fa_aggragated_actinactin_angle(actin_data, actin_angle_data, fa_data, voxel_size_xyz, self.dist_threshold)
        angle_savepd = pd.DataFrame({'Angle actin':anglelst1, 'Angle mem':anglelst2})

        # filament2filament_dist_lst_savepd.loc['distance'] = filament2filament_dist_lst
        angle_savepd.to_csv(f'{self.datanewpath}/{self.dataid}_{self.focal_aggregated_actin_angles_pair_filename}', index=False)






class actin_to_vesicle_vis(process_organelles):
    def __init__(self, dataid):
        super().__init__(dataid)
    
    def check_actin_vesicle_type(self):
        # ISG	Actin	MT	Mito ER
        return self.check_type(self.actin_type, self.vesicle_type)

    def step(self):

        if not self.check_actin_vesicle_type():
            print(f'{self.dataid} actin_to_vesicle_vis data not exist.')

        else:
            organelletypelst = [self.actin_type, self.vesicle_type]
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()

            self.actin_to_vesicle_visimg(overlap = self.plotdata_overlap)

    def actin_to_vesicle_visimg(self, overlap):
        self.get_datanewpath()
        self.actin_to_vesicle_filepath = f'{self.datanewpath}/{self.dataid}_{self.actin_to_vesicle_filename}'

        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)
        
        with open(f'{self.actinfilepath}', 'r') as f:
            actin_data = json.load(f)

        if len(actin_data) == 0:
            print('no actin')
            vesicle_mem_mask = mrcfile.open(self.vesiclefilepath, permissive=True).data
            vis_img = np.zeros_like(vesicle_mem_mask)
        else:
            vesicle_mem_mask = mrcfile.open(self.vesiclefilepath, permissive=True).data

            actin_data_new = distances_generator.shift_bias(actin_data, shift_xml_xyz)

            # dist = distances_generator.coords_to_mem_distance_generator(actin_data_new, vesicle_mem_mask, voxel_size_xyz)

            vis_img = distances_generator.filadistance2mem_vis(actin_data_new, vesicle_mem_mask, voxel_size_xyz)

        vis_img3 = vis_img.astype(np.int16)
        outputname = f'{self.singledatavis_outplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.actin_vesicle_vis_filename}'
        with mrcfile.new(outputname, overwrite=True) as mrc:
            mrc.set_data(vis_img3)







class microtube_to_actin(process_organelles):
    def __init__(self, dataid):
        super().__init__(dataid)
    
    def check_microtube_actin_type(self):
        # ISG	Actin	MT	Mito ER
        return self.check_type(self.microtube_type, self.actin_type)

    def step(self):

        if not self.check_microtube_actin_type():
            print(f'{self.dataid} microtube_to_actin data not exist.')

        else:
            organelletypelst = [self.microtube_type, self.actin_type]
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()
            self.microtube_to_actin_map_plot(overlap = self.plotdata_overlap)


    def microtube_to_actin_map_plot(self, overlap):
        self.get_datanewpath()
        self.microtube_surround_actin_vect_map_filepath= f'{self.datanewpath}/{self.dataid}_{self.microtube_surround_actin_vect_map_filename}'
        self.microtube_surround_actin_vect_dist_map_filepath = f'{self.datanewpath}/{self.dataid}_{self.microtube_surround_actin_vect_dist_map_filename}'
        if not os.path.exists(self.microtube_surround_actin_vect_map_filepath) or not os.path.exists(self.microtube_surround_actin_vect_dist_map_filepath) or overlap == True:
            self.microtube_to_actin_map()


        with open(self.microtube_surround_actin_vect_map_filepath, 'r', encoding='utf-8') as f:
            filament_vect_2_MT_map_lst = json.load(f)
        print(self.dataid, 'MT num,', len(filament_vect_2_MT_map_lst))  ## list[x,y,z]

        filament_vect_2_MT_map_lst_new = []
        vect_count_n = 0
        for filament in filament_vect_2_MT_map_lst:
            curcoord = [ np.array(vect) for vect in filament ]
            filament_vect_2_MT_map_lst_new.append(curcoord)
            vect_count_n += len(curcoord)

        filament_vect_2_MT_map_lst = filament_vect_2_MT_map_lst_new

        outplot.vis_microtube_surround_actin(filament_vect_2_MT_map_lst, self.dataid, len(filament_vect_2_MT_map_lst), vect_count_n) 
        plt.savefig(f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.microtube_surround_actin_plotname}')
        # plt.show()
        plt.close()

        # ----
        filament_vect_2_MT_dist_map_lst = pd.read_csv(self.microtube_surround_actin_vect_dist_map_filepath)

        print(self.dataid, 'filament 2 MT num,', len(filament_vect_2_MT_dist_map_lst))  ## list[x,y,z]

        outplot.vis_microtube_surround_actin_dist(filament_vect_2_MT_dist_map_lst, self.dataid, len(filament_vect_2_MT_dist_map_lst), vect_count_n)
        plt.savefig(f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.microtube_surround_actin_dist_plotname}')
        plt.show()
        plt.close()




    def microtube_to_actin_map(self):
        # shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(arg, self.dataid)

        with open(f'{self.actinfilepath}', 'r') as f:
            actin_data = json.load(f)

        with open(f'{self.microtubefilepath}', 'r') as f:
            MT_data = json.load(f)


        filament_vect_2_MT_map_lst = distances_generator.actin_to_microtube_distance(actin_data, MT_data)
        filament_vect_2_MT_map_lst_save =[]
        
        for filament in filament_vect_2_MT_map_lst:
            cur_vect = [ [vect[0], vect[1], vect[2]]  for vect in filament ]
            filament_vect_2_MT_map_lst_save.append(cur_vect)

        with open(f'{self.datanewpath}/{self.dataid}_{self.microtube_surround_actin_vect_map_filename}', 'w', encoding='utf-8') as f:
                json.dump(filament_vect_2_MT_map_lst_save, f)

        self.microtube_to_actin_dist_map()


    def microtube_to_actin_dist_map(self):
        # shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(arg, self.dataid)

        with open(f'{self.actinfilepath}', 'r') as f:
            actin_data = json.load(f)

        with open(f'{self.microtubefilepath}', 'r') as f:
            MT_data = json.load(f)


        filament_vect_2_MT_map_lst = distances_generator.actin_to_microtube_distance(actin_data, MT_data)
        filament_vect_2_MT_map_lst_save =[]
        
        for filament in filament_vect_2_MT_map_lst:
            cur_vect = [ [vect[0], vect[1], vect[2]]  for vect in filament ]
            filament_vect_2_MT_map_lst_save.append(cur_vect)

        with open(f'{self.datanewpath}/{self.dataid}_{self.microtube_surround_actin_vect_map_filename}', 'w', encoding='utf-8') as f:
                json.dump(filament_vect_2_MT_map_lst_save, f)



        dist = distances_generator.actin_to_microtube_distance2(actin_data, MT_data)

        dist_savepd = pd.DataFrame(columns=['Distance'], data= dist)
        dist_savepd.to_csv(f'{self.datanewpath}/{self.dataid}_{self.microtube_surround_actin_vect_dist_map_filename}', index=False)





class actin_to_mito(process_organelles):
    def __init__(self, dataid):
        super().__init__(dataid)
    
    def check_organelle_type(self,organelletype):
        # ISG	Actin	MT	Mito ER
        return self.check_type(organelletype[0], organelletype[1])

    def step(self):
        organelletypelst = [self.mito_type, self.actin_type]
        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} actin_to_mito data not exist.')

        else:
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()
            self.actin_to_mito_plot(overlap = self.plotdata_overlap)
        print('Done with actin_to_mito.')


    def actin_to_mito_plot(self, overlap):
        self.actin_to_mito_filepath = f'{self.datanewpath}/{self.dataid}_{self.actin_to_mito_filename}'
        if not os.path.exists(self.actin_to_mito_filepath) or overlap == True:
            self.actin_to_mito_distance()

        with open(f'{self.actinfilepath}', 'r') as f:
            actin_data = json.load(f)
        dist2 = pd.read_csv(self.actin_to_mito_filepath)
        outplot.actindist2mito_hist(dist2, self.dataid, len(actin_data), len(dist2))
        plotname = f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.actin_to_mito_plotname}'
        plt.savefig(plotname)
        # plt.show()
        plt.close()



    def actin_to_mito_distance(self):

        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)
        with open(f'{self.actinfilepath}', 'r') as f:
            filament_data = json.load(f)

        mito_mem_mask = mrcfile.open(self.mitofilepath, permissive=True).data
        filament_data_new = distances_generator.shift_bias(filament_data, shift_xml_xyz)

        dist = distances_generator.coords_to_mem_distance_generator(filament_data_new, mito_mem_mask, voxel_size_xyz)
        dist_savepd = pd.DataFrame(columns=['Distance'], data= dist)

        dist_savepd.to_csv(self.actin_to_mito_filepath, index=False)


class actin_to_endoreticulum(process_organelles):
    def __init__(self, dataid):
        super().__init__(dataid)
    
    def check_organelle_type(self,organelletype):
        # ISG	Actin	MT	Mito ER
        return self.check_type(organelletype[0], organelletype[1])

    def step(self):
        organelletypelst = [self.endo_reticulum_type, self.actin_type]
        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} actin_to_endoreticulum data not exist.')

        else:
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()
            self.actin_to_endoreticulum_plot(overlap = self.plotdata_overlap)
        print('Done with actin_to_endoreticulum.')


    def actin_to_endoreticulum_plot(self, overlap):

        self.actin_to_endoreticulum_filepath = f'{self.datanewpath}/{self.dataid}_{self.actin_to_endoreticulum_filename}'
        if not os.path.exists(self.actin_to_endoreticulum_filepath) or overlap == True:
            self.actin_to_endoreticulum_distance()

        with open(f'{self.actinfilepath}', 'r') as f:
            actin_data = json.load(f)
        dist2 = pd.read_csv(self.actin_to_endoreticulum_filepath)
        outplot.actindist2er_hist(dist2, self.dataid, len(actin_data), len(dist2))
        plotname = f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.actin_to_endoreticulum_plotname}'
        plt.savefig(plotname)
        # plt.show()
        plt.close()


    def actin_to_endoreticulum_distance(self):

        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)
        with open(f'{self.actinfilepath}', 'r') as f:
            filament_data = json.load(f)

        er_mem_mask = mrcfile.open(self.endoreticulumfilepath, permissive=True).data
        filament_data_new = distances_generator.shift_bias(filament_data, shift_xml_xyz)

        dist = distances_generator.coords_to_mem_distance_generator(filament_data_new, er_mem_mask, voxel_size_xyz)
        dist_savepd = pd.DataFrame(columns=['Distance'], data= dist)

        dist_savepd.to_csv(self.actin_to_endoreticulum_filepath, index=False)


class microtube_to_vesicle(process_organelles):
    def __init__(self, dataid):
        super().__init__(dataid)
    

    def step(self):
        organelletypelst = [self.microtube_type, self.vesicle_type]
        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} microtube_to_vesicle data not exist.')

        else:
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()

            self.microtube_to_vesicle_plot(overlap = self.plotdata_overlap)
        print('Done with microtube_to_vesicle.')


    def microtube_to_vesicle_plot(self, overlap):
        self.microtube_to_vesicle_filepath = f'{self.datanewpath}/{self.dataid}_{self.microtube_to_vesicle_filename}'
        if not os.path.exists(self.microtube_to_vesicle_filepath) or overlap == True:
            self.microtube_to_vesicle_distance()

        with open(f'{self.microtubefilepath}', 'r') as f:
            MT_data = json.load(f)

        dist2 = pd.read_csv(self.microtube_to_vesicle_filepath)
        outplot.dist_mt2isghist(dist2, self.dataid, len(MT_data), len(dist2))
        plotname = f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.microtube_to_vesicle_plotname}'
        plt.savefig(plotname)
        # plt.show()
        plt.close()


    def microtube_to_vesicle_distance(self):
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)

        with open(f'{self.microtubefilepath}', 'r') as f:
            MT_data = json.load(f)

        vesicle_mem_mask = mrcfile.open(self.vesiclefilepath, permissive=True).data
        MT_data_new = distances_generator.shift_bias(MT_data, shift_xml_xyz)

        dist = distances_generator.coords_to_mem_distance_generator(MT_data_new, vesicle_mem_mask, voxel_size_xyz)
        dist_savepd = pd.DataFrame(columns=['Distance'], data= dist)
        dist_savepd.to_csv(self.microtube_to_vesicle_filepath, index=False)


class microtube_to_vesicle_vis(process_organelles):
    def __init__(self, dataid):
        super().__init__(dataid)
    
    def check_actin_vesicle_type(self):
        # ISG	Actin	MT	Mito ER
        return self.check_type(self.microtube_type, self.vesicle_type)

    def step(self):

        organelletypelst = [self.microtube_type, self.vesicle_type]
        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} microtube_vesicle_vis data not exist.')

        else:
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()

            self.microtube_to_vesicle_visimg(overlap = self.plotdata_overlap)


    def microtube_to_vesicle_visimg(self, overlap):
        self.get_datanewpath()
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)

        with open(f'{self.microtubefilepath}', 'r') as f:
            MT_data = json.load(f)

        if len(MT_data) == 0:
            print('no microtube data')
            vesicle_mem_mask = mrcfile.open(self.vesiclefilepath, permissive=True).data
            vis_img = np.zeros_like(vesicle_mem_mask)

        else:
            vesicle_mem_mask = mrcfile.open(self.vesiclefilepath, permissive=True).data
            MT_data_new = distances_generator.shift_bias(MT_data, shift_xml_xyz)



        # dist = distances_generator.coords_to_mem_distance_generator(actin_data_new, vesicle_mem_mask, voxel_size_xyz)

            vis_img = distances_generator.filadistance2mem_vis(MT_data_new, vesicle_mem_mask, voxel_size_xyz)

        vis_img3 = vis_img.astype(np.int16)
        outputname = f'{self.singledatavis_outplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.microtube_vesicle_vis_filename}'
        with mrcfile.new(outputname, overwrite=True) as mrc:
            mrc.set_data(vis_img3)



class microtube_to_mito(process_organelles):
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'microtube_to_vesicle'

    def step(self):
        organelletypelst = [self.microtube_type, self.mito_type]
        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} {self.name} data not exist.')

        else:
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()

            self.microtube_to_mito_plot(overlap = self.plotdata_overlap)
        print(f'Done with {self.name}.')

    def microtube_to_mito_plot(self, overlap):
        
        self.microtube_to_mito_filepath = f'{self.datanewpath}/{self.dataid}_{self.microtube_to_mito_filename}'
        self.curdatapath = self.microtube_to_mito_filepath
        if not os.path.exists(self.curdatapath) or overlap == True:
            self.microtube_to_mito_distance()

        with open(f'{self.microtubefilepath}', 'r') as f:
            MT_data = json.load(f)

        dist2 = pd.read_csv(self.curdatapath)
        outplot.dist_mt2mitohist(dist2, self.dataid, len(MT_data), len(dist2))

        plotname = f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.microtube_to_mito_plotname}'
        plt.savefig(plotname)
        # plt.show()
        plt.close()

    def microtube_to_mito_distance(self):
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)

        with open(f'{self.microtubefilepath}', 'r') as f:
            MT_data = json.load(f)
        mito_mem_mask = mrcfile.open(self.mitofilepath, permissive=True).data
        MT_data_new = distances_generator.shift_bias(MT_data, shift_xml_xyz)

        dist = distances_generator.coords_to_mem_distance_generator(MT_data_new, mito_mem_mask, voxel_size_xyz)
        dist_savepd = pd.DataFrame(columns=['Distance'], data= dist)
        dist_savepd.to_csv(self.curdatapath, index=False)


class microtube_to_endoreticulum(process_organelles):
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'microtube_to_endoreticulum'

    def step(self):
        organelletypelst = [self.microtube_type, self.endo_reticulum_type]
        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} {self.name} data not exist.')

        else:
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()

            self.microtube_to_endoreticulum_plot(overlap = self.plotdata_overlap)
        print(f'Done with {self.name}.')


    def microtube_to_endoreticulum_plot(self, overlap):
        
        self.microtube_to_endoreticulum_filepath = f'{self.datanewpath}/{self.dataid}_{self.microtube_to_endoreticulum_filename}'
        self.curdatapath = self.microtube_to_endoreticulum_filepath
        if not os.path.exists(self.curdatapath) or overlap == True:
            self.microtube_to_endoreticulum_distance()

        with open(f'{self.microtubefilepath}', 'r') as f:
            MT_data = json.load(f)

        dist2 = pd.read_csv(self.curdatapath)
        outplot.dist_mt2erhist(dist2, self.dataid, len(MT_data), len(dist2))

        plotname = f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.microtube_to_endoreticulum_plotname}'
        plt.savefig(plotname)
        # plt.show()
        plt.close()

    def microtube_to_endoreticulum_distance(self):
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)

        with open(f'{self.microtubefilepath}', 'r') as f:
            MT_data = json.load(f)
        er_mem_mask = mrcfile.open(self.endoreticulumfilepath, permissive=True).data
        MT_data_new = distances_generator.shift_bias(MT_data, shift_xml_xyz)

        dist = distances_generator.coords_to_mem_distance_generator(MT_data_new, er_mem_mask, voxel_size_xyz)
        dist_savepd = pd.DataFrame(columns=['Distance'], data= dist)
        dist_savepd.to_csv(self.curdatapath, index=False)


class vesicle_to_mito(process_organelles):
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'vesicle_to_mito'

    def step(self):
        organelletypelst = [self.vesicle_type, self.mito_type]
        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} {self.name} data not exist.')

        else:
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()

            self.vesicle_to_mito_plot(overlap = self.plotdata_overlap)
        print(f'Done with {self.name}.')


    def vesicle_to_mito_plot(self, overlap):

        self.vesicle_to_mito_filepath = f'{self.datanewpath}/{self.dataid}_{self.vesicle_to_mito_filename}'
        self.curdatapath = self.vesicle_to_mito_filepath
        if not os.path.exists(self.curdatapath) or overlap == True:
            self.vesicle_to_mito_distance()
   
        dist3 = pd.read_csv(self.curdatapath)
        outplot.dist_isg2mito_hist(dist3, self.dataid, len(dist3))
        plotname = f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.vesicle_to_mito_plotname}'
        plt.savefig(plotname)
        plt.close()


    def vesicle_to_mito_distance(self):

        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)

        vesicle_mem_mask = mrcfile.open(self.vesiclefilepath, permissive=True).data
        mito_mem_mask = mrcfile.open(self.mitofilepath, permissive=True).data

        # resize
        zoom_degree = 0.10
        vesicle_mem_mask = scipy.ndimage.zoom(vesicle_mem_mask, zoom=zoom_degree, order=0)
        mito_mem_mask = scipy.ndimage.zoom(mito_mem_mask, zoom=zoom_degree, order=0)

        vesicle_edt_mask = scipy.ndimage.distance_transform_edt(vesicle_mem_mask)
        vesicle_edge_coords =np.array(np.where(vesicle_edt_mask == 1)).T #[[xyz], [xyz]]
        # vesicle_cent_coord = [np.mean(np.where(vesicle_edt_mask == 1)[0]), np.mean(np.where(vesicle_edt_mask == 1)[1]), np.mean(np.where(vesicle_edt_mask == 1)[2])]

        mito_edt_mask = scipy.ndimage.distance_transform_edt(mito_mem_mask)
        mito_edge_coords = np.array(np.where(mito_edt_mask == 1)).T #[[xyz], [xyz]]

        # vesicle_edge_coords = vesicle_edge_coords * np.array([1/zoom_degree, 1/zoom_degree, 1/zoom_degree,])
        # mito_edge_coords = mito_edge_coords * np.array([1/zoom_degree, 1/zoom_degree, 1/zoom_degree,])

        print('zoom ', zoom_degree)

        batchsize = 1024
        nn = len(vesicle_edge_coords) // batchsize +1
        distt = []

        vesicle_edge_coords, mito_edge_coords = vesicle_edge_coords.astype(np.float), mito_edge_coords.astype(np.float)
        vesicle_edge_coords, mito_edge_coords =torch.from_numpy(vesicle_edge_coords), torch.from_numpy(mito_edge_coords)
        vesicle_edge_coords, mito_edge_coords = vesicle_edge_coords.to(device), mito_edge_coords.to(device)


        for ii in range(nn):
            curisgcoord = vesicle_edge_coords[batchsize * ii : batchsize * (ii+1)]
            
            tempdistt = torch.cdist(curisgcoord, mito_edge_coords, compute_mode='use_mm_for_euclid_dist')
            tempdistt = tempdistt.cpu().numpy()
            for jj in range(len(tempdistt)):
                distt.extend([value for value in tempdistt[jj]])
            del tempdistt   # len = 375670784 *52 ??  # 46803968 * 26 

        print('len distt', len(distt)) 
        # 1211921105 for zoom 0.5   
        # 73498191 for zoom 0.25  
        # 29657790 for zoom 0.2
        # 12004180 for zoom 0.16
        # 1624105 for zoom 0.1


        # to nm
        avesize = np.average(voxel_size_xyz)
        distt = [dist * (avesize/10) * (1/zoom_degree) for dist in distt ] # nm

        distpd = pd.DataFrame(columns=['Distance'], data = distt)
        distpd.to_csv(self.curdatapath, index=False)


class vesicle_to_endoreticulum(process_organelles):
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'vesicle_to_endoreticulum'

    def step(self):
        organelletypelst = [self.vesicle_type, self.endo_reticulum_type]
        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} {self.name} data not exist.')

        else:
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()

            self.vesicle_to_endoreticulum_plot(overlap = self.plotdata_overlap)
        print(f'Done with {self.name}.')

    def vesicle_to_endoreticulum_plot(self, overlap):
        
        self.vesicle_to_endoreticulum_filepath = f'{self.datanewpath}/{self.dataid}_{self.vesicle_to_endoreticulum_filename}'
        self.curdatapath = self.vesicle_to_endoreticulum_filepath
        if not os.path.exists(self.curdatapath) or overlap == True:
            self.vesicle_to_endoreticulum_distance()
   
        dist3 = pd.read_csv(self.curdatapath)
        outplot.dist_isg2er_hist(dist3, self.dataid, len(dist3))
        plotname = f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.vesicle_to_endoreticulum_plotname}'
        plt.savefig(plotname)
        plt.close()


    def vesicle_to_endoreticulum_distance(self):
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)



        vesicle_mem_mask = mrcfile.open(self.vesiclefilepath, permissive=True).data
        er_mem_mask = mrcfile.open(self.endoreticulumfilepath, permissive=True).data


        # resize
        zoom_degree = 0.10

        vesicle_mem_mask = scipy.ndimage.zoom(vesicle_mem_mask, zoom=zoom_degree, order=0)
        er_mem_mask = scipy.ndimage.zoom(er_mem_mask, zoom=zoom_degree, order=0)


        vesicle_edt_mask = scipy.ndimage.distance_transform_edt(vesicle_mem_mask)
        vesicle_edge_coords =np.array(np.where(vesicle_edt_mask == 1)).T #[[xyz], [xyz]]
        # vesicle_cent_coord = [np.mean(np.where(vesicle_edt_mask == 1)[0]), np.mean(np.where(vesicle_edt_mask == 1)[1]), np.mean(np.where(vesicle_edt_mask == 1)[2])]

        # er_edtmask = scipy.ndimage.distance_transform_edt(er_mem_mask)
        er_coords = np.array(np.where(er_mem_mask != 0)).T #[[xyz], [xyz]]
        er_edge_coords = er_coords

        print('zoom ', zoom_degree)

        batchsize = 1024
        nn = len(vesicle_edge_coords) // batchsize +1
        distt = []

        vesicle_edge_coords, er_edge_coords = vesicle_edge_coords.astype(np.float), er_edge_coords.astype(np.float)
        vesicle_edge_coords, er_edge_coords =torch.from_numpy(vesicle_edge_coords), torch.from_numpy(er_edge_coords)
        vesicle_edge_coords, er_edge_coords = vesicle_edge_coords.to(device), er_edge_coords.to(device)

        for ii in range(nn):
            curisgcoord = vesicle_edge_coords[batchsize * ii : batchsize * (ii+1)]
            
            tempdistt = torch.cdist(curisgcoord, er_edge_coords, compute_mode='use_mm_for_euclid_dist')
            tempdistt = tempdistt.cpu().numpy()
            for jj in range(len(tempdistt)):
                distt.extend([value for value in tempdistt[jj]])
            del tempdistt   # len = 375670784 *52 ??  # 46803968 * 26 

        # print('len distt', len(distt)) 

        # to nm
        avesize = np.average(voxel_size_xyz)
        distt = [dist * (avesize/10) * (1/zoom_degree) for dist in distt ] # nm

        distpd = pd.DataFrame(columns=['Distance'], data = distt)
        distpd.to_csv(self.curdatapath, index=False)


class mito_to_endoreticulum(process_organelles):
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'mito_to_endoreticulum'

    def step(self):
        organelletypelst = [self.mito_type, self.endo_reticulum_type]
        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} {self.name} data not exist.')

        else:
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()

            self.mito_to_endoreticulum_plot(overlap = self.plotdata_overlap)
        print(f'Done with {self.name}.')


    def mito_to_endoreticulum_plot(self, overlap):
        
        self.mito_to_endoreticulum_filepath = f'{self.datanewpath}/{self.dataid}_{self.mito_to_endoreticulum_filename}'
        self.curdatapath = self.mito_to_endoreticulum_filepath
        if not os.path.exists(self.curdatapath) or overlap == True:
            self.mito_to_endoreticulum_distance()
   
        dist3 = pd.read_csv(self.curdatapath)
        outplot.dist_mito2er_hist(dist3, self.dataid, len(dist3))
        plotname = f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.mito_to_endoreticulum_plotname}'
        plt.savefig(plotname)
        plt.close()


    def mito_to_endoreticulum_distance(self):
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)

        mito_mem_mask = mrcfile.open(self.mitofilepath, permissive=True).data
        er_mem_mask = mrcfile.open(self.endoreticulumfilepath, permissive=True).data

        # resize
        zoom_degree = 0.10

        mito_mem_mask = scipy.ndimage.zoom(mito_mem_mask, zoom=zoom_degree, order=0)
        er_mem_mask = scipy.ndimage.zoom(er_mem_mask, zoom=zoom_degree, order=0)


        mito_edt_mask = scipy.ndimage.distance_transform_edt(mito_mem_mask)
        mito_edge_coords =np.array(np.where(mito_edt_mask == 1)).T #[[xyz], [xyz]]
        # vesicle_cent_coord = [np.mean(np.where(vesicle_edt_mask == 1)[0]), np.mean(np.where(vesicle_edt_mask == 1)[1]), np.mean(np.where(vesicle_edt_mask == 1)[2])]

        # er_edtmask = scipy.ndimage.distance_transform_edt(er_mem_mask)
        er_coords = np.array(np.where(er_mem_mask != 0)).T #[[xyz], [xyz]]
        er_edge_coords = er_coords

        print('zoom ', zoom_degree)

        batchsize = 1024
        nn = len(mito_edge_coords) // batchsize +1
        distt = []

        mito_edge_coords, er_edge_coords = mito_edge_coords.astype(np.float), er_edge_coords.astype(np.float)
        mito_edge_coords, er_edge_coords =torch.from_numpy(mito_edge_coords), torch.from_numpy(er_edge_coords)
        mito_edge_coords, er_edge_coords = mito_edge_coords.to(device), er_edge_coords.to(device)

        for ii in range(nn):
            curisgcoord = mito_edge_coords[batchsize * ii : batchsize * (ii+1)]
            
            tempdistt = torch.cdist(curisgcoord, er_edge_coords, compute_mode='use_mm_for_euclid_dist')
            tempdistt = tempdistt.cpu().numpy()
            for jj in range(len(tempdistt)):
                distt.extend([value for value in tempdistt[jj]])
            del tempdistt   # len = 375670784 *52 ??  # 46803968 * 26 


        # to nm
        avesize = np.average(voxel_size_xyz)
        distt = [dist * (avesize/10) * (1/zoom_degree) for dist in distt ] # nm

        distpd = pd.DataFrame(columns=['Distance'], data = distt)
        distpd.to_csv(self.curdatapath, index=False)


class actin_angle_distance_pair(process_organelles):
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'actin_angle_distance_pair'

    def step(self):
        organelletypelst = [self.actin_type, self.vesicle_type]
        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} {self.name} data not exist.')

        else:
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()
            self.actin_angle_distance_pair_plot(overlap = self.plotdata_overlap)
        print(f'Done with {self.name}.')


    def actin_angle_distance_pair_plot(self, overlap):
        
        self.actin_angle_distance_pair_filepath = f'{self.datanewpath}/{self.dataid}_{self.actin_angle_distance_pair_filename}'
        self.curdatapath = self.actin_angle_distance_pair_filepath
        if not os.path.exists(self.curdatapath) or overlap == True:
            self.actin_angle_distance_pair_map()


        with open(f'{self.actinfilepath}', 'r') as f:
            filament_data = json.load(f)

        dist3 = pd.read_csv(self.curdatapath)
        # print(dist3)
        filaments_min_d_lst, filaments_angle_lst =  dist3['Distances'], dist3['Angles']
        if dist3.shape[0] == 0:
            print(f'{self.dataid} only have one filament.')
        else:
            outplot.vis_actindistangle(filaments_min_d_lst, filaments_angle_lst, self.dataid, len(filament_data), len(filaments_min_d_lst))
            plotname = f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.actin_angle_distance_pair_plotname}'
            plt.savefig(plotname)
            # plt.show()
            plt.close()


    def actin_angle_distance_pair_map(self):

        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)
        with open(f'{self.actinfilepath}', 'r') as f:
            filament_data = json.load(f)
        vesicle_mem_mask = mrcfile.open(self.vesiclefilepath, permissive=True).data
        imgsize = vesicle_mem_mask.shape
        filament_data_new = distances_generator.shift_bias(filament_data, shift_xml_xyz)


        # filaments_min_d_lst, filaments_angle_lst, = distances_generator.actinangle2dist(filament_data_new, imgsize, voxel_size_xyz)
        filaments_min_d_lst, filaments_angle_lst, = self.actinangle2dist(filament_data_new, imgsize, voxel_size_xyz)
        
        distangle_savepd = pd.DataFrame(data= {'Distances':filaments_min_d_lst,
                                          'Angles': filaments_angle_lst })
        distangle_savepd.to_csv(self.curdatapath, index=False)



    def actinangle2dist(self, filament_coord_extend_dict, imgsize, voxel_size_xyz):
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
                            '''
                            update min_dist, cur_angle
                            '''
                            temp_d = min_d
                            vect2 = rest_filament_vect_extend_lst[kk][min_d_idx]
                            if (vect1 == vect2).all():
                                print(vect1, vect2)
                                continue
                            else:
                                vector_dot_product = np.dot(vect1, vect2)  # when nan, 1600
                                arccos = np.arccos(vector_dot_product / (np.linalg.norm(vect1) * np.linalg.norm(vect2)))
                                angle = np.degrees(arccos)
                                if angle >= 180:
                                    temp_angle = 0
                                elif angle > 90:
                                    temp_angle = 180 - angle
                                else:
                                    temp_angle = angle
                                    if angle >0 and angle <=90:
                                        pass
                                    # else: 
                                        # print(temp_angle == float(nan))

                                        # print(temp_angle>=0 or temp_angle<=180)
                                        # print(temp_angle)
                                        # print(temp_angle == nan)
                                        # print(temp_angle == np.nan)
                                        # print(type(temp_angle))

                                if temp_angle>=0 or temp_angle<=180:  # exclude nan data
                                    min_dist = temp_d
                                    cur_angle = temp_angle
                                # else:
                                #     print(11, temp_d, temp_angle)

                    filament_dist_lst.append(min_dist)
                    filament_angle_lst.append(cur_angle)

                filaments_min_d_lst.extend(filament_dist_lst)
                filaments_angle_lst.extend(filament_angle_lst)


        return  filaments_min_d_lst, filaments_angle_lst



class focal_to_focal(process_organelles):
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'focal_to_focal'
        

    def step(self):
        print(f'dataid {self.dataid}')
        organelletypelst = [self.focal_adhesion_type, self.focal_adhesion_type]
        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} {self.name} data not exist.')

        else:
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()

            self.focal_to_focal_plot(overlap = self.plotdata_overlap)
        print(f'Done with {self.name}.')


    def focal_to_focal_plot(self,  overlap):

        self.focal_to_focal_distance_filepath = f'{self.datanewpath}/{self.dataid}_{self.focal_to_focal_distance_filename}'
        self.curdatapath = self.focal_to_focal_distance_filepath
        if not os.path.exists(self.curdatapath) or overlap == True:
            self.focal_to_focal_distance()


        fa_data = pd.read_csv(self.focaladhesionfilepath)

        dist3 = pd.read_csv(self.focal_to_focal_distance_filepath)

        outplot.dist_focal2focal_hist(dist3, self.dataid, len(fa_data),)
        plotname = f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.focal_to_focal_distance_plotname}'
        plt.savefig(plotname)
        # plt.show()
        plt.close()
     

    def focal_to_focal_distance(self):
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)
        
        fa_data = pd.read_csv(self.focaladhesionfilepath)
        fa_data = [[fa_data['X'][i], fa_data['Y'][i], fa_data['Z'][i] ] for i in range(fa_data.shape[0]) ]

        Fa_data_new = fa_data
        dist = distances_generator.FA_pairdist(Fa_data_new, voxel_size_xyz)  # nm
        dist_savepd = pd.DataFrame(columns=['Distance'], data= dist)

        # filament2filament_dist_lst_savepd.loc['distance'] = filament2filament_dist_lst
        dist_savepd.to_csv(f'{self.datanewpath}/{self.dataid}_{self.focal_to_focal_distance_filename}', index=False)


    def reflash_inputdata(self):  # for updated false coord to refresh input
        organelletypelst = [self.focal_adhesion_type, self.focal_adhesion_type]
        self.organelletypelist = organelletypelst
        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} {self.name} data not exist.')
        else:
            self.generate_processed_data()


class focal_aggregated_actin_angle(process_organelles):
    '''
    Actins aggregated within one focal adhesion are regareded as one cluster.
    calculate the angle distribution of the actin end to focal adhesion angle.
    distributions.
    '''
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'focal_aggregated_actin_angle'
        self.dist_threshold = 100 # nm
    def step(self):
        organelletypelst = [self.focal_adhesion_type, self.actin_type]
        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} {self.name} data not exist.')
        else:
            if not self.check_data_updated(organelletypelst):
                # print(11)
                self.generate_processed_data()

            # self.focal_aggregated_actin_angle_check()
            self.focal_aggregated_actin_angle_plot(overlap = self.plotdata_overlap)
        print(f'Done with {self.name}.')

    def focal_aggregated_actin_angle_plot(self, overlap):
        # print(1111)
        self.focal_aggregated_actin_angle_map_filepath = f'{self.datanewpath}/{self.dataid}_{self.focal_aggregated_actin_angle_map_filename}'
        self.curdatapath = self.focal_aggregated_actin_angle_map_filepath
        if not os.path.exists(self.curdatapath) or overlap == True:
            self.focal_aggregated_actin_angle_map()


        # focal_aggregated_actin_map_filename

        with open(f'{self.actinfilepath}', 'r') as f:
            actin_data = json.load(f)
        fa_data = pd.read_csv(self.focaladhesionfilepath)


        angle_dist_pairs = pd.read_csv(self.curdatapath) ## angles

        angle, distance =  angle_dist_pairs['Angle'], angle_dist_pairs['Distance']
        

        outplot.fa_actin_angle_dist_map_plot(angle_dist_pairs, self.dataid,)
        plotname = f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.focal_aggregated_actin_angle_map_plotname}'
        plt.savefig(plotname)
        plt.show()
        plt.close()
     

    def focal_aggregated_actin_angle_map(self):
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)
        
        with open(self.actinfilepath, 'r') as f:
            actin_data = json.load(f)
        fa_data = pd.read_csv(self.focaladhesionfilepath)

        actin_data = distances_generator.shift_bias(actin_data, shift_xml_xyz)

        # anglelst, distlst = distances_generator.actin2Fa_angle(actin_data, fa_data, voxel_size_xyz, self.dist_threshold)
        anglelst, distlst = self.actin2Fa_angle(actin_data, fa_data, voxel_size_xyz, self.dist_threshold)
        angle_savepd = pd.DataFrame({'Angle':anglelst, 'Distance':distlst})

        # filament2filament_dist_lst_savepd.loc['distance'] = filament2filament_dist_lst
        angle_savepd.to_csv(f'{self.datanewpath}/{self.dataid}_{self.focal_aggregated_actin_angle_map_filename}', index=False)
        

    def focal_aggregated_actin_angle_check(self):
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)
        with open(self.actinfilepath, 'r') as f:
            actin_data = json.load(f)
        fa_data = pd.read_csv(self.focaladhesionfilepath)

        actin_data = distances_generator.shift_bias(actin_data, shift_xml_xyz)

        anglelst, distlst = distances_generator.actin2Fa_angle_check(actin_data, fa_data, voxel_size_xyz, self.dist_threshold)
        angle_savepd = pd.DataFrame({'Angle':anglelst, 'Distance':distlst})  

        # filament2filament_dist_lst_savepd.loc['distance'] = filament2filament_dist_lst
        angle_savepd.to_csv(f'{self.datanewpath}/{self.dataid}_{self.focal_aggregated_actin_angle_map_filename}', index=False)



    def actin2Fa_angle(self,filament_coords, fa_data, voxelsize, dist_threshold):
        
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
            # print('matrix shape',dists_matrix.shape)
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
                    fa_actin_dist_lst.append(min_d)
                    fa_actin_actincoords.append([curfilaments_coords[0], curfilaments_coords[-1]][min_d_idx])
                    fa_actin_actinidx.append(fila_idx)
                    fa_actin_faidx.append(i)
                # currespond_angle_for_filament[i] = angle



        return fa_actin_angle_lst, fa_actin_dist_lst


class focal_aggregated_actin_angles_pair(process_organelles):
    '''
    Actins aggregated within one focal adhesion are regareded as one cluster.
    calculate the angle distribution of the actin to other actin and the angle to mem
    output: 2d plot angle to other angle Vs angle to mem
    
    '''

    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'focal_aggregated_actin_angles_pair'
        self.dist_threshold = 100 # nm
    def step(self):
        organelletypelst = [self.focal_adhesion_type, self.actin_type]
        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} {self.name} data not exist.')
        else:
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()

            self.focal_aggregated_actin_angles_pair_plot(overlap = self.plotdata_overlap)
        print(f'Done with {self.name}.')



    def focal_aggregated_actin_angles_pair_plot(self, overlap):

        self.focal_aggregated_actin_angles_pair_filepath = f'{self.datanewpath}/{self.dataid}_{self.focal_aggregated_actin_angles_pair_filename}'
        self.curdatapath = self.focal_aggregated_actin_angles_pair_filepath
        if not os.path.exists(self.curdatapath) or overlap == True:
            self.focal_aggregated_actin_angles_pair_map()

        with open(f'{self.actinfilepath}', 'r') as f:
            actin_data = json.load(f)
        fa_data = pd.read_csv(self.focaladhesionfilepath)

        angle_angle_pairs = pd.read_csv(self.curdatapath) ## angles
        # angle, distance =  angle_angle_pairs['Angle actin'], angle_angle_pairs['Angle mem']

        
        outplot.fa_actin_angle_angle_pair_plot(angle_angle_pairs, self.dataid,)
        plotname = f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.focal_aggregated_actin_angles_pair_plotname}'
        plt.savefig(plotname)
        # plt.show()
        plt.close()
     
        self.angle_to_mem_plot()

     
    def focal_aggregated_actin_angles_pair_map(self):
        
        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)
        
        with open(self.actinfilepath, 'r') as f:
            actin_data = json.load(f)
        with open(self.actinanglefilepath, 'r') as f:
            actin_angle_data = json.load(f)        
        fa_data = pd.read_csv(self.focaladhesionfilepath)





        anglelst1, anglelst2 = distances_generator.Fa_aggragated_actinactin_angle(actin_data, actin_angle_data, fa_data, voxel_size_xyz, self.dist_threshold)
        angle_savepd = pd.DataFrame({'Angle actin':anglelst1, 'Angle mem':anglelst2})

        # filament2filament_dist_lst_savepd.loc['distance'] = filament2filament_dist_lst
        angle_savepd.to_csv(f'{self.datanewpath}/{self.dataid}_{self.focal_aggregated_actin_angles_pair_filename}', index=False)



    def angle_to_mem_plot(self):

        with open(self.actinanglefilepath, 'r') as f:
            actin_angle_data = json.load(f)    

        anglelst = list( actin_angle_data.values() )

        anglepd = pd.DataFrame({'angle1':anglelst, 'angle2':anglelst})

        df = anglepd
        plt.hist(anglelst, bins=20, rwidth=0.8, range=(0,90))


        # sns.jointplot(x=df['angle1'], y=df['angle2'],
        #     data=df,
        #     color='k',
        #     kind='scatter',             #reg添加线性回归线
        #     height=16,
        #     ratio=5,
        #     marginal_kws=dict(bins=20, rug=True,))

        plt.title(f'{self.dataid} fa_actin_angle_angle_pair_plot plot \n point num {df.shape[0]}')

    # plt.show()

        plotname = f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_anglecheck.png'
        plt.savefig(plotname)
        # plt.show()
        plt.close()


class focal_to_ca(process_organelles):

    def __init__(self, dataid) -> None:
        super().__init__(dataid)
        self.name = f'focal_to_ca'
        
    def step(self):
        print(f'dataid {self.dataid}')
        organelletypelst = [self.ca_type, self.focal_adhesion_type]
        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} {self.name} data not exist.')

        else:
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()
                # print('line1559')
            self.focal_to_ca_plot(overlap = self.plotdata_overlap)
        print(f'Done with {self.name}.')


    def focal_to_ca_plot(self, overlap):
        
        self.focal_to_ca_distance_filepath = f'{self.datanewpath}/{self.dataid}_{self.focal_to_ca_distance_filename}'
        self.curdatapath = self.focal_to_ca_distance_filepath
        if not os.path.exists(self.curdatapath) or overlap == True:
            self.focal_to_ca_distance()

        ca_data = pd.read_csv(self.cafilepath)
        fa_data = pd.read_csv(self.focaladhesionfilepath)

        dist3 = pd.read_csv(self.focal_to_ca_distance_filepath)

        outplot.dist_focal2ca_hist(dist3, self.dataid, len(fa_data), len(ca_data))
        plotname = f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.focal_to_ca_distance_plotname}'
        plt.savefig(plotname)
        # plt.show()
        plt.close()
     

    def focal_to_ca_distance(self):

        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)
        
        fa_data = pd.read_csv(self.focaladhesionfilepath)
        fa_data = [[fa_data['X'][i], fa_data['Y'][i], fa_data['Z'][i] ] for i in range(fa_data.shape[0]) ]
        ca_data = pd.read_csv(self.cafilepath)
        ca_data  = [[ca_data['X'][i], ca_data['Y'][i], ca_data['Z'][i] ] for i in range(ca_data.shape[0]) ]


        Fa_data_new = fa_data
        dist = distances_generator.FA2ca_dist(Fa_data_new, ca_data, voxel_size_xyz)  # nm
        dist_savepd = pd.DataFrame(columns=['Distance'], data= dist)

        # filament2filament_dist_lst_savepd.loc['distance'] = filament2filament_dist_lst
        dist_savepd.to_csv(f'{self.datanewpath}/{self.dataid}_{self.focal_to_ca_distance_filename}', index=False)





def main():
    '''
    get dataid
    get catagories for calculation
    get list if catagory available
    check if filtered data exist
    generate figure 
    '''
    dataidlst = preprocess.read_dataid()
    print('The dataidlst is: ', dataidlst)
    # dataidlst = dataidlst[:1]
    # 20220325_2.8_16
    # dataidlst =  [ '20220918_30_004',]
    #     # dataidlst = ['20210517_30_009', '20210517_30_010', '20210517_30_011', '20210517_30_012', '20210517_30_013', '20210517_30_014', '20220918_30_001', '20220918_30_004', '20220219_5_002', '20220219_5_003', '20220305_5_001', '20220305_5_002', '20220305_5_003', '20220305_5_004', '20220305_5_005', '20220717_5_001', '20210323_2.8_006', '20210323_2.8_007', '20210323_2.8_008', '20210323_2.8_010', '20210323_2.8_011', '20210323_2.8_012', '20210323_2.8_014', '20210323_2.8_016', '20220325_30_10', '20220325_30_11', '20220325_30_3', '20220325_30_4', '20220325_5_1', '20220325_5_2', '20220325_5_3', '20220325_5_4', '20220325_5_5', '20220325_5_6', '20220325_2.8_16', '20220325_2.8_20', '20220326_2.8_1', '20220326_2.8_2']
    # dataidlst = dataidlst[:1]

    for dataid in dataidlst: # check all data 
        print(dataid)
        aaa = check_data(dataid)
        aaa.step()


    # for dataid in dataidlst:
    #     print(dataid)

    #     a2vv = actin_to_vesicle_vis(dataid)
    #     a2vv.step()
    #     # mt 8nm
    #     # actin 4 nm

    #     m2vv = microtube_to_vesicle_vis(dataid)
    #     m2vv.step()


    # Actin-Vesicle_dist-angle pair: X轴是Actin-vesicle dist，Y轴是对应点的actin所成角度，角度数值要用90减一下才是从XY平面的角度。画一个散点或者热力图。
    # 另外在之前actin-vesicle_dist中距离近的那个点是哪个数据的哪一根actin，编号就用input Actin.xml中的编号就行，我可以直接去tomo里面画图找出来。


    for dataid in dataidlst:
        print(dataid)

        # vesicle actin microtube mito endo_reticulum
        # print(catagory)
        
        a2a = actin_to_actin(dataid)
        a2a.step()

        # a2v = actin_to_vesicle(dataid)
        # a2v.step()

        a2v_d2a = actin_to_vesicle_dist_angle_pair(dataid)
        a2v_d2a.step()

        mt2a = microtube_to_actin(dataid)
        mt2a.step()

        # a2mi = actin_to_mito(dataid) 
        # a2mi.step()

        # a2e = actin_to_endoreticulum(dataid)
        # a2e.step()

        # mt2i = microtube_to_vesicle(dataid)
        # mt2i.step()

        # mt2mi = microtube_to_mito(dataid)
        # mt2mi.step()



        # mt2er = microtube_to_endoreticulum(dataid)
        # mt2er.step()

        # i2mi =vesicle_to_mito(dataid)
        # i2mi.step()

        # i2er = vesicle_to_endoreticulum(dataid)
        # i2er.step()

        # mi2er = mito_to_endoreticulum(dataid)
        # mi2er.step()

        # acta2dp = actin_angle_distance_pair(dataid)
        # # acta2dp.plotdata_overlap = True
        # acta2dp.step()

        # fa2fa = focal_to_focal(dataid)
        # # fa2fa.plotdata_overlap = True
        # ## fa2fa.reflash_inputdata()
        # fa2fa.step()

        # faaa = focal_aggregated_actin_angle(dataid)
        # # faaa.plotdata_overlap = True
        # faaa.step()

        # faaap = focal_aggregated_actin_angles_pair(dataid)
        # # faaap.plotdata_overlap = True
        # faaap.step()


        # fa2ca = focal_to_ca(dataid)
        # # fa2ca.plotdata_overlap = True
        # fa2ca.step()


    #     # vis actin/microtube  

    print('Done.')


if __name__ == '__main__':
    main()
    # single_data()
 


# %%
# dataidlst = ['20210517_30_009', '20210517_30_010', '20210517_30_011', '20210517_30_012', '20210517_30_013', '20210517_30_014', '20220918_30_001', '20220918_30_004', '20220219_5_002', '20220219_5_003', '20220305_5_001', '20220305_5_002', '20220305_5_003', '20220305_5_004', '20220305_5_005', '20220717_5_001', '20210323_2.8_006', '20210323_2.8_007', '20210323_2.8_008', '20210323_2.8_010', '20210323_2.8_011', '20210323_2.8_012', '20210323_2.8_014', '20210323_2.8_016', '20220325_30_10', '20220325_30_11', '20220325_30_3', '20220325_30_4', '20220325_5_1', '20220325_5_2', '20220325_5_3', '20220325_5_4', '20220325_5_5', '20220325_5_6', '20220325_2.8_16', '20220325_2.8_20', '20220326_2.8_1', '20220326_2.8_2']

# print(len(dataidlst))
# %%
