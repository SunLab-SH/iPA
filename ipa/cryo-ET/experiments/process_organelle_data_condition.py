# %%
from cmath import nan
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

print(torch.__version__)
print(torch.cuda.is_available())

device = torch.device('cuda')

#%%

def obtain_dir():
    for maindir, subdir, file_name_list in os.walk(os.path.join(arg.root_dir, arg.data_root_dir), topdown=False):
        locate = np.array(subdir)

    subdir1 = [  os.path.join(arg.root_dir, arg.data_root_dir, subdir) for subdir in locate]
    # print(subdir1)
    subdir2 = []
    for dir in subdir1:
        for maindir, subdir, file_name_list in os.walk(dir, topdown=False):
            condition = np.array(subdir)
        for singlecondition in condition:
            subdir2.append( os.path.join(dir, singlecondition) )
    # print(subdir2)
    subdir3 = []
    for dir in subdir2:
        for maindir, subdir, file_name_list in os.walk(dir, topdown=False):
            replica = np.array(subdir)
        for singlereplica in replica:

            subdir3.append( os.path.join(dir, singlereplica))
    # print(subdir3)

    return subdir2, subdir3  # subdir2 condition path   subdir3 data path 



class process_organelles_condition(process_organelles):
    def __init__(self,dataid):
        super(process_organelles_condition, self).__init__(dataid)
        self.obtain_conditionlst()
        self.set_conditionplotname()
    def obtain_conditionlst(self):
        self.conditionpathlst, self.filepathlst = preprocess.obtain_dir()
        self.conditionlst = [ path.split('\\')[-1] for path in self.conditionpathlst ]
        tempconditionlst = []
        for condition in self.conditionlst:
            if condition not in tempconditionlst:
                tempconditionlst.append(condition)
        self.conditionlst = tempconditionlst
        del tempconditionlst
        self.filelst  = [ path.split('\\')[-1] for path in self.filepathlst ]

    def obtain_location(self, conditionname):
        'self.condtion  to find location'
        for pathname in self.conditionpathlst:
            if conditionname in pathname:
                self.location =  pathname.split('\\')[-2]


        
    def obtain_condition_data(self, conditionword, filetypeword):
        filelist = self.get_filelist(os.path.join(arg.root_dir, arg.processed_data_dir),[])
        
        filelist_new = list()

        for filepath in filelist:
            if conditionword in filepath and filetypeword in filepath:
                filelist_new.append(filepath)


        self.filteredfilelst = filelist_new

        self.curconditionfilenamelst = [ path.split('\\')[-2] for path in self.filteredfilelst ]
        #print('The filelist of focal_ca is:', filelist_new)
       
        self.obtain_location(conditionword)


    def count_actinnum(self, conditionname):
        conditionword = conditionname
        filetypeword = self.actinfilename

        filelist = self.get_filelist(os.path.join(arg.root_dir, arg.processed_data_dir),[])
        filelist_new = list()

        for filepath in filelist:
            if conditionword in filepath and filetypeword in filepath:
                filelist_new.append(filepath)
        
        filteredfilelst = filelist_new

        num = 0
        for name in filteredfilelst:
            with open(f'{name}', 'r') as f:
                filament_data = json.load(f)
    
            num += len(filament_data)

        self.curcondition_actinnum = num



    
    def set_conditionplotname(self):

        self.actin_surround_actin_condition_plotname = f'actin_surrounding_actin_probability_density_condition.pdf'

        self.actin_to_vesicle_condition_plotname = f'actin_distance_to_vesicle_distance_condition.pdf' 

        self.microtube_surround_actin_condition_plotname = f'microtube_surrounding_actin_probability_density_condition.pdf'

        self.actin_to_mito_condition_plotname = f'actin_to_mito_distance_condition.pdf'

        self.actin_to_endoreticulum_condition_plotname = f'actin_to_endoreticulum_distance_condition.pdf'

        self.microtube_to_vesicle_condition_plotname = f'microtube_to_vesicle_distance_condition.pdf'

        self.microtube_to_mito_condition_plotname = f'microtube_to_mito_distance_condition.pdf'

        self.microtube_to_endoreticulum_condition_plotname = f'microtube_to_endoreticulum_distance_condition.pdf'

        self.vesicle_to_mito_condition_plotname = f'vesicle_to_mito_distance_condition.pdf'

        self.vesicle_to_endoreticulum_condition_plotname = f'vesicle_to_endoreticulum_distances_condition.pdf'

        self.mito_to_endoreticulum_condition_plotname = f'mito_to_endoreticulum_distances_condition.pdf'

        self.actin_angle_distance_pair_condition_filename = f'actin_angle_distance_pair_condition.pdf'

        self.focal_to_focal_distance_condition_plotname = f'focal_to_focal_distance_condition.pdf'

        self.focal_aggregated_actin_angle_map_condition_plotname = f'focal_aggregated_actin_angle_map_condition.pdf'

        self.focal_aggregated_actin_angles_pair_condition_plotname = f'focal_aggregated_actin_angles_pair_condition.pdf'

        self.focal_to_ca_distance_condition_plotname = f'focal_to_ca_distance_condition.pdf'



class actin_to_actin_condition(process_organelles_condition):
    def __init__(self,dataid):
        super(actin_to_actin_condition,self).__init__(dataid)
        self.name = f'actin_to_actin_condition'


    def step(self):
        print('conditions', self.conditionlst)

        for conditionname in self.conditionlst:
            print(conditionname)
            self.conditionname = conditionname
            self.obtain_condition_data(conditionname, self.actin_surround_actin_vect_map_filename)
            print(self.location)
            # print(1001,conditionname, self.filteredfilelst)
            if len (self.filteredfilelst) != 0:
                self.actin_surround_actin_map_condition_plot()

        print(f'Done with {self.name}.')


    def actin_surround_actin_map_condition_plot(self):
        filament_vect_map_condition = []
        for dataname in self.filteredfilelst:
            with open(f'{dataname}', 'r', encoding='utf-8') as f:
                filament_vect_map_lst = json.load(f)
            print('filament num,', len(filament_vect_map_lst))  ## list[x,y,z]

            filament_vect_map_lst_new = []
            vect_count_n = 0
            for filament in filament_vect_map_lst:
                curcoord = [ np.array(vect) for vect in filament ]
                filament_vect_map_lst_new.append(curcoord)
                vect_count_n += len(curcoord)

            filament_vect_map_lst = filament_vect_map_lst_new
            filament_vect_map_condition.append(filament_vect_map_lst)


        outplot.dist2fila_condition_plot(filament_vect_map_condition, self.conditionname, self.curconditionfilenamelst)
        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_{self.actin_surround_actin_condition_plotname}'
        print(plotname)
        plt.tight_layout()
        # plt.clim()
        plt.clim(vmin = 0, vmax = 11000)
        plt.savefig(plotname)
        plt.show()
        plt.close()



class actin_to_vesicle_condition(process_organelles_condition):
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'actin_to_vesicle_condition'


    def step(self):
        print('conditions', self.conditionlst)

        for conditionname in self.conditionlst:
            print(conditionname)
            self.conditionname = conditionname
            self.obtain_condition_data(conditionname, self.actin_to_vesicle_filename)
            # print(1001,conditionname, self.filteredfilelst)
            if len (self.filteredfilelst) != 0:
                self.actin_to_vesicle_condition_plot()

        print(f'Done with {self.name}.')



    def actin_to_vesicle_condition_plot(self):
        dist_condition_collection = []
        distss_in_one_condition = []
        distss_prob_in_one_condition = []
        bins_ = 50

        for dataname in self.filteredfilelst:
            distances = pd.read_csv(f'{dataname}')
            distances = distances['Distance']
            distances_new =distances.tolist()
            print(dataname)
            print('max distance', np.max(distances))
            print('num', len(distances_new))
            dist_condition_collection.extend(distances_new)

        
            disttt, bins = np.histogram(distances,bins=bins_,range=[0,1200])
            distss_in_one_condition.append(disttt)
            distss_prob_in_one_condition.append(disttt / np.sum(disttt))
            print(np.sum(distss_prob_in_one_condition[-1]))
        

        disttt, _ = np.histogram(dist_condition_collection,bins=bins_,range=[0,1200])
        distss_prob_in_one_condition = [dist / np.sum(disttt) for dist in disttt ]
        
        outplot.dist2mem_condition_plot(distss_prob_in_one_condition, self.conditionname, self.curconditionfilenamelst, bins=bins_, x_lim = 1200)
        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_{self.actin_to_vesicle_condition_plotname}'
        plt.tight_layout()
        plt.savefig(plotname)
        plt.show()
        plt.close()

        outplot.dist2mem_condition_plot2(dist_condition_collection, self.conditionname, self.curconditionfilenamelst, bins=bins_, x_lim = 1200)
        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_2_{self.actin_to_vesicle_condition_plotname}'
        # add plot here, and change name above avoid overwriting
        plt.tight_layout()
        plt.savefig(plotname)
        plt.show()
        plt.close()



class microtube_to_actin_condition(process_organelles_condition):
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'microtube_to_actin_condition'


    def step(self):
        print('conditions', self.conditionlst)

        for conditionname in self.conditionlst:
            print(conditionname)
            self.conditionname = conditionname
            self.obtain_condition_data(conditionname, self.microtube_surround_actin_vect_map_filename)
            # print(1001,conditionname, self.filteredfilelst)
            if len (self.filteredfilelst) != 0:
                self.microtube_to_actin_map_condition_plot()

        print(f'Done with {self.name}.')


    def microtube_to_actin_map_condition_plot(self):
        filament_vect_map_condition = []
        for dataname in self.filteredfilelst:
            with open(f'{dataname}', 'r', encoding='utf-8') as f:
                filament_vect_map_lst = json.load(f)
            print('MT num,', len(filament_vect_map_lst))  ## list[x,y,z]

            filament_vect_map_lst_new = []
            vect_count_n = 0
            for filament in filament_vect_map_lst:
                curcoord = [ np.array(vect) for vect in filament ]
                filament_vect_map_lst_new.append(curcoord)
                vect_count_n += len(curcoord)

            filament_vect_map_lst = filament_vect_map_lst_new
            filament_vect_map_condition.append(filament_vect_map_lst)


        figname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_{self.microtube_surround_actin_condition_plotname}'
        outplot.dist_fila2MT_condition_plot(filament_vect_map_condition, self.conditionname, self.curconditionfilenamelst)
        plt.tight_layout()
        plt.clim(vmin = 0, vmax = 100)
        plt.savefig(figname)
        plt.show()
        plt.close()



class actin_to_mito_condition(process_organelles_condition):
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'actin_to_mito_condition'


    def step(self):
        print('conditions', self.conditionlst)

        for conditionname in self.conditionlst:
            print(conditionname)
            self.conditionname = conditionname
            self.obtain_condition_data(conditionname, self.actin_to_mito_filename)
            # print(1001,conditionname, self.filteredfilelst)
            if len (self.filteredfilelst) != 0:
                self.actin_to_mito_condition_plot()

        print(f'Done with {self.name}.')



    def actin_to_mito_condition_plot(self):

        dist_condition_collection = []
        distss_in_one_condition = []
        distss_prob_in_one_condition = []
        
        bins_ = 12
        for dataname in self.filteredfilelst:
            distances = pd.read_csv(f'{dataname}')
            distances = distances['Distance']
            distances_new =distances.tolist()
            print(dataname)
            print('max distance', np.max(distances))
            print('num', len(distances_new))
            dist_condition_collection.extend(distances_new)

        
            disttt, bins = np.histogram(distances,bins=bins_, range=[0,1200])
            distss_in_one_condition.append(disttt)
            distss_prob_in_one_condition.append(disttt / np.sum(disttt))
            print(np.sum(distss_prob_in_one_condition[-1]))

        disttt, _ = np.histogram(dist_condition_collection,bins=bins_,range=[0,1200])
        distss_prob_in_one_condition = [dist / np.sum(disttt) for dist in disttt ]
        distss_prob_in_one_condition2 = copy.deepcopy(distss_prob_in_one_condition)

        outplot.dist2mitomem_condition_plot(distss_prob_in_one_condition, self.conditionname, self.curconditionfilenamelst, bins=bins_, x_lim = 1200)
        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_{self.actin_to_mito_condition_plotname}'
        plt.tight_layout()
        plt.savefig(plotname)
        plt.show()
        plt.close()

        outplot.dist2mitomem_condition_plot2(dist_condition_collection, self.conditionname, self.curconditionfilenamelst, bins=bins_, x_lim = 1200)
        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_2_{self.actin_to_mito_condition_plotname}'
        plt.tight_layout()
        plt.savefig(plotname)
        # plt.show()

        plt.close()


class actin_to_endoreticulum_condition(process_organelles_condition):
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'actin_to_endoreticulum_condition'


    def step(self):
        print('conditions', self.conditionlst)
        for conditionname in self.conditionlst:
            print(conditionname)
            self.conditionname = conditionname
            self.obtain_condition_data(conditionname, self.actin_to_endoreticulum_filename)
            # print(1001,conditionname, self.filteredfilelst)
            if len (self.filteredfilelst) != 0:
                self.actin_to_endoreticulum_condition_plot()

        print(f'Done with {self.name}.')



    def actin_to_endoreticulum_condition_plot(self):

        dist_condition_collection = []
        distss_in_one_condition = []
        distss_prob_in_one_condition = []
        bins_ = 12


        for dataname in self.filteredfilelst:
            distances = pd.read_csv(f'{dataname}')
            distances = distances['Distance']
            distances_new =distances.tolist()
            print(dataname)
            print('max distance', np.max(distances))
            print('num', len(distances_new))
            dist_condition_collection.extend(distances_new)

        
            disttt, bins = np.histogram(distances,bins=bins_, range=[0,1200])
            distss_in_one_condition.append(disttt)
            distss_prob_in_one_condition.append(disttt / np.sum(disttt))
            print(np.sum(distss_prob_in_one_condition[-1]))
        

        disttt, _ = np.histogram(dist_condition_collection,bins=bins_,range=[0,1200])
        distss_prob_in_one_condition = [dist / np.sum(disttt) for dist in disttt ]

        outplot.dist2ermem_condition_plot(distss_prob_in_one_condition, self.conditionname, self.curconditionfilenamelst, bins=bins_, x_lim = 1200)
        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_{self.actin_to_endoreticulum_condition_plotname}'
        plt.tight_layout()
        plt.savefig(plotname)
        plt.show()
        plt.close()

        outplot.dist2ermem_condition_plot2(dist_condition_collection, self.conditionname, self.curconditionfilenamelst, bins=bins_, x_lim = 1200)
        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_2_{self.actin_to_endoreticulum_condition_plotname}'
        plt.tight_layout()
        plt.savefig(plotname)
        plt.show()
        plt.close()

class microtube_to_vesicle_condition(process_organelles_condition):
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'microtube_to_vesicle_condition'


    def step(self):
        print('conditions', self.conditionlst)

        for conditionname in self.conditionlst:
            print(conditionname)
            self.conditionname = conditionname
            self.obtain_condition_data(conditionname, self.microtube_to_vesicle_filename)
            # print(1001,conditionname, self.filteredfilelst)
            if len (self.filteredfilelst) != 0:
                self.microtube_to_vesicle_condition_plot()

        print(f'Done with {self.name}.')


    def microtube_to_vesicle_condition_plot(self):
        dist_condition_collection = []
        distss_in_one_condition = []
        distss_prob_in_one_condition = []
        bins_ = 50

        for dataname in self.filteredfilelst:
            distances = pd.read_csv(f'{dataname}')
            distances = distances['Distance']
            distances_new =distances.tolist()
            print(dataname)
            print('max distance', np.max(distances))
            print('num', len(distances_new))
            dist_condition_collection.extend(distances_new)

        
            disttt, bins = np.histogram(distances,bins=bins_,range=[0,1200])
            distss_in_one_condition.append(disttt)
            distss_prob_in_one_condition.append(disttt / np.sum(disttt))
            print(np.sum(distss_prob_in_one_condition[-1]))
        

        if len(dist_condition_collection) != 0:
            max_= np.ceil(np.max(dist_condition_collection)/100) * 100
        else:
            max_ = 1200
        disttt, _ = np.histogram(dist_condition_collection,bins=bins_,range=[0,max_])

        distss_prob_in_one_condition = [dist / np.sum(disttt) for dist in disttt ]
        
        outplot.dist2mem_condition_plot_notitle(distss_prob_in_one_condition, self.conditionname, self.curconditionfilenamelst, bins=bins_, x_lim = max_)
        plt.title(f'{self.conditionname}_condition_MT2isg_distance_density_distribution bin = {bins_}\n  {[i for i in self.curconditionfilenamelst]}')
        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_{self.microtube_to_vesicle_condition_plotname}'
        plt.tight_layout()
        plt.savefig(plotname)
        plt.show()
        plt.close()

        outplot.dist2mem_condition_plot2_notitle(dist_condition_collection, self.conditionname, self.curconditionfilenamelst, bins=bins_, x_lim = max_)
        plt.title(f'{self.conditionname}_condition_MT2isg_distance_density_distribution bin = {bins_}\n  {[i for i in self.curconditionfilenamelst]}')
        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_2_{self.microtube_to_vesicle_condition_plotname}'
        plt.tight_layout()
        plt.savefig(plotname)
        # plt.show()
        plt.close()

class microtube_to_mito_condition(process_organelles_condition):
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'microtube_to_mito_condition'

    def step(self):
        print('conditions', self.conditionlst)

        for conditionname in self.conditionlst:
            print(conditionname)
            self.conditionname = conditionname
            self.obtain_condition_data(conditionname, self.microtube_to_mito_filename)
            # print(1001,conditionname, self.filteredfilelst)
            if len (self.filteredfilelst) != 0:
                self.microtube_to_mito_condition_plot()

        print(f'Done with {self.name}.')


    def microtube_to_mito_condition_plot(self):
        
        dist_condition_collection = []
        distss_in_one_condition = []
        distss_prob_in_one_condition = []
        bins_ = 12

        for dataname in self.filteredfilelst:
            distances = pd.read_csv(f'{dataname}')
            distances = distances['Distance']
            distances_new =distances.tolist()
            print(dataname)
            print('max distance', np.max(distances))
            print('num', len(distances_new))
            dist_condition_collection.extend(distances_new)

        
            disttt, bins = np.histogram(distances,bins=bins_,range=[0,1200])
            distss_in_one_condition.append(disttt)
            distss_prob_in_one_condition.append(disttt / np.sum(disttt))
            print(np.sum(distss_prob_in_one_condition[-1]))
        

        if len(dist_condition_collection) != 0:
            max_= np.ceil(np.max(dist_condition_collection)/100) * 100
        else:
            max_ = 1200
        disttt, _ = np.histogram(dist_condition_collection,bins=bins_,range=[0,max_])

        distss_prob_in_one_condition = [dist / np.sum(disttt) for dist in disttt ]
        
        outplot.dist2mem_condition_plot_notitle(distss_prob_in_one_condition, self.conditionname, self.curconditionfilenamelst, bins=bins_, x_lim = max_)
        plt.title(f'{self.conditionname}_condition_MT to mito mem distance_density_distribution bin = {bins_}\n  {[i for i in self.curconditionfilenamelst]}')
        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_{self.microtube_to_mito_condition_plotname}'
        plt.tight_layout()
        plt.savefig(plotname)       
        plt.show()
        plt.close()

        outplot.dist2mem_condition_plot2_notitle(dist_condition_collection, self.conditionname, self.curconditionfilenamelst, bins=bins_, x_lim = max_)
        plt.title(f'{self.conditionname}_condition_MT to mito mem distance_density_distribution bin = {bins_}\n  {[i for i in self.curconditionfilenamelst]}')
        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_2_{self.microtube_to_mito_condition_plotname}'
        plt.tight_layout()
        plt.savefig(plotname)       
        # plt.show()
        plt.close()

class microtube_to_endoreticulum_condition(process_organelles_condition):
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'microtube_to_endoreticulum_condition'

    def step(self):
        print('conditions', self.conditionlst)

        for conditionname in self.conditionlst:
            print(conditionname)
            self.conditionname = conditionname
            self.obtain_condition_data(conditionname, self.microtube_to_endoreticulum_filename)
            # print(1001,conditionname, self.filteredfilelst)
            if len (self.filteredfilelst) != 0:
                self.microtube_to_endoreticulum_condition_plot()

        print(f'Done with {self.name}.')


    def microtube_to_endoreticulum_condition_plot(self):
        
        dist_condition_collection = []
        distss_in_one_condition = []
        distss_prob_in_one_condition = []
        bins_ = 12

        for dataname in self.filteredfilelst:
            distances = pd.read_csv(f'{dataname}')
            distances = distances['Distance']
            distances_new =distances.tolist()
            print(dataname)
            print('max distance', np.max(distances))
            print('num', len(distances_new))
            dist_condition_collection.extend(distances_new)

        
            disttt, bins = np.histogram(distances,bins=bins_,range=[0,1200])
            distss_in_one_condition.append(disttt)
            distss_prob_in_one_condition.append(disttt / np.sum(disttt))
            print(np.sum(distss_prob_in_one_condition[-1]))
        

        if len(dist_condition_collection) != 0:
            max_= np.ceil(np.max(dist_condition_collection)/100) * 100
        else:
            max_ = 1200
        disttt, _ = np.histogram(dist_condition_collection,bins=bins_,range=[0,max_])

        distss_prob_in_one_condition = [dist / np.sum(disttt) for dist in disttt ]
        

        outplot.dist2mem_condition_plot_notitle(distss_prob_in_one_condition, self.conditionname, self.curconditionfilenamelst, bins=bins_, x_lim = max_)
        plt.title(f'{self.conditionname}_condition_MT to er mem distance_density_distribution bin = {bins_}\n  {[i for i in self.curconditionfilenamelst]}')
        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_{self.microtube_to_endoreticulum_condition_plotname}'
        plt.tight_layout()
        plt.savefig(plotname)            
        plt.show()
        plt.close()

        outplot.dist2mem_condition_plot2_notitle(dist_condition_collection, self.conditionname, self.curconditionfilenamelst, bins=bins_, x_lim = max_)
        plt.title(f'{self.conditionname}_condition_MT to er mem distance_density_distribution bin = {bins_}\n  {[i for i in self.curconditionfilenamelst]}')
        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_2_{self.microtube_to_endoreticulum_condition_plotname}'
        plt.tight_layout()
        plt.savefig(plotname)            
        # plt.show()
        plt.close()

class vesicle_to_mito_condition(process_organelles_condition):
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'vesicle_to_mito_condition'

    def step(self):
        print('conditions', self.conditionlst)

        for conditionname in self.conditionlst:
            print(conditionname)
            self.conditionname = conditionname
            self.obtain_condition_data(conditionname, self.vesicle_to_mito_filename)
            # print(1001,conditionname, self.filteredfilelst)
            if len (self.filteredfilelst) != 0:
                self.vesicle_to_mito_condition_plot()

        print(f'Done with {self.name}.')


    def vesicle_to_mito_condition_plot(self,):

        dist_condition_collection = []
        distss_in_one_condition = []
        distss_prob_in_one_condition = []
        bins_ = 12


        for dataname in self.filteredfilelst:
            distances = pd.read_csv(f'{dataname}')
            distances = distances['Distance']
            distances_new =distances.tolist()
            print(dataname)
            print('max distance', np.max(distances))
            print('num', len(distances_new))
            dist_condition_collection.extend(distances_new)

        
            disttt, bins = np.histogram(distances,bins=bins_)
            distss_in_one_condition.append(disttt)
            distss_prob_in_one_condition.append(disttt / np.sum(disttt))
            print(np.sum(distss_prob_in_one_condition[-1]))
        
        print('len(dist_condition_collection)', len(dist_condition_collection))
        if len(dist_condition_collection) != 0:
            max_= np.ceil(np.max(dist_condition_collection)/100) * 100
        else:
            max_ = 1200
        disttt, _ = np.histogram(dist_condition_collection,bins=bins_,range=[0,max_])

        distss_prob_in_one_condition = [dist / np.sum(disttt) for dist in disttt ]
        
        outplot.dist2mem_condition_plot_notitle(distss_prob_in_one_condition, self.conditionname, self.curconditionfilenamelst, bins=bins_, x_lim = max_)
        plt.title(f'{self.conditionname}_condition_isg2mito_distance_density_distribution bin = {bins_}\n  {[i for i in self.curconditionfilenamelst]}')
        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_{self.vesicle_to_mito_condition_plotname}'
        plt.tight_layout()
        plt.savefig(plotname)        
        plt.show()
        plt.close()


        # outplot.dist2mem_condition_plot2_notitle(distss_prob_in_one_condition, self.conditionname, self.curconditionfilenamelst, bins=bins_, x_lim = max_)
        # plt.title(f'{self.conditionname}_condition_isg2mito_distance_density_distribution bin = {bins_}\n  {[i for i in self.curconditionfilenamelst]}')
        # plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_2_{self.vesicle_to_mito_condition_plotname}'
        # plt.tight_layout()
        # plt.savefig(plotname)            
        # # plt.show()
        # plt.close()

        pass

class vesicle_to_endoreticulum_condition(process_organelles_condition):
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'vesicle_to_endoreticulum_condition'


    def step(self):
        print('conditions', self.conditionlst)

        for conditionname in self.conditionlst:
            print(conditionname)
            self.conditionname = conditionname
            self.obtain_condition_data(conditionname, self.vesicle_to_endoreticulum_filename)
            # print(1001,conditionname, self.filteredfilelst)
            if len (self.filteredfilelst) != 0:
                self.vesicle_to_endoreticulum_condition_plot()

        print(f'Done with {self.name}.')


    def vesicle_to_endoreticulum_condition_plot(self):


        dist_condition_collection = []
        distss_in_one_condition = []
        distss_prob_in_one_condition = []
        bins_ = 12

        for dataname in self.filteredfilelst:
            distances = pd.read_csv(f'{dataname}')
            distances = distances['Distance']
            distances_new =distances.tolist()
            print(dataname)
            print('max distance', np.max(distances))
            print('num', len(distances_new))
            dist_condition_collection.extend(distances_new)

        
            disttt, bins = np.histogram(distances,bins=bins_)
            distss_in_one_condition.append(disttt)
            distss_prob_in_one_condition.append(disttt / np.sum(disttt))
            print(np.sum(distss_prob_in_one_condition[-1]))
        

        if len(dist_condition_collection) != 0:
            max_= np.ceil(np.max(dist_condition_collection)/100) * 100
        else:
            max_ = 1200
        disttt, _ = np.histogram(dist_condition_collection,bins=bins_,range=[0,max_])

        distss_prob_in_one_condition = [dist / np.sum(disttt) for dist in disttt ]
        
        outplot.dist2mem_condition_plot_notitle(distss_prob_in_one_condition, self.conditionname, self.curconditionfilenamelst, bins=bins_, x_lim = max_)
        plt.title(f'{self.conditionname}_condition_isg2er_distance_density_distribution bin = {bins_}\n  {[i for i in self.curconditionfilenamelst]}')
        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_{self.vesicle_to_endoreticulum_condition_plotname}'
        plt.tight_layout()
        plt.savefig(plotname)
        plt.show()
        plt.close()

        # outplot.dist2mem_condition_plot2_notitle(distss_prob_in_one_condition, self.conditionname, self.curconditionfilenamelst, bins=bins_, x_lim = max_)
        # plt.title(f'{self.conditionname}_condition_isg2er_distance_density_distribution bin = {bins_}\n  {[i for i in self.curconditionfilenamelst]}')
        # plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_2_{self.vesicle_to_endoreticulum_condition_plotname}'
        # plt.tight_layout()
        # plt.savefig(plotname)
        # # plt.show()
        # plt.close()

        pass

class mito_to_endoreticulum_condition(process_organelles_condition):
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'mito_to_endoreticulum_condition'

    def step(self):
        print('conditions', self.conditionlst)

        for conditionname in self.conditionlst:
            print(conditionname)
            self.conditionname = conditionname
            self.obtain_condition_data(conditionname, self.mito_to_endoreticulum_filename)
            # print(1001,conditionname, self.filteredfilelst)
            if len (self.filteredfilelst) != 0:
                self.mito_to_endoreticulum_condition_plot()

        print(f'Done with {self.name}.')


    def mito_to_endoreticulum_condition_plot(self):

        dist_condition_collection = []
        distss_in_one_condition = []
        distss_prob_in_one_condition = []

        bins_ = 12


        for dataname in self.filteredfilelst:
            distances = pd.read_csv(f'{dataname}')
            distances = distances['Distance']
            distances_new =distances.tolist()
            print(dataname)
            print('max distance', np.max(distances))
            print('num', len(distances_new))
            dist_condition_collection.extend(distances_new)

        
            disttt, bins = np.histogram(distances,bins=bins_)
            distss_in_one_condition.append(disttt)
            distss_prob_in_one_condition.append(disttt / np.sum(disttt))
            print(np.sum(distss_prob_in_one_condition[-1]))
        

        if len(dist_condition_collection) != 0:
            max_= np.ceil(np.max(dist_condition_collection)/100) * 100
        else:
            max_ = 1200
        disttt, _ = np.histogram(dist_condition_collection,bins=bins_,range=[0,max_])

        distss_prob_in_one_condition = [dist / np.sum(disttt) for dist in disttt ]
        
        outplot.dist2mem_condition_plot_notitle(distss_prob_in_one_condition, self.conditionname, self.curconditionfilenamelst, bins=bins_, x_lim = max_)
        plt.title(f'{self.conditionname}_condition_mito2er_distance_density_distribution bin = {bins_}\n  {[i for i in self.curconditionfilenamelst]}')
        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_{self.mito_to_endoreticulum_condition_plotname}'
        plt.tight_layout()
        plt.savefig(plotname)
        plt.show()
        plt.close()

        # outplot.dist2mem_condition_plot2_notitle(distss_prob_in_one_condition, self.conditionname, self.curconditionfilenamelst, bins=bins_, x_lim = max_)
        # plt.title(f'{self.conditionname}_condition_mito2er_distance_density_distribution bin = {bins_}\n  {[i for i in self.curconditionfilenamelst]}')
        # plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_2_{self.mito_to_endoreticulum_condition_plotname}'
        # plt.tight_layout()
        # plt.savefig(plotname)
        # plt.show()
        # plt.close()

        pass

class actin_angle_distance_pair_condition(process_organelles_condition):
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'actin_angle_distance_pair_condition'

    def step(self):
        print('conditions', self.conditionlst)

        for conditionname in self.conditionlst:
            print(conditionname)
            self.conditionname = conditionname
            self.obtain_condition_data(conditionname, self.actin_angle_distance_pair_filename)
            # print(1001,conditionname, self.filteredfilelst)
            if len (self.filteredfilelst) != 0:
                self.actin_angle_distance_pair_condition_plot()

        print(f'Done with {self.name}.')


    def actin_angle_distance_pair_condition_plot(self):

        dist_condition_collection = []
        distss_in_one_condition = []
        distss_prob_in_one_condition = []
        angle_condition_collection = []

        bins_ = 12

        for dataname in self.filteredfilelst:
            data = pd.read_csv(f'{dataname}')
            distances = data['Distances']
            distances_new =distances.tolist()
            angles = data['Angles']
            angles_new =angles.tolist()

            print(dataname)
            if len(distances_new) != 0:
                print('max distance', np.max(distances_new))
                print('angle range', np.min(angles_new), np.max(angles_new))
                print('num', len(distances_new))


            dist_condition_collection.extend(distances_new)
            angle_condition_collection.extend(angles_new)
        
        # if len(dist_condition_collection) == 0:
        #     dist_condition_collection = [0]
        #     angle_condition_collection = [0]

        if len(dist_condition_collection) != 0:
            max_= np.ceil(np.max(dist_condition_collection)/100) * 100
        else:
            max_ = 1200

        disttt, _ = np.histogram(dist_condition_collection,bins=bins_,range=[0,max_])
        distss_prob_in_one_condition = [dist / np.sum(disttt) for dist in disttt ]
        
        print('num', len(dist_condition_collection))
        # print('max', np.max(dist_condition_collection))
        outplot.vis_actindistangle_conditionplot(dist_condition_collection, angle_condition_collection,)
        plt.title(f'{self.conditionname} condition filament angle2dist distribution  \n {[i for i in self.curconditionfilenamelst]}' )
        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_{self.actin_angle_distance_pair_condition_filename}'
        plt.tight_layout()
        plt.clim(vmin = 0, vmax = 2500)
        plt.savefig(plotname)
        plt.show()
        plt.close()

        outplot.vis_actindistangle_conditionplot2(dist_condition_collection, angle_condition_collection,)
        plt.title(f'{self.conditionname} condition filament angle2dist distribution  \n {[i for i in self.curconditionfilenamelst]}' )
        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_2_{self.actin_angle_distance_pair_condition_filename}'
        plt.tight_layout()
        plt.clim(vmin = 0, vmax = 250)
        plt.savefig(plotname)
        plt.show()
        plt.close()


class focal_to_focal_condition(process_organelles_condition):
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'focal_to_focal'
        
    def step(self):
        print('conditions', self.conditionlst)
        for conditionname in self.conditionlst:
            print('focal_focal conditionname:', conditionname)
            self.conditionname = conditionname
            self.obtain_condition_data(conditionname, self.focal_to_focal_distance_filename)
            print(1001,conditionname, self.filteredfilelst)
            if len (self.filteredfilelst) != 0:                
                self.focal_to_focal_condition_plot()
        print(f'Done with {self.name}.')

    def focal_to_focal_condition_plot(self):
        merged_fa_distance =  []
        for filename in self.filteredfilelst:
            dist3 = pd.read_csv(filename)
            distance =  list(dist3['Distance'])
            merged_fa_distance.extend(distance)
            # print(len(angle1), len(merged_data_angle1))
        merged_fa_distance_pd = pd.DataFrame({'Distance':merged_fa_distance})

        outplot.dist_focal2focal_condition_hist(merged_fa_distance_pd, self.conditionname, self.curconditionfilenamelst)
        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_{self.focal_to_focal_distance_condition_plotname}'
        plt.xlim([-500,2500])
        plt.ylim([0, 0.0012])
        plt.tight_layout()
        plt.savefig(plotname)
        plt.show()
        plt.close()

class focal_aggregated_actin_angle_condition(process_organelles_condition):
    '''
    Actins aggregated within one focal adhesion are regareded as one cluster.
    calculate the angle distribution of the actin end to focal adhesion angle.
    distributions.
    
    '''
    def __init__(self, dataid):

        super().__init__(dataid)
        self.name = f'focal_aggregated_actin_angle'
        self.dist_threshold = 200 # nm
    def step(self):
        print('conditions', self.conditionlst)
        for conditionname in self.conditionlst:
            print(conditionname)
            self.conditionname = conditionname
            self.obtain_condition_data(conditionname, self.focal_aggregated_actin_angle_map_filename)
            # print(1001,conditionname, self.filteredfilelst)
            if len (self.filteredfilelst) != 0:
                self.focal_aggregated_actin_angle_condition_plot()

        print(f'Done with {self.name}.')


    def focal_aggregated_actin_angle_condition_plot(self):
        merged_data_angle =  []
        merged_data_distance =  []
        for filename in self.filteredfilelst:
            angle_dist_pairs = pd.read_csv(filename)
            angle, distance =  list(angle_dist_pairs['Angle']), list(angle_dist_pairs['Distance'])
            merged_data_angle.extend(angle)
            merged_data_distance.extend(distance)
            # print(len(angle1), len(merged_data_angle1))

        angle_angle_pairs_pd = pd.DataFrame({'Angle':merged_data_angle, 'Distance':merged_data_distance})

        outplot.fa_actin_angle_dist_map_condition_plot(angle_angle_pairs_pd, self.conditionname, self.curconditionfilenamelst)
        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_{self.focal_aggregated_actin_angle_map_condition_plotname}'
        plt.tight_layout()
        plt.savefig(plotname)
        plt.show()
        plt.close()

class focal_aggregated_actin_angles_pair_condition(process_organelles_condition):
    '''
    Actins aggregated within one focal adhesion are regareded as one cluster.
    calculate the angle distribution of the actin to other actin and the angle to mem
    output: 2d plot angle to other angle Vs angle to mem
    '''

    def __init__(self, dataid):
        super().__init__(dataid)
        self.dataid = None
        self.name = f'focal_aggregated_actin_angles_pair_condition'
        self.dist_threshold = 200 # nm
    def step(self):
        
        print('conditions', self.conditionlst)
        for conditionname in self.conditionlst:
            print(conditionname)
            self.conditionname = conditionname
            self.obtain_condition_data(conditionname, self.focal_aggregated_actin_angles_pair_filename)
            
            # print(1001,conditionname, self.filteredfilelst)
            if len (self.filteredfilelst) != 0:
                self.focal_aggregated_actin_angles_pair_condition_plot()


        print(f'Done with {self.name}.')



    def focal_aggregated_actin_angles_pair_condition_plot(self):

        merged_data_angle1 =  []
        merged_data_angle2 =  []
        for filename in self.filteredfilelst:
            angle_angle_pairs = pd.read_csv(filename)
            angle1, angle2 =  list(angle_angle_pairs['Angle actin']), list(angle_angle_pairs['Angle mem'])
            merged_data_angle1.extend(angle1)
            merged_data_angle2.extend(angle2)

            # print(len(angle1), len(merged_data_angle1))

        angle_angle_pairs_pd = pd.DataFrame({'Angle actin':merged_data_angle1, 'Angle mem':merged_data_angle2})

        
        outplot.fa_actin_angle_angle_pair_condition_plot(angle_angle_pairs_pd, self.conditionname,self.curconditionfilenamelst)
        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_{self.focal_aggregated_actin_angles_pair_condition_plotname}'
        plt.tight_layout()
        plt.savefig(plotname)
        plt.show()
        plt.close()
     


class focal_to_ca_condition(process_organelles_condition):
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'focal_to_ca'
        
    def step(self):
        print('conditions', self.conditionlst)
        for conditionname in self.conditionlst:
            self.conditionname = conditionname
            self.obtain_condition_data(conditionname, self.focal_to_ca_distance_filename)
            # print(1001,conditionname, self.filteredfilelst)
            if len (self.filteredfilelst) != 0:
                self.focal_to_ca_condition_plot()
        print(f'Done with {self.name}.')


    def focal_to_ca_condition_plot(self):
        merged_ca_distance =  []
        for filename in self.filteredfilelst:
            dist3 = pd.read_csv(filename)
            distance =  list(dist3['Distance'])
            merged_ca_distance.extend(distance)
            # print(len(angle1), len(merged_data_angle1))
        merged_ca_distance_pd = pd.DataFrame({'Distance': merged_ca_distance})
        merged_ca_distance_pd.to_csv('D:/filament_project_weimin/data/220814/results/merged_ca_distance_pd_test.csv',index=False)

        outplot.dist_focal2ca_condition_hist(merged_ca_distance_pd, self.conditionname, self.curconditionfilenamelst)
        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_{self.focal_to_ca_distance_condition_plotname}'
        plt.xlim([-500,2500])
        plt.ylim([0, 0.0012])
        plt.tight_layout()
        plt.savefig(plotname)
        plt.show()
        plt.close()

#%%

def main():
    '''
    get dataid
    get catagories for calculation
    get list if catagory available
    check if filtered data exist
    generate figure 
    '''
    dataidlst = preprocess.read_dataid()
    print(dataidlst)
    dataid = dataidlst[0]

    # dataid = '20210517_30_010' 
    # # vesicle actin microtube mito endo_reticulum
    # # print(catagory)

    # for dataid in dataidlst:

    if True:
        a2a = actin_to_actin_condition(dataid)
        a2a.step()  
        # retu

        # a2v = actin_to_vesicle_condition(dataid)
        # a2v.step()
        # # bin = 50    / 100

        # mt2a = microtube_to_actin_condition(dataid)
        # mt2a.step()
        # # dist distribution

        # # a2mi = actin_to_mito_condition(dataid) 
        # # a2mi.step()

        # # a2e = actin_to_endoreticulum_condition(dataid)
        # # a2e.step()

        # mt2i = microtube_to_vesicle_condition(dataid)
        # mt2i.step()
        # # bin = 50 / 100

        # mt2mi = microtube_to_mito_condition(dataid)
        # mt2mi.step()

        # mt2er = microtube_to_endoreticulum_condition(dataid)
        # mt2er.step()

        # i2mi =vesicle_to_mito_condition(dataid)
        # i2mi.step()

        # # isg - mito  distance (Y) / isg - mt distance (x)

        # i2er = vesicle_to_endoreticulum_condition(dataid)
        # i2er.step()

        # mi2er = mito_to_endoreticulum_condition(dataid)
        # mi2er.step()

        # acta2dp = actin_angle_distance_pair_condition(dataid)
        # acta2dp.step()

        # fa2fa = focal_to_focal_condition(dataid)
        # fa2fa.step()

        # faaa = focal_aggregated_actin_angle_condition(dataid)
        # faaa.step()

        # faaap = focal_aggregated_actin_angles_pair_condition(dataid)
        # faaap.step()

        # fc = focal_to_ca_condition(dataid)
        # fc.step()

# plt.clim(0, 250) #colorbar_weimin, set the lim scale of the colorbar

if __name__ == '__main__':
    main()
    # single_data()
    print('Done.')

# %%
1
# %%
