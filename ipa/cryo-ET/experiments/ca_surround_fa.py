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

import seaborn as sns
print(torch.__version__)
print(torch.cuda.is_available())

device = torch.device('cuda')





# %%

def getsingledatafilename():
    filename =  'ca_surround_fa_num.csv' 
    plotname =  'ca_surround_fa_num.png'

    return filename, plotname
def ca_surround_fa_num_generator(FAcoords, cacoords, voxel_size_xyz, cutoff):
    
    '''input cutoff in nm'''
    cutoff_invoxel = cutoff / np.average(voxel_size_xyz) * 10
    avesize = np.average(voxel_size_xyz)

    Dists = cdist(FAcoords, cacoords, metric='euclidean')
    ca_surround_fa_lst = list()
    for ii, singleca in enumerate(cacoords):
        cur_ca_curround_fa_num = 0
        temp_ca_surround_fa_lst =[]
        for jj, singlefa in enumerate(FAcoords):
            dist = np.linalg.norm(np.array(singleca) - np.array(singlefa))
            if dist <= cutoff_invoxel:
                cur_ca_curround_fa_num += 1
                temp_ca_surround_fa_lst.append(singlefa)
        ca_surround_fa_lst.append(cur_ca_curround_fa_num)

        if ca_surround_fa_lst[-1] > 10:
            print('cur ca coord:', singleca)
            print(cur_ca_curround_fa_num)
            print('cur fa', temp_ca_surround_fa_lst)

    return ca_surround_fa_lst




class ca_surround_fa(process_organelles):
    def __init__(self, dataid) -> None:
        super().__init__(dataid)
        self.name = f'ca_surround_fa_num'
        self.ca_surround_fa_num_filename , self.ca_surround_fa_num_plotname = getsingledatafilename()

    def step(self):
        print(f'dataid {self.dataid}')
        organelletypelst = [self.ca_type, self.focal_adhesion_type]
        print(organelletypelst)
        if not self.check_organelle_type(organelletypelst):
            print(f'{self.dataid} {self.name} data not exist.')

        else:
            if not self.check_data_updated(organelletypelst):
                self.generate_processed_data()
                # print('line1559')
            self.ca_surround_fa_num_plot(overlap = self.plotdata_overlap)
        print(f'Done with {self.name}.')


    def ca_surround_fa_num_plot(self, overlap):
        
        self.ca_surround_fa_num_filepath = f'{self.datanewpath}/{self.dataid}_{self.ca_surround_fa_num_filename}'
        self.curdatapath = self.ca_surround_fa_num_filepath
        if not os.path.exists(self.curdatapath) or overlap == True:
            self.ca_surround_fa_num_distance()

        ca_data = pd.read_csv(self.cafilepath)
        fa_data = pd.read_csv(self.focaladhesionfilepath)

        num3 = pd.read_csv(self.ca_surround_fa_num_filepath)
        print(np.array(num3).reshape(1,-1))

        sns.distplot(a=num3, kde=True)
        # sns.kdeplot(data=dist, shade=True)
        plt.title(f'{self.dataid} ca_surround_fa_num distribution \n datanum {len(num3)} fa num {len(fa_data)} ca num {len(ca_data)}' )
        # plt.ylabel('probability density')
        # plt.xlabel('focal - calcium distance (nm)')  


        plotname = f'{self.singledataoutplotpath}/{self.location}_{self.condition}_{self.dataid}_{self.ca_surround_fa_num_plotname}'
        # plt.savefig(plotname)
        plt.show()
        plt.close()
     

    def ca_surround_fa_num_distance(self):

        shift_xml_xyz, voxel_size_xyz = preprocess.load_parameters(self.dataid)
        cutoff = 100 # nm
        fa_data = pd.read_csv(self.focaladhesionfilepath)
        fa_data = [[fa_data['X'][i], fa_data['Y'][i], fa_data['Z'][i] ] for i in range(fa_data.shape[0]) ]
        ca_data = pd.read_csv(self.cafilepath)
        ca_data  = [[ca_data['X'][i], ca_data['Y'][i], ca_data['Z'][i] ] for i in range(ca_data.shape[0]) ]


        Fa_data_new = fa_data
        ca_surround_fa_num = ca_surround_fa_num_generator(Fa_data_new, ca_data, voxel_size_xyz, cutoff)  # nm
        numpd = pd.DataFrame(columns=['Number'], data= ca_surround_fa_num)
        
        # filament2filament_dist_lst_savepd.loc['distance'] = filament2filament_dist_lst
        numpd.to_csv(f'{self.datanewpath}/{self.dataid}_{self.ca_surround_fa_num_filename}', index=False)




def main_single():

    dataidlst = preprocess.read_dataid()
    print(dataidlst)
    # dataidlst = dataidlst[:1]

    # dataidlst = ['20220325_2.8_20']
    # dataidlst = ['20210323_2.8_011']

    for dataid in dataidlst:
        print(dataid)

        casfa = ca_surround_fa(dataid)
        casfa.plotdata_overlap = True
        casfa.step()






#%%  conditions


from process_organelle_data_condition import process_organelles_condition



class ca_surround_fa_condition(process_organelles_condition):
    def __init__(self, dataid):
        super().__init__(dataid)
        self.name = f'ca_surround_fa_condition'
        self.ca_surround_fa_num_filename, _ = getsingledatafilename()
        self.ca_surround_fa_num_condition_filename = 'ca_surround_fa_num_condition.png'
    def step(self):
        print('conditions', self.conditionlst)
        for conditionname in self.conditionlst:
            print(conditionname)
            self.conditionname = conditionname
            self.obtain_condition_data(conditionname, self.ca_surround_fa_num_filename)
            # print(1001,conditionname, self.filteredfilelst)
            if len (self.filteredfilelst) != 0:
                self.ca_surround_fa_condition_plot()

        print(f'Done with {self.name}.')



    def ca_surround_fa_condition_plot(self):


        merged_fa_distance =  []

        for filename in self.filteredfilelst:
            dist3 = pd.read_csv(filename)

            distance =  list(dist3['Number'])

            merged_fa_distance.extend(distance)

            # print(len(angle1), len(merged_data_angle1))


        plt.hist(merged_fa_distance, bins=20,range=[0,20], density=True)
        plt.title(f'{self.location}_{self.conditionname}_ca_surround_fa_num_condition \n num{len(merged_fa_distance)}')
        plt.xlabel('Length')
        plt.ylabel('Probability')
        plt.xlim([0,20])
        plt.ylim([0,1])

        plotname = f'{self.conditiondataoutplotpath}/{self.location}_{self.conditionname}_{self.ca_surround_fa_num_condition_filename}'
        plt.tight_layout()
        # plt.savefig(plotname)
        plt.show()
        plt.close()









def main_condition():

    dataidlst = preprocess.read_dataid()
    dataid = '1'

    casfa = ca_surround_fa_condition(dataid)
    casfa.plotdata_overlap = True
    casfa.step()

if __name__ == '__main__':
    # main_single()
    main_condition()


#%%








