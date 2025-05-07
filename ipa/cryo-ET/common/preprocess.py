from matplotlib.cbook import maxdict
import numpy as np
import scipy
import scipy.ndimage
import skimage.feature
import skimage.morphology
import pandas as pd
import matplotlib.pyplot as plt
import mrcfile, tifffile
import xml
import xml.etree.ElementTree as ET
import json
import os, sys
import copy
from common.parser import arg
from scipy.spatial.distance import cdist

def read_dataid():

    idlst = pd.read_csv(os.path.join(arg.root_dir, arg.parameter_dir, 'dataid.csv'))

    idlst = [ idlst['ID'].iloc[i] for i in range(0,idlst.shape[0])  ]
    # print(idlst)

    return idlst # ['20210517_30_011', '20210517_30_012',]


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
            subdir3.append( os.path.join(dir, singlereplica) )
    # print(subdir3)

    return subdir2, subdir3  # subdir2 condition location   subdir3 data location 



def load_files(data_dir, check_ = False):

    mrcname = ''
    mrcname_fullfill = ''
    filamentdata_xmlname = ''
    filamentdata_jsonname_full = ''
    mtdata_xmlname = ''
    MT_jsonname_full = ''
    mito_mrcname = ''
    er_mrcname = ''
    for maindir, subdir, file_name_list in os.walk(data_dir, topdown=False):
        filelist = np.array(file_name_list)
    # print(filelist)

    for name in filelist:
        if '.mrc' in name and 'mito' not in name.lower() and 'er' not in name.lower():
            if 'fullfill' not in name:
                mrcname = os.path.join(data_dir, name)   
            elif 'fullfill' in name:
                mrcname_fullfill = os.path.join(data_dir, name)
            # file_idx = mrcname.split('.')[0]       
        elif '.xml' in name and ('filament' in name.lower() or 'actin' in name.lower()):
            if 'fill' not in name:
                filamentdata_xmlname = os.path.join(data_dir, name)
        elif '.json' in name and 'filament' in name and 'filled' in name:
                filamentdata_jsonname_full = os.path.join(data_dir, name)
        
        elif '.xml' in name and ('mt' in name.lower() or 'microtube' in name.lower()):
            if 'fill' not in name:
                mtdata_xmlname = os.path.join(data_dir, name)
        elif '.json' in name and 'mt' in name.lower() and 'filled' in name:
                MT_jsonname_full = os.path.join(data_dir, name)
        elif '.mrc' in name and 'mito_filled' in name.lower():
                mito_mrcname = os.path.join(data_dir, name)
        elif '.mrc' in name and 'er' in name.lower():
                er_mrcname = os.path.join(data_dir, name)
        else:
            pass

    if check_:
        if os.path.exists(mrcname): print('isg mrcname:',mrcname) #'mrcname_fullfill', mrcname_fullfill)
        if os.path.exists(mrcname_fullfill): print('isg mrcname processed:', mrcname_fullfill)
        if os.path.exists(filamentdata_xmlname): print('actin xmlname:', filamentdata_xmlname)
        if os.path.exists(filamentdata_jsonname_full): print('actin xml processed:', filamentdata_jsonname_full)
        if os.path.exists(mtdata_xmlname): print('MT name', mtdata_xmlname)
        if os.path.exists(MT_jsonname_full): print('MT processed:', MT_jsonname_full)
        if os.path.exists(mito_mrcname): print('mito data:', mito_mrcname)
        if os.path.exists(er_mrcname): print('er data:', er_mrcname)


    filenamelst = [mrcname, mrcname_fullfill, filamentdata_xmlname, 
                    filamentdata_jsonname_full, 
                    mtdata_xmlname, MT_jsonname_full,
                    mito_mrcname,
                    er_mrcname]

    

    return filenamelst



def get_filelist(dir, Filelist):
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
            get_filelist(newDir, Filelist)
    return Filelist

def load_focaladhesion_data(data_dir, dataid, check_ = False):
    '''
    path F:\\salilab\\projects\\filament_project_weimin\\folcon\\FA
    '''

    fa_xmlname = ''

    filelist = get_filelist(data_dir, [])
    filelist = [filename for filename in filelist if dataid in filename]

    return filelist


def load_ca_data(data_dir, dataid, check_ = False):
    '''
    path F:\\salilab\\projects\\filament_project_weimin\\folcon\\ca
    '''

    return load_focaladhesion_data(data_dir, dataid,check_)




def load_parameters(file_idx):
    parpath_ = os.path.join(arg.root_dir, arg.parameter_dir)
    df = pd.read_excel(f'{parpath_}/Mrc_offset.xlsx')
    shift_data=df.values
    # print(shift_data)

    shift_xml_xyz = [0,0,0]
    testn = 0
    for line in shift_data:
        if file_idx in line:
            testn = 1
            shift_xml_xyz[0] = line[3]
            shift_xml_xyz[1] = line[2]
            shift_xml_xyz[2] = line[1]

    assert testn == 1
    # print('shift_xml_xyz:', shift_xml_xyz)

    ## voxelsize 

    df2 = pd.read_excel(f'{parpath_}/Voxel_Size.xlsx')
    voxel_data=df2.values
    voxel_size_xyz = [0,0,0]
    for line in voxel_data:
        if file_idx in line:
            voxel_size_xyz[0] = line[1] * 10
            voxel_size_xyz[1] = line[2] * 10
            voxel_size_xyz[2] = line[3] * 10
            
    # print('voxel size:', voxel_size_xyz)

    return shift_xml_xyz, voxel_size_xyz



# def get_filelist(dir, Filelist):
#     newDir = dir
#     if os.path.isfile(dir):
#         Filelist.append(dir)
#         # # 若只是要返回文件文，使用这个
#         # Filelist.append(os.path.basename(dir))
#     elif os.path.isdir(dir):
#         for s in os.listdir(dir):
#             # 如果需要忽略某些文件夹，使用以下代码
#             #if s == "xxx":
#                 #continue
#             newDir=os.path.join(dir,s)
#             get_filelist(newDir, Filelist)

#     return Filelist


# lists = get_filelist(os.path.join(arg.root_dir, arg.data_root_dir), [])
# actin_file_lst = [file for file in lists if 'actin.xml' in file.lower()]





def fill_actin_point_gap(filename, resolution = 40):
    '''
    fill pints gap to generate representative points for filaments
    reorganize coordinates zyx to xyz
    
    '''
    print(f'resolution: { resolution} A')

    xmlname = f'{filename}'
    f = open(xmlname, encoding='utf-8') 
    xml_txt = f.read()
    root = ET.fromstring(xml_txt)

    # tree = ET.parse(xmlname)
    # root = tree.getroot() # Workbook
    # tag = Element.tag
    # attrib = Element.attrib
    # value = Element.text

    pagenamelst = ['nodes', 'points', 'segments']
    pageidxlst = []
    for ii, pagename in enumerate(pagenamelst):
        for jj, child in enumerate(root):
            # print(child)
            # print(child.tag,":", child.attrib)
            if pagename in str(child.attrib).lower():
                pageidxlst.append(jj)
                # print(pagename)

    assert len(pageidxlst) == 3

    Nodes = root[int(pageidxlst[0])]
    Points = root[int(pageidxlst[1])]
    Segments = root[int(pageidxlst[2])]



    # for points 
    # df.loc[0]=['cat', 3] 
    # Points = root[2]
    # for tabblee in Points:
    #     print(tabblee.tag)
    # print(list(Points)[0].tag)
    assert 'table' in str(list(Points)[0].tag).lower()
    tablee = list(Points)[0]

    rowws = list(tablee) #tablee.getchildren()
    # roww_eg = tablee.getchildren()[0:20]
    # rowws = roww_eg
    # set column name 

    start_n = 0
    for row_i, roww in enumerate(rowws):
        if 'row' not in str(roww.tag).lower():
            start_n = row_i+1
        else: break

    # print(start_n)



    columnname = []
    for roww in rowws[:start_n + 1]:
        # print(roww.tag)
        if 'row' in str(roww.tag).lower():
            cells = list(roww)  #roww.getchildren()

            for cell in cells:
                datas = list(cell)  #cell.getchildren()
                for data in datas:
                    columnname.append(data.text)
    # print('103', columnname)

    Points_pd =pd.DataFrame(columns=columnname)
    # print(Points_pd)


    # obtain data
    for i, roww in enumerate(rowws[start_n + 1:]):
        cells = roww.getchildren()
        templst = []
        for cell in cells:
            datas = cell.getchildren()
            for data in datas:
                templst.append(data.text)
        Points_pd.loc[i]= templst

    print(Points_pd.head(5))


    Points_pd_x = Points_pd['X Coord']
    Points_pd_y = Points_pd['Y Coord']
    Points_pd_z = Points_pd['Z Coord']

    # print(Points_pd_x.head(5), Points_pd_y.head(5), Points_pd_z.head(5))


    #
    # for Segments 
    # df.loc[0]=['cat', 3] 

    # Segments = root[3]
    assert 'table' in str(list(Points)[0].tag).lower()

    tablee = Segments.getchildren()[0]
    rowws = tablee.getchildren()
    # roww_eg = tablee.getchildren()[0:20]
    # rowws = roww_eg
    # set column name 

    start_n = 0
    for row_i, roww in enumerate(rowws):
        if 'row' not in str(roww.tag).lower():
            start_n = row_i+1
        else: break


    for roww in rowws[:start_n+1]:
        # print(roww.tag)
        if 'row' in str(roww.tag).lower():
            cells = roww.getchildren()
            columnname = []
            for cell in cells:
                datas = cell.getchildren()
                for data in datas:
                    columnname.append(data.text)
    # print('140', columnname)

    segments_pd =pd.DataFrame(columns=columnname)

    # segments_pd_points =segments_pd['Point IDs']


    # print(Points_pd.head(5))

    # obtain data
    for i, roww in enumerate(rowws[start_n+1:]):
        cells = roww.getchildren()
        templst = []
        for cell in cells:
            datas = cell.getchildren()
            for data in datas:
                templst.append(data.text)
        segments_pd.loc[i]= templst

    print(segments_pd.head(5))


    # for test
    # segments_pd = segments_pd[:10]


    points = Points_pd

    #
    # print(points.iloc[1])
    coords_lst = []
    for idx in range(len(points)):
        coords_lst.append([float(Points_pd_z.iloc[idx]), float(Points_pd_y.iloc[idx]), float(Points_pd_x.iloc[idx])])  #zyx from csv

    coords = np.array(coords_lst).reshape(-1,3)
    # coords = coords - np.array([shift_xml_xyz[2], shift_xml_xyz[1], shift_xml_xyz[0]])
    # print(coords)  # [x,y,z]



    #  
    # input relations 
    filament_idx_lst = []
    filament_coord_lst = []
    segments = segments_pd
    # print(segments.columns)
    # ['Segment ID', 'CurvedLength', 'MeanRadius', 'Volume',
    #        'OrientationTheta', 'OrientationPhi', 'SubgraphID', 'ChordLength',
    #        'Tortuosity', 'TensorXX', 'TensorYY', 'TensorZZ', 'TensorXY',
    #        'TensorXZ', 'TensorYZ', 'Node ID #1', 'Node ID #2', 'Point IDs']


    # print(segments['Segment ID'])
    # print(segments['Point IDs'])
    correspond_id = segments['Segment ID']
    correspond_points_id = segments['Point IDs']

    #   angle =  ['OrientationTheta']  # angle to z axis 

    # print(correspond_id)

    for idx in range(len(segments)):
    # for idx in range(5):
        assert idx ==  int(correspond_id[idx]) ## match 
        coord_idx_lst = [ int(pointidx) for pointidx in correspond_points_id[idx].split(',') ]# str
        # print(coord_idx_lst)
        # print(len(coord_idx_lst))
        temp_coordlst = []
        for pointidx in coord_idx_lst:
            temp_coordlst.append(coords[pointidx])  #   coords[pointidx] unit A  , voxel_coords

        filament_idx_lst.append(coord_idx_lst)
        filament_coord_lst.append(temp_coordlst)

    # print(len(filament_coord_lst))
    # print(len(filament_coord_lst[0]))



    # fill gap

    
    def get_line(point1, point2):
        # point1 to point2
        # resolution 40=4nm
        vect = np.array([point2[0] - point1[0], point2[1] - point1[1], point2[2] - point1[2]])
        unit_vect = vect / np.linalg.norm(vect) * resolution  # 40 = 4 nm
        thresh_len = resolution * 0.2
        points_sum = np.linalg.norm(vect) // resolution + 1
        pointslst = [ point1 + unit_vect * i for i in range(int(points_sum))]
        while True:
            if np.linalg.norm(pointslst[-1] - point1) > np.linalg.norm(vect):
                pointslst = pointslst[:-1]
            else:
                break
        # if np.linalg.norm(pointslst[-1] - np.array(point2)) < thresh_len:
        #     pointslst = pointslst[:-1]

        # print(unit_vect, pointslst)
        # voxel_loc = np.around(pointslst)
        # print(voxel_loc)
        return pointslst


    # extend coord to narrow interval
    filament_coord_extend_lst = [] 

    for filament_coord in filament_coord_lst:
        single_filament_coord_extend_lst = []
        for i in range(1, len(filament_coord)):
            coord1, coord2 = filament_coord[i-1], filament_coord[i]  # array
            temp_coord_extend_lst = get_line(coord1, coord2)
            # print('temp_coord_extend_lst', len(temp_coord_extend_lst))
            single_filament_coord_extend_lst.extend(temp_coord_extend_lst)
        single_filament_coord_extend_lst.append(filament_coord[-1])
        # print(len(single_filament_coord_extend_lst))
        filament_coord_extend_lst.append(single_filament_coord_extend_lst)
            # print(coord1, coord2)
            # print(filament_coord_extend_lst)  

    # filament_coord_lst
    coords_intotal1 = 0
    for single_filament_coord_extend in filament_coord_lst:
        coords_intotal1 += len(single_filament_coord_extend)
    print(f'coord intotal before: {coords_intotal1}')
    coords_intotal = 0
    for single_filament_coord_extend in filament_coord_extend_lst:
        coords_intotal += len(single_filament_coord_extend)
    print(f'coord intotal after: {coords_intotal}')


    ## save file
    img_columnname = [f'filament_{i}' for i in range(len(filament_coord_extend_lst))]
    # print(img_columnname)
    # df = pd.DataFrame(columns=[img_columnname])
    # print(df)
    points_all_dict = {}

    for idx, filament_vects in enumerate(filament_coord_extend_lst):
        points_all_dict[f'{img_columnname[idx]}'] =[]
        for coord in filament_vects:
            points_all_dict[f'{img_columnname[idx]}'] .append(list(coord))
            # df.loc[f'{img_columnname[idx]}'] = temp_vectlst

    # print(vect_all_dict)
    #%%
    # print(vect_all_dict.keys())

    return points_all_dict




def get_actin_to_pm_angle(filename):
    '''
    obtain actin angle to plasma membrane from segmentation page
    '''


    xmlname = f'{filename}'
    f = open(xmlname, encoding='utf-8') 
    xml_txt = f.read()
    root = ET.fromstring(xml_txt)

    # tree = ET.parse(xmlname)
    # root = tree.getroot() # Workbook
    # tag = Element.tag
    # attrib = Element.attrib
    # value = Element.text

    pagenamelst = ['nodes', 'points', 'segments']
    pageidxlst = []
    for ii, pagename in enumerate(pagenamelst):
        for jj, child in enumerate(root):
            # print(child)
            # print(child.tag,":", child.attrib)
            if pagename in str(child.attrib).lower():
                pageidxlst.append(jj)
                # print(pagename)

    assert len(pageidxlst) == 3

    Nodes = root[int(pageidxlst[0])]
    Points = root[int(pageidxlst[1])]
    Segments = root[int(pageidxlst[2])]



    # for Segments 
    # df.loc[0]=['cat', 3] 

    # Segments = root[3]
    assert 'table' in str(list(Points)[0].tag).lower()

    tablee = Segments.getchildren()[0]
    rowws = tablee.getchildren()
    # roww_eg = tablee.getchildren()[0:20]
    # rowws = roww_eg
    # set column name 

    start_n = 0
    for row_i, roww in enumerate(rowws):
        if 'row' not in str(roww.tag).lower():
            start_n = row_i+1
        else: break


    for roww in rowws[:start_n+1]:
        # print(roww.tag)
        if 'row' in str(roww.tag).lower():
            cells = roww.getchildren()
            columnname = []
            for cell in cells:
                datas = cell.getchildren()
                for data in datas:
                    columnname.append(data.text)
    # print('140', columnname)

    segments_pd =pd.DataFrame(columns=columnname)

    # segments_pd_points =segments_pd['Point IDs']


    # print(Points_pd.head(5))

    # obtain data
    for i, roww in enumerate(rowws[start_n+1:]):
        cells = roww.getchildren()
        templst = []
        for cell in cells:
            datas = cell.getchildren()
            for data in datas:
                templst.append(data.text)
        segments_pd.loc[i]= templst

    print(segments_pd.head(5))



    # input relations 
    filament_idx_lst = []
    filament_coord_lst = []
    segments = segments_pd
    # print(segments.columns)
    # ['Segment ID', 'CurvedLength', 'MeanRadius', 'Volume',
    #        'OrientationTheta', 'OrientationPhi', 'SubgraphID', 'ChordLength',
    #        'Tortuosity', 'TensorXX', 'TensorYY', 'TensorZZ', 'TensorXY',
    #        'TensorXZ', 'TensorYZ', 'Node ID #1', 'Node ID #2', 'Point IDs']


    # print(segments['Segment ID'])
    # print(segments['Point IDs'])
    correspond_id = segments['Segment ID']
    correspond_points_id = segments['Point IDs']
    correspond_angle = segments['OrientationTheta']


    ## save file
    img_columnname = [f'filament_{i}' for i in range(correspond_angle.shape[0])]

    # print(img_columnname)
    # df = pd.DataFrame(columns=[img_columnname])
    # print(df)
    angle_all_dict = {}

    for idx, name in enumerate(img_columnname):
        angle_all_dict[f'{name}'] = correspond_angle.iloc[idx]


    return angle_all_dict



def fill_microtube_point_gap(filename,  resolution = 40):
    '''
    fill pints gap to generate representative points for filaments
    reorganize coordinates zyx to xyz
    
    '''
    print(f'resolution: { resolution} A')

    xmlname = filename
    f = open(xmlname, encoding='utf-8') 
    xml_txt = f.read()
    root = ET.fromstring(xml_txt)

    # tree = ET.parse(xmlname)
    # root = tree.getroot() # Workbook
    # tag = Element.tag
    # attrib = Element.attrib
    # value = Element.text

    pagenamelst = ['nodes', 'points', 'segments']
    pageidxlst = []
    for ii, pagename in enumerate(pagenamelst):
        for jj, child in enumerate(root):
            # print(child)
            # print(child.tag,":", child.attrib)
            if pagename in str(child.attrib).lower():
                pageidxlst.append(jj)
                # print(pagename)

    assert len(pageidxlst) == 3

    Nodes = root[int(pageidxlst[0])]
    Points = root[int(pageidxlst[1])]
    Segments = root[int(pageidxlst[2])]



    # for points 
    # df.loc[0]=['cat', 3] 
    # Points = root[2]
    # for tabblee in Points:
    #     print(tabblee.tag)
    # print(list(Points)[0].tag)
    assert 'table' in str(list(Points)[0].tag).lower()
    tablee = list(Points)[0]

    rowws = list(tablee) #tablee.getchildren()
    # roww_eg = tablee.getchildren()[0:20]
    # rowws = roww_eg
    # set column name 

    start_n = 0
    for row_i, roww in enumerate(rowws):
        if 'row' not in str(roww.tag).lower():
            start_n = row_i+1
        else: break

    # print(start_n)



    columnname = []
    for roww in rowws[:start_n + 1]:
        # print(roww.tag)
        if 'row' in str(roww.tag).lower():
            cells = list(roww)  #roww.getchildren()

            for cell in cells:
                datas = list(cell)  #cell.getchildren()
                for data in datas:
                    columnname.append(data.text)
    # print('103', columnname)

    Points_pd =pd.DataFrame(columns=columnname)
    # print(Points_pd)


    # obtain data
    for i, roww in enumerate(rowws[start_n + 1:]):
        cells = roww.getchildren()
        templst = []
        for cell in cells:
            datas = cell.getchildren()
            for data in datas:
                templst.append(data.text)
        Points_pd.loc[i]= templst

    print(Points_pd.head(5))


    Points_pd_x = Points_pd['X Coord']
    Points_pd_y = Points_pd['Y Coord']
    Points_pd_z = Points_pd['Z Coord']

    # print(Points_pd_x.head(5), Points_pd_y.head(5), Points_pd_z.head(5))


    #
    # for Segments 
    # df.loc[0]=['cat', 3] 

    # Segments = root[3]
    assert 'table' in str(list(Points)[0].tag).lower()

    tablee = Segments.getchildren()[0]
    rowws = tablee.getchildren()
    # roww_eg = tablee.getchildren()[0:20]
    # rowws = roww_eg
    # set column name 

    start_n = 0
    for row_i, roww in enumerate(rowws):
        if 'row' not in str(roww.tag).lower():
            start_n = row_i+1
        else: break


    for roww in rowws[:start_n+1]:
        print(roww.tag)
        if 'row' in str(roww.tag).lower():
            cells = roww.getchildren()
            columnname = []
            for cell in cells:
                datas = cell.getchildren()
                for data in datas:
                    columnname.append(data.text)
    # print('140', columnname)

    segments_pd =pd.DataFrame(columns=columnname)

    # segments_pd_points =segments_pd['Point IDs']




    # print(Points_pd.head(5))

    # obtain data
    for i, roww in enumerate(rowws[start_n+1:]):
        cells = roww.getchildren()
        templst = []
        for cell in cells:
            datas = cell.getchildren()
            for data in datas:
                templst.append(data.text)
        segments_pd.loc[i]= templst

    print(segments_pd.head(5))


    # for test
    # segments_pd = segments_pd[:10]


    points = Points_pd


    #
    # print(points.iloc[1])
    coords_lst = []
    for idx in range(len(points)):
        coords_lst.append([float(Points_pd_z.iloc[idx]), float(Points_pd_y.iloc[idx]), float(Points_pd_x.iloc[idx])])  #zyx from s

    coords = np.array(coords_lst).reshape(-1,3)
    # coords = coords - np.array([shift_xml_xyz[2], shift_xml_xyz[1], shift_xml_xyz[0]])
    # print(coords)  # [x,y,z]



    #  
    # input relations 
    filament_idx_lst = []
    filament_coord_lst = []
    segments = segments_pd
    # print(segments.columns)
    # ['Segment ID', 'CurvedLength', 'MeanRadius', 'Volume',
    #        'OrientationTheta', 'OrientationPhi', 'SubgraphID', 'ChordLength',
    #        'Tortuosity', 'TensorXX', 'TensorYY', 'TensorZZ', 'TensorXY',
    #        'TensorXZ', 'TensorYZ', 'Node ID #1', 'Node ID #2', 'Point IDs']


    # print(segments['Segment ID'])
    # print(segments['Point IDs'])
    correspond_id = segments['Segment ID']
    correspond_points_id = segments['Point IDs']
    # print(correspond_id)

    for idx in range(len(segments)):
    # for idx in range(5):
        assert idx ==  int(correspond_id[idx]) ## match 
        coord_idx_lst = [ int(pointidx) for pointidx in correspond_points_id[idx].split(',') ]# str
        # print(coord_idx_lst)
        # print(len(coord_idx_lst))
        temp_coordlst = []
        for pointidx in coord_idx_lst:
            temp_coordlst.append(coords[pointidx])  #   coords[pointidx] unit A  , voxel_coords

        filament_idx_lst.append(coord_idx_lst)
        filament_coord_lst.append(temp_coordlst)

    # print(len(filament_coord_lst))
    # print(len(filament_coord_lst[0]))



    # fill gap

    
    def get_line(point1, point2):
        # point1 to point2
        # resolution 40=4nm
        vect = np.array([point2[0] - point1[0], point2[1] - point1[1], point2[2] - point1[2]])
        unit_vect = vect / np.linalg.norm(vect) * resolution  # 4 nm
        thresh_len = resolution * 0.3 
        points_sum = np.linalg.norm(vect) // resolution + 1
        pointslst = [ point1 + unit_vect * i for i in range(int(points_sum))]
        while True:
            if np.linalg.norm(pointslst[-1] - point1) > np.linalg.norm(vect):
                pointslst = pointslst[:-1]
            else:
                break
        # if np.linalg.norm(pointslst[-1] - point2) < thresh_len:
        #     pointslst = pointslst[:-1]
        # print(unit_vect, pointslst)
        # voxel_loc = np.around(pointslst)
        # print(voxel_loc)
        return pointslst


    # extend coord to narrow interval
    filament_coord_extend_lst = [] 

    for filament_coord in filament_coord_lst:
        single_filament_coord_extend_lst = []
        for i in range(1, len(filament_coord)):
            coord1, coord2 = filament_coord[i-1], filament_coord[i]  # array
            temp_coord_extend_lst = get_line(coord1, coord2)
            # print('temp_coord_extend_lst', len(temp_coord_extend_lst))
            single_filament_coord_extend_lst.extend(temp_coord_extend_lst)
        single_filament_coord_extend_lst.append(filament_coord[-1])
        # print(len(single_filament_coord_extend_lst))
        filament_coord_extend_lst.append(single_filament_coord_extend_lst)
            # print(coord1, coord2)
            # print(filament_coord_extend_lst)  

    # filament_coord_lst
    coords_intotal1 = 0
    for single_filament_coord_extend in filament_coord_lst:
        coords_intotal1 += len(single_filament_coord_extend)
    print(f'coord intotal before: {coords_intotal1}')
    coords_intotal = 0
    for single_filament_coord_extend in filament_coord_extend_lst:
        coords_intotal += len(single_filament_coord_extend)
    print(f'coord intotal after: {coords_intotal}')


    ## save file
    img_columnname = [f'MT_{i}' for i in range(len(filament_coord_extend_lst))]
    # print(img_columnname)
    # df = pd.DataFrame(columns=[img_columnname])
    # print(df)
    points_all_dict = {}

    for idx, filament_vects in enumerate(filament_coord_extend_lst):
        points_all_dict[f'{img_columnname[idx]}'] =[]
        for coord in filament_vects:
            points_all_dict[f'{img_columnname[idx]}'] .append(list(coord))
            # df.loc[f'{img_columnname[idx]}'] = temp_vectlst

    # print(vect_all_dict)
    #
    # print(vect_all_dict.keys())
    
    for idx, filament_vects in enumerate(filament_coord_extend_lst):
        for coord in filament_vects:
            if len(points_all_dict[f'{img_columnname[idx]}']) <= 2:
                del points_all_dict[f'{img_columnname[idx]}']


    return points_all_dict
    # curpath = data_dir
    # outputname = f'{curpath}/{file_idx}_MT_filled_points.json'
    # with open(outputname, 'w') as f:
    #     json.dump(points_all_dict, f)

    # # print('Points gap filled.')
    # print('MT xml processed:', outputname)

    # return outputname




def fill_memholes(raw_img, start_idx):

    start_idx = 0

    img02 = raw_img
    img02_fill = np.zeros_like(img02)

    for i in range(img02.shape[0]):

        temp_img = skimage.morphology.closing(img02[i])
        img02_fill[i] = scipy.ndimage.morphology.binary_fill_holes(temp_img)
        

    img02_fill = img02_fill.astype(int)
    img02_fill_bkp = copy.deepcopy(img02_fill)
    # check overlap 
    for i in range(int(start_idx) + 1, img02.shape[0]):

        temp_img = copy.deepcopy(img02_fill[i])
        temp_img_lst = copy.deepcopy(img02_fill[i-1])
        temp_img_raw_lst = copy.deepcopy(img02_fill_bkp[i-1])

        # diff = np.sum(temp_img ^ temp_img_lst) / (np.sum(temp_img_lst) + 1e-6)


        # if True:
        temp_img_lst_edt = scipy.ndimage.distance_transform_edt(temp_img_lst)
        temp_img_revert = temp_img * (-1) + 1
        temp_cent_lst = skimage.feature.peak_local_max(temp_img_lst_edt, footprint=np.ones((21,) * (temp_img_lst_edt.ndim)), exclude_border=True)

        check_edt_lst = []
        filtered_temp_cent_lst = []
        for num_ in range(len(temp_cent_lst)):
            coord = temp_cent_lst[num_]
            check_edt_lst.append(temp_img_lst_edt[coord[0]][coord[1]])

            if temp_img[coord[0]][coord[1]] < 1 and temp_img_lst_edt[coord[0]][coord[1]] > 20 :
                filtered_temp_cent_lst.append(coord)


        ## judge if on the same circle, see if on same label
        # if len(filtered_temp_cent_lst) > 1 :

        temp_img_lst_label, _ = scipy.ndimage.label(temp_img_lst)
        temp_label = [ temp_img_lst_label[coord[0]][coord[1]] for coord in filtered_temp_cent_lst   ]

        cur_labellst = list(set(temp_label))



        # filtered cur labellst
        cur_labellst_new = []
        for label in cur_labellst:
            temp_label_img = np.zeros_like(temp_img_lst_label)
            temp_label_img[np.where(temp_img_lst_label == label)] = 1
            temp_label_img = scipy.ndimage.morphology.binary_dilation(temp_label_img, iterations=3)

            temp_label_img2 = np.zeros_like(temp_img_lst_label)
            temp_label_img2[np.where(temp_img_lst_label != 0)] = 1
            temp_label_img2[np.where(temp_img_lst_label == label)] = 0
            temp_label_img2 = scipy.ndimage.morphology.binary_dilation(temp_label_img2, iterations=1)

            overlap1_img = temp_img * temp_label_img
            overlap2_img = temp_img * temp_label_img2
            overlap1_img[np.where(overlap2_img == 1)]= 0

            if np.sum(overlap1_img) > 0:
                cur_labellst_new.append(label)

        cur_labellst = cur_labellst_new

        
        for label in cur_labellst:
            
            temp_seqlst = np.where(np.array(temp_label) == label )[0]
            cur_filtered_temp_cent_lst = [ filtered_temp_cent_lst[i] for i in temp_seqlst ]


            cent = [int(np.average([coord[0] for coord in cur_filtered_temp_cent_lst ])), int(np.average([coord[1] for coord in cur_filtered_temp_cent_lst ]))]

            centimg = np.ones_like(temp_img_lst_edt)
            centimg[cent[0]][cent[1]] = 0
            cent_edt = scipy.ndimage.distance_transform_edt(centimg)

            filtered_label_img_0 = np.zeros_like(temp_img_lst_edt)
            filtered_label_img_0[np.where(temp_img_lst_label == label)] = 1
            filtered_label_img = filtered_label_img_0 * temp_img

            filtered_label_img_0_dilation = scipy.ndimage.morphology.binary_dilation(filtered_label_img_0,iterations=5)
            filtered_label_lst_img = filtered_label_img_0_dilation * temp_img_raw_lst

            edge2cent_dist = cent_edt[np.where(filtered_label_img == 1.)]
            edge2cent_lst_dist = cent_edt[np.where(filtered_label_lst_img == 1.)]

            dist_change = (np.mean(edge2cent_dist) - np.mean(edge2cent_lst_dist))/5

            # for all coords in label
            coords = np.where(temp_img_lst_label == label)


            if dist_change < 0.4: # shrink
                modified_label_img = scipy.ndimage.morphology.binary_erosion(filtered_label_img_0, iterations=max(1, int((abs(dist_change)))))
            elif dist_change < 0.8:
            # else:
                modified_label_img = copy.deepcopy(filtered_label_img_0)
            else: 
                modified_label_img = scipy.ndimage.morphology.binary_dilation(filtered_label_img_0, iterations=1)


            modified_label_img = modified_label_img.astype(int)


            tmp_img02fill_slice = img02_fill[i] + modified_label_img
            tmp_img02fill_slice[np.where(tmp_img02fill_slice !=0)] = 1
            img02_fill[i] = tmp_img02fill_slice

    raw_img = img02_fill
    return raw_img




def fullfill_mem(mrcname, file_idx, data_dir):
    img01 = mrcfile.open( mrcname, permissive=True).data
    print(img01.shape)

    img02 = scipy.ndimage.morphology.binary_closing(img01) 
    img02_fill = np.zeros_like(img02)

    for i in range(img02.shape[0]):
        temp_img = skimage.morphology.closing(img02[i])
        img02_fill[i] = scipy.ndimage.morphology.binary_fill_holes(temp_img)

    # max_slide_sum = 0
    # max_idx = 0
    # for i in range(img02_fill.shape[0]):
        
    #     if np.sum(img02_fill[i]) > max_slide_sum:
    #         max_slide_sum = np.sum(img02_fill[i])
    #         max_idx = i
    # print('start idx:', max_idx)

    # start_idx = 0
    # one direction
    img02 = fill_memholes(img02, 0)
    img02_re = img02[::-1,:,:]
    img02_re = fill_memholes(img02_re, 0) #(img02_re.shape[0]-max_idx-1))
    img02 = img02_re[::-1,:,:]

    # another direction
    # img02 = fill_memholes(img02, max_idx)
    # img02_re = img02[::-1,:,:]
    # img02_re = fill_memholes(img02_re, (img02_re.shape[0]-max_idx-1))
    # img02 = img02_re[::-1,:,:]

    # save mrc
    curpath = data_dir
    img02 = img02.astype(np.int8)
    outputname = f'{curpath}/{file_idx}_fullfill_vesicle_mask.mrc'
    with mrcfile.new(outputname, overwrite=True) as mrc:
        mrc.set_data(img02 )
    print('mrc processed:', outputname)

    return outputname





#----------------------

class FA_cleaner():
    '''
    clean False FA from images
    
    '''
    def __init__(self, ) -> None:

        pass

    def getdata(self,FA_image, FA_coords, false_coords):
        self.img = FA_image
        self.coords = FA_coords
        self.false_coords = false_coords


    def clean(self):
        '''
        if dist < 14 ,remove
        '''
        self.newimg = copy.deepcopy(self.img)
        self.replace_box()



    def replace_box(self):
        
        self.recognized_edge = 12
        recognized_edge = self.recognized_edge
        
        self.replace_edge = 9


        edge = recognized_edge

        self.empty_box = np.zeros([2*edge+1,2*edge+1,2*edge+1])   

        dist_array = cdist(self.false_coords, self.coords)
        falsecoord_pair = np.where(dist_array <= 10) # 14* sqrt(3)
        self.removed_coordlst = []
        self.confirmed_falsecoord = []
        for i in range(len(falsecoord_pair[0])):
            false_coord1 = self.false_coords[falsecoord_pair[0][i]]
            match_coord2 = self.coords[falsecoord_pair[1][i]]
            # print(147, false_coord1, match_coord2)

            # if True:
            # print(false_coord1[0], match_coord2[0],  (false_coord1[0] > (match_coord2[0]-edge)) and (false_coord1[0] < (match_coord2[0]+edge)))
            if (false_coord1[0] > (match_coord2[0]-edge)) and (false_coord1[0] < (match_coord2[0]+edge)) and\
                (false_coord1[1] > (match_coord2[1]-edge)) and (false_coord1[1] < (match_coord2[1]+edge) )and \
                (false_coord1[2] > (match_coord2[2]-edge)) and (false_coord1[2] < (match_coord2[2]+edge)):
                # print(111)
                # replace image
                self.confirmed_falsecoord.append(false_coord1)
                self.removed_coordlst.append(match_coord2)
                # x,y,z = match_coord2[0], match_coord2[1], match_coord2[2]
                self.replace_local_box(match_coord2, replace_edge = self.replace_edge)


        # print(self.false_coords, self.confirmed_falsecoord)
        
        self.confirmed_facoord = []
        for coord in self.coords:
            if coord not in self.removed_coordlst:
                self.confirmed_facoord.append(coord)

        self.revised_coord = []
        for coord in self.false_coords:
            if coord not in self.confirmed_falsecoord:
                self.revised_coord.append(coord)

    def replace_local_box(self, coord, replace_edge):
        x_min, x_max = max(coord[0]-replace_edge, 0), min(coord[0]+replace_edge, self.img.shape[0])
        y_min, y_max = max(coord[1]-replace_edge, 0), min(coord[1]+replace_edge, self.img.shape[1])
        z_min, z_max = max(coord[2]-replace_edge, 0), min(coord[2]+replace_edge, self.img.shape[2])

        cur_box = copy.deepcopy(self.img[x_min:x_max+1, y_min:y_max+1,z_min:z_max+1])
        # cent_coord = [math.floor((x_min + x_max)/2),math.floor((y_min + y_max)/2),math.floor((z_min + z_max)/2)]
        cent_coord = coord

        for x in range(x_min,x_max+1):
            for y in range(y_min, y_max+1):
                for z in range(z_min, z_max+1):
                    if (x - cent_coord[0])**2 + (y - cent_coord[1])**2 + (z - cent_coord[2])**2 <= replace_edge**2:                    
                        cur_box[x-x_min][y-y_min][z-z_min] = 0
        
        self.newimg[x_min:x_max+1, y_min:y_max+1,z_min:z_max+1] = cur_box





    def check_removed_coord(self):
        for coord in self.removed_coordlst:
            print(coord)
            temp_img = self.img[coord[0], coord[1]-25:coord[1]+26, coord[2]-25:coord[2]+26]
            plt.imshow(temp_img)
            plt.title(f'{coord}')
            plt.show()
            plt.close()

    def check_box(self,coords):
        for coord in coords:
            print(coord)
            temp_img = self.img[max(coord[0]-10,0):min(coord[0]+10,self.img.shape[0]), max((coord[1]-25), 0):min((coord[1]+26),self.img.shape[1]), max((coord[2]-25),0):min((coord[2]+26),self.img.shape[2])]
            for i in range(3):
                cent = int(temp_img.shape[0]/2) 
                plt.imshow(temp_img[cent-5 + 5*i])
                plt.title(f'{coord} {i}')
                plt.show()
                plt.close()

    def check_revisedbox(self,coords):
        for coord in coords:
            print(coord)
            temp_img = self.newimg[max(coord[0]-10,0):min(coord[0]+10,self.img.shape[0]), max((coord[1]-25), 0):min((coord[1]+26),self.img.shape[1]), max((coord[2]-25),0):min((coord[2]+26),self.img.shape[2])]
            for i in range(3):
                cent = int(temp_img.shape[0]/2) 
                plt.imshow(temp_img[cent-5 + 5*i])
                plt.title(f'{coord} {i}')
                plt.show()
                plt.close()




def filter_fa_point(fa_raw_points_name, false_fa_points_name, image_name):
    def read_mrcimg(mrcname):
        img = mrcfile.open(mrcname, permissive=True).data
        print(mrcname, 'shape',img.shape)
        return img[:,::-1,:]
    
    def read_xml(xmlname, imgdata_shape):
        '''
        read coordinates from FA xml
        
        xyz in img
        #(235, 1440, 1024)

        xyz in xml
        # 1024, Y = 1440-y, z 

        '''
        xmlname = f'{xmlname}'
        f = open(xmlname, encoding='utf-8') 
        xml_txt = f.read()
        root = ET.fromstring(xml_txt)
        FAs_lst = []
        for ii, child in enumerate(root):
            # if ii > 1: break
            # print(child.tag, child.attrib)
            single_fa = child
            for jj, child in enumerate(single_fa):
                if child.tag == 'PickPosition': 
                    # print(child.tag)
                    FA_lst =[int(child.attrib['Z']), imgdata_shape[1] - int(child.attrib['Y']), int(child.attrib['X'])]
                    FAs_lst.append(FA_lst)
                # print(child.tag, child.attrib)
        # (Z_xml, 1440 - Y_xml, X_xml)
    
        # print(FAs_lst)
        return FAs_lst

    def read_csv(csvname):
        false_FA = pd.read_csv(csvname)
        false_FAlst = []
        for i in range(false_FA.shape[0]):
            coord = false_FA.iloc[i]
            false_FAlst.append([coord[2], coord[1], coord[0]])


        return false_FAlst
    
    image = read_mrcimg(image_name)
    imgdata_shape = image.shape

    fa_raw_points = read_xml(fa_raw_points_name, imgdata_shape)
    false_fa_points = read_csv(false_fa_points_name)


    fac = FA_cleaner()  
    fac.getdata(image, fa_raw_points, false_fa_points)
    fac.clean()
    faimg = fac.newimg 

    return fac.confirmed_facoord, fac.newimg 


