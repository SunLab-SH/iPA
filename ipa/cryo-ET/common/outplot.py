import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.ndimage
import copy
from matplotlib.image import NonUniformImage
from scipy.interpolate import griddata #引入scipy中的二维插值库
import time


def vis_actin_surround_actin(actin_vect_map_lst, file_idx, actin_num, point_num, ):
    
    # vis_img = np.zeros([603, 603])
    # n =600
    points_x = []
    points_y = []
    # points = []
    for idx, actin_vects in enumerate(actin_vect_map_lst):
        # print(vis_img.shape)
        if len(actin_vects) >0:
            for vect in actin_vects:
                points_x.append(vect[0])
                points_y.append(vect[1])
                # points.append([vect[0],vect[1]])

    bins_ = 50
    fig = plt.figure(figsize=(8,6), dpi=1200)
    hist, xedges, yedges = np.histogram2d(points_x, points_y, bins=bins_, range =[[-1000,1000],[-1000,1000]] )
    # print(hist)
    values = hist.reshape(-1, 1)
    points = []
    for i in range(hist.shape[0]):
        for j in range(hist.shape[1]):
            points.append([i,j])
    # points_x, points_y = np.mgrid[1:100:100j, 1:100:100j]
    # points = np.array([ [points_x[i], points_y[i]] for i in range(len(points_x)) ])

    grid_x, grid_y = np.mgrid[1:hist.shape[0]:1000j, 1:hist.shape[1]:1000j]
    # points = np.array(points).reshape(-1,2)
    # values = hist.
    # # points = np.mgrid[1:100:100j, 1:100:100j]
    # grid_x, grid_y = np.mgrid[1:100:1000j, 1:100:1000j]
    grid_z0 = griddata(points, values, (grid_x, grid_y), method='linear',fill_value=0)

    grid_z0[np.where(grid_z0 < (np.max(grid_z0) * 0.3) )] = np.nan
    grid_z0[np.where(grid_z0 < (np.max(grid_z0) * 0.45) )] = 0
    
    # print(grid_z0)
    plt.imshow(grid_z0,interpolation='bilinear', extent=(-100,100,-100,100),cmap= 'RdYlBu_r')
    plt.title(f'{file_idx} Actin surrounding probability density ditribution testClim \n climmax = {np.max(grid_z0)}, actin num {actin_num}, points num {point_num}' )
    
    # plt.grid(True)
    plt.colorbar()
    #plt.clim(vmin=0, vmax=1100)
    
    # plt.show()
    # plt.close()

    # hist[np.where(hist < (np.max(hist) * 0.2) )] = np.nan
    # hist[np.where(hist < (np.max(hist) * 0.3) )] = 0
    # # X,Y = np.meshgrid(xedges,yedges)
    # plt.imshow(grid_z0, interpolation='bilinear', cmap=  'RdYlBu_r')

    # plt.grid(True)
    # plt.colorbar()
    # plt.show()


def vis_actin_surround_actin_dist(dist, file_idx, actin_num, point_num, ):
    
    sns.distplot(a=dist, kde=True)
    # sns.kdeplot(data=dist, shade=True)
    plt.title(f'{file_idx} actin to actin distance distribution \n actin num {actin_num}, points num {point_num}' )
    plt.ylabel('probability density')
    plt.xlabel('actin - actin distance (nm)')
    # plt.xlim(0,1500)
    # plt.ylim(0, 0.002)
    # plt.savefig(f'{file_idx}_distance_distribution.png')







def actin_to_vesicle_dist_hist(dist, file_idx, actin_num, point_num):
    sns.distplot(a=dist, kde=True)
    # sns.kdeplot(data=dist, shade=True)
    plt.title(f'{file_idx} actin to vesicle distance distribution \n actin num {actin_num}, points num {point_num}' )
    plt.ylabel('probability density')
    plt.xlabel('actin - vesicle distance (nm)')
    # plt.xlim(0,1500)
    # plt.ylim(0, 0.002)
    # plt.savefig(f'{file_idx}_distance_distribution.png')

def vis_microtube_surround_actin(actin_vect_map_lst, file_idx, actin_num, point_num):
    
    vis_img = np.zeros([201, 201])
    for idx, actin_vects in enumerate(actin_vect_map_lst):
        # print(vis_img.shape)
        if len(actin_vects) >0:
            for vect in actin_vects:
                temp_vect = np.around(vect)
                temp_vect = temp_vect / np.array([10,10,10])
                temp_vect = [int(temp_vect[0]+vis_img.shape[0]/2), int(temp_vect[1]+100)]
                vis_img[temp_vect[0]][temp_vect[1]] = vis_img[temp_vect[0]][temp_vect[1]] + 1
                # vis_img[temp_vect[0]+1][temp_vect[1]] = vis_img[temp_vect[0]][temp_vect[1]] + 1
                # vis_img[temp_vect[0]-1][temp_vect[1]] = vis_img[temp_vect[0]][temp_vect[1]] + 1
                # vis_img[temp_vect[0]][temp_vect[1]+1] = vis_img[temp_vect[0]][temp_vect[1]] + 1
                # vis_img[temp_vect[0]][temp_vect[1]-1] = vis_img[temp_vect[0]][temp_vect[1]] + 1

    print(np.max(vis_img))
   
    # vis_img2 = copy.deepcopy(vis_img)
    # vis_img2[np.where(vis_img2 < (np.max(vis_img2) * 0.3) )] = 0
    # vis_img_filtered = scipy.ndimage.gaussian_filter(vis_img2, sigma=0.5)

    # print('filteredmax', np.max(vis_img_filtered))

    fig = plt.figure(figsize=(8,6), dpi=1200)
    plt.imshow(vis_img, cmap=  'RdYlBu_r'  ) #'gray_r')  #RdYlBu
    plt.title(f'{file_idx} microtube surrounding probability density ditribution \n clim {np.max(vis_img)}, actin num {actin_num}, points num {point_num}')
    
    # plt.axis('off')
    plt.colorbar() 
    plt.clim(0, 250) #colorbar_weimin
    plt.savefig(f'{file_idx}_microtube_surrounding_probability_density_ditribution.pdf')



def vis_microtube_surround_actin_dist(dist, file_idx, actin_num, point_num, ):
    
    sns.distplot(a=dist, kde=True)
    # sns.kdeplot(data=dist, shade=True)
    plt.title(f'{file_idx} actin to MT distance distribution \n actin num {actin_num}, points num {point_num}' )
    plt.ylabel('probability density')
    plt.xlabel('actin - MT distance (nm)')
    # plt.xlim(0,1500)
    # plt.ylim(0, 0.002)
    # plt.savefig(f'{file_idx}_distance_distribution.png')





def actindist2mito_hist(dist, file_idx, actin_num, point_num):

    sns.distplot(a=dist, kde=True)
    # sns.kdeplot(data=dist, shade=True)
    plt.title(f'{file_idx} actin to mito mem distance distribution \n actin num {actin_num}, points num {point_num}' )
    plt.ylabel('probability density')
    plt.xlabel('actin - mito membrane distance (nm)')
    # plt.xlim(0,1500)
    # plt.ylim(0, 0.002)
    # plt.savefig(f'{file_idx}_distance_distribution.png')


def actindist2er_hist(dist, file_idx, actin_num, point_num):
    sns.distplot(a=dist, kde=True)
    # sns.kdeplot(data=dist, shade=True)
    plt.title(f'{file_idx} actin to er mem distance distribution \n actin num {actin_num}, points num {point_num}' )
    plt.ylabel('probability density')
    plt.xlabel('actin - er membrane distance (nm)')
    # plt.xlim(0,1500)
    # plt.ylim(0, 0.002)
    # plt.savefig(f'{file_idx}_distance_distribution.png')


def dist_mt2isghist(dist, file_idx, actin_num, point_num):
    sns.distplot(a=dist, kde=True)
    # sns.kdeplot(data=dist, shade=True)
    plt.title(f'{file_idx} MT to isg mem distance distribution \n actin num {actin_num}, points num {point_num}' )
    plt.ylabel('probability density')
    plt.xlabel('MT - isg membrane distance (nm)')
    # plt.xlim(0,1500)
    # plt.ylim(0, 0.002)
    # plt.savefig(f'{file_idx}_distance_distribution.png')

def dist_mt2mitohist(dist, file_idx, actin_num, point_num):
    sns.distplot(a=dist, kde=True)
    # sns.kdeplot(data=dist, shade=True)
    plt.title(f'{file_idx} MT to mito mem distance distribution \n actin num {actin_num}, points num {point_num}' )
    plt.ylabel('probability density')
    plt.xlabel('MT - mito membrane distance (nm)')
    # plt.xlim(0,1500)
    # plt.ylim(0, 0.002)
    # plt.savefig(f'{file_idx}_distance_distribution.png')

def dist_mt2erhist(dist, file_idx, actin_num, point_num):
    sns.distplot(a=dist, kde=True)
    # sns.kdeplot(data=dist, shade=True)
    plt.title(f'{file_idx} MT to ER mem distance distribution \n actin num {actin_num}, points num {point_num}' )
    plt.ylabel('probability density')
    plt.xlabel('MT -ER membrane distance (nm)')
    # plt.xlim(0,1500)
    # plt.ylim(0, 0.002)
    # plt.savefig(f'{file_idx}_distance_distribution.png')

def dist_isg2mito_hist(dist, file_idx, point_num):
    sns.distplot(a=dist, kde=True)
    # sns.kdeplot(data=dist, shade=True)
    plt.title(f'{file_idx} isgmem to mitomem distance distribution \n points num {point_num}' )
    plt.ylabel('probability density')
    plt.xlabel('isgedge - mitoedge distance (nm)')
    # plt.xlim(0,1500)
    # plt.ylim(0, 0.002)
    # plt.savefig(f'{file_idx}_distance_distribution.png')

def dist_isg2er_hist(dist, file_idx, point_num):
    sns.distplot(a=dist, kde=True)
    # sns.kdeplot(data=dist, shade=True)
    plt.title(f'{file_idx} isgmem to ermem distance distribution \n points num {point_num}' )
    plt.ylabel('probability density')
    plt.xlabel('isgedge - eredge distance (nm)')
    # plt.xlim(0,1500)
    # plt.ylim(0, 0.002)
    # plt.savefig(f'{file_idx}_distance_distribution.png')



def dist_focal2focal_hist(dist, file_idx, point_num):
    sns.distplot(a=dist, kde=True)
    # sns.kdeplot(data=dist, shade=True)
    plt.title(f'{file_idx} focal to focal distance distribution \n points num {point_num}' )
    plt.ylabel('probability density')
    plt.xlabel('focal - focal distance (nm)')  


def dist_focal2focal_condition_hist(dist, file_idx, point_num):
    sns.distplot(a=dist, kde=True)
    # sns.kdeplot(data=dist, shade=True)
    plt.title(f'{file_idx} focal to focal distance condition distribution \n points num {point_num}' )
    plt.ylabel('probability density')
    plt.xlabel('focal - focal distance (nm)')  



def dist_focal2ca_hist(dist, file_idx, point_num, ca_num):
    sns.distplot(a=dist, kde=True)
    # sns.kdeplot(data=dist, shade=True)
    plt.title(f'{file_idx} focal to calcium distance distribution \n fa num {point_num} ca num {ca_num}' )
    plt.ylabel('probability density')
    plt.xlabel('focal - calcium distance (nm)')  


def dist_focal2ca_condition_hist(dist, file_idx, filelst):
    sns.distplot(a=dist, kde=True)
    # sns.kdeplot(data=dist, shade=True)
    plt.title(f'{file_idx} focal to calcium distance condition distribution \n {filelst}' )
    plt.ylabel('probability density')
    plt.xlabel('focal - calcium distance (nm)')  


#________________________



def dist_mito2er_hist(dist, file_idx, point_num):
    sns.distplot(a=dist, kde=True)
    # sns.kdeplot(data=dist, shade=True)
    plt.title(f'{file_idx} mitomem to ermem distance distribution \n points num {point_num}' )
    plt.ylabel('probability density')
    plt.xlabel('mitoedge - eredge distance (nm)')
    # plt.xlim(0,1500)
    # plt.ylim(0, 0.002)
    # plt.savefig(f'{file_idx}_distance_distribution.png')


def dist_mito2er_condition_hist(dist, file_idx, point_num, conditionnamelst):
    sns.distplot(a=dist, kde=True)
    # sns.kdeplot(data=dist, shade=True)
    plt.title(f'{file_idx} mitomem to ermem distance distribution \n points num {point_num} \n {conditionnamelst}' )
    plt.ylabel('probability density')
    plt.xlabel('mitoedge - eredge distance (nm)')
    # plt.xlim(0,1500)
    # plt.ylim(0, 0.002)
    # plt.savefig(f'{file_idx}_distance_distribution.png')


def vis_actin_surround(actin_vect_map_lst, file_idx, actin_num, point_num):
    
    vis_img = np.zeros([203, 203])
    # n =600
    for idx, actin_vects in enumerate(actin_vect_map_lst):
        # print(vis_img.shape)
        if len(actin_vects) >0:
            for vect in actin_vects:
                temp_vect = np.around(vect)
                temp_vect = temp_vect / np.array([10,10,10]) * np.array([1,1,1])
                temp_vect = [int(temp_vect[0]+vis_img.shape[0]/2), int(temp_vect[1]+vis_img.shape[1]/2)]
                vis_img[temp_vect[0]][temp_vect[1]] = vis_img[temp_vect[0]][temp_vect[1]] + 1
                vis_img[temp_vect[0]+1][temp_vect[1]] = vis_img[temp_vect[0]+1][temp_vect[1]] + 1
                vis_img[temp_vect[0]-1][temp_vect[1]] = vis_img[temp_vect[0]-1][temp_vect[1]] + 1
                vis_img[temp_vect[0]][temp_vect[1]+1] = vis_img[temp_vect[0]][temp_vect[1]+1] + 1
                vis_img[temp_vect[0]][temp_vect[1]-1] = vis_img[temp_vect[0]][temp_vect[1]-1] + 1

    print(np.max(vis_img))
    # vis_img = scipy.ndimage.binary_dilation(vis_img,iterations=3).astype(int)
    # intensity_threshold = [0 + 0.1 * x for x in range(11)]
    # for thresh in intensity_threshold:
    #     vis_img2 = copy.deepcopy(vis_img)
    #     vis_img2[np.where(vis_img2 < (np.max(vis_img2)*thresh) )] = 0
    #     vis_img_filtered = scipy.ndimage.gaussian_filter(vis_img2, sigma=0.5)
    #     fig = plt.figure(figsize=(8,6), dpi=1200)
    #     plt.imshow(vis_img_filtered, cmap=  'RdYlBu_r'  ) #'gray_r')  #RdYlBu
    #     plt.title(f'{file_idx} Actin surrounding probability density ditribution\n sigma = {0.5}, thresh = {thresh}, max {np.max(vis_img_filtered)} \n actin num {actin_num}, points num {point_num}')
        
    #     # plt.axis('off')
    #     plt.colorbar()
    #     plt.savefig(f'./temp/{file_idx}_actin_surrounding_sigma0.5_thresh_{thresh}.png')
    #     plt.show()
    #     plt.close()


    # sigmalst = [0 + 0.1 * x for x in range(11)]
    # for sigma_ in sigmalst:
    #     vis_img_filtered = scipy.ndimage.gaussian_filter(vis_img, sigma=sigma_)
    #     fig = plt.figure(figsize=(8,6), dpi=1200)
    #     plt.imshow(vis_img_filtered, cmap=  'RdYlBu_r'  ) #'gray_r')  #RdYlBu
    #     plt.title(f'{file_idx} Actin surrounding probability density ditribution\n sigma = {sigma_}, max {np.max(vis_img_filtered)} \n actin num {actin_num}, points num {point_num}')
        
    #     # plt.axis('off')
    #     plt.colorbar()
    #     plt.show()
    #     plt.close()

    # vis_img2 = copy.deepcopy(vis_img)
    # vis_img2[np.where(vis_img2 < (np.max(vis_img2) * 0.3) )] = 0
    # vis_img_filtered = scipy.ndimage.gaussian_filter(vis_img2, sigma=0.5)

    # print('filteredmax', np.max(vis_img_filtered))

    fig = plt.figure(figsize=(8,6), dpi=1200)
    plt.imshow(vis_img, cmap=  'RdYlBu_r'  ) #'gray_r')  #RdYlBu
    plt.title(f'{file_idx} Actin surrounding probability density ditribution \n actin num {actin_num}, points num {point_num}')
    
    # plt.axis('off')
    plt.colorbar()
    # plt.savefig(f'{file_idx}_actin_surrounding_probability_density_ditribution.pdf')

def vis_actin_surround2(actin_vect_map_lst, file_idx, actin_num, point_num):
    
    vis_img = np.zeros([403, 403])
    # n =600
    for idx, actin_vects in enumerate(actin_vect_map_lst):
        # print(vis_img.shape)
        if len(actin_vects) >0:
            for vect in actin_vects:
                temp_vect = np.around(vect)
                temp_vect = temp_vect / np.array([10,10,10]) * np.array([2,2,2])
                temp_vect = [int(temp_vect[0]+vis_img.shape[0]/2), int(temp_vect[1]+vis_img.shape[1]/2)]
                vis_img[temp_vect[0]][temp_vect[1]] = vis_img[temp_vect[0]][temp_vect[1]] + 1
                vis_img[temp_vect[0]+1][temp_vect[1]] = vis_img[temp_vect[0]+1][temp_vect[1]] + 1
                vis_img[temp_vect[0]-1][temp_vect[1]] = vis_img[temp_vect[0]-1][temp_vect[1]] + 1
                vis_img[temp_vect[0]][temp_vect[1]+1] = vis_img[temp_vect[0]][temp_vect[1]+1] + 1
                vis_img[temp_vect[0]][temp_vect[1]-1] = vis_img[temp_vect[0]][temp_vect[1]-1] + 1

    print(np.max(vis_img))
    # vis_img = scipy.ndimage.binary_dilation(vis_img,iterations=3).astype(int)
    # intensity_threshold = [0 + 0.1 * x for x in range(11)]
    # for thresh in intensity_threshold:
    #     vis_img2 = copy.deepcopy(vis_img)
    #     vis_img2[np.where(vis_img2 < (np.max(vis_img2)*thresh) )] = 0
    #     vis_img_filtered = scipy.ndimage.gaussian_filter(vis_img2, sigma=0.5)
    #     fig = plt.figure(figsize=(8,6), dpi=1200)
    #     plt.imshow(vis_img_filtered, cmap=  'RdYlBu_r'  ) #'gray_r')  #RdYlBu
    #     plt.title(f'{file_idx} Actin surrounding probability density ditribution\n sigma = {0.5}, thresh = {thresh}, max {np.max(vis_img_filtered)} \n actin num {actin_num}, points num {point_num}')
        
    #     # plt.axis('off')
    #     plt.colorbar()
    #     plt.savefig(f'./temp/{file_idx}_actin_surrounding_sigma0.5_thresh_{thresh}.png')
    #     plt.show()
    #     plt.close()


    # sigmalst = [0 + 0.1 * x for x in range(11)]
    # for sigma_ in sigmalst:
    #     vis_img_filtered = scipy.ndimage.gaussian_filter(vis_img, sigma=sigma_)
    #     fig = plt.figure(figsize=(8,6), dpi=1200)
    #     plt.imshow(vis_img_filtered, cmap=  'RdYlBu_r'  ) #'gray_r')  #RdYlBu
    #     plt.title(f'{file_idx} Actin surrounding probability density ditribution\n sigma = {sigma_}, max {np.max(vis_img_filtered)} \n actin num {actin_num}, points num {point_num}')
        
    #     # plt.axis('off')
    #     plt.colorbar()
    #     plt.show()
    #     plt.close()

    vis_img2 = copy.deepcopy(vis_img)
    vis_img2[np.where(vis_img2 < (np.max(vis_img2) * 0.3) )] = np.nan
    vis_img2[np.where(vis_img2 < (np.max(vis_img2) * 0.5) )] = 0.1
    vis_img_filtered = scipy.ndimage.gaussian_filter(vis_img2, sigma=0.5)

    # print('filteredmax', np.max(vis_img_filtered))

    fig = plt.figure(figsize=(8,6), dpi=1200)
    plt.imshow(vis_img2, cmap=  'RdYlBu_r'  ) #'gray_r')  #RdYlBu
    plt.title(f'{file_idx} Actin surrounding probability density ditribution \n actin num {actin_num}, points num {point_num}')
    
    # plt.axis('off')
    plt.colorbar()
    # plt.savefig(f'{file_idx}_actin_surrounding_probability_density_ditribution.pdf')



def vis_actindistangle(actins_min_d_lst, actins_angle_lst, file_idx, actin_num, point_num):

    max_ = np.max(actins_min_d_lst)
    points_x = actins_min_d_lst
    points_y = actins_angle_lst

    if len(actins_min_d_lst)>300:    
        bins_ = [50,50]
    elif len(actins_min_d_lst)>100:
        bins_ = [50,50]
    else:
        bins_ = [10, 10]

    fig = plt.figure(figsize=(8,6), dpi=1200)
    hist, xedges, yedges = np.histogram2d(points_x, points_y, bins=bins_, range =[[0,np.ceil(max_/100)*100],[0,90]] )
    # print(hist)
    
    values = hist.reshape(-1, 1)
    points = []
    for i in range(hist.shape[0]):
        for j in range(hist.shape[1]):
            points.append([i,j])

    # points_x, points_y = np.mgrid[1:100:100j, 1:100:100j]
    # points = np.array([ [points_x[i], points_y[i]] for i in range(len(points_x)) ])

    grid_x, grid_y = np.mgrid[0:hist.shape[0]:1000j, 0:hist.shape[1]:1000j]


    grid_z0 = griddata(points, values, (grid_x, grid_y), method='linear',fill_value = 0)


    # grid_z0[np.where(grid_z0 < (np.max(grid_z0) * 0.4) )] = np.nan
    # grid_z0[np.where(grid_z0 < (np.max(grid_z0) * 0.45) )] = 0

    # print(grid_z0)
    plt.imshow(grid_z0,interpolation='bilinear', cmap= 'RdYlBu_r', origin='lower')
    # plt.imshow(grid_z0,interpolation='bilinear', 
    #             extent=( 0,max_, 0,90,),
    #             cmap= 'RdYlBu_r', 
    #             origin='lower')
    plt.title(f'{file_idx} Actin dist VS angle \n actin num {actin_num}, points num {point_num}')
    plt.xticks(np.arange(0,1000,200),np.arange(0, int(np.ceil(max_/100)*100), int(np.ceil(max_/100)*10*2)))
    plt.yticks(np.arange(0,990,110),np.arange(0,90,10))
    plt.xlabel('Distance (nm)')
    plt.ylabel('Angle')
    # plt.grid(True)
    plt.colorbar()
    





def vis_actindistangle_condition(actins_min_d_lst, actins_angle_lst, file_idx,actinnum, point_num):

    max_ = np.max(actins_min_d_lst)
    points_x = actins_min_d_lst
    points_y = actins_angle_lst

    if len(actins_min_d_lst)>300:    
        bins_ = [50,50]
    elif len(actins_min_d_lst)>100:
        bins_ = [50,50]
    else:
        bins_ = [10, 10]

    fig = plt.figure(figsize=(8,6), dpi=1200)
    hist, xedges, yedges = np.histogram2d(points_x, points_y, bins=bins_, range =[[0,np.ceil(max_/100)*100],[0,90]] )
    # print(hist)
    
    values = hist.reshape(-1, 1)
    points = []
    for i in range(hist.shape[0]):
        for j in range(hist.shape[1]):
            points.append([i,j])

    # points_x, points_y = np.mgrid[1:100:100j, 1:100:100j]
    # points = np.array([ [points_x[i], points_y[i]] for i in range(len(points_x)) ])

    grid_x, grid_y = np.mgrid[0:hist.shape[0]:1000j, 0:hist.shape[1]:1000j]


    grid_z0 = griddata(points, values, (grid_x, grid_y), method='linear',fill_value = 0)


    # grid_z0[np.where(grid_z0 < (np.max(grid_z0) * 0.4) )] = np.nan
    # grid_z0[np.where(grid_z0 < (np.max(grid_z0) * 0.45) )] = 0

    # print(grid_z0)
    plt.imshow(grid_z0,interpolation='bilinear', cmap= 'RdYlBu_r', origin='lower')
    # plt.imshow(grid_z0,interpolation='bilinear', 
    #             extent=( 0,max_, 0,90,),
    #             cmap= 'RdYlBu_r', 
    #             origin='lower')
    plt.title(f'{file_idx} Actin dist VS angle \n actin num {actinnum}, points num {point_num}')
    plt.xticks(np.arange(0,1000,200),np.arange(0, int(np.ceil(max_/100)*100), int(np.ceil(max_/100)*10*2)))
    plt.yticks(np.arange(0,990,110),np.arange(0,90,10))
    plt.xlabel('Distance (nm)')
    plt.ylabel('Angle')
    # plt.grid(True)
    plt.colorbar()
    plt.clim(vmin=0, vmax=2500) #colorbar_weimin
    



def fa_actin_angle_dist_map_plot(angle_dist_pd, dataid):

    df = angle_dist_pd

    sns.jointplot(x=df['Distance'], y=df['Angle'],
                data=df,
                color='k',
                kind='scatter',             #reg添加线性回归线
                height=16,
                ratio=5,
                marginal_kws=dict(bins=20, rug=True,))

   # plt.xlabel('Distance (nm)')
    # plt.ylabel('Angle')
    plt.title(f'{dataid} fa_actin_angle_dist_pair plot in nm \n point num {angle_dist_pd.shape[0]}')
    # plt.show()




def fa_actin_angle_dist_map_condition_plot(angle_dist_pd, conditionid, dataidlst ):

    df = angle_dist_pd

    sns.jointplot(x=df['Distance'], y=df['Angle'],
                data=df,
                color='k',
                kind='scatter',             #reg添加线性回归线
                height=16,
                ratio=5,
                marginal_kws=dict(bins=20, rug=True,))

    # plt.xlabel('Distance (nm)')
    # plt.ylabel('Angle')
    
    plt.title(f'{conditionid} fa_actin_angle_dist_pair plot in nm \n point num {angle_dist_pd.shape[0]} \n {dataidlst}')
    # plt.show()


def fa_actin_angle_angle_pair_plot(angle_angle_pd, dataid):

    df = angle_angle_pd

    sns.jointplot(x=df['Angle mem'], y=df['Angle actin'],
                data=df,
                color='k',
                kind='scatter',             #reg添加线性回归线
                height=16,
                ratio=5,
                marginal_kws=dict(bins=20, rug=True,))
    plt.xlim([0,90])
    plt.ylim([0,90])
    plt.title(f'{dataid} fa_actin_angle_angle_pair_plot plot \n point num {angle_angle_pd.shape[0]}')
    # plt.show()



def fa_actin_angle_angle_pair_condition_plot(angle_angle_pd, conditionid, dataidlst):

    df = angle_angle_pd

    sns.jointplot(x=df['Angle mem'], y=df['Angle actin'],
                data=df,
                color='k',
                kind='scatter',             #reg添加线性回归线
                height=16,
                ratio=5,
                marginal_kws=dict(bins=20, rug=True,))
    plt.xlim([0,90])
    plt.ylim([0,90])
    plt.title(f'{conditionid} fa_actin_angle_angle_pair_plot plot \n point num {angle_angle_pd.shape[0]} \n {dataidlst}')
    # plt.show()







def dist2mem_condition_plot_old(distss_prob_in_one_condition, condition, validpath, bins = 100, x_lim = 1200):


    ave_distss_prob_in_one_condition = np.average(distss_prob_in_one_condition, axis=0)
    std_distss_prob_in_one_condition = np.std(distss_prob_in_one_condition, axis=0, ddof = 1)
    se_distss_prob_in_one_condition = std_distss_prob_in_one_condition / math.sqrt(len(validpath))

    print(len(ave_distss_prob_in_one_condition))

    x = [x_lim/bins * (0.5+i) for i in range(len(ave_distss_prob_in_one_condition))]

    fig =plt.figure(figsize=(8,6),  dpi=1200)
    lower_bound = [max(0, ave_distss_prob_in_one_condition[i] - se_distss_prob_in_one_condition[i]) for i in range(len(ave_distss_prob_in_one_condition))]
    upper_bound = [ave_distss_prob_in_one_condition[i] + se_distss_prob_in_one_condition[i] for i in range(len(ave_distss_prob_in_one_condition))]
    plt.plot(x, ave_distss_prob_in_one_condition)
    plt.fill_between(x, lower_bound, upper_bound, alpha = 0.3)

    plt.xlabel('Distance to membrane (nm)')
    plt.ylabel('Probability density')
    plt.title(f'{condition}_condition_distance_density_distribution bin = {bins}\n  {[i for i in validpath]}')
    plt.tight_layout()

    # plt.legend()

    # plt.savefig(f'{sign}_condition_distance_density_distribution.pdf')
    # plt.show()
    # plt.close()


def dist2mem_condition_plot(distss_prob_in_one_condition, condition, validpath, bins = 100, x_lim = 1200):


    # ave_distss_prob_in_one_condition = np.average(distss_prob_in_one_condition, axis=0)
    # std_distss_prob_in_one_condition = np.std(distss_prob_in_one_condition, axis=0, ddof = 1)
    # se_distss_prob_in_one_condition = std_distss_prob_in_one_condition / math.sqrt(len(validpath))

    # print(len(ave_distss_prob_in_one_condition))

    x = [x_lim/bins * (0.5+i) for i in range(len(distss_prob_in_one_condition))]
    distss_prob_in_one_condition = [distss_prob / (x_lim / bins) for distss_prob in distss_prob_in_one_condition]
    fig =plt.figure(figsize=(8,6),  dpi=1200)
    # lower_bound = [max(0, ave_distss_prob_in_one_condition[i] - se_distss_prob_in_one_condition[i]) for i in range(len(ave_distss_prob_in_one_condition))]
    # upper_bound = [ave_distss_prob_in_one_condition[i] + se_distss_prob_in_one_condition[i] for i in range(len(ave_distss_prob_in_one_condition))]
    plt.plot(x, distss_prob_in_one_condition)
    # plt.fill_between(x, lower_bound, upper_bound, alpha = 0.3)

    plt.xlabel('Distance to membrane (nm)')
    plt.ylabel('Probability density')
    plt.title(f'{condition}_condition_distance_density_distribution bin = {bins}\n  {[i for i in validpath]}')
    plt.tight_layout()


    # plt.legend()

    # plt.savefig(f'{sign}_condition_distance_density_distribution.pdf')
    # plt.show()
    # plt.close()

def dist2mem_condition_plot_notitle(distss_prob_in_one_condition, condition, validpath, bins = 100, x_lim = 1200):


    # ave_distss_prob_in_one_condition = np.average(distss_prob_in_one_condition, axis=0)
    # std_distss_prob_in_one_condition = np.std(distss_prob_in_one_condition, axis=0, ddof = 1)
    # se_distss_prob_in_one_condition = std_distss_prob_in_one_condition / math.sqrt(len(validpath))

    # print(len(ave_distss_prob_in_one_condition))

    x = [x_lim/bins * (0.5+i) for i in range(len(distss_prob_in_one_condition))]
    distss_prob_in_one_condition = [distss_prob / (x_lim / bins) for distss_prob in distss_prob_in_one_condition]
    fig =plt.figure(figsize=(8,6),  dpi=1200)
    # lower_bound = [max(0, ave_distss_prob_in_one_condition[i] - se_distss_prob_in_one_condition[i]) for i in range(len(ave_distss_prob_in_one_condition))]
    # upper_bound = [ave_distss_prob_in_one_condition[i] + se_distss_prob_in_one_condition[i] for i in range(len(ave_distss_prob_in_one_condition))]
    plt.plot(x, distss_prob_in_one_condition)
    # plt.fill_between(x, lower_bound, upper_bound, alpha = 0.3)

    plt.xlabel('Distance to membrane (nm)')
    plt.ylabel('Probability density')
    # plt.title(f'{condition}_condition_distance_density_distribution bin = {bins}\n  {[i for i in validpath]}')
    plt.tight_layout()






def dist2mem_condition_plot2(dist, file_idx, validpath, bins = 100, x_lim = 1200 ):
    fig =plt.figure(figsize=(8,6),  dpi=1200)
    # sns.displot(dist, kde=True, bins=bins ) #,rug=True, bins=bins, )
    sns.distplot(a=dist, kde=True, bins=bins,) #rug=True, bins=bins, )
    # sns.kdeplot(data=dist, shade=True)
    plt.title(f'{file_idx} condition actin to mem distance distribution  \n {[i for i in validpath]}' )
    plt.ylabel('probability density')
    plt.xlabel('actin - membrane distance (nm)')
    plt.xlim(0, x_lim)


def dist2mem_condition_plot2_notitle(dist, file_idx, validpath, bins = 100, x_lim = 1200 ):
    fig =plt.figure(figsize=(8,6),  dpi=1200)
    # sns.displot(dist, kde=True, bins=bins ) #,rug=True, bins=bins, )
    sns.distplot(a=dist, kde=True, bins=bins,) #rug=True, bins=bins, )
    # sns.kdeplot(data=dist, shade=True)
    # plt.title(f'{file_idx} condition actin to mem distance distribution  \n {[i for i in validpath]}' )
    plt.ylabel('probability density')
    plt.xlabel('actin - membrane distance (nm)')
    plt.xlim(0, x_lim)


def dist2mitomem_condition_plot(distss_prob_in_one_condition, condition, validpath, bins = 100, x_lim = 1200):

    x = [x_lim/bins * (0.5+i) for i in range(len(distss_prob_in_one_condition))]
    distss_prob_in_one_condition = [distss_prob / (x_lim / bins) for distss_prob in distss_prob_in_one_condition]
    fig =plt.figure(figsize=(8,6),  dpi=1200)
    # lower_bound = [max(0, ave_distss_prob_in_one_condition[i] - se_distss_prob_in_one_condition[i]) for i in range(len(ave_distss_prob_in_one_condition))]
    # upper_bound = [ave_distss_prob_in_one_condition[i] + se_distss_prob_in_one_condition[i] for i in range(len(ave_distss_prob_in_one_condition))]
    plt.plot(x, distss_prob_in_one_condition)
    # plt.fill_between(x, lower_bound, upper_bound, alpha = 0.3)

    plt.xlabel('Distance to membrane (nm)')
    plt.ylabel('Probability density')
    plt.title(f'{condition}_condition_actin distance to mito _density_distribution bin = {bins}\n  {[i for i in validpath]}')
    plt.tight_layout()



def dist2mitomem_condition_plot2(dist, file_idx, validpath, bins = 100, x_lim = 1200 ):
    fig =plt.figure(figsize=(8,6),  dpi=1200)
    # sns.displot(dist, kde=True, bins=bins ) #,rug=True, bins=bins, )
    sns.distplot(a=dist, kde=True, bins=bins,) #rug=True, bins=bins, )
    # sns.kdeplot(data=dist, shade=True)
    plt.title(f'{file_idx} condition actin to mito mem distance distribution  \n {[i for i in validpath]}' )
    plt.ylabel('probability density')
    plt.xlabel('actin - membrane distance (nm)')
    plt.xlim(0, x_lim)

def dist2ermem_condition_plot(distss_prob_in_one_condition, condition, validpath, bins = 100, x_lim = 1200):

    x = [x_lim/bins * (0.5+i) for i in range(len(distss_prob_in_one_condition))]
    distss_prob_in_one_condition = [distss_prob / (x_lim / bins) for distss_prob in distss_prob_in_one_condition]
    fig =plt.figure(figsize=(8,6),  dpi=1200)
    # lower_bound = [max(0, ave_distss_prob_in_one_condition[i] - se_distss_prob_in_one_condition[i]) for i in range(len(ave_distss_prob_in_one_condition))]
    # upper_bound = [ave_distss_prob_in_one_condition[i] + se_distss_prob_in_one_condition[i] for i in range(len(ave_distss_prob_in_one_condition))]
    plt.plot(x, distss_prob_in_one_condition)
    # plt.fill_between(x, lower_bound, upper_bound, alpha = 0.3)

    plt.xlabel('Distance to membrane (nm)')
    plt.ylabel('Probability density')
    plt.title(f'{condition}_condition_actin distance to er _density_distribution bin = {bins}\n  {[i for i in validpath]}')
    plt.tight_layout()



def dist2ermem_condition_plot2(dist, file_idx, validpath, bins = 100, x_lim = 1200 ):
    fig =plt.figure(figsize=(8,6),  dpi=1200)
    # sns.displot(dist, kde=True, bins=bins ) #,rug=True, bins=bins, )
    sns.distplot(a=dist, kde=True, bins=bins,) #rug=True, bins=bins, )
    # sns.kdeplot(data=dist, shade=True)
    plt.title(f'{file_idx} condition actin to er mem distance distribution  \n {[i for i in validpath]}' )
    plt.ylabel('probability density')
    plt.xlabel('actin - membrane distance (nm)')
    plt.xlim(0, x_lim)


def dist2fila_condition_plot(actin_vect_map_condition, condition, validpath):
    
    points_x = []
    points_y = []
    # points = []
    for actin_vect_map_lst in actin_vect_map_condition:
        for idx, actin_vects in enumerate(actin_vect_map_lst):
            # print(vis_img.shape)
                for vect in actin_vects:
                    if len(actin_vects) >0:
                        points_x.append(vect[0])
                        points_y.append(vect[1])
                    # points.append([vect[0],vect[1]])

    bins_ = 96  ## works
    fig = plt.figure(figsize=(8,6), dpi=1200)
    hist, xedges, yedges = np.histogram2d(points_x, points_y, bins=bins_, range =[[-500,500],[-500,500]] )
    # print(hist)
    values = hist.reshape(-1, 1)
    points = []
    for i in range(hist.shape[0]):
        for j in range(hist.shape[1]):
            points.append([i,j])
    grid_x, grid_y = np.mgrid[1:hist.shape[0]:5000j, 1:hist.shape[1]:5000j]

    grid_z0 = griddata(points, values, (grid_x, grid_y), method='linear',fill_value=0)
    print(' griddata')
    # alphalst = [i * 0.1 for i in range(1,11)]
    # for i, alpha in enumerate(alphalst):
    #     grid_z0_test = copy.deepcopy(grid_z0)
    #     grid_z0_test[np.where(grid_z0_test < (np.max(grid_z0_test) * alpha) )] = np.nan
    #     plt.imshow(grid_z0_test, extent=(-100,100,-100,100),cmap= 'RdYlBu_r')
    #     plt.title(f'{condition} Actin surrounding probability density ditribution\n alpha ={alpha},\n {validpath}  ')
    #     plt.colorbar()
    #     plt.savefig(f'{condition}_condition_actin_surrounding_alpha={alpha}.png')
    #     plt.show()
    #     plt.close()
    

    grid_z0[np.where(grid_z0 < (np.max(grid_z0) * 0.3) )] = np.nan
    grid_z0[np.where(grid_z0 < (np.max(grid_z0) * 0.45) )] = 0

    plt.imshow(grid_z0, extent=(-50,50,-50,50),cmap= 'RdYlBu_r')
    plt.title(f'{condition} Actin surrounding probability density ditribution \n {validpath}')
    
    # plt.grid(True)
    plt.colorbar()
    #plt.savefig()
    #plt.savefig(f'{file_idx}_actin_surrounding_probability_density_ditribution.pdf')
    # plt.show()



def dist2fila_condition_plot_old(actin_vect_map_condition, condition, validpath):
    
    vis_img = np.zeros([201, 201])
    for actin_vect_map_lst in actin_vect_map_condition:
        for idx, actin_vects in enumerate(actin_vect_map_lst):
            # print(vis_img.shape)
            if len(actin_vects) >0:
                for vect in actin_vects:
                    temp_vect = np.around(vect)
                    temp_vect = temp_vect / np.array([10,10,10])
                    temp_vect = [int(temp_vect[0]+vis_img.shape[0]/2), int(temp_vect[1]+100)]
                    vis_img[temp_vect[0]][temp_vect[1]] = vis_img[temp_vect[0]][temp_vect[1]] + 1
                    # vis_img[temp_vect[0]+1][temp_vect[1]] = vis_img[temp_vect[0]][temp_vect[1]] + 1
                    # vis_img[temp_vect[0]-1][temp_vect[1]] = vis_img[temp_vect[0]][temp_vect[1]] + 1
                    # vis_img[temp_vect[0]][temp_vect[1]+1] = vis_img[temp_vect[0]][temp_vect[1]] + 1
                    # vis_img[temp_vect[0]][temp_vect[1]-1] = vis_img[temp_vect[0]][temp_vect[1]] + 1

    print(np.max(vis_img))
    vis_img[np.where(vis_img == 0) ] = np.nan
    # vis_img = scipy.ndimage.binary_dilation(vis_img,iterations=3).astype(int)

    # vis_img2 = copy.deepcopy(vis_img)
    # vis_img2[np.where(vis_img2 < (np.max(vis_img2) * 0.3) )] = 0
    # vis_img_filtered = scipy.ndimage.gaussian_filter(vis_img2, sigma=0.5)
    
    # print('filteredmax', np.max(vis_img_filtered))

    fig = plt.figure(figsize=(8,6), dpi=1200)
    plt.imshow(vis_img, cmap=  'RdYlBu_r'  ) #'gray_r')  #RdYlBu
    plt.title(f'{condition} Actin surrounding probability density ditribution \n {validpath}')
    
    # plt.axis('off')
    plt.colorbar()
    # plt.savefig(f'{file_idx}_actin_surrounding_probability_density_ditribution.pdf')


def dist_fila2MT_condition_plot(actin_vect_map_condition, condition, validpath):
    '''
    Description: this function is used to filter the microtube_surround_actin data, yielding the filtered figure
    '''
    points_x = []
    points_y = []
    # points = []
    for actin_vect_map_lst in actin_vect_map_condition:
        for idx, actin_vects in enumerate(actin_vect_map_lst):
            # print(vis_img.shape)
            if len(actin_vects) >0:
                for vect in actin_vects:
                    points_x.append(vect[0])
                    points_y.append(vect[1])
                    # points.append([vect[0],vect[1]])

    bins_ = 96  ## works
    fig = plt.figure(figsize=(8,6), dpi=1200)
    hist, xedges, yedges = np.histogram2d(points_x, points_y, bins=bins_, range =[[-500,500],[-500,500]] )
    # print(hist)
    values = hist.reshape(-1, 1)
    points = []
    for i in range(hist.shape[0]):
        for j in range(hist.shape[1]):
            points.append([i,j])
    grid_x, grid_y = np.mgrid[1:hist.shape[0]:5000j, 1:hist.shape[1]:5000j]

    grid_z0 = griddata(points, values, (grid_x, grid_y), method='linear',fill_value=0)
    print('linear griddata')
    # alphalst = [i * 0.1 for i in range(1,11)]
    # for i, alpha in enumerate(alphalst):
    #     grid_z0_test = copy.deepcopy(grid_z0)
    #     grid_z0_test[np.where(grid_z0_test < (np.max(grid_z0_test) * alpha) )] = np.nan
    #     plt.imshow(grid_z0_test, extent=(-100,100,-100,100),cmap= 'RdYlBu_r')
    #     plt.title(f'{condition} Actin surrounding probability density ditribution\n alpha ={alpha},\n {validpath}  ')
    #     plt.colorbar()
    #     plt.savefig(f'{condition}_condition_actin_surrounding_alpha={alpha}.png')
    #     plt.show()
    #     plt.close()
    

    grid_z0[np.where(grid_z0 < (np.max(grid_z0) * 0.02) )] = np.nan #0.015
    #grid_z0[np.where(grid_z0 < (np.max(grid_z0) * 0.2) )] = 0 #0.2

    plt.imshow(grid_z0, extent=(-50,50,-50,50),cmap= 'RdYlBu_r')
    plt.title(f'{condition} MT surrounding actin probability density ditribution \n {validpath}')
    
    # plt.grid(True)
    plt.colorbar()
     
    # plt.show()

    # plt.savefig(f'{file_idx}_actin_surrounding_probability_density_ditribution.pdf')



def dist_fila2MT_condition_plot_old(actin_vect_map_condition, condition, validpath):
    
    vis_img = np.zeros([201, 201])
    for actin_vect_map_lst in actin_vect_map_condition:
        for idx, actin_vects in enumerate(actin_vect_map_lst):
            # print(vis_img.shape)
            if len(actin_vects) >0:
                for vect in actin_vects:
                    temp_vect = np.around(vect)
                    temp_vect = temp_vect / np.array([10,10,10])
                    temp_vect = [int(temp_vect[0]+vis_img.shape[0]/2), int(temp_vect[1]+100)]
                    vis_img[temp_vect[0]][temp_vect[1]] = vis_img[temp_vect[0]][temp_vect[1]] + 1
                    # vis_img[temp_vect[0]+1][temp_vect[1]] = vis_img[temp_vect[0]][temp_vect[1]] + 1
                    # vis_img[temp_vect[0]-1][temp_vect[1]] = vis_img[temp_vect[0]][temp_vect[1]] + 1
                    # vis_img[temp_vect[0]][temp_vect[1]+1] = vis_img[temp_vect[0]][temp_vect[1]] + 1
                    # vis_img[temp_vect[0]][temp_vect[1]-1] = vis_img[temp_vect[0]][temp_vect[1]] + 1

    print(np.max(vis_img))
    # vis_img = scipy.ndimage.binary_dilation(vis_img,iterations=3).astype(int)

    fig = plt.figure(figsize=(8,6), dpi=1200)
    plt.imshow(vis_img, cmap=  'RdYlBu_r'  ) #'gray_r')  #RdYlBu
    plt.title(f'{condition} MT surrounding actin probability density ditribution \n {validpath}')
    
    # plt.axis('off')
    plt.colorbar()
    # plt.savefig(f'{file_idx}_actin_surrounding_probability_density_ditribution.pdf')





def vis_actindistangle_conditionplot(actins_min_d_lst, actins_angle_lst, ):
    if len(actins_angle_lst) >0:
        max_ = np.max(actins_min_d_lst)
    else: 
        max_ = 1200
    points_x = actins_min_d_lst
    points_y = actins_angle_lst

    if len(actins_min_d_lst)>300:    
        bins_ = [50,50]
    elif len(actins_min_d_lst)>100:
        bins_ = [50,50]
    else:
        bins_ = [10, 10]

    fig = plt.figure(figsize=(8,6), dpi=1200)
    hist, xedges, yedges = np.histogram2d(points_x, points_y, bins=bins_, range =[[0,np.ceil(max_/100)*100],[0,90]] )
    # print(hist)
    
    values = hist.reshape(-1, 1)
    points = []
    for i in range(hist.shape[0]):
        for j in range(hist.shape[1]):
            points.append([i,j])

    # points_x, points_y = np.mgrid[1:100:100j, 1:100:100j]
    # points = np.array([ [points_x[i], points_y[i]] for i in range(len(points_x)) ])

    grid_x, grid_y = np.mgrid[0:hist.shape[0]:1000j, 0:hist.shape[1]:1000j]


    grid_z0 = griddata(points, values, (grid_x, grid_y), method='linear',fill_value = 0)


    # grid_z0[np.where(grid_z0 < (np.max(grid_z0) * 0.4) )] = np.nan
    # grid_z0[np.where(grid_z0 < (np.max(grid_z0) * 0.45) )] = 0

    # print(grid_z0)
    plt.imshow(grid_z0,interpolation='bilinear', cmap= 'RdYlBu_r', origin='lower')
    # plt.imshow(grid_z0,interpolation='bilinear', 
    #             extent=( 0,max_, 0,90,),
    #             cmap= 'RdYlBu_r', 
    #             origin='lower')
    # plt.title(f'{file_idx} Actin dist VS angle \n actin num {actin_num}, points num {point_num}')
    plt.xticks(np.arange(0,1000,200),np.arange(0, int(np.ceil(max_/100)*100), int(np.ceil(max_/100)*10*2)))
    plt.yticks(np.arange(0,990,110),np.arange(0,90,10))
    plt.xlabel('Distance (nm)')
    plt.ylabel('Angle')
    # plt.grid(True)
    plt.colorbar()
    

def vis_actindistangle_conditionplot2(filaments_min_d_lst, filaments_angle_lst):

    zoomsize = 5 # zoom_the_figure_weimin

    if len(filaments_min_d_lst) > 0:
        max_ = np.max(filaments_min_d_lst)
    else:
        max_ = 1200
    points_x = filaments_min_d_lst
    points_y = filaments_angle_lst

    if len(filaments_min_d_lst)>300:    
        bins_ = [int(50),int(50)]
    elif len(filaments_min_d_lst)>100:
        bins_ = [int(50),int(50)]
    else:
        bins_ = [int(10), int(10)]

    # if len(filaments_min_d_lst)>300:    
    #     bins_ = [int(50/zoomsize),int(50/zoomsize)]
    # elif len(filaments_min_d_lst)>100:
    #     bins_ = [int(50/zoomsize),int(50/zoomsize)]
    # else:
    #     bins_ = [int(10/zoomsize), int(10/zoomsize)]




    fig = plt.figure(figsize=(8,6), dpi=1200)
    hist, xedges, yedges = np.histogram2d(points_x, points_y, bins=bins_, range =[[0,np.ceil(max_/100)*100/zoomsize],[0,90/zoomsize]] )
    # print(hist)
    

    plt.imshow(hist)

    
    values = hist.reshape(-1, 1)
    points = []
    for i in range(hist.shape[0]):
        for j in range(hist.shape[1]):
            points.append([i,j])

    # points_x, points_y = np.mgrid[1:100:100j, 1:100:100j]
    # points = np.array([ [points_x[i], points_y[i]] for i in range(len(points_x)) ])


    grid_x, grid_y = np.mgrid[0:hist.shape[0]-1:100j, 0:hist.shape[1]-1:100j]
    # x_stick = 100
    # y_stick = 100
    

    grid_z0 = griddata(points, values, (grid_x, grid_y), method='linear',fill_value = 0)


    # grid_z0[np.where(grid_z0 < (np.max(grid_z0) * 0.4) )] = np.nan
    # grid_z0[np.where(grid_z0 < (np.max(grid_z0) * 0.45) )] = 0

    # print(grid_z0)
    plt.imshow(grid_z0, cmap= 'RdYlBu_r', origin='lower')
    # plt.imshow(grid_z0,interpolation='bilinear', 
    #             extent=( 0,max_, 0,90,),
    #             cmap= 'RdYlBu_r', 
    #             origin='lower')
    # plt.title(f'{file_idx} Actin dist VS angle \n filament num {filament_num}, points num {point_num}')
    plt.xticks(np.arange(0,100+20,20),np.arange(0, int(np.ceil(max_/100)*100/zoomsize)+int(np.ceil(max_/100)*10*2/zoomsize), int(np.ceil(max_/100)*10*2/zoomsize)))
    plt.yticks(np.arange(0,100+10,10),np.arange(0,100/zoomsize+10/zoomsize,10/zoomsize))
    plt.xlabel('Distance (nm)')
    plt.ylabel('Angle')
    # plt.grid(True)
    plt.colorbar()

