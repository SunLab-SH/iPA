#%%

import numpy as np
import mrcfile
import matplotlib.pyplot as plt
import sys, os
import tifffile
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from scipy.ndimage import gaussian_filter
from scipy.spatial.distance import cdist
import copy
import concurrent.futures


dirpath = os.path.dirname(os.path.abspath(__file__))
os.chdir(dirpath)
sys.path.append(os.pardir)
print('curpath', os.getcwd())

from isg_seg_on_sim_experiment.parser import args

#%%
# load data from mrc file

class ImageProcessor:
    def __init__(self, path):
        self.path = path

    def load_mrc(self, path):
        with mrcfile.open(path) as mrc:
            data = mrc.data
        return data

    def normalize(self, image):
        min_val = np.min(image)
        max_val = np.max(image)
        print(min_val, max_val)
        return (image - min_val) / (max_val - min_val)


    def __call__(self, rescale = True):
        # Load the image
        image = self.load_mrc(self.path)
        image = image[:, ::-1, :] # to match images on ImageJ
        # 使得image的value为绝对值
        image = np.abs(image)
        image = self.normalize(image)
        print(np.min(image), np.max(image))
        self.original_shape = image.shape
        print("Original shape:", self.original_shape)

        if not rescale:
            return image
        else:# return rescaled_image
            # Calculate the current shape
            z_dim, x_dim, y_dim = image.shape

            # assert image.shape[1] == image.shape[2]
            new_shape = [image.shape[1], image.shape[1], image.shape[2]]

            # Calculate the ratio of old dimension to new dimension
            # [z,y,x] 31.3，31.3，100.3nm voxel representation
            z_ratio = 3.2 # new_shape[0] / z_dim 
            x_ratio = new_shape[1] / x_dim
            y_ratio = new_shape[2] / y_dim

            # Perform the interpolation
            rescaled_image = ndi.interpolation.zoom(image, (z_ratio, x_ratio, y_ratio), order=1)

            print("New shape:", rescaled_image.shape)  # should print (800, 800, 800)
            return rescaled_image


def check_mrc(img): 
    # check img py plot slide by slide
    for i in range(0,img.shape[0],16):
        plt.imshow(img[i,:,:])
        plt.title('slice {}, intensity {}'.format(i, [np.min(img[i]), np.max(img[i])]))
        plt.colorbar()
        plt.show()
    


# %%

def seperate_mask(label_mask, mrc, seqlst):
    labellst = seqlst

    def return_maskrange(label):        
        coords = np.where(label_mask == label)
        shape_ = label_mask.shape

        x_min, x_max = np.clip([np.min(coords[0])-5, np.max(coords[0])+6], 0, shape_[0])
        y_min, y_max = np.clip([np.min(coords[1])-5, np.max(coords[1])+6], 0, shape_[1])
        z_min, z_max = np.clip([np.min(coords[2])-5, np.max(coords[2])+6], 0, shape_[2])

        return (slice(x_min, x_max), slice(y_min, y_max), slice(z_min, z_max))

    rangelst = [return_maskrange(label) for label in labellst]

    def return_tinymask(idx):
        label = labellst[idx]
        temp_mask = np.zeros_like(label_mask)
        temp_mask[label_mask == label] = label
        return [temp_mask[rangelst[idx]], mrc[rangelst[idx]]]
        
    masklst = [return_tinymask(i) for i in range(len(labellst))]

    return masklst, rangelst



def discrete_sphere(max_r, min_r):
    # Create a cube of side length 2*max_r + 3
    cube = np.ones((int(2*max_r + 3),)*3)

    # Get the coordinates of all points in the cube
    coords = np.array(np.where(cube)).T

    # Compute the distance of each point from the center of the cube
    center = np.array([max_r + 1, max_r + 1, max_r + 1])
    dists = np.linalg.norm(coords - center, axis=1)

    # Get the unique distances and their counts (i.e., the volumes)
    r_lst, counts = np.unique(dists, return_counts=True)

    # Filter out distances less than min_r
    valid_indices = r_lst > min_r
    r_lst = r_lst[valid_indices]
    r_volume_lst = counts[valid_indices]

    return np.array([r_lst, r_volume_lst])




def get_rad(centcoord, mask, r_seq, theta):
    all_coords = np.array(np.where(mask != 0)).T
    dist = cdist(np.array([centcoord]), all_coords)[0]
    
    r_lst, r_volume = r_seq

    # Compute mask volumes for each radius
    mask_volumes = np.array([np.sum(dist <= r) for r in r_lst])

    # Compute ratio of mask volume to theoretical volume
    ratiolst = mask_volumes / r_volume

    # Find the first radius where the ratio exceeds the threshold
    r = next((r for r, ratio in zip(r_lst, ratiolst) if ratio >= theta), 0)

    return r

def fit_spheres(blobs, edt_values, mask, r_seq, theta):
    final_cent_coords = []
    final_cent_rads = []
    iterated_mask = mask.copy()
    
    for i, blob in enumerate(blobs):
        edt_value = edt_values[i]
        
        # Skip if edt_value is 0 or blob is within any existing sphere
        if edt_value == 0:
            continue
        if final_cent_coords and any(np.linalg.norm(np.array(final_cent_coords) - blob, axis=1) < final_cent_rads):
            continue
        
        # Compute radius
        r = get_rad(blob, mask, r_seq, theta)
        
        if r > 0:
            final_cent_coords.append(blob)
            final_cent_rads.append(r)
            
            # Update iterated_mask
            temp_edt = np.ones_like(iterated_mask)
            temp_edt[blob[0], blob[1], blob[2]] = 0
            temp_edt = ndi.distance_transform_edt(temp_edt)
            iterated_mask[temp_edt < r] = 0
    
    return final_cent_coords, final_cent_rads


def assign_labels(mask, final_cent_coords, final_cent_rads):
    instance_mask = np.zeros_like(mask)
    
    if not final_cent_coords:
        return instance_mask

    coords = np.array(np.where(mask > 0)).T

    def set_label(coord):
        dists = cdist([coord], final_cent_coords)[0]
        ratios = dists / final_cent_rads
        label = np.argmin(ratios)
        return label + 1

    labels = list(map(set_label, coords))
    
    for coord, label in zip(coords, labels):
        instance_mask[coord[0], coord[1], coord[2]] = label

    return instance_mask



def blob_fit(mask, mrc, min_r=1.5, theta=0.8, check_=False):
    mask = (mask != 0).astype(int)
    n = mrc.shape[0] // 2

    # Find local maxima in the filtered images
    sigmas = np.arange(1, 10, 1)
    blobs = [peak_local_max(gaussian_filter(mrc, sigma),
                            footprint=np.ones((3,) * mask.ndim),
                            exclude_border=True)
             for sigma in sigmas]
    blobs = np.concatenate(blobs)


    if check_:
        fig = plt.figure(figsize=(18,12))
        ax = plt.subplot(111)
        ax.imshow(mrc[n])
        ax.axis('off')
        blobs_lst = list(blobs)
        for num_ in range(len(blobs_lst)):
            circle1 = plt.Circle((blobs_lst[num_][2], blobs_lst[num_][1]), 0.7, color = 'r', linewidth=2, fill = False )
            # plt.gcf().gca().add_artist(circle1) 
            plt.gcf().gca().add_patch(circle1) 
            # add_patch  


    if not blobs.size:
        return mask if mask.sum() > 15 else np.zeros_like(mask)

    # Compute intensity and EDT of each blob
    intensities = mrc[tuple(blobs.T)]
    edt = ndi.distance_transform_edt(ndi.binary_dilation(mask))
    edt_values = edt[tuple(blobs.T)]

    # Sort blobs by intensity
    order = np.argsort(intensities)[::-1]
    blobs, intensities, edt_values = blobs[order], intensities[order], edt_values[order]

    # Compute radii of spheres for each possible radius
    max_r = edt_values.max() // 1 + 3
    r_seq = discrete_sphere(max_r, min_r)


    # Fit spheres to blobs
    final_cent_coords, final_cent_rads = fit_spheres(blobs, edt_values, mask, r_seq, theta)



    # Assign labels to mask
    instance_mask = assign_labels(mask, final_cent_coords, final_cent_rads)

    return instance_mask

#%%
def rescale2origion(image, original_shape):
    cur_shape = image.shape
    zoom_factors = [original_shape[0]/cur_shape[0], original_shape[1]/cur_shape[1], original_shape[2]/cur_shape[2]]
    rescaled_image = ndi.zoom(image, zoom_factors, order=1)
    print(rescaled_image.shape)
    return rescaled_image

#%%


def process_file(filename):
    print('processing {}'.format(filename))
    # load data
    processor = ImageProcessor(os.path.join(args.datapath, filename))
    interpolated_mrc = processor(rescale = True)
    original_shape = processor.original_shape

    # load mask
    interpolated_mask = np.zeros_like(interpolated_mrc)
    interpolated_mask[np.where(interpolated_mrc > 0)] = 1

    mrc, mask = interpolated_mrc, interpolated_mask
    label_mask,_ = ndi.label(mask)
    index_seq = [i for i in range(1, np.max(label_mask)+1)]
    print(index_seq[-1])

    batch_accumulate_seq = 0
    instance_mask = np.zeros_like(mask)

    for index, seqlst in enumerate(index_seq):
        masklst, rangelst = seperate_mask(label_mask, mrc,  [seqlst])
        idx_lst = [i for i in range(len(masklst))]

        def set_fitmask(idx):
            test_mask, test_mrc = masklst[idx][0], masklst[idx][1]
            instance_tinymask = blob_fit(test_mask, test_mrc, check_=False)
            return instance_tinymask

        instance_tinymask_lst = list(map(set_fitmask, idx_lst)) 
        label_seq = [ np.max(mask) for mask in instance_tinymask_lst ]
        accumulate_seq = [np.sum(label_seq[:i]) + batch_accumulate_seq for i in range(len(label_seq))]

        for idx in idx_lst:
            temp_mask = np.zeros_like(mask)
            accumu_label_mask = np.zeros_like(instance_tinymask_lst[idx])
            accumu_label_mask[np.where(instance_tinymask_lst[idx] != 0)] = accumulate_seq[idx]
            temp_mask[rangelst[idx]] = instance_tinymask_lst[idx] + accumu_label_mask
            instance_mask += temp_mask
            del temp_mask, accumu_label_mask
        
        batch_accumulate_seq += np.sum(label_seq)

        del masklst, rangelst, instance_tinymask_lst, label_seq, accumulate_seq
    
    instance_mask = rescale2origion(instance_mask, original_shape)
    instance_mask = instance_mask.astype(np.int16)
    tifffile.imsave(f'{args.resultpath}/instance_{filename}', instance_mask)
    plt.imshow(instance_mask[int(instance_mask.shape[0]/2)])
    plt.title('file {}, slice {}'.format(filename, int(instance_mask.shape[0]/2)))
    plt.show()
    plt.close()
#%%
file_list=[]
with open(f'../file_list_MT_ex4.txt','r') as f:     # Modified
	for line in f:
		file_list.append(line.strip('\n'))

# mrc_files = []
# for root, dirs, files in os.walk(args.datapath):
#     for file in files:
#         if os.path.splitext(file)[1] == '.mrc':
#             # mrc_files.append(os.path.join(root, file))
#             mrc_files.append(file)

# print(mrc_files)
# mrc_files = ['0+10-1-1_SIM_volumn_bin2_ISG.mrc'] # for test

#%%
mrc_files = []
for file_name in file_list:
    mrc_file = f'{file_name}_volumn_bin2_ISG_threshold.mrc'
    print(mrc_file)
    mrc_files.append(mrc_file) 

# Use a ThreadPoolExecutor to process the files in parallel
with concurrent.futures.ThreadPoolExecutor(max_workers=30) as executor:
    executor.map(process_file, mrc_files)

#%%