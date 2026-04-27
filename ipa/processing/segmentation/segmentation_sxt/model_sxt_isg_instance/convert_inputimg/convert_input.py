#%%

import numpy as np

def convert_img(img_3d):
    img_3d_new = np.zeros_like(img_3d)

    for i in range(img_3d.shape[0]):
        curslice = img_3d[i]
        curslice = (curslice - np.min(curslice) )/ np.max(curslice) * 255
        img_3d_new[i] = curslice

    return img_3d_new.astype(np.uint8)

