#%%
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
import cv2
#%%


# imagepath = f'F:\\salilab\\projects\\auto-segmentation_Aneesh\\model_for_insulin_vesicle\\Pytorch-UNet-master\\data\\imgs\\0cdf5b5d0ce1_02.jpg'
imagepath = f'F:\\salilab\\projects\\auto-segmentation_Aneesh\\model_for_insulin_vesicle\\Pytorch-UNet\\data\\\imgs\\766_8_211_img.jpg'

aa = Image.open(imagepath)
print(np.array(aa).shape)
print(np.array(aa))

img =  Image.open(imagepath)
print(img)
print(np.max(img))
print(img.size)
#%%
img_gray = img.convert("L")
img = np.array(img)[:,:,:3]
img_gray = np.asarray(img_gray)
img[:,:,0], img[:,:,1], img[:,:,2] = img_gray,img_gray,img_gray
print(img.shape)
#%%
print(img)
abcc = Image.fromarray(img)
plt.imshow(abcc)
plt.show()
#%%
aa = aa.convert("L")
bb = np.array(aa)
print(np.array(aa))
print(bb.shape)
plt.imshow(aa)
plt.show()


image = cv2.imread(imagepath, 0)
print(image.shape)
print(image)
# %%
aaa = np.empty([5,5])
bb = Image.fromarray(aaa)
print(np.array(bb))
plt.imshow(bb)
plt.show()
# %%
