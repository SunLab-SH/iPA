'''
=================
3D wireframe plot
=================

A very basic demonstration of a wireframe plot.
'''
#%%
from mpl_toolkits.mplot3d import axes3d
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys, os

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#%%
os.chdir("F:\\Bing\\Features\\Results\\")   

os.getcwd() 

seq = '2.8-30_P2_PM_NE_mask'
#seq = 'HighGlu_t20'
#%%
# Grab some test data.
#X, Y, Z = axes3d.get_test_data(0.05)
#print(Z)
raw_mask = np.load('masks\\' + seq + '.npy') # T, z, C(PM/NE), x, y


k1 = np.loadtxt('01_index\\' + str(seq) + '_ne_index_t0.xvg')
xs = [item[0] for item in k1]
ys = [item[1] for item in k1]
zs = [item[2] for item in k1]
k2 = np.loadtxt('01_index\\' + str(seq) + '_pm_index_t0.xvg')
xs2 = [item[0] for item in k2]
ys2 = [item[1] for item in k2]
zs2 = [item[2] for item in k2]

#print(xs)

sep_tiff01 = np.ones((raw_mask.shape[1], raw_mask.shape[3], raw_mask.shape[4], 3))
#sep_tiff01 = np.zeros(500, 500, 500)
#[-19712, -13312, -13312]
#%%
for i in range(0, len(k1)):
    sep_tiff01[int(xs[i]),int(ys[i]),int(zs[i]),:] = 0

for i in range(0, len(k2)):
    sep_tiff01[int(xs2[i]),int(ys2[i]),int(zs2[i]),:] = 0.5
#%%
fig, (ax1) = plt.subplots(1,1, figsize=(7,7))

#cmap = matplotlib.cm.spring  # Can be any colormap that you want after the cm
#cmap.set_bad(color='white')

ax1.imshow(sep_tiff01[15, :, :])
#ax1.imshow(sep_tiff01[300, :, :])
#ax1.imshow(sep_tiff02[232, :, :])
#sys.exit()

# Plot a basic wireframe.
#ax.scatter(xs, ys, zs, c='k', marker='.')
#ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)

plt.show()
fig.savefig('01_index\\' + 'check_' + str(seq) + 't0_ne.png', transparent=True, dpi=300)

# %%
