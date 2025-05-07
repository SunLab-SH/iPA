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
# file = 'I:\\Bing\\fluorescence\\Results\\SIM\\0min_3_5_SIM_PM_NE_mask.npy'
os.chdir("I:\\Bing\\fluorescence\\Results\\SIM\\")   

os.getcwd() 

seq = '0min_3_5_SIM'
#seq = 'HighGlu_t20'
#%%
# Grab some test data.
#X, Y, Z = axes3d.get_test_data(0.05)
#print(Z)
raw_mask = np.load(seq + '_PM_NE_mask.npy') # T, z, C(PM/NE), x, y


k1 = np.loadtxt(str(seq) + '_ne_index.xvg')
xs = [item[0] for item in k1]
ys = [item[1] for item in k1]
zs = [item[2] for item in k1]
k2 = np.loadtxt(str(seq) + '_pm_index.xvg')
xs2 = [item[0] for item in k2]
ys2 = [item[1] for item in k2]
zs2 = [item[2] for item in k2]

#print(xs)

sep_tiff01 = np.ones((raw_mask.shape[0], raw_mask.shape[2], raw_mask.shape[3], 3))
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

ax1.imshow(sep_tiff01[1, :, :])
#ax1.imshow(sep_tiff01[300, :, :])
#ax1.imshow(sep_tiff02[232, :, :])
#sys.exit()

# Plot a basic wireframe.
#ax.scatter(xs, ys, zs, c='k', marker='.')
#ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)

plt.show()
fig.savefig('check_' + str(seq) + '_ne.png', transparent=True, dpi=300)

# %%
