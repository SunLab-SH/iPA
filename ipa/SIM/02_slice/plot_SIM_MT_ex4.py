'''
=================
3D wireframe plot
=================

A very basic demonstration of a wireframe plot.
'''
#%%
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import sys, os

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

os.chdir("i:\\Bing\\fluorescence\\Results\\SIM\\")   

os.getcwd() 

# Grab some test data.
#X, Y, Z = axes3d.get_test_data(0.05)
#print(Z)

# seqs = ["2.8-30_P3S3-1", "2.8-30_P3S3-2", 
#         "2.8-30_P17S1-1", "16.7-5_P4", 
#         "16.7-30_P14S1-2", "16.7-30_P22"]


# condition = "LowGlu"
# timeslice = "t8"

#%%
file_dir = 'I:\\Bing\\fluorescence\\3D\\00_index\\'
file_dir2 = 'I:\\Bing\\fluorescence\\3D\\02_slice\\'

file_list=[]
with open(f'../file_list.txt','r') as f:
	for line in f:
		file_list.append(line.strip('\n'))

for file_name in file_list:
    print(file_name)


    ne = np.loadtxt(f'{file_dir}{file_name}_bin2_ne_index.xvg')
    zne = [item[0] for item in ne]
    yne = [item[1] for item in ne]
    xne = [item[2] for item in ne]
    pm = np.loadtxt(f'{file_dir}{file_name}_bin2_pm_index.xvg')
    zpm = [item[0] for item in pm]
    ypm = [item[1] for item in pm]
    xpm = [item[2] for item in pm]

    k0 = np.loadtxt(f'{file_dir2}{file_name}_bin2_slice_shell0.xvg')
    zs0 = [item[0] for item in k0]
    ys0 = [item[1] for item in k0]
    xs0 = [item[2] for item in k0]
    k1 = np.loadtxt(f'{file_dir2}{file_name}_bin2_slice_shell1.xvg')
    zs1 = [item[0] for item in k1]
    ys1 = [item[1] for item in k1]
    xs1 = [item[2] for item in k1]
    k2 = np.loadtxt(f'{file_dir2}{file_name}_bin2_slice_shell2.xvg')
    zs2 = [item[0] for item in k2]
    ys2 = [item[1] for item in k2]
    xs2 = [item[2] for item in k2]
    k3 = np.loadtxt(f'{file_dir2}{file_name}_bin2_slice_shell3.xvg')
    zs3 = [item[0] for item in k3]
    ys3 = [item[1] for item in k3]
    xs3 = [item[2] for item in k3]
    k4 = np.loadtxt(f'{file_dir2}{file_name}_bin2_slice_shell4.xvg')
    zs4 = [item[0] for item in k4]
    ys4 = [item[1] for item in k4]
    xs4 = [item[2] for item in k4]
    k5 = np.loadtxt(f'{file_dir2}{file_name}_bin2_slice_shell5.xvg')
    zs5 = [item[0] for item in k5]
    ys5 = [item[1] for item in k5]
    xs5 = [item[2] for item in k5]
    k6 = np.loadtxt(f'{file_dir2}{file_name}_bin2_slice_shell6.xvg')
    zs6 = [item[0] for item in k6]
    ys6 = [item[1] for item in k6]
    xs6 = [item[2] for item in k6]
    k7 = np.loadtxt(f'{file_dir2}{file_name}_bin2_slice_shell7.xvg')
    zs7 = [item[0] for item in k7]
    ys7 = [item[1] for item in k7]
    xs7 = [item[2] for item in k7]
    k8 = np.loadtxt(f'{file_dir2}{file_name}_bin2_slice_shell8.xvg')
    zs8 = [item[0] for item in k8]
    ys8 = [item[1] for item in k8]
    xs8 = [item[2] for item in k8]
    # k9 = np.loadtxt(f'03_slice\\{seq}_slice_shell9.xvg')
    # xs9 = [item[0] for item in k9]
    # ys9 = [item[1] for item in k9]
    # zs9 = [item[2] for item in k9]
    # k10 = np.loadtxt(f'03_slice\\{seq}_slice_shell10.xvg')
    # xs10 = [item[0] for item in k10]
    # ys10 = [item[1] for item in k10]
    # zs10 = [item[2] for item in k10]
    # k11 = np.loadtxt(f'03_slice\\{seq}_slice_shell11.xvg')
    # xs11 = [item[0] for item in k11]
    # ys11 = [item[1] for item in k11]
    # zs11 = [item[2] for item in k11]
    # k12 = np.loadtxt(f'03_slice\\{seq}_slice_shell12.xvg')
    # xs12 = [item[0] for item in k12]
    # ys12 = [item[1] for item in k12]
    # zs12 = [item[2] for item in k12]
    # k13 = np.loadtxt(f'03_slice\\{seq}_slice_shell13.xvg')
    # xs13 = [item[0] for item in k13]
    # ys13 = [item[1] for item in k13]
    # zs13 = [item[2] for item in k13]
    # k14 = np.loadtxt(f'03_slice\\{seq}_slice_shell14.xvg')
    # xs14 = [item[0] for item in k14]
    # ys14 = [item[1] for item in k14]
    # zs14 = [item[2] for item in k14]
    # k15 = np.loadtxt(f'03_slice\\{seq}_slice_shell15.xvg')
    # xs15 = [item[0] for item in k15]
    # ys15 = [item[1] for item in k15]
    # zs15 = [item[2] for item in k15]
    # k16 = np.loadtxt(f'03_slice\\{seq}_slice_shell16.xvg')
    # xs16 = [item[0] for item in k16]
    # ys16 = [item[1] for item in k16]
    # zs16 = [item[2] for item in k16]

    # kk3 = np.loadtxt("../../../vesicle/" + str(condition) + "_" + str(timeslice) + "_coord_voxel_all.xvg")
    # kxs3 = [item[0] for item in kk3]
    # kys3 = [item[1] for item in kk3]
    # kzs3 = [item[2] for item in kk3]


    # print(len(kk3))
    sep_tiff01 = np.ones((800,800,64,3))

    # # rgb - darkorange [255/255, 140/255, 0]
    # for i in range(0, len(kk3)):
    #     sep_tiff01[int(kxs3[i]),int(kys3[i]),int(kzs3[i]),0] = 1
    #     sep_tiff01[int(kxs3[i]),int(kys3[i]),int(kzs3[i]),1] = 0
    #     sep_tiff01[int(kxs3[i]),int(kys3[i]),int(kzs3[i]),2] = 0

    # ne    
    for i in range(0, len(k0)):
        sep_tiff01[int(xs0[i]),int(ys0[i]),int(zs0[i]),:] = 0

    for i in range(0, len(k1)):
        sep_tiff01[int(xs1[i]),int(ys1[i]),int(zs1[i]),0] = 0
        sep_tiff01[int(xs1[i]),int(ys1[i]),int(zs1[i]),1] = 0.2
        sep_tiff01[int(xs1[i]),int(ys1[i]),int(zs1[i]),2] = 1

    for i in range(0, len(k2)):
        sep_tiff01[int(xs2[i]),int(ys2[i]),int(zs2[i]),0] = 0
        sep_tiff01[int(xs2[i]),int(ys2[i]),int(zs2[i]),1] = 0.1
        sep_tiff01[int(xs2[i]),int(ys2[i]),int(zs2[i]),2] = 0

    for i in range(0, len(k3)):
        sep_tiff01[int(xs3[i]),int(ys3[i]),int(zs3[i]),0] = 0
        sep_tiff01[int(xs3[i]),int(ys3[i]),int(zs3[i]),1] = 0.2
        sep_tiff01[int(xs3[i]),int(ys3[i]),int(zs3[i]),2] = 1

    for i in range(0, len(k4)):
        sep_tiff01[int(xs4[i]),int(ys4[i]),int(zs4[i]),0] = 0
        sep_tiff01[int(xs4[i]),int(ys4[i]),int(zs4[i]),1] = 0.1
        sep_tiff01[int(xs4[i]),int(ys4[i]),int(zs4[i]),2] = 0

    for i in range(0, len(k5)):
        sep_tiff01[int(xs5[i]),int(ys5[i]),int(zs5[i]),0] = 0
        sep_tiff01[int(xs5[i]),int(ys5[i]),int(zs5[i]),1] = 0.2
        sep_tiff01[int(xs5[i]),int(ys5[i]),int(zs5[i]),2] = 1

    for i in range(0, len(k6)):
        sep_tiff01[int(xs6[i]),int(ys6[i]),int(zs6[i]),0] = 0.1
        sep_tiff01[int(xs6[i]),int(ys6[i]),int(zs6[i]),1] = 0.1
        sep_tiff01[int(xs6[i]),int(ys6[i]),int(zs6[i]),2] = 0.1

    for i in range(0, len(k7)):
        sep_tiff01[int(xs7[i]),int(ys7[i]),int(zs7[i]),0] = 0
        sep_tiff01[int(xs7[i]),int(ys7[i]),int(zs7[i]),1] = 0.2
        sep_tiff01[int(xs7[i]),int(ys7[i]),int(zs7[i]),2] = 1

    for i in range(0, len(k8)):
        sep_tiff01[int(xs8[i]),int(ys8[i]),int(zs8[i]),0] = 0
        sep_tiff01[int(xs8[i]),int(ys8[i]),int(zs8[i]),1] = 1
        sep_tiff01[int(xs8[i]),int(ys8[i]),int(zs8[i]),2] = 0

    # for i in range(0, len(k9)):
    #     sep_tiff01[int(xs9[i]),int(ys9[i]),int(zs9[i]),0] = 0
    #     sep_tiff01[int(xs9[i]),int(ys9[i]),int(zs9[i]),1] = 0
    #     sep_tiff01[int(xs9[i]),int(ys9[i]),int(zs9[i]),2] = 1

    # for i in range(0, len(k10)):
    #     sep_tiff01[int(xs10[i]),int(ys10[i]),int(zs10[i]),0] = 0
    #     sep_tiff01[int(xs10[i]),int(ys10[i]),int(zs10[i]),1] = 0
    #     sep_tiff01[int(xs10[i]),int(ys10[i]),int(zs10[i]),2] = 0

    # for i in range(0, len(k11)):
    #     sep_tiff01[int(xs11[i]),int(ys11[i]),int(zs11[i]),0] = 0
    #     sep_tiff01[int(xs11[i]),int(ys11[i]),int(zs11[i]),1] = 0
    #     sep_tiff01[int(xs11[i]),int(ys11[i]),int(zs11[i]),2] = 1

    # for i in range(0, len(k12)):
    #     sep_tiff01[int(xs12[i]),int(ys12[i]),int(zs12[i]),0] = 0
    #     sep_tiff01[int(xs12[i]),int(ys12[i]),int(zs12[i]),1] = 0
    #     sep_tiff01[int(xs12[i]),int(ys12[i]),int(zs12[i]),2] = 0

    # for i in range(0, len(k13)):
    #     sep_tiff01[int(xs13[i]),int(ys13[i]),int(zs13[i]),0] = 0
    #     sep_tiff01[int(xs13[i]),int(ys13[i]),int(zs13[i]),1] = 0
    #     sep_tiff01[int(xs13[i]),int(ys13[i]),int(zs13[i]),2] = 1

    # for i in range(0, len(k14)):
    #     sep_tiff01[int(xs14[i]),int(ys14[i]),int(zs14[i]),0] = 0
    #     sep_tiff01[int(xs14[i]),int(ys14[i]),int(zs14[i]),1] = 0
    #     sep_tiff01[int(xs14[i]),int(ys14[i]),int(zs14[i]),2] = 0

    # for i in range(0, len(k15)):
    #     sep_tiff01[int(xs15[i]),int(ys15[i]),int(zs15[i]),0] = 0
    #     sep_tiff01[int(xs15[i]),int(ys15[i]),int(zs15[i]),1] = 0
    #     sep_tiff01[int(xs15[i]),int(ys15[i]),int(zs15[i]),2] = 1

    # for i in range(0, len(k16)):
    #     sep_tiff01[int(xs16[i]),int(ys16[i]),int(zs16[i]),0] = 0
    #     sep_tiff01[int(xs16[i]),int(ys16[i]),int(zs16[i]),1] = 0
    #     sep_tiff01[int(xs16[i]),int(ys16[i]),int(zs16[i]),2] = 0

    for i in range(0, len(ne)):
        sep_tiff01[int(xne[i]),int(yne[i]),int(zne[i]),:] = [254/255,67/255,101/255]

    for i in range(0, len(pm)):
        sep_tiff01[int(xpm[i]),int(ypm[i]),int(zpm[i]),:] = [0,0,0]

    fig, (ax1) = plt.subplots(1,1, figsize=(7,7))
    #ax1.imshow(sep_tiff01[48, :, :])
    ax1.imshow(sep_tiff01[:, :, 11])
    #ax1.imshow(sep_tiff02[232, :, :])
    #sys.exit()

    # Plot a basic wireframe.
    #ax.scatter(xs, ys, zs, c='k', marker='.')
    #ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)
    plt.title(file_name)
    plt.show()
    fig.savefig("plot_"+ file_name + "_bin2_z1_" + "8slice.png", transparent=True, dpi=300)

    break
# %%
