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

os.chdir("i:\\Bing\\fluorescence\\Results\\")   

os.getcwd() 

# Grab some test data.
#X, Y, Z = axes3d.get_test_data(0.05)
#print(Z)

# seqs = ["2.8-30_P3S3-1", "2.8-30_P3S3-2", 
#         "2.8-30_P17S1-1", "16.7-5_P4", 
#         "16.7-30_P14S1-2", "16.7-30_P22"]

seqs = ['2.8-30_P3S3-1',
        '2.8-30_P2',
        '2.8-30_P3S3-2',
        '2.8-30_P10S2-1',
        '2.8-30_P10S2-2',
        '2.8-30_P16',
        '2.8-30_P17S1-1',
        '2.8-30_test',
        '16.7-5_P4',
        '16.7-5_P7-1',
        '16.7-5_P7-2',
        '16.7-5_P10S3-1',
        '16.7-5_P10S3-2',
        '16.7-5_P20',
        '16.7-30_P14S1-1',
        '16.7-30_P14S1-2',
        '16.7-30_P14S1-3',
        '16.7-30_P21',
        '16.7-30_P22-1',
        '16.7-30_P22']
# condition = "LowGlu"
timeslice = "t8"

#%%
for idx in range(len(seqs)):
    seq = seqs[idx]

    ne = np.loadtxt(f'01_index\\{seq}_PM_NE_mask_ne_index_{timeslice}.xvg')
    xne = [item[0] for item in ne]
    yne = [item[1] for item in ne]
    zne = [item[2] for item in ne]
    pm = np.loadtxt(f'01_index\\{seq}_PM_NE_mask_pm_index_{timeslice}.xvg')
    xpm = [item[0] for item in pm]
    ypm = [item[1] for item in pm]
    zpm = [item[2] for item in pm]

    k0 = np.loadtxt(f'03_slice\\{seq}_slice_shell0_{timeslice}.xvg')
    xs0 = [item[0] for item in k0]
    ys0 = [item[1] for item in k0]
    zs0 = [item[2] for item in k0]
    k1 = np.loadtxt(f'03_slice\\{seq}_slice_shell1_{timeslice}.xvg')
    xs1 = [item[0] for item in k1]
    ys1 = [item[1] for item in k1]
    zs1 = [item[2] for item in k1]
    k2 = np.loadtxt(f'03_slice\\{seq}_slice_shell2_{timeslice}.xvg')
    xs2 = [item[0] for item in k2]
    ys2 = [item[1] for item in k2]
    zs2 = [item[2] for item in k2]
    k3 = np.loadtxt(f'03_slice\\{seq}_slice_shell3_{timeslice}.xvg')
    xs3 = [item[0] for item in k3]
    ys3 = [item[1] for item in k3]
    zs3 = [item[2] for item in k3]
    k4 = np.loadtxt(f'03_slice\\{seq}_slice_shell4_{timeslice}.xvg')
    xs4 = [item[0] for item in k4]
    ys4 = [item[1] for item in k4]
    zs4 = [item[2] for item in k4]
    k5 = np.loadtxt(f'03_slice\\{seq}_slice_shell5_{timeslice}.xvg')
    xs5 = [item[0] for item in k5]
    ys5 = [item[1] for item in k5]
    zs5 = [item[2] for item in k5]
    k6 = np.loadtxt(f'03_slice\\{seq}_slice_shell6_{timeslice}.xvg')
    xs6 = [item[0] for item in k6]
    ys6 = [item[1] for item in k6]
    zs6 = [item[2] for item in k6]
    k7 = np.loadtxt(f'03_slice\\{seq}_slice_shell7_{timeslice}.xvg')
    xs7 = [item[0] for item in k7]
    ys7 = [item[1] for item in k7]
    zs7 = [item[2] for item in k7]
    k8 = np.loadtxt(f'03_slice\\{seq}_slice_shell8_{timeslice}.xvg')
    xs8 = [item[0] for item in k8]
    ys8 = [item[1] for item in k8]
    zs8 = [item[2] for item in k8]
    k9 = np.loadtxt(f'03_slice\\{seq}_slice_shell9_{timeslice}.xvg')
    xs9 = [item[0] for item in k9]
    ys9 = [item[1] for item in k9]
    zs9 = [item[2] for item in k9]
    k10 = np.loadtxt(f'03_slice\\{seq}_slice_shell10_{timeslice}.xvg')
    xs10 = [item[0] for item in k10]
    ys10 = [item[1] for item in k10]
    zs10 = [item[2] for item in k10]
    k11 = np.loadtxt(f'03_slice\\{seq}_slice_shell11_{timeslice}.xvg')
    xs11 = [item[0] for item in k11]
    ys11 = [item[1] for item in k11]
    zs11 = [item[2] for item in k11]
    k12 = np.loadtxt(f'03_slice\\{seq}_slice_shell12_{timeslice}.xvg')
    xs12 = [item[0] for item in k12]
    ys12 = [item[1] for item in k12]
    zs12 = [item[2] for item in k12]
    k13 = np.loadtxt(f'03_slice\\{seq}_slice_shell13_{timeslice}.xvg')
    xs13 = [item[0] for item in k13]
    ys13 = [item[1] for item in k13]
    zs13 = [item[2] for item in k13]
    k14 = np.loadtxt(f'03_slice\\{seq}_slice_shell14_{timeslice}.xvg')
    xs14 = [item[0] for item in k14]
    ys14 = [item[1] for item in k14]
    zs14 = [item[2] for item in k14]
    k15 = np.loadtxt(f'03_slice\\{seq}_slice_shell15_{timeslice}.xvg')
    xs15 = [item[0] for item in k15]
    ys15 = [item[1] for item in k15]
    zs15 = [item[2] for item in k15]
    k16 = np.loadtxt(f'03_slice\\{seq}_slice_shell16_{timeslice}.xvg')
    xs16 = [item[0] for item in k16]
    ys16 = [item[1] for item in k16]
    zs16 = [item[2] for item in k16]

    # kk3 = np.loadtxt("../../../vesicle/" + str(condition) + "_" + str(timeslice) + "_coord_voxel_all.xvg")
    # kxs3 = [item[0] for item in kk3]
    # kys3 = [item[1] for item in kk3]
    # kzs3 = [item[2] for item in kk3]


    # print(len(kk3))
    mask_t = np.load(f'masks\\{seq}_PM_NE_mask.npy')
    sep_tiff01 = np.ones((mask_t.shape[1], mask_t.shape[3], mask_t.shape[4],3))

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
        sep_tiff01[int(xs1[i]),int(ys1[i]),int(zs1[i]),1] = 0
        sep_tiff01[int(xs1[i]),int(ys1[i]),int(zs1[i]),2] = 1

    for i in range(0, len(k2)):
        sep_tiff01[int(xs2[i]),int(ys2[i]),int(zs2[i]),0] = 0
        sep_tiff01[int(xs2[i]),int(ys2[i]),int(zs2[i]),1] = 0
        sep_tiff01[int(xs2[i]),int(ys2[i]),int(zs2[i]),2] = 0

    for i in range(0, len(k3)):
        sep_tiff01[int(xs3[i]),int(ys3[i]),int(zs3[i]),0] = 0
        sep_tiff01[int(xs3[i]),int(ys3[i]),int(zs3[i]),1] = 0
        sep_tiff01[int(xs3[i]),int(ys3[i]),int(zs3[i]),2] = 1

    for i in range(0, len(k4)):
        sep_tiff01[int(xs4[i]),int(ys4[i]),int(zs4[i]),0] = 0
        sep_tiff01[int(xs4[i]),int(ys4[i]),int(zs4[i]),1] = 0
        sep_tiff01[int(xs4[i]),int(ys4[i]),int(zs4[i]),2] = 0

    for i in range(0, len(k5)):
        sep_tiff01[int(xs5[i]),int(ys5[i]),int(zs5[i]),0] = 0
        sep_tiff01[int(xs5[i]),int(ys5[i]),int(zs5[i]),1] = 0
        sep_tiff01[int(xs5[i]),int(ys5[i]),int(zs5[i]),2] = 1

    for i in range(0, len(k6)):
        sep_tiff01[int(xs6[i]),int(ys6[i]),int(zs6[i]),0] = 0
        sep_tiff01[int(xs6[i]),int(ys6[i]),int(zs6[i]),1] = 0
        sep_tiff01[int(xs6[i]),int(ys6[i]),int(zs6[i]),2] = 0

    for i in range(0, len(k7)):
        sep_tiff01[int(xs7[i]),int(ys7[i]),int(zs7[i]),0] = 0
        sep_tiff01[int(xs7[i]),int(ys7[i]),int(zs7[i]),1] = 0
        sep_tiff01[int(xs7[i]),int(ys7[i]),int(zs7[i]),2] = 1

    for i in range(0, len(k8)):
        sep_tiff01[int(xs8[i]),int(ys8[i]),int(zs8[i]),0] = 0
        sep_tiff01[int(xs8[i]),int(ys8[i]),int(zs8[i]),1] = 0
        sep_tiff01[int(xs8[i]),int(ys8[i]),int(zs8[i]),2] = 0

    for i in range(0, len(k9)):
        sep_tiff01[int(xs9[i]),int(ys9[i]),int(zs9[i]),0] = 0
        sep_tiff01[int(xs9[i]),int(ys9[i]),int(zs9[i]),1] = 0
        sep_tiff01[int(xs9[i]),int(ys9[i]),int(zs9[i]),2] = 1

    for i in range(0, len(k10)):
        sep_tiff01[int(xs10[i]),int(ys10[i]),int(zs10[i]),0] = 0
        sep_tiff01[int(xs10[i]),int(ys10[i]),int(zs10[i]),1] = 0
        sep_tiff01[int(xs10[i]),int(ys10[i]),int(zs10[i]),2] = 0

    for i in range(0, len(k11)):
        sep_tiff01[int(xs11[i]),int(ys11[i]),int(zs11[i]),0] = 0
        sep_tiff01[int(xs11[i]),int(ys11[i]),int(zs11[i]),1] = 0
        sep_tiff01[int(xs11[i]),int(ys11[i]),int(zs11[i]),2] = 1

    for i in range(0, len(k12)):
        sep_tiff01[int(xs12[i]),int(ys12[i]),int(zs12[i]),0] = 0
        sep_tiff01[int(xs12[i]),int(ys12[i]),int(zs12[i]),1] = 0
        sep_tiff01[int(xs12[i]),int(ys12[i]),int(zs12[i]),2] = 0

    for i in range(0, len(k13)):
        sep_tiff01[int(xs13[i]),int(ys13[i]),int(zs13[i]),0] = 0
        sep_tiff01[int(xs13[i]),int(ys13[i]),int(zs13[i]),1] = 0
        sep_tiff01[int(xs13[i]),int(ys13[i]),int(zs13[i]),2] = 1

    for i in range(0, len(k14)):
        sep_tiff01[int(xs14[i]),int(ys14[i]),int(zs14[i]),0] = 0
        sep_tiff01[int(xs14[i]),int(ys14[i]),int(zs14[i]),1] = 0
        sep_tiff01[int(xs14[i]),int(ys14[i]),int(zs14[i]),2] = 0

    for i in range(0, len(k15)):
        sep_tiff01[int(xs15[i]),int(ys15[i]),int(zs15[i]),0] = 0
        sep_tiff01[int(xs15[i]),int(ys15[i]),int(zs15[i]),1] = 0
        sep_tiff01[int(xs15[i]),int(ys15[i]),int(zs15[i]),2] = 1

    for i in range(0, len(k16)):
        sep_tiff01[int(xs16[i]),int(ys16[i]),int(zs16[i]),0] = 0
        sep_tiff01[int(xs16[i]),int(ys16[i]),int(zs16[i]),1] = 0
        sep_tiff01[int(xs16[i]),int(ys16[i]),int(zs16[i]),2] = 0

    for i in range(0, len(ne)):
        sep_tiff01[int(xne[i]),int(yne[i]),int(zne[i]),:] = [1,0,0]

    for i in range(0, len(pm)):
        sep_tiff01[int(xpm[i]),int(ypm[i]),int(zpm[i]),:] = [0,1,0]

    fig, (ax1) = plt.subplots(1,1, figsize=(7,7))
    #ax1.imshow(sep_tiff01[48, :, :])
    ax1.imshow(sep_tiff01[15, :, :])
    #ax1.imshow(sep_tiff02[232, :, :])
    #sys.exit()

    # Plot a basic wireframe.
    #ax.scatter(xs, ys, zs, c='k', marker='.')
    #ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)
    plt.title(seq)
    plt.show()
    fig.savefig("03_slice\\plot\\plot_"+ seq + "z15_t8_" + "16slice.png", transparent=True, dpi=300)

# %%
