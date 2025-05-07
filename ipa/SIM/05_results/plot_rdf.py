#!/usr/bin/python

#%%
# Import module for plot "matplotlib"
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import rc
import matplotlib
import matplotlib.ticker as mtick
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import rcParams
from pylab import *
import pylab as pylab
import matplotlib.patches as patches
from matplotlib.pyplot import *

#%%
def poly5_func(xdata, pars):
    '''
    calculate a fifth order polynomial equation with xdata and a list of parameter 
    '''
    return pars[0]*xdata**5 +pars[1]*xdata**4 + pars[2]*xdata**3 + pars[3]*xdata**2 + pars[4]*xdata +pars[5]


seqs_basal = ['2.8-30_P3S3-1', '2.8-30_P3S3-2',
        '2.8-30_P10S2-1', '2.8-30_P10S2-2', '2.8-30_P16', 
        '2.8-30_P17S1-1', '2.8-30_test']
seqs_5 = ['16.7-5_P4', '16.7-5_P7-1', '16.7-5_P7-2', '16.7-5_P10S3-1', 
        '16.7-5_P10S3-2', '16.7-5_P20']
seqs_30 = ['16.7-30_P14S1-1', '16.7-30_P14S1-2', '16.7-30_P14S1-3', '16.7-30_P21',
        '16.7-30_P22', '16.7-30_P22-1'] # '2.8-30_P2', 

N_shell = 17

# cellnumber = '2.8-30_P3S3-1'

k1_sum = []
for seq in seqs_basal:
    # Read data file
    k1 = np.loadtxt('../../05_rdf/results/'+ str(seq) + '_isgrdf_ne-pm.xvg')
    k1_sum.append(k1)

k1_mean = np.mean(np.array(k1_sum), axis=0)
k1_std = np.std(np.array(k1_sum), axis=0)

k2_sum = []
for seq in seqs_5:
    # Read data file
    k2 = np.loadtxt('../../05_rdf/results/'+ str(seq) + '_isgrdf_ne-pm.xvg')
    k2_sum.append(k2)

k2_mean = np.mean(np.array(k2_sum), axis=0)
k2_std = np.std(np.array(k2_sum), axis=0)

k3_sum = []
for seq in seqs_30:
    # Read data file
    k3 = np.loadtxt('../../05_rdf/results/'+ str(seq) + '_isgrdf_ne-pm.xvg')
    k3_sum.append(k3)

k3_mean = np.mean(np.array(k3_sum), axis=0)
k3_std = np.std(np.array(k3_sum), axis=0)

Rcell = 53946
Rne = 34939
Risg = 1289
dr = (Rcell - Rne)/16
xdata = np.linspace(1, 16, len(k1_mean))*dr

'''
# Plots
'''

# Specifiy environmental parameter
rc('font',**{'family':'serif','serif':['Arial Narrow']})

# Create axes 
fig = plt.figure(figsize=(6,4)) #inches
fig.subplots_adjust(left=0.3, right=0.7, bottom=0.3, top=0.7)

# Main figure
ax1 = plt.subplot2grid((1,1), (0, 0))

#ax1.set_title("Plot title...")    
ax1.set_xlabel('distance ($\,$$\mu$m)',fontweight="normal",fontsize="8")
ax1.set_ylabel('g(r)',fontweight="normal",fontsize="8",color='k')
ax1.tick_params(direction='out', pad=1)
ax1.xaxis.set_label_coords(0.5, -0.15)
ax1.yaxis.set_label_coords(-0.14, 0.5)
ax1.set_xlim([0,2])
ax1.set_ylim([0,0.2])
xmajorLocator   = MultipleLocator(0.4)
xmajorFormatter = FormatStrFormatter('%.1f')
xminorLocator   = MultipleLocator(0.2)
ax1.xaxis.set_major_locator(xmajorLocator)
ax1.xaxis.set_major_formatter(xmajorFormatter)
ax1.xaxis.set_minor_locator(xminorLocator)
ymajorLocator   = MultipleLocator(0.1)
ymajorFormatter = FormatStrFormatter('%.1f')
yminorLocator   = MultipleLocator(0.05)
ax1.yaxis.set_major_locator(ymajorLocator)
ax1.yaxis.set_major_formatter(ymajorFormatter)
ax1.yaxis.set_minor_locator(yminorLocator)
#rcParams['xtick.direction'] = 'in'
#rcParams['ytick.direction'] = 'in'
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(8)
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(8)
    tick.label.set_color('k')
for axis in ['bottom','left']:ax1.spines[axis].set_linewidth(0.5)
for axis in ['top','right']:ax1.spines[axis].set_linewidth(0)
for line in ax1.xaxis.get_ticklines():
    line.set_markersize(2)
    line.set_markeredgewidth(0.5)
for line in ax1.yaxis.get_ticklines():
    line.set_markersize(2)
    line.set_markeredgewidth(0.5)
for line in ax1.xaxis.get_minorticklines():
    line.set_markersize(1)
    line.set_markeredgewidth(0.5)
for line in ax1.yaxis.get_minorticklines():
    line.set_markersize(1)
    line.set_markeredgewidth(0.5)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')

ax1.plot((0,100),(1/16,1/16),'k',linestyle="--",linewidth = 0.5)
ydata = k1_mean[:]
#yerr= k1[:,1]
for i in range(0,len(ydata)):
    ax1.plot(((xdata[i]-dr)/10000, xdata[i]/10000),(ydata[i],ydata[i]),linestyle='-',c='k',linewidth=0.5,alpha=1)
    ax1.plot((xdata[i]-0.5*dr)/10000,ydata[i],marker='o',linestyle='none',c='k',markersize=0.75)
    ax1.errorbar((xdata[i]-0.5*dr)/10000,ydata[i],k1_std[i],c='k',linewidth=0.5,capsize=0.8)
    ax1.plot(((xdata[i]-dr)/10000, xdata[i]/10000),(k2_mean[i],k2_mean[i]),linestyle='-',c='r',linewidth=0.5,alpha=1)
    ax1.plot((xdata[i]-0.5*dr)/10000,k2_mean[i],marker='o',linestyle='none',c='r',markersize=0.75)
    ax1.errorbar((xdata[i]-0.5*dr)/10000,k2_mean[i],k2_std[i],c='r',linewidth=0.5,capsize=0.8)
    ax1.plot(((xdata[i]-dr)/10000, xdata[i]/10000),(k3_mean[i],k3_mean[i]),linestyle='-',c='g',linewidth=0.5,alpha=1)
    ax1.plot((xdata[i]-0.5*dr)/10000,k3_mean[i],marker='o',linestyle='none',c='g',markersize=0.75)
    ax1.errorbar((xdata[i]-0.5*dr)/10000,k3_mean[i],k3_std[i],c='g',linewidth=0.5,capsize=0.8)
    #ax1.errorbar((xdata[i]-0.5*dr)/10000,ydata[i], yerr = yerr[i],marker='.',linestyle='none',markersize=0.75, fmt='-',capsize=0.75, capthick=0.25, elinewidth=0.5,linewidth=0.5,c='K',alpha=0.6)
    if i < len(k1_mean)-1:
        ax1.plot((xdata[i]/10000, xdata[i]/10000),(ydata[i],ydata[i+1]),linestyle='-',c='k',linewidth=0.5,alpha=1)
        ax1.plot((xdata[i]/10000, xdata[i]/10000),(k2_mean[i],k2_mean[i+1]),linestyle='-',c='r',linewidth=0.5,alpha=1)
        ax1.plot((xdata[i]/10000, xdata[i]/10000),(k3_mean[i],k3_mean[i+1]),linestyle='-',c='g',linewidth=0.5,alpha=1)

# new = [1.537e-20, -7.306e-16, 1.283e-11, -1.004e-07, 0.0004144, -0.4888]
### Fitting the output rdf

#----------simulation parameters for 1-1 ----------
R = 65000 # PBC radius, A
R_NUCLEUS = 43000 # NE radius, A
N_GRANULES = 306 # Number of granules
R_GRANULES = 1200 # Radius of granules, A
R_PATCHES = 100 # Radius of the binding patch,
R_GLUCOSE = 100 # A, Radius of glucose, A
ISOS_CONTACT_RANGE = R_GLUCOSE*3 # contact range of granule surface to PM under which a collision is counted, A

# xdata = np.linspace(R_GRANULES + R_PATCHES, R - R_NUCLEUS - (R_GRANULES + R_PATCHES + ISOS_CONTACT_RANGE), 100)
# ydata2_out_fit = poly5_func(xdata, new)
#ax1.plot(xdata/10000 -0.1,ydata2_out_fit,linestyle='-', color='red',linewidth=0.5,alpha=1) 
        
#ax1.plot((0.9, 1),(2.42,2.42),'k',linestyle="-",linewidth = 0.5)
#ax1.plot((0.9, 1),(2.21,2.21),'red',linestyle="-",linewidth = 0.5)
#ax1.plot((0.06, 0.16),(2.0,2.0),'green',linestyle="-",linewidth = 0.5)

# build a rectangle in axes coords
left, width = .25, .5
bottom, height = .25, .5
right = left + width
top = bottom + height
#ax1.text(0.02*(left+right), 0.96*(bottom+top), '0$\,$mM Glu',horizontalalignment='left',verticalalignment='center',fontsize=8,fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.52*(left+right), 0.96*(bottom+top), 'input RDF',horizontalalignment='left',verticalalignment='center',fontsize=8,fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.52*(left+right), 0.87*(bottom+top), 'output RDF',horizontalalignment='left',verticalalignment='center',fontsize=8,fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.1*(left+right), 0.78*(bottom+top), '25$\,$mM Glu + Ex-4',horizontalalignment='left',verticalalignment='center',fontsize=8,fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.02*(left+right), 0.89*(bottom+top), 'anodic layer',horizontalalignment='left',verticalalignment='center',fontsize=20,fontname="Arial",fontweight="normal",color='chocolate',transform=ax1.transAxes)
#ax1.text(0.02*(left+right), 0.83*(bottom+top), 'bilayer',horizontalalignment='left',verticalalignment='center',fontsize=20,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)

plt.show()

# Save figures, (PNG,JPG,EPS,SVG,PGF and PDF supported, PDF is recommanded or PGF)
# fig.savefig('../../05_rdf/results/'+ str(seqs[0]) + '_rdf_ne_pm.png', transparent=True, dpi=300)
fig.savefig('../../05_rdf/results/' + '_total_rdf_ne_pm.png', transparent=True, dpi=300)

#print("--- %s seconds ---" % (time.time() - start_time))

# %%
