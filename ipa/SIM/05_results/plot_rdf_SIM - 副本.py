#!/usr/bin/python

# Import module for plot "matplotlib"
#%%
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
import scipy

def poly5_func(xdata, pars):
    '''
    calculate a fifth order polynomial equation with xdata and a list of parameter 
    '''
    return pars[0]*xdata**5 +pars[1]*xdata**4 + pars[2]*xdata**3 + pars[3]*xdata**2 + pars[4]*xdata +pars[5]

file_dir = 'I:\\Bing\\fluorescence\\3D\\04_rdf\\' 
N_shell = 10

file_list=[]
with open(f'../file_list_MT_ex4.txt','r') as f:                           # Modified
	for line in f:
		file_list.append(line.strip('\n'))

# Read data file
basal_fp = []
for file_name in file_list[:10]:   # 2.8 mM
    print('b:'+file_name)
    k1 = np.loadtxt(f'{file_dir}/results/{file_name}_bin2_isgrdf_ne-pm.xvg')
    basal_fp.append(k1)

glu_5_fp = []
for file_name in file_list[10:19]:   # 16.7 mM - 5 min
    print('5:'+file_name)
    k1 = np.loadtxt(f'{file_dir}/results/{file_name}_bin2_isgrdf_ne-pm.xvg')
    glu_5_fp.append(k1)

glu_30_fp = []
for file_name in file_list[19:]:   # 16.7 mM - 30 min
    print('30:'+file_name)
    k1 = np.loadtxt(f'{file_dir}/results/{file_name}_bin2_isgrdf_ne-pm.xvg')
    glu_30_fp.append(k1)   

basal_mean = np.mean(np.array(basal_fp), axis=0)
basal_sem = np.std(np.array(basal_fp), axis=0) / np.sqrt(len(basal_fp)-1)
glu_5_mean = np.mean(np.array(glu_5_fp), axis=0)
glu_5_sem = np.std(np.array(glu_5_fp), axis=0) / np.sqrt(len(glu_5_fp)-1)
glu_30_mean = np.mean(np.array(glu_30_fp), axis=0)
glu_30_sem = np.std(np.array(glu_30_fp), axis=0) / np.sqrt(len(glu_30_fp)-1)

Rcell = 53946
Rne = 34939
Risg = 1289
dr = (Rcell - Rne)/9

xdata = np.linspace(1, 9, len(k1))*dr


#%%
'''
# Plots
'''

# Specifiy environmental parameter
rc('font',**{'family':'serif','serif':['Arial Narrow']})

# Create axes 
fig = plt.figure(figsize=(8,6)) #inches
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
ax1.set_ylim([0,0.3])
xmajorLocator   = MultipleLocator(0.4)
xmajorFormatter = FormatStrFormatter('%.1f')
xminorLocator   = MultipleLocator(0.2)
ax1.xaxis.set_major_locator(xmajorLocator)
ax1.xaxis.set_major_formatter(xmajorFormatter)
ax1.xaxis.set_minor_locator(xminorLocator)
ymajorLocator   = MultipleLocator(0.1)
ymajorFormatter = FormatStrFormatter('%.1f')
yminorLocator   = MultipleLocator(0.25)
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

# ax1.plot((0,100),(1,1),'k',linestyle="--",linewidth = 0.5)

#yerr= k1[:,1]
for i in range(0,len(basal_mean)):
    line1, = ax1.plot(((xdata[i]-dr)/10000, xdata[i]/10000),(basal_mean[i],basal_mean[i]),linestyle='-',c='k',linewidth=0.5,alpha=1)
    ax1.plot((xdata[i]-0.5*dr)/10000,basal_mean[i],marker='o',linestyle='none',c='k',markersize=0.75)
    ax1.errorbar((xdata[i]-0.5*dr)/10000,basal_mean[i], yerr = basal_sem[i],marker='.',linestyle='none',markersize=0.75, fmt='-',capsize=0.75, capthick=0.25, elinewidth=0.5,linewidth=0.5,c='k',alpha=0.6)

    line2, = ax1.plot(((xdata[i]-dr)/10000, xdata[i]/10000),(glu_5_mean[i],glu_5_mean[i]),linestyle='-',c='r',linewidth=0.5,alpha=1)
    ax1.plot((xdata[i]-0.5*dr)/10000,glu_5_mean[i],marker='o',linestyle='none',c='r',markersize=0.75)
    ax1.errorbar((xdata[i]-0.5*dr)/10000,glu_5_mean[i], yerr = glu_5_sem[i],marker='.',linestyle='none',markersize=0.75, fmt='-',capsize=0.75, capthick=0.25, elinewidth=0.5,linewidth=0.5,c='r',alpha=0.6)

    line3, = ax1.plot(((xdata[i]-dr)/10000, xdata[i]/10000),(glu_30_mean[i],glu_30_mean[i]),linestyle='-',c='g',linewidth=0.5,alpha=1)
    ax1.plot((xdata[i]-0.5*dr)/10000,glu_30_mean[i],marker='o',linestyle='none',c='g',markersize=0.75)
    ax1.errorbar((xdata[i]-0.5*dr)/10000,glu_30_mean[i], yerr = glu_30_sem[i],marker='.',linestyle='none',markersize=0.75, fmt='-',capsize=0.75, capthick=0.25, elinewidth=0.5,linewidth=0.5,c='g',alpha=0.6)
    if i < len(basal_mean)-1:
        ax1.plot((xdata[i]/10000, xdata[i]/10000),(basal_mean[i],basal_mean[i+1]),linestyle='-',c='k',linewidth=0.5,alpha=1)
        ax1.plot((xdata[i]/10000, xdata[i]/10000),(glu_5_mean[i],glu_5_mean[i+1]),linestyle='-',c='r',linewidth=0.5,alpha=1)
        ax1.plot((xdata[i]/10000, xdata[i]/10000),(glu_30_mean[i],glu_30_mean[i+1]),linestyle='-',c='g',linewidth=0.5,alpha=1)


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
# plt.legend([p2, p1], ["line 2", "line 1"], loc='upper left')
# ax1.legend([line1, line2, line3], ['2,8mM-30min', '16,7mM-5min', '16,7mM-30min'])
plt.show()

# Save figures, (PNG,JPG,EPS,SVG,PGF and PDF supported, PDF is recommanded or PGF)
fig.savefig('_rdf_ne_pm.png', transparent=True, dpi=300)

#print("--- %s seconds ---" % (time.time() - start_time))

# %%
