
import numpy as np
import sys, math
from timeit import default_timer as timer
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from matplotlib.font_manager import FontProperties
#import dask.dataframe as dd
from matplotlib import rc
from tqdm import tqdm
from scipy.ndimage.filters import gaussian_filter

#### FUNCTION ####

def cal_2d(x,y,temp,R,pmf_max):
    H, xedges, yedges = np.histogram2d(x,y,density=True,bins=40,range=[[0,50],[0,y.max()]])
    stepx = xedges[1]-xedges[0]
    stepy = yedges[1]-yedges[0]
    xx, yy = np.mgrid[0:50:stepx,yedges.min():yedges.max():stepy]
    pos = np.dstack((xx, yy))
    pmax = 0
    for i in H:
        p = i.sum()
        if p >=pmax:
            pmax = p
    print("Found pmax = ",pmax)

    for i in range(len(H)):
        for j in range(len(H.T)):
            if H[i,j]!=0:
                H[i,j]=-R*temp*np.log(H[i,j]/pmax)
            else:
                H[i,j]=pmf_max
    return pos,H

### READING INPUT ###

file1 = input("Enter colvars: ")
ofile1 = input("Enter output: ")

### Making INPUT files ###
df = pd.read_csv(file1,header=0,delim_whitespace=True)

### INPUT ###
TEMP = 300
R = 0.002 # kcal/mol
PMF_max = 6 # kcal / mol

X = df['sl']
y_list = ['si','ac','k2']

print("Finished reading input")

fig,ax_list = plt.subplots(1,3, sharex=False, sharey=False, squeeze=True,tight_layout=True)
ax_list = ax_list.flatten()

y_labels = ['Stacking Index, Xi_pi','Acylindricity, c','Relative Anisotropy, K2']

for a in range(ax_list.size):
    ofile = open("system-0-sl-"+y_list[a]+".dat", "w")
    ofile.write("sl  "+y_list[a]+"   pmf"+'\n')
    
    pos,H = cal_2d(X,df["%s"%(y_list[a])],TEMP,R,PMF_max)
    ofile.write("%s\t%s\t%s\n" % (pos[0],pos[1],H))
    
    H = gaussian_filter(H,0.6)
    
    plot = ax_list[a].contourf(pos[:,:,0],pos[:,:,1],H,cmap="viridis_r",extend='max')
    cbar = fig.colorbar(plot,ax=ax_list[a],shrink=0.9,spacing='proportional')
    cbar.set_label('kcal/mol',rotation=270,labelpad=15)
    cbar.minorticks_on()
    ax_list[a].set_xlabel('SL')
    ax_list[a].set_ylabel("%s"%(y_labels[a]))
    ofile.close()

# ### MISCELLANEOUS ###
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(14,4)
plt.tight_layout(pad=0.05)

# ### SAVE FIGURE ### 
# plt.savefig(ofile2,dpi=300,orientation='lanscape', papertype=None, format=None,
#          transparent=False, bbox_inches=None, pad_inches=None,metadata=None)

plt.show()
