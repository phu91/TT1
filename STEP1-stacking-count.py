import MDAnalysis as md
import matplotlib.pyplot as plt 
import numpy as np 
from MDAnalysis.analysis import distances
import pandas as pd
import scipy.stats
import seaborn as sns 
import itertools
import time
from tqdm import tqdm
from numpy import linalg

def count_si(x):
    s = 0
    if x[1]<=3.5:
        if x[2]<=30:
            s = 1.0
        if 30< x[2] <=60:
            s = 0.5
        if 60<x[2]<=90:
            s = 0.25
        if 90<x[2]<=120:
            s = 0.25
        if 120<x[2]<=150:
            s = 0.5
        if 150<x[2]<=180:
            s = 1.0
    else:
        s = 0
    return s

u = md.Universe('dry.tt1_12mol_cubic.prmtop','tt1-12-test2-100ns.dcd')

csvfile = "system-B-si.dat"

distance = []
angle = []

for ts in tqdm(u.trajectory[::1]):
    C6  = u.select_atoms('name C6')
    O19 = u.select_atoms('name O19')
    vec1 = C6.positions-O19.positions
    C31 = u.select_atoms('name C31') 
    O7  = u.select_atoms('name O7')
    vec2 = C31.positions-O7.positions
    for x1,y1,z1 in vec1:
        for x2,y2,z2 in vec1:
            vec1_mag = np.sqrt(x1**2+y1**2+z1**2)
            vec2_mag = np.sqrt(x2**2+y2**2+z2**2)
            d11 = np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
            dot = x1*x2+y1*y2+z1*z2
            theta11 = float(np.arccos(dot/(vec1_mag*vec2_mag))*180/3.14)
            angle.append((int(ts.frame)+1,d11,theta11))
    for x1,y1,z1 in vec1:
        for x2,y2,z2 in vec2:
            vec1_mag = np.sqrt(x1**2+y1**2+z1**2)
            vec2_mag = np.sqrt(x2**2+y2**2+z2**2)
            d12 = np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
            dot = x1*x2+y1*y2+z1*z2
            theta12 = float(np.arccos(dot/(vec1_mag*vec2_mag))*180/3.14)
            angle.append((int(ts.frame)+1,d12,theta12))
    for x1,y1,z1 in vec2:
        for x2,y2,z2 in vec2:
            vec1_mag = np.sqrt(x1**2+y1**2+z1**2)
            vec2_mag = np.sqrt(x2**2+y2**2+z2**2)
            d22 = np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
            dot = x1*x2+y1*y2+z1*z2
            theta22 = float(np.arccos(dot/(vec1_mag*vec2_mag))*180/3.14)
            angle.append((int(ts.frame)+1,d22,theta22))
#    time.sleep(10)

df = pd.DataFrame(angle,columns=['f','d','a'])
df = df[df['d']!=0]
df = df.dropna().drop_duplicates(keep='first',ignore_index=True)
df['si'] = df.apply(count_si,axis=1)
df = df.groupby(['f']).sum()
df = df.drop(['d','a'],axis=1)
df.reset_index().to_csv(csvfile,index=False,index_label=False,sep='\t')
