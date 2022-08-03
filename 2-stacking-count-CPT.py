import MDAnalysis as md
import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd
import seaborn as sns 
from tqdm import tqdm
import math
#import sys 


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

u = md.Universe('dry.cpt_24.prmtop','cpt_24_100.dcd')

ftrfile = "system-0.ftr"
csvfile = "system-0.dat"

raw_data = []

for ts in tqdm(u.trajectory[::1]):
    C16  = u.select_atoms('name C16')
    O2 = u.select_atoms('name O2')
    vec1 = C16.positions-O2.positions

    for x1,y1,z1 in vec1:
        vec1_mag = np.sqrt(x1**2+y1**2+z1**2)
        for x2,y2,z2 in vec1:
            vec2_mag = np.sqrt(x2**2+y2**2+z2**2)
            d11 = (np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2))
            dot = x1*x2+y1*y2+z1*z2
            theta11 = (np.arccos(dot/(vec1_mag*vec2_mag))*180/3.14).round(1)
            raw_data.append((ts.frame+1,round(d11,1),round(theta11,1)))
            dot = 0
            d11 = 0

df = pd.DataFrame(raw_data,columns=['f','d','a'])

df = df[df['d']!=0]
df = df.dropna().drop_duplicates(keep='first',ignore_index=True)
df['si'] = df.apply(count_si,axis=1)

df = df.groupby(['f']).sum()
# print(df)
df = df.drop(['d','a'],axis=1)
# df.reset_index().to_feather(ftrfile)
df.reset_index().to_csv(csvfile,index=False,index_label=False,sep='\t')
