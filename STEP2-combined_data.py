
import numpy as np
import sys, math
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 

### READING INPUT ###

file1 = input("Enter colvars: ")
file2 = input("Enter SI: ")
ofile1 = "system-0-combined.dat"

### Making INPUT files ###
df1 = pd.read_csv(file1,header=0,delim_whitespace=True,comment='#',names=['frame','sl','k2','ac'])
df2 = pd.read_csv(file2,header=0,delim_whitespace=True,comment='#',names=['frame','si'])

data = {'f':df2['frame'],'si':df2['si'],'sl':df1['sl'],'ac':df1['ac'],'k2':df1['k2']}
df = pd.DataFrame.from_dict(data=data)

print(df)
df.to_csv(ofile1,sep=' ',index=False)
