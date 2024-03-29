#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import sys, math
from timeit import default_timer as timer
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from matplotlib.font_manager import FontProperties
import MDAnalysis as md
from matplotlib import rc
from tqdm import tqdm


def find_peak(data,sl_min,sl_max,rc,rc_min,rc_max):
    sl_max = data['sl']<=sl_max
    sl_min = data['sl']>=sl_min
    if rc =="si":
        si_max = data['si']<=rc_max
        si_min = data['si']>=rc_min
        frame = data.loc[sl_min & sl_max & si_min & si_max,['f']]-1
    elif rc =="ac":
        ac_max = data['ac']<=rc_max
        ac_min = data['ac']>=rc_min
        frame = data.loc[sl_min & sl_max & ac_min & ac_max,['f']]-1
    elif rc =="k2":
        k2_max = df['k2']<=rc_max
        k2_min = df['k2']>=rc_min
        frame = data.loc[sl_min & sl_max & k2_min & k2_max,['f']]-1
    frames = frame.values
    nf = +len(frames)
    # print("There are %s to write." % nf)
    return frames

### I/O ###

### INPUT EXAMPLES ###

# These steps can be skipped if all parameters are concatnated. 
# df1 = pd.read_csv("si-example.dat",header=0,delim_whitespace=True)
# df2 = pd.read_csv("si-ac-k2-example.dat",header=0,delim_whitespace=True)
# df = pd.concat([df1,df2],axis=1)

df = pd.read_csv("system-A-combined.dat",header=0,delim_whitespace=True)

#This choose S_L (PARAM1) as the reference range, and then select feature (SI, AC, K2) or PARAM2 with its min-max range 
# Need to be rewitten to select correct parameters 

peak = str(input("Peak name: "))
min_sl = float(input("min S_L: "))
max_sl = float(input("max S_L: "))
feature = str(input("Choose si,ac,k2: "))
param_min = float(input("min value: "))
param_max = float(input("max value: "))

frames = find_peak(df,min_sl,max_sl,feature,param_min,param_max)
nf = +len(frames)

print("There are %s to write." % nf)
 
u = md.Universe('PRMTOP','DCD')
selected_Atoms = u.select_atoms('all')
with md.Writer("OUTPUT.DCD",selected_Atoms.n_atoms) as W:
    for ts in u.trajectory[::1]:
        current = ts.frame
        for i in frames:
            if current == i:
                W.write(selected_Atoms)
