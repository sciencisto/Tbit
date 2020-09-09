#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 17:53:15 2020

@author: kristine
"""

import glob as glob 
import pandas as pd
import numpy as np
"""
path = '/home/kristine/Desktop/new_VDJ/Test/Distances'
all_files = sorted(glob.glob(path + "/*.tsv"))
print(all_files)
li = []

for filename in all_files:
    df = pd.read_csv(filename, sep = '\t', index_col=None, header=None)
    li.append(df)

frame = pd.concat(li, axis=0, ignore_index=True)
t_dmat = frame.to_numpy()

path = '/home/kristine/Desktop/new_VDJ/Test/binary'
all_files = sorted(glob.glob(path + "/*.tsv"))
print(all_files)
li = []

for filename in all_files:
    df = pd.read_csv(filename, sep = '\t', index_col=None, header=None)
    li.append(df)

frame = pd.concat(li, axis=0, ignore_index=True)
b_dmat = frame.to_numpy()


t = pd.read_csv('/home/kristine/Desktop/new_VDJ/Test/test_set.tsv', sep="\t", header = None)
t.columns = ["CDR3", "V", "Ep"]
t = list(t.Ep)

t_Epi = []
for i in range(len(t)):
    E = t[i].split(":")[0]
    t_Epi.append(E)
    
"""

# flatten each, and make three columns, with 830 in each 

#flatten
#transpose

from sklearn.utils import resample
import pandas as pd, seaborn as sns
sns.set_style("white")

#bin_ = np.array(b_dmat).flatten()  
#dist = np.array(t_dmat).flatten()  

#y = np.repeat(t_Epi,830)

#d = {"a": dist, "b": y, "c": bin_}
#df = pd.DataFrame(d)

#df1 = df[df.c == 1]
#df1 = df1.drop(c)

g = sns.pairplot(df1, hue="b", palette="dark")




































