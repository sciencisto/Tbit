#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 14:08:37 2020

@author: kristine
"""

"""
./wraparound.py > run_calc.sh
"""
    
filepath = input("Enter path to the kristine_thesis_runthrough directory") #'/home/kristine/Documents/Thesis/github_vers4/kristine_thesis_runthrough'
filename = input("Enter filename")

Binary_method_list = ["exact", "group"]

print("echo '########################################################'")      
print("echo 'Starting the TCRdistance program, Kristines Thesis 2020'")
print("echo '########################################################'")      


import pandas as pd
import math
TSV = pd.read_csv(filepath+"/"+filename, sep="\t", header=None)
length = math.ceil(len(TSV)/50)
remain = len(TSV)-(length*50)
del TSV
      
print("mkdir TCR_validate")
for i in range(length):
    print("python3 validate_fold1.py '%s/' '%s' 50 %s" % (filepath, filename, i))
    print("mv *dmat.tsv TCR_validate")

print("python3 eval_distribution_plot.py 'TCR_validate' 'TCR_validate'")
      
print("mv *.png distribution_plots")

for i in range(length):
    print("python3 subsample_validate.py '%s' '%s' '%s/Binary_group' %s 50 %s" % (filepath, filename, filepath, i, remain))


#print("echo '########################################################'")      
#print("echo 'Subsampling'")
#print("echo '########################################################'")    
#print("python3 subsample_all.py '%s' '%s' '%s/Binary_group' %s 50 %s" % (filepath, filename, filepath, length, remain))

print("chmod 777 *")