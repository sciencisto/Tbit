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
      
print("mkdir CDR1dist_volume")
for i in range(length):
    print("python3 calc_TCR_distance.py '%s/' '%s' 50 %s 'none' 'volume' 'none'" % (filepath, filename, i))
    print("mv *dmat.tsv CDR1dist_volume")

print("mkdir CDR2dist_volume")
for i in range(length):
    print("python3 calc_TCR_distance.py '%s/' '%s' 50 %s 'none' 'none' 'volume'" % (filepath, filename, i))
    print("mv *dmat.tsv CDR2dist_volume")
    
print("mkdir CDR3dist_volume")
for i in range(length):
    print("python3 calc_TCR_distance.py '%s/' '%s' 50 %s 'volume' 'none' 'none'" % (filepath, filename, i))
    print("mv *dmat.tsv CDR3dist_volume")
    
print("echo '########################################################'")      
print("echo 'Making distribution plots for all'")
print("echo '########################################################'")      

print("python3 eval_distribution_plot.py 'CDR1dist_volume' 'CDR1_volume'")
print("python3 eval_distribution_plot.py 'CDR2dist_volume' 'CDR2_volume'")
print("python3 eval_distribution_plot.py 'CDR3dist_volume' 'CDR3_volume'")     
      
print("mv *.png distribution_plots")

print("echo '########################################################'")      
print("echo 'Subsampling'")
print("echo '########################################################'")    
print("python3 subsample_all.py '%s' '%s' '%s/Binary_group' %s 50 %s" % (filepath, filename, filepath, length, remain))

print("chmod 777 *")