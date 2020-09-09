#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 08:56:32 2020

@author: kristine
"""

filepath = input("Enter path directory")

filename = input("Enter filename")

print("echo '########################################################'")      
print("echo 'Starting the TBit'")
print("echo '########################################################'")      

      
import pandas as pd
import math

#Import the file and calculate the number of chunks needed
TSV = pd.read_csv(filepath+"/"+filename, sep="\t", header=None)
length = math.ceil(len(TSV)/50)
#count how many TCRs are present in the last chunk
remain = len(TSV)-((length-1)*50)
del TSV
#Impart the file and calculate the number of chunks needed
TSV = pd.read_csv(filepath+"/"+filename, sep="\t", header=None)
length = math.ceil(len(TSV)/50)
#count how many TCRs are present in the last chunk
remain = len(TSV)-((length-1)*50)
del TSV
      
print("mkdir Tbit_dist")
for i in range(length):
    print("python3 calc_TCR_distance_model.py '%s/' '%s' 50 %s" % (filepath, 
                                                                   filename, i))
    print("mv *dmat.tsv Tbit_dist")


##############################################################################
# This part is only if evaluation of the distance metric is desired
##############################################################################

#Distribution of the distance
print("python3 eval_distribution_plot.py 'BetaDistance' 'BetaDistance'")

#Generating a Binary file, notice the option "group" require a library file. 
#Alternatively choose option "exact"
print("mkdir Binary")
for i in range(length):
    print("python3 eval_TCR_distance_binary.py '%s/' '%s' 50 %s 'group'" 
          % (filepath, filename, i))
    print("mv *dmat.tsv Binary")

#Subsampling the file outouts one csv file with two columns, distance and 
#binary information here i get some problems
for i in range(length):
    print("python3 subsample_validate.py '%s/BetaDistance' '%s' '%s/Binary' %s 50 %s" % 
          (filepath, filename, filepath, i, remain))

print("python3 ROC.py '%s/BetaDist_subsampled.csv'" % (filepath))


print("chmod 777 *")