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

CDR12_method_list = ["TM", "RMSD", "Blosum62", "hydrophobicity"]
CDR3_method_list = ["Hamming", "Blosum62", "Blosum45", "Pam10", "K-Mer_3", "K-Mer_4", "K-Mer_5", "Atchley_Factors", "hydrophobicity", "Immunomap"]
Binary_method_list = ["exact", "group"]

print("echo '########################################################'")      
print("echo 'Starting the TCRdistance program, Kristines Thesis 2020'")
print("echo '########################################################'")      

      
##########################################################################
#writing out the chunks, i wish to have 
      
import pandas as pd
import math
TSV = pd.read_csv(filepath+"/"+filename, sep="\t", header=None)
length = math.ceil(len(TSV)/50)
remain = len(TSV)-(length*50)
del TSV
      
for j in range(len(CDR12_method_list)):
    print("echo 'Calculating dmats for the %s method for CDR1'" % (CDR12_method_list[j]))
    print("mkdir CDR1dist_%s" % (CDR12_method_list[j]))
    for i in range(length):
        print("python3 TCR_distance.py '%s/' '%s' 50 %s 'none' '%s' 'none'" % (filepath, filename, i, CDR12_method_list[j]))
    print("mv *dmat.tsv CDR1dist_%s" % (CDR12_method_list[j]))
for j in range(len(CDR12_method_list)):
    print("echo 'Calculating dmats for the %s method for CDR2'" % (CDR12_method_list[j]))
    print("mkdir CDR2dist_%s" % (CDR12_method_list[j]))
    for i in range(length):
        print("python3 TCR_distance.py '%s/' '%s' 50 %s 'none' 'none' '%s'" % (filepath, filename, i, CDR12_method_list[j]))
    print("mv *dmat.tsv CDR2dist_%s" % (CDR12_method_list[j]))
for j in range(len(CDR3_method_list)):
    print("echo 'Calculating dmats for the %s method for CDR3'" % (CDR3_method_list[j]))
    print("mkdir CDR3dist_%s" % (CDR3_method_list[j]))
    for i in range(length):
        print("python3 TCR_distance.py '%s/' '%s' 50 %s '%s' 'none' 'none'" % (filepath, filename, i, CDR3_method_list[j]))
    print("mv *dmat.tsv CDR3dist_%s" % (CDR3_method_list[j]))
    
print("echo '########################################################'")      
print("echo 'all dmats have been calculated, now the binary files are generated'")
print("echo '########################################################'")      

for j in range(len(Binary_method_list)):
    print("echo 'Calculating binary dmat for the %s match'" % (Binary_method_list[j]))
    print("mkdir Binary_%s" % (Binary_method_list[j]))
    for i in range(length):
        print("python3 TCR_distance_binary.py '%s/' '%s' 50 %s '%s'" % (filepath, filename, i, Binary_method_list[j]))
    print("mv *dmat.tsv Binary_%s" % (Binary_method_list[j]))

print("echo '########################################################'")      
print("echo 'Making distribution plots for all'")
print("echo '########################################################'")      
print("mkdir distribution_plots")
for j in range(len(CDR12_method_list)):
    print("python3 distribution_plot.py 'CDR1dist_%s' 'CDR1_%s'" % (CDR12_method_list[j], CDR12_method_list[j]))
    print("python3 distribution_plot.py 'CDR2dist_%s' 'CDR2_%s'" % (CDR12_method_list[j], CDR12_method_list[j]))
for j in range(len(CDR3_method_list)):
    print("python3 distribution_plot.py 'CDR3dist_%s' 'CDR3_%s'" % (CDR3_method_list[j], CDR3_method_list[j]))
for j in range(len(Binary_method_list)):
    print("python3 distribution_plot.py 'Binary_%s' 'binary_%s'" % (Binary_method_list[j], Binary_method_list[j]))
print("mv *.png distribution_plots")

print("echo '########################################################'")      
print("echo 'Subsampling'")
print("echo '########################################################'")    
print("python3 generate_sub_cdr1.py '%s' '%s' '%s/Binary_group' %s 50 %s" % (filepath, filename, filepath, length, remain))
print("python3 generate_sub_cdr2.py '%s' '%s' '%s/Binary_group' %s 50 %s" % (filepath, filename, filepath, length, remain))
print("python3 generate_sub_cdr3.py '%s' '%s' '%s/Binary_group' %s 50 %s" % (filepath, filename, filepath, length, remain))

print("mkdir distance_matrices")
print("chmod 777 *")
print("mv CDR1dist_* distance_matrices")
print("mv CDR2dist_* distance_matrices")
print("mv CDR3dist_* distance_matrices")
