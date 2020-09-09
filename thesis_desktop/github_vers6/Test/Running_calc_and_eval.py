#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 14:08:37 2020

@author: kristine
"""

"""
./Running_calc_and_eval.py > calculate_and_evaluate.sh
"""
    
filepath = input("Enter path directory")
filename = input("Enter filename")

CDR12_method_list = ["TM", "RMSD", "Blosum62", "hydrophobicity"]
CDR3_method_list = ["Hamming", "Blosum62", "Blosum45", "Pam10", "K-Mer_3", "K-Mer_4", "K-Mer_5", "Atchley_Factors", "hydrophobicity", "Immunomap"]
Binary_method_list = ["group"]

print("echo '########################################################'")      
print("echo 'Starting the TCRdistance program, Kristines Thesis 2020'")
print("echo '########################################################'")      

      
##########################################################################
#firstly i count the number of chunks i want.    
import pandas as pd
import math
TSV = pd.read_csv(filepath+"/"+filename, sep="\t", header=None)
length = math.ceil(len(TSV)/50)
rest = len(TSV)-(length-1)*50
del TSV

print("mkdir all_methods")
for i in range(length):
    print("python3 calc_TCR_distance_all_methods.py '%s/' '%s' 50 %s" % (filepath, filename, i))
    print("mv *dmat.tsv all_methods")
#calculate the dmats for all methdods and moving them into their designated directories      
print("mkdir Immunomap")
for i in range(length):
    print("python3 calc_TCR_distance.py '%s/' '%s' 50 %s 'Immunomap' 'none' 'none'" % (filepath, filename, i))
    print("mv *dmat.tsv Immunomap")
    
print("echo '########################################################'")      
print("echo 'all dmats have been calculated, now the binary files are generated'")
print("echo '########################################################'")      

for j in range(len(Binary_method_list)):
    print("echo 'Calculating binary dmat for the %s match'" % (Binary_method_list[j]))
    print("mkdir Binary_%s" % (Binary_method_list[j]))
    for i in range(length):
        print("python3 eval_TCR_distance_binary.py '%s/' '%s' 20 %s '%s'" % (filepath, filename, i, Binary_method_list[j]))
    print("mv *dmat.tsv Binary_%s" % (Binary_method_list[j]))

print("echo '########################################################'")      
print("echo 'Generation of confusion matrices and TPR & FPR'")
print("echo '########################################################'")     
print("mkdir ROC_group_matches")

#a little piece of hard coding for shittyness      
all_folders = ['all_methods','Immunomap']

for i in range(len(all_folders)):
    print("python3 eval_confus_matrix.py '%s/%s/' '%s/Binary_group' 'ROC_group_matches'"% (filepath, all_folders[i], filepath))            

print("echo '########################################################'")      
print("echo 'Generation of TPR and ROC curves'")
print("echo '########################################################'")     

print("python3 eval_ROC.py '%s/ROC_group_matches'" % (filepath))
print("mv *.png ROC_group_matches")

print("echo '########################################################'")      
print("echo 'Starting logarithmic regression calculations to find coeficients'")
print("echo '########################################################'")     
      
      
#print("echo 'AFFINITY PROPERGATION'")
#print("echo '########################################################'")  
      
#print("mkdir Affinity_propergation")
#for j in range(len(CDR12_method_list)):
#    print("python3 Affinity_propergation.py 'CDR1dist_%s' 'CDR1_%s' '%s' '%s'" % (CDR12_method_list[j], CDR12_method_list[j], filepath, filename))
#    print("python3 Affinity_propergation.py 'CDR2dist_%s' 'CDR2_%s' '%s' '%s'" % (CDR12_method_list[j], CDR12_method_list[j], filepath, filename))
#for j in range(len(CDR3_method_list)):
#    print("python3 Affinity_propergation.py 'CDR3dist_%s' 'CDR3_%s' '%s' '%s'" % (CDR3_method_list[j], CDR3_method_list[j], filepath, filename))
#print("mv *AP.png Affinity_propergation")
#print("mv *AP_stats.txt Affinity_propergation")

#print("echo 'Networks'")
#print("echo '########################################################'")  
#print("mkdir Networks")
#print("cd Networks")
#print("mkdir _0.6")
#print("mkdir _0.8")
#print("chmod 777 *")
#print("cd ..")
#for j in range(len(CDR12_method_list)):
#    print("python3 network.py 'CDR1dist_%s' 'CDR1_%s' '%s' '%s' 0.6" % (CDR12_method_list[j], CDR12_method_list[j], filepath, filename))
#    print("python3 network.py 'CDR2dist_%s' 'CDR2_%s' '%s' '%s' 0.6" % (CDR12_method_list[j], CDR12_method_list[j], filepath, filename))
#for j in range(len(CDR3_method_list)):
#    print("python3 network.py 'CDR3dist_%s' 'CDR3_%s' '%s' '%s' 0.6" % (CDR3_method_list[j], CDR3_method_list[j], filepath, filename))
#print("mv *network.png Networks/_0.6")
#for j in range(len(CDR12_method_list)):
#    print("python3 network.py 'CDR1dist_%s' 'CDR1_%s' '%s' '%s' 0.8" % (CDR12_method_list[j], CDR12_method_list[j], filepath, filename))
#    print("python3 network.py 'CDR2dist_%s' 'CDR2_%s' '%s' '%s' 0.8" % (CDR12_method_list[j], CDR12_method_list[j], filepath, filename))
#for j in range(len(CDR3_method_list)):
#    print("python3 network.py 'CDR3dist_%s' 'CDR3_%s' '%s' '%s' 0.8" % (CDR3_method_list[j], CDR3_method_list[j], filepath, filename))
#print("mv *network.png Networks/_0.8")

#print("mkdir distance_matrices")
#print("chmod 777 *")
#print("mv CDR1dist_* distance_matrices")
#print("mv CDR2dist_* distance_matrices")
#print("mv CDR3dist_* distance_matrices")