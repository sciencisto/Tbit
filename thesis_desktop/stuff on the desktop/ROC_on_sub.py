#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 13:17:24 2020

@author: kristine
"""
import pandas as pd
import glob as glob
import numpy as np
from sklearn import metrics
from sklearn.metrics import roc_auc_score
from sys import argv 
import time 

#original ROC
def parse_commandline():  
    """
    Rembember to login to databar by ssh -X.
	commandline call: python3 distribution_plot.py "concat_test/blos62/" "CDR1_blos" 
	python3 distribution_plot.py "/work1/s111518/new_VDJ/CDR2dist_TM/" ""/work1/s111518/new_VDJ/Binary_group/"
	"""
    filepath = argv[1] #path to distance matrix

    return filepath

def import_files(path):   
    d = {"Hamming": [1], "Blosum62": [1], "Blosum45": [1], "Pam10": [1], "K-Mer_3": [1], "K-Mer_4": [1], 
         "K-Mer_5": [1], "Atchley_Factors": [1], "hydrophobicity": [1], "Bin": [0]}
    df = pd.DataFrame(data=d)    
    file = pd.read_csv(path, index_col=None, header = None)
    file.columns = df.columns
    frame = pd.concat([df, file])
    return frame

def ROC_in(distance_m, distance_b):
    y_true = np.array(distance_b).flatten()
    y_scores = np.array(distance_m).flatten()
    fpr, tpr, thresholds = metrics.roc_curve(y_true, y_scores, pos_label = 1)
    AUC = roc_auc_score(y_true, y_scores)
    print(AUC)
    return 1-fpr, 1-tpr, thresholds, AUC

def send_to_roc(df):
    AUC_list = []

    methods = list(df.columns)

    for i in range(len(methods)):
        ROC_out = ROC_in(1-np.array(df[methods[1]]), list(df.Bin))
        np.savetxt(methods[i] +"_FPR.csv", ROC_out[0], delimiter=',')
        np.savetxt(methods[i] +"_TPR.csv", ROC_out[1], delimiter=',')
        np.savetxt(methods[i] +"_threshold.csv", ROC_out[2], delimiter=',')
        AUC = ROC_out[3]
    AUC_list.append(AUC)

    return AUC

if __name__ == '__main__':
    start_overall_time = time.time()   
    filepath = parse_commandline()
    df = import_files(filepath)
    outputs = send_to_roc(df)
    np.savetxt("CDR3_AUC.csv", outputs, delimiter=',')
    print("--- %s seconds ---" % (time.time()- start_overall_time))