#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 11:53:55 2020

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
	python3 ROC_with_package.py "/work1/s111518/new_VDJ/CDR1dist_RMSD/" "/work1/s111518/new_VDJ/Binary_group/"
	"""
    path_dist = argv[1] #path to distance matrix
    path_bin = argv[2] #type of calculation

    return path_dist, path_bin


def import_files(path):
    """Function that imports all the chuncks and return them as one big numpy array """
    all_files = glob.glob(path + "/*.tsv")
    li = []
    for i in range(len(all_files)):
        print(all_files[i])
        file = pd.read_csv(all_files[i], index_col=None, sep="\t", header = None)
        li.append(file)
    frame = np.concatenate((li), axis=0)
    del li 
    del file 
    return frame

def ROC_in(distance_m, distance_b):
    y_true = np.array(distance_b).flatten()
    y_scores = np.array(distance_m).flatten()
    fpr, tpr, thresholds = metrics.roc_curve(y_true, y_scores, pos_label = 1)
    AUC = roc_auc_score(y_true, y_scores)
    print(AUC)
    return 1-fpr, 1-tpr, thresholds, AUC

if __name__ == '__main__':
    start_overall_time = time.time()   
    path_dist, path_bin = parse_commandline()
    distance_m = import_files(path_dist)
    distance_b = import_files(path_bin)
    ROC_out = ROC_in(1-distance_m, distance_b)
    name = path_dist.split("/")[-2]
    np.savetxt(name +"_FPR.csv", ROC_out[0], delimiter=',')
    np.savetxt(name +"_TPR.csv", ROC_out[1], delimiter=',')
    np.savetxt(name +"_threshold.csv", ROC_out[2], delimiter=',')
    np.savetxt(name +"_AUC.csv", ROC_out[3], delimiter=',')
    print("--- %s seconds ---" % (time.time()- start_overall_time))
