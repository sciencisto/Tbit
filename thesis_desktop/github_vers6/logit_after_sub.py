#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 07:42:21 2020

@author: kristine
"""

#make downsampled file

"""
I will import each line of each matrix and concatenate them into a dataframe + downsample
"""
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.utils import resample
import statsmodels.api as sm
from sys import argv 
import time 
from sklearn.metrics import accuracy_score

def parse_commandline():  
    """
	commandline call: python3 Evaluation_2.py "path to distance chunks" "path to binary chunks" 
	
    EX:
    python3 log_reg.py '/home/kristine/Documents/Thesis/github_vers4/kristine_thesis_runthrough/distance_matrices' 'train_VDJ_mini.tsv' "/home/kristine/Documents/Thesis/github_vers4/kristine_thesis_runthrough/Binary_group" 3 20 2

	"""
    path_dist = argv[1]

    return path_dist

def import_files(path_dist):
    d = {"Hamming": [1], "Blosum62": [1], "Blosum45": [1], "Pam10": [1], "K-Mer_3": [1], "K-Mer_4": [1], 
         "K-Mer_5": [1], "Atchley_Factors": [1], "hydrophobicity": [1], "Bin": [0]}
    
    df = pd.DataFrame(data=d)
    
    for i in range(5):
        fd = pd.read_csv(path_dist+"/"+i+"_CDR3_subsampled.csv")
        
        df_final = pd.concat([df, fd], ignore_index = True)
    
    return df_final

def log_Regress(total):
    
    y = total.Bin
    X = total.drop("Bin", axis=1)
    logit_model=sm.Logit(y,X)
    result=logit_model.fit()
    print("COEFICIENTS AND RESULTS")
    print(result.summary2())
    clf_2 = LogisticRegression().fit(X, y)
    pred_y_2 = clf_2.predict(X)
    print("Classes predicted")
    print( np.unique( pred_y_2 ) )
    print("accuray score:")
    print( accuracy_score(y, pred_y_2) )
    
    return result

if __name__ == '__main__':
    start_overall_time = time.time()   
    path_dist = parse_commandline()
    
    import_ = import_files(path_dist)
    log_data = log_Regress(import_)
    print("--------------------------------")
    print("--- %s seconds ---" % (time.time()- start_overall_time))
    

