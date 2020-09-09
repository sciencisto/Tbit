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
    python3 generate_sub.py 'work1/s111518/new_VDJ' 'Currated_VDJdb.tsv'  'work1/s111518/new_VDJ/Binary_group' 87 50 46
	"""
    path_dist = argv[1]
    filename = argv[2]
    path_binary = argv[3]
    chunk_num = int(float(argv[4]))
    chunk_size = int(float(argv[5]))
    chunk_remain = int(float(argv[6]))
    return path_dist, filename, path_binary, chunk_num, chunk_size, chunk_remain

def import_lines(filepath, i):
    with open(filepath, 'r') as fp:
        lines = [line.strip("\n").split("\t") for line in fp][i]
    the_line = np.array(lines)
    the_line = the_line.astype(np.float)
    return the_line

def getting_one_line_from_Each(filepath, filename, line, binarypath):
    empty = []
    CDR3_method_list = ["Hamming", "Blosum62", "Blosum45", "Pam10", "K-Mer_3", "K-Mer_4", "K-Mer_5", "Atchley_Factors", "hydrophobicity"]
    for j in range(len(CDR3_method_list)):
        y = import_lines("%s/CDR3dist_%s%s" % (filepath, CDR3_method_list[j], filename), line)
        empty.append(y)
    print("reading lines...")
    x = import_lines(binarypath+"/"+filename, line)
    empty.append(x)

    df = pd.DataFrame(empty)
    df = df.transpose()
    df.columns = CDR3_method_list + ["Bin"]
    return df

def reample_df(df):
    df_majority = df[df.Bin==0]
    df_minority = df[df.Bin==1]
    
    n = len(df_minority)
    
    #df_majority_downsampled = df_majority.sample(n=n, random_state = 1)
    df_majority_downsampled = resample(df_majority, 
                                 replace=False,    # sample without replacement
                                 n_samples=n,     # to match minority class
                                 random_state=123) # reproducible results
    df_downsampled = pd.concat([df_majority_downsampled, df_minority])
    return df_downsampled

def over_all_lines(filepath, filename, binarypath, chunk_len):
    d = {"Hamming": [1], "Blosum62": [1], "Blosum45": [1], "Pam10": [1], "K-Mer_3": [1], "K-Mer_4": [1], 
         "K-Mer_5": [1], "Atchley_Factors": [1], "hydrophobicity": [1], "Bin": [0]}
    df = pd.DataFrame(data=d)
    
    for i in range(chunk_len):
        total_df = getting_one_line_from_Each(filepath, filename, i, binarypath)    
        sub_sampled_df = reample_df(total_df)
        print("subsampling data")
        df = pd.concat([df, sub_sampled_df], ignore_index = True)
    
    return df

def input_all_chunks(path, filename, binarypath, chunks, chunk_len, remain):
    #skal starte med same fil i alle foldere
    d = {"Hamming": [1], "Blosum62": [1], "Blosum45": [1], "Pam10": [1], "K-Mer_3": [1], "K-Mer_4": [1], 
         "K-Mer_5": [1], "Atchley_Factors": [1], "hydrophobicity": [1], "Bin": [0]}
    df = pd.DataFrame(data=d)
    
    file_n = filename.split(".")[-2]
        
    for i in range(chunks):
        if i in range(chunks-1):
            new = over_all_lines(path,"/%s_"%(i)+file_n+"_dmat.tsv", binarypath, chunk_len)
            df_new = pd.concat([df, new], ignore_index = True)
        else:
            new1 = over_all_lines(path,"/%s_"%(i)+file_n+"_dmat.tsv", binarypath, remain)
            df_new1 = pd.concat([df, new1], ignore_index = True)
    df_final = pd.concat([df_new, df_new1], ignore_index = True)
    print(len(df_final))
    return df_final        


#def log_Regress(total):
    
#    y = total.Bin
#    X = total.drop("Bin", axis=1)
#    logit_model=sm.Logit(y,X)
#    result=logit_model.fit()
#    print("COEFICIENTS AND RESULTS")
#    print(result.summary2())
#    clf_2 = LogisticRegression().fit(X, y)
#    pred_y_2 = clf_2.predict(X)
#    print("Classes predicted")
#    print( np.unique( pred_y_2 ) )
#    print("accuray score:")
#    print( accuracy_score(y, pred_y_2) )
    
#    return result


if __name__ == '__main__':
    start_overall_time = time.time()   
    path_dist, filename, path_binary, chunk_num, chunk_size, chunk_remain = parse_commandline()
    
    total = input_all_chunks(path_dist, filename, path_binary, chunk_num, chunk_size, chunk_remain)
#    log_data = log_Regress(total)    
        
    np.savetxt("CDR3_subsampled.csv", total, delimiter=',')
    print("--------------------------------")
    print("--- %s seconds ---" % (time.time()- start_overall_time))
    

