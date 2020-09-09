#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 19:02:52 2020

@author: kristine
"""

#make downsampled file

"""
I will import each line of each matrix and concatenate them into a dataframe + downsample
"""
import pandas as pd
import numpy as np
from sklearn.utils import resample
from sys import argv 
import time 

def parse_commandline():  
    """
	commandline call: python3 Evaluation_2.py "path to distance chunks" "path to binary chunks" 
	
    EX:
    python3 log_reg.py '/home/kristine/Documents/Thesis/github_vers4/kristine_thesis_runthrough/distance_matrices' 'train_VDJ_mini.tsv' "/home/kristine/Documents/Thesis/github_vers4/kristine_thesis_runthrough/Binary_group" 3 20 2
    python3 generate_sub_cdr2.py 'work1/s111518/new_VDJ' 'Currated_VDJdb.tsv'  'work1/s111518/new_VDJ/Binary_group' 87 50 46
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
    CDR2_method_list = ["TM", "RMSD", "Blosum62", "hydrophobicity"]
    for j in range(len(CDR2_method_list)):
        y = import_lines("%s/CDR2dist_%s%s" % (filepath, CDR2_method_list[j], filename), line)
        empty.append(y)
    print("reading lines...")
    x = import_lines(binarypath+"/"+filename, line)
    empty.append(x)

    df = pd.DataFrame(empty)
    df = df.transpose()
    df.columns = CDR2_method_list + ["Bin"]
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
    d = {"TM": [1], "RMSD": [1], "Blosum62": [1], "hydrophobicity": [1], "Bin": [0]}
    df = pd.DataFrame(data=d)
    
    for i in range(chunk_len):
        total_df = getting_one_line_from_Each(filepath, filename, i, binarypath)    
        sub_sampled_df = reample_df(total_df)
        print("subsampling data")
        df = pd.concat([df, sub_sampled_df], ignore_index = True)
    
    return df

def input_all_chunks(path, filename, binarypath, chunks, chunk_len, remain):
    #skal starte med same fil i alle foldere
    d = {"TM": [1], "RMSD": [1], "Blosum62": [1], "hydrophobicity": [1], "Bin": [0]}
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

if __name__ == '__main__':
    start_overall_time = time.time()   
    path_dist, filename, path_binary, chunk_num, chunk_size, chunk_remain = parse_commandline()
    
    total = input_all_chunks(path_dist, filename, path_binary, chunk_num, chunk_size, chunk_remain)
#    log_data = log_Regress(total)    
        
    np.savetxt("CDR2_subsampled.csv", total, delimiter=',')
    print("--------------------------------")
    print("--- %s seconds ---" % (time.time()- start_overall_time))
    

