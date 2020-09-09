#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 15:16:11 2020

@author: kristine
"""

#Distribution script

import pandas as pd
import glob as glob
import matplotlib.pyplot as plt
import numpy as np
#import seaborn as sns
from sys import argv 
import time 


def parse_commandline():  
    """
    Rembember to login to databar by ssh -X.
	commandline call: python3 distribution_plot.py "concat_test/blos62/" "CDR1_blos" 
	python3 distribution_plot.py "/work1/s111518/CDR2dist_TM/Train/" "CDR2_TM"
	"""
    path_dist = argv[1] #path to distance matrix
    name = argv[2] #type of calculation

    return path_dist, name

    
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

def distribution(distance_matrix, name):
    """Function that print two immages desrcibing the distrbution of the distance matrix"""
    print("Making the plain distance distribution plot")
    plt.figure(figsize=(8,6))
    plt.hist(np.array(distance_matrix).flatten(), bins=25, color='grey',  edgecolor='black', linewidth=1, alpha=0.7)
    plt.ylabel('Number of TCRs')
    plt.xlabel('Pairwise %s Distance' % (name))
    plt.savefig(name+"_distance_distribution.png", bbox_inches = 'tight')
  
if __name__ == '__main__':
    start_overall_time = time.time()   
    path_dist, name = parse_commandline()
    distance_m = import_files(path_dist)
    distribution(distance_m, name)    
    print("--- %s seconds ---" % (time.time()- start_overall_time))
