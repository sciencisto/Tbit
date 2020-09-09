import pandas as pd
import glob as glob
import numpy as np
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

def quanties(distance_matrix, name):
    a = np.quantile(distance_matrix, 0.05)    
    b = np.quantile(distance_matrix, 0.10)
    c = np.quantile(distance_matrix, 0.20)
    d = np.quantile(distance_matrix, 0.25)
    e = np.quantile(distance_matrix, 0.30)
    f = np.quantile(distance_matrix, 0.40)
    g = np.quantile(distance_matrix, 0.50)
    h = np.quantile(distance_matrix, 0.60)
    i = np.quantile(distance_matrix, 0.70)
    j = np.quantile(distance_matrix, 0.75)
    k = np.quantile(distance_matrix, 0.80)
    l = np.quantile(distance_matrix, 0.90)
    m = np.quantile(distance_matrix, 1.0)
    
    threshold_list = np.array([a, b, c, d, e, f, g, h, i, j, k, l, m])
    np.savetxt(("%s_Threshold_list.csv"%(name)), threshold_list, delimiter=',')
    return threshold_list 
    
if __name__ == '__main__':
    start_overall_time = time.time()   
    path_dist, name = parse_commandline()
    distance_m = import_files(path_dist)
    quanties(distance_m, name)  
    print("--- %s seconds ---" % (time.time()- start_overall_time))
