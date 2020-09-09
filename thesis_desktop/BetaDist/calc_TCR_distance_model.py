#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 13:52:30 2020

@author: kristine
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 20:46:24 2020

@author: kristine
"""

"""
This script contains functions to compare TCR beta chains. 
It imports a file with the following format and output a complete distance 
matrix.

Input file format:
------------------------------
    INPUT:  .tsv file containing In frame cdr3 seq, V-gene and Epitope:HLA gene 
    
    CDR3 sequence \t V-gene \t Epitope:HLA gene
    
    Example: 
        
    CATSNDRDLDEQFF	TRBV15_01	RLRPGGKKK:HLA-A_03:01
    CASSYGQAYQPQHF	TRBV06-1_01	ARMILMTHF:HLA-B_27
    CASSLAVSSEQFF	TRBV11-2_01	ARMILMTHF:HLA-B_27
    CSVGGTGTYGYTF	TRBV29-1_01	ARMILMTHF:HLA-B_27
    CASSIVGHMNTEAFF	TRBV19_01	ARMILMTHF:HLA-B_27
    CAISDGAGTETQYF	TRBV10-3_01	GTSGSPIIDK:HLA-A_11:01

Theis script can run without the epitopes. They are used at a later stage and 
user can settle for only imputing the file once by including the epitopes.

Dependencies:
------------------------------
The scripts: calc_CDR12_distance.py, calc_CDR3_distance.py and 
calc_substitution_matrices

Modules: Biopython, operator & multiprocessing.

Functionality
---------------------------------
This script take a chunk of a file and generate a distance matrix based on the 
chosen methods. Both are administered in the script: Running_calc_and_eval.py

Output:
------------------------------
cc_"sample_name"_distance_matrix.tsv file distances of 50 TCRs against 
all TCRs

USAGE: python3 TCR_distance.py 'filepath to input file' 'input file name' 
chunk_size chunk_number "CDR3 distance method" "CDR1 distance method" 
"CDR2 distance method"
"""
  
#Calculate distance matrix
from calc_CDR12_distance import create_CDR12_dict
from calc_CDR3_distance import dist_CDR3_pairwise_ham, CDR3_diagonal_dict, dist_CDR3_alignment, kmers_dist, Atchley_euclidean_dist, make_AA_to_Atchley_dict, Immunomap_Dist, physiocehmical_prop, vol_distance
from calc_substitution_matrices import Blosum62, pam10, Blosum45

#import packages:
import numpy as np
from sys import argv 
import os
import multiprocessing as mp
import time 
import gc
import operator

##############################################################################
# 1 Parse commandline 
##############################################################################

def parse_commandline():  
    """
	commandline call: python3 TCR_distance.py 'filepath to input file' 
    'input file name' chunk_size chunk_number "CDR3 distance method" 
    "CDR1 distance method" "CDR2 distance method"
    ---------------------------------------------------------------------------    
    'filepath to input file': '//'
    
    'input file name': '..'
    
    chunk_size: a float, no ''
    chunk_number: a float no ''
    
    CDR3 distance method options:   "Hamming", "Blosum62", "Blosum45", 
                                    "Pam10", "K-Mer_3", "K-Mer_4", "K-Mer_5", 
                                    "Atchley_Factors", "hydrophobicity", 
                                    "Immunomap", "none"
    CDR1 distance method options:   "TM", "RMSD", "Blosum62", 
                                    "hydrophobicity", "none"
    CDR2 distance method options:   "TM", "RMSD", "Blosum62", "hydrophobicity"
                                    "none"
	"""
    wdir_path = argv[1] 
    filtered_input_filename = argv[2] 
    chunk_size = int(float(argv[3])) 
    chunk_num = int(float(argv[4])) 
    return wdir_path, filtered_input_filename, chunk_size, chunk_num

###############################################################################
#Read in data
###############################################################################
def read_TCR_list(filtered_input_filename):
	"""
	Reads in the sample line by line; appends each TCR to a list 
	"""
	with open(filtered_input_filename,'r') as f:
		lines = [line.strip('\n').split('\t') for line in f]
	# generating a list only of CDR3s to be used to calculate the CDR3 distance
	CDR3 = list()
	for line in lines: 
		CDR3.append(line[0])
        
	# generating a list with only V-genes to calulcate CDR1 and CDR2 distance   
	V = list()
	for line in lines: 
		V.append(line[1])
	#    
	return lines, CDR3, V

def dist_i_list_parallel(TCR_i, CDR3_list, V_list, path):
    """return one line of the distance matrix for each i for the chosen 
    dist_option. For the moment there is a bunch of options. When i am finished
    There will be one or ONE combination of CDR1, 2 and CDR3"""
    
    CDR1_TM = create_CDR12_dict(path+"Libraries/dist_CDR1_TM_1.csv")
    dist_CDR1_TM = np.array([CDR1_TM[TCR_i[1], v] for v in V_list])*0.17
    
    CDR1_B62 = create_CDR12_dict(path+"Libraries/CDR1_blos62.csv")
    dist_CDR1_B62 = np.array([CDR1_B62[TCR_i[1], v] for v in V_list])*0.15
            
    CDR2_RMSD = create_CDR12_dict(path+"Libraries/CDR2_RMSD.csv")
    dist_CDR2_RMSD = np.array([CDR2_RMSD[TCR_i[1], v] for v in V_list])*0.21
    
    CDR2_TM = create_CDR12_dict(path+"Libraries/dist_CDR2_TM_1.csv")
    dist_CDR2_TM = np.array([CDR2_TM[TCR_i[1], v] for v in V_list])*0.28
           
    diagonal10 = CDR3_diagonal_dict(CDR3_list, pam10, -30, -0.5)
    dist_CDR3_pam10 = np.array([dist_CDR3_alignment(TCR_i[0], cdr3, pam10, -30, -0.5, 
                        diagonal10, operator.ge)/1.5 for cdr3 in CDR3_list])*0.19
    
    distance = (dist_CDR1_TM + dist_CDR1_B62 + dist_CDR2_RMSD +
                dist_CDR2_TM + dist_CDR3_pam10)

    return distance
    
def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def write_out_dist_i_chunk(chunk_enumerator, sample_name, distance_line):
	"""
	Adds a distance row for each TCR_i in the format X X X d(TCR_i, TCR_j)
	"""
	with open(chunk_enumerator + "_" + sample_name + "_dmat.tsv", 'a+b') as f:
		np.savetxt(f, distance_line, delimiter = "", newline = "\t")
		last_tab_in_line_size_in_bytes = -1
		f.seek(last_tab_in_line_size_in_bytes, 2)
		f.truncate()
		f.write(b"\n")
	del distance_line

################################################################################

if __name__ == '__main__':
    start_overall_time = time.time()
    
    #what will be parsed from comandline in time:
    wdir_path, filtered_input, chunk_size, chunk_num = parse_commandline()
    print("#################################################################################")
    
    # importing the sample file and making lists
    print("Reading in the sample file...")
    sample_name = os.path.splitext(os.path.basename(filtered_input))[0]
    input_list = read_TCR_list(wdir_path + filtered_input)
    TCR_list = input_list[0]
    CDR3_list = input_list[1]
    V_list = input_list[2]

    print("Administering the number of chunks...")
    TCR_chunk = list(chunks(TCR_list, chunk_size))[chunk_num]
    CDR3_chunk = list(chunks(CDR3_list, chunk_size))[chunk_num]
    V_chunk = list(chunks(V_list, chunk_size))[chunk_num]
	#
    print("...initiating parralellization...")
    start_parallel_time = time.time()
	#
    cores = 8
    Que = list()
    print("deleting current files at location")
    if os.path.isfile(str(chunk_num) + "_" + sample_name + "_dmat.tsv"):
        os.remove(str(chunk_num) + "_" + sample_name + "_dmat.tsv")
	#	

    pool = mp.Pool(processes = cores, maxtasksperchild = 1000)
    print("Initilizing the loop that calculate each line of the chunk")
    for ii, TCR_i in enumerate(TCR_chunk):
        print(ii, TCR_i)
        process = pool.apply_async(dist_i_list_parallel, 
                              (TCR_i, CDR3_list, V_list, wdir_path))
    
        Que.append(process)

    print("printing the chunk into a file")
    for oo in range(len(Que)):
        function_output = Que[oo].get()
        write_out_dist_i_chunk(str(chunk_num), sample_name, function_output)
        del function_output
    pool.terminate()
    del Que
    gc.collect()
    print("parallel calculating time")
    print("--- %s seconds ---" % (time.time() - start_parallel_time))
    print("overall used time")
    print("--- %s seconds ---" % (time.time()- start_overall_time))
