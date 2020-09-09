 
#Calculate "truth" matrix
#import packages:
import numpy as np
from sys import argv 
import os
import multiprocessing as mp
import time 
import gc
from calc_CDR12_distance import create_CDR12_dict


##############################################################################
# 1 Parse commandline 
##############################################################################

def parse_commandline():  
    """
	commandline call: python3 eval_TCR_distance_binary.py input_filepath 
    input_filename chunk_size chunk_num "option"
    
    option is either exact or group. You can always choose exact but have to
    have a library to choose group. 
	
	"""
    wdir_path = argv[1]
    filtered_input_filename = argv[2] 
    chunk_size = int(float(argv[3])) 
    chunk_num = int(float(argv[4]))
    option = argv[5]

    return wdir_path, filtered_input_filename, chunk_size, chunk_num, option

###############################################################################
#Read in data
###############################################################################
def read_TCR_list(filtered_input_filename):
	"""
	Reads in the sample line by line; appends each TCR to a list 
	"""
	with open(filtered_input_filename,'r') as f:
		lines = [line.strip('\n').split('\t') for line in f]
	#
	Epitope = list()
	for line in lines: 
		Epitope.append(line[2])
	#    
	return lines, Epitope

def dist_i_list_parallel(TCR_i, Epitope, path, option):
    """return one line of the distance matrix for each i for the chosen 
    dist_option. For the moment there is a bunch of options. When i am finished
    There will be one or ONE combination of CDR1, 2 and CDR3"""
    
    if option == "exact":
        dist_i = np.array([dist_binary(TCR_i, cdr3) for cdr3 in Epitope])
    elif option == "group":
        lib = create_CDR12_dict(path+"/Libraries/Binary_epitope_lib.csv")
        dist_i = np.array([lib[TCR_i, cdr3] for cdr3 in Epitope])
    return dist_i

def dist_binary(TCR1_CDR3, TCR2_CDR3):
    if TCR1_CDR3 == TCR2_CDR3:
        distance = 1
    else:
        distance = 0
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
    wdir_path, filtered_input, chunk_size, chunk_num, option = parse_commandline()

    sample_name = os.path.splitext(os.path.basename(filtered_input))[0]
#	distance_function = choose_distance_function(distance_option)
    input_list = read_TCR_list(wdir_path + filtered_input)
    TCR_list = input_list[0]
    Epitope = input_list[1]
    
    TCR_chunk = list(chunks(TCR_list, chunk_size))[chunk_num]
    Epitope_chunk = list(chunks(Epitope, chunk_size))[chunk_num]
	#
    start_parallel_time = time.time()
	#
    cores = 8
    Que = list()
    if os.path.isfile(str(chunk_num) + "_" + sample_name + "_dmat.tsv"):
        os.remove(str(chunk_num) + "_" + sample_name + "_dmat.tsv")
	#	

    pool = mp.Pool(processes = cores, maxtasksperchild = 1000)
    for ii, TCR_i in enumerate(Epitope_chunk):
        print(ii, TCR_i)
        process = pool.apply_async(dist_i_list_parallel, 
                              (TCR_i, Epitope, wdir_path, option))
    
        Que.append(process)
		#del process
	
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
