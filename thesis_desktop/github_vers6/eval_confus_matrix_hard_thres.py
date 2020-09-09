import glob as glob
import numpy as np
from sklearn.metrics import confusion_matrix
from sys import argv 
import time 


def parse_commandline():  
    """
	commandline call: python3 Evaluation_2.py "path to distance chunks" "path to binary chunks" 
	
    EX:
    python3 Evaluation_2.py '/home/kristine/Documents/Thesis/github_vers4/concat_test/blos62/' '/home/kristine/Documents/Thesis/github_vers4/concat_test/binary/'
    python3 Evaluation_2.py '/work1/s111518/CDR3dist_kmer5/Train/' '/work1/s111518/Binary/Train_group/'
	"""
    path_dist = argv[1] #path to distance matrix
    path_bin = argv[2] #path to binary, truth matrix
    path_output = argv[3]

    return path_dist, path_bin, path_output

def import_matrix_new(dist_matrix, bin_path, recgon):
    """Function that import the distance matrix chunk and the binary distance matrix chunk"""
    #open tsv file of a specific sorted name
    print("Finding distance and matching binany chunks.....")
    distance_matrix_chunck = np.genfromtxt(dist_matrix, delimiter="\t") 
    binary_matrix_chunck = np.genfromtxt(bin_path+"/%s"%(recgon), delimiter="\t")
    #return both files
    return distance_matrix_chunck, binary_matrix_chunck

def _confusion_matrix(distance_matrix, binary_matrix, theshold):
    """
    Function that create a confusion matrix for each distance file
    
              |        True Condition           |
    ----------|---------------------------------|
    Predicted | True postive   | False positive |
    condition |----------------|----------------|
              | False positive | True Negative  |
      
    """
    
    _true = np.array(binary_matrix).flatten()
    distance_array = np.array(distance_matrix).flatten()
    _predicted = np.zeros(len(distance_array))

    for i in range(len(distance_array)):
        if distance_array[i] <= theshold:
            _predicted[i] = 1
    print("Calculating intermidiate confusion matrix for threshold: %s" % (theshold))
    tn, fp, fn, tp = confusion_matrix(_true, _predicted).ravel()
    confus_matrix =  np.array([[tp, fp],[fn, tn]])        
    return confus_matrix

def TPR_curve(distance_matrix, binary_matrix, theshold_list):
    """Function that outputs a csv file with all the TPRs of the distance matrix"""
    confus_mat_list = []
    
    for i in range(len(theshold_list)):
        cm  = _confusion_matrix(distance_matrix, binary_matrix, theshold_list[i])
        confus_mat_list.append(cm)
    
    return confus_mat_list

def make_confus(path_dist, path_bin, theshold_list):
    """Function that loop all the data into the confusion matrix generator"""
    empty_list = []    
    all_files = sorted(glob.glob(path_dist + "/*.tsv"))
    chunks = (len(all_files)) 
    
    for i in range(len(all_files)):
        print("-------------------------------------------------------------------------------")
        print("Chunk %s/%s" % (i, len(all_files)))
        print("The confusion matrices along the threshold vector is being calculated for %s" % (all_files[i]))
        print("-------------------------------------------------------------------------------")
        codon = all_files[i].split("/")[-1]
        imported_files = import_matrix_new(all_files[i], path_bin, codon)
        _confus_matrix = TPR_curve(imported_files[0], imported_files[1], theshold_list)
        empty_list.append(_confus_matrix)

    return empty_list, chunks

def confus_concat(total_confus, thes_num, chunks):
    """Function that concatenate the confusion matrix for each threshold"""
    x = []
    for i in range(chunks):
        mat = total_confus[i][thes_num]
        x.append(mat)
    confus_matrix = sum(x)
    print("Confusion matrix: [[%s, %s], [%s, %s]]" % (confus_matrix[0][0], confus_matrix[0][1], confus_matrix[1][0], confus_matrix[1][1]))
    #TPR defines how many correct positive results occur among all positive samples available during the test.
    TPR = confus_matrix[0][0]/(confus_matrix[0][0]+confus_matrix[1][0])
    #FPR defines how many incorrect positive results occur among all negative samples available during the test.
    FPR = confus_matrix[0][1]/(confus_matrix[0][1]+confus_matrix[1][1])
    print("The true positive rate become: %s" % (TPR))
    print("The true neagtive rate become: %s" % (FPR))
    return TPR, FPR

def confus_concat1(total_confus, theshold_list, chunks):
    tpr =[]
    fpr =[]
    for j in range(len(theshold_list)-1):
        print("--------------------------------")
        print("For threshold: %s" % (theshold_list[j]))
        print("--------------------------------")
        y = confus_concat(total_confus, j, int(chunks-1))
        tpr.append(y[0])
        fpr.append(y[1])
    return tpr, fpr
    
if __name__ == '__main__':
    start_overall_time = time.time()   
    path_dist, path_bin, path_output = parse_commandline()
    print("Calculating the confusion matrixes and True Positive rates")
    print("#################################################################################")
    
    #
    theshold_list = [0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.01]
    #
    #name = path_dist.split("/")[-2]

    #theshold_list = np.genfromtxt(name+"_Threshold_list.csv", delimiter=",")
    #theshold_list = list(theshold_list)
    
    print("Theshold list is: %s" % (theshold_list))
    total_confus = make_confus(path_dist, path_bin, theshold_list)
    print("#################################################################################")
    print("The final matrixes along the theshold vector along all chunks are:")
    print("#################################################################################")
    TPR_list = confus_concat1(total_confus[0], theshold_list, total_confus[1])
    print("--------------------------------")
    print("Therefore the TPR list is: %s" % (TPR_list[0]))
    print("And the FPR list is: %s" % (TPR_list[1]))
    name_ = path_dist.split("/")[-2]
    np.savetxt(("%s/%s_TPR.csv"%(path_output, name_)), TPR_list[0], delimiter=',')
    np.savetxt(("%s/%s_FPR.csv"%(path_output, name_)), TPR_list[1], delimiter=',')
    print("--------------------------------")
    print("--- %s seconds ---" % (time.time()- start_overall_time))
    
