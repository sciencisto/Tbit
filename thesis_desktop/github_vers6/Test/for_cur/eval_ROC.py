#Evaluation and quantification script

import matplotlib.pyplot as plt
import numpy as np
from sys import argv 
import time 

def parse_commandline():  
    """
    commandline call: python3 distribution_plot.py "concat_test/blos62/" "CDR1_blos" 
	python3 distribution_plot.py "/work1/s111518/CDR2dist_TM/Train/" "CDR2_TM"
	"""
    path = argv[1] #path to distance matrix
    
    
    return path
##############################################################################
#   CDR1 and CDR2
##############################################################################

def CDR1_TPR_graph(loop, type_, TM, Blosum62, RMSD, hydro): #, tpr_Hydrophobicity_list):
    #x axis  
    plt.figure(figsize=(10,7))
    plt.plot(np.genfromtxt('CDR1dist_RMSD_Threshold_list.csv', delimiter=","), np.append(RMSD, 1.0), '-o', c = "black", alpha = 0.7, label = "RMSD")    
    plt.plot(np.genfromtxt('CDR1dist_TM_Threshold_list.csv', delimiter=","), np.append(TM, 1.0), '-o', c = "darkred", alpha = 0.7, label = "TM")    
    plt.plot(np.genfromtxt('CDR1dist_Blosum62_Threshold_list.csv', delimiter=","), np.append(Blosum62, 1.0), '-o', c = "grey", alpha = 0.7, label = "Blosum62")
    plt.plot(np.genfromtxt('CDR1dist_hydrophobicity_Threshold_list.csv', delimiter=","), np.append(hydro, 1.0), '-o', c = "red", alpha = 0.7, label = "Hydrophobicity")
    plt.ylabel('True Positive Rate')
    plt.xlabel('Threshold')
    plt.title("%s True positive rate, Sensitivity of predictors" % (loop))
    plt.legend(loc='lower right');
    plt.savefig(type_+"_"+loop+"_TRP_graph.png", bbox_inches = 'tight')
    
def CDR1_specifi_graph(loop, type_, TM, Blosum62, RMSD, hydro): #, tpr_Hydrophobicity_list):
    #x axis
    plt.figure(figsize=(10,7))
    plt.plot(np.genfromtxt('CDR1dist_RMSD_Threshold_list.csv', delimiter=","), 1-np.append(RMSD, 1.0), '-o', c = "black", alpha = 0.7, label = "RMSD")    
    plt.plot(np.genfromtxt('CDR1dist_TM_Threshold_list.csv', delimiter=","), 1-np.append(TM, 1.0), '-o', c = "darkred", alpha = 0.7, label = "TM")    
    plt.plot(np.genfromtxt('CDR1dist_Blosum62_Threshold_list.csv', delimiter=","), 1-np.append(Blosum62, 1.0), '-o', c = "grey", alpha = 0.7, label = "Blosum62")
    plt.plot(np.genfromtxt('CDR1dist_hydrophobicity_Threshold_list.csv', delimiter=","), 1-np.append(hydro, 1.0), '-o', c = "red", alpha = 0.7, label = "Hydrophobicity")
    plt.ylabel('1-False Positive rate (Specificity)')
    plt.xlabel('Threshold')
    plt.title("%s Specificity of predictors, 1-False Positive rate" % (loop))
    plt.legend(loc='lower right');
    plt.savefig(type_+"_"+loop+"_Specificity_graph.png", bbox_inches = 'tight')

def CDR2_TPR_graph(loop, type_, TM, Blosum62, RMSD, hydro): #, tpr_Hydrophobicity_list):
    #x axis  
    plt.figure(figsize=(10,7))
    plt.plot(np.genfromtxt('CDR2dist_RMSD_Threshold_list.csv', delimiter=","), np.append(RMSD, 1.0), '-o', c = "black", alpha = 0.7, label = "RMSD")    
    plt.plot(np.genfromtxt('CDR2dist_TM_Threshold_list.csv', delimiter=","), np.append(TM, 1.0), '-o', c = "darkred", alpha = 0.7, label = "TM")    
    plt.plot(np.genfromtxt('CDR2dist_Blosum62_Threshold_list.csv', delimiter=","), np.append(Blosum62, 1.0), '-o', c = "grey", alpha = 0.7, label = "Blosum62")
    plt.plot(np.genfromtxt('CDR2dist_hydrophobicity_Threshold_list.csv', delimiter=","), np.append(hydro, 1.0), '-o', c = "red", alpha = 0.7, label = "Hydrophobicity")
    plt.ylabel('True Positive Rate')
    plt.xlabel('Threshold')
    plt.title("%s True positive rate, Sensitivity of predictors" % (loop))
    plt.legend(loc='lower right');
    plt.savefig(type_+"_"+loop+"_TRP_graph.png", bbox_inches = 'tight')
    
def CDR2_specifi_graph(loop, type_, TM, Blosum62, RMSD, hydro): #, tpr_Hydrophobicity_list):
    #x axis
    plt.figure(figsize=(10,7))
    plt.plot(np.genfromtxt('CDR2dist_RMSD_Threshold_list.csv', delimiter=","), 1-np.append(RMSD, 1.0), '-o', c = "black", alpha = 0.7, label = "RMSD")    
    plt.plot(np.genfromtxt('CDR2dist_TM_Threshold_list.csv', delimiter=","), 1-np.append(TM, 1.0), '-o', c = "darkred", alpha = 0.7, label = "TM")    
    plt.plot(np.genfromtxt('CDR2dist_Blosum62_Threshold_list.csv', delimiter=","), 1-np.append(Blosum62, 1.0), '-o', c = "grey", alpha = 0.7, label = "Blosum62")
    plt.plot(np.genfromtxt('CDR2dist_hydrophobicity_Threshold_list.csv', delimiter=","), 1-np.append(hydro, 1.0), '-o', c = "red", alpha = 0.7, label = "Hydrophobicity")
    plt.ylabel('1-False Positive rate (Specificity)')
    plt.xlabel('Threshold')
    plt.title("%s Specificity of predictors, 1-False Positive rate" % (loop))
    plt.legend(loc='lower right');
    plt.savefig(type_+"_"+loop+"_Specificity_graph.png", bbox_inches = 'tight')

   
def CDR1_CDR2_ROC_graph(loop, type_, TM_TPR, TM_FPR, blos_TPR, blos_FPR, RMSD_TPR, RMSD_FPR, 
                        hydro_TPR, hydro_FPR): 
    
    plt.figure(figsize=(10,7))
    plt.plot(np.array([0,1]), np.array([0,1]), '--', c = "blue", alpha = 0.6, label ="Random Guess")
    
    plt.plot(TM_FPR, TM_TPR, '-o', c = "darkred", alpha = 0.7, label = "TM")    
    plt.plot(blos_FPR, blos_TPR, '-o', c = "grey", alpha = 0.7, label = "Blosum62")    
    plt.plot(RMSD_FPR, RMSD_TPR, '-o', c = "black", alpha = 0.7, label = "RMSD")
    plt.plot(hydro_FPR, hydro_TPR, '-o', c = "red", alpha = 0.7, label = "Hydrophobicity")
    
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive rate')
    plt.title("%s Reciever operating characteristics curves" % (loop))
    plt.legend(loc='lower right');
    plt.savefig(type_+"_"+loop+"_ROC_graph.png", bbox_inches = 'tight')
    
##############################################################################
#   CDR3
##############################################################################    
    
def CDR3_TPR_graph(type_, k3, k4, k5, ham, atchley, blos62, blos45, pam10, hydro, immuno): #, tpr_Hydrophobicity_list):
    #x axis
    plt.figure(figsize=(10,7))
    
    plt.plot(np.genfromtxt('CDR3dist_K-Mer_3_Threshold_list.csv', delimiter=","), np.append(k3, 1.0), '-o', c = "mediumseagreen", alpha = 0.7, label = "K-mer, k=3") 
    plt.plot(np.genfromtxt('CDR3dist_K-Mer_4_Threshold_list.csv', delimiter=","), np.append(k4, 1.0), '-o', c = "green", alpha = 0.7, label = "K-mer, k=4")    
    plt.plot(np.genfromtxt('CDR3dist_Hamming_Threshold_list.csv', delimiter=","), np.append(ham, 1.0), '-o', c = "black", alpha = 0.7, label = "Hamming")
    plt.plot(np.genfromtxt('CDR3dist_hydrophobicity_Threshold_list.csv', delimiter=","), np.append(hydro, 1.0), '-o', c = "red", alpha = 0.7, label = "Hydrophobicity")
    plt.plot(np.genfromtxt('CDR3dist_Atchley_Factors_Threshold_list.csv', delimiter=","), np.append(atchley, 1.0), '-o', c = "grey", alpha = 0.7, label = "Atchley Factors")
    plt.plot(np.genfromtxt('CDR3dist_Blosum62_Threshold_list.csv', delimiter=","), np.append(blos62, 1.0), '-o', c = "darkred", alpha = 0.7, label = "Blosum 62")
    plt.plot(np.genfromtxt('CDR3dist_Blosum45_Threshold_list.csv', delimiter=","), np.append(blos45, 1.0), '-o', c = "teal", alpha = 0.7, label = "Blosum 45")
    plt.plot(np.genfromtxt('CDR3dist_Pam10_Threshold_list.csv', delimiter=","), np.append(pam10, 1.0), '-o', c = "darkturquoise", alpha = 0.7, label = "PAM 10")
    plt.plot(np.genfromtxt('CDR3dist_K-Mer_5_Threshold_list.csv', delimiter=","), np.append(k5, 1.0), c = 'darkgreen', alpha = 0.7, label = "K-mer, k=5")
    plt.plot(np.genfromtxt('CDR3dist_Immunomap_Threshold_list.csv', delimiter=","), np.append(immuno, 1.0), '-', c = "blue", alpha = 0.7, label = "Immunomap (Benchmark)")

    plt.ylabel('True Positive Rate (Sensitivity)')
    plt.xlabel('Distance Threshold')
    plt.title("CDR3 True positive rate, Sensitivity of predictors")
    plt.legend(loc='lower left');
    plt.savefig(type_+"_CDR3_TRP_graph.png", bbox_inches = 'tight')
    
def CDR3_speci_graph(type_, k3, k4, k5, ham, atchley, blos62, blos45, pam10, hydro, immuno): #, tpr_Hydrophobicity_list):
    #x axis
    plt.figure(figsize=(10,7))
    plt.plot(np.genfromtxt('CDR3dist_K-Mer_3_Threshold_list.csv', delimiter=","), 1-np.append(k3, 1.0), '-o', c = "mediumseagreen", alpha = 0.7, label = "K-mer, k=3") 
    plt.plot(np.genfromtxt('CDR3dist_K-Mer_4_Threshold_list.csv', delimiter=","), 1-np.append(k4, 1.0), '-o', c = "green", alpha = 0.7, label = "K-mer, k=4")    
    plt.plot(np.genfromtxt('CDR3dist_Hamming_Threshold_list.csv', delimiter=","), 1-np.append(ham, 1.0), '-o', c = "black", alpha = 0.7, label = "Hamming")
    plt.plot(np.genfromtxt('CDR3dist_hydrophobicity_Threshold_list.csv', delimiter=","), 1-np.append(hydro, 1.0), '-o', c = "red", alpha = 0.7, label = "Hydrophobicity")
    plt.plot(np.genfromtxt('CDR3dist_Atchley_Factors_Threshold_list.csv', delimiter=","), 1-np.append(atchley, 1.0), '-o', c = "grey", alpha = 0.7, label = "Atchley Factors")
    plt.plot(np.genfromtxt('CDR3dist_Blosum62_Threshold_list.csv', delimiter=","), 1-np.append(blos62, 1.0), '-o', c = "darkred", alpha = 0.7, label = "Blosum 62")
    plt.plot(np.genfromtxt('CDR3dist_Blosum45_Threshold_list.csv', delimiter=","), 1-np.append(blos45, 1.0), '-o', c = "teal", alpha = 0.7, label = "Blosum 45")
    plt.plot(np.genfromtxt('CDR3dist_Pam10_Threshold_list.csv', delimiter=","), 1-np.append(pam10, 1.0), '-o', c = "darkturquoise", alpha = 0.7, label = "PAM 10")
    plt.plot(np.genfromtxt('CDR3dist_K-Mer_5_Threshold_list.csv', delimiter=","), 1-np.append(k5, 1.0), c = 'darkgreen', alpha = 0.7, label = "K-mer, k=5")
    plt.plot(np.genfromtxt('CDR3dist_Immunomap_Threshold_list.csv', delimiter=","), 1-np.append(immuno, 1.0), '-', c = "blue", alpha = 0.7, label = "Immunomap (Benchmark)")

    plt.ylabel('1-False Positive rate (Specificity)')
    plt.xlabel('Distance Threshold')
    plt.title("CDR3 Specificity of predictors, 1-False Positive rate")
    plt.legend(loc='upper left');
    plt.savefig(type_+"_CDR3_specification_graph.png", bbox_inches = 'tight')

def ROC_graph(type_, CDR3_kmer3_TPR, CDR3_kmer3_FPR, CDR3_kmer4_TPR, CDR3_kmer4_FPR, CDR3_Hamming_TPR, 
              CDR3_Hamming_FPR, CDR3_Atchley_TPR, CDR3_Atchley_FPR, CDR3_blos62_FPR, CDR3_blos62_TPR, 
              CDR3_blos45_FPR, CDR3_blos45_TPR, CDR3_pam10_FPR, CDR3_pam10_TPR, CDR3_Hydro_FPR, 
              CDR3_Hydro_TPR, CDR3_Immunomap_grouped_TPR, CDR3_Immunomap_grouped_FPR, 
              CDR3_kmer5_grouped_TPR, CDR3_kmer5_grouped_FPR):    
    #graph
    plt.figure(figsize=(10,7))

    plt.plot(np.array([0,1]), np.array([0,1]), '--', c = "blue", alpha = 0.6, label ="Random Guess")
    
    plt.plot(CDR3_kmer3_FPR, CDR3_kmer3_TPR, '-o', c = "mediumseagreen", alpha = 0.7, label = "K-mer, k=3")    
    plt.plot(CDR3_kmer4_FPR, CDR3_kmer4_TPR, '-o', c = "green", alpha = 0.7, label = "K-mer, k=4")    
    plt.plot(CDR3_Hamming_FPR, CDR3_Hamming_TPR, '-o', c = "black", alpha = 0.7, label = "Hamming")
    plt.plot(CDR3_Hydro_FPR, CDR3_Hydro_TPR, '-o', c = "red", alpha = 0.7, label = "Hydrophobicity")
    plt.plot(CDR3_Atchley_FPR, CDR3_Atchley_TPR, '-o', c = "grey", alpha = 0.7, label = "Atchley")
    plt.plot(CDR3_blos62_FPR, CDR3_blos62_TPR, '-o', c = "darkred", alpha = 0.7, label = "Blosum 62")
    plt.plot(CDR3_blos45_FPR, CDR3_blos45_TPR, '-o', c = "teal", alpha = 0.7, label = "Blosum 45")
    plt.plot(CDR3_pam10_FPR, CDR3_pam10_TPR, '-o', c = "darkturquoise", alpha = 0.7, label = "PAM 10")
    plt.plot(CDR3_kmer5_grouped_FPR, CDR3_kmer5_grouped_TPR, '-o', c = 'darkgreen', alpha = 0.7, label = "K-mer, k=5")
    plt.plot(CDR3_Immunomap_grouped_FPR, CDR3_Immunomap_grouped_TPR, '-', c = "blue", alpha = 0.7, label = "Immunomap (Benchmark)")
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive rate')
    plt.title("Reciver Operator Characteristics curves for CDR3")
    plt.legend(loc='lower right');
    plt.savefig(type_+"_ROC_graph.png", bbox_inches = 'tight')

#####################################################################################################################
#   PATH TO THE GROUPED DATA    
#####################################################################################################################

if __name__ == '__main__':
    start_overall_time = time.time()   
    path = parse_commandline()
    type1 = path.split("/")[-1]
    type1 = type1.split("_")[-2]
    CDR1dist_Blosum62_FPR = np.genfromtxt(path+'/CDR1dist_Blosum62_FPR.csv', delimiter=",")
    CDR1dist_Blosum62_TPR = np.genfromtxt(path+'/CDR1dist_Blosum62_TPR.csv', delimiter=",")
    CDR1dist_RMSD_FPR = np.genfromtxt(path+'/CDR1dist_RMSD_FPR.csv', delimiter=",")
    CDR1dist_RMSD_TPR = np.genfromtxt(path+'/CDR1dist_RMSD_TPR.csv', delimiter=",")
    CDR1dist_TM_FPR = np.genfromtxt(path+'/CDR1dist_TM_FPR.csv', delimiter=",")
    CDR1dist_TM_TPR = np.genfromtxt(path+'/CDR1dist_TM_TPR.csv', delimiter=",")
    CDR1dist_hydrophobicity_FPR = np.genfromtxt(path+'/CDR1dist_hydrophobicity_FPR.csv', delimiter=",")
    CDR1dist_hydrophobicity_TPR = np.genfromtxt(path+'/CDR1dist_hydrophobicity_TPR.csv', delimiter=",")
    CDR2dist_Blosum62_FPR = np.genfromtxt(path+'/CDR2dist_Blosum62_FPR.csv', delimiter=",")
    CDR2dist_Blosum62_TPR = np.genfromtxt(path+'/CDR2dist_Blosum62_TPR.csv', delimiter=",")
    CDR2dist_RMSD_FPR = np.genfromtxt(path+'/CDR2dist_RMSD_FPR.csv', delimiter=",")
    CDR2dist_RMSD_TPR = np.genfromtxt(path+'/CDR2dist_RMSD_TPR.csv', delimiter=",")
    CDR2dist_TM_FPR = np.genfromtxt(path+'/CDR2dist_TM_FPR.csv', delimiter=",")
    CDR2dist_TM_TPR = np.genfromtxt(path+'/CDR2dist_TM_TPR.csv', delimiter=",")
    CDR2dist_hydrophobicity_FPR = np.genfromtxt(path+'/CDR2dist_hydrophobicity_FPR.csv', delimiter=",")
    CDR2dist_hydrophobicity_TPR = np.genfromtxt(path+'/CDR2dist_hydrophobicity_TPR.csv', delimiter=",")
    CDR3dist_Atchley_Factors_FPR = np.genfromtxt(path+'/CDR3dist_Atchley_Factors_FPR.csv', delimiter=",")
    CDR3dist_Atchley_Factors_TPR = np.genfromtxt(path+'/CDR3dist_Atchley_Factors_TPR.csv', delimiter=",")
    CDR3dist_Blosum45_FPR = np.genfromtxt(path+'/CDR3dist_Blosum45_FPR.csv', delimiter=",")
    CDR3dist_Blosum45_TPR = np.genfromtxt(path+'/CDR3dist_Blosum45_TPR.csv', delimiter=",")
    CDR3dist_Blosum62_FPR = np.genfromtxt(path+'/CDR3dist_Blosum62_FPR.csv', delimiter=",")
    CDR3dist_Blosum62_TPR = np.genfromtxt(path+'/CDR3dist_Blosum62_FPR.csv', delimiter=",")
    CDR3dist_Hamming_FPR = np.genfromtxt(path+'/CDR3dist_Hamming_FPR.csv', delimiter=",")
    CDR3dist_Hamming_TPR = np.genfromtxt(path+'/CDR3dist_Hamming_TPR.csv', delimiter=",")
    CDR3dist_Immunomap_FPR = np.genfromtxt(path+'/CDR3dist_Immunomap_FPR.csv', delimiter=",")
    CDR3dist_Immunomap_TPR = np.genfromtxt(path+'/CDR3dist_Immunomap_TPR.csv', delimiter=",")
    CDR3dist_KMer_3_FPR = np.genfromtxt(path+'/CDR3dist_K-Mer_3_FPR.csv', delimiter=",")
    CDR3dist_KMer_3_TPR = np.genfromtxt(path+'/CDR3dist_K-Mer_3_TPR.csv', delimiter=",")
    CDR3dist_KMer_4_FPR = np.genfromtxt(path+'/CDR3dist_K-Mer_4_FPR.csv', delimiter=",")
    CDR3dist_KMer_4_TPR = np.genfromtxt(path+'/CDR3dist_K-Mer_4_TPR.csv', delimiter=",")
    CDR3dist_KMer_5_FPR = np.genfromtxt(path+'/CDR3dist_K-Mer_5_FPR.csv', delimiter=",")
    CDR3dist_KMer_5_TPR = np.genfromtxt(path+'/CDR3dist_K-Mer_5_TPR.csv', delimiter=",")
    CDR3dist_Pam10_FPR = np.genfromtxt(path+'/CDR3dist_Pam10_FPR.csv', delimiter=",")
    CDR3dist_Pam10_TPR = np.genfromtxt(path+'/CDR3dist_Pam10_TPR.csv', delimiter=",")
    CDR3dist_hydrophobicity_FPR = np.genfromtxt(path+'/CDR3dist_hydrophobicity_FPR.csv', delimiter=",")
    CDR3dist_hydrophobicity_TPR = np.genfromtxt(path+'/CDR3dist_hydrophobicity_TPR.csv', delimiter=",")
    
    CDR1_TPR_graph("CDR1", type1, CDR1dist_TM_TPR, CDR1dist_Blosum62_TPR, CDR1dist_RMSD_TPR, CDR1dist_hydrophobicity_TPR)
    CDR1_specifi_graph("CDR1", type1, CDR1dist_TM_FPR, CDR1dist_Blosum62_FPR, CDR1dist_RMSD_FPR, CDR1dist_hydrophobicity_FPR)
    CDR1_CDR2_ROC_graph("CDR1", type1, CDR1dist_TM_TPR, CDR1dist_TM_FPR, CDR1dist_Blosum62_TPR, CDR1dist_Blosum62_FPR, 
                        CDR1dist_RMSD_TPR, CDR1dist_RMSD_FPR, CDR1dist_hydrophobicity_TPR, CDR1dist_hydrophobicity_FPR)

    CDR2_TPR_graph("CDR2", type1, CDR2dist_TM_TPR, CDR2dist_Blosum62_TPR, CDR2dist_RMSD_TPR, CDR2dist_hydrophobicity_TPR)
    CDR2_specifi_graph("CDR2", type1, CDR2dist_TM_FPR, CDR2dist_Blosum62_FPR, CDR2dist_RMSD_FPR, CDR2dist_hydrophobicity_FPR)
    CDR1_CDR2_ROC_graph("CDR2", type1, CDR2dist_TM_TPR, CDR2dist_TM_FPR, CDR2dist_Blosum62_TPR, CDR2dist_Blosum62_FPR, 
                        CDR2dist_RMSD_TPR, CDR2dist_RMSD_FPR, CDR2dist_hydrophobicity_TPR, CDR2dist_hydrophobicity_FPR)

    CDR3_TPR_graph(type1, CDR3dist_KMer_3_TPR, CDR3dist_KMer_4_TPR, CDR3dist_KMer_5_TPR, CDR3dist_Hamming_TPR, 
                   CDR3dist_Atchley_Factors_TPR, CDR3dist_Blosum62_TPR, CDR3dist_Blosum45_TPR, CDR3dist_Pam10_TPR, 
                   CDR3dist_hydrophobicity_TPR, CDR3dist_Immunomap_TPR)
    CDR3_speci_graph(type1, CDR3dist_KMer_3_FPR, CDR3dist_KMer_4_FPR, CDR3dist_KMer_5_FPR, CDR3dist_Hamming_FPR, 
                   CDR3dist_Atchley_Factors_FPR, CDR3dist_Blosum62_FPR, CDR3dist_Blosum45_FPR, CDR3dist_Pam10_FPR, 
                   CDR3dist_hydrophobicity_FPR, CDR3dist_Immunomap_FPR)
    ROC_graph(type1, CDR3dist_KMer_3_TPR, CDR3dist_KMer_3_FPR, CDR3dist_KMer_4_TPR, CDR3dist_KMer_4_FPR, CDR3dist_Hamming_TPR, 
              CDR3dist_Hamming_FPR, CDR3dist_Atchley_Factors_TPR, CDR3dist_Atchley_Factors_FPR, CDR3dist_Blosum62_TPR, 
              CDR3dist_Blosum62_TPR, CDR3dist_Blosum45_FPR, CDR3dist_Blosum45_TPR, CDR3dist_Pam10_FPR, CDR3dist_Pam10_TPR, 
              CDR3dist_hydrophobicity_FPR, CDR3dist_hydrophobicity_TPR, CDR3dist_Immunomap_TPR, CDR3dist_Immunomap_FPR, 
              CDR3dist_KMer_5_TPR, CDR3dist_KMer_5_FPR)
    print("--- %s seconds ---" % (time.time()- start_overall_time))

