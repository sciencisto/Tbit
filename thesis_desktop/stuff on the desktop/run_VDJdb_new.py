filepath = "/work1/s111518/new_VDJ"
filename = "Currated_VDJdb.tsv"

CDR12_method_list = ["TM", "RMSD", "Blosum62", "hydrophobicity"]
CDR3_method_list = ["Hamming", "Blosum62", "Blosum45", "Pam10", "K-Mer_3", "K-Mer_4", "K-Mer_5", "Atchley_Factors", "hydrophobicity", "Immunomap"]
Binary_method_list = ["exact", "group"]

print("echo '########################################################'")      
print("echo 'Starting the TCRdistance program, Kristines Thesis 2020'")
print("echo '########################################################'")      

      
##########################################################################
#firstly i count the number of chunks i want.    
#import pandas as pd
#import math
#TSV = pd.read_csv(filepath+"/"+filename, sep="\t", header=None)
length = 87
#length = math.ceil(len(TSV)/50)
rest = 46
#rest = len(TSV)-(length-1)*50
#del TSV

#calculate the dmats for all methdods and moving them into their designated directories      
for j in range(len(CDR12_method_list)):
    print("echo 'Calculating dmats for the %s method for CDR1'" % (CDR12_method_list[j]))
    print("mkdir CDR1dist_%s" % (CDR12_method_list[j]))
    for i in range(length):
        print("python3 calc_TCR_distance.py '%s/' '%s' 50 %s 'none' '%s' 'none'" % (filepath, filename, i, CDR12_method_list[j]))
    print("mv *dmat.tsv CDR1dist_%s" % (CDR12_method_list[j]))
for j in range(len(CDR12_method_list)):
    print("echo 'Calculating dmats for the %s method for CDR2'" % (CDR12_method_list[j]))
    print("mkdir CDR2dist_%s" % (CDR12_method_list[j]))
    for i in range(length):
        print("python3 calc_TCR_distance.py '%s/' '%s' 50 %s 'none' 'none' '%s'" % (filepath, filename, i, CDR12_method_list[j]))
    print("mv *dmat.tsv CDR2dist_%s" % (CDR12_method_list[j]))
"""
for j in range(len(CDR3_method_list)):
    print("echo 'Calculating dmats for the %s method for CDR3'" % (CDR3_method_list[j]))
    print("mkdir CDR3dist_%s" % (CDR3_method_list[j]))
    for i in range(length):
        print("python3 calc_TCR_distance.py '%s/' '%s' 50 %s '%s' 'none' 'none'" % (filepath, filename, i, CDR3_method_list[j]))
    print("mv *dmat.tsv CDR3dist_%s" % (CDR3_method_list[j]))
    
print("echo '########################################################'")      
print("echo 'all dmats have been calculated, now the binary files are generated'")
print("echo '########################################################'")      

for j in range(len(Binary_method_list)):
    print("echo 'Calculating binary dmat for the %s match'" % (Binary_method_list[j]))
    print("mkdir Binary_%s" % (Binary_method_list[j]))
    for i in range(length):
        print("python3 eval_TCR_distance_binary.py '%s/' '%s' 50 %s '%s'" % (filepath, filename, i, Binary_method_list[j]))
    print("mv *dmat.tsv Binary_%s" % (Binary_method_list[j]))

print("echo '########################################################'")      
print("echo 'Making distribution plots for all'")
print("echo '########################################################'")      
print("mkdir distribution_plots")
for j in range(len(CDR12_method_list)):
    print("python3 eval_distribution_plot.py 'CDR1dist_%s' 'CDR1_%s'" % (CDR12_method_list[j], CDR12_method_list[j]))
    print("python3 eval_distribution_plot.py 'CDR2dist_%s' 'CDR2_%s'" % (CDR12_method_list[j], CDR12_method_list[j]))
for j in range(len(CDR3_method_list)):
    print("python3 eval_distribution_plot.py 'CDR3dist_%s' 'CDR3_%s'" % (CDR3_method_list[j], CDR3_method_list[j]))
for j in range(len(Binary_method_list)):
    print("python3 eval_distribution_plot.py 'Binary_%s' 'binary_%s'" % (Binary_method_list[j], Binary_method_list[j]))
print("mv *.png distribution_plots")

print("echo '########################################################'")      
print("echo 'Generation of confusion matrices and TPR & FPR'")
print("echo '########################################################'")     
print("mkdir ROC_exact_matches")
print("mkdir ROC_group_matches")

#a little piece of hard coding for shittyness      
all_folders = ['CDR1dist_TM','CDR1dist_RMSD','CDR1dist_Blosum62','CDR1dist_hydrophobicity','CDR2dist_TM',
               'CDR2dist_RMSD','CDR2dist_Blosum62','CDR2dist_hydrophobicity', 'CDR3dist_Hamming',
               'CDR3dist_Blosum62','CDR3dist_Blosum45','CDR3dist_Pam10','CDR3dist_K-Mer_3','CDR3dist_K-Mer_4',
               'CDR3dist_K-Mer_5','CDR3dist_Atchley_Factors','CDR3dist_hydrophobicity','CDR3dist_Immunomap']

for i in range(len(all_folders)):
    print("python3 eval_confus_matrix.py '%s/%s/' '%s/Binary_group' 'ROC_group_matches'"% (filepath, all_folders[i], filepath))            
for j in range(len(all_folders)):
    print("python3 eval_confus_matrix.py '%s/%s/' '%s/Binary_exact' 'ROC_exact_matches'"% (filepath, all_folders[j], filepath))            

print("echo '########################################################'")      
print("echo 'Generation of TPR and ROC curves'")
print("echo '########################################################'")     
print("python3 eval_ROC.py '%s/ROC_exact_matches'" % (filepath))
print("mv *.png ROC_exact_matches")
print("python3 eval_ROC.py '%s/ROC_group_matches'" % (filepath))
print("mv *.png ROC_group_matches")

print("echo '########################################################'")      
print("echo 'The distance matrices and curves are ready'")
print("echo '########################################################'")     
      
print("echo '########################################################'")      
print("echo 'Starting logarithmic regression calculations to find coeficients'")
print("echo '########################################################'")     
 
print("echo 'CDR1_group'")
print("echo '########################################################'")  
print("python3 eval_log_reg_CDR12.py '%s' '%s' '%s/Binary_group' %s 50 %s 1 "% (filepath, filename, filepath, length, rest))        
print("echo 'CDR1_Exact'")
print("echo '########################################################'")  
print("python3 eval_log_reg_CDR12.py '%s' '%s' '%s/Binary_exact' %s 50 %s 1 "% (filepath, filename, filepath, length, rest))        
print("echo 'CDR2_group'")
print("echo '########################################################'")  
print("python3 eval_log_reg_CDR12.py '%s' '%s' '%s/Binary_group' %s 50 %s 2 "% (filepath, filename, filepath, length, rest))        
print("echo 'CDR2_exact'")
print("echo '########################################################'")  
print("python3 eval_log_reg_CDR12.py '%s' '%s' '%s/Binary_exact' %s 50 %s 2 "% (filepath, filename, filepath, length, rest))        
print("echo 'CDR3_group'")
print("echo '########################################################'")  
print("python3 eval_log_reg_CDR3.py '%s' '%s' '%s/Binary_group' %s 50 %s "% (filepath, filename, filepath, length, rest))        
print("echo 'CDR3_exact'")
print("echo '########################################################'")  
print("python3 eval_log_reg_CDR3.py '%s' '%s' '%s/Binary_exact' %s 50 %s "% (filepath, filename, filepath, length, rest))        
"""