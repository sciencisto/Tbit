
NB! the ./Running_calc_and_eval.py has not been updated. So it does not run 
as whished at the moment. Further, i have seperated out the calculations and 
evaluation in two folders, which is not in compliance with the code. 


Distance Metric for T-cell receptor beta subunit (TCRß) based on 
3D conformational predictions and sequence comparisons

###############################################################################
M.Sc. Thesis, Master of Science in Engineering, Technical University of Denmark
Kristine Degn, 2020
###############################################################################

This program computes the pairwise distance between all input TCRs. The 
program is run by downloading all the scripts and libraries into one folder. 

The program is executed in terminal:


./Running_calc_and_eval.py > calculate_and_evaluate.sh
/work1/s111518/Train
train_VDJ.tsv

Now you write the filepath to the folder without ' or ". Eg. /work1/s111518
without the seemingly last / press enter and write the input file name, 
agian without ' or ". 

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
    

in terminal write ./calculate_and_evaluate.sh (remeber to make it executable)


Now the calculations are running and you will be able to acces the following 
output: 
----------------------------------
1. ROC_exact_matches: here is the TPR and ROC curves. 
2. distribution_plots: distribution of the distance measure
3. distance_matrice: the distance matrices calculated. 


 
