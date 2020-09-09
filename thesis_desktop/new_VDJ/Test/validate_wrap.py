"""
./BetaDistanceStart.py > run_calc.sh
"""
    
filepath = input("Enter path to the directory") 
filename = input("Enter filename")

print("echo '########################################################'")      
print("echo 'Starting the BetaDist, Kristines Thesis 2020'")
print("echo '########################################################'")      

#Import packages
import pandas as pd
import math

#Impart the file and calculate the number of chunks needed
TSV = pd.read_csv(filepath+"/"+filename, sep="\t", header=None)
length = math.ceil(len(TSV)/50)
#count how many TCRs are present in the last chunk
remain = len(TSV)-((length-1)*50)
del TSV
      
print("mkdir BetaDistance")
for i in range(length):
    print("python3 calc_TCR_distance_model.py '%s/' '%s' 50 %s" % (filepath, 
                                                                   filename, i))
    print("mv *dmat.tsv BetaDistance")


##############################################################################
# This part is only if evaluation of the distance metric is desired
##############################################################################

#Distribution of the distance
print("python3 eval_distribution_plot.py 'BetaDistance' 'BetaDistance'")

#Generating a Binary file, notice the option "group" require a library file. 
#Alternatively choose option "exact"
print("mkdir Binary")
for i in range(length):
    print("python3 eval_TCR_distance_binary.py '%s/' '%s' 50 %s 'group'" 
          % (filepath, filename, i))
    print("mv *dmat.tsv Binary")

#Subsampling the file outouts one csv file with two columns, distance and 
#binary information
for i in range(length):
    print("python3 subsample_validate.py '%s' '%s' '%s/Binary' %s 50 %s" % 
          (filepath, filename, filepath, i, remain))

print("chmod 777 *")