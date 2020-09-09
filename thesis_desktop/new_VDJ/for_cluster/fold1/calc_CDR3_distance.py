"""
THIS SCRIPT IS USED TO COMPARE TWO CDR3-SEQUENCES WHEN COMPARING TCRs
---------------------------------------------------------------------------
Function content of this script:
    
    Alignments
    ---------------------------------------------
    1. Hamming distance
    2. Alignment comparison diagonal dictionary
    3. Pairwise alignment tools
    4. k-mer comparisons
    5. Atchley Factors
    6. Hydrophobicity
    
    Then there is also the Immunomap method.
    
    Timing:
    ----------------------------------------------
    For two CDR3s, TCR1_CDR3 = "CASSYSRTGSYEQYF", TCR2_CDR3 = "CASSVEGPGELFF",
    the calulation time is as follows:
        
                     |  the first | the subsequent | the distance
                     --------------------------------------------
    Hamming distance | 1.26e-05 s |   1.26e-05 s   |    0.32
    Blosum 62        | 0.0016 s*  |   0.0003 s     |    0.68
    Blosum 45        | 0.0007 s*  |   0.0003 s     |    0.62
    Pam10            | 0.0007 s*  |   0.0003 s     |    0.78
    K-mer, k = 3     | 5.65e-05 s |   5.65e-05 s   |    0.91
    K-mer, k = 4     | 9.66e-05 s |   9.66e-05 s   |    0.95
    K-mer, k = 5     | 6.34e-05 s |   6.34e-05 s   |    1.00
    Atchley factor   | 0.0055 s   |   0.0055 s     |    0.16
    Hydrophobicity   | 0.0009 s   |   0.0009 s     |    0.55
    
    *Depend on the length of the total list of TCRs to be compared. Here only 
    two in the list.


Input:
-------------------------------------------------------------------------------    
TCR1_CDR3, TCR2_CDR3 strings of Amino Acids as well as A TCR list, all CDR3s 
to be pairwise calculated, and a k-value. All inputs are administered by he
calc_TCR_distance.py 

Dependencies:
------------------------------
the script rely on the Biopython module, pandas and the 
calc_substitution_matrices  script. 

Output
-------------------------------------------------------------------------------
A pairwise distance, one float.

"""
import numpy as np
from Bio import pairwise2
##############################################################################
#Sequence comparison: Alignment
##############################################################################

#Hamming
def dist_CDR3_pairwise_ham(TCR1_CDR3, TCR2_CDR3):
    """This function return a distance between two CDR3 sequences based 
    on their hamming distance, and scaled based on length"""
    
    distance = len(list(filter(lambda x : ord(x[0])^ord(x[1]), 
                               zip(TCR1_CDR3, TCR2_CDR3))))
    
    #scale by dividing by the maximum CDR3 sequence length
    distance = distance/(25)
    
    return distance

#Substitution matrix based   
def CDR3_diagonal_dict(TCR_list, metric, gappenalty, gap_extend):
    """This function return a global substitution matrix score on a CDR3 
    sequence towards itself. It is used in dist_CDR3_alignment to scale to 1-0    
    The input is a list of CDR3 sequences and the output is a dictionary
    with each sequence as key at the score as value
    
        TCR_list: a list of CDR3 sequences 
        metric: The chosen substiturion matrix from calc_substitution_matrices
        gappenalty: Chosen gap-penalty, needs to be a negative number        
    """
    
    #generate an empty list
    diagonal = []
    
    # loop over each CDR3 sequence in the TCR_list
    for ii in range(len(TCR_list)):
        #for each CDR3 sequence align globally to itself
        for a in pairwise2.align.globalds(TCR_list[ii], TCR_list[ii], metric, 
                                          gappenalty, gap_extend):
            al1,al2, score, begin, end = a
            #score to itself is found
            score = a[2]
        #each score is appended to the empty list
        diagonal.append(score)
    #the scores and the CDR3 sequences is converted to a dictionary with 
    #CDR3 sequence as keys and their score as values
    zipbObj = zip(TCR_list, diagonal)
    diagonal_dict = dict(zipbObj)

    return diagonal_dict

def dist_CDR3_alignment(TCR1_CDR3, TCR2_CDR3, metric, gappenalty, gap_extend, 
                        diagonal, scale_operator): 
    """This function returns a global alignment distance between two CDR3 
    sequences based on the input parameters:
    
        TCR1_CDR3, TCR2_CDR3: 
            the two CDR3 sequences to compare.  
        metric: 
            The chosen substiturion matrix from substitution_matrices
        gappenalty: 
            Chosen gap-penalty, needs to be a negative number
        diagonal: 
            the diagonal score dictionary calulated for scaling, hence,
            a score of zero means completly identical and a score of 1 is 
            completly different.
        scale_op: the choice of operator to scale, the choices are: 
            operator.ge if the score should be scaled on the maximum diagonal 
            value or operator.le if the score should be scaled on the minimum
            diagonal value.
        
        Example:
            dist_CDR3_alignment('CASSYSRTGSYEQYF', 'CASSVEGPGELFF', Blosum45, 
                                -10, -0.5, Blos45diagonal*, operator.ge)
            
            * Calculated with: 
                Blos45diagonal = CDR3_diagonal_dict(TCR_list, Blosum45, -10, -0.5)
    """ 
    #The two CDR3 sequences are aligned globally
    for a in pairwise2.align.globalds(TCR1_CDR3, TCR2_CDR3, metric, 
                                      gappenalty, gap_extend):
        al1,al2, score, begin, end = a
        dist = a[2] 
    
    #scaled based on the chosen diagonal
        if scale_operator(diagonal[TCR1_CDR3], diagonal[TCR2_CDR3]):
            distance = 1-(dist/diagonal[TCR1_CDR3])
        else:
            distance = 1-(dist/diagonal[TCR2_CDR3])
    
    return distance

##############################################################################
#Sequence comparison: k-mer comparison
##############################################################################

def kmers_dist(TCR1_CDR3, TCR2_CDR3, k):
    """Returns a jaccard dissimilarity score based on reuccuring k-mers, 
    0 = completely identical and 1 = completly dissimilar"""
    
    #k-mers in CDR3 for TCR1 
    #################################
    
    #Empty list
    list_of_TCR1_kmers = []    
    #Start the sliding window to find AA sequences of length k
    for ii in range(len(TCR1_CDR3)):
        kmer = TCR1_CDR3[ii:k+ii]
        #Append the k-mers to the emplty list
        list_of_TCR1_kmers.append(kmer)        
        
        #At the end of the CDR3 sequence, the last entries in the sliding 
        #windown is shorter than k, so they are removed:
        TCR1_long_kmers = [s for s in list_of_TCR1_kmers if len(s) == k]
    
    #The unique kmers for the CDR3 sequencxe is saved in a new list
    TCR1_unique_kmers = set(TCR1_long_kmers)

    #k-mers in CDR3 for TCR2
    #################################
    
    #For TCR2 the same logic is used:
    list_of_TCR2_kmers = []
    for ii in range(len(TCR2_CDR3)):
        kmer = TCR2_CDR3[ii:k+ii]
        list_of_TCR2_kmers.append(kmer)        
        
        TCR2_long_kmers = [s for s in list_of_TCR2_kmers if len(s) == k]
    
    TCR2_unique_kmers = set(TCR2_long_kmers)
    
    #calculating the elements for a jaccard similarity
    intersection = len(TCR1_unique_kmers.intersection(TCR2_unique_kmers)) 
    union = len(TCR1_unique_kmers.union(TCR2_unique_kmers)) 
    
    #Calulating their distance (jaccard dissimilarity)
    distance = 1-(intersection / union)

    return distance

##############################################################################
#Sequence comparison: Atchley factors, Code written by Milena
##############################################################################

    
def make_AA_to_Atchley_dict():
    """
    Atchley_factor_data order of the rows follows: 
    list_of_AAs=['A','C', 'D','E','F','G','H','I','K','L','M','N','P','Q','R',
    'S','T','V','W','Y','Z']
    each element of the list follows the format: ['AF1 AF2 AF3 AF4 AF5', ...]
    # AA_to_Atchley dict  has (letter symbol of AminoAcid  as sting, AF index 
    as integer) as keys, and the corresponding AF as float as values
    """

    Atchley_factor_data = ['-0.591 -1.302 -0.733 1.570 -0.146', 
                           '-1.343 0.465 -0.862 -1.020 -0.255',
                           '1.050 0.302 -3.656 -0.259 -3.242',
                           '1.357 -1.453 1.477 0.113 -0.837',
                           '-1.006 -0.590 1.891 -0.397 0.412',
                           '-0.384 1.652 1.330 1.045 2.064',
                           '0.336 -0.417 -1.673 -1.474 -0.078',
                           '-1.239 -0.547 2.131 0.393 0.816',
                           '1.831 -0.561 0.533 -0.277 1.648',
                           '-1.019 -0.987 -1.505 1.266 -0.912',
                           '-0.663 -1.524 2.219 -1.005 1.212',
                           '0.945 0.828 1.299 -0.169 0.993',
                           '0.189 2.081 -1.628 0.421 -1.392',
                           '0.931 -0.179 -3.005 -0.503 -1.853',
                           '1.538 -0.055 1.502 0.440 2.897',
                           '-0.228 1.339 -4.760 0.670 -2.647',
                           '-0.032 0.326 2.213 0.908 1.313',
                           '-1.337 -0.279 -0.544 1.242 -1.262',
                           '-0.595 0.009 0.672 -2.128 -0.184',
                           '0.260 0.830 3.097 -0.838 1.512',
                           '0.000 0.000 0.000 0.000 0.000']

    list_of_AAs=['A','C', 'D','E','F','G','H','I','K','L','M','N','P','Q','R',
                 'S','T','V','W','Y','Z']

    AA_to_Atchley = dict()
    for (AA, row) in zip(list_of_AAs, Atchley_factor_data): 
        for (ii, entry) in enumerate(row.split(" "), start =1):
            AA_to_Atchley[AA, ii] = float(entry)

    return AA_to_Atchley

def Atchley_euclidean_dist(s1,s2, AA_to_Atch = make_AA_to_Atchley_dict(), 
                           AF_list = [1,2,3,4,5]):
    """
    Returns the distance calculated as euclidean distance of average AFs for 
    each CDR3 E.g. each AF is calculated as average AF value for the CDR3 
    resulting in a 5-tuple corresponding to each AF. For 2 CDR3s the distance 
    between them is calculated as euleadian distance between the two 5-tuples.
    """
    s1_list = np.array([sum(AA_to_Atch[AA, AF] for AA in 
                            list(s1))/float(len(s1)) for AF in AF_list])
    s2_list = np.array([sum(AA_to_Atch[AA, AF] for AA in 
                            list(s2))/float(len(s2)) for AF in AF_list])
    distance = np.sqrt(np.sum((np.array(s1_list)-np.array(s2_list))**2))
    del s1_list
    del s2_list
    
    #scaleing based on emperical distribution evidence for 20.000 CDR3s 
    distance = distance/4.5
    
    return distance

# Kyte & Doolittle index of hydrophobicity 
# J. Mol. Biol. 157:105-132(1982). positivity mean hydrophobic
kd = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
      "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
      "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
      "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2} 

def physiocehmical_prop(TCR1_CDR3, TCR2_CDR3):
    """Function that returns a distance based on the difference in 
    hydrophobicity between two CDR3 sequences"""
    
    #Empty lists
    TCR1_letters = []
    TCR2_letters = []
    
    #A hydrophobicity score is picked up for each letter
    #from the kd dictionary
    for i in range(len(TCR1_CDR3)):
        letter = kd[TCR1_CDR3[i]]
        TCR1_letters.append(letter)
    
    #The average hydrophibicity for the CDR3 sequence is computed
    TCR1_hydrophobicity = np.average(TCR1_letters)
    #* note: avaerage was chosen to keep the score relatively low. Summing 
    #resulted in the same distribution around a higer mean.
    
    #The same is computed for the other CDR3
    for i in range(len(TCR2_CDR3)):
        letter = kd[TCR2_CDR3[i]]
        TCR2_letters.append(letter)
    TCR2_hydrophobicity = np.average(TCR2_letters)
    
    #The smaller hydrophobicity is subtracted from the lager 
    #hydrophobicity. Hence giving a positive number that inform the 
    #difference between the CDR3s, here interpreted as distance.
    
    # Logic: A large difference == a large distance, unlikely to bind the 
    #same epitope.
    
    if TCR1_hydrophobicity > TCR2_hydrophobicity:
        distance = TCR1_hydrophobicity-TCR2_hydrophobicity
    else:
        distance = TCR2_hydrophobicity-TCR1_hydrophobicity
    
    #scaleing based on emperical distribution evidence for 20.000 CDR3s 
    distance = distance/3
    return distance

#Zamyatnin, A.A., Protein volume in solution, Prog. Biophys. Mol. Biol., 24:107-123 (1972), PMID: 4566650.

vol = {"A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5,
      "Q": 143.8, "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7,
      "L": 166.7, "K": 168.6, "M": 162.9, "F": 189.9, "P": 112.7,
      "S": 89.0, "T": 116.1, "W": 227.7, "Y": 193.6, "V": 140.0}

def vol_distance(TCR1_CDR3, TCR2_CDR3):
    """Function that returns a distance based on the difference in 
    volumne in solution between two CDR3 sequences"""
    
    #Empty lists
    TCR1_letters = []
    TCR2_letters = []
    
    #A volume score is picked up for each letter
    #from the Zamyatnin dictionary
    for i in range(len(TCR1_CDR3)):
        letter = vol[TCR1_CDR3[i]]
        TCR1_letters.append(letter)
    
    #The sum of the volume for the CDR3 sequence is computed
    TCR1_volume = np.sum(TCR1_letters)
    
    #The same is computed for the other CDR3
    for i in range(len(TCR2_CDR3)):
        letter = vol[TCR2_CDR3[i]]
        TCR2_letters.append(letter)
    TCR2_volume = np.sum(TCR2_letters)
    
    #The smaller volume is subtracted from the lager 
    #volume. Hence giving a positive number that inform the 
    #difference between the CDR3s, here interpreted as distance.
    
    # Logic: A large difference == a large distance, unlikely to bind the 
    #same epitope.
    
    if TCR1_volume > TCR2_volume:
        distance = TCR1_volume-TCR2_volume
    else:
        distance = TCR2_volume-TCR1_volume
    
    distance = distance/1300
    #scaleing based on emperical distribution evidence for 20.000 CDR3s 
    return distance
##############################################################################
#Competitor methods
##############################################################################
def Immunomap_Dist(TCR1_CDR3, TCR2_CDR3, metric, gappenalty, gap_extend, 
                        diagonal): 
 
    #The two CDR3 sequences are aligned globally
    for a in pairwise2.align.globalds(TCR1_CDR3, TCR2_CDR3, metric, 
                                      gappenalty, gap_extend):
        al1,al2, score, begin, end = a
        dist = a[2] 
    
    #calculating the Immunomapdistance for each pair of sequences. 
    Immunomap_dist = (1-(dist/diagonal[TCR1_CDR3]))*(1-(dist/diagonal[TCR2_CDR3]))
    
    return Immunomap_dist

