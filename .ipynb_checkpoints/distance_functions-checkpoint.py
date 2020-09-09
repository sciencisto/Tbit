def levenshteinDistance(s1, s2):
	if len(s1) > len(s2):
		s1, s2 = s2, s1
	#
	distances = range(len(s1) + 1)
	for i2, c2 in enumerate(s2):
		distances_ = [i2+1]
		for i1, c1 in enumerate(s1):
			if c1 == c2:
				distances_.append(distances[i1])
			else:
				distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
		distances = distances_
	return distances[-1]


def make_AA_to_Atchley_dict():
	"""
	REFERENCE: 
	Atchley et al. "Solving the protein sequence metric problem", 2005, PNAS, vol. 102, no. 18, pp: 6395-6400
	
	Atchley_factor_data order of the rows follows: 
	list_of_AAs=['A','C', 'D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','Z']
	each element of the list follows the format: ['AF1 AF2 AF3 AF4 AF5', ...]

	# AA_to_Atchley dict  has (letter symbol of AminoAcid  as sting, AF index as integer) as keys, and the corresponding AF as   
    float as values
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
	#
	list_of_AAs=['A','C', 'D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','Z']
	#
	AA_to_Atchley = dict()
	for (AA, row) in zip(list_of_AAs, Atchley_factor_data): 
		for (ii, entry) in enumerate(row.split(" "), start =1):
			AA_to_Atchley[AA, ii] = float(entry)
	#
	return AA_to_Atchley


def Atchley_euclidean_dist(s1,s2, AA_to_Atch = make_AA_to_Atchley_dict(), AF_list = [1,2,3,4,5]):
    """
    Returns the distance calculated as euclidean distance of average AFs for each CDR3
    E.g. each AF is calculated as average AF value for the CDR3 resulting in a 5-tuple 
    corresponding to each AF. For 2 CDR3s the distance between them is calculated as 
    euleadian distance between the two 5-tuples.
    """
    import numpy as np
    s1_list = np.array([sum(AA_to_Atch[AA, AF] for AA in list(s1))/float(len(s1)) for AF in AF_list])
    s2_list = np.array([sum(AA_to_Atch[AA, AF] for AA in list(s2))/float(len(s2)) for AF in AF_list])
    distance = np.sqrt(np.sum((np.array(s1_list)-np.array(s2_list))**2))
    del s1_list
    del s2_list
    return distance


def Hamming_pairwise_dist(s1, s2):
    """This function return a hamming distance between two aminoacid sequences, 
    scaled based on length"""
    
    distance = len(list(filter(lambda x : ord(x[0])^ord(x[1]), 
                               zip(s1, s2))))
    #scale by dividing by the maximum CDR3 sequence length
    distance = distance/(max([len(s1), len(s2)]))
    return distance


def BLOSUM45_score_dist(s1,s2):
    from Bio.Align import substitution_matrices
    from Bio import Align
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -10
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM45")
    aligner.mode = "local"
    score_s12 = aligner.score(s1,s2)
    score11 = aligner.score(s1,s1)
    score22 = aligner.score(s2,s2)
    distance = 1- score_s12/max(score11,score22)
    return distance

def BLOSUM62_score_dist(s1,s2):
    from Bio.Align import substitution_matrices
    from Bio import Align
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -10
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.mode = "local"
    score_s12 = aligner.score(s1,s2)
    score11 = aligner.score(s1,s1)
    score22 = aligner.score(s2,s2)
    distance = 1- score_s12/max(score11,score22)
    return distance


def PAM30_score_dist(s1,s2):
    from Bio.Align import substitution_matrices
    from Bio import Align
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -30
    aligner.substitution_matrix = substitution_matrices.load("PAM30")
    aligner.mode = "local"
    score_s12 = aligner.score(s1,s2)
    score11 = aligner.score(s1,s1)
    score22 = aligner.score(s2,s2)
    distance = 1- score_s12/max(score11,score22)
    return distance

def kmers_dist(s1, s2, k):
    """Returns a jaccard dissimilarity score based on reuccuring k-mers, 
    0 = completely identical and 1 = completly dissimilar"""
    #k-mers in CDR3 for TCR1 
    #################################
    s1_kmers = [s1[ii:k+ii] for ii in range(len(s1))]
    s1_kmers = [kmer for kmer in s1_kmers if len(kmer) == k]
    s1_kmers = list(set(s1_kmers))
    
    #k-mers in CDR3 for TCR2
    #################################  
    #For TCR2 the same logic is used:
    s2_kmers = [s2[ii:k+ii] for ii in range(len(s2))]
    s2_kmers = [kmer for kmer in s2_kmers if len(kmer) == k]
    s2_kmers = list(set(s2_kmers))
    #calculating the elements for a jaccard similarity
    intersection = len(set(s1_kmers).intersection(set(s2_kmers))) 
    union = len(set(s1_kmers).union(set(s2_kmers))) 
    #Calulating the distance (jaccard dissimilarity)
    distance = 1-(intersection / union)
    return distance


def Kyte_Doolittle_hydrophobicity_dist(s1, s2):
    """
    returns a distance based on the hydrophobicity index 
    calculated based on Kyte & Doolittle index of hydrophobicity 
    J. Mol. Biol. 157:105-132(1982). positivity mean hydrophobic
    The values are scaled by the maximum possible difference. The 
    most hydrophobic AA is Isoleucine with a score of 4.5, and the 
    most hydrophilic AA is Arginine R with a score of -4.5. The 
    maximum distance is therefore 9. 
    """
    import numpy as np
    #
    KD_dict = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
          "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
          "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
          "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2} 
    max_dist = 9.0
    score_s1 = np.average([KD_dict[aa] for aa in list(s1)])
    score_s2 = np.average([KD_dict[aa] for aa in list(s2)])
    if score_s2 > score_s1: 
        distance = (score_s2 - score_s1)/max_dist
    else: 
        distance = (score_s1 - score_s2)/max_dist    
    return distance



def CDR12_dist(Vgene1, Vgene2, distance_measure, organism):
    """
    NOTE has to be run outside the libraries
    Input: 
    - V gene sequence 1 
    - CDR3 sequence 2
    - V gene sequence 2 
    - distance measure: "TM": 1- TM score of canonical structures
                        "RMSD": RMSD between canonical structures
                        "Levenstein": Levenstein distance
                        "Hamming": Hamming distance
                        "AF": distance calculated based on the 5 AF distance 
                        "BLOSUM45": BLOSUM45 alignment score scaled by the max value.
                        "kmer_N": k mer distance of length of k-mer = N
    - organism:         "human" or "murine"
    There are two main types of distance matrices:
        1. The structural based libraies: 
            The files are absolute files, that was generated using the LYRA 
            framework to pairwisely compare mainchain conformations of 
            all V-genes in their predicted conformations. Comparison frameworks 
            can be altered, however 1-TM and normalized RMSD score files exist
            today.  

    Output: - CDR1 distance value 
            - CDR2 distance value
    """
    import pandas as pd 
    import pickle
    if distance_measure == "TM": 
        if organism == "human":
            CDR1_TM = pd.read_pickle("./libraries/CDR1_TM_human.pickle")
            CDR2_TM = pd.read_pickle("./libraries/CDR2_TM_human.pickle")
        elif organism == "murine":
            CDR1_TM = pd.read_pickle("./libraries/CDR1_TM_murine.pickle")
            CDR2_TM = pd.read_pickle("./libraries/CDR2_TM_murine.pickle")
        distance_CDR1 = CDR1_TM.loc[Vgene1, Vgene2]
        distance_CDR2 = CDR2_TM.loc[Vgene1, Vgene2]
    elif distance_measure == "RMSD":
        if organism == "human":
            CDR1_RMSD = pd.read_pickle("./libraries/CDR1_RMSD_human.pickle")
            CDR2_RMSD = pd.read_pickle("./libraries/CDR2_RMSD_human.pickle")
        elif organism == "murine":
            CDR1_RMSD = pd.read_pickle("./libraries/CDR1_RMSD_murine.pickle")
            CDR2_RMSD = pd.read_pickle("./libraries/CDR2_RMSD_murine.pickle")
        distance_CDR1 = CDR1_RMSD.loc[Vgene1, Vgene2]
        distance_CDR2 = CDR2_RMSD.loc[Vgene1, Vgene2]
    elif distance_measure in ["Levenstein", "Hamming", "AF", "BLOSUM45", "BLOSUM62", "PAM30", "Kyte-Doolittle"]:
        if organism == "human":
            Vgene_to_CDR_dict =  pickle.load(open("./libraries/Vgene_CDR1_CDR2_seq_human.pickle", "rb"))
        elif organism == "murine":
            Vgene_to_CDR_dict =  pickle.load(open("./libraries/Vgene_CDR1_CDR2_seq_murine.pickle", "rb"))
        CDR1_seq_1, CDR2_seq_1 = Vgene_to_CDR_dict[Vgene1]
        CDR1_seq_2, CDR2_seq_2 = Vgene_to_CDR_dict[Vgene2]
        if distance_measure == "Levenstein":
            distance_CDR1 = levenshteinDistance(CDR1_seq_1, CDR1_seq_2)
            distance_CDR2 = levenshteinDistance(CDR2_seq_1, CDR2_seq_2)
        elif distance_measure == "Hamming":
            distance_CDR1 = Hamming_pairwise_dist(CDR1_seq_1, CDR1_seq_2)
            distance_CDR2 = Hamming_pairwise_dist(CDR2_seq_1, CDR2_seq_2)
        elif distance_measure == "AF":
            distance_CDR1 = Atchley_euclidean_dist(CDR1_seq_1, CDR1_seq_2)
            distance_CDR2 = Atchley_euclidean_dist(CDR2_seq_1, CDR2_seq_2)
        elif distance_measure == "BLOSUM45":
            distance_CDR1 = BLOSUM45_score_dist(CDR1_seq_1, CDR1_seq_2)
            distance_CDR2 = BLOSUM45_score_dist(CDR2_seq_1, CDR2_seq_2)
        elif distance_measure == "BLOSUM62":
            distance_CDR1 = BLOSUM62_score_dist(CDR1_seq_1, CDR1_seq_2)
            distance_CDR2 = BLOSUM62_score_dist(CDR2_seq_1, CDR2_seq_2)
        elif distance_measure == "PAM30":
            distance_CDR1 = PAM30_score_dist(CDR1_seq_1, CDR1_seq_2)
            distance_CDR2 = PAM30_score_dist(CDR2_seq_1, CDR2_seq_2)
        elif distance_measure == "Kyte-Doolittle":
            distance_CDR1 = Kyte_Doolittle_hydrophobicity_dist(CDR1_seq_1, CDR1_seq_2)
            distance_CDR2 = Kyte_Doolittle_hydrophobicity_dist(CDR2_seq_1, CDR2_seq_2)
    elif distance_measure.startswith("kmer_"):
        if organism == "human":
            Vgene_to_CDR_dict =  pickle.load(open("./libraries/Vgene_CDR1_CDR2_seq_human.pickle", "rb"))
        elif organism == "murine":
            Vgene_to_CDR_dict =  pickle.load(open("./libraries/Vgene_CDR1_CDR2_seq_murine.pickle", "rb"))
        CDR1_seq_1, CDR2_seq_1 = Vgene_to_CDR_dict[Vgene1]
        CDR1_seq_2, CDR2_seq_2 = Vgene_to_CDR_dict[Vgene2]
        k = int(distance_measure.split("_")[-1])
        distance_CDR1 = kmers_dist(CDR1_seq_1, CDR1_seq_2, k)
        distance_CDR2 = kmers_dist(CDR2_seq_1, CDR2_seq_2, k)
    else: 
        print("wrong distance measure imputed. Please specify again.")
    return distance_CDR1, distance_CDR2




def CDR3_dist(CDR3_1, CDR3_2, distance_measure):
    """
    Input: 
    - CDR3 sequence 1
    - CDR3 sequence 2
    - distance measure: "Levenstein": Levenstein distance
                        "Hamming": Hamming distance
                        "AF": distance calculated based on the 5 AF distance 
                        "BLOSUM45": BLOSUM45 alignment score scaled by the max value.
                        "kmer_N": k mer distance of length of k-mer = N
    Output: - CDR3 distance value 
    """
    if distance_measure in ["Levenstein", "Hamming", "AF", "BLOSUM45", "BLOSUM62", "PAM30", "Kyte-Doolittle"]: 
        if distance_measure == "Levenstein":
            distance_CDR3 = levenshteinDistance(CDR3_1, CDR3_2)
        elif distance_measure == "Hamming":
            distance_CDR3 = Hamming_pairwise_dist(CDR3_1, CDR3_2)
        elif distance_measure == "AF":
            distance_CDR3 = Atchley_euclidean_dist(CDR3_1, CDR3_2)
        elif distance_measure == "BLOSUM45":
            distance_CDR3 = BLOSUM45_score_dist(CDR3_1, CDR3_2)
        elif distance_measure == "BLOSUM62":
            distance_CDR3 = BLOSUM62_score_dist(CDR3_1, CDR3_2)
        elif distance_measure == "PAM30":
            distance_CDR3 = PAM30_score_dist(CDR3_1, CDR3_2)
        elif distance_measure == "Kyte-Doolittle":
            distance_CDR3 = Kyte_Doolittle_hydrophobicity_dist(CDR3_1, CDR3_2)
    elif distance_measure.startswith("kmer_"):
        k = int(distance_measure.split("_")[-1])
        distance_CDR3 = kmers_dist(CDR3_1, CDR3_2, k)
    else: 
        print("wrong distance measure imputed. Please specify again.")
    return distance_CDR3



def Tbit_dist(Vgene_1, CDR3_1, Vgene_2, CDR3_2, model):
    """
    NB: Has to be run outside the libraries directory
    Input: 
            Vgene_1 - name of the V gene of the 1st TCR
            CDR3_1  - AA sequence of CDR3 region of 1st TCR
            Vgene_2 - name of the V gene of the 2nd TCR
            CDR3_2  - AA sequence of CDR3 region of the 2nd TCR
            model_coeff - which model should be used (after testing redundant)
    Output: 
            distance value calculated with the model 
    """
    import pandas as pd 
    import pickle
    from collections import defaultdict
    # import model coefficients: 
    if model  == "Logistic regression":
        model_coeff =  pickle.load(open("./Tbit_log_regression_coeff.pickle", "rb"))
    distance_measures = list(set([method.split(" ")[-1] for method in model_coeff.keys()]))
    Tbit_distance_dict = defaultdict()
    for distance_measure in distance_measures:
        if distance_measure in ["TM", "RMSD"]:
            dist_CDR1, dist_CDR2 = CDR12_dist(Vgene_1, Vgene_2, distance_measure)
            Tbit_distance_dict["CDR1 " + distance_measure] = dist_CDR1*model_coeff["CDR1 " + distance_measure]
            Tbit_distance_dict["CDR2 " + distance_measure] = dist_CDR2*model_coeff["CDR2 " + distance_measure]
        else: 
            dist_CDR1, dist_CDR2 = CDR12_dist(Vgene_1, Vgene_2, distance_measure)
            dist_CDR3 = CDR3_dist(CDR3_1, CDR3_2, distance_measure)
            Tbit_distance_dict["CDR1 " + distance_measure] = dist_CDR1*model_coeff["CDR1 " + distance_measure]
            Tbit_distance_dict["CDR2 " + distance_measure] = dist_CDR2*model_coeff["CDR2 " + distance_measure]
            Tbit_distance_dict["CDR3 " + distance_measure] = dist_CDR3*model_coeff["CDR3 " + distance_measure]
    Tbit_dist = sum(Tbit_distance_dict.values())
    return Tbit_dist


#def TCRdist(Vgene_1, CDR3_1): 
#    import pickle
#    CDR3_1 = CDR3_1[2:]
#    Vgene_to_CDR_dict =  pickle.load(open("./libraries/Vgene_CDR1_CDR2_seq_murine.pickle", "rb"))
#    Vgene_to_CDR25_dict =  pickle.load(open("./libraries/Vgene_CDR2_5_seq_murine.pickle", "rb"))


