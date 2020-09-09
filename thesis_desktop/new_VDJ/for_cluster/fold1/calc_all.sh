Enter path to the kristine_thesis_runthrough directoryEnter filenameecho '########################################################'
echo 'Starting the TCRdistance program, Kristines Thesis 2020'
echo '########################################################'
echo 'Calculating dmats for the TM method for CDR1'
mkdir CDR1dist_TM
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 0 'none' 'TM' 'none'
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 1 'none' 'TM' 'none'
mv *dmat.tsv CDR1dist_TM
echo 'Calculating dmats for the RMSD method for CDR1'
mkdir CDR1dist_RMSD
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 0 'none' 'RMSD' 'none'
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 1 'none' 'RMSD' 'none'
mv *dmat.tsv CDR1dist_RMSD
echo 'Calculating dmats for the Blosum62 method for CDR1'
mkdir CDR1dist_Blosum62
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 0 'none' 'Blosum62' 'none'
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 1 'none' 'Blosum62' 'none'
mv *dmat.tsv CDR1dist_Blosum62
echo 'Calculating dmats for the hydrophobicity method for CDR1'
mkdir CDR1dist_hydrophobicity
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 0 'none' 'hydrophobicity' 'none'
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 1 'none' 'hydrophobicity' 'none'
mv *dmat.tsv CDR1dist_hydrophobicity
echo 'Calculating dmats for the TM method for CDR2'
mkdir CDR2dist_TM
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 0 'none' 'none' 'TM'
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 1 'none' 'none' 'TM'
mv *dmat.tsv CDR2dist_TM
echo 'Calculating dmats for the RMSD method for CDR2'
mkdir CDR2dist_RMSD
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 0 'none' 'none' 'RMSD'
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 1 'none' 'none' 'RMSD'
mv *dmat.tsv CDR2dist_RMSD
echo 'Calculating dmats for the Blosum62 method for CDR2'
mkdir CDR2dist_Blosum62
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 0 'none' 'none' 'Blosum62'
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 1 'none' 'none' 'Blosum62'
mv *dmat.tsv CDR2dist_Blosum62
echo 'Calculating dmats for the hydrophobicity method for CDR2'
mkdir CDR2dist_hydrophobicity
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 0 'none' 'none' 'hydrophobicity'
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 1 'none' 'none' 'hydrophobicity'
mv *dmat.tsv CDR2dist_hydrophobicity
echo 'Calculating dmats for the Hamming method for CDR3'
mkdir CDR3dist_Hamming
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 0 'Hamming' 'none' 'none'
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 1 'Hamming' 'none' 'none'
mv *dmat.tsv CDR3dist_Hamming
echo 'Calculating dmats for the Blosum62 method for CDR3'
mkdir CDR3dist_Blosum62
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 0 'Blosum62' 'none' 'none'
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 1 'Blosum62' 'none' 'none'
mv *dmat.tsv CDR3dist_Blosum62
echo 'Calculating dmats for the Blosum45 method for CDR3'
mkdir CDR3dist_Blosum45
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 0 'Blosum45' 'none' 'none'
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 1 'Blosum45' 'none' 'none'
mv *dmat.tsv CDR3dist_Blosum45
echo 'Calculating dmats for the Pam10 method for CDR3'
mkdir CDR3dist_Pam10
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 0 'Pam10' 'none' 'none'
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 1 'Pam10' 'none' 'none'
mv *dmat.tsv CDR3dist_Pam10
echo 'Calculating dmats for the K-Mer_3 method for CDR3'
mkdir CDR3dist_K-Mer_3
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 0 'K-Mer_3' 'none' 'none'
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 1 'K-Mer_3' 'none' 'none'
mv *dmat.tsv CDR3dist_K-Mer_3
echo 'Calculating dmats for the K-Mer_4 method for CDR3'
mkdir CDR3dist_K-Mer_4
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 0 'K-Mer_4' 'none' 'none'
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 1 'K-Mer_4' 'none' 'none'
mv *dmat.tsv CDR3dist_K-Mer_4
echo 'Calculating dmats for the K-Mer_5 method for CDR3'
mkdir CDR3dist_K-Mer_5
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 0 'K-Mer_5' 'none' 'none'
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 1 'K-Mer_5' 'none' 'none'
mv *dmat.tsv CDR3dist_K-Mer_5
echo 'Calculating dmats for the Atchley_Factors method for CDR3'
mkdir CDR3dist_Atchley_Factors
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 0 'Atchley_Factors' 'none' 'none'
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 1 'Atchley_Factors' 'none' 'none'
mv *dmat.tsv CDR3dist_Atchley_Factors
echo 'Calculating dmats for the hydrophobicity method for CDR3'
mkdir CDR3dist_hydrophobicity
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 0 'hydrophobicity' 'none' 'none'
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 1 'hydrophobicity' 'none' 'none'
mv *dmat.tsv CDR3dist_hydrophobicity
echo 'Calculating dmats for the Immunomap method for CDR3'
mkdir CDR3dist_Immunomap
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 0 'Immunomap' 'none' 'none'
python3 calc_TCR_distance.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 1 'Immunomap' 'none' 'none'
mv *dmat.tsv CDR3dist_Immunomap
echo '########################################################'
echo 'all dmats have been calculated, now the binary files are generated'
echo '########################################################'
echo 'Calculating binary dmat for the exact match'
mkdir Binary_exact
python3 eval_TCR_distance_binary.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 0 'exact'
python3 eval_TCR_distance_binary.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 1 'exact'
mv *dmat.tsv Binary_exact
echo 'Calculating binary dmat for the group match'
mkdir Binary_group
python3 eval_TCR_distance_binary.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 0 'group'
python3 eval_TCR_distance_binary.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/' 'small_set.tsv' 50 1 'group'
mv *dmat.tsv Binary_group
echo '########################################################'
echo 'Making distribution plots for all'
echo '########################################################'
mkdir distribution_plots
python3 eval_distribution_plot.py 'CDR1dist_TM' 'CDR1_TM'
python3 eval_distribution_plot.py 'CDR2dist_TM' 'CDR2_TM'
python3 eval_distribution_plot.py 'CDR1dist_RMSD' 'CDR1_RMSD'
python3 eval_distribution_plot.py 'CDR2dist_RMSD' 'CDR2_RMSD'
python3 eval_distribution_plot.py 'CDR1dist_Blosum62' 'CDR1_Blosum62'
python3 eval_distribution_plot.py 'CDR2dist_Blosum62' 'CDR2_Blosum62'
python3 eval_distribution_plot.py 'CDR1dist_hydrophobicity' 'CDR1_hydrophobicity'
python3 eval_distribution_plot.py 'CDR2dist_hydrophobicity' 'CDR2_hydrophobicity'
python3 eval_distribution_plot.py 'CDR3dist_Hamming' 'CDR3_Hamming'
python3 eval_distribution_plot.py 'CDR3dist_Blosum62' 'CDR3_Blosum62'
python3 eval_distribution_plot.py 'CDR3dist_Blosum45' 'CDR3_Blosum45'
python3 eval_distribution_plot.py 'CDR3dist_Pam10' 'CDR3_Pam10'
python3 eval_distribution_plot.py 'CDR3dist_K-Mer_3' 'CDR3_K-Mer_3'
python3 eval_distribution_plot.py 'CDR3dist_K-Mer_4' 'CDR3_K-Mer_4'
python3 eval_distribution_plot.py 'CDR3dist_K-Mer_5' 'CDR3_K-Mer_5'
python3 eval_distribution_plot.py 'CDR3dist_Atchley_Factors' 'CDR3_Atchley_Factors'
python3 eval_distribution_plot.py 'CDR3dist_hydrophobicity' 'CDR3_hydrophobicity'
python3 eval_distribution_plot.py 'CDR3dist_Immunomap' 'CDR3_Immunomap'
python3 eval_distribution_plot.py 'Binary_exact' 'binary_exact'
python3 eval_distribution_plot.py 'Binary_group' 'binary_group'
mv *.png distribution_plots
echo '########################################################'
echo 'Subsampling'
echo '########################################################'
python3 generate_sub_cdr1.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1' 'small_set.tsv' '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/Binary_group' 2 50 -46
python3 generate_sub_cdr2.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1' 'small_set.tsv' '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/Binary_group' 2 50 -46
python3 generate_sub_cdr3.py '/home/kristine/Desktop/new_VDJ/for_cluster/fold1' 'small_set.tsv' '/home/kristine/Desktop/new_VDJ/for_cluster/fold1/Binary_group' 2 50 -46
mkdir distance_matrices
chmod 777 *
mv CDR1dist_* distance_matrices
mv CDR2dist_* distance_matrices
mv CDR3dist_* distance_matrices
