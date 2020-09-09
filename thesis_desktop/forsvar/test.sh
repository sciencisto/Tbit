Enter path directoryEnter filenameecho '########################################################'
echo 'Starting the BetaDist, Kristines Thesis 2020'
echo '########################################################'
mkdir BetaDistance
python3 calc_TCR_distance_model.py '/home/kristine/Desktop/BetaDist/' 'defence_set.tsv' 50 0
mv *dmat.tsv BetaDistance
python3 eval_distribution_plot.py 'BetaDistance' 'BetaDistance'
mkdir Binary
python3 eval_TCR_distance_binary.py '/home/kristine/Desktop/BetaDist/' 'defence_set.tsv' 50 0 'group'
mv *dmat.tsv Binary
python3 subsample_validate.py '/home/kristine/Desktop/BetaDist/BetaDistance' 'defence_set.tsv' '/home/kristine/Desktop/BetaDist/Binary' 0 50 5
python3 ROC.py '/home/kristine/Desktop/BetaDist/BetaDist_subsampled.csv'
chmod 777 *
