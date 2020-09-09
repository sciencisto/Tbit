Enter path directoryEnter filenameecho '########################################################'
echo 'Starting the BetaDist, Kristines Thesis 2020'
echo '########################################################'
mkdir BetaDistance
python3 calc_TCR_distance_model.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 0
mv *dmat.tsv BetaDistance
python3 calc_TCR_distance_model.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 1
mv *dmat.tsv BetaDistance
python3 calc_TCR_distance_model.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 2
mv *dmat.tsv BetaDistance
python3 calc_TCR_distance_model.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 3
mv *dmat.tsv BetaDistance
python3 calc_TCR_distance_model.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 4
mv *dmat.tsv BetaDistance
python3 calc_TCR_distance_model.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 5
mv *dmat.tsv BetaDistance
python3 calc_TCR_distance_model.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 6
mv *dmat.tsv BetaDistance
python3 calc_TCR_distance_model.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 7
mv *dmat.tsv BetaDistance
python3 calc_TCR_distance_model.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 8
mv *dmat.tsv BetaDistance
python3 calc_TCR_distance_model.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 9
mv *dmat.tsv BetaDistance
python3 calc_TCR_distance_model.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 10
mv *dmat.tsv BetaDistance
python3 calc_TCR_distance_model.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 11
mv *dmat.tsv BetaDistance
python3 calc_TCR_distance_model.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 12
mv *dmat.tsv BetaDistance
python3 calc_TCR_distance_model.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 13
mv *dmat.tsv BetaDistance
python3 calc_TCR_distance_model.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 14
mv *dmat.tsv BetaDistance
python3 calc_TCR_distance_model.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 15
mv *dmat.tsv BetaDistance
python3 calc_TCR_distance_model.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 16
mv *dmat.tsv BetaDistance
python3 eval_distribution_plot.py 'BetaDistance' 'BetaDistance'
mkdir Binary
python3 eval_TCR_distance_binary.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 0 'group'
mv *dmat.tsv Binary
python3 eval_TCR_distance_binary.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 1 'group'
mv *dmat.tsv Binary
python3 eval_TCR_distance_binary.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 2 'group'
mv *dmat.tsv Binary
python3 eval_TCR_distance_binary.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 3 'group'
mv *dmat.tsv Binary
python3 eval_TCR_distance_binary.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 4 'group'
mv *dmat.tsv Binary
python3 eval_TCR_distance_binary.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 5 'group'
mv *dmat.tsv Binary
python3 eval_TCR_distance_binary.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 6 'group'
mv *dmat.tsv Binary
python3 eval_TCR_distance_binary.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 7 'group'
mv *dmat.tsv Binary
python3 eval_TCR_distance_binary.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 8 'group'
mv *dmat.tsv Binary
python3 eval_TCR_distance_binary.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 9 'group'
mv *dmat.tsv Binary
python3 eval_TCR_distance_binary.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 10 'group'
mv *dmat.tsv Binary
python3 eval_TCR_distance_binary.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 11 'group'
mv *dmat.tsv Binary
python3 eval_TCR_distance_binary.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 12 'group'
mv *dmat.tsv Binary
python3 eval_TCR_distance_binary.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 13 'group'
mv *dmat.tsv Binary
python3 eval_TCR_distance_binary.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 14 'group'
mv *dmat.tsv Binary
python3 eval_TCR_distance_binary.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 15 'group'
mv *dmat.tsv Binary
python3 eval_TCR_distance_binary.py '/home/kristine/Desktop/BetaDist/' 'example_TCR.tsv' 50 16 'group'
mv *dmat.tsv Binary
python3 subsample_validate.py '/home/kristine/Desktop/BetaDist/BetaDistance' 'example_TCR.tsv' '/home/kristine/Desktop/BetaDist/Binary' 0 50 30
python3 subsample_validate.py '/home/kristine/Desktop/BetaDist/BetaDistance' 'example_TCR.tsv' '/home/kristine/Desktop/BetaDist/Binary' 1 50 30
python3 subsample_validate.py '/home/kristine/Desktop/BetaDist/BetaDistance' 'example_TCR.tsv' '/home/kristine/Desktop/BetaDist/Binary' 2 50 30
python3 subsample_validate.py '/home/kristine/Desktop/BetaDist/BetaDistance' 'example_TCR.tsv' '/home/kristine/Desktop/BetaDist/Binary' 3 50 30
python3 subsample_validate.py '/home/kristine/Desktop/BetaDist/BetaDistance' 'example_TCR.tsv' '/home/kristine/Desktop/BetaDist/Binary' 4 50 30
python3 subsample_validate.py '/home/kristine/Desktop/BetaDist/BetaDistance' 'example_TCR.tsv' '/home/kristine/Desktop/BetaDist/Binary' 5 50 30
python3 subsample_validate.py '/home/kristine/Desktop/BetaDist/BetaDistance' 'example_TCR.tsv' '/home/kristine/Desktop/BetaDist/Binary' 6 50 30
python3 subsample_validate.py '/home/kristine/Desktop/BetaDist/BetaDistance' 'example_TCR.tsv' '/home/kristine/Desktop/BetaDist/Binary' 7 50 30
python3 subsample_validate.py '/home/kristine/Desktop/BetaDist/BetaDistance' 'example_TCR.tsv' '/home/kristine/Desktop/BetaDist/Binary' 8 50 30
python3 subsample_validate.py '/home/kristine/Desktop/BetaDist/BetaDistance' 'example_TCR.tsv' '/home/kristine/Desktop/BetaDist/Binary' 9 50 30
python3 subsample_validate.py '/home/kristine/Desktop/BetaDist/BetaDistance' 'example_TCR.tsv' '/home/kristine/Desktop/BetaDist/Binary' 10 50 30
python3 subsample_validate.py '/home/kristine/Desktop/BetaDist/BetaDistance' 'example_TCR.tsv' '/home/kristine/Desktop/BetaDist/Binary' 11 50 30
python3 subsample_validate.py '/home/kristine/Desktop/BetaDist/BetaDistance' 'example_TCR.tsv' '/home/kristine/Desktop/BetaDist/Binary' 12 50 30
python3 subsample_validate.py '/home/kristine/Desktop/BetaDist/BetaDistance' 'example_TCR.tsv' '/home/kristine/Desktop/BetaDist/Binary' 13 50 30
python3 subsample_validate.py '/home/kristine/Desktop/BetaDist/BetaDistance' 'example_TCR.tsv' '/home/kristine/Desktop/BetaDist/Binary' 14 50 30
python3 subsample_validate.py '/home/kristine/Desktop/BetaDist/BetaDistance' 'example_TCR.tsv' '/home/kristine/Desktop/BetaDist/Binary' 15 50 30
python3 subsample_validate.py '/home/kristine/Desktop/BetaDist/BetaDistance' 'example_TCR.tsv' '/home/kristine/Desktop/BetaDist/Binary' 16 50 30
python3 ROC.py '/home/kristine/Desktop/BetaDist/BetaDist_subsampled.csv'
chmod 777 *
