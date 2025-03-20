#!/bin/bash
###DSSO
##simulated
directories=("1clv_5" "1dfj_10" "1dfj_11" "1kxp_10" "1r0r_5" "2ayo_10" "2ayo_15" "2b42_10" "2hle_10" "2hle_15")

for dir in "${directories[@]}"; do
  filename=$(echo "$dir" | cut -d'_' -f1)
  number=$(echo "$dir" | cut -d'_' -f2)

  csv_file="/home/muskaan/projects/easal_related/Integrative_docking_benchmark/benchmark/simulated/DSSO/test_FP/${filename}_DSSO_${number}.csv"

  ~/imp-clean/build/setup_environment.sh python ~/Integrative_docking_benchmark/scripts/compare_ensembles/crosslink_distance_perc_satisfied_calculation/imp_output/calc_xlink_dist_perc_rmf.py "filtered_rmfs/${dir}_filtered.rmf3" "$csv_file" "32" &

done
