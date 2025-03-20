#!/bin/bash

# run it in easal_output/

###DSSO
#simulated
directories=("1clv_cl5" "1dfj_cl10" "1dfj_cl11" "1kxp_cl10" "1r0r_cl5" "2ayo_cl10" "2ayo_cl15" "2b42_cl10" "2hle_cl10" "2hle_cl15")

for dir in "${directories[@]}"; do
    cd "$dir"

    filename=$(echo "$dir" | cut -d'_' -f1)
    number=$(echo "$dir" | sed 's/.*_cl\([0-9]*\)/\1/')

    #chainA and chainB based on the directory name
    if [[ "$filename" == 1clv* ]]; then
        chainA="A"
        chainB="I"
    elif [[ "$filename" == 1dfj* || "$filename" == 1r0r* ]]; then
        chainA="E"
        chainB="I"
    elif [[ "$filename" == 1kxp* ]]; then
        chainA="A"
        chainB="D"
    else
        chainA="A"
        chainB="B"
    fi

    csv_file="~/Integrative_docking_benchmark/benchmark/simulated/DSSO/test_FP/${filename}_DSSO_${number}.csv"
    python ~/Integrative_docking_benchmark/scripts/compare_ensembles/crosslink_distance_perc_satisfied_calculation/easal_output/calc_xlink_dist_perc_pdb.py "$chainA" "$chainB" "$csv_file" "og"

    cd ../
done
