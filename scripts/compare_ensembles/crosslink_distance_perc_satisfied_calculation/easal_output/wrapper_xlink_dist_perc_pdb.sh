#!/bin/bash

# run it in easal_output/

###DSSO
#simulated
directories=("1clv_cl2" "1clv_cl6" "1dfj_cl3" "1dfj_cl9" "1dfj_cl12" "1kxp_cl4" "1kxp_cl7" "1kxp_cl11" "1r0r_cl3" "1r0r_cl7" "2ayo_cl4" "2ayo_cl8" "2ayo_cl13" "2b42_cl5" "2b42_cl10" "2hle_cl5" "2hle_cl10" "2hle_cl14")

for dir in "${directories[@]}"; do
    cd "DSSO_pdbs/$dir"

    filename=$(echo "$dir" | cut -d'_' -f1)
    number=$(echo "$dir" | sed 's/.*_cl\([0-9]*\)/\1/')


    #directory based on the directory name
    if (( number <= 5 )); then
        direc2="less_than_5"
    elif (( number >= 6 && number <= 10 )); then
        direc2="6-10"
    else
        direc2="more_than_10"
    fi

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

    csv_file="$HOME/Integrative_docking_benchmark/benchmark/simulated/DSSO/$direc2/${filename}_DSSO_${number}.csv"
    python ~/Integrative_docking_benchmark/scripts/compare_ensembles/crosslink_distance_perc_satisfied_calculation/easal_output/calc_xlink_dist_perc_pdb.py "$chainA" "$chainB" "$csv_file" "og"

    cd ../../..
done

# #experimental
directories=("gata_gatc_cl3" "gcvpa_gcvpb_cl5" "phes_phet_cl8" "roca_putc_cl2" "sucd_succ_cl4")

for dir in "${directories[@]}"; do
    cd "DSSO_pdbs/$dir"

    filename=$(echo "$dir" | cut -d'_' -f1)
    filename+="_"
    filename+=$(echo "$dir" | cut -d'_' -f2)
    number=$(echo "$dir" | sed 's/.*_cl\([0-9]*\)/\1/')

    csv_file="$HOME/Integrative_docking_benchmark/benchmark/experimental/crosslinks/${filename}_DSSO_${number}.csv"
    python ~/Integrative_docking_benchmark/scripts/compare_ensembles/crosslink_distance_perc_satisfied_calculation/easal_output/calc_xlink_dist_perc_pdb.py "A" "B" "$csv_file" "og"

    cd ../../..
done

###DMTMM
directories=("1clv_cl8" "1dfj_cl4" "1kxp_cl7" "1r0r_cl6" "2ayo_cl5" "2b42_cl10" "2hle_cl9")

for dir in "${directories[@]}"; do
    cd "DMTMM_pdbs/$dir"

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

    csv_file="$HOME/Integrative_docking_benchmark/benchmark/simulated/DMTMM/${filename}_DMTMM_${number}.csv"
    python ~/Integrative_docking_benchmark/scripts/compare_ensembles/crosslink_distance_perc_satisfied_calculation/easal_output/calc_xlink_dist_pdb.py "$chainA" "$chainB" "$csv_file" "og"

    cd ../../
done
