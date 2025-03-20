#!/bin/bash
# ###DSSO
cd DSSO
for protein in 1clv_5 1dfj_11 1dfj_10 1kxp_10 1r0r_5 2ayo_10 2ayo_15 2b42_10 2hle_10 2hle_15; do
    cd "$protein"

    if [[ "$protein" == 1dfj* ||  "$protein" == 1r0r* ]]; then
        chainA="E"
        chainB="I"
    elif [[ "$protein" == 1kxp* ]]; then
        chainA="A"
        chainB="D"
    else
        chainA="A"
        chainB="B"
    fi

    for i in $(seq 20); do
        folder_name="$i"
        mpirun -np 4 ~/imp-clean/build/setup_environment.sh python ~/Integrative_docking_benchmark/scripts/imp/modeling/sample_imp_distance_restraint.py "prod" "$chainA" "$chainB" "$folder_name" "DSSO" "32" 2> "err_$folder_name.log" &
    done

    cd ..
done

wait
