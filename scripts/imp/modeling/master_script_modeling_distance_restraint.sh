#!/bin/bash
# ###DSSO
cd DSSO
for protein in 1clv_2 1clv_6; do
    cd "$protein"

    if [[ "$protein" == 1clv* ]]; then
        chainA="A"
        chainB="I"
    elif [[ "$protein" == 1dfj* ||  "$protein" == 1r0r* ]]; then
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
        mpirun -np 4 /home/muskaan/imp-clean/build/setup_environment.sh python ~/projects/easal_related_stuff/Integrative_docking_benchmark/scripts/imp/modeling/sample_imp_distance_restraint.py "prod" "$chainA" "$chainB" "$folder_name" "DSSO" "32" 2> "err_$folder_name.log" &
    done

    cd ..
done

sleep 600;

cd DSSO
for protein in 1dfj_3 1dfj_9 1dfj_12 ; do
    cd "$protein"

    if [[ "$protein" == 1clv* ]]; then
        chainA="A"
        chainB="I"
    elif [[ "$protein" == 1dfj* ||  "$protein" == 1r0r* ]]; then
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
        mpirun -np 4 /home/muskaan/imp-clean/build/setup_environment.sh python ~/projects/easal_related_stuff/Integrative_docking_benchmark/scripts/imp/modeling/sample_imp_distance_restraint.py "prod" "$chainA" "$chainB" "$folder_name" "DSSO" "32" 2> "err_$folder_name.log" &
    done

    cd ..
done

sleep 600;

cd DSSO
for protein in 1kxp_4 1kxp_7 1kxp_11 ; do
    cd "$protein"

    if [[ "$protein" == 1clv* ]]; then
        chainA="A"
        chainB="I"
    elif [[ "$protein" == 1dfj* ||  "$protein" == 1r0r* ]]; then
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
        mpirun -np 4 /home/muskaan/imp-clean/build/setup_environment.sh python ~/projects/easal_related_stuff/Integrative_docking_benchmark/scripts/imp/modeling/sample_imp_distance_restraint.py "prod" "$chainA" "$chainB" "$folder_name" "DSSO" "32" 2> "err_$folder_name.log" &
    done

    cd ..
done

sleep 600;


cd DSSO
for protein in 1r0r_3 1r0r_7 2ayo_4 ; do
    cd "$protein"

    if [[ "$protein" == 1clv* ]]; then
        chainA="A"
        chainB="I"
    elif [[ "$protein" == 1dfj* ||  "$protein" == 1r0r* ]]; then
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
        mpirun -np 4 /home/muskaan/imp-clean/build/setup_environment.sh python ~/projects/easal_related_stuff/Integrative_docking_benchmark/scripts/imp/modeling/sample_imp_distance_restraint.py "prod" "$chainA" "$chainB" "$folder_name" "DSSO" "32" 2> "err_$folder_name.log" &
    done

    cd ..
done

sleep 600;




 cd DSSO
 for protein in 2ayo_8 2ayo_13 2b42_5 2b42_10 2hle_5 2hle_10 2hle_14 gata_gatc gcvpa_gcvpb phes_phet roca_putc sucd_succ ; do
     cd "$protein"

     if [[ "$protein" == 1clv* ]]; then
         chainA="A"
         chainB="I"
     elif [[ "$protein" == 1dfj* ||  "$protein" == 1r0r* ]]; then
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
         mpirun -np 4 /home/muskaan/imp-clean/build/setup_environment.sh python ~/projects/easal_related_stuff/Integrative_docking_benchmark/scripts/imp/modeling/sample_imp_distance_restraint.py "prod" "$chainA" "$chainB" "$folder_name" "DSSO" "32" 2> "err_$folder_name.log" &
     done

     cd ..
 done

 sleep 600;


###DMTMM
cd ../DMTMM/
for protein in 1clv_8 1dfj_4 1kxp_7 1r0r_6 2ayo_5 2b42_10 2hle_9; do
    cd "$protein"

    if [[ "$protein" == 1clv* ]]; thenzz
        chainA="A"
        chainB="I"
    elif [[ "$protein" == 1dfj* ||  "$protein" == 1r0r* ]]; then
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
        mpirun -np 4 /home/muskaan/imp-clean/build/setup_environment.sh python ~/projects/easal_related_stuff/Integrative_docking_benchmark/scripts/imp/modeling/sample_imp_distance_restraint.py "prod" "$chainA" "$chainB" "$folder_name" "DMTMM" "20" 2> "err_$folder_name.log" &
    done

    cd ..
done
