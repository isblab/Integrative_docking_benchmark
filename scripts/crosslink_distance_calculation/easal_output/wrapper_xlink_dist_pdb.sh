#TODO consistent naming of all wrapper scripts.Consider adding the name wrapper to the script file name to disambiguate from the script that
# calculates actual distances in the model
# crosslink_distances_easal_wrapper.sh


#!/bin/bash
###DSSO
#simulated
directories=("1clv_cl2" "1clv_cl6" "1dfj_cl3" "1dfj_cl9" "1dfj_cl12" "1kxp_cl4" "1kxp_cl7" "1kxp_cl11" "1r0r_cl3" "1r0r_cl7" "2ayo_cl4" "2ayo_cl8" "2ayo_cl13" "2b42_cl5" "2b42_cl10" "2hle_cl5" "2hle_cl10" "2hle_cl14")

for dir in "${directories[@]}"; do
    cd "DSSO/simulated/$dir"

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

    csv_file="$HOME/EASAL/benchmark/simulated/DSSO/$direc2/${filename}_DSSO_${number}.csv"
    python /home/muskaan/EASAL/scripts/crosslink_distance_calculation/easal_output/calc_xlink_distance_pdb.py "$chainA" "$chainB" "$csv_file"

    cd ../../..
done

# #experimental
directories=("gata_gatc_cl3" "gcvpa_gcvpb_cl5" "phes_phet_cl8" "roca_putc_cl2" "sucd_succ_cl4")

for dir in "${directories[@]}"; do
    cd "DSSO/experimental/$dir"

    filename=$(echo "$dir" | cut -d'_' -f1)
    filename+="_"
    filename+=$(echo "$dir" | cut -d'_' -f2)
    number=$(echo "$dir" | sed 's/.*_cl\([0-9]*\)/\1/')

    csv_file="$HOME/EASAL/benchmark/experimental/crosslinks/${filename}_DSSO_${number}.csv"
    python /home/muskaan/EASAL/scripts/crosslink_distance_calculation/easal_output/calc_xlink_distance_pdb.py "A" "B" "$csv_file"

    cd ../../..
done

###EDC
directories=("1clv_cl8" "1dfj_cl4" "1kxp_cl7" "1r0r_cl6" "2ayo_cl5" "2b42_cl10" "2hle_cl9")

for dir in "${directories[@]}"; do
    cd "EDC/$dir"

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

    csv_file="$HOME/EASAL/benchmark/simulated/EDC/${filename}_EDC_${number}.csv"
    python /home/muskaan/EASAL/scripts/crosslink_distance_calculation/easal_output/calc_xlink_distance_pdb.py "$chainA" "$chainB" "$csv_file"

    cd ../../
done
