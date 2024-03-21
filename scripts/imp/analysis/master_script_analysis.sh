#!/bin/bash
# ###DSSO
mkdir DSSO_analysis
cd DSSO_analysis
#
for protein in 1clv_2 1clv_6 1dfj_3 1dfj_9 1dfj_12 1kxp_4 1kxp_7 1kxp_11 1r0r_3 1r0r_7 2ayo_4 2ayo_8 2ayo_13 2b42_5 2b42_10 2hle_5 2hle_10 2hle_14 gata_gatc_3 gcvpa_gcvpb_5 phes_phet_8 roca_putc_2 sucd_succ_4; do
    mkdir "$protein"
    cd "$protein"
    if [[ "$protein" == 1clv* ]]; then
        cp -r ../../density_A_I.txt .
    elif [[ "$protein" == 1dfj* ||  "$protein" == 1r0r* ]]; then
        cp -r ../../density_E_I.txt .
    elif [[ "$protein" == 1kxp* ]]; then
        cp -r ../../density_A_D.txt .
    else
        cp -r ../../density_A_B.txt .
    fi

    echo "$protein....."
    python ~/EASAL/scripts/imp/analysis/end_to_end_analysis.py "$protein" "DSSO"

    cd ..
done

###EDC
mkdir ../EDC_analysis
cd ../EDC_analysis

for protein in 1clv_8 1dfj_4 1kxp_7 1r0r_6 2ayo_5 2b42_10 2hle_9; do
    mkdir "$protein"
    cd "$protein"

    if [[ "$protein" == 1clv* ]]; then
        cp -r ../../density_A_I.txt .
    elif [[ "$protein" == 1dfj* ||  "$protein" == 1r0r* ]]; then
        cp -r ../../density_E_I.txt .
    elif [[ "$protein" == 1kxp* ]]; then
        cp -r ../../density_A_D.txt .
    else
        cp -r ../../density_A_B.txt .
    fi
    echo "$protein....."
    python ~/EASAL/scripts/imp/analysis/end_to_end_analysis.py "$protein" "EDC"

    cd ..
done
