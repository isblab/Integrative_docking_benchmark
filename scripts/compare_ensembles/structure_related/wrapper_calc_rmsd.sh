#!/bin/bash


for case in "1r0r_DSSO_3" "1clv_DSSO_2" "1kxp_DSSO_4" "2ayo_DSSO_4" "2b42_DSSO_5" "1dfj_DSSO_3" "2hle_DSSO_5" "1dfj_DMTMM_4" "1clv_DMTMM_8" "1kxp_DMTMM_7" "1r0r_DMTMM_6" "2ayo_DMTMM_5" "2b42_DMTMM_10" "2hle_DMTMM_9" "1dfj_DSSO_9" "1clv_DSSO_6" "1kxp_DSSO_7" "1r0r_DSSO_7" "2ayo_DSSO_8" "2b42_DSSO_10" "2hle_DSSO_10" "1dfj_DSSO_12" "1kxp_DSSO_11" "2ayo_DSSO_13" "2hle_DSSO_14" "gata_gatc_DSSO_3" "gcvpa_gcvpb_DSSO_5" "roca_putc_DSSO_2" "sucd_succ_DSSO_4" "phes_phet_DSSO_8"; do
    if [ ${#case} -lt 15 ]; then
        if [[ $case == *"DSSO"* ]]; then
            rmf_file="~/easal/imp_output/DSSO_analysis/${case%DSSO*}${case##*_}/sampcon_0_extracted.rmf3"
            easal_output_direc="~/easal/easal_output/DSSO/simulated/${case%DSSO*}cl${case##*_}/"
        elif [[ $case == *"DMTMM"* ]]; then
            rmf_file="~/easal/imp_output/DMTMM_analysis/${case%DMTMM*}${case##*_}/sampcon_0_extracted.rmf3"
            easal_output_direc="~/easal/easal_output/DMTMM/${case%DMTMM*}cl${case##*_}/"
        fi
    else
        rmf_file="~/easal/imp_output/DSSO_analysis/${case%DSSO*}${case##*_}/sampcon_0_extracted.rmf3"
        easal_output_direc="~/easal/easal_output/DSSO/experimental/${case%DSSO*}cl${case##*_}/"
    fi

    if [[ $case == *"1dfj"* || $case == *"1r0r"* ]]; then
        native_pdbfile="~/easal-dev/files/${case:0:4}.pdb"
        chain_A="E"
        chain_B="I"
    elif [[ $case == *"1kxp"* ]]; then
        native_pdbfile="~/easal-dev/files/1kxp.pdb"
        chain_A="A"
        chain_B="D"
    elif [[ $case == *"2hle"* || $case == *"2b42"* || $case == *"2ayo"* ]]; then
        native_pdbfile="~/easal-dev/files/${case:0:4}.pdb"
        chain_A="A"
        chain_B="B"
    elif [[ $case == *"1clv"* ]]; then
        native_pdbfile="~/easal-dev/files/${case:0:4}.pdb"
        chain_A="A"
        chain_B="I"
    else
        native_pdbfile="~/easal-dev/files/${case%_DSSO*}.pdb"
        chain_A="A"
        chain_B="B"
    fi

    outputfile="$case.txt"

    ~/imp-clean/build/setup_environment.sh python ~/Integrative_docking_benchmark/scripts/analysis/structure_related/calc_rmsd.py "$native_pdbfile" "$chain_A" "$chain_B" "$rmf_file" "$easal_output_direc" "$outputfile" &
done
