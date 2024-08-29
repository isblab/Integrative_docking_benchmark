#!/bin/bash
###DSSO
##simulated
directories=("1clv_2" "1clv_6" "1dfj_3" "1dfj_9" "1dfj_12" "1kxp_4" "1kxp_7" "1kxp_11" "1r0r_3" "1r0r_7" "2ayo_4" "2ayo_8" "2ayo_13" "2b42_5" "2b42_10" "2hle_5" "2hle_10" "2hle_14")

for dir in "${directories[@]}"; do
  filename=$(echo "$dir" | cut -d'_' -f1)

  if [[ "$filename" == 1clv* ]]; then
      chainA="A"
      chainB="I"
  elif [[ "$filename" == 1dfj* ||  "$filename" == 1r0r* ]]; then
      chainA="E"
      chainB="I"
  elif [[ "$filename" == 1kxp* ]]; then
      chainA="A"
      chainB="D"
  else
      chainA="A"
      chainB="B"
  fi

  ~/imp-clean/build/setup_environment.sh python ~/Integrative_docking_benchmark/scripts/imp/filtering_by_clashes.py "imp_output/DSSO_analysis/$dir/sampcon_0_extracted.rmf3" "imp_output/DSSO/$dir/${filename}.pdb" "$chainA" "$chainB" "filtered_rmfs/DSSO/${dir}_filtered.rmf3" &

done

##experimental
directories=("gata_gatc_3" "gcvpa_gcvpb_5" "phes_phet_8" "roca_putc_2" "sucd_succ_4")

for dir in "${directories[@]}"; do
    filename=$(echo "$dir" | cut -d'_' -f1)
    filename+="_"
    filename+=$(echo "$dir" | cut -d'_' -f2)
    number=$(echo "$dir" | cut -d'_' -f3)

    ~/imp-clean/build/setup_environment.sh python ~/Integrative_docking_benchmark/scripts/imp/filtering_by_clashes.py "imp_output/DSSO_analysis/$dir/sampcon_0_extracted.rmf3" "imp_output/DSSO/$dir/${filename}.pdb" "A" "B" "filtered_rmfs/DSSO/${dir}_filtered.rmf3" &

done

###DMTMM
directories=("1clv_8" "1dfj_4" "1kxp_7" "1r0r_6" "2ayo_5" "2b42_10" "2hle_9")

for dir in "${directories[@]}"; do
  filename=$(echo "$dir" | cut -d'_' -f1)
  number=$(echo "$dir" | cut -d'_' -f2)

  if [[ "$filename" == 1clv* ]]; then
      chainA="A"
      chainB="I"
  elif [[ "$filename" == 1dfj* ||  "$filename" == 1r0r* ]]; then
      chainA="E"
      chainB="I"
  elif [[ "$filename" == 1kxp* ]]; then
      chainA="A"
      chainB="D"
  else
      chainA="A"
      chainB="B"
  fi

  ~/imp-clean/build/setup_environment.sh python ~/Integrative_docking_benchmark/scripts/imp/filtering_by_clashes.py "imp_output/DMTMM_analysis/$dir/sampcon_0_extracted.rmf3" "imp_output/DMTMM/$dir/${filename}.pdb" "$chainA" "$chainB" "filtered_rmfs/DMTMM/${dir}_filtered.rmf3" &

done
