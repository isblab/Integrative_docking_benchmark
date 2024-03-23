#TODO rename file to be more descriptive and unambigious. master_script.sh means what? 


#!/bin/bash
###DSSO
##simulated
directories=("1clv_2" "1clv_6" "1dfj_3" "1dfj_9" "1dfj_12" "1kxp_4" "1kxp_7" "1kxp_11" "1r0r_3" "1r0r_7" "2ayo_4" "2ayo_8" "2ayo_13" "2b42_5" "2b42_10" "2hle_5" "2hle_10" "2hle_14")

for dir in "${directories[@]}"; do
  filename=$(echo "$dir" | cut -d'_' -f1)
  number=$(echo "$dir" | cut -d'_' -f2)

  if (( number <= 5 )); then
      direc2="less_than_5"
  elif (( number >= 6 && number <= 10 )); then
      direc2="6-10"
  else
      direc2="more_than_10"
  fi
  csv_file="$HOME/EASAL/benchmark/simulated/DSSO/$direc2/${filename}_DSSO_${number}.csv"

  ~/imp-clean/build/setup_environment.sh python ~/EASAL/scripts/crosslink_distance_calculation/imp_output/get_xlviol_val_set_v2_easal.py "../DSSO_analysis/$dir/sampcon_0_extracted.rmf3" "$csv_file" "32"

done

##experimental
directories=("gata_gatc_3" "gcvpa_gcvpb_5" "phes_phet_8" "roca_putc_2" "sucd_succ_4")

for dir in "${directories[@]}"; do
    filename=$(echo "$dir" | cut -d'_' -f1)
    filename+="_"
    filename+=$(echo "$dir" | cut -d'_' -f2)
    number=$(echo "$dir" | cut -d'_' -f3)

    csv_file="$HOME/EASAL/benchmark/experimental/crosslinks/${filename}_DSSO_${number}.csv"
    ~/imp-clean/build/setup_environment.sh python ~/EASAL/scripts/crosslink_distance_calculation/imp_output/get_xlviol_val_set_v2_easal.py "../DSSO_analysis/$dir/sampcon_0_extracted.rmf3" "$csv_file" "32"
done

###EDC
directories=("1clv_8" "1dfj_4" "1kxp_7" "1r0r_6" "2ayo_5" "2b42_10" "2hle_9")

for dir in "${directories[@]}"; do
  filename=$(echo "$dir" | cut -d'_' -f1)
  number=$(echo "$dir" | cut -d'_' -f2)
  csv_file="$HOME/EASAL/benchmark/simulated/EDC/${filename}_EDC_${number}.csv"

  ~/imp-clean/build/setup_environment.sh python ~/EASAL/scripts/crosslink_distance_calculation/imp_output/get_xlviol_val_set_v2_easal.py "../EDC_analysis/$dir/sampcon_0_extracted.rmf3" "$csv_file" "20"

done
