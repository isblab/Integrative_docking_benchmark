#!/bin/bash

# for EDC
#go to /home/muskaan/easal_benchmark_JWALK_results/JWALK_results_EDC/ directory and run the following:
patterns=(
    "1clv 8 EDC"
    "1r0r 6 EDC"
    "2ayo 5 EDC"
    "1kxp 7 EDC"
    "2hle 9 EDC"
    "1dfj 4 EDC"
    "2b42 10 EDC"
)

# Iterate over the patterns and execute the Python script
for pattern in "${patterns[@]}"; do
    set -- $pattern
    pdb_file="${1}_interprotein.csv"
    python ~/EASAL/scripts/2_getting_random_xlinks.py "$pdb_file" "$2" "$3"
done

# for DSSO
#go to /home/muskaan/easal_benchmark_JWALK_results/JWALK_results_DSSO/ directory and run the following:
# patterns=(
#     "1clv 2 DSSO"
#     "1clv 6 DSSO"
#     "1r0r 3 DSSO"
#     "1r0r 7 DSSO"
#     "2ayo 4 DSSO"
#     "2ayo 8 DSSO"
#     "2ayo 13 DSSO"
#     "1kxp 4 DSSO"
#     "1kxp 7 DSSO"
#     "1kxp 11 DSSO"
#     "2hle 5 DSSO"
#     "2hle 10 DSSO"
#     "2hle 14 DSSO"
#     "1dfj 3 DSSO"
#     "1dfj 9 DSSO"
#     "1dfj 12 DSSO"
#     "2b42 5 DSSO"
#     "2b42 10 DSSO"
# )
#
# # Iterate over the patterns and execute the Python script
# for pattern in "${patterns[@]}"; do
#     set -- $pattern
#     pdb_file="${1}_crosslink_list_interprotein_crosslinks.csv"
#     python ~/EASAL/scripts/2_getting_random_xlinks.py "$pdb_file" "$2" "$3"
# done
