#!/bin/bash

patterns=(
    "1clv_A A 1 471"
    "1clv_B I 1 32"
    "1r0r_A E 106 379"
    "1r0r_B I 135 185"
    "2ayo_A A 91 494"
    "2ayo_B B 609 684"
    "1kxp_A A 3 377"
    "1kxp_B D 17 474"
    "2hle_A A 17 196"
    "2hle_B B 28 165"
    "1dfj_A E 27 150"
    "1dfj_B I 1 456"
    "2b42_A A 22 402"
    "2b42_B B 29 213"
)

# Iterate over the patterns and execute the Python script
for pattern in "${patterns[@]}"; do
    set -- $pattern
    pdb_file="${1}.pdb"
    python ../../../../scripts/get_selected_regions_from_AF2_structure.py "$pdb_file" "$2" "$3" "$4"
done
