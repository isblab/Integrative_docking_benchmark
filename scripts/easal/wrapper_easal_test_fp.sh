#!/bin/bash

#Defaults. Note: Can change these parameters based on the last settings in settings.ini
def_pdbfile='phes_phet'
def_chainA='chain_A = "A"'
def_chainB='chain_B = "B"'
def_activeUpperDelta='32'
def_crossLinkCount='crossLinkCount = 8'
def_crossLinkSatisfyThres='crossLinkSatisfyThres = 6'
def_crossLinks='crossLinks = {78, 650, 85, 650, 87, 679, 321, 502, 342, 610, 181, 593, 75, 650, 237, 502}'
log_file="logfile.txt"

###### DSSO simulated ######
mkdir DSSO
for protein in 1clv_5 1dfj_10 1dfj_11 1kxp_10 1r0r_5 2ayo_10 2ayo_15 2b42_10 2hle_10 2hle_15; do
    crossLinkCount="crossLinkCount = $(echo "$protein" | cut -d'_' -f2)"
    proteinName=$(echo "$protein" | cut -d'_' -f1)
    number=$(echo "$protein" | cut -d'_' -f2)

    # if [[ $number -lt 3 ]]; then
    #   crossLinkSatisfyThres='crossLinkSatisfyThres = 2'
    # else
    #   crossLinkSatisfyThres="crossLinkSatisfyThres = $((number-2))"
    # fi

    if [[ "$protein" == 1clv* ]]; then
      chainA='chain_A = "A"'
      chainB='chain_B = "I"'

      if [[ "$protein" == 1clv_5 ]]; then
        crossLinks='crossLinks = {188, 504, 233, 504, 300, 504, 303, 504, 209, 519}'
        crossLinkSatisfyThres='crossLinkSatisfyThres = 3'
      fi

    elif [[ "$protein" == 1dfj* ||  "$protein" == 1r0r* ]]; then
      chainA='chain_A = "E"'
      chainB='chain_B = "I"'

      if [[ "$protein" == 1dfj_11 ]]; then
        crossLinks='crossLinks = {61, 347, 98, 347, 1, 449, 104, 347, 31, 81, 37, 141, 104, 316, 31, 449, 122, 302, 9, 189, 1, 226}'
        crossLinkSatisfyThres='crossLinkSatisfyThres = 7'

      elif [[ "$protein" == 1dfj_10 ]]; then
        crossLinks='crossLinks = {91, 278, 61, 347, 31, 141, 37, 84, 104, 347, 31, 84, 98, 316, 37, 336, 16, 190, 115, 334}'
        crossLinkSatisfyThres='crossLinkSatisfyThres = 7'

      elif [[ "$protein" == 1r0r_5 ]]; then
        crossLinks='crossLinks = {94, 13, 136, 13, 170, 13, 170, 29, 275, 56}'
        crossLinkSatisfyThres='crossLinkSatisfyThres = 3'
      fi

    elif [[ "$protein" == 1kxp* ]]; then
      chainA='chain_A = "A"'
      chainB='chain_B = "D"'

      if [[ "$protein" == 1kxp_10 ]]; then
        crossLinks='crossLinks = {315, 421, 328, 207, 291, 436, 315, 420, 326, 150, 284, 303, 113, 295, 291, 150, 206, 295, 263, 264}'
        crossLinkSatisfyThres='crossLinkSatisfyThres = 7'
      fi

    else
      chainA='chain_A = "A"'
      chainB='chain_B = "B"'

      if [[ "$protein" == 2ayo_10 ]]; then
        crossLinks='crossLinks = {315, 63, 335, 27, 335, 33, 188, 48, 239, 48, 299, 48, 255, 33, 237, 6, 368, 1, 443, 38}'
        crossLinkSatisfyThres='crossLinkSatisfyThres = 7'

      elif [[ "$protein" == 2ayo_15 ]]; then
        crossLinks='crossLinks = {213, 48, 335, 29, 315, 6, 312, 6, 335, 27, 315, 63, 239, 6, 255, 48, 333, 63, 315, 11, 341, 11, 333, 29, 464, 20, 419, 29, 454, 11}'
        crossLinkSatisfyThres='crossLinkSatisfyThres = 11'

      elif [[ "$protein" == 2b42_10 ]]; then
        crossLinks='crossLinks = {286, 40, 9, 95, 111, 95, 115, 95, 315, 95, 111, 99, 111, 135, 9, 154, 144, 52, 35, 27}'
        crossLinkSatisfyThres='crossLinkSatisfyThres = 7'

      elif [[ "$protein" == 2hle_10 ]]; then
        crossLinks='crossLinks = {31, 71, 24, 71, 149, 116, 158, 112, 158, 80, 158, 86, 162, 80, 149, 166, 133, 69, 96, 78}'
        crossLinkSatisfyThres='crossLinkSatisfyThres = 7'

      elif [[ "$protein" == 2hle_15 ]]; then
        crossLinks='crossLinks = {149, 166, 158, 80, 162, 80, 149, 80, 24, 95, 31, 116, 24, 71, 31, 44, 31, 71, 149, 60, 24, 112, 24, 67, 174, 52, 179, 39, 136, 46}'
        crossLinkSatisfyThres='crossLinkSatisfyThres = 11'
      fi

    fi

    mkdir DSSO/$protein
    sed -i "s/$def_pdbfile/${protein:0:4}/g;
            s/$def_chainA/$chainA/g;
            s/$def_chainB/$chainB/g;
            s/$def_crossLinkCount/$crossLinkCount/g;
            s/$def_crossLinkSatisfyThres/$crossLinkSatisfyThres/g;
            s/$def_crossLinks/$crossLinks/g" "settings.ini"

    cp -r settings.ini DSSO/$protein
    build/easal >> DSSO/$protein/"$log_file" 2>&1
    mv $proteinName*.txt DSSO/$protein/
    def_pdbfile=${protein:0:4}
    def_chainA=$chainA
    def_chainB=$chainB
    def_crossLinkCount=$crossLinkCount
    def_crossLinkSatisfyThres=$crossLinkSatisfyThres
    def_crossLinks=$crossLinks

done
