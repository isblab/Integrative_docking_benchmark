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
for protein in 1clv_2 1clv_6 1dfj_3 1dfj_9 1dfj_12 1kxp_4 1kxp_7 1kxp_11 1r0r_3 1r0r_7 2ayo_4 2ayo_8 2ayo_13 2b42_5 2b42_10 2hle_5 2hle_10 2hle_14; do
    crossLinkCount="crossLinkCount = $(echo "$protein" | cut -d'_' -f2)"
    proteinName=$(echo "$protein" | cut -d'_' -f1)
    number=$(echo "$protein" | cut -d'_' -f2)

    if [[ $number -lt 3 ]]; then
      crossLinkSatisfyThres='crossLinkSatisfyThres = 2'
    else
      crossLinkSatisfyThres="crossLinkSatisfyThres = $((number-2))"
    fi

    if [[ "$protein" == 1clv* ]]; then
      chainA='chain_A = "A"'
      chainB='chain_B = "I"'

      if [[ "$protein" == 1clv_2 ]]; then
        crossLinks='crossLinks = {233, 504, 262, 504}'

      elif [[ "$protein" == 1clv_6 ]]; then
        crossLinks='crossLinks = {188, 504, 233, 504, 300, 504, 303, 504, 300, 511, 262, 504}'
      fi

    elif [[ "$protein" == 1dfj* ||  "$protein" == 1r0r* ]]; then
      chainA='chain_A = "E"'
      chainB='chain_B = "I"'

      if [[ "$protein" == 1dfj_3 ]]; then
        crossLinks='crossLinks = {31, 84, 31, 81, 91, 290}'

      elif [[ "$protein" == 1dfj_9 ]]; then
        crossLinks='crossLinks = {91, 278, 61, 347, 31, 141, 37, 84, 104, 347, 31, 84, 98, 316, 37, 336, 31, 449}'

      elif [[ "$protein" == 1dfj_12 ]]; then
        crossLinks='crossLinks = {31, 141, 37, 254, 98, 316, 98, 290, 61, 347, 98, 347, 1, 449, 104, 347, 31, 81, 37, 141, 104, 316, 31, 449}'

      elif [[ "$protein" == 1r0r_3 ]]; then
        crossLinks='crossLinks = {94, 34, 170, 55, 170, 29}'

      elif [[ "$protein" == 1r0r_7 ]]; then
        crossLinks='crossLinks = {94, 13, 136, 13, 170, 13, 136, 29, 170, 29, 94, 34, 170, 55}'
      fi

    elif [[ "$protein" == 1kxp* ]]; then
      chainA='chain_A = "A"'
      chainB='chain_B = "D"'
      if [[ "$protein" == 1kxp_4 ]]; then
        crossLinks='crossLinks = {326, 419, 113, 284, 328, 419, 328, 150}'

      elif [[ "$protein" == 1kxp_7 ]]; then
        crossLinks='crossLinks = {315, 421, 328, 207, 291, 436, 315, 420, 326, 150, 284, 303, 113, 295}'

      elif [[ "$protein" == 1kxp_11 ]]; then
        crossLinks='crossLinks = {291, 419, 291, 303, 291, 153, 291, 420, 315, 444, 328, 444, 326, 420, 284, 419, 118, 292, 328, 420, 328, 207}'
      fi

    else
      chainA='chain_A = "A"'
      chainB='chain_B = "B"'

      if [[ "$protein" == 2ayo_4 ]]; then
        crossLinks='crossLinks = {188, 48, 299, 11, 333, 63, 335, 29}'

      elif [[ "$protein" == 2ayo_8 ]]; then
        crossLinks='crossLinks = {315, 63, 335, 27, 335, 33, 188, 48, 239, 48, 299, 48, 255, 33, 237, 6}'

      elif [[ "$protein" == 2ayo_13 ]]; then
        crossLinks='crossLinks = {213, 48, 335, 29, 315, 6, 312, 6, 335, 27, 315, 63, 239, 6, 255, 48, 333, 63, 315, 11, 341, 11, 333, 29, 333, 33}'

      elif [[ "$protein" == 2b42_5 ]]; then
        crossLinks='crossLinks = {9, 95, 111, 135, 111, 95, 111, 154, 111, 99}'

      elif [[ "$protein" == 2b42_10 ]]; then
        crossLinks='crossLinks = {286, 40, 9, 95, 111, 95, 115, 95, 315, 95, 111, 99, 111, 135, 9, 154, 111, 154, 286, 154}'

      elif [[ "$protein" == 2hle_5 ]]; then
        crossLinks='crossLinks = {31, 60, 162, 80, 149, 86, 31, 44, 31, 131}'

      elif [[ "$protein" == 2hle_10 ]]; then
        crossLinks='crossLinks = {31, 71, 24, 71, 149, 116, 158, 112, 158, 80, 158, 86, 162, 80, 149, 166, 149, 60, 158, 133}'

      elif [[ "$protein" == 2hle_14 ]]; then
        crossLinks='crossLinks = {149, 166, 158, 80, 162, 80, 149, 80, 24, 95, 31, 116, 24, 71, 31, 44, 31, 71, 149, 60, 24, 112, 24, 67, 24, 80, 149, 131}'
      fi

    fi

    mkdir DSSO/$protein
    sed -i "s/$def_pdbfile/${protein:0:4}/g;
            s/$def_chainA/$chainA/g;
            s/$def_chainB/$chainB/g;
            s/$def_crossLinkCount/$crossLinkCount/g;
            s/$def_crossLinkSatisfyThres/$crossLinkSatisfyThres/g;
            s/$def_crossLinks/$crossLinks/g" "/home/muskaan/easal-dev/settings.ini"

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

###### DSSO experimental ######

for protein in gata_gatc_3 gcvpa_gcvpb_5 phes_phet_8 roca_putc_2 sucd_succ_4; do
    chainA='chain_A = "A"'
    chainB='chain_B = "B"'

    crossLinkCount="crossLinkCount = $(echo "$protein" | cut -d'_' -f3)"
    number=$(echo "$protein" | cut -d'_' -f3)

    proteinName=$(echo "$protein" | cut -d'_' -f1)
    proteinName+="_"
    proteinName+=$(echo "$protein" | cut -d'_' -f2)

    if [[ $number -lt 3 ]]; then
      crossLinkSatisfyThres='crossLinkSatisfyThres = 2'
    else
      crossLinkSatisfyThres="crossLinkSatisfyThres = $((number-2))"
    fi

    if [[ "$protein" == gata_gatc ]]; then
      crossLinks='crossLinks = {85, 70, 85, 79, 344, 10}'

    elif [[ "$protein" == gcvpa_gcvpb ]]; then
      crossLinks='crossLinks = {305, 2, 177, 301, 173, 301, 305, 468, 45, 306}'

    elif [[ "$protein" == phes_phet ]]; then
      crossLinks='crossLinks = {78, 650, 85, 650, 87, 679, 321, 502, 342, 610, 181, 593, 75, 650, 237, 502}'

    elif [[ "$protein" == roca_putc ]]; then
      crossLinks='crossLinks = {272, 272, 6, 158}'

    elif [[ "$protein" == sucd_succ ]]; then
      crossLinks='crossLinks = {256, 146, 256, 56, 256, 1, 243, 1}'
    fi

    mkdir DSSO/$protein
    sed -i "s/$def_pdbfile/$proteinName/g;
            s/$def_chainA/$chainA/g;
            s/$def_chainB/$chainB/g;
            s/$def_crossLinkCount/$crossLinkCount/g;
            s/$def_crossLinkSatisfyThres/$crossLinkSatisfyThres/g;
            s/$def_crossLinks/$crossLinks/g" "/home/muskaan/easal-dev/settings.ini"

    cp -r settings.ini DSSO/$protein
    build/easal >> DSSO/$protein/"$log_file" 2>&1
    mv $proteinName*.txt DSSO/$protein/
    def_pdbfile=$proteinName
    def_chainA=$chainA
    def_chainB=$chainB
    def_crossLinkCount=$crossLinkCount
    def_crossLinkSatisfyThres=$crossLinkSatisfyThres
    def_crossLinks=$crossLinks

done

###### EDC simulated ######
mkdir EDC
for protein in 1clv_8 1dfj_4 1kxp_7 1r0r_6 2ayo_5 2b42_10 2hle_9; do
  crossLinkCount="crossLinkCount = $(echo "$protein" | cut -d'_' -f2)"
  number=$(echo "$protein" | cut -d'_' -f2)
  crossLinkSatisfyThres="crossLinkSatisfyThres = $((number-2))"
  proteinName=$(echo "$protein" | cut -d'_' -f1)

  if [[ "$protein" == 1clv* ]]; then
    chainA='chain_A = "A"'
    chainB='chain_B = "I"'
    crossLinks='crossLinks = {141, 519, 135, 526, 141, 513, 114, 513, 143, 513, 116, 513, 149, 513, 135, 519}'

  elif [[ "$protein" == 1dfj* ||  "$protein" == 1r0r* ]]; then
    chainA='chain_A = "E"'
    chainB='chain_B = "I"'

    if [[ "$protein" == 1dfj_4 ]]; then
      crossLinks='crossLinks = {111, 431, 121, 439, 121, 374, 2, 437}'

    elif [[ "$protein" == 1r0r_6 ]]; then
      crossLinks='crossLinks = {54, 7, 60, 10, 181, 19, 60, 7, 60, 19, 54, 10}'
    fi

  elif [[ "$protein" == 1kxp* ]]; then
    chainA='chain_A = "A"'
    chainB='chain_B = "D"'
    crossLinks='crossLinks = {361, 298, 316, 417, 292, 132, 311, 142, 286, 412, 270, 424, 167, 138}'

  else
    chainA='chain_A = "A"'
    chainB='chain_B = "B"'

    if [[ "$protein" == 2ayo_5 ]]; then
      crossLinks='crossLinks = {296, 34, 334, 32, 144, 51, 201, 39, 295, 16}'

    elif [[ "$protein" == 2hle_9 ]]; then
      crossLinks='crossLinks = {39, 119, 29, 62, 29, 69, 154, 87, 62, 59, 62, 128, 29, 34, 43, 59, 51, 128}'

    elif [[ "$protein" == 2b42_10 ]]; then
      crossLinks='crossLinks = {355, 11, 217, 11, 355, 78, 284, 121, 320, 121, 354, 11, 284, 119, 305, 11, 309, 11, 10, 121}'
    fi

  fi

  mkdir EDC/$protein
  sed -i "s/$def_pdbfile/${protein:0:4}/g;
          s/$def_chainA/$chainA/g;
          s/$def_chainB/$chainB/g;
          s/$def_crossLinkCount/$crossLinkCount/g;
          s/$def_crossLinkSatisfyThres/$crossLinkSatisfyThres/g;
          s/$def_activeUpperDelta/20/g;
          s/$def_crossLinks/$crossLinks/g" "/home/muskaan/easal-dev/settings.ini"

  cp -r settings.ini EDC/$protein
  build/easal >> EDC/$protein/"$log_file" 2>&1
  mv $proteinName*.txt EDC/$protein/
  def_pdbfile=${protein:0:4}
  def_chainA=$chainA
  def_chainB=$chainB
  def_crossLinkCount=$crossLinkCount
  def_crossLinkSatisfyThres=$crossLinkSatisfyThres
  def_crossLinks=$crossLinks

done
