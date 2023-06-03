#!/bin/bash
# This file is part of EASAL.
#
# EASAL is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# EASAL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


ncores=(2 4) 
stepsizes=(1.2)
maxSamplers=(10 100)
molA_files=("files/3AtomsA.pdb")
molB_files=("files/3AtomsB.pdb")

#molA_files=("files/T2B.pdb")
#molB_files=("files/T2B.pdb")

#molA_files=("files/3AtomsA.pdb" "files/6AtomsA.pdb")
#molB_files=("files/3AtomsB.pdb" "files/6AtomsB.pdb")
if [ ${#molA_files[@]} -ne ${#molB_files[@]} ]; then
    echo "Molecule lists of different lenghts!";
    exit 1
fi

if [ ! -d "SettingsFiles" ]; then
    echo "Directory SettingsFiles does not exist."
    exit 1
fi

if [ ! -d "DriverData" ]; then
    echo "Directory DriverData does not exist."
    exit 1
fi

base_settings_file="settings.ini"

uname=$(whoami)
fileName='easal_benchmark'
fileName+=_`date +%F_%T`
fileName+='.csv'
touch $fileName

for core in "${ncores[@]}" #Run it for different number of 
do
    for step in "${stepsizes[@]}" #Run it for different step sizes
    do  
        for (( i=0; i<${#molA_files[@]}; i++ ));
        do
			for samplers in "${maxSamplers[@]}"
			do
            #echo $i ${molA_files[$i]} ${molB_files[$i]}
            molA="$(cut -d'/' -f2<<<${molA_files[$i]} | cut -d'.' -f1)"
            molB="$(cut -d'/' -f2<<<${molB_files[$i]} | cut -d'.' -f1)"
            #echo $molA $molB
            
            jobName="${molA}_${molB}_${step}_${core}cores_${samplers}Samplers"
            jobName+=_`date +%F_%T`

            if [ -d "DriverData/$jobName" ]; then
                echo "Filename clash: DriverData/$jobName"
                exit 1
            fi

            #mkdir "DriverData/$jobName"

            settings_file=SettingsFiles/settings_${jobName}.ini
            mol_A=$(echo \"${molA_files[$i]}\" | sed 's/\//\\\//g')
            mol_B=$(echo \"${molB_files[$i]}\" | sed 's/\//\\\//g')

            sed "/stepSize/s/= .*/= ${step}/
				 /MaxSamplers/s/= .*/= ${samplers}/
                /dataDirectory/s/= .*/= \"DriverData\/${jobName}\"/
                9 s/= .*/= ${mol_A}/
                20 s/= .*/= ${mol_B}/
                /sessionId/s/= .*/= \"${jobName}\"/" $base_settings_file > $settings_file

    		output=$(build/easal -settings ${settings_file} -max-threads ${core}) 

			done
        done
	done
done
