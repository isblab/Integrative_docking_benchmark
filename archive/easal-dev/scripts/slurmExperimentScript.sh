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


ncores=(10 9 8 7 6 5 4 3 2 1) 
#ncores=(10 8 6 4 2) 
#ncores=(8) 

#stepsizes=(1.6)
stepsizes=(1.6 0.8)

#maxSamplers=(100)
maxSamplers=(30 100 300 1000 10000)
#molA_files=("files/6AtomsA.pdb")
#molB_files=("files/6AtomsB.pdb")

molA_files=("files/10AtomsA.pdb")
molB_files=("files/10AtomsB.pdb")

mts=(1)

#molA_files=("files/3AtomsA.pdb" "files/6AtomsA.pdb" "files/6pocketedA.pdb" "files/10AtomsA.pdb")
#molB_files=("files/3AtomsB.pdb" "files/6AtomsB.pdb" "files/6pocketedB.pdb" "files/10AtomsB.pdb")

#molA_files=("files/3AtomsA.pdb" "files/6AtomsA.pdb" "files/6pocketedA.pdb") 
#molB_files=("files/3AtomsB.pdb" "files/6AtomsB.pdb" "files/6pocketedB.pdb")
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
				
				for mt in "${mts[@]}"
				do
					#max_threads=$mt
					max_threads=$((mt * core))
					jobName="${molA}_${molB}_${step}_${core}cores_${samplers}Samplers_${max_threads}Threads_"
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

					mem=20gb
					echo $jobName
				
					#Run the program
					#echo "sbatch --cpus-per-task=$core slurm.sh ${settings_file}"

					if [ $max_threads -eq 0 ]
					then
						echo "sbatch --cpus-per-task=$core --mem=$mem --job-name=$jobName MPI.sh ${settings_file} with auto option"
						output=$(sbatch --cpus-per-task=$core --mem=$mem --job-name=$jobName MPI.sh ${settings_file} ${jobName}) 
					else
						echo "sbatch --cpus-per-task=$core --mem=$mem --job-name=$jobName MPI.sh ${settings_file} $max_threads"
						output=$(sbatch --cpus-per-task=$core --mem=$mem --job-name=$jobName MPI.sh ${settings_file} ${jobName} $max_threads) 
					fi
					
					#Monitor it till its done
					jobID="$(cut -d' ' -f4<<<$output)"
					echo "JobId: $jobID"
					running=1
					while [ $running -eq 1 ]
					do
						squeueOutput=$(squeue -j $jobID)
						#echo "$squeueOutput" | awk '{print $13}' | grep R && running=1 || running=0
						echo "$squeueOutput" | grep $jobID && running=1 || running=0
						sleep 1m
					done
					jobStatus=$(sacct -j $jobID --format=TotalCPU,Elapsed,MaxRSS |tail -1)
					cpuTime=$(echo $jobStatus|cut -d' ' -f1)
					wallTime=$(echo $jobStatus|cut -d' ' -f2)
					memory=$(echo $jobStatus|cut -d' ' -f3)
					numNodes=$(wc -l DriverData\/${jobName}/RoadMap.txt | cut -d' ' -f 1)

					echo $jobName,$molA-$molB,$step,$core,$samplers,$wallTime,$cpuTime,$memory,$numNodes,$max_threads >> $fileName
				done
			done
        done
	done
done
