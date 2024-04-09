################ Perform end-to-end analysis for IMP models ################
############## ------>"May the Force serve u well..." <------###############
############################################################################

############# One above all #############
##-------------------------------------##
import numpy as np
import sys
import os
import shutil
import subprocess
import glob

############ This is the way ############
##-------------------------------------##
# Check all the path variables before proceeding.

# Define the protein pairs and protein lengths for contact map script.
imp_path = "/home/muskaan/imp-clean/build/setup_environment.sh"
pdb_id = sys.argv[1]
xlinker = sys.argv[2]
modeling_dir_path = f"../../{xlinker}/{pdb_id}/"
density_file =glob.glob('density*.txt')[0]

def return_major_cluster():
	file = "model_analysis/summary_hdbscan_clustering.dat"

	with open (file, "r") as f:
		cluster_summary = f.readlines()

	models_count = {}
	#Get the index for the N_models header
	index = cluster_summary[0].split(",").index("N_models")

	# Get cluster no. and model counts for all clusters.
	for i in range(1,len(cluster_summary)):
		line = cluster_summary[i].split(",")
		if line[0] != "-1":
			models_count[int(line[0])] = int(line[index])
		else:
			continue                              #"Jaisa chal raha h, chalne do!!!"

	clust = [(i,models_count[i]) for i in range(len(models_count)) if models_count[i] == max(models_count.values())][0]
	print("Max models in cluster: ", clust[0], " no. of models = ", clust[1])
	return clust


############# PMI Analysis ##############
##-------------------------------------##
print("\n<-----------PMI Analysis----------->")
subprocess.call( [f"{imp_path}", "python", "/home/muskaan/EASAL/scripts/imp/analysis/run_analysis_trajectories.py", f"{modeling_dir_path}", f"{xlinker}"] )
cluster = return_major_cluster()


############ Model Extraction ###########
##-------------------------------------##
print("\n<-----------Model Extraction----------->")
if int( cluster[1] ) >= 30000:
	print(">30000 models....\n using Variable filter \n")

	subprocess.call( ["python", "/home/muskaan/EASAL/scripts/imp/analysis/variable_filter_v1.py", "-c", f"{cluster[0]}", "-g", "model_analysis"] )
	subprocess.call( [f"{imp_path}", "python", "/home/muskaan/EASAL/scripts/imp/analysis/run_extract_models.py", f"{modeling_dir_path}", f"{cluster[0]}"] )
	
else:
	subprocess.call( [f"{imp_path}", "python", "/home/muskaan/EASAL/scripts/imp/analysis/run_extract_models.py", f"{modeling_dir_path}", f"{cluster[0]}"] )


########### Sampcon analysis ############
##-------------------------------------##
print("\n<-----------Sampcon----------->")
sampcon_path = "/home/muskaan/imp-clean/imp/modules/sampcon/pyext/src/exhaust.py"
subprocess.call([
	f"{imp_path}",
	"python",
	f"{sampcon_path}",
	"-n", f"{pdb_id}",
	"-m", "cpu_omp",
	"-c", "2",
	"-d", f"{density_file}",
	"-gp",
	"-g", "2",
	"-sa", f"model_analysis/A_models_clust{cluster[0]}.txt",
	"-sb", f"model_analysis/B_models_clust{cluster[0]}.txt",
	"-ra", f"model_analysis/A_models_clust{cluster[0]}.rmf3",
	"-rb", f"model_analysis/B_models_clust{cluster[0]}.rmf3",
	"--align",
	# "--prism"
	])
