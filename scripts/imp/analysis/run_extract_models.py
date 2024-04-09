import numpy as np
import pandas as pd
import math
import glob
import sys
import os

sys.path.append('/home/muskaan/imp-clean/PMI_analysis/pyext/src/')
from analysis_trajectories import *

#################################
########### MAIN ################
#################################
nproc = 30
top_dir = sys.argv[1]
analys_dir = './model_analysis/'

# Check if analysis dir exists
# if not os.path.isdir(analys_dir):
#     os.makedirs(analys_dir)

# How are the trajectories dir names
dir_head = 'run_'
out_dirs = glob.glob(top_dir+'/'+dir_head+'*/')
print(out_dirs)
c = sys.argv[2]
################################
# Extract frames
################################
# Load module
AT = AnalysisTrajectories(
    out_dirs, dir_name=dir_head, analysis_dir=analys_dir, nproc=nproc
)

HA = AT.get_models_to_extract(
        analys_dir + "good_scoring_models_A_cluster" + str(c) + "_detailed.csv")

HB = AT.get_models_to_extract(
        analys_dir + "good_scoring_models_B_cluster" + str(c) + "_detailed.csv")

rmf_file_out_A = "A_models_clust" + str(c) + ".rmf3"
rmf_file_out_B = "B_models_clust" + str(c) + ".rmf3"

# Arguments for do_extract_models_single_rmf:
# HA :: Dataframe object from AT.get_models_to_extract()
# file_out :: The output rmf file
AT.do_extract_models_single_rmf(
    HA,
    rmf_file_out_A,  # RMF file outputted for Sample A
    top_dir,  # Top directory containing the PMI output folders
    analys_dir,  # Analysis directory to write RMF and score files
    scores_prefix="A_models_clust" + str(c),
)  # Prefix for the scores file

AT.do_extract_models_single_rmf(
    HB, rmf_file_out_B, top_dir, analys_dir, scores_prefix="B_models_clust" + str(c)
)
exit()
