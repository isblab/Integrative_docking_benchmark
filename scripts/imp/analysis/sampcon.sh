#! /bin/bash

# # $imp python ~/Documents/modeling_micos_complex/micos/scripts/analysis/run_analysis_trajectories.py ../../prod_run/mic19_1full_1N_1C_independent
#
# cd model_analysis/
#
# # read summary_hdbscan_clustering.dat and save second row, first value as a cluster
# cluster=$(awk 'NR==2 {print $1}' summary_hdbscan_clustering.dat | cut -d',' -f1)
# echo "The value is: $cluster"
#
# num_models=$(awk 'NR==2 {print $1}' summary_hdbscan_clustering.dat | cut -d',' -f11)
# echo "The value is: $num_models"
#
# if [ "$num_models" -gt 30000 ]; then
#     python ~/Documents/modeling_micos_complex/micos/scripts/analysis/variable_filter_v1.py -c "$cluster" -g model_analysis/
# fi
#

$imp  python  ~/imp-clean/imp/modules/sampcon/pyext/src/exhaust.py -n 1clv -a -m cpu_omp -c 0 -cc 2 -pr -d model_analysis/density.txt -gp -g 2.0  -sa model_analysis/A_models_clust4.txt -sb model_analysis/B_models_clust4.txt  -ra model_analysis/A_models_clust4.rmf3 -rb model_analysis/B_models_clust4.rmf3

echo "Done!"
