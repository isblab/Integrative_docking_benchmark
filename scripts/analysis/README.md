## Distance between crosslinked residues and percentage of crosslinks satisfied in the models

### EASAL models
Run `crosslink_distance_perc_satisfied_calculation/imp_output/calc_xlink_dist_perc_pdb.py chainA chainB crosslinks.csv`, where crosslinks.csv is the crosslink file, for each case in the easal output directories containing all the pdbs.

Run the following wrapper script for 30 cases:
```
wrapper_xlink_dist_perc_pdb.sh  
```
### IMP models
Run `crosslink_distance_perc_satisfied_calculation/imp_output/calc_xlink_dist_perc_rmf.py sampcon_0_extracted.rmf3 crosslinks.csv 32`, where sampcon_0_extracted.rmf3 is generated after analysis scripts, crosslinks.csv is the crosslink file and 32 is the max violation length for DSSO (20 for EDC), for each case.

Run the following wrapper script for 30 cases:
```
wrapper_xlink_dist_perc_rmf.sh  
```
Run scripts in `crosslink_distance_perc_satisfied_calculation/plots/` for summary (#TODO which scripts) and complex-wise plots (#TODO which scripts). 

## Comparison of crosslink distances in native structure and models

#TODO say how distances in native were calculated? did we have wrapper for it? 

Run `crosslink_distance_perc_satisfied_calculation/plots/plot_native_model_xlink_dist.py complexwise` with 'complexwise' flag to plot violin plots of the selected cases and 'summary' flag to plot the average crosslink distances in models and native structure.

## RMSD in IMP models, EASAL models and native structure
Run `structure_related/calc_and_plot_rmsd.py` to plot the minimum RMSD in IMP and EASAL models with respect to the native structure.

## Comparison of time taken in IMP and EASAL sampling
Run `time_related/compare_efficiency.py` to plot the time taken (in minutes) to sample models in each case.
