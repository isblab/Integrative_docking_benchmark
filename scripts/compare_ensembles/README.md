## Distance between crosslinked residues and percentage of crosslinks satisfied in the models

### EASAL models
Run `crosslink_distance_perc_satisfied_calculation/imp_output/calc_xlink_dist_perc_pdb.py chainA chainB crosslinks.csv`, where crosslinks.csv is the crosslink file, for each case in the easal output directories containing all the PDBs.

Run the following wrapper script for 30 cases:
```
wrapper_xlink_dist_perc_pdb.sh  
```
### IMP models
Run `crosslink_distance_perc_satisfied_calculation/imp_output/calc_xlink_dist_perc_rmf.py sampcon_0_extracted.rmf3 crosslinks.csv 32`, where sampcon_0_extracted.rmf3 is generated after analysis scripts, crosslinks.csv is the crosslink file and 32 is the max violation length for DSSO (20 for DMTMM), for each case.

Run the following wrapper script for 30 cases:
```
wrapper_xlink_dist_perc_rmf.sh  
```
Run `crosslink_distance_perc_satisfied_calculation/plots/plot_summary_dist.py` and `crosslink_distance_perc_satisfied_calculation/plots/plot_summary_perc_satisfied` for summary.

Run `crosslink_distance_perc_satisfied_calculation/plots/plot_crosslink_distance.py` and `crosslink_distance_perc_satisfied_calculation/plots/plot_crosslink_perc_satisfied.py`for complex-wise plots. 

## Comparison of crosslink distances in native structure and models
To calculate the difference in the crosslink distance in the native structure and the models and plot it,
Run `crosslink_distance_perc_satisfied_calculation/plots/calc_and_plot_native_model_xlink_dist.py complexwise` with 'complexwise' flag to plot violin plots of the selected cases and 'summary' flag to plot the average crosslink distances in models and native structure.

## RMSD in IMP models, EASAL models and native structure
1. Run `structure_related/calc_rmsd.py native_pdb, chain_A, chain_B, rmf_file, easal_ouput_direc, outputfile` for calculating RMSD of the models to the native structure and saving it into a text file. Here, native_pdb file, names of the Chain A and Chain B, RMF file from IMP, EASAL output directory with all the models, and output file name are given as arguments.

Use `structure_related/wrapper_calc_rmsd.sh` to do the above task for all the benchmark cases.

2. Run `structure_related/plot_rmsd.py min` to plot the minimum RMSD using 'min' flag and plot RMSD in all the models using 'all' flag. 

### Superposing structures of the native structure, IMP and EASAL model

Run `structure_related/align_ccm_pdbs.py`

This will align and save the CA atoms of the specified IMP and EASAL model which has the minimum RMSD and the native structure.

## Comparison of time taken in IMP and EASAL sampling
Run `time_related/compare_efficiency.py` to plot the time taken (in minutes and CPU hours) to sample models using both methods.
