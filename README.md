# EASAL vs IMP

To compare methods for integrative docking of complexes based on crosslinking. We compare EASAL and IMP.

## To run EASAL
1. Download EASAL from `https://bitbucket.org/geoplexity/easal-dev/branch/shruthi`
2. Store all the input pdb files in `easal-dev/scripts/pdb_files`.
3. Change the following parameters for each case in `easal-dev/settings.ini`:

A. [Input]
file: name of pdb file containing two proteins
chain_A, chain_B: chain name in pdb file.

B. [Sampling]
stepSize: step size of sampling, smaller means finer (and slower) sample procedure. We have set it to 5

C. [Constraint]
activeUpperDelta: change the upper limit to 32 or 20 here based on the crosslinker length
crossLinkCount: total count of cross links in input file
crossLinkSatisfyThres: threshold of taking a configuration into account. If no configuration has this many cross links feasible, no configuration will be stored. Left to 3 for all the benchmark cases. 
crossLinks: list of all cross links in the form of {A1, B1, A2, B2, ...} corresponding to cross link file you sent me before.
 
 Run the following wrapper script to run EASAL for 30 benchmark cases:
 
```
 scripts/easal/wrapper_easal.sh
```
 
3. After you run EASAL, you'll see a result file named "A_clB_ssC.txt" where A is the input file name, B count of cross links, and C step size, corresponding to what you set in settings.ini. Move this A_clB_ssC.txt file in `easal-dev/scripts` and run `python3 result2pdb.py A B D` where A is the name of pdb, B count of cross links, D chain B

## To run IMP
1. Make directory for each complex containing the pdb file and crosslink file. For example, `imp/DSSO/1clv_2/` has 1clv.pdb and 1clv_DSSO_2.csv. Similary, for cases with EDC crosslinks.
2. Run `scripts/imp/modeling/sample_imp.py prod chainA chainB 1 DSSO 25`, where prod is for production run (10000 frames), chain A and B names, 1 is run number, crosslinker name and crosslinker length. 

Run the following wrapper script for 30 benchmark cases:

```
scripts/imp/modeling/master_script_modeling.sh
```
3. Save the density files for each complex in `imp` directory and make the direc for each case to save the analysis results. For example, `imp/DSSO_analysis/1clv_2`.
4. Run `scripts/imp/analysis/end_to_end_analysis.py` in each directory.

Run the following wrapper script for 30 benchmark cases:

```
scripts/imp/analysis/master_script_analysis.sh
```

5. We used -s -ct commands to cluster the models at a precision higher than the sampling precision.
6. Run the following command in each directory and change the cluster number based on each case.
```
 ~/imp-clean/build/setup_environment.sh python scripts/imp/analysis/extract_sampcon.py sampcon_0_extracted.rmf3 model_analysis/A_models_clust1.rmf3 sampcon/cluster.0.sample_A.txt model_analysis/B_models_clust1.rmf3 sampcon/cluster.0.sample_B.txt``` 
 
 
 
 
 
 
 
 
 
 
