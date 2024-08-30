## **Installation of integrative modeling software**
### EASAL:
* Download EASAL from `https://bitbucket.org/geoplexity/easal-dev/branch/shruthi`

### IMP:
* Compile IMP from source code. See [IMP installation](https://github.com/salilab/imp)

## **Runnning EASAL:**

### Inputs
1. Store all the input PDB files in `easal-dev/scripts/pdb_files`.
2. Change the following parameters for each case in `easal-dev/settings.ini`:

[Input]
* `file`: name of PDB file containing two proteins. e.g., 1clv.pdb
* `chain_A, chain_B`: chain name in PDB file. e.g., A,I for 1clv

[Sampling]
* `stepSize`: step size of sampling, smaller corresponds to finer sampling (and slower sampling procedure). We have set it to 5 A. If you dont get enough models satisfying crosslinks, you can consider halving the step size. 

[Constraint]
* `activeUpperDelta`: parameter to set the upper bound on the crosslink distance. We have set it to 32 or 20, based on the crosslinker length.
* `activeUpperLambda`: parameter to set the upper bound on the crosslink distance. We have set it to 0.
* `activeLowerDelta`: parameter to set the lower bound on the crosslink distance. We have set it to 0.
* `activeLowerLambda`: parameter to set the lower bound on the crosslink distance. We have set it to 2.
* `crossLinkCount`: total count of cross links in input file.
* `crossLinkSatisfyThres`:  threshold for the number of crosslinks to be satisfied by a configuration. We have set it to `n-2` where n is the number of crosslinks. Decrease this number if no configurations were found satisfying these many crosslinks. 
* `crossLinks`: list of all cross links in the form of {A1, B1, A2, B2, ...} corresponding to cross link file you sent me before.
* `smartCrossLinkMode`: number of walls between the maximum and minimum crosslink distance. We have set it to 2 (default), i.e., two walls corresponding to the upper and lower bound on crosslink distance. You can increase this value to introduce intermediate walls. 

### Run command
In `easal-dev` directory, run the following command for each case:
`'build/easal' `

You can run the following wrapper script to run EASAL for 30 benchmark cases:
 
```
 scripts/easal/wrapper_easal.sh
```

The above command for a complex will generate a text file `A_clB_ssC.txt` where A is the name of the input PDB file, B is the number of cross links, and C is the step size, corresponding to what you set in `settings.ini`. This text file contains the translations and rotations on the second protein in the complex. 

### Outputs
To write the EASAL output as PDB files:

Run `python3 easal-dev/scripts/result2pdb.py F A B D` where F is path of `A_clB_ssC.txt`, A is name of the input PDB file, B is the number of crosslinks and D is Chain B name.

This will return in `A_clB` directory containing pdb files based on the translations and rotations in the `A_clB_ssC.txt` file.

## **Runnning IMP:**

### Inputs
Make a directory for each complex containing the PDB file and crosslink file. For example, `imp/DSSO/1clv_2/` has 1clv.pdb and 1clv_DSSO_2.csv. Similarly, for cases with DMTMM crosslinks.

### Run command
1. Run the modeling script for each case as follows:
   
`~/imp-clean/build/setup_environment.sh python scripts/imp/modeling/sample_imp.py prod chainA chainB 1 DSSO 25`

where prod is for production run (10000 frames), chainA and chainB names are chain names in the PDB, 1 is the run index, DSSO is the crosslinker type and 25 is the average crosslinker length. 

You can run the following wrapper script for 30 benchmark cases:
```
scripts/imp/modeling/master_script_modeling.sh
```

### Outputs
1. Run the analysis script as `python scripts/imp/analysis/end_to_end_analysis.py` in each directory.

You can run the following wrapper script for 30 benchmark cases:
```
sh scripts/imp/analysis/master_script_analysis.sh
```

2. Subsequently, we used -s -ct commands to cluster the models at a precision worse than the sampling precision.
Run the following command in the directory containing the `sampcon` results for reclustering models at a worse precision.

```
~/imp-clean/build/setup_environment.sh python  ~/imp-clean/imp/modules/sampcon/pyext/src/exhaust.py -n 1clv_2 -a -m cpu_omp -c 0 -s -ct 12 -cc 4 -pr -d density_A_I.txt -gp -g 2.0  -sa model_analysis/A_models_clust1.txt -sb model_analysis/B_models_clust1.txt  -ra model_analysis/A_models_clust1.rmf3 -rb model_analysis/B_models_clust1.rmf3
```

3. Run the following command to get the RMF file containing all the models in the largest cluster. 

```
 ~/imp-clean/build/setup_environment.sh python scripts/imp/analysis/extract_sampcon.py sampcon_0_extracted.rmf3 model_analysis/A_models_clust1.rmf3 sampcon/cluster.0.sample_A.txt model_analysis/B_models_clust1.rmf3 sampcon/cluster.0.sample_B.txt
```

4. Finally, filter the model for clashes by removing the models having overlapping beads. Run the following script for each case:
`~/imp-clean/build/setup_environment.sh scripts/imp/analysis/filtering_by_clashes.py sampcon_0_extracted.rmf3 pdbfile chainA chainB ouputfile`
where sampcon_0_extracted.rmf3 is the RMF file conatining all the models in the largest cluser, pdbfile is the path to pdb file, chainA and chainB are chain names in the PDB, outputfile is the name of the output rmf3 file. 

The filtered RMF will be used for comparing IMP and EASAL model ensembles. 

You can run the following wrapper script for 30 benchmark cases:
``` 
scripts/imp/analysis/wrapper_filtering_models.sh
```
