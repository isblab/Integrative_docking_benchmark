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
* `stepSize`: step size of sampling, smaller means finer sampling (and slower sampling procedure). We have set it to 5 A. If you dont get enough models satisfying crosslinks, you can consider halving the step size. 

[Constraint]
* `activeUpperDelta`: change the upper limit to 32 or 20 here based on the crosslinker length.
* `crossLinkCount`: total count of cross links in input file.
* `crossLinkSatisfyThres`: threshold of taking a configuration into account. If no configuration has this many cross links feasible, no configuration will be stored. Left to 3 for all the benchmark cases. 
* `crossLinks`: list of all cross links in the form of {A1, B1, A2, B2, ...} corresponding to cross link file you sent me before.
* `activeLowerLambda`: lower bound on the crosslink distance.
* `smartCrossLinkMode`: number of intermediate walls between the maximum and minimum crosslink distance.

### Run command
In `easal-dev` directory, run the following command for each case:
`'build/easal' `

You can run the following wrapper script to run EASAL for 30 benchmark cases:
 
```
 scripts/easal/wrapper_easal.sh
```

### Outputs
1. The above command will generate a text file names as `A_clB_ssC.txt` where A is the name of the input PDB file, B is the number of cross links, and C is the step size, corresponding to what you set in `settings.ini`. 
2. Run `python3 result2pdb.py F A B D` where F is path of `A_clB_ssC.txt`, A is name of the input PDB file, B is the number of crosslinks and D is Chain B name.

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

2. Run the following command to get the RMF file containing all the models in the largest cluster. This will be used for comparing IMP and EASAL model ensembles. 

```
 ~/imp-clean/build/setup_environment.sh python scripts/imp/analysis/extract_sampcon.py sampcon_0_extracted.rmf3 model_analysis/A_models_clust1.rmf3 sampcon/cluster.0.sample_A.txt model_analysis/B_models_clust1.rmf3 sampcon/cluster.0.sample_B.txt
```

**Note** 
We used -s -ct commands to cluster the models at a precision worse than the sampling precision.
Run the following command in the directory containing the `sampcon` results for reclustering models at a worse precision.

```
~/imp-clean/build/setup_environment.sh python  ~/imp-clean/imp/modules/sampcon/pyext/src/exhaust.py -n 1clv_2 -a -m cpu_omp -c 0 -s -ct 12 -cc 4 -pr -d density_A_I.txt -gp -g 2.0  -sa model_analysis/A_models_clust1.txt -sb model_analysis/B_models_clust1.txt  -ra model_analysis/A_models_clust1.rmf3 -rb model_analysis/B_models_clust1.rmf3
```
