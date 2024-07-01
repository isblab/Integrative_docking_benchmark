# Integrative docking benchmark 

This repository contains a benchmark of hetero-dimeric complexes, with component atomic structures and chemical crosslinks. It is used to compare integrative modeling methods: EASAL and IMP.
**
Where is the benchmark? Creation scripts? **

## **Installation of integrative modeling software**
### EASAL:
* Download EASAL from `https://bitbucket.org/geoplexity/easal-dev/branch/shruthi`

### IMP:
* Compile IMP from source code. See [IMP installation](https://github.com/salilab/imp)

## **Runnning EASAL:**

### Inputs
1. Store all the input pdb files in `easal-dev/scripts/pdb_files`.
2. Change the following parameters for each case in `easal-dev/settings.ini`:

[Input]
* file: name of pdb file containing two proteins. e.g., 1clv.pdb
* chain_A, chain_B: chain name in pdb file. e.g., A,I for 1clv

[Sampling]
* stepSize: step size of sampling, smaller means finer (and slower) sample procedure. We have set it to 5

[Constraint]
* activeUpperDelta: change the upper limit to 32 or 20 here based on the crosslinker length
* crossLinkCount: total count of cross links in input file
* crossLinkSatisfyThres: threshold of taking a configuration into account. If no configuration has this many cross links feasible, no configuration will be stored. Left to 3 for all the benchmark cases. 
* crossLinks: list of all cross links in the form of {A1, B1, A2, B2, ...} corresponding to cross link file you sent me before.

### Run command
In `easal-dev` directory, run the following command for each case:
`'build/easal' `

You can run the following wrapper script to run EASAL for 30 benchmark cases:
 
```
 scripts/easal/wrapper_easal.sh
```

### Outputs
1. The above command will generate a text file names as `A_clB_ssC.txt` where A is the name of the input PDB file, B count of cross links, and C step size, corresponding to what you set in `settings.ini`. 
2. Move this `A_clB_ssC.txt` file to `easal-dev/scripts` 
3. Run `python3 result2pdb.py A B D` where A is the name of the input PDB file, B is the count of crosslinks, D is chain B.

## **Runnning IMP:**

### Inputs
Make a directory for each complex containing the pdb file and crosslink file. For example, `imp/DSSO/1clv_2/` has 1clv.pdb and 1clv_DSSO_2.csv. Similarly, for cases with EDC crosslinks.

### Run command
1. Run the modeling script for each case as follows:
`~/imp-clean/build/setup_environment.sh python scripts/imp/modeling/sample_imp.py prod chainA chainB 1 DSSO 25`
where prod is for production run (10000 frames), chain A and B names, 1 is run number, crosslinker name and crosslinker length. 

You can run the following wrapper script for 30 benchmark cases:
```
scripts/imp/modeling/master_script_modeling.sh
```

### Outputs
Run the analysis script as `python scripts/imp/analysis/end_to_end_analysis.py` in each directory.

You can run the following wrapper script for 30 benchmark cases:
```
sh scripts/imp/analysis/master_script_analysis.sh
```
**Note** 
We used -s -ct commands to cluster the models at a precision worse than the sampling precision.

Run the following command in the directory containing the `sampcon` results for reclustering models at a worse precision.

```
 ~/imp-clean/build/setup_environment.sh python scripts/imp/analysis/extract_sampcon.py sampcon_0_extracted.rmf3 model_analysis/A_models_clust1.rmf3 sampcon/cluster.0.sample_A.txt model_analysis/B_models_clust1.rmf3 sampcon/cluster.0.sample_B.txt
```


## **Information**
**Author(s):** Yichi Zhang, Muskaan Jindal, Shruthi Viswanath, Meera Sitharam  
**License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0
International License.  

 
 
 
 
 
 
 
 
