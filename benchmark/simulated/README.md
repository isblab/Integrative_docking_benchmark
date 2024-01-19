# Simulated crosslinks on experimental structures of protein complexes

## Structures
Obtained from cMNXL dataset (DOI: 10.1016/j.str.2018.04.016 ) which is a subset of ZDOCK benchmark (https://zlab.umassmed.edu/benchmark/). 

## Simulating crosslinks on protein complexes

The crosslinks are simulated on the protein complexes using JWALK, downloaded from http://jwalk.ismb.lon.ac.uk/jwalk/

The pdbs are stored in ~/Jwalk_2.0/Jwalk/

Run `python ~/Jwalk_2.0/Jwalk/Jwalk.v1.1.py -max_dist 20 -aa1 ASP -aa2 ASP`

 -max_dist is the maximum distance between the crosslinked residues.  We used 32 and 20 for DSSO and EDC crosslinker, respectively.

 -aa1 and -aa2 are crosslinked amino acids. For EDC, aa1 and aa2 were set to all four combinations of ASP and GLU. For DSSO aa1 = aa2 = LYS (which is the default). 

This will give both inter- and intra-protein crosslinks in `Jwalk_results/` directory. Save the results for EDC and DSSO in separate directories and run the following commands in them.

To filter the interprotein crosslinks:
Run `python ~/EASAL/scripts/1_processing_jwalk_output.py Jwalk_results/`

To get a set of `n` random crosslinks:
Run `python ~/EASAL/scripts/2_getting_random_xlinks.py 1dfj_interprotein_crosslinks.csv DSSO 4` where n=4

You can use master script `master_script.sh` to get different combination of number of crosslinks and crosslinker. 

