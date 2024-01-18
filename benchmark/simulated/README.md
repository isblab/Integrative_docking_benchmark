#Simulating crosslinks on protein complexes:

The crosslinks are simulated on the protein complexes using JWALK, downloaded from http://jwalk.ismb.lon.ac.uk/jwalk/. 
The pdbs are stored in /home/muskaan/Jwalk_2.0/Jwalk/
Run `python /home/muskaan/Jwalk_2.0/Jwalk/Jwalk.v1.1.py -max_dist 32` where -max_dist is the maximum distance between the crosslinked residues. We used 32 and 20 for DSSO and EDC crosslinker, respectively.
This will give both inter- and intra-protein crosslinks in Jwalk_results/ directory. Saved the results for EDC and DSSO in separate directories and run the following commands in them.

To filter the interprotein crosslinks:
Run `python ~/EASAL/scripts/1_processing_jwalk_output.py Jwalk_results/`

To get n crosslinks randomly:
Run `python ~/EASAL/scripts/2_getting_random_xlinks.py 1dfj_crosslink_list_interprotein_crosslinks.csv 4` where n=4

You can use master script to get different combination of number of crosslinks and crosslinker:
```
master_script.sh
```

