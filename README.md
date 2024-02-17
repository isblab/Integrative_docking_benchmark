# EASAL

To compare methods for integrative docking of complexes based on crosslinking. We compare EASAL and IMP.

Download EASAL from 

[Input]
file: name of pdb file containing two proteins
chain_A, chain_B: chain name in pdb file.
[Sampling]
stepSize: step size of sampling, smaller means finer (and slower) sample procedure.
[Constraint]
activeUpperDelta: change the upper limit to 32 or 20 here. 
crossLinkCount: total count of cross links in input file
crossLinkSatisfyThres: threshold of taking a configuration into account. If no configuration has this many cross links feasible, no configuration will be stored. Left to 3 for all the benchmark cases. 
crossLinks: list of all cross links in the form of {A1, B1, A2, B2, ...} corresponding to cross link file you sent me before.
 
After you run EASAL, you'll see a result file named "A_clB_ssC.txt" where A is the input file name, B count of cross links, and C step size, corresponding to what you set in settings.ini. Please run the attached Python script with the output in the same folder. For the script, run "python3 result2pdb.py A B D" where A is the name of pdb, B count of cross links, D name for 2nd chain (i.e. the one to be applied transformations stored in .txt file) 
