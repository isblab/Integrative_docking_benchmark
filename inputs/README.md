Running EASAL on this benchmark

### Input PDBs
 PDBs are in  `input_to_easal/`. We are taking the chain A, keeping it fixed, and docking chain B on it.
 Take the CA coordinates of all residues for these chains from these PDBs.
 Assume each residue is a sphere of 2.8 A radius.
 Make sure your code is robust to handle alternate atom locations etc.

### Crosslinks
Crosslinks for corresponding receptor-ligand complexes are in `input_to_easal/crosslinks/

### To generate a new set of crosslinks
Use the `generate_xl.py` script to get any set of crosslinks where the CA's of the residues are within a certain distance cutoff.

E.g.


`python generate_xl.py --pdb_file 1AVX.pdb -dt 6 --receptor_chains A --ligand_chains B --total_xl 10
`
