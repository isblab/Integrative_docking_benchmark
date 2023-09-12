import os, sys
import string
from Bio.PDB import PDBParser, PDBIO
import math
import random
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--pdb_file", type=str, required=True, help="Input PDB file of complex for which need crosslinks between receptor and ligand chains")
parser.add_argument("--receptor_chains",type=str, default="A",help="String of receptor chains, e.g. AB")
parser.add_argument("--ligand_chains",type=str, default="B",help="String of ligand chains, e.g. CD")
parser.add_argument("--dt", type=float, required=True, help="Max distance between crosslinked residues defined by linker length")
parser.add_argument("--fp", type=int, default=0, help="Number of false positives in the total set of crosslinks")
parser.add_argument("--total_xl", type=int, required=True, help="Total number of crosslinks")
args = parser.parse_args()


def get_distance(c1, c2):
    raw_d = (
        (c1[0] - c2[0]) ** 2
        + (c1[1] - c2[1]) ** 2
        + (c1[2] - c2[2]) ** 2
    )
    d = math.sqrt(raw_d)
    return d


def select_random_entries(input_list, num_out):
    random.shuffle(input_list)
    outlist = []
    while len(outlist) < num_out:
        choice = random.choice(input_list)
        if choice not in outlist:
            outlist.append(choice)
    return outlist

def create_pseudobond_file(fname, xl_list, color):
    model_num = 0
    atom_type = "ca"
    with open(fname, "w") as outf:
        for link in xl_list:
            chain1 = link.split(",")[0][-1].strip()
            res1 = int(link.split(",")[1])
            chain2 = link.split(",")[2][-1].strip()
            res2 = int(link.split(",")[3])
            pseudobond = f"#{model_num}:{res1}.{chain1}@{atom_type} #{model_num}:{res2}.{chain2}@{atom_type} {color}"
            outf.write(pseudobond + "\n")


##########################################################################################################################################
################################################################### Main #################################################################
##########################################################################################################################################


parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", args.pdb_file)

rchains = args.receptor_chains
lchains = args.ligand_chains

rcoords = [] # list of atomic coordinates of CA's of receptor chain.
lcoords =[] # list of atomic coordinates of CA's of ligand chain.

# Format: each of the above lists is a list of tuples. Each tuple has (residue number, chain, xyz coords of CA of residue)

for model in structure[0]: # Assuming only one model in the PDB
    for chain in model:

        if chain not in rchains and chain not in lchains: #TODO test if string membership is recorded this way
            continue

        for residue in chain:
            if residue.has_id('CA'):
                atom = residue['CA']

                # get residue number
                # get chain name
                # get coordinates
                atom_entry = (atom.get_full_id()[3][1],atom.get_full_id()[2],atom.get_coord()) #TODO test

                # If that does not work try the below for the chain name
                # str(atom.get_parent().get_parent().get_id())

                if chain in rchains:
                    rcoords.append(atom_entry)
                elif chain in lchains:
                    lcoords.append(atom_entry)

crosslinks =[]
noncrosslinks = []

for ri,rres in enumerate(rcoords):
    for li,lres in enumerate(lcoords):

        distance = get_distance(rres[2], lres[2]) # TODO test

        if distance <= args.dt:
            crosslinks.append((ri,li))    # enough to store the indices of the lists
        else:
            noncrosslinks.append((ri,li))

## STEP 3. Select a random subset of entries and save in file
num_false_xl = args.fp
num_true_xl = int(args.total_xl - args.fp)

all_xls = select_random_entries(crosslinks, num_true_xl)
false_xls = select_random_entries(noncrosslinks, num_false_xl)
all_xls.extend(false_xls)
random.shuffle(all_xls)

with open('xl_'+args.pdb_file+'_d'+int(args.dt)+'.txt', "w") as outf:

        for lnk in all_xls:
            ri,li=lnk[0],lnk[1]
            outf.write(f"{rcoords[ri][1]},{rcoords[ri][0]},{lcoords[li][1]},{lcoords[li][0]}" + "\n")
