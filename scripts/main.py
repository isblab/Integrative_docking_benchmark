from Bio.PDB import PDBParser, Chain
from Bio.PDB.PDBIO import PDBIO
import math
import random
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--mode", type=int, default=0, help='0 for creating input, 1 for handling EASAL output')
parser.add_argument(
    "--pdb_file",
    type=str,
    required=True,
    help="Input PDB file of complex for which need crosslinks between receptor and ligand chains",
)
parser.add_argument(
    "--receptor_chains",
    type=str,
    default="A",
    help="String of receptor chains, e.g. AB",
)
parser.add_argument(
    "--ligand_chains",
    type=str,
    default="B",
    help="String of ligand chains, e.g. CD"
)
parser.add_argument(
    "--dt",
    type=float,
    required=True,
    help="Max distance between crosslinked residues defined by linker length",
)
parser.add_argument(
    "--fp",
    type=int,
    default=0,
    help="Number of false positives in the total set of crosslinks",
)
parser.add_argument(
    "--total_xl", type=int, required=True, help="Total number of crosslinks"
)
args = parser.parse_args()


def get_distance(c1, c2):
    raw_d = (c1[0] - c2[0]) ** 2 + (c1[1] - c2[1]) ** 2 + (c1[2] - c2[2]) ** 2
    d = math.sqrt(raw_d)
    return d


def select_random_entries(input_list, num_out):
    if num_out > len(input_list):
        raise Exception(
            f"Not enough crosslinks within the distance threshold found: {num_out}>{len(input_list)}"
        )
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


def get_calphas_for_a_chain(chain):
    calphas = []
    residues = [residue for residue in chain.get_residues()]
    for res in residues:
        calphas.append(CAlphaAtom(chain=chain, residue=res))
    return calphas


class CAlphaAtom:
    def __init__(self, chain, residue) -> None:
        self.chain = chain
        self.residue = residue
        self.coords = self.get_ca_coords()

    def get_ca_coords(self):
        coords = None
        for atom in self.residue.get_atoms():
            if atom.get_name() == "CA":
                coords = atom.get_coord()
        return coords


##########################################################################################################################################
################################################################### Main #################################################################
##########################################################################################################################################


if __name__ == '__main__':
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", args.pdb_file)

    if args.mode == 0:
        rchains = [Chain.Chain(ch_name) for ch_name in args.receptor_chains]
        lchains = [Chain.Chain(ch_name) for ch_name in args.ligand_chains]

        rchain_calphas = []
        lchain_calphas = []

        for chain in structure[0].get_chains():  # Assuming only one model in the PDB
            if chain in rchains:
                print(f"Receptor chain found {chain}")
                rchain_calphas.append(get_calphas_for_a_chain(chain))
            elif chain in lchains:
                print(f"Ligand chain found {chain}")
                lchain_calphas.append(get_calphas_for_a_chain(chain))
        with open(args.pdb_file[:-4] + '_A.pdb', 'w') as fileA:
            for chain in rchain_calphas:
                i = 0
                for res in chain:
                    fileA.write('ATOM ' + str(res.residue.id[1]) + ' CA A ' + res.residue.resname + ' ' + str(
                        res.coords[0]) + ' ' + str(
                        res.coords[1]) + ' ' + str(res.coords[2]) + ' 2.8\n')
                    i += 1

        with open(args.pdb_file[:-4] + '_B.pdb', 'w') as fileB:
            for chain in lchain_calphas:
                i = 0
                for res in chain:
                    fileB.write(
                        'ATOM ' + str(i) + ' CA B ' + res.residue.resname + ' ' + str(res.coords[0]) + ' ' + str(
                            res.coords[1]) + ' ' + str(res.coords[2]) + ' 2.8\n')
                    i += 1

        crosslinks = []
        noncrosslinks = []

        for rchain in rchain_calphas:
            for lchain in lchain_calphas:
                for rres in rchain:
                    for lres in lchain:
                        distance = get_distance(rres.coords, lres.coords)  # TODO test

                        if distance <= args.dt:
                            crosslinks.append((rres, lres, distance))
                        else:
                            noncrosslinks.append((rres, lres, distance))

        # STEP 3. Select a random subset of entries and save in file
        num_false_xl = args.fp
        num_true_xl = int(args.total_xl - args.fp)

        all_xls = select_random_entries(crosslinks, num_true_xl)
        false_xls = select_random_entries(noncrosslinks, num_false_xl)
        all_xls.extend(false_xls)
        random.shuffle(all_xls)

        with open(f"xl_{args.pdb_file.split('/')[-1]}_d{args.dt}.txt", "w") as outf:
            outf.write("num1, res1, prot1, num2, res2, prot2\n")
            for lnk in all_xls:
                rca, lca = lnk[0], lnk[1]

                outf.write(
                    f"{rchain_calphas[0].index(rca)}, {rca.residue.get_id()[1]}, {rca.chain.get_id()}, "
                    f"{lchain_calphas[0].index(lca)}, {lca.residue.get_id()[1]}, {lca.chain.get_id()}, {lnk[2]}\n"
                )

    elif args.mode == 1:
        for chain in structure.get_chains():
            if chain.get_id() == args.receptor_chains:
                molA = chain
            elif chain.get_id() == args.ligand_chains:
                molB = chain
        io = PDBIO()

        with open('allOris.txt') as allOris:
            for modelcount,line in enumerate(allOris):
                cur_line = list(map(float, line.split()))
                flip_no = int(cur_line[0])
                ori_no = int(cur_line[1])
                trans = np.array(cur_line[2:5])
                rot = np.array(cur_line[5:14]).reshape((3, 3))
                print(modelcount)
                for atom in molB.get_atoms():
                    # print(atom.coord)
                    # print(rot)
                    atom.set_coord(rot @ atom.coord + trans) # multiply rotation matrix and add translation
                    # print(atom.get_coord())
                    # exit()
                output_pdb_filename = f'{args.pdb_file}_model_{modelcount}.pdb'
                with open(output_pdb_filename, 'w') as output_pdb:
                    io.set_structure(structure)
                    io.save(output_pdb)

                print(f'Saved PDB file: {output_pdb_filename}')

                # pack up molA and molB and output as 1avx_modelcount.pdb file
                # save_pdb(structure)

                #TODO make sure residue numbers maintained, all atoms are there
