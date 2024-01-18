# To extract the residues from AF2 predicted monomer strutcture which are in original pdb of the complex

from Bio import PDB
import os, sys

def extract_residues(input_pdb, output_pdb, chain_id, start_res, end_res):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', input_pdb)

    extracted_structure = PDB.Model.Model(0)

    for model in structure:
        for chain in model:
            new_chain = PDB.Chain.Chain(chain_id) # to change the chain id in output pdb
            extracted_structure.add(new_chain)
            for residue in chain:
                if start_res <= residue.id[1] <= end_res:
                    new_residue = residue.copy()
                    new_chain.add(new_residue)
                    new_residue.child_dict = {}

    io = PDB.PDBIO()
    io.set_structure(extracted_structure)
    io.save(output_pdb)


input_pdb = sys.argv[1]
chain_id = sys.argv[2]
start_res = int(sys.argv[3])
end_res = int(sys.argv[4])
output_pdb = f"{os.path.splitext(input_pdb)[0].split('_')[0]}_{chain_id}.pdb"

extract_residues(input_pdb, output_pdb, chain_id, start_res, end_res)
