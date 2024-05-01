from Bio.PDB import PDBParser
import os, sys
import IMP
import RMF
import IMP.rmf
import IMP.core
import IMP.atom
import IMP.pmi.analysis
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

def get_rmsd_imp_easal(native_pdbfile, chain_A, chain_B, rmf_file, easal_output_direc):

    #Native structure
    pdb_ca_mdl = IMP.Model()
    pdb_ca = IMP.atom.read_pdb(native_pdbfile,pdb_ca_mdl,IMP.atom.CAlphaPDBSelector())
    pdb_ca_mdl.update()

    # Selecting chain A and B in pdb
    pdb_native_A = IMP.atom.Selection(pdb_ca,resolution=1,chain_id=chain_A).get_selected_particles()
    pdb_native_B = IMP.atom.Selection(pdb_ca,resolution=1,chain_id=chain_B).get_selected_particles()

    ##IMP models
    mdl = IMP.Model()
    rmf_fh = RMF.open_rmf_file_read_only(rmf_file)
    hier = IMP.rmf.create_hierarchies(rmf_fh, mdl)
    rmsd_imp = []

    for frame in tqdm(range(rmf_fh.get_number_of_frames())):
        IMP.rmf.load_frame(rmf_fh, frame)
        mdl.update()

        # Selecting chain A and B in rmf
        rmf_A = IMP.atom.Selection(hier,resolution=1,molecule=chain_A).get_selected_particles()
        rmf_B = IMP.atom.Selection(hier,resolution=1,molecule=chain_B).get_selected_particles()

        #Get transformations for each rmf and calculate rmsd
        rmsd_imp.append(IMP.atom.get_pairwise_rmsd_score(pdb_native_A, rmf_B ,rmf_A, pdb_native_B))

    print('imp done')
    # ##EASAL models
    rmsd_easal = []
    for pdb_file in os.listdir(easal_output_direc):
        if pdb_file.endswith(".pdb"):
            # print(pdb_file)
            pdbfile_path = os.path.join(os.path.expanduser(easal_output_direc), pdb_file)

            pdb = IMP.Model()
            pdbfile = IMP.atom.read_pdb(pdbfile_path,pdb,IMP.atom.CAlphaPDBSelector())
            pdb.update()

            pdb_easal_A = IMP.atom.Selection(pdbfile,resolution=1,chain_id=chain_A).get_selected_particles()
            pdb_easal_B = IMP.atom.Selection(pdbfile,resolution=1,chain_id=chain_B).get_selected_particles()

            rmsd_easal.append(IMP.atom.get_pairwise_rmsd_score(pdb_native_A, pdb_easal_B, pdb_easal_A, pdb_native_B))

    print('easal done')

    return rmsd_imp, rmsd_easal


def main():
    best_rmsd_imp, best_rmsd_easal= [],[]
    native_pdbfile, chain_A, chain_B, rmf_file, easal_output_direc, outputfile = sys.argv[1:7]
    all_rmsd_imp, all_rmsd_easal = get_rmsd_imp_easal(native_pdbfile, chain_A, chain_B, rmf_file, easal_output_direc)
    best_rmsd_imp.append(min(all_rmsd_imp))
    best_rmsd_easal.append(min(all_rmsd_easal))
    # with open(outputfile, 'w') as outf:
    #     outf.write(f'min imp {float(best_rmsd_imp [0])}\n min easal {float(best_rmsd_easal[0])}\n rmsd imp {list(all_rmsd_imp)}\n rmsd easal {list(all_rmsd_easal)}')
    # print(best_rmsd_imp, best_rmsd_easal)
    print('Frame with min rmsd', all_rmsd_imp.index(min(all_rmsd_imp)))
    print('PDB file with min rmsd',all_rmsd_easal.index(min(all_rmsd_easal)))
if __name__ == "__main__":
    main()
