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
    pdb_A = IMP.atom.Selection(pdb_ca,resolution=1,chain_id=chain_A).get_selected_particles()
    pdb_B = IMP.atom.Selection(pdb_ca,resolution=1,chain_id=chain_B).get_selected_particles()

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

        # for sel1 in rmf_A:
        #     hier_A = IMP.atom.Hierarchy.setup_particle(sel1)
        #     # print(IMP.atom.get_representation(hier_A))
        #     rb_A = IMP.core.RigidBody.add_member(sel1)
        #     print(type(rb_A))
        # exit()
        #
        #
        # for sel2 in rmf_B:
        #     hier_B = IMP.atom.Hierarchy.setup_particle(sel2)
        #
        # coords_pdb[chain_A] = [IMP.core.XYZ(i).get_coordinates() for i in pdb_A]
        # coords_pdb[chain_B] = [IMP.core.XYZ(i).get_coordinates() for i in pdb_B]
        #
        # coords_rmf[chain_A] = [IMP.core.XYZ(i).get_coordinates() for i in rmf_A]
        # coords_rmf[chain_B] = [IMP.core.XYZ(i).get_coordinates() for i in rmf_B]

        #Get transformations for each rmf and calculate rmsd
        rmsd_imp.append(IMP.atom.get_pairwise_rmsd_score(pdb_A, pdb_B,rmf_A, rmf_B))
        #
        # # Aligning chain A
        # _, transformation = IMP.pmi.analysis.Alignment(query=coords_rmf, template=coords_pdb).align()
        # # print(transformation)
        #
        # # Transformation on both chains in rmf file
        # IMP.atom.transform(hier_A,transformation)
        # IMP.atom.transform(hier_B,transformation)
        #
        # # exit()
        # transformed_rmf_B = IMP.atom.Selection(hier_B,resolution=1,chain_id=chain_B).get_selected_particles()
        # print(IMP.core.RigidBody.get_coordinates(hier_B))
        #
        #
        # # RMSD calculation
        # rmsd_imp.append(IMP.atom.get_rmsd(transformed_rmf_B,pdb_B))
        # print(rmsd_imp)

        # exit()

    print('imp done')
    # ##EASAL models
    rmsd_easal = []
    for pdb_file in os.listdir(easal_output_direc):
        if pdb_file.endswith(".pdb"):
            pdbfile_path = os.path.join(os.path.expanduser(easal_output_direc), pdb_file)

            pdb = IMP.Model()
            pdbfile = IMP.atom.read_pdb(pdbfile_path,pdb,IMP.atom.CAlphaPDBSelector())
            pdb.update()

            pdb_easal_A = IMP.atom.Selection(pdbfile,resolution=1,chain_id=chain_A).get_selected_particles()
            pdb_easal_B = IMP.atom.Selection(pdbfile,resolution=1,chain_id=chain_B).get_selected_particles()

            rmsd_easal.append(IMP.atom.get_pairwise_rmsd_score(pdb_A, pdb_B, pdb_easal_A, pdb_easal_B))

    print('easal done')

    return rmsd_imp, rmsd_easal


def main():
    best_rmsd_imp, best_rmsd_easal= [],[]
    native_pdbfile, chain_A, chain_B, rmf_file, easal_output_direc, outputfile = sys.argv[1:7]
    all_rmsd_imp, all_rmsd_easal = get_rmsd_imp_easal(native_pdbfile, chain_A, chain_B, rmf_file, easal_output_direc)
    best_rmsd_imp.append(min(all_rmsd_imp))
    best_rmsd_easal.append(min(all_rmsd_easal))
    with open(outputfile, 'w') as outf:
        outf.write(f'min imp {float(best_rmsd_imp [0])}\n min easal {float(best_rmsd_easal[0])}\n rmsd imp {list(all_rmsd_imp)}\n rmsd easal {list(all_rmsd_easal)}')
    # print(best_rmsd_imp, best_rmsd_easal)

if __name__ == "__main__":
    main()
