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


def get_coord_from_pdb(pdb_file, chain_B):

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    coord = []
    for model in structure:
        for chain in model:
            if chain.id == chain_B:
                for residue in chain:
                    if 'CA' in residue and residue.get_id()[0] == ' ':
                        coord.append(residue['CA'].get_coord())
    return coord


def calc_rmsd(mdl, native):
    coord1 = np.array(mdl)
    coord2 = np.array(native)
    rmsd = np.sqrt(np.mean(np.sum((coord1-coord2) ** 2, axis =1)))
    return rmsd

def get_rmsd_imp_easal(native_pdbfile, chain_A, chain_B, rmf_file, easal_output_direc):
    ## Native structure
    ts = get_coord_from_pdb(native_pdbfile, chain_B)

    ##IMP models
    mdl = IMP.Model()
    rmf_fh = RMF.open_rmf_file_read_only(rmf_file)
    hier = IMP.rmf.create_hierarchies(rmf_fh, mdl)
    rmsd_imp = []

    for frame in tqdm(range(rmf_fh.get_number_of_frames())):
        IMP.rmf.load_frame(rmf_fh, frame)
        mdl.update()

        #Get transformations for each rmf

        pdb_ca_mdl = IMP.Model()
        pdb_ca = IMP.atom.read_pdb(native_pdbfile,pdb_ca_mdl,IMP.atom.CAlphaPDBSelector())
        pdb_ca_mdl.update()

        coords_pdb_ca = {}
        coords_ccm = {}

        sel_ca_pdb = IMP.atom.Selection(pdb_ca,resolution=1,chain_id=chain_A).get_selected_particles()
        sel_ccm = IMP.atom.Selection(hier,resolution=1,molecule=chain_A).get_selected_particles()

        coords_pdb_ca[chain_A] = [IMP.core.XYZ(i).get_coordinates() for i in sel_ca_pdb]
        # print(coords_pdb_ca[chain_A])
        coords_ccm[chain_A] = [IMP.core.XYZ(i).get_coordinates() for i in sel_ccm]
        # print(coords_ccm[chain_A])
        # print(len(coords_pdb_ca), len(coords_ccm))

        _, transformation = IMP.pmi.analysis.Alignment(query=coords_pdb_ca, template=coords_ccm).align()
        # print(transformation)
        IMP.atom.transform(pdb_ca, transformation)
        # IMP.atom.write_pdb(pdb_ca,"transformed.pdb")
        transformed_pdb_B = IMP.atom.Selection(pdb_ca,resolution=1,chain_id=chain_B).get_selected_particles()
        transformed_pdb_coord_B = [IMP.core.XYZ(i).get_coordinates() for i in transformed_pdb_B]

        #Calculate rmsd with tranformed pdb
        mdl_coord = []
        c1 = IMP.atom.Selection(hierarchy=hier,
                                molecule=chain_B,
                                resolution=1).get_selected_particles()
        for x in range(len(c1)):
            mdl_coord.append(IMP.core.XYZ.get_coordinates(IMP.core.XYZ(c1[x])))
        rmsd_imp.append(calc_rmsd(mdl_coord, transformed_pdb_coord_B))
        # exit()


    # ##EASAL models
    rmsd_easal = []
    for pdb_file in os.listdir(easal_output_direc):
        if pdb_file.endswith(".pdb"):
            pdbfile_path = os.path.join(os.path.expanduser(easal_output_direc), pdb_file)
            rmsd_easal.append(calc_rmsd(get_coord_from_pdb(pdbfile_path, chain_B),ts))
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
