from Bio.PDB import PDBParser
import os, sys
import IMP
import RMF
import IMP.rmf
import IMP.core
import IMP.atom
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
                        # print(residue)
                        coord.append(residue['CA'].get_coord())
    return coord


def calc_rmsd(mdl, native):
    coord1 = np.array(mdl)
    coord2 = np.array(native)
    rmsd = np.sqrt(np.mean(np.sum((coord1-coord2) ** 2, axis =1)))
    return rmsd

def get_rmsd_imp_easal(native_pdbfile, chain_B, rmf_file, easal_output_direc):
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
        mdl_coord = []
        c1 = IMP.atom.Selection(hierarchy=hier,
                                molecule=chain_B,
                                resolution=1).get_selected_particles()
        for x in range(len(c1)):
            mdl_coord.append(IMP.core.XYZ.get_coordinates(IMP.core.XYZ(c1[x])))
        rmsd_imp.append(calc_rmsd(mdl_coord, ts))

    ##EASAL models
    rmsd_easal = []
    for pdb_file in os.listdir(easal_output_direc):
        if pdb_file.endswith(".pdb"):
            pdbfile_path = os.path.join(os.path.expanduser(easal_output_direc), pdb_file)
            rmsd_easal.append(calc_rmsd(get_coord_from_pdb(pdbfile_path, chain_B),ts))

    return rmsd_imp, rmsd_easal


def main():
    best_rmsd_imp, best_rmsd_easal= [],[]

    input_cases = [ "1dfj_DSSO_3", "1clv_DSSO_2", "1kxp_DSSO_4", "1r0r_DSSO_3", "2ayo_DSSO_4", "2b42_DSSO_5", "2hle_DSSO_5",
        "1dfj_EDC_4", "1clv_EDC_8", "1kxp_EDC_7", "1r0r_EDC_6", "2ayo_EDC_5", "2b42_EDC_10", "2hle_EDC_9",
        "1dfj_DSSO_9", "1clv_DSSO_6", "1kxp_DSSO_7", "1r0r_DSSO_7", "2ayo_DSSO_8", "2b42_DSSO_10", "2hle_DSSO_10",
        "1dfj_DSSO_12", "1kxp_DSSO_11", "2ayo_DSSO_13", "2hle_DSSO_14",
        "gata_gatc_DSSO_3", "gcvpa_gcvpb_DSSO_5","roca_putc_DSSO_2", "sucd_succ_DSSO_4", "phes_phet_DSSO_8"]

    for case in input_cases:
        # print(len(case), case)
        if len(case) < 15:
            if 'DSSO' in case:
                rmf_file = '/home/muskaan/easal/imp_output/DSSO_analysis/' + case.split('DSSO')[0] + case.split('_')[-1]+ '/sampcon_0_extracted.rmf3'
                easal_output_direc = '/home/muskaan/easal/easal_output/DSSO/simulated/'+case.split('DSSO')[0]+ 'cl'+ case.split('_')[-1]+'/'
            elif 'EDC' in case:
                rmf_file = '/home/muskaan/easal/imp_output/EDC_analysis/' + case.split('EDC')[0] + case.split('_')[-1]+ '/sampcon_0_extracted.rmf3'
                easal_output_direc = '/home/muskaan/easal/easal_output/EDC/'+case.split('EDC')[0]+ 'cl'+ case.split('_')[-1]+'/'
        else:
            rmf_file = '/home/muskaan/easal/imp_output/DSSO_analysis/' + case.split('DSSO')[0] + case.split('_')[-1]+ '/sampcon_0_extracted.rmf3'
            easal_output_direc = '/home/muskaan/easal/easal_output/DSSO/experimental/'+case.split('DSSO')[0]+ 'cl'+ case.split('_')[-1]+'/'

        if '1dfj' in case or '1clv' in case or '1r0r' in case:
            native_pdbfile = f'/home/muskaan/easal-dev/scripts/pdbfile/{case[0:4]}.pdb'
            chain_B = 'I'

        elif '1kxp' in case:
            native_pdbfile = '/home/muskaan/easal-dev/scripts/pdbfile/1kxp.pdb'
            chain_B = 'D'

        elif '2hle' in case or '2b42' in case or '2ayo' in case:
            native_pdbfile = f'/home/muskaan/easal-dev/scripts/pdbfile/{case[0:4]}.pdb'
            chain_B = 'B'
        else:
            native_pdbfile = '/home/muskaan/easal-dev/scripts/pdbfile/'+case.split('_DSSO')[0]+'.pdb'
            chain_B = 'B'

        # print(native_pdbfile, rmf_file, easal_output_direc)
        all_rmsd_imp, all_rmsd_easal = get_rmsd_imp_easal(native_pdbfile, chain_B, rmf_file, easal_output_direc)
        best_rmsd_imp.append(min(all_rmsd_imp))
        best_rmsd_easal.append(min(all_rmsd_easal))

    print(best_rmsd_imp, best_rmsd_easal)
    plt.scatter(best_rmsd_imp, best_rmsd_easal)
    plt.xlabel('Minimum RMSD in IMP Ensemble (Å)', fontsize=14)
    plt.ylabel('Minimum RMSD in EASAL Ensemble (Å)', fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.xlim(0, 80)
    plt.ylim(0, 80)
    plt.savefig('/home/muskaan/easal/plots/best_rmsd.png',dpi=600)
    plt.show()
    # plot_rmsd_in_all_models(case, all_rmsd_imp, all_rmsd_easal)


if __name__ == "__main__":
    main()
