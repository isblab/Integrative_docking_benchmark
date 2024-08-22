import os, sys
import IMP
import RMF
import IMP.core
import IMP.rmf
import IMP.pmi.analysis
import glob
from Bio.PDB import *


###################################################################################################
############################################# Inputs ##############################################
###################################################################################################

all_pdb_files = ['/home/muskaan/easal/pdbfile/2b42.pdb', '/home/muskaan/easal/pdbfile/1dfj.pdb', '/home/muskaan/easal/pdbfile/roca_putc.pdb', '/home/muskaan/easal/pdbfile/phes_phet.pdb']

easal_files = ['/home/muskaan/easal/easal_output/EDC/2b42_cl10/node71_flip0_1111111110_0.pdb',
                '/home/muskaan/easal/easal_output/DSSO/simulated/1dfj_cl3/node1_flip5_111_0.pdb',
                '/home/muskaan/easal/easal_output/DSSO/experimental/roca_putc_cl2/node6_flip5_11_0.pdb',
                '/home/muskaan/easal/easal_output/DSSO/experimental/phes_phet_cl8/node226_flip2_11110111_0.pdb']

for pdb in easal_files:
    mdl = IMP.Model()
    pdb_ca_easal = IMP.atom.read_pdb(pdb,mdl,IMP.atom.CAlphaPDBSelector())
    mdl.update()

    IMP.atom.write_pdb(pdb_ca_easal, f"./ca_easal_{pdb.split('/')[-1]}")

# The all proteins list has the following architecture:
# [{protein:{chain_id,residue range}}, {protein:{chain_id,residue range}]
# The order of entries in the offset list must be the same as that in the pdb_files list

for file in all_pdb_files:
    pdb_files = [file]

    if '2b42.pdb' in file:
        all_proteins = [{'A':{'A':[0,range(1,42)]},'A':{'A':[0,range(47,69)]},'A':{'A':[0,range(78,263)]},'A':{'A':[0,range(266,382)]}}]
        frame_num = 8
        input_file = '/home/muskaan/easal/imp_output/EDC_analysis/2b42_10/sampcon_0_extracted.rmf3'

    elif '1dfj.pdb' in file:
        all_proteins = [{'E':{'E':[0,range(1,124)]}, 'I':{'I':[0,range(1,456)]}}]
        frame_num = 115
        input_file = '/home/muskaan/easal/imp_output/DSSO_analysis/1dfj_3/sampcon_0_extracted.rmf3'

    elif 'roca_putc.pdb' in file:
        all_proteins = [{'A':{'A':[0,range(1,516)]}}]
        frame_num = 26
        input_file = '/home/muskaan/easal/imp_output/DSSO_analysis/roca_putc_2/sampcon_0_extracted.rmf3'

    elif 'phes_phet.pdb' in file:
        all_proteins = [{'A':{'A':[0,range(1,345)]}}]
        frame_num = 8
        input_file = '/home/muskaan/easal/imp_output/DSSO_analysis/phes_phet_8/sampcon_0_extracted.rmf3'

###################################################################################################
##################################### Get transformations #########################################
###################################################################################################

    for file_index in range(len(pdb_files)):
        print(f"Aligning: {pdb_files[file_index]}")
        pdb_file = pdb_files[file_index]
        proteins = all_proteins[file_index]
        ccm_mdl = IMP.Model()
        ccm = RMF.open_rmf_file_read_only(input_file)
        hier = IMP.rmf.create_hierarchies(ccm, ccm_mdl)[0]
        IMP.rmf.load_frame(ccm, frame_num)
        ccm_mdl.update()
        pdb_ca_mdl = IMP.Model()
        pdb_ca = IMP.atom.read_pdb(pdb_file,pdb_ca_mdl,IMP.atom.CAlphaPDBSelector())
        pdb_ca_mdl.update()

        coords_pdb_ca = {}
        coords_ccm = {}

        for prot in proteins.keys():
            print(proteins)
            for chain_id in proteins[prot]:
                protein_name = prot
                print(type(IMP.atom.get_leaves(hier)))
                sel_ca_pdb = IMP.atom.Selection(pdb_ca,resolution=1,chain_id=chain_id,residue_indexes=[i for i in proteins[prot][chain_id][1]]).get_selected_particles()
                sel_ccm = IMP.atom.Selection(hier,resolution=1,molecule=protein_name,copy_index=proteins[prot][chain_id][0],residue_indexes=[i for i in proteins[prot][chain_id][1]]).get_selected_particles()

                # Remove coarse grained beads
                new_ccm_sel = []
                for selection in sel_ccm:
                    if not IMP.atom.Fragment.get_is_setup(selection):
                        new_ccm_sel.append(selection)

                coords_pdb_ca[protein_name] = [IMP.core.XYZ(i).get_coordinates() for i in sel_ca_pdb]
                coords_ccm[protein_name] = [IMP.core.XYZ(i).get_coordinates() for i in new_ccm_sel]
                # print(len(coords_pdb_ca[protein_name]),len(coords_ccm[protein_name]))
        _, transformation = IMP.pmi.analysis.Alignment(query=coords_ccm, template=coords_pdb_ca).align()
        print(transformation)


        ###################################################################################################
        #################################### Transform and write PDB ######################################
        ###################################################################################################
        IMP.atom.transform(IMP.atom.Hierarchy(hier), transformation)
        IMP.atom.write_pdb(IMP.atom.Hierarchy(hier), f"./aligned_imp_mdl_{pdb_file.split('/')[-1]}")
        IMP.atom.write_pdb(pdb_ca, f"./ca_native_{pdb_file.split('/')[-1]}")
