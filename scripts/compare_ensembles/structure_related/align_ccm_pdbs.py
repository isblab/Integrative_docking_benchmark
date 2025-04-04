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

all_pdb_files = ['~/easal/pdbfile/2b42.pdb', '~/easal/pdbfile/1dfj.pdb', '~/easal/pdbfile/gcvpa_gcvpb.pdb', '~/easal/pdbfile/phes_phet.pdb']

easal_files = ['~/easal/easal_output/DMTMM/2b42_cl10/node116_flip0_1111011111_0.pdb',
                '~/easal/easal_output/DSSO/simulated/1dfj_cl3/node1_flip5_111_0.pdb',
                '~/easal/easal_output/DSSO/experimental/gcvpa_gcvpb_cl5/node30_flip6_11110_0.pdb',
                '~/easal/easal_output/DSSO/experimental/phes_phet_cl8/node271_flip4_11110111_0.pdb']

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
        frame_num = 939
        input_file = '~/easal/imp_output/DMTMM_analysis/2b42_10/sampcon_0_extracted.rmf3'

    elif '1dfj.pdb' in file:
        all_proteins = [{'E':{'E':[0,range(1,124)]}}]
        frame_num = 6941
        input_file = '~/easal/imp_output/DSSO_analysis/1dfj_3/sampcon_0_extracted.rmf3'

    elif 'gcvpa_gcvpb.pdb' in file:
        all_proteins = [{'A':{'A':[0,range(1,448)]}}]
        frame_num = 14435
        input_file = '~/easal/imp_output/DSSO_analysis/gcvpa_gcvpb_5/sampcon_0_extracted.rmf3'

    elif 'phes_phet.pdb' in file:
        all_proteins = [{'A':{'A':[0,range(1,345)]}}]
        frame_num = 17953
        input_file = '~/easal/imp_output/DSSO_analysis/phes_phet_8/sampcon_0_extracted.rmf3'

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
