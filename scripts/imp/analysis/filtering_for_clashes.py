import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.basic
import IMP.pmi.io.crosslink
import IMP.pmi.restraints.crosslinking
import os,sys,math,numpy
import ihm
import glob
from tqdm import tqdm

all_rmf = sys.argv[1]
threshold = 2.51577 * 2
pdb = sys.argv[2]
chain_A = sys.argv[3]
chain_B = sys.argv[4]

mdl = IMP.Model()
rmf_fh = RMF.open_rmf_file_read_only(all_rmf)
hier = IMP.rmf.create_hierarchies(rmf_fh, mdl)[0]

fh_out = RMF.create_rmf_file(sys.argv[5])
IMP.rmf.add_hierarchy(fh_out, hier)

# Selecting chain A and B in rmf
chain_A_in_rmf = IMP.atom.Selection(hier,resolution=1,molecule=chain_A).get_selected_particles()
chain_B_in_rmf = IMP.atom.Selection(hier,resolution=1,molecule=chain_B).get_selected_particles()


for frame in tqdm(range(rmf_fh.get_number_of_frames())):
    IMP.rmf.load_frame(rmf_fh, frame)
    mdl.update()

    stop_iteration = False  # Flag to stop outer loop

    for i in chain_A_in_rmf:
        for j in chain_B_in_rmf:
            dist = IMP.core.get_distance(IMP.core.XYZ(i), IMP.core.XYZ(j))
            if dist < threshold:
                print(dist)
                stop_iteration = True
                break

        if stop_iteration:
            break

    if not stop_iteration:
        IMP.rmf.save_frame(fh_out, str(frame))
