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
import numpy as np

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

total_counts = []
for frame in tqdm(range(rmf_fh.get_number_of_frames())):
    IMP.rmf.load_frame(rmf_fh, frame)
    mdl.update()

    stop_iteration = False  # Flag to stop outer loop

    chain_A_coords = np.array([IMP.core.XYZ(i).get_coordinates() for i in chain_A_in_rmf])
    chain_B_coords = np.array([IMP.core.XYZ(j).get_coordinates() for j in chain_B_in_rmf])

    from scipy.spatial import distance_matrix

    dist_matrix = distance_matrix(chain_A_coords, chain_B_coords)

    count_not_in_threshold = np.sum(dist_matrix < (threshold))
    total_counts.append(count_not_in_threshold)

    if count_not_in_threshold > 0:
        stop_iteration = True

    if not stop_iteration:
        print(f'saved_{frame}')
        IMP.rmf.save_frame(fh_out, str(frame))
